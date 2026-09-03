/*********************************************************************/
/* File:   device_sparsecholesky.cpp                                 */
/* Date:   4. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_SPARSECHOLESKY_CPP

#include <la.hpp>
#include <gpukernel.hpp>
#include <climits>

namespace ngla
{
  using namespace ngs_gpu;

  namespace
  {
    string SubstituteScal (string src, const string & value)
    {
      const string token = "SCAL";
      for (size_t p = src.find(token); p != string::npos; p = src.find(token, p))
        src.replace (p, token.size(), value);
      return src;
    }

    // MicroTask::BT
    constexpr int L_BLOCK = 0, B_BLOCK = 1, LB_BLOCK = 2;

    /*
      The persistent solves: a group is GROUP_SIZE_X lanes (the simd
      width) times GROUP_SIZE_Y independent rows. Each row takes jobs
      from cnt, waits until incomingdep[job] has been counted down to
      zero by the jobs it depends on, does the work of the micro-task
      and decrements the counters of the dependent jobs.

      hy is passed twice: plainly for the ordinary loads and stores,
      and as an atomic buffer for the accumulations into rows that
      other jobs update at the same time.
    */
    const char * kernel_source = R"RAW(

      KERNEL(chol_reorder, GLOBAL_IN(SCAL,src), GLOBAL(SCAL,dst), GLOBAL_IN(int,order), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n && order[i] != -1) dst[order[i]] = src[i];
      }

      KERNEL(chol_reorder_add, GLOBAL_IN(SCAL,src), GLOBAL(SCAL,dst), GLOBAL_IN(int,order),
                               VALUE(SCAL,s), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n && order[i] != -1) dst[i] += s*src[order[i]];
      }

      KERNEL(chol_multdiag, GLOBAL_IN(SCAL,diag), GLOBAL(SCAL,vec), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n) vec[i] *= diag[i];
      }

      // working counters back to their initial values, job counters to zero
      KERNEL(chol_reset, GLOBAL_IN(int,dep0), GLOBAL(int,dep), GLOBAL_IN(int,dep0t), GLOBAL(int,dept),
                         GLOBAL(int,cnt), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n) { dep[i] = dep0[i]; dept[i] = dep0t[i]; }
        if (i < 2) cnt[i] = 0;
      }


      #define CHOL_TASK_SETUP                                                  \
        int blocknr  = microtasks[4*myjob];                                    \
        int type     = microtasks[4*myjob+1];                                  \
        int bblock   = microtasks[4*myjob+2];                                  \
        int nbblocks = microtasks[4*myjob+3];                                  \
        int rfirst = blocks[blocknr], rnext = blocks[blocknr+1];               \
        int rsize = rnext - rfirst;                                            \
        int base = firstinrow_ri[rfirst] + rsize - 1;                          \
        int ext_size = firstinrow[rfirst+1] - firstinrow[rfirst] - rsize + 1;  \
        int mfirst = int((long)ext_size * bblock / nbblocks);                  \
        int mnext  = int((long)ext_size * (bblock+1) / nbblocks);

      KERNEL(chol_solveL, GLOBAL_IN(int,depidx), GLOBAL_IN(int,depdata), GLOBAL_ATOMIC(int,incomingdep),
                          GLOBAL(SCAL,hy), GLOBAL_ATOMIC(SCAL,hy_atomic), GLOBAL_ATOMIC(int,cnt),
                          GLOBAL_IN(int,microtasks), GLOBAL_IN(int,blocks), GLOBAL_IN(int,rowindex2),
                          GLOBAL_IN(int,firstinrow_ri), GLOBAL_IN(int,firstinrow), GLOBAL_IN(SCAL,lfact),
                          VALUE(int,njobs))
      {
        SHARED(int, myjobs, 16);
        int lane = int(LOCAL_ID_X), row = int(LOCAL_ID_Y);
        int lanes = int(GROUP_SIZE_X);

        while (true)
          {
            if (lane == 0) myjobs[row] = ATOMIC_ADD(&cnt[0], 1);
            SIMD_BARRIER();
            int myjob = myjobs[row];
            if (myjob >= njobs) break;

            if (lane == 0)
              while (ATOMIC_LOAD(&incomingdep[myjob]) != 0) { }
            SIMD_BARRIER();
            DEVICE_FENCE();

            CHOL_TASK_SETUP

            if (type == 0 || type == 2)     // L_BLOCK, LB_BLOCK
              for (int i = rfirst; i < rnext-1; i++)
                {
                  // lane reading hy[i+1] next iteration waits for the lane updating it now
                  SIMD_BARRIER();
                  SCAL hyi = hy[i];
                  int off = firstinrow[i] - i - 1;
                  for (int j = lane + rfirst; j < rnext; j += lanes)
                    if (j > i) hy[j] -= lfact[off+j] * hyi;
                  SIMD_BARRIER();
                }
            SIMD_BARRIER();   // B part reads the hy written by the L part

            if ((type == 1 || type == 2) && ext_size != 0)    // B_BLOCK, LB_BLOCK
              for (int j = lane; j < mnext-mfirst; j += lanes)
                {
                  SCAL temp = 0;
                  for (int i = rfirst; i < rnext; i++)
                    temp += lfact[firstinrow[i] + rnext-i-1 + mfirst + j] * hy[i];
                  ATOMIC_ADD(&hy_atomic[rowindex2[base + mfirst + j]], -temp);
                }

            SIMD_BARRIER();
            DEVICE_FENCE();
            if (lane == 0)
              for (int k = depidx[myjob]; k < depidx[myjob+1]; k++)
                ATOMIC_ADD(&incomingdep[depdata[k]], -1);
          }
      }


      KERNEL(chol_solveLT, GLOBAL_IN(int,depidx), GLOBAL_IN(int,depdata), GLOBAL_ATOMIC(int,incomingdep),
                           GLOBAL(SCAL,hy), GLOBAL_ATOMIC(SCAL,hy_atomic), GLOBAL_ATOMIC(int,cnt),
                           GLOBAL_IN(int,microtasks), GLOBAL_IN(int,blocks), GLOBAL_IN(int,rowindex2),
                           GLOBAL_IN(int,firstinrow_ri), GLOBAL_IN(int,firstinrow), GLOBAL_IN(SCAL,lfact),
                           VALUE(int,njobs))
      {
        SHARED(int, myjobs, 16);
        int lane = int(LOCAL_ID_X), row = int(LOCAL_ID_Y);
        int lanes = int(GROUP_SIZE_X);

        while (true)
          {
            if (lane == 0) myjobs[row] = ATOMIC_ADD(&cnt[1], 1);
            SIMD_BARRIER();
            int myjob = njobs-1 - myjobs[row];
            if (myjob < 0) break;

            if (lane == 0)
              while (ATOMIC_LOAD(&incomingdep[myjob]) != 0) { }
            SIMD_BARRIER();
            DEVICE_FENCE();

            CHOL_TASK_SETUP

            if ((type == 1 || type == 2) && ext_size != 0)    // B_BLOCK, LB_BLOCK
              for (int i = rfirst + lane; i < rnext; i += lanes)
                {
                  int first = firstinrow[i] + rnext-i-1 + mfirst;
                  SCAL val = 0;
                  for (int j = 0; j < mnext-mfirst; j++)
                    val += lfact[first + j] * hy[rowindex2[base + mfirst + j]];
                  ATOMIC_ADD(&hy_atomic[i], -val);
                }
            SIMD_BARRIER();   // L part reads the hy updated by the B part

            if (type == 0 || type == 2)     // L_BLOCK, LB_BLOCK
              for (int j = rsize-1; j >= 0; j--)
                {
                  SCAL hyj = hy[rfirst+j];
                  for (int i = lane; i < j; i += lanes)
                    hy[rfirst+i] -= lfact[firstinrow[rfirst+i] - i - 1 + j] * hyj;
                  SIMD_BARRIER();   // hy[rfirst+j-1] is read next iteration
                }

            SIMD_BARRIER();
            DEVICE_FENCE();
            if (lane == 0)
              for (int k = depidx[myjob]; k < depidx[myjob+1]; k++)
                ATOMIC_ADD(&incomingdep[depdata[k]], -1);
          }
      }

    )RAW";


    template <typename T>
    class DeviceCholeskyKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> reorder, reorder_add, multdiag, reset, solveL, solveLT;
      unsigned groupsize;      // elementwise kernels
      unsigned lanes;          // simd width, x-size of the solve groups
      unsigned solvegroups;    // persistent groups

      DeviceCholeskyKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        if constexpr (is_same_v<T,double>)
          if (!device->HasFloat64())
            throw Exception("DeviceSparseCholesky<double> on "+device->Name()+
                            ", which has no fp64 - use DeviceSparseCholesky<float>");

        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         SubstituteScal (kernel_source, scal));
        reorder     = library->GetKernel ("chol_reorder");
        reorder_add = library->GetKernel ("chol_reorder_add");
        multdiag    = library->GetKernel ("chol_multdiag");
        reset       = library->GetKernel ("chol_reset");
        solveL      = library->GetKernel ("chol_solveL");
        solveLT     = library->GetKernel ("chol_solveLT");
        queue = device->DefaultQueue();

        bool gpu = device->SimdWidth() > 1;
        groupsize = gpu ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
        // the host backend runs groups one after another, so all jobs
        // must be taken by the rows of a single group
        lanes = gpu ? device->SimdWidth() : 8;
        solvegroups = gpu ? 512 : 1;
      }

      void LaunchElementwise (const shared_ptr<Kernel> & kernel, size_t n,
                              std::initializer_list<KernelArg> args) const
      {
        if (n == 0) return;
        unsigned groups = (n + groupsize-1) / groupsize;
        queue->Launch (*kernel, Dim3(groups), Dim3(groupsize), args);
      }

      static const DeviceCholeskyKernels & Get()
      {
        static mutex mtx;
        static shared_ptr<DeviceCholeskyKernels> cached;

        auto dev = GetGpuDevice();
        auto lock = lock_guard<mutex>(mtx);
        if (!cached || cached->device != dev)
          cached = make_shared<DeviceCholeskyKernels> (dev);
        return *cached;
      }
    };


    template <typename T>
    TypedBuffer<T> Upload (Device & device, FlatArray<T> data)
    {
      auto buf = device.NewBuffer<T> (max<size_t>(data.Size(),1), MemType::Device);
      buf.H2D (data.Data(), data.Size());
      return buf;
    }

    // dependency table as csr
    void TableToCSR (FlatTable<int> table, Array<int> & index, Array<int> & data)
    {
      index.SetSize (table.Size()+1);
      data.SetSize (table.AsArray().Size());
      index[0] = 0;
      for (size_t i = 0; i < table.Size(); i++)
        index[i+1] = index[i] + table[i].Size();
      data = table.AsArray();
    }
  }



  template <typename T>
  template <typename TM>
  DeviceSparseCholesky<T> :: DeviceSparseCholesky (const SparseCholeskyTM<TM> & mat)
  {
    height = mat.Height();
    width = mat.Width();
    nused = mat.GetNUsed();

    const auto & kern = DeviceCholeskyKernels<T>::Get();
    device = kern.device;
    queue = kern.queue;
    memtype = PreferredMemType();

    // micro-tasks and their dependency graphs
    auto tasks = mat.GetMicroTasks();
    njobs = tasks.Size();
    Array<int> hosttasks(4*njobs);
    for (size_t i = 0; i < njobs; i++)
      {
        hosttasks[4*i]   = tasks[i].blocknr;
        hosttasks[4*i+1] = int(tasks[i].type);
        hosttasks[4*i+2] = tasks[i].bblock;
        hosttasks[4*i+3] = tasks[i].nbblocks;
      }

    auto dep = mat.GetMicroDependency();
    auto dept = mat.GetMicroDependencyTranspose();
    Array<int> depidx, depdata, deptidx, deptdata;
    TableToCSR (dep, depidx, depdata);
    TableToCSR (dept, deptidx, deptdata);

    // the forward solve waits for every predecessor, the backward one
    // for every successor
    Array<int> incoming(njobs), incoming_trans(njobs);
    incoming = 0;
    for (size_t i = 0; i < njobs; i++)
      for (int d : dep[i])
        {
          if (d <= int(i)) throw Exception("DeviceSparseCholesky: dependency graph must be directional");
          incoming[d]++;
        }
    for (size_t i = 0; i < njobs; i++)
      incoming_trans[i] = dep[i].Size();

    // 64-bit offsets to int32
    auto firstinrow = mat.GetFirstInRow();
    auto firstinrow_ri = mat.GetFirstInRowRI();
    if (firstinrow[nused] > size_t(INT_MAX) || firstinrow_ri[nused] > size_t(INT_MAX))
      throw Exception("DeviceSparseCholesky: factor too large for 32-bit offsets");
    Array<int> hfirstinrow(nused+1), hfirstinrow_ri(nused+1);
    for (size_t i = 0; i <= nused; i++)
      {
        hfirstinrow[i] = int(firstinrow[i]);
        hfirstinrow_ri[i] = int(firstinrow_ri[i]);
      }

    auto lfact = mat.GetLFact();
    auto diag = mat.GetDiag();
    Array<T> hlfact(lfact.Size()), hdiag(diag.Size());
    for (size_t i = 0; i < lfact.Size(); i++) hlfact[i] = T(lfact[i]);
    for (size_t i = 0; i < diag.Size(); i++) hdiag[i] = T(diag[i]);

    cout << IM(7) << "DeviceSparseCholesky<" << (is_same_v<T,double> ? "double" : "float")
         << "> height = " << height << ", nused = " << nused << ", lfact = " << lfact.Size()
         << ", microtasks = " << njobs << endl;

    dev_microtasks = Upload<int> (*device, hosttasks);
    dev_depidx = Upload<int> (*device, depidx);
    dev_depdata = Upload<int> (*device, depdata);
    dev_depidx_trans = Upload<int> (*device, deptidx);
    dev_depdata_trans = Upload<int> (*device, deptdata);
    dev_incomingdep0 = Upload<int> (*device, incoming);
    dev_incomingdep0_trans = Upload<int> (*device, incoming_trans);
    dev_incomingdep = device->NewBuffer<int> (max<size_t>(njobs,1), MemType::Device);
    dev_incomingdep_trans = device->NewBuffer<int> (max<size_t>(njobs,1), MemType::Device);
    dev_cnt = device->NewBuffer<int> (2, MemType::Device);

    dev_blocks = Upload<int> (*device, mat.GetBlocks());
    dev_rowindex2 = Upload<int> (*device, mat.GetRowIndex2());
    dev_firstinrow_ri = Upload<int> (*device, hfirstinrow_ri);
    dev_firstinrow = Upload<int> (*device, hfirstinrow);
    dev_lfact = Upload<T> (*device, hlfact);
    dev_diag = Upload<T> (*device, hdiag);
    dev_order = Upload<int> (*device, mat.GetOrder());
    dev_hx = device->NewBuffer<T> (max<size_t>(nused,1), MemType::Device);
  }


  template <typename T>
  void DeviceSparseCholesky<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceSparseCholesky::MultAdd"); RegionTimer reg(t);
    if (x.Size() != width || y.Size() != height)
      throw Exception("DeviceSparseCholesky::MultAdd - size mismatch");
    if (nused == 0) return;

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    const auto & kern = DeviceCholeskyKernels<T>::Get();

    kern.LaunchElementwise (kern.reorder, height,
                            { ux.DevArgRO(), KernelArg(dev_hx), KernelArg(dev_order), KernelArg(int(height)) });
    kern.LaunchElementwise (kern.reset, max<size_t>(njobs,2),
                            { KernelArg(dev_incomingdep0), KernelArg(dev_incomingdep),
                              KernelArg(dev_incomingdep0_trans), KernelArg(dev_incomingdep_trans),
                              KernelArg(dev_cnt), KernelArg(int(njobs)) });

    // one row of lanes per job in flight; the backward solve has more
    // parallelism per task and runs 8 rows per group
    queue->Launch (*kern.solveL, Dim3(kern.solvegroups), Dim3(kern.lanes, 1),
                   { KernelArg(dev_depidx), KernelArg(dev_depdata), KernelArg(dev_incomingdep),
                     KernelArg(dev_hx), KernelArg(dev_hx), KernelArg(dev_cnt),
                     KernelArg(dev_microtasks), KernelArg(dev_blocks), KernelArg(dev_rowindex2),
                     KernelArg(dev_firstinrow_ri), KernelArg(dev_firstinrow), KernelArg(dev_lfact),
                     KernelArg(int(njobs)) });

    kern.LaunchElementwise (kern.multdiag, nused,
                            { KernelArg(dev_diag), KernelArg(dev_hx), KernelArg(int(nused)) });

    queue->Launch (*kern.solveLT, Dim3(kern.solvegroups), Dim3(kern.lanes, 8),
                   { KernelArg(dev_depidx_trans), KernelArg(dev_depdata_trans), KernelArg(dev_incomingdep_trans),
                     KernelArg(dev_hx), KernelArg(dev_hx), KernelArg(dev_cnt),
                     KernelArg(dev_microtasks), KernelArg(dev_blocks), KernelArg(dev_rowindex2),
                     KernelArg(dev_firstinrow_ri), KernelArg(dev_firstinrow), KernelArg(dev_lfact),
                     KernelArg(int(njobs)) });

    kern.LaunchElementwise (kern.reorder_add, height,
                            { KernelArg(dev_hx), uy.DevArgRW(), KernelArg(dev_order),
                              KernelArg(T(s)), KernelArg(int(height)) });
  }


  template <typename T>
  AutoVector DeviceSparseCholesky<T> :: CreateRowVector () const
  {
    return make_unique<DeviceVector<T>> (width, memtype);
  }

  template <typename T>
  AutoVector DeviceSparseCholesky<T> :: CreateColVector () const
  {
    return make_unique<DeviceVector<T>> (height, memtype);
  }

  template <typename T>
  BaseMatrix::OperatorInfo DeviceSparseCholesky<T> :: GetOperatorInfo () const
  {
    return { string("DeviceSparseCholesky<") + (is_same_v<T,double> ? "double" : "float") + ">",
             height, width };
  }

  template <typename T>
  ostream & DeviceSparseCholesky<T> :: Print (ostream & ost) const
  {
    ost << "DeviceSparseCholesky<" << (is_same_v<T,double> ? "double" : "float")
        << ">, height = " << height << ", nused = " << nused
        << ", microtasks = " << njobs << ", on " << device->Name() << endl;
    return ost;
  }


  template class DeviceSparseCholesky<double>;
  template class DeviceSparseCholesky<float>;
  template DeviceSparseCholesky<double>::DeviceSparseCholesky (const SparseCholeskyTM<double> &);
  template DeviceSparseCholesky<double>::DeviceSparseCholesky (const SparseCholeskyTM<float> &);
  template DeviceSparseCholesky<float>::DeviceSparseCholesky (const SparseCholeskyTM<double> &);
  template DeviceSparseCholesky<float>::DeviceSparseCholesky (const SparseCholeskyTM<float> &);
}
