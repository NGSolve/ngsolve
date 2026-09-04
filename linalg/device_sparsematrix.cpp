/*********************************************************************/
/* File:   device_sparsematrix.cpp                                   */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   3. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_SPARSEMATRIX_CPP

#include <la.hpp>
#include <gpukernel.hpp>
#include <climits>
#include <chrono>

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

    const char * kernel_source = R"RAW(

      // y[row] += s * sum_j val[j]*x[colnr[j]], one work-item per row
      KERNEL(spmv_row, GLOBAL_IN(int,firsti), GLOBAL_IN(int,colnr), GLOBAL_IN(SCAL,val),
                       GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y), VALUE(SCAL,s), VALUE(int,h))
      {
        int row = int(GLOBAL_ID_X);
        if (row >= h) return;
        SCAL sum = 0;
        int last = firsti[row+1];
        for (int j = firsti[row]; j < last; j++)
          sum += val[j]*x[colnr[j]];
        y[row] += s*sum;
      }

      // lanes consecutive work-items share a row (lanes a power of two
      // dividing the group size), partial sums reduced in group memory
      KERNEL(spmv_lanes, GLOBAL_IN(int,firsti), GLOBAL_IN(int,colnr), GLOBAL_IN(SCAL,val),
                         GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y), VALUE(SCAL,s),
                         VALUE(int,lanes), VALUE(int,h))
      {
        SHARED(SCAL, tmp, 1024);
        int lid = int(LOCAL_ID_X);
        int lane = lid & (lanes-1);
        int row = int(GLOBAL_ID_X) / lanes;
        SCAL sum = 0;
        if (row < h)
          {
            int last = firsti[row+1];
            for (int j = firsti[row]+lane; j < last; j += lanes)
              sum += val[j]*x[colnr[j]];
          }
        tmp[lid] = sum;
        BARRIER();
        for (int d = lanes/2; d > 0; d /= 2)
          {
            if (lane < d) tmp[lid] += tmp[lid+d];
            BARRIER();
          }
        if (lane == 0 && row < h) y[row] += s*tmp[lid];
      }

      // two rows per lane group, two steps unrolled: four gathers in
      // flight per lane, which is what hides the gather latency
      KERNEL(spmv_rows2u2, GLOBAL_IN(int,firsti), GLOBAL_IN(int,colnr), GLOBAL_IN(SCAL,val),
                                 GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y), VALUE(SCAL,s),
                                 VALUE(int,lanes), VALUE(int,h))
      {
        SHARED(SCAL, tmp, 1024);
        SHARED(SCAL, tmp2, 1024);
        int lid = int(LOCAL_ID_X);
        int lane = lid & (lanes-1);
        int row = 2*(int(GLOBAL_ID_X) / lanes);
        SCAL sum0 = 0, sum1 = 0;
        if (row < h)
          {
            int j0 = firsti[row] + lane, l0 = firsti[row+1];
            int j1 = (row+1 < h) ? firsti[row+1] + lane : 0, l1 = (row+1 < h) ? firsti[row+2] : 0;
            while (j0 < l0 || j1 < l1)
              {
                int j0b = j0 + lanes, j1b = j1 + lanes;
                int c0 = (j0 < l0) ? colnr[j0] : 0, c0b = (j0b < l0) ? colnr[j0b] : 0;
                int c1 = (j1 < l1) ? colnr[j1] : 0, c1b = (j1b < l1) ? colnr[j1b] : 0;
                SCAL v0 = (j0 < l0) ? val[j0] : SCAL(0), v0b = (j0b < l0) ? val[j0b] : SCAL(0);
                SCAL v1 = (j1 < l1) ? val[j1] : SCAL(0), v1b = (j1b < l1) ? val[j1b] : SCAL(0);
                SCAL x0 = x[c0], x0b = x[c0b], x1 = x[c1], x1b = x[c1b];
                sum0 += v0*x0 + v0b*x0b;
                sum1 += v1*x1 + v1b*x1b;
                j0 += 2*lanes; j1 += 2*lanes;
              }
          }
        tmp[lid] = sum0; tmp2[lid] = sum1;
        BARRIER();
        for (int d = lanes/2; d > 0; d /= 2)
          {
            if (lane < d) { tmp[lid] += tmp[lid+d]; tmp2[lid] += tmp2[lid+d]; }
            BARRIER();
          }
        if (lane == 0 && row < h) y[row] += s*tmp[lid];
        if (lane == 0 && row+1 < h) y[row+1] += s*tmp2[lid];
      }

    )RAW";


    // the spmv kernels, compiled once per device and scalar type
    template <typename T>
    class DeviceSparseKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> spmv_row, spmv_lanes, spmv_rows2u2;
      unsigned groupsize;
      string forced_variant;     // NGS_GPU_SPMV=lanes|rows2u2 skips the timing
      int forced_lanes = 0;      // NGS_GPU_SPMV_LANES

      DeviceSparseKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        if constexpr (is_same_v<T,double>)
          if (!device->HasFloat64())
            throw Exception("DeviceSparseMatrix<double> on "+device->Name()+
                            ", which has no fp64 - use DeviceSparseMatrix<float>");

        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         SubstituteScal (kernel_source, scal));
        spmv_row     = library->GetKernel ("spmv_row");
        spmv_lanes   = library->GetKernel ("spmv_lanes");
        spmv_rows2u2 = library->GetKernel ("spmv_rows2u2");
        queue = device->DefaultQueue();
        // the cpu reference backend runs one OS thread per work-item
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        if (getenv("NGS_GPU_SPMV_GROUP")) groupsize = atoi (getenv("NGS_GPU_SPMV_GROUP"));
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
        forced_variant = getenv("NGS_GPU_SPMV") ? getenv("NGS_GPU_SPMV") : "";
        forced_lanes = getenv("NGS_GPU_SPMV_LANES") ? atoi (getenv("NGS_GPU_SPMV_LANES")) : 0;
      }

      static const DeviceSparseKernels & Get()
      {
        static mutex mtx;
        static shared_ptr<DeviceSparseKernels> cached;

        auto dev = GetGpuDevice();
        auto lock = lock_guard<mutex>(mtx);
        if (!cached || cached->device != dev)
          cached = make_shared<DeviceSparseKernels> (dev);
        return *cached;
      }
    };
  }



  template <typename T>
  template <typename TM>
  DeviceSparseMatrix<T> :: DeviceSparseMatrix (const SparseMatrixTM<TM> & mat)
  {
    height = mat.Height();
    width = mat.Width();
    nze = mat.NZE();
    if (nze > size_t(INT_MAX) || height > size_t(INT_MAX) || width > size_t(INT_MAX))
      throw Exception("DeviceSparseMatrix: matrix too large for 32-bit indices");

    const auto & kern = DeviceSparseKernels<T>::Get();
    device = kern.device;
    queue = kern.queue;
    memtype = PreferredMemType();

    // row starts to int32, values to T
    Array<int> firsti (height+1);
    auto hfirsti = mat.GetFirstArray();
    for (size_t i = 0; i <= height; i++) firsti[i] = int(hfirsti[i]);

    Array<T> values (nze);
    auto hvalues = mat.GetValues();
    for (size_t j = 0; j < nze; j++) values[j] = T(hvalues(j));

    // the index buffers are never written by a kernel, keep them off the host
    dev_firsti = device->NewBuffer<int> (height+1, MemType::Device);
    dev_colnr  = device->NewBuffer<int> (max<size_t>(nze,1), MemType::Device);
    dev_values = device->NewBuffer<T> (max<size_t>(nze,1), MemType::Device);

    dev_firsti.H2D (firsti.Data(), height+1);
    dev_colnr.H2D (mat.GetColIndices().Data(), nze);
    dev_values.H2D (values.Data(), nze);

    choice = AutoTune (dev_firsti, dev_colnr, dev_values, height, width);
    cout << IM(7) << "DeviceSparseMatrix<" << (is_same_v<T,double> ? "double" : "float")
         << "> height = " << height << ", width = " << width << ", nze = " << nze
         << ", kernel " << choice.kernel->Name() << ", lanes = " << choice.lanes << endl;
  }


  /*
    The best kernel and lane count differ between gpus: a discrete card
    hides the gather latency by occupancy and wants more lanes per row,
    unified-memory gpus want few lanes and several rows in flight per
    lane. Timing the candidates on the matrix itself costs a few solves
    once and decides it for this matrix.
  */
  template <typename T>
  typename DeviceSparseMatrix<T>::SpMVChoice
  DeviceSparseMatrix<T> :: AutoTune (const TypedBuffer<int> & firsti, const TypedBuffer<int> & colnr,
                                     const TypedBuffer<T> & values, size_t rows, size_t cols) const
  {
    const auto & kern = DeviceSparseKernels<T>::Get();
    SpMVChoice best;
    best.kernel = kern.spmv_row;
    if (device->SimdWidth() <= 1 || rows == 0 || nze == 0) return best;

    std::vector<SpMVChoice> candidates;
    std::vector<int> lanes_list;
    for (int l = 4; l <= 32 && size_t(l) <= device->SimdWidth() && size_t(l) <= kern.groupsize; l *= 2)
      lanes_list.push_back (l);
    if (kern.forced_lanes) lanes_list = { kern.forced_lanes };
    for (int l : lanes_list)
      {
        if (kern.forced_variant != "rows2u2") candidates.push_back ({ kern.spmv_lanes, l, 1 });
        if (kern.forced_variant != "lanes")   candidates.push_back ({ kern.spmv_rows2u2, l, 2 });
      }
    if (candidates.size() == 1) return candidates[0];

    auto x = device->NewBuffer<T> (cols, MemType::Device);
    auto y = device->NewBuffer<T> (rows, MemType::Device);
    {
      Array<T> zeros(max(rows, cols)); zeros = T(0);
      x.H2D (zeros.Data(), cols); y.H2D (zeros.Data(), rows);
    }
    double best_time = 1e300;
    for (auto & c : candidates)
      {
        LaunchSpMV (firsti, colnr, values, c, rows, KernelArg(x), KernelArg(y), T(1));
        queue->Finish();
        auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < 3; i++)
          LaunchSpMV (firsti, colnr, values, c, rows, KernelArg(x), KernelArg(y), T(1));
        queue->Finish();
        double t = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
        if (t < best_time) { best_time = t; best = c; }
      }
    return best;
  }


  template <typename T>
  void DeviceSparseMatrix<T> :: BuildTranspose() const
  {
    lock_guard<mutex> lock(trans_mutex);
    if (devt_firsti) return;

    static Timer t("DeviceSparseMatrix::BuildTranspose"); RegionTimer reg(t);

    Array<int> firsti(height+1), colnr(nze);
    Array<T> values(nze);
    queue->Finish();
    dev_firsti.D2H (firsti.Data(), height+1);
    dev_colnr.D2H (colnr.Data(), nze);
    dev_values.D2H (values.Data(), nze);

    // csr of the transpose: count per column, prefix sum, scatter
    Array<int> tfirsti(width+1);
    tfirsti = 0;
    for (size_t j = 0; j < nze; j++) tfirsti[colnr[j]+1]++;
    for (size_t c = 0; c < width; c++) tfirsti[c+1] += tfirsti[c];

    Array<int> pos(width), tcolnr(nze);
    Array<T> tvalues(nze);
    for (size_t c = 0; c < width; c++) pos[c] = tfirsti[c];
    for (size_t i = 0; i < height; i++)
      for (int j = firsti[i]; j < firsti[i+1]; j++)
        {
          int k = pos[colnr[j]]++;
          tcolnr[k] = int(i);
          tvalues[k] = values[j];
        }

    auto tf = device->NewBuffer<int> (width+1, MemType::Device);
    auto tc = device->NewBuffer<int> (max<size_t>(nze,1), MemType::Device);
    auto tv = device->NewBuffer<T> (max<size_t>(nze,1), MemType::Device);
    tf.H2D (tfirsti.Data(), width+1);
    tc.H2D (tcolnr.Data(), nze);
    tv.H2D (tvalues.Data(), nze);

    choice_trans = AutoTune (tf, tc, tv, width, height);
    devt_colnr = tc;
    devt_values = tv;
    devt_firsti = tf;    // last: its presence marks the transpose as ready
  }


  template <typename T>
  void DeviceSparseMatrix<T> :: LaunchSpMV (const TypedBuffer<int> & firsti,
                                            const TypedBuffer<int> & colnr,
                                            const TypedBuffer<T> & values,
                                            const SpMVChoice & ch, size_t rows,
                                            KernelArg x, KernelArg y, T s) const
  {
    if (rows == 0) return;
    const auto & kern = DeviceSparseKernels<T>::Get();
    if (ch.lanes == 1)
      {
        unsigned groups = (rows + kern.groupsize-1) / kern.groupsize;
        queue->Launch (*ch.kernel, Dim3(groups), Dim3(kern.groupsize),
                       { KernelArg(firsti), KernelArg(colnr), KernelArg(values),
                         x, y, KernelArg(s), KernelArg(int(rows)) });
      }
    else
      {
        size_t items = ((rows + ch.rows_per_group-1) / ch.rows_per_group) * ch.lanes;
        unsigned groups = (items + kern.groupsize-1) / kern.groupsize;
        queue->Launch (*ch.kernel, Dim3(groups), Dim3(kern.groupsize),
                       { KernelArg(firsti), KernelArg(colnr), KernelArg(values),
                         x, y, KernelArg(s), KernelArg(int(ch.lanes)), KernelArg(int(rows)) });
      }
  }


  template <typename T>
  void DeviceSparseMatrix<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceSparseMatrix::MultAdd"); RegionTimer reg(t);
    if (x.Size() != width || y.Size() != height)
      throw Exception("DeviceSparseMatrix::MultAdd - size mismatch");

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    LaunchSpMV (dev_firsti, dev_colnr, dev_values, choice, height,
                ux.DevArgRO(), uy.DevArgRW(), T(s));
  }


  template <typename T>
  void DeviceSparseMatrix<T> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceSparseMatrix::MultTransAdd"); RegionTimer reg(t);
    if (x.Size() != height || y.Size() != width)
      throw Exception("DeviceSparseMatrix::MultTransAdd - size mismatch");

    if (!devt_firsti) BuildTranspose();

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    LaunchSpMV (devt_firsti, devt_colnr, devt_values, choice_trans, width,
                ux.DevArgRO(), uy.DevArgRW(), T(s));
  }


  template <typename T>
  AutoVector DeviceSparseMatrix<T> :: CreateRowVector () const
  {
    return make_unique<DeviceVector<T>> (width, memtype);
  }

  template <typename T>
  AutoVector DeviceSparseMatrix<T> :: CreateColVector () const
  {
    return make_unique<DeviceVector<T>> (height, memtype);
  }

  template <typename T>
  BaseMatrix::OperatorInfo DeviceSparseMatrix<T> :: GetOperatorInfo () const
  {
    return { string("DeviceSparseMatrix<") + (is_same_v<T,double> ? "double" : "float")
             + "> (nze=" + ToString(nze) + ")", height, width };
  }

  template <typename T>
  ostream & DeviceSparseMatrix<T> :: Print (ostream & ost) const
  {
    ost << "DeviceSparseMatrix<" << (is_same_v<T,double> ? "double" : "float")
        << ">, height = " << height << ", width = " << width
        << ", nze = " << nze << ", on " << device->Name() << endl;
    return ost;
  }


  template class DeviceSparseMatrix<double>;
  template class DeviceSparseMatrix<float>;
  template DeviceSparseMatrix<double>::DeviceSparseMatrix (const SparseMatrixTM<double> &);
  template DeviceSparseMatrix<double>::DeviceSparseMatrix (const SparseMatrixTM<float> &);
  template DeviceSparseMatrix<float>::DeviceSparseMatrix (const SparseMatrixTM<double> &);
  template DeviceSparseMatrix<float>::DeviceSparseMatrix (const SparseMatrixTM<float> &);
}
