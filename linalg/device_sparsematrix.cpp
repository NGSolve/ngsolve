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

      // y[row] = beta*y[row] + s * sum_j val[j]*x[colnr[j]]; beta 0 overwrites
      #define STORE_ROW(r, v) { if (beta == SCAL(0)) y[r] = s*(v); else y[r] = beta*y[r] + s*(v); }

      KERNEL(spmv_row, GLOBAL_IN(int,firsti), GLOBAL_IN(int,colnr), GLOBAL_IN(SCAL,val),
                       GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y), VALUE(SCAL,s), VALUE(SCAL,beta), VALUE(int,h))
      {
        int row = int(GLOBAL_ID_X);
        if (row >= h) return;
        SCAL sum = 0;
        int last = firsti[row+1];
        for (int j = firsti[row]; j < last; j++)
          sum += val[j]*x[colnr[j]];
        STORE_ROW(row, sum);
      }

      KERNEL(scale_vec, GLOBAL(SCAL,y), VALUE(SCAL,beta), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i >= n) return;
        y[i] = (beta == SCAL(0)) ? SCAL(0) : beta*y[i];
      }

      KERNEL(spmv_sym, GLOBAL_IN(int,firsti), GLOBAL_IN(int,colnr), GLOBAL_IN(SCAL,val),
                       GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y), VALUE(SCAL,s), VALUE(int,h))
      {
        int row = int(GLOBAL_ID_X);
        if (row >= h) return;
        SCAL xr = x[row];
        SCAL sum = 0;
        int last = firsti[row+1];
        for (int j = firsti[row]; j < last; j++)
          {
            int c = colnr[j];
            SCAL v = val[j];
            sum += v*x[c];
            if (c != row) ATOMIC_ADD (&y[c], s*v*xr);
          }
        ATOMIC_ADD (&y[row], s*sum);
      }

      // transposed product by scattering: y[colnr[j]] += s*val[j]*x[row],
      // lanes consecutive work-items share a row
      KERNEL(spmvT_lanes, GLOBAL_IN(int,firsti), GLOBAL_IN(int,colnr), GLOBAL_IN(SCAL,val),
                          GLOBAL_IN(SCAL,x), GLOBAL_ATOMIC(SCAL,y), VALUE(SCAL,s),
                          VALUE(int,lanes), VALUE(int,h))
      {
        int lane = int(LOCAL_ID_X) & (lanes-1);
        int row = int(GLOBAL_ID_X) / lanes;
        if (row >= h) return;
        SCAL xr = s*x[row];
        int last = firsti[row+1];
        for (int j = firsti[row]+lane; j < last; j += lanes)
          ATOMIC_ADD (&y[colnr[j]], val[j]*xr);
      }

      // lanes consecutive work-items share a row (lanes a power of two
      // dividing the group size), partial sums reduced in group memory
      KERNEL(spmv_lanes, GLOBAL_IN(int,firsti), GLOBAL_IN(int,colnr), GLOBAL_IN(SCAL,val),
                         GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y), VALUE(SCAL,s), VALUE(SCAL,beta),
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
        if (lane == 0 && row < h) STORE_ROW(row, tmp[lid]);
      }

      // two rows per lane group, two steps unrolled: four gathers in
      // flight per lane, which is what hides the gather latency
      KERNEL(spmv_rows2u2, GLOBAL_IN(int,firsti), GLOBAL_IN(int,colnr), GLOBAL_IN(SCAL,val),
                                 GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y), VALUE(SCAL,s), VALUE(SCAL,beta),
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
        if (lane == 0 && row < h) STORE_ROW(row, tmp[lid]);
        if (lane == 0 && row+1 < h) STORE_ROW(row+1, tmp2[lid]);
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
      shared_ptr<Kernel> spmv_row, spmv_lanes, spmv_rows2u2, spmv_sym, spmvT_lanes, scale_vec;
      unsigned groupsize;

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
        spmv_sym     = library->GetKernel ("spmv_sym");
        spmvT_lanes  = library->GetKernel ("spmvT_lanes");
        scale_vec    = library->GetKernel ("scale_vec");
        queue = device->DefaultQueue();
        // the cpu reference backend runs one OS thread per work-item
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
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
  DeviceSparseMatrix<T> :: DeviceSparseMatrix (const SparseMatrixTM<TM> & mat, bool asymmetric)
  {
    symmetric = asymmetric;
    height = mat.Height();
    width = mat.Width();
    nze = mat.NZE();
    if (nze > size_t(INT_MAX) || height > size_t(INT_MAX) || width > size_t(INT_MAX))
      throw Exception("DeviceSparseMatrix: matrix too large for 32-bit indices");

    const auto & kern = DeviceSparseKernels<T>::Get();
    device = kern.device;
    queue = kern.queue;
    memtype = PreferredMemType();

    static Timer tup("DeviceSparseMatrix ctor upload");
    auto hfirsti = mat.GetFirstArray();
    auto hcolnr = mat.GetColIndices();
    auto hvalues = mat.GetValues();

    // converted straight into the buffer (unified memory) or into the
    // backend's pinned staging chunks, no host temporary
    tup.Start();
    dev_firsti = device->NewBuffer<int> (height+1, memtype);
    dev_colnr  = device->NewBuffer<int> (max<size_t>(nze,1), memtype);
    dev_values = device->NewBuffer<T> (max<size_t>(nze,1), memtype);

    dev_firsti.Fill (height+1, [&] (int * dst, size_t off, size_t n)
      { ParallelFor (n, [&] (size_t i) { dst[i] = int(hfirsti[off+i]); }); });
    dev_colnr.H2D (hcolnr.Data(), nze);
    dev_values.Fill (nze, [&] (T * dst, size_t off, size_t n)
      { ParallelFor (n, [&] (size_t j) { dst[j] = T(hvalues(off+j)); }); });
    tup.Stop();

    if (symmetric)
      {
        choice.kernel = kern.spmv_sym;
        cout << IM(7) << "DeviceSparseMatrix<" << (is_same_v<T,double> ? "double" : "float")
             << "> symmetric storage, height = " << height << ", nze = " << nze << endl;
        return;
      }
    choice = ChooseKernel (height);
    lanes_trans = ChooseLanesTrans (height);
    cout << IM(7) << "DeviceSparseMatrix<" << (is_same_v<T,double> ? "double" : "float")
         << "> height = " << height << ", width = " << width << ", nze = " << nze
         << ", kernel " << choice.kernel->Name() << ", lanes = " << choice.lanes
         << ", transposed lanes = " << lanes_trans << endl;
  }


  template <typename T>
  void DeviceSparseMatrix<T> :: LaunchSym (KernelArg x, KernelArg y, T s, T beta) const
  {
    const auto & kern = DeviceSparseKernels<T>::Get();
    unsigned groups = (height + kern.groupsize-1) / kern.groupsize;
    if (beta != T(1))
      queue->Launch (*kern.scale_vec, Dim3(groups), Dim3(kern.groupsize),
                     { y, KernelArg(beta), KernelArg(int(height)) });
    queue->Launch (*kern.spmv_sym, Dim3(groups), Dim3(kern.groupsize),
                   { KernelArg(dev_firsti), KernelArg(dev_colnr), KernelArg(dev_values),
                     x, y, KernelArg(s), KernelArg(int(height)) });
  }


  /*
    Fixed choice, from sweeps on an M4 Pro and an RTX 5090 (2026-09-10):
    a discrete card hides the gather latency by occupancy and wants
    spmv_lanes with lanes growing with the row length (4 / 8 / 16 for
    average rows below 16 / 48 / above); unified-memory gpus want few
    lanes and two rows in flight per lane (spmv_rows2u2, 4 lanes, 8 for
    rows of 48+). Timing the candidates per matrix was tried and dropped:
    it cost a few solves per matrix and mis-picked on a cold or busy gpu.
  */
  template <typename T>
  typename DeviceSparseMatrix<T>::SpMVChoice
  DeviceSparseMatrix<T> :: ChooseKernel (size_t rows) const
  {
    const auto & kern = DeviceSparseKernels<T>::Get();
    SpMVChoice ch;
    ch.kernel = kern.spmv_row;
    if (device->SimdWidth() <= 1 || rows == 0 || nze == 0) return ch;

    double avg = double(nze) / rows;
    bool unified = device->IsUnifiedMemory();
    int lanes;
    if (unified)
      lanes = avg < 48 ? 4 : 8;
    else
      lanes = avg < 16 ? 4 : avg < 48 ? 8 : 16;
    lanes = int(min<size_t> (lanes, min (device->SimdWidth(), size_t(kern.groupsize))));

    ch.kernel = unified ? kern.spmv_rows2u2 : kern.spmv_lanes;
    ch.lanes = lanes;
    ch.rows_per_group = unified ? 2 : 1;
    return ch;
  }


  /*
    Scatter lanes from the same sweep (2026-09-10): the atomic traffic
    favours more lanes per row on a discrete card (4 / 16 / 32 for
    average rows below 16 / 48 / above), few on unified memory (4, 8 for
    rows of 48+). The scatter costs 1.4-2x a forward product; a caller
    who needs many transposed products uploads the host transpose.
  */
  template <typename T>
  int DeviceSparseMatrix<T> :: ChooseLanesTrans (size_t rows) const
  {
    const auto & kern = DeviceSparseKernels<T>::Get();
    if (device->SimdWidth() <= 1 || rows == 0 || nze == 0) return 1;
    double avg = double(nze) / rows;
    int lanes;
    if (device->IsUnifiedMemory())
      lanes = avg < 48 ? 4 : 8;
    else
      lanes = avg < 16 ? 4 : avg < 48 ? 16 : 32;
    return int(min<size_t> (lanes, min (device->SimdWidth(), size_t(kern.groupsize))));
  }


  // y = beta*y + s*A^T x
  template <typename T>
  void DeviceSparseMatrix<T> :: LaunchSpMVT (KernelArg x, KernelArg y, T s, T beta) const
  {
    const auto & kern = DeviceSparseKernels<T>::Get();
    if (beta != T(1))
      {
        unsigned groups = (width + kern.groupsize-1) / kern.groupsize;
        queue->Launch (*kern.scale_vec, Dim3(groups), Dim3(kern.groupsize),
                       { y, KernelArg(beta), KernelArg(int(width)) });
      }
    if (height == 0) return;
    size_t items = height * lanes_trans;
    unsigned groups = (items + kern.groupsize-1) / kern.groupsize;
    queue->Launch (*kern.spmvT_lanes, Dim3(groups), Dim3(kern.groupsize),
                   { KernelArg(dev_firsti), KernelArg(dev_colnr), KernelArg(dev_values),
                     x, y, KernelArg(s), KernelArg(int(lanes_trans)), KernelArg(int(height)) });
  }


  template <typename T>
  void DeviceSparseMatrix<T> :: LaunchSpMV (const TypedBuffer<int> & firsti,
                                            const TypedBuffer<int> & colnr,
                                            const TypedBuffer<T> & values,
                                            const SpMVChoice & ch, size_t rows,
                                            KernelArg x, KernelArg y, T s, T beta) const
  {
    if (rows == 0) return;
    const auto & kern = DeviceSparseKernels<T>::Get();
    if (ch.lanes == 1)
      {
        unsigned groups = (rows + kern.groupsize-1) / kern.groupsize;
        queue->Launch (*ch.kernel, Dim3(groups), Dim3(kern.groupsize),
                       { KernelArg(firsti), KernelArg(colnr), KernelArg(values),
                         x, y, KernelArg(s), KernelArg(beta), KernelArg(int(rows)) });
      }
    else
      {
        size_t items = ((rows + ch.rows_per_group-1) / ch.rows_per_group) * ch.lanes;
        unsigned groups = (items + kern.groupsize-1) / kern.groupsize;
        queue->Launch (*ch.kernel, Dim3(groups), Dim3(kern.groupsize),
                       { KernelArg(firsti), KernelArg(colnr), KernelArg(values),
                         x, y, KernelArg(s), KernelArg(beta), KernelArg(int(ch.lanes)), KernelArg(int(rows)) });
      }
  }

  template <typename T>
  void DeviceSparseMatrix<T> :: Mult (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceSparseMatrix::Mult"); RegionTimer reg(t);
    if (x.Size() != width || y.Size() != height)
      throw Exception("DeviceSparseMatrix::Mult - size mismatch");
    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    if (symmetric) { LaunchSym (ux.DevArgRO(), uy.DevArgW(), T(1), T(0)); return; }
    LaunchSpMV (dev_firsti, dev_colnr, dev_values, choice, height,
                ux.DevArgRO(), uy.DevArgW(), T(1), T(0));
  }

  template <typename T>
  void DeviceSparseMatrix<T> :: MultTrans (const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceSparseMatrix::MultTrans"); RegionTimer reg(t);
    if (x.Size() != height || y.Size() != width)
      throw Exception("DeviceSparseMatrix::MultTrans - size mismatch");
    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    if (symmetric) { LaunchSym (ux.DevArgRO(), uy.DevArgW(), T(1), T(0)); return; }
    LaunchSpMVT (ux.DevArgRO(), uy.DevArgW(), T(1), T(0));
  }


  template <typename T>
  void DeviceSparseMatrix<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceSparseMatrix::MultAdd"); RegionTimer reg(t);
    if (x.Size() != width || y.Size() != height)
      throw Exception("DeviceSparseMatrix::MultAdd - size mismatch");

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    if (symmetric) { LaunchSym (ux.DevArgRO(), uy.DevArgRW(), T(s), T(1)); return; }
    LaunchSpMV (dev_firsti, dev_colnr, dev_values, choice, height,
                ux.DevArgRO(), uy.DevArgRW(), T(s), T(1));
  }


  template <typename T>
  void DeviceSparseMatrix<T> :: MultTransAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceSparseMatrix::MultTransAdd"); RegionTimer reg(t);
    if (x.Size() != height || y.Size() != width)
      throw Exception("DeviceSparseMatrix::MultTransAdd - size mismatch");

    DeviceVectorWrapper<T> ux(x, memtype);
    DeviceVectorWrapper<T> uy(y, memtype);
    if (symmetric) { LaunchSym (ux.DevArgRO(), uy.DevArgRW(), T(s), T(1)); return; }
    LaunchSpMVT (ux.DevArgRO(), uy.DevArgRW(), T(s), T(1));
  }




  template <typename T>
  BaseMatrix::OperatorInfo DeviceSparseMatrix<T> :: GetOperatorInfo () const
  {
    return { string("DeviceSparseMatrix<") + (is_same_v<T,double> ? "double" : "float")
             + (symmetric ? "> symmetric (nze=" : "> (nze=") + ToString(nze) + ")", height, width };
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
  template DeviceSparseMatrix<double>::DeviceSparseMatrix (const SparseMatrixTM<double> &, bool);
  template DeviceSparseMatrix<double>::DeviceSparseMatrix (const SparseMatrixTM<float> &, bool);
  template DeviceSparseMatrix<float>::DeviceSparseMatrix (const SparseMatrixTM<double> &, bool);
  template DeviceSparseMatrix<float>::DeviceSparseMatrix (const SparseMatrixTM<float> &, bool);
}
