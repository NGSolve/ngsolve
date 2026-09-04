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

    )RAW";


    // the spmv kernels, compiled once per device and scalar type
    template <typename T>
    class DeviceSparseKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> spmv_row, spmv_lanes;
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
        spmv_row   = library->GetKernel ("spmv_row");
        spmv_lanes = library->GetKernel ("spmv_lanes");
        queue = device->DefaultQueue();
        // the cpu reference backend runs one OS thread per work-item
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
      }

      // work-items per row: a power of two up to the simd width,
      // matching the average row length
      int ChooseLanes (size_t nze, size_t rows) const
      {
        size_t simd = device->SimdWidth();
        if (simd <= 1 || rows == 0) return 1;
        double avg = double(nze) / rows;
        int lanes = 1;
        while (2*lanes <= avg && size_t(2*lanes) <= simd && size_t(2*lanes) <= groupsize)
          lanes *= 2;
        return lanes;
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
    lanes = kern.ChooseLanes (nze, height);

    cout << IM(7) << "DeviceSparseMatrix<" << (is_same_v<T,double> ? "double" : "float")
         << "> height = " << height << ", width = " << width << ", nze = " << nze
         << ", lanes = " << lanes << endl;

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

    lanes_trans = DeviceSparseKernels<T>::Get().ChooseLanes (nze, width);
    devt_colnr = tc;
    devt_values = tv;
    devt_firsti = tf;    // last: its presence marks the transpose as ready
  }


  template <typename T>
  void DeviceSparseMatrix<T> :: LaunchSpMV (const TypedBuffer<int> & firsti,
                                            const TypedBuffer<int> & colnr,
                                            const TypedBuffer<T> & values,
                                            int alanes, size_t rows,
                                            KernelArg x, KernelArg y, T s) const
  {
    if (rows == 0) return;
    const auto & kern = DeviceSparseKernels<T>::Get();
    if (alanes == 1)
      {
        unsigned groups = (rows + kern.groupsize-1) / kern.groupsize;
        queue->Launch (*kern.spmv_row, Dim3(groups), Dim3(kern.groupsize),
                       { KernelArg(firsti), KernelArg(colnr), KernelArg(values),
                         x, y, KernelArg(s), KernelArg(int(rows)) });
      }
    else
      {
        size_t items = rows * alanes;
        unsigned groups = (items + kern.groupsize-1) / kern.groupsize;
        queue->Launch (*kern.spmv_lanes, Dim3(groups), Dim3(kern.groupsize),
                       { KernelArg(firsti), KernelArg(colnr), KernelArg(values),
                         x, y, KernelArg(s), KernelArg(int(alanes)), KernelArg(int(rows)) });
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
    LaunchSpMV (dev_firsti, dev_colnr, dev_values, lanes, height,
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
    LaunchSpMV (devt_firsti, devt_colnr, devt_values, lanes_trans, width,
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
