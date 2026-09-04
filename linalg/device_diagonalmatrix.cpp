/*********************************************************************/
/* File:   device_diagonalmatrix.cpp                                 */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   4. Sep. 2026                                              */
/*********************************************************************/

#define FILE_DEVICE_DIAGONALMATRIX_CPP

#include <la.hpp>
#include <gpukernel.hpp>

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

      // y[i] += s * d[i] * x[i]
      KERNEL(diag_mult, GLOBAL_IN(SCAL,d), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                        VALUE(SCAL,s), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n) y[i] += s * d[i] * x[i];
      }

    )RAW";


    template <typename T>
    class DeviceDiagonalKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> mult;
      unsigned groupsize;

      DeviceDiagonalKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         SubstituteScal (kernel_source, scal));
        mult = library->GetKernel ("diag_mult");
        queue = device->DefaultQueue();
        groupsize = (device->SimdWidth() > 1) ? 256 : 64;
        groupsize = min<size_t> (groupsize, device->MaxThreadsPerGroup());
      }

      static const DeviceDiagonalKernels & Get()
      {
        static mutex mtx;
        static shared_ptr<DeviceDiagonalKernels> cached;

        auto dev = GetGpuDevice();
        auto lock = lock_guard<mutex>(mtx);
        if (!cached || cached->device != dev)
          cached = make_shared<DeviceDiagonalKernels> (dev);
        return *cached;
      }
    };
  }



  template <typename T>
  template <typename TS>
  DeviceDiagonalMatrix<T> :: DeviceDiagonalMatrix (FlatVector<TS> adiag)
    : diag (adiag.Size(), PreferredMemType())
  {
    // fills the host side, uploaded on first device use
    diag.FVScal() = adiag;
  }


  template <typename T>
  void DeviceDiagonalMatrix<T> :: MultAdd (double s, const BaseVector & x, BaseVector & y) const
  {
    static Timer t("DeviceDiagonalMatrix::MultAdd"); RegionTimer reg(t);
    size_t n = diag.Size();
    if (x.Size() != n || y.Size() != n)
      throw Exception("DeviceDiagonalMatrix::MultAdd - size mismatch");
    if (n == 0) return;

    DeviceVectorWrapper<T> ux(x, diag.GetMemType());
    DeviceVectorWrapper<T> uy(y, diag.GetMemType());

    const auto & kern = DeviceDiagonalKernels<T>::Get();
    unsigned groups = (n + kern.groupsize-1) / kern.groupsize;
    kern.queue->Launch (*kern.mult, Dim3(groups), Dim3(kern.groupsize),
                        { diag.DevArgRO(), ux.DevArgRO(), uy.DevArgRW(),
                          KernelArg(T(s)), KernelArg(int(n)) });
  }


  template <typename T>
  AutoVector DeviceDiagonalMatrix<T> :: CreateRowVector () const
  {
    return make_unique<DeviceVector<T>> (diag.Size(), diag.GetMemType());
  }

  template <typename T>
  AutoVector DeviceDiagonalMatrix<T> :: CreateColVector () const
  {
    return make_unique<DeviceVector<T>> (diag.Size(), diag.GetMemType());
  }

  template <typename T>
  BaseMatrix::OperatorInfo DeviceDiagonalMatrix<T> :: GetOperatorInfo () const
  {
    return { string("DeviceDiagonalMatrix<") + (is_same_v<T,double> ? "double" : "float") + ">",
             diag.Size(), diag.Size() };
  }

  template <typename T>
  ostream & DeviceDiagonalMatrix<T> :: Print (ostream & ost) const
  {
    ost << "DeviceDiagonalMatrix<" << (is_same_v<T,double> ? "double" : "float")
        << ">, size = " << diag.Size() << endl;
    return ost;
  }


  template class DeviceDiagonalMatrix<double>;
  template class DeviceDiagonalMatrix<float>;
  template DeviceDiagonalMatrix<double>::DeviceDiagonalMatrix (FlatVector<double>);
  template DeviceDiagonalMatrix<double>::DeviceDiagonalMatrix (FlatVector<float>);
  template DeviceDiagonalMatrix<float>::DeviceDiagonalMatrix (FlatVector<double>);
  template DeviceDiagonalMatrix<float>::DeviceDiagonalMatrix (FlatVector<float>);
}
