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


      // the projector: a 0/1 diagonal kept as a bit mask, packed into 32 bit words.
      #define NGS_PROJ_SET(i)  ((bits[(i) >> 5] >> ((i) & 31)) & 1)

      KERNEL(proj_mult, GLOBAL_IN(int,bits), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                        VALUE(int,keep), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n) y[i] = (NGS_PROJ_SET(i) == keep) ? x[i] : SCAL(0);
      }

      KERNEL(proj_multadd, GLOBAL_IN(int,bits), GLOBAL_IN(SCAL,x), GLOBAL(SCAL,y),
                           VALUE(SCAL,s), VALUE(int,keep), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n && NGS_PROJ_SET(i) == keep) y[i] += s * x[i];
      }

      // in place: drop everything the projector does not keep
      KERNEL(proj_project, GLOBAL_IN(int,bits), GLOBAL(SCAL,x),
                           VALUE(int,keep), VALUE(int,n))
      {
        int i = int(GLOBAL_ID_X);
        if (i < n && NGS_PROJ_SET(i) != keep) x[i] = SCAL(0);
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
      shared_ptr<Kernel> proj_mult, proj_multadd, proj_project;
      unsigned groupsize;

      DeviceDiagonalKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         SubstituteScal (kernel_source, scal));
        mult = library->GetKernel ("diag_mult");
        proj_mult    = library->GetKernel ("proj_mult");
        proj_multadd = library->GetKernel ("proj_multadd");
        proj_project = library->GetKernel ("proj_project");
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





  class DeviceProjector : public Projector
  {
    MemType memtype;
    TypedBuffer<int> dev_bits;

    void Launch (Kernel & kernel, size_t n,
                 const std::vector<KernelArg> & args) const
    {
      const auto & kern = DeviceDiagonalKernels<double>::Get();
      unsigned groups = (n + kern.groupsize-1) / kern.groupsize;
      kern.queue->Launch (kernel, Dim3(groups), Dim3(kern.groupsize), args);
    }

  public:
    DeviceProjector (const Projector & proj)
      : Projector (proj.Mask(), proj.KeepValues()), memtype (PreferredMemType())
    {
      auto & ba = *Mask();
      size_t nwords = (ba.Size() + 31) / 32;
      std::vector<int> words (max<size_t>(nwords,1), 0);
      for (size_t i = 0; i < ba.Size(); i++)
        if (ba.Test(i)) words[i >> 5] |= (1 << (i & 31));
      dev_bits = DeviceDiagonalKernels<double>::Get().device->NewBuffer<int> (words.size(), MemType::Device);
      dev_bits.H2D (words.data(), words.size());
    }

    void Mult (const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("DeviceProjector::Mult"); RegionTimer reg(t);
      size_t n = Mask()->Size();
      if (n == 0) return;
      DeviceVectorWrapper<double> ux(x, memtype);
      DeviceVectorWrapper<double> uy(y, memtype);
      Launch (*DeviceDiagonalKernels<double>::Get().proj_mult, n,
              { dev_bits, ux.DevArgRO(), uy.DevArgW(),
                KernelArg(int(KeepValues())), KernelArg(int(n)) });
    }

    void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("DeviceProjector::MultAdd"); RegionTimer reg(t);
      size_t n = Mask()->Size();
      if (n == 0) return;
      DeviceVectorWrapper<double> ux(x, memtype);
      DeviceVectorWrapper<double> uy(y, memtype);
      Launch (*DeviceDiagonalKernels<double>::Get().proj_multadd, n,
              { dev_bits, ux.DevArgRO(), uy.DevArgRW(),
                KernelArg(s), KernelArg(int(KeepValues())), KernelArg(int(n)) });
    }

    void Project (BaseVector & x) const override
    {
      static Timer t("DeviceProjector::Project"); RegionTimer reg(t);
      size_t n = Mask()->Size();
      if (n == 0) return;
      DeviceVectorWrapper<double> ux(x, memtype);
      Launch (*DeviceDiagonalKernels<double>::Get().proj_project, n,
              { dev_bits, ux.DevArgRW(),
                KernelArg(int(KeepValues())), KernelArg(int(n)) });
    }

    AutoVector CreateRowVector () const override
    { return make_unique<DeviceVector<double>> (Mask()->Size(), memtype); }
    AutoVector CreateColVector () const override
    { return make_unique<DeviceVector<double>> (Mask()->Size(), memtype); }

    BaseMatrix::OperatorInfo GetOperatorInfo () const override
    { return { "DeviceProjector", Mask()->Size(), Mask()->Size() }; }
  };


  shared_ptr<BaseMatrix> Projector :: CreateDeviceMatrix () const
  {
    if (ngs_gpu::HasDevice())
      return make_shared<DeviceProjector> (*this);
    return BaseMatrix::CreateDeviceMatrix();
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
