/*********************************************************************/
/* File:   devicevector.cpp                                          */
/* Date:   29. Aug. 2026                                             */
/*********************************************************************/

#define FILE_DEVICEVECTOR_CPP

#include <la.hpp>
#include <gpukernel.hpp>

namespace ngla
{
  using namespace ngs_gpu;

  shared_ptr<Device> GetGpuDevice()
  {
    auto dev = ngs_gpu::GetDevice();
    if (!dev)
      throw Exception("no gpu backend registered - import ngsolve.ngsmetal or "
                      "ngsolve.ngscuda, or SetGPUDevice(GetCPUDevice())");
    return dev;
  }

  MemType PreferredMemType()
  {
    return GetGpuDevice()->IsUnifiedMemory() ? MemType::Shared : MemType::Device;
  }


  namespace
  {
    string Substitute (string src, const string & token, const string & value)
    {
      for (size_t p = src.find(token); p != string::npos; p = src.find(token, p))
        src.replace (p, token.size(), value);
      return src;
    }

    // elements handled by one work-item - the memory bound kernels want
    // more than one in flight, and the tail is taken by the last one only
    constexpr int DV_PER_THREAD = 4;

    const char * kernel_source = R"RAW(

      KERNEL(dv_setscalar, GLOBAL(SCAL,x), VALUE(SCAL,a), VALUE(int,n))
      {
        int i = 4*int(GLOBAL_ID_X);
        if (i+3 < n)
          { x[i] = a; x[i+1] = a; x[i+2] = a; x[i+3] = a; }
        else
          for ( ; i < n; i++) x[i] = a;
      }

      KERNEL(dv_scale, GLOBAL(SCAL,x), VALUE(SCAL,a), VALUE(int,n))
      {
        int i = 4*int(GLOBAL_ID_X);
        if (i+3 < n)
          { x[i] *= a; x[i+1] *= a; x[i+2] *= a; x[i+3] *= a; }
        else
          for ( ; i < n; i++) x[i] *= a;
      }

      // scale by a value living in gpu memory (a DeviceScalar) -
      // no host round-trip, so sequences using it can be captured
      KERNEL(dv_scale_dev, GLOBAL(SCAL,x), GLOBAL_IN(SCAL,a), VALUE(int,n))
      {
        SCAL s = a[0];
        int i = 4*int(GLOBAL_ID_X);
        if (i+3 < n)
          { x[i] *= s; x[i+1] *= s; x[i+2] *= s; x[i+3] *= s; }
        else
          for ( ; i < n; i++) x[i] *= s;
      }

      KERNEL(dv_set, GLOBAL(SCAL,y), GLOBAL_IN(SCAL,x), VALUE(SCAL,a), VALUE(int,n))
      {
        int i = 4*int(GLOBAL_ID_X);
        if (i+3 < n)
          { y[i] = a*x[i]; y[i+1] = a*x[i+1]; y[i+2] = a*x[i+2]; y[i+3] = a*x[i+3]; }
        else
          for ( ; i < n; i++) y[i] = a*x[i];
      }

      KERNEL(dv_axpy, GLOBAL(SCAL,y), GLOBAL_IN(SCAL,x), VALUE(SCAL,a), VALUE(int,n))
      {
        int i = 4*int(GLOBAL_ID_X);
        if (i+3 < n)
          { y[i] += a*x[i]; y[i+1] += a*x[i+1]; y[i+2] += a*x[i+2]; y[i+3] += a*x[i+3]; }
        else
          for ( ; i < n; i++) y[i] += a*x[i];
      }

      // axpy with the coefficient living in gpu memory (a DeviceScalar)
      KERNEL(dv_axpy_dev, GLOBAL(SCAL,y), GLOBAL_IN(SCAL,x), GLOBAL_IN(SCAL,a), VALUE(int,n))
      {
        SCAL s = a[0];
        int i = 4*int(GLOBAL_ID_X);
        if (i+3 < n)
          { y[i] += s*x[i]; y[i+1] += s*x[i+1]; y[i+2] += s*x[i+2]; y[i+3] += s*x[i+3]; }
        else
          for ( ; i < n; i++) y[i] += s*x[i];
      }

    )RAW";


    /*
      The elementwise kernels, compiled once per device and scalar type.
      Rebuilt if the registered backend changed underneath us.
    */
    template <typename T>
    class DeviceVectorKernels
    {
      shared_ptr<Library> library;
    public:
      shared_ptr<Device> device;
      shared_ptr<ngs_gpu::Queue> queue;
      shared_ptr<Kernel> setscalar, scale, scale_dev, set, axpy, axpy_dev;
      unsigned groupsize;

      DeviceVectorKernels (shared_ptr<Device> adevice)
        : device(adevice)
      {
        if constexpr (is_same_v<T,double>)
          if (!device->HasFloat64())
            throw Exception("DeviceVector<double> on "+device->Name()+
                            ", which has no fp64 - use DeviceVector<float>");

        string scal = is_same_v<T,double> ? "double" : "float";
        library = device->CompileSource (string(code_gpukernel) +
                                         Substitute (kernel_source, "SCAL", scal));
        setscalar = library->GetKernel ("dv_setscalar");
        scale     = library->GetKernel ("dv_scale");
        scale_dev = library->GetKernel ("dv_scale_dev");
        set       = library->GetKernel ("dv_set");
        axpy      = library->GetKernel ("dv_axpy");
        axpy_dev  = library->GetKernel ("dv_axpy_dev");

        queue = device->DefaultQueue();
        // a wide group is what saturates memory on a gpu, but the cpu
        // reference backend runs one OS thread per work-item
        groupsize = (device->SimdWidth() > 1) ? device->MaxThreadsPerGroup() : 64;
      }

      static const DeviceVectorKernels & Get()
      {
        static mutex mtx;
        static shared_ptr<DeviceVectorKernels> cached;

        auto dev = GetGpuDevice();
        auto lock = lock_guard<mutex>(mtx);
        if (!cached || cached->device != dev)
          cached = make_shared<DeviceVectorKernels> (dev);
        return *cached;
      }

      void Launch (const shared_ptr<Kernel> & kernel, size_t n,
                   std::initializer_list<KernelArg> args) const
      {
        if (n == 0) return;
        size_t items = (n + DV_PER_THREAD-1) / DV_PER_THREAD;
        unsigned groups = (items + groupsize - 1) / groupsize;
        queue->Launch (*kernel, Dim3(groups), Dim3(groupsize), args);
      }
    };


    // v -> dst, converting between the scalar types
    template <typename T>
    void CopyFromBaseVector (const BaseVector & v, FlatVector<T> dst)
    {
      std::visit ([&] (auto proto)
      {
        auto src = v.FV<decltype(proto)>();
        if constexpr (requires { dst(0) = src(0); })
          dst = src;
        else
          throw Exception("DeviceVector: cannot convert from complex vector");
      }, v.GetScalarType());
    }

    template <typename T>
    void CopyToBaseVector (FlatVector<T> src, const BaseVector & v)
    {
      std::visit ([&] (auto proto)
      {
        auto dst = v.FV<decltype(proto)>();
        if constexpr (requires { dst(0) = src(0); })
          dst = src;
        else
          throw Exception("DeviceVector: cannot convert to complex vector");
      }, v.GetScalarType());
    }

    template <typename T>
    bool HasScalarType (const BaseVector & v)
    {
      return std::visit ([] (auto proto)
                         { return is_same_v<decltype(proto),T>; },
                         v.GetScalarType());
    }
  }



  template <typename T>
  void DeviceVector<T> :: AllocBuffer (size_t asize)
  {
    const auto & kern = DeviceVectorKernels<T>::Get();
    queue = kern.queue;

    size_t padded = align > 1 ? ((asize+align-1)/align)*align : asize;
    devbuffer = kern.device->NewBuffer (padded*sizeof(T), memtype);
    devoffset = 0;

    if (memtype == MemType::Shared)
      {
        host_data = devbuffer->HostData<T>();
        host_uptodate = dev_uptodate = true;
      }
    else
      {
        hostmem.SetSize (asize);
        host_data = hostmem.Data();
        host_uptodate = dev_uptodate = false;
      }

    if (padded > asize)
      {
        Array<T> zeros(padded-asize);
        zeros = T(0);
        devbuffer->H2DArray<T> (zeros.Data(), padded-asize, asize);
      }
  }

  template <typename T>
  DeviceVector<T> :: DeviceVector (size_t asize, MemType amemtype, size_t aalign)
    : align(aalign), memtype(amemtype)
  {
    this->size = asize;
    AllocBuffer (asize);
  }

  template <typename T>
  DeviceVector<T> :: DeviceVector (const BaseVector & v, MemType amemtype, size_t aalign)
    : align(aalign), memtype(amemtype)
  {
    this->size = v.Size();
    AllocBuffer (this->size);

    CopyFromBaseVector<T> (v, FlatVector<T>(this->size, host_data));
    host_uptodate = true;
    InvalidateDevice();
    UpdateDevice();
  }

  template <typename T>
  DeviceVector<T> :: ~DeviceVector() { }


  template <typename T>
  KernelArg DeviceVector<T> :: DevArgRO() const
  {
    UpdateDevice();
    return KernelArg (*devbuffer, devoffset*sizeof(T));
  }

  template <typename T>
  KernelArg DeviceVector<T> :: DevArgRW() const
  {
    UpdateDevice();
    InvalidateHost();
    return KernelArg (*devbuffer, devoffset*sizeof(T));
  }

  template <typename T>
  KernelArg DeviceVector<T> :: DevArgW() const
  {
    dev_uptodate = true;
    InvalidateHost();
    return KernelArg (*devbuffer, devoffset*sizeof(T));
  }


  template <typename T>
  void DeviceVector<T> :: UpdateHost () const
  {
    // a kernel may still be writing the buffer we are about to read
    if (queue) queue->Finish();

    if (host_uptodate) return;
    if (dev_uptodate)
      devbuffer->D2HArray<T> (host_data, this->size, devoffset);
    host_uptodate = true;
  }

  template <typename T>
  void DeviceVector<T> :: UpdateDevice () const
  {
    if (dev_uptodate) return;
    if (host_uptodate)
      devbuffer->H2DArray<T> (host_data, this->size, devoffset);
    dev_uptodate = true;
  }


  template <typename T>
  BaseVector & DeviceVector<T> :: Scale (double scal)
  {
    if (scal == 1.0) return *this;
    const auto & kern = DeviceVectorKernels<T>::Get();
    kern.Launch (kern.scale, this->size,
                 { DevArgRW(), KernelArg(T(scal)), KernelArg(int(this->size)) });
    return *this;
  }

  template <typename T>
  BaseVector & DeviceVector<T> :: Scale (BaseScalar & scal)
  {
    // the DeviceScalar holds a double, so the buffer types match only
    // for T=double; float vectors read the value over the host
    if constexpr (is_same_v<T,double>)
      if (auto devscal = dynamic_cast<DeviceScalar*> (&scal))
        {
          const auto & kern = DeviceVectorKernels<T>::Get();
          kern.Launch (kern.scale_dev, this->size,
                       { DevArgRW(), devscal->DevArg(), KernelArg(int(this->size)) });
          return *this;
        }
    return Scale (scal.GetD());
  }

  template <typename T>
  BaseVector & DeviceVector<T> :: SetScalar (double scal)
  {
    const auto & kern = DeviceVectorKernels<T>::Get();
    kern.Launch (kern.setscalar, this->size,
                 { DevArgW(), KernelArg(T(scal)), KernelArg(int(this->size)) });
    return *this;
  }

  template <typename T>
  BaseVector & DeviceVector<T> :: Set (double scal, const BaseVector & v)
  {
    if (v.Size() != this->size)
      throw Exception("DeviceVector::Set - size mismatch");
    if (&v == this) return Scale(scal);
    DeviceVectorWrapper<T> dv(v, memtype);
    const auto & kern = DeviceVectorKernels<T>::Get();
    kern.Launch (kern.set, this->size,
                 { DevArgW(), dv.DevArgRO(),
                   KernelArg(T(scal)), KernelArg(int(this->size)) });
    return *this;
  }

  template <typename T>
  BaseVector & DeviceVector<T> :: Add (double scal, const BaseVector & v)
  {
    if (v.Size() != this->size)
      throw Exception("DeviceVector::Add - size mismatch");
    if (&v == this) return Scale(1.0+scal);
    DeviceVectorWrapper<T> dv(v, memtype);
    const auto & kern = DeviceVectorKernels<T>::Get();
    kern.Launch (kern.axpy, this->size,
                 { DevArgRW(), dv.DevArgRO(),
                   KernelArg(T(scal)), KernelArg(int(this->size)) });
    return *this;
  }


  template <typename T>
  BaseVector & DeviceVector<T> :: Add (BaseScalar & scal, const BaseVector & v)
  {
    // buffer types match only for T=double, see Scale
    if constexpr (is_same_v<T,double>)
      if (auto devscal = dynamic_cast<DeviceScalar*> (&scal))
        {
          if (v.Size() != this->size)
            throw Exception("DeviceVector::Add - size mismatch");
          DeviceVectorWrapper<T> dv(v, memtype);
          const auto & kern = DeviceVectorKernels<T>::Get();
          kern.Launch (kern.axpy_dev, this->size,
                       { DevArgRW(), dv.DevArgRO(),
                         devscal->DevArg(), KernelArg(int(this->size)) });
          return *this;
        }
    return Add (scal.GetD(), v);
  }


  template <typename T>
  void * DeviceVector<T> :: Memory () const
  {
    // the caller gets a writable pointer, so the device copy may go stale
    UpdateHost();
    InvalidateDevice();
    return host_data;
  }

  template <typename T>
  FlatVector<T> DeviceVector<T> :: FVScal () const
  {
    UpdateHost();
    InvalidateDevice();
    return FlatVector<T> (this->size, host_data);
  }

  template <typename T>
  AutoVector DeviceVector<T> :: CreateVector () const
  {
    return make_unique<DeviceVector<T>> (this->size, memtype, align);
  }

  template <typename T>
  ostream & DeviceVector<T> :: Print (ostream & ost) const
  {
    const T * data = HostData();
    ost << "DeviceVector<" << (is_same_v<T,double> ? "double" : "float")
        << ">, size = " << this->size << endl;
    for (size_t i = 0; i < this->size; i++)
      ost << data[i] << "\n";
    return ost;
  }



  template <typename T>
  DeviceVectorWrapper<T> :: DeviceVectorWrapper (const BaseVector & avec, MemType amemtype)
    : vec(avec)
  {
    this->size = vec.Size();

    if (auto p = dynamic_cast<const DeviceVector<T>*> (&vec))
      {
        // alias the storage, take over its flags and hand them back later
        alias_of = p;
        this->devbuffer = p->devbuffer;
        this->queue = p->queue;
        this->devoffset = p->devoffset;
        this->align = p->align;
        this->memtype = p->memtype;
        this->host_data = p->host_data;
        this->host_uptodate = p->host_uptodate;
        this->dev_uptodate = p->dev_uptodate;
        return;
      }

    // a device buffer of our own, the host side stays with vec
    this->memtype = MemType::Device;
    this->align = 1;
    this->AllocBuffer (this->size);

    if (HasScalarType<T> (vec))
      {
        this->hostmem.SetSize(0);
        this->host_data = vec.FV<T>().Data();
      }
    else
      {
        converting = true;
        CopyFromBaseVector<T> (vec, FlatVector<T>(this->size, this->host_data));
      }

    initial_host_uptodate = true;
    this->host_uptodate = true;
    this->dev_uptodate = false;
  }

  template <typename T>
  DeviceVectorWrapper<T> :: ~DeviceVectorWrapper()
  {
    if (alias_of)
      {
        alias_of->host_uptodate = this->host_uptodate;
        alias_of->dev_uptodate = this->dev_uptodate;
        this->host_data = nullptr;
        return;
      }

    // only pay for the copy back if a kernel actually wrote
    if (initial_host_uptodate && !this->host_uptodate)
      {
        this->UpdateHost();
        if (converting)
          CopyToBaseVector<T> (FlatVector<T>(this->size, this->host_data), vec);
      }
    this->host_data = nullptr;
  }


  template <typename T>
  shared_ptr<BaseScalar> DeviceVector<T> :: CreateScalar () const
  {
    return make_shared<DeviceScalar>();
  }


  DeviceScalar :: DeviceScalar (double d)
  {
    auto dev = GetGpuDevice();
    devbuffer = dev->NewBuffer (sizeof(double), PreferredMemType());
    queue = dev->DefaultQueue();
    devbuffer->H2D (&d, sizeof(double));
  }

  void DeviceScalar :: Set (double d)
  {
    devbuffer->H2D (&d, sizeof(double));
  }

  void DeviceScalar :: Set (Complex c)
  {
    throw Exception ("DeviceScalar is real-valued");
  }

  double DeviceScalar :: GetD () const
  {
    queue->Finish();   // kernels writing the value may still be queued
    double d;
    devbuffer->D2H (&d, sizeof(double));
    return d;
  }

  Complex DeviceScalar :: GetC () const
  {
    throw Exception ("DeviceScalar is real-valued");
  }


  template class DeviceVector<double>;
  template class DeviceVector<float>;
  template class DeviceVectorWrapper<double>;
  template class DeviceVectorWrapper<float>;
}
