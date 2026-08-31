/*********************************************************************/
/* File:   devicevector.cpp                                          */
/* Date:   29. Aug. 2026                                             */
/*********************************************************************/

#define FILE_DEVICEVECTOR_CPP

#include <la.hpp>
#include <gpukernel.hpp>
#include <map>

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


      /*
        two-stage dot product: dv_dot1 writes one partial sum per group,
        dv_dot2 reduces the partials with a single group into res[0].
        Both stages stay on the device, so with res pointing into a
        DeviceScalar the reduction is graph-capturable.
      */
      KERNEL(dv_dot1, GLOBAL_IN(SCAL,x), GLOBAL_IN(SCAL,y), GLOBAL(SCAL,partial), VALUE(int,n))
      {
        SHARED(SCAL, tmp, 1024);
        SCAL s = 0;
        for (int i = int(GLOBAL_ID_X); i < n; i += int(GROUP_SIZE_X*NUM_GROUPS_X))
          s += x[i]*y[i];
        tmp[LOCAL_ID_X] = s;
        BARRIER();
        for (uint d = GROUP_SIZE_X/2; d > 0; d /= 2)
          {
            if (LOCAL_ID_X < d) tmp[LOCAL_ID_X] += tmp[LOCAL_ID_X+d];
            BARRIER();
          }
        if (LOCAL_ID_X == 0) partial[GROUP_ID_X] = tmp[0];
      }

      KERNEL(dv_dot2, GLOBAL_IN(SCAL,partial), GLOBAL(SCAL,res), VALUE(int,m))
      {
        SHARED(SCAL, tmp, 1024);
        SCAL s = 0;
        for (int i = int(LOCAL_ID_X); i < m; i += int(GROUP_SIZE_X))
          s += partial[i];
        tmp[LOCAL_ID_X] = s;
        BARRIER();
        for (uint d = GROUP_SIZE_X/2; d > 0; d /= 2)
          {
            if (LOCAL_ID_X < d) tmp[LOCAL_ID_X] += tmp[LOCAL_ID_X+d];
            BARRIER();
          }
        if (LOCAL_ID_X == 0) res[0] = tmp[0];
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
      shared_ptr<Kernel> dot1, dot2;
      shared_ptr<Buffer> dot_scratch;   // dotgroups partials + one result slot
      unsigned dotgroups, dotgroupsize;
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
        dot1      = library->GetKernel ("dv_dot1");
        dot2      = library->GetKernel ("dv_dot2");

        queue = device->DefaultQueue();
        // a wide group is what saturates memory on a gpu, but the cpu
        // reference backend runs one OS thread per work-item
        groupsize = (device->SimdWidth() > 1) ? device->MaxThreadsPerGroup() : 64;

        // the reduction tree wants a power of two, at most the shared size
        dotgroupsize = 1;
        while (2*dotgroupsize <= groupsize && 2*dotgroupsize <= 1024)
          dotgroupsize *= 2;
        // enough groups to saturate a gpu, few for the cpu reference backend
        dotgroups = (device->SimdWidth() > 1) ? 256 : 4;
        dot_scratch = device->NewBuffer ((dotgroups+1)*sizeof(T), MemType::Device);
      }

      // res: where dv_dot2 puts the result (buffer + byte offset)
      void LaunchDot (KernelArg x, KernelArg y, KernelArg res, int n) const
      {
        queue->Launch (*dot1, Dim3(dotgroups), Dim3(dotgroupsize),
                       { x, y, KernelArg(*dot_scratch), KernelArg(int(n)) });
        queue->Launch (*dot2, Dim3(1), Dim3(dotgroupsize),
                       { KernelArg(*dot_scratch), res, KernelArg(int(dotgroups)) });
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


    // v[offset : offset+dst.Size()) -> dst, converting the scalar types
    template <typename T>
    void CopyFromBaseVector (const BaseVector & v, FlatVector<T> dst, size_t offset = 0)
    {
      std::visit ([&] (auto proto)
      {
        auto src = v.FV<decltype(proto)>().Range(offset, offset+dst.Size());
        if constexpr (requires { dst(0) = src(0); })
          dst = src;
        else
          throw Exception("DeviceVector: cannot convert from complex vector");
      }, v.GetScalarType());
    }

    template <typename T>
    void CopyToBaseVector (FlatVector<T> src, const BaseVector & v, size_t offset = 0)
    {
      std::visit ([&] (auto proto)
      {
        auto dst = v.FV<decltype(proto)>().Range(offset, offset+src.Size());
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
  double DeviceVector<T> :: InnerProductD (const BaseVector & v2) const
  {
    if (v2.Size() != this->size)
      throw Exception("DeviceVector::InnerProductD - size mismatch");
    DeviceVectorWrapper<T> dv(v2, memtype);
    const auto & kern = DeviceVectorKernels<T>::Get();
    // the host reads back anyway, so skip the second reduction kernel:
    // fetch the partials (2KB cost the same as 8 bytes) and sum here
    kern.queue->Launch (*kern.dot1, Dim3(kern.dotgroups), Dim3(kern.dotgroupsize),
                        { DevArgRO(), dv.DevArgRO(),
                          KernelArg(*kern.dot_scratch), KernelArg(int(this->size)) });
    kern.queue->Finish();
    T partials[256];
    kern.dot_scratch->D2HArray (partials, kern.dotgroups);
    T res = 0;
    for (unsigned i = 0; i < kern.dotgroups; i++)
      res += partials[i];
    return res;
  }

  template <typename T>
  void DeviceVector<T> :: InnerProduct (const BaseVector & v2, BaseScalar & scal, bool conjugate) const
  {
    // result straight into the DeviceScalar - no host round-trip, capturable
    if constexpr (is_same_v<T,double>)
      if (auto devscal = dynamic_cast<DeviceScalar*> (&scal))
        {
          if (v2.Size() != this->size)
            throw Exception("DeviceVector::InnerProduct - size mismatch");
          DeviceVectorWrapper<T> dv(v2, memtype);
          const auto & kern = DeviceVectorKernels<T>::Get();
          kern.LaunchDot (DevArgRO(), dv.DevArgRO(), devscal->DevArg(), int(this->size));
          return;
        }
    BaseVector::InnerProduct (v2, scal, conjugate);
  }

  template <typename T>
  double DeviceVector<T> :: L2Norm () const
  {
    return sqrt (fabs (InnerProductD (*this)));
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
  DeviceVectorWrapper<T> :: DeviceVectorWrapper (const BaseVector & avec, MemType amemtype,
                                                 optional<IntRange> opt_range)
    : vec(avec)
  {
    IntRange range = opt_range.value_or (IntRange(0, avec.Size()));
    this->size = range.Size();
    subrange = (range.Size() != avec.Size());

    if (auto p = dynamic_cast<const DeviceVector<T>*> (&vec))
      {
        // alias the storage, take over its flags and hand them back later
        alias_of = p;
        // a device write to a sub-range must not leave the rest of the
        // vector stale, so bring the whole vector to the device first
        if (subrange)
          p->UpdateDevice();
        this->devbuffer = p->devbuffer;
        this->queue = p->queue;
        this->devoffset = p->devoffset + range.First();
        this->align = p->align;
        this->memtype = p->memtype;
        this->host_data = p->host_data ? p->host_data + range.First() : nullptr;
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
        this->host_data = vec.FV<T>().Data() + range.First();
      }
    else
      {
        converting = true;
        convert_offset = range.First();
        CopyFromBaseVector<T> (vec, FlatVector<T>(this->size, this->host_data), convert_offset);
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
        if (subrange)
          {
            // an update of the sub-range cannot validate the full
            // vector, but staleness propagates
            alias_of->host_uptodate &= this->host_uptodate;
            alias_of->dev_uptodate &= this->dev_uptodate;
          }
        else
          {
            // full alias: same storage, the flags simply carry over
            alias_of->host_uptodate = this->host_uptodate;
            alias_of->dev_uptodate = this->dev_uptodate;
          }
        this->host_data = nullptr;
        return;
      }

    // only pay for the copy back if a kernel actually wrote
    if (initial_host_uptodate && !this->host_uptodate)
      {
        this->UpdateHost();
        if (converting)
          CopyToBaseVector<T> (FlatVector<T>(this->size, this->host_data), vec, convert_offset);
      }
    this->host_data = nullptr;
  }


  template <typename T>
  AutoVector DeviceVector<T> :: Range (T_Range<size_t> range) const
  {
    return make_unique<DeviceVectorWrapper<T>>
      (*this, memtype, IntRange(range.First(), range.Next()));
  }


  template <typename T>
  shared_ptr<BaseScalar> DeviceVector<T> :: CreateScalar () const
  {
    return make_shared<DeviceScalar>();
  }


  void EvalScalarExpr (const std::string & params, const std::string & expr,
                       const std::vector<ngs_gpu::KernelArg> & args)
  {
    struct Entry
    {
      shared_ptr<Device> device;
      shared_ptr<Library> library;
      shared_ptr<Kernel> kernel;
    };
    static mutex mtx;
    static std::map<string, Entry> cache;

    auto dev = GetGpuDevice();
    string key = params + "|" + expr;

    shared_ptr<Kernel> kernel;
    {
      auto lock = lock_guard<mutex>(mtx);
      auto & e = cache[key];
      if (!e.kernel || e.device != dev)
        {
          string src = string(code_gpukernel) +
            "\nKERNEL(sc_expr, GLOBAL(double,r)" + params + ")\n" +
            "{ r[0] = " + expr + "; }\n";
          e.device = dev;
          e.library = dev->CompileSource (src);
          e.kernel = e.library->GetKernel ("sc_expr");
        }
      kernel = e.kernel;
    }
    dev->DefaultQueue()->Launch (*kernel, Dim3(1), Dim3(1), args);
  }


  DeviceScalar :: DeviceScalar (double d)
  {
    auto dev = GetGpuDevice();
    devbuffer = dev->NewBuffer (sizeof(double), PreferredMemType());
    queue = dev->DefaultQueue();
    devbuffer->H2D (&d, sizeof(double));
  }

  DeviceScalar :: DeviceScalar (const DeviceScalar & s2)
    : DeviceScalar ()
  {
    (*this) = s2;
  }

  DeviceScalar & DeviceScalar :: operator= (const DeviceScalar & s2)
  {
    std::vector<KernelArg> args = { DevArg(), KernelArg(*s2.devbuffer) };
    EvalScalarExpr (", GLOBAL_IN(double,s1)", "s1[0]", args);
    return *this;
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
