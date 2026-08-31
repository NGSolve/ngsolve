#ifndef FILE_DEVICEVECTOR_HPP
#define FILE_DEVICEVECTOR_HPP

/*********************************************************************/
/* File:   devicevector.hpp                                          */
/* Date:   29. Aug. 2026                                             */
/*********************************************************************/

/*
  Backend-independent gpu vector, built on ngs_gpu (ngstd/gpuwrapper.hpp).

  Storage is one ngs_gpu::Buffer. With MemType::Shared the buffer is host
  addressable and host_data points into it, with MemType::Device a host
  mirror is kept and transferred on demand. Which side is current is
  tracked by host_uptodate / dev_uptodate; on shared memory both stay true
  and UpdateHost only waits for the queue.

  The elementwise kernels are written in the common syntax of
  ngstd/gpukernel.hpp and compiled at runtime, so the same source runs on
  metal, cuda and the cpu reference backend.
*/

#include <gpuwrapper.hpp>
#include "basevector.hpp"

namespace ngla
{
  using ngs_gpu::MemType;

  // the device gpu vectors live on - the registered backend, throws if none
  NGS_DLL_HEADER shared_ptr<ngs_gpu::Device> GetGpuDevice();

  // shared where it is free (unified memory), a device buffer otherwise
  NGS_DLL_HEADER MemType PreferredMemType();

  template <typename T> class DeviceVectorWrapper;
  class UnifiedVectorWrapper;   // cuda wrapper (ngscuda), aliases like DeviceVectorWrapper


  /*
    A scalar living in gpu memory, at a fixed device address. Kernels
    and backend libraries can read and write it without a host
    round-trip, which is what makes operator sequences using it
    recordable into a device graph.
  */
  class NGS_DLL_HEADER DeviceScalar : public BaseScalar
  {
  protected:
    shared_ptr<ngs_gpu::Buffer> devbuffer;   // one double
    shared_ptr<ngs_gpu::Queue> queue;

  public:
    DeviceScalar (double d = 0.0);

    void Set (double d) override;
    void Set (Complex c) override;
    double GetD () const override;
    Complex GetC () const override;

    // buffer of the value, for kernel arguments and backend access
    const shared_ptr<ngs_gpu::Buffer> & DevBuffer() const { return devbuffer; }
    ngs_gpu::KernelArg DevArg() const { return ngs_gpu::KernelArg(*devbuffer); }
  };


  template <typename T>
  class NGS_DLL_HEADER DeviceVector : public S_BaseVector<T>
  {
  protected:
    shared_ptr<ngs_gpu::Buffer> devbuffer;
    shared_ptr<ngs_gpu::Queue> queue;
    size_t devoffset = 0;              // elements into devbuffer
    size_t align = 1;
    MemType memtype = MemType::Shared;

    mutable T * host_data = nullptr;   // into devbuffer, into hostmem, or foreign
    Array<T> hostmem;                  // only if the host copy is ours and separate
    mutable bool host_uptodate = false, dev_uptodate = false;

    DeviceVector () { }
    void AllocBuffer (size_t asize);

    friend class DeviceVectorWrapper<T>;
    friend class UnifiedVectorWrapper;

  public:
    DeviceVector (size_t asize, MemType amemtype = MemType::Shared, size_t aalign = 1);
    DeviceVector (const BaseVector & v, MemType amemtype = MemType::Shared, size_t aalign = 1);
    virtual ~DeviceVector();

    using S_BaseVector<T>::Size;

    // buffer + offset for a kernel launch, with the transfers it implies
    ngs_gpu::KernelArg DevArgRO() const;   // kernel reads
    ngs_gpu::KernelArg DevArgRW() const;   // kernel reads and writes
    ngs_gpu::KernelArg DevArgW() const;    // kernel overwrites, no upload needed

    const shared_ptr<ngs_gpu::Buffer> & DevBufferRO() const
    { UpdateDevice(); return devbuffer; }
    const shared_ptr<ngs_gpu::Buffer> & DevBufferRW() const
    { UpdateDevice(); InvalidateHost(); return devbuffer; }
    size_t DevOffset() const { return devoffset; }
    MemType GetMemType() const { return memtype; }
    const shared_ptr<ngs_gpu::Queue> & GetQueue() const { return queue; }

    void UpdateHost () const;
    void UpdateDevice () const;
    // on shared memory host and device are the same storage, never stale
    void InvalidateHost() const { if (memtype != MemType::Shared) host_uptodate = false; }
    void InvalidateDevice() const { if (memtype != MemType::Shared) dev_uptodate = false; }
    bool IsHostUptodate() const { return host_uptodate; }
    bool IsDevUptodate() const { return dev_uptodate; }

    // read-only host view, does not invalidate the device copy
    const T * HostData() const { UpdateHost(); return host_data; }

    BaseVector & operator= (double d) { return SetScalar(d); }
    BaseVector & operator= (const BaseVector & v) { return Set(1.0, v); }
    DeviceVector & operator= (const DeviceVector & v) { Set(1.0, v); return *this; }

    template <typename T2>
    DeviceVector & operator= (const VVecExpr<T2> & v)
    { BaseVector::operator= (v); return *this; }

    virtual BaseVector & Scale (double scal) override;
    virtual BaseVector & Scale (BaseScalar & scal) override;
    virtual BaseVector & SetScalar (double scal) override;
    virtual BaseVector & Set (double scal, const BaseVector & v) override;
    virtual BaseVector & Add (double scal, const BaseVector & v) override;
    virtual BaseVector & Add (BaseScalar & scal, const BaseVector & v) override;

    virtual void * Memory () const override;
    virtual FlatVector<T> FVScal () const override;
    virtual AutoVector CreateVector () const override;
    virtual shared_ptr<BaseScalar> CreateScalar () const override;
    virtual ostream & Print (ostream & ost) const override;
  };


  /*
    Makes any BaseVector usable as a DeviceVector<T> for the lifetime of the
    wrapper. If the argument already is one its storage is aliased and its
    up-to-date flags are handed back on destruction, otherwise a device
    buffer is allocated and the values are transferred back - but only if
    the device actually wrote to them.
  */
  template <typename T>
  class NGS_DLL_HEADER DeviceVectorWrapper : public DeviceVector<T>
  {
    const BaseVector & vec;
    const DeviceVector<T> * alias_of = nullptr;
    bool initial_host_uptodate = false;
    bool converting = false;         // hostmem holds a converted copy of vec
  public:
    DeviceVectorWrapper (const BaseVector & avec, MemType amemtype = MemType::Shared);
    ~DeviceVectorWrapper();

    using DeviceVector<T>::operator=;
  };


#if !defined(FILE_DEVICEVECTOR_CPP)
  extern template class DeviceVector<double>;
  extern template class DeviceVector<float>;
  extern template class DeviceVectorWrapper<double>;
  extern template class DeviceVectorWrapper<float>;
#endif
}

#endif
