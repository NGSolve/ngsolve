#ifndef FILE_GPUWRAPPER
#define FILE_GPUWRAPPER

/*********************************************************************/
/* File:   gpuwrapper.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   29. Aug. 2025                                             */
/*********************************************************************/

/*
  Backend-independent GPU interface (CUDA / Metal).

  Kernels are generated as source (see tinybla), compiled at runtime,
  looked up by name, bound to positional arguments and dispatched.

  No backend headers here - ngscuda and ngsmetal implement the interface.
*/

#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <type_traits>
#include <initializer_list>
#include <ostream>

namespace ngs_gpu
{
  using std::string;
  using std::shared_ptr;

  enum class MemType { Device,    // not host visible
                       Shared };  // host + device (unified / managed)


  struct Dim3
  {
    unsigned x, y, z;
    constexpr Dim3 (unsigned ax = 1, unsigned ay = 1, unsigned az = 1)
      : x(ax), y(ay), z(az) { }
  };


  // public interface is non-virtual, backends override the Do* hooks,
  // otherwise the typed overloads get hidden in derived classes
  class Buffer
  {
  protected:
    size_t size;   // bytes
    MemType memtype;
    Buffer (size_t asize, MemType amemtype) : size(asize), memtype(amemtype) { }

    virtual void * DoHostPtr() const = 0;
    virtual void DoH2D (const void * src, size_t bytes, size_t offset) = 0;
    virtual void DoD2H (void * dst, size_t bytes, size_t offset) const = 0;

  public:
    virtual ~Buffer() = default;

    size_t Size() const { return size; }
    MemType GetMemType() const { return memtype; }

    // nullptr if not host visible
    void * HostPtr() const { return DoHostPtr(); }
    template <typename T> T * HostData() const { return static_cast<T*>(DoHostPtr()); }

    void H2D (const void * src, size_t bytes, size_t offset = 0)
    { DoH2D (src, bytes, offset); }
    void D2H (void * dst, size_t bytes, size_t offset = 0) const
    { DoD2H (dst, bytes, offset); }
  };


  /*
    Typed view of a buffer: sizes, offsets and transfers count elements of T.
    A value type sharing ownership of the raw buffer, copies alias.
  */
  template <typename T>
  class TypedBuffer
  {
    shared_ptr<Buffer> buf;
  public:
    typedef T value_type;

    TypedBuffer () = default;
    explicit TypedBuffer (shared_ptr<Buffer> abuf) : buf(std::move(abuf)) { }

    explicit operator bool() const { return bool(buf); }
    const shared_ptr<Buffer> & Raw() const { return buf; }

    size_t Size() const { return buf->Size() / sizeof(T); }
    size_t Bytes() const { return buf->Size(); }
    MemType GetMemType() const { return buf->GetMemType(); }

    // nullptr if not host visible
    T * HostData() const { return buf->HostData<T>(); }
    T & operator[] (size_t i) const { return HostData()[i]; }

    // const like a shared_ptr: the handle, not the memory
    void H2D (const T * src, size_t n, size_t offset = 0) const
    { buf->H2D (src, n*sizeof(T), offset*sizeof(T)); }
    void D2H (T * dst, size_t n, size_t offset = 0) const
    { buf->D2H (dst, n*sizeof(T), offset*sizeof(T)); }
  };


  // one positional kernel argument: a buffer binding, or bytes by value
  class KernelArg
  {
  public:
    enum class Kind { Buffer, Value };
    static constexpr size_t MAX_INLINE = 64;

  private:
    Kind kind;
    Buffer * buf = nullptr;
    size_t offset = 0;      // byte offset into buf
    size_t nbytes = 0;
    alignas(16) unsigned char inlinedata[MAX_INLINE];

  public:
    KernelArg (Buffer & b, size_t off = 0)
      : kind(Kind::Buffer), buf(&b), offset(off) { }

    // templated so shared_ptr<Derived> needs only one user-defined conversion
    template <typename TB, typename = std::enable_if_t<std::is_base_of_v<Buffer, TB>>>
    KernelArg (const shared_ptr<TB> & b, size_t off = 0)
      : KernelArg (static_cast<Buffer&>(*b), off) { }

    // offset in elements
    template <typename T>
    KernelArg (const TypedBuffer<T> & b, size_t elemoff = 0)
      : KernelArg (*b.Raw(), elemoff*sizeof(T)) { }

    template <typename T,
              typename = std::enable_if_t<std::is_trivially_copyable_v<T> &&
                                          !std::is_base_of_v<Buffer, T>>>
    KernelArg (const T & val)
      : kind(Kind::Value), nbytes(sizeof(T))
    {
      static_assert (sizeof(T) <= MAX_INLINE, "kernel argument too large");
      std::memcpy (inlinedata, &val, sizeof(T));
    }

    // raw bytes by value
    KernelArg (const void * data, size_t bytes)
      : kind(Kind::Value), nbytes(bytes)
    {
      if (bytes > MAX_INLINE) throw std::runtime_error ("kernel argument too large");
      std::memcpy (inlinedata, data, bytes);
    }

    Kind GetKind() const { return kind; }
    Buffer * GetBuffer() const { return buf; }
    size_t Offset() const { return offset; }
    size_t NBytes() const { return nbytes; }
    const void * Data() const { return inlinedata; }
  };


  class Library;

  // resource usage of a compiled kernel, 0 where the backend has no query
  struct KernelInfo
  {
    size_t registers = 0;            // per thread
    size_t local_bytes = 0;          // stack / spills per thread
    size_t shared_bytes = 0;         // static group memory
    size_t max_threads_per_group = 0;
    size_t max_groups_per_unit = 0;  // resident groups per SM/core at groupsize
  };

  inline std::ostream & operator<< (std::ostream & ost, const KernelInfo & ki)
  {
    return ost << "regs=" << ki.registers << " local=" << ki.local_bytes
               << "B shared=" << ki.shared_bytes << "B maxthreads=" << ki.max_threads_per_group
               << " groups/unit=" << ki.max_groups_per_unit;
  }

  class Kernel
  {
    friend class Library;
    shared_ptr<Library> library;   // the raw handle points into the module
  public:
    virtual ~Kernel() = default;
    virtual string Name() const = 0;
    const shared_ptr<Library> & GetLibrary() const { return library; }
    virtual KernelInfo Info (size_t groupsize = 0) const { return {}; }
  };


  class Library : public std::enable_shared_from_this<Library>
  {
  protected:
    virtual shared_ptr<Kernel> DoGetKernel (const string & name) = 0;
  public:
    virtual ~Library() = default;

    shared_ptr<Kernel> GetKernel (const string & name)
    {
      auto kernel = DoGetKernel (name);
      kernel->library = shared_from_this();
      return kernel;
    }
  };


  class Queue
  {
  protected:
    virtual void DoLaunch (Kernel & kernel, Dim3 groups, Dim3 groupsize,
                           const std::vector<KernelArg> & args,
                           size_t dynamic_group_memory) = 0;
  public:
    virtual ~Queue() = default;

    void Launch (Kernel & kernel, Dim3 groups, Dim3 groupsize,
                 const std::vector<KernelArg> & args,
                 size_t dynamic_group_memory = 0)
    { DoLaunch (kernel, groups, groupsize, args, dynamic_group_memory); }

    void Launch (Kernel & kernel, Dim3 groups, Dim3 groupsize,
                 std::initializer_list<KernelArg> args,
                 size_t dynamic_group_memory = 0)
    { DoLaunch (kernel, groups, groupsize, args, dynamic_group_memory); }

    void Launch (const shared_ptr<Kernel> & kernel, Dim3 groups, Dim3 groupsize,
                 std::initializer_list<KernelArg> args,
                 size_t dynamic_group_memory = 0)
    { DoLaunch (*kernel, groups, groupsize, args, dynamic_group_memory); }

    virtual void Finish() = 0;    // submit and wait
  };


  class Device
  {
  public:
    virtual ~Device() = default;

    virtual string Name() const = 0;
    virtual bool HasFloat64() const = 0;
    virtual bool IsUnifiedMemory() const = 0;
    virtual size_t MaxThreadsPerGroup() const = 0;
    virtual size_t SimdWidth() const = 0;      // warp / simdgroup
    // parallel processors (cuda SMs, apple gpu cores, cpu threads);
    // metal has no query for it, the default suits apple silicon
    virtual size_t ComputeUnits() const { return 16; }

    // raw bytes
    shared_ptr<Buffer> NewBuffer (size_t bytes, MemType mt = MemType::Device)
    { return DoNewBuffer (bytes, mt); }
    // n elements of T
    template <typename T>
    TypedBuffer<T> NewBuffer (size_t n, MemType mt = MemType::Device)
    { return TypedBuffer<T> (DoNewBuffer (n*sizeof(T), mt)); }

    virtual shared_ptr<Library> CompileSource (const string & source) = 0;
    virtual shared_ptr<Queue> DefaultQueue() = 0;

  protected:
    virtual shared_ptr<Buffer> DoNewBuffer (size_t bytes, MemType mt) = 0;
  };


  /*
    The backend (ngsmetal / ngscuda) registers its creator, ngstd hands out
    the device. Storage is in gpuwrapper.cpp: the backends are pybind modules
    built with hidden visibility, an inline static would not be shared.
  */
  typedef std::function<shared_ptr<Device>()> DeviceCreator;

  // reference implementation on the host, always available.
  // kernels are compiled with the host c++ compiler, see gpukernel.hpp
  shared_ptr<Device> GetCpuDevice();

  void SetDeviceCreator (DeviceCreator creator);
  bool HasDevice();
  shared_ptr<Device> GetDevice();     // nullptr if no backend registered
}

#endif
