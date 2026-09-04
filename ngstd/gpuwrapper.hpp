#ifndef FILE_GPUWRAPPER
#define FILE_GPUWRAPPER

/*********************************************************************/
/* File:   gpuwrapper.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   29. Aug. 2026                                             */
/*********************************************************************/

/*
  Backend-independent GPU interface (CUDA / Metal).

  Kernels are generated as source (see tinybla), compiled at runtime,
  looked up by name, bound to positional arguments and dispatched.

  No backend headers here - ngscuda and ngsmetal implement the interface.
*/

#include <complex>
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
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


  /*
    Scalar type of a kernel argument (element type for buffers), for
    checking launches against the kernel declaration. size 0 = unknown,
    unknown types are not checked.
  */
  template <typename T> struct IsStdComplex : std::false_type { };
  template <typename T> struct IsStdComplex<std::complex<T>> : std::true_type { };

  struct ArgType
  {
    enum Kind : unsigned char { Unknown, Float, Int, UInt, Bool, Complex };
    Kind kind = Unknown;
    unsigned size = 0;

    constexpr ArgType () = default;
    constexpr ArgType (Kind k, unsigned sz) : kind(k), size(sz) { }

    template <typename T>
    static constexpr ArgType Of()
    {
      if constexpr (std::is_same_v<T,bool>) return { Bool, 1 };
      else if constexpr (std::is_floating_point_v<T>) return { Float, sizeof(T) };
      else if constexpr (IsStdComplex<T>::value) return { Complex, sizeof(T) };   // size of both parts
      else if constexpr (std::is_integral_v<T>)
        return { std::is_signed_v<T> ? Int : UInt, sizeof(T) };
      else return { Unknown, sizeof(T) };   // size is still checked
    }
    // c type name as written in a kernel source, e.g. "unsigned int"
    static ArgType FromName (const string & name);

    bool Known() const { return kind != Unknown && size != 0; }
    bool operator== (const ArgType & t2) const { return kind == t2.kind && size == t2.size; }
    bool operator!= (const ArgType & t2) const { return !(*this == t2); }
    string ToString() const;    // "float32", "int64", "complex64", "8 bytes", "?"
  };


  // declared arguments of a kernel, from KERNEL(name, GLOBAL_IN(T,x), VALUE(T,a), ...)
  struct KernelSignature
  {
    enum Kind : unsigned char { Unknown, Buffer, Value };
    struct Arg { Kind kind; ArgType type; string name; };
    std::vector<Arg> args;
  };

  // all KERNEL declarations in a source, keyed by kernel name
  std::map<string, KernelSignature> ParseKernelSignatures (const string & source);


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
    ArgType type;           // element type / value type, unknown for raw buffers
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
      : KernelArg (*b.Raw(), elemoff*sizeof(T)) { type = ArgType::Of<T>(); }

    template <typename T,
              typename = std::enable_if_t<std::is_trivially_copyable_v<T> &&
                                          !std::is_base_of_v<Buffer, T>>>
    KernelArg (const T & val)
      : kind(Kind::Value), nbytes(sizeof(T)), type(ArgType::Of<T>())
    {
      static_assert (sizeof(T) <= MAX_INLINE, "kernel argument too large");
      std::memcpy (inlinedata, &val, sizeof(T));
    }

    // raw bytes by value
    KernelArg (const void * data, size_t bytes, ArgType atype = {})
      : kind(Kind::Value), nbytes(bytes), type(atype)
    {
      if (!type.size) type.size = unsigned(bytes);
      if (bytes > MAX_INLINE) throw std::runtime_error ("kernel argument too large");
      std::memcpy (inlinedata, data, bytes);
    }

    // e.g. a raw buffer whose type the caller knows
    KernelArg & WithType (ArgType t) { type = t; return *this; }

    Kind GetKind() const { return kind; }
    ArgType GetType() const { return type; }
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
    // nullptr if the source had no KERNEL(...) declaration for it
    const KernelSignature * Signature() const;
  };


  class Library : public std::enable_shared_from_this<Library>
  {
    friend class Device;
    std::map<string, KernelSignature> signatures;   // parsed from the source
  protected:
    virtual shared_ptr<Kernel> DoGetKernel (const string & name) = 0;
  public:
    virtual ~Library() = default;

    const KernelSignature * FindSignature (const string & name) const
    {
      auto it = signatures.find(name);
      return (it != signatures.end()) ? &it->second : nullptr;
    }

    shared_ptr<Kernel> GetKernel (const string & name)
    {
      auto kernel = DoGetKernel (name);
      kernel->library = shared_from_this();
      return kernel;
    }
  };


  inline const KernelSignature * Kernel :: Signature() const
  { return library ? library->FindSignature(Name()) : nullptr; }

  // throws if the arguments do not match the kernel's declaration
  void CheckKernelArgs (const Kernel & kernel, const std::vector<KernelArg> & args);


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
    {
      CheckKernelArgs (kernel, args);
      DoLaunch (kernel, groups, groupsize, args, dynamic_group_memory);
    }

    void Launch (Kernel & kernel, Dim3 groups, Dim3 groupsize,
                 std::initializer_list<KernelArg> args,
                 size_t dynamic_group_memory = 0)
    { Launch (kernel, groups, groupsize, std::vector<KernelArg>(args), dynamic_group_memory); }

    void Launch (const shared_ptr<Kernel> & kernel, Dim3 groups, Dim3 groupsize,
                 std::initializer_list<KernelArg> args,
                 size_t dynamic_group_memory = 0)
    { Launch (*kernel, groups, groupsize, std::vector<KernelArg>(args), dynamic_group_memory); }

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

    // the source's KERNEL declarations are kept for checking launches
    shared_ptr<Library> CompileSource (const string & source)
    {
      auto lib = DoCompileSource (source);
      lib->signatures = ParseKernelSignatures (source);
      return lib;
    }

    virtual shared_ptr<Queue> DefaultQueue() = 0;

  protected:
    virtual shared_ptr<Buffer> DoNewBuffer (size_t bytes, MemType mt) = 0;
    virtual shared_ptr<Library> DoCompileSource (const string & source) = 0;
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
