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

  T is double, float or Complex (= std::complex<double>, the kernels see
  it as Complex<double>). The complex vector needs an fp64 device.
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
  template <typename T> class DeviceScalar;


  // compiles (once per shape) and launches the one-thread kernel
  // { <body> } - the backend of the scalar expressions below
  NGS_DLL_HEADER void EvalScalarExpr (const std::string & params, const std::string & body,
                                      const std::vector<ngs_gpu::KernelArg> & args);

  // the precision as spelled in a kernel
  template <typename T> constexpr const char * DevScalTypeName ()
  { return std::is_same_v<T,double> ? "double" : "float"; }

  template <typename T>
  concept IsDevScalarExpr = requires (const T & e, std::string & c, std::string & p,
                                      std::vector<ngs_gpu::KernelArg> & a)
  { e.Emit (c, p, a); typename T::value_type; };


  /*
    A scalar living in gpu memory, at a fixed device address, in the
    precision T of the vectors it belongs to (double or float). Kernels
    and backend libraries can read and write it without a host
    round-trip, which is what makes operator sequences using it
    recordable into a device graph.
  */
  template <typename T>
  class NGS_DLL_HEADER DeviceScalar : public BaseScalar
  {
  protected:
    ngs_gpu::TypedBuffer<T> devbuffer;   // one value
    shared_ptr<ngs_gpu::Queue> queue;

  public:
    typedef T value_type;

    DeviceScalar (double d = 0.0);

    // value copy on the device, no host round-trip
    DeviceScalar (const DeviceScalar & s2);
    DeviceScalar & operator= (const DeviceScalar & s2);
    DeviceScalar & operator= (double d) { Set(d); return *this; }

    void Set (double d) override;
    void Set (Complex c) override;
    double GetD () const override;
    Complex GetC () const override;

    // buffer of the value, for kernel arguments and backend access
    const shared_ptr<ngs_gpu::Buffer> & DevBuffer() const { return devbuffer.Raw(); }
    ngs_gpu::KernelArg DevArg() const { return ngs_gpu::KernelArg(devbuffer); }
    const shared_ptr<ngs_gpu::Queue> & GetQueue() const { return queue; }

    // evaluate a scalar expression on the device, e.g.
    //   alpha = rz / pq;
    // one single-thread kernel, no host round-trip, capturable
    template <IsDevScalarExpr Expr>
    DeviceScalar & operator= (const Expr & e)
    {
      static_assert (std::is_same_v<typename Expr::value_type, T>,
                     "DeviceScalar: expression of another precision");
      std::vector<ngs_gpu::KernelArg> args;
      std::string body = "s0[0] = ";
      std::string params = std::string(", GLOBAL(") + DevScalTypeName<T>() + ",s0)";
      args.push_back (DevArg());
      e.Emit (body, params, args);
      body += ";";
      EvalScalarExpr (params, body, args);
      return *this;
    }
  };


  /*
    Scalar expression templates. A node contributes its sub-expression
    to the kernel source and its buffers/values to the argument list;
    the assignment above compiles the collected source once per shape.
    All nodes of one expression share the precision T.
  */

  template <typename T>
  struct DevScalExpr        // leaf: a DeviceScalar
  {
    typedef T value_type;
    const DeviceScalar<T> * s;
    void Emit (std::string & code, std::string & params,
               std::vector<ngs_gpu::KernelArg> & args) const
    {
      auto name = "s" + std::to_string(args.size());
      params += std::string(", GLOBAL_IN(") + DevScalTypeName<T>() + "," + name + ")";
      code += name + "[0]";
      args.push_back (s->DevArg());
    }
  };

  template <typename T>
  struct DevConstExpr       // leaf: a value, passed as kernel argument
  {
    typedef T value_type;
    T val;
    void Emit (std::string & code, std::string & params,
               std::vector<ngs_gpu::KernelArg> & args) const
    {
      auto name = "s" + std::to_string(args.size());
      params += std::string(", VALUE(") + DevScalTypeName<T>() + "," + name + ")";
      code += name;
      args.push_back (ngs_gpu::KernelArg(val));
    }
  };

  template <IsDevScalarExpr L, IsDevScalarExpr R>
  struct DevBinExpr
  {
    typedef typename L::value_type value_type;
    L l; R r; char op;
    void Emit (std::string & code, std::string & params,
               std::vector<ngs_gpu::KernelArg> & args) const
    {
      code += "(";
      l.Emit (code, params, args);
      code += op;
      r.Emit (code, params, args);
      code += ")";
    }
  };

  template <IsDevScalarExpr E>
  struct DevFuncExpr
  {
    typedef typename E::value_type value_type;
    E e; const char * func;
    void Emit (std::string & code, std::string & params,
               std::vector<ngs_gpu::KernelArg> & args) const
    {
      code += std::string(func) + "(";
      e.Emit (code, params, args);
      code += ")";
    }
  };

  /*
    Operator syntax:  alpha = rz / pq;  neg_alpha = -alpha;
                      tau = Sqrt(rho);  x = 0.5*a + b;
    Operands are expressions, DeviceScalars or plain numbers; at least
    one gpu-flavoured operand is required, so the templates never touch
    ordinary arithmetic. The precision is that of the gpu-flavoured
    operands, which must agree; plain numbers are converted to it.
  */

  template <typename T>
  concept IsDevScalarLike = requires { typename std::remove_cvref_t<T>::value_type; } &&
    std::is_base_of_v<DeviceScalar<typename std::remove_cvref_t<T>::value_type>, std::remove_cvref_t<T>>;
  template <typename T>
  concept IsDevScalarOperand = IsDevScalarExpr<std::remove_cvref_t<T>> ||
    IsDevScalarLike<T> || std::is_arithmetic_v<std::remove_cvref_t<T>>;
  template <typename L, typename R>
  concept IsDevScalarOp = IsDevScalarOperand<L> && IsDevScalarOperand<R> &&
    !(std::is_arithmetic_v<std::remove_cvref_t<L>> && std::is_arithmetic_v<std::remove_cvref_t<R>>);

  // the precision of a binary operation: from the gpu-flavoured side
  template <typename L, typename R>
  using DevOpType = std::conditional_t<std::is_arithmetic_v<std::remove_cvref_t<L>>,
                                       typename std::remove_cvref_t<R>::value_type,
                                       typename std::remove_cvref_t<L>::value_type>;

  template <typename T, typename X>
  auto ToDevExpr (const X & x)
  {
    if constexpr (std::is_arithmetic_v<X>)
      return DevConstExpr<T>{T(x)};
    else if constexpr (IsDevScalarLike<X>)
      {
        static_assert (std::is_same_v<typename X::value_type, T>,
                       "DeviceScalar expression mixes precisions");
        return DevScalExpr<T>{&x};
      }
    else
      {
        static_assert (std::is_same_v<typename X::value_type, T>,
                       "DeviceScalar expression mixes precisions");
        return x;
      }
  }

#define NGS_DEVSCAL_BINOP(OP, CH)                                                   \
  template <typename L, typename R> requires IsDevScalarOp<L,R>                     \
  auto operator OP (const L & l, const R & r)                                       \
  {                                                                                 \
    typedef DevOpType<L,R> T;                                                       \
    auto el = ToDevExpr<T>(l); auto er = ToDevExpr<T>(r);                           \
    return DevBinExpr<decltype(el), decltype(er)>{el, er, CH};                      \
  }
  NGS_DEVSCAL_BINOP(+, '+')
  NGS_DEVSCAL_BINOP(-, '-')
  NGS_DEVSCAL_BINOP(*, '*')
  NGS_DEVSCAL_BINOP(/, '/')
#undef NGS_DEVSCAL_BINOP

  template <typename X> requires (IsDevScalarExpr<std::remove_cvref_t<X>> || IsDevScalarLike<X>)
  auto operator- (const X & x)
  {
    typedef typename std::remove_cvref_t<X>::value_type T;
    auto e = ToDevExpr<T>(x);
    return DevFuncExpr<decltype(e)>{e, "-"};
  }

  template <typename X> requires (IsDevScalarExpr<std::remove_cvref_t<X>> || IsDevScalarLike<X>)
  auto Sqrt (const X & x)
  {
    typedef typename std::remove_cvref_t<X>::value_type T;
    auto e = ToDevExpr<T>(x);
    return DevFuncExpr<decltype(e)>{e, "sqrt"};
  }


  /*
    Several assignments fused into one kernel launch:
      Eval (Assign(alpha, rz/pq), Assign(neg_alpha, -alpha));
    The statements run in order in a single thread, so a later one may
    read a scalar a previous one wrote.
  */
  template <IsDevScalarExpr E>
  struct DevScalStmt
  {
    typedef typename E::value_type T;
    DeviceScalar<T> * dest;
    E e;
    void EmitStmt (std::string & body, std::string & params,
                   std::vector<ngs_gpu::KernelArg> & args) const
    {
      auto name = "s" + std::to_string(args.size());
      params += std::string(", GLOBAL(") + DevScalTypeName<T>() + "," + name + ")";
      args.push_back (dest->DevArg());
      body += name + "[0] = ";
      e.Emit (body, params, args);
      body += "; ";
    }
  };

  template <typename T, typename E> requires IsDevScalarOperand<E>
  auto Assign (DeviceScalar<T> & d, const E & e)
  {
    auto ee = ToDevExpr<T>(e);
    return DevScalStmt<decltype(ee)>{&d, ee};
  }

  template <typename... E>
  void Eval (const DevScalStmt<E> &... stmts)
  {
    std::string body, params;
    std::vector<ngs_gpu::KernelArg> args;
    (stmts.EmitStmt (body, params, args), ...);
    EvalScalarExpr (params, body, args);
  }


  template <typename T>
  class NGS_DLL_HEADER DeviceVector : public S_BaseVector<T>
  {
  protected:
    ngs_gpu::TypedBuffer<T> devbuffer;
    shared_ptr<ngs_gpu::Queue> queue;
    size_t devoffset = 0;              // elements into devbuffer

    // the scaled operations, with the coefficient in the vector's type
    BaseVector & ScaleT (T scal);
    BaseVector & SetScalarT (T scal);
    BaseVector & SetT (T scal, const BaseVector & v);
    BaseVector & AddT (T scal, const BaseVector & v);
    T DotT (const BaseVector & v2, bool conjugate) const;   // result over the host
    size_t align = 1;
    MemType memtype = MemType::Shared;

    mutable T * host_data = nullptr;   // into devbuffer, into hostmem, or foreign
    Array<T> hostmem;                  // only if the host copy is ours and separate
    mutable bool host_uptodate = false, dev_uptodate = false;

    DeviceVector () { }
    void AllocBuffer (size_t asize);

    friend class DeviceVectorWrapper<T>;

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
    { UpdateDevice(); return devbuffer.Raw(); }
    const shared_ptr<ngs_gpu::Buffer> & DevBufferRW() const
    { UpdateDevice(); InvalidateHost(); return devbuffer.Raw(); }
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
    virtual BaseVector & Scale (Complex scal) override;
    virtual BaseVector & Scale (BaseScalar & scal) override;
    virtual BaseVector & SetScalar (double scal) override;
    virtual BaseVector & SetScalar (Complex scal) override;
    virtual BaseVector & Set (double scal, const BaseVector & v) override;
    virtual BaseVector & Set (Complex scal, const BaseVector & v) override;
    virtual BaseVector & Add (double scal, const BaseVector & v) override;
    virtual BaseVector & Add (Complex scal, const BaseVector & v) override;
    virtual BaseVector & Add (BaseScalar & scal, const BaseVector & v) override;

    // using BaseVector::InnerProduct;
    using S_BaseVector<T>::InnerProduct;
    virtual double InnerProductD (const BaseVector & v2) const override;
    virtual Complex InnerProductC (const BaseVector & v2, bool conjugate = false) const override;
    virtual void InnerProduct (const BaseVector & v2, BaseScalar & scal, bool conjugate = false) const override;
    virtual double L2Norm () const override;

    virtual void * Memory () const override;
    virtual FlatVector<T> FVScal () const override;
    virtual AutoVector CreateVector () const override;
    virtual AutoVector Range (T_Range<size_t> range) const override;
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
    bool subrange = false;           // wraps a sub-range of vec
    bool initial_host_uptodate = false;
    bool converting = false;         // hostmem holds a converted copy of vec
    size_t convert_offset = 0;       // into vec, for the converted sub-range
  public:
    DeviceVectorWrapper (const BaseVector & avec, MemType amemtype = MemType::Shared,
                         optional<IntRange> opt_range = nullopt);
    ~DeviceVectorWrapper();

    using DeviceVector<T>::operator=;
  };


#if !defined(FILE_DEVICEVECTOR_CPP)
  extern template class DeviceScalar<double>;
  extern template class DeviceScalar<float>;
  extern template class DeviceVector<double>;
  extern template class DeviceVector<float>;
  extern template class DeviceVector<Complex>;
  extern template class DeviceVectorWrapper<double>;
  extern template class DeviceVectorWrapper<float>;
  extern template class DeviceVectorWrapper<Complex>;
#endif
}

#endif
