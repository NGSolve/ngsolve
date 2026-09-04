/*********************************************************************/
/* File:   gpukernel.h                                               */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   4. Sep. 2026                                              */
/*********************************************************************/

#ifdef __CUDACC__

  #define NGS_INDICES
  #define KERNEL(name, ...)   extern "C" __global__ void name(__VA_ARGS__ NGS_INDICES)
  // group size and resident groups per SM the kernel is compiled for
  #define KERNEL_BOUNDS(name, maxthreads, mingroups, ...) \
    extern "C" __global__ void __launch_bounds__(maxthreads, mingroups) name(__VA_ARGS__ NGS_INDICES)

  #define GLOBAL(T,name)      T * name
  #define GLOBAL_IN(T,name)   const T * name
  #define VALUE(T,name)       T name

  #define GLOBAL_PTR(T)       T *
  #define LOCAL_PTR(T)        T *
  #define GPU_FUNC            __device__

  typedef unsigned int   uint;
  typedef unsigned short ushort;

  #define GLOBAL_ID_X   (blockIdx.x*blockDim.x + threadIdx.x)
  #define GLOBAL_ID_Y   (blockIdx.y*blockDim.y + threadIdx.y)
  #define GLOBAL_ID_Z   (blockIdx.z*blockDim.z + threadIdx.z)
  #define LOCAL_ID_X    threadIdx.x
  #define LOCAL_ID_Y    threadIdx.y
  #define LOCAL_ID_Z    threadIdx.z
  #define GROUP_ID_X    blockIdx.x
  #define GROUP_ID_Y    blockIdx.y
  #define GROUP_ID_Z    blockIdx.z
  #define GROUP_SIZE_X  blockDim.x
  #define GROUP_SIZE_Y  blockDim.y
  #define GROUP_SIZE_Z  blockDim.z
  #define NUM_GROUPS_X  gridDim.x
  #define NUM_GROUPS_Y  gridDim.y
  #define NUM_GROUPS_Z  gridDim.z

  // program-scope table, may be indexed by a runtime value
  #define CONSTANT_ARRAY(T,name,N)  __device__ const T name[N]

  #define SHARED(T,name,N)       alignas(16) __shared__ T name[N]
  #define SHARED_2D(T,name,N,M)  alignas(16) __shared__ T name[N][M]
  #define BARRIER()              __syncthreads()

  #define GLOBAL_ATOMIC(T,name)  T * name
  #define ATOMIC_PTR(T)          T *
  #define ATOMIC_ADD(p,val)      atomicAdd(p, val)
  template <typename T> __device__ inline T ngs_atomic_load (T * p)
  { return *(volatile T*)p; }
  #define ATOMIC_LOAD(p)         ngs_atomic_load(p)

  // synchronisation below the group: the lanes of one warp / simdgroup
  // (the x-lanes of a group whose GROUP_SIZE_X is the simd width), and
  // a fence making this thread's global writes visible device-wide
  #define SIMD_BARRIER()         __syncwarp()
  #define DEVICE_FENCE()         __threadfence()

  #define NGS_CPLX_FUNC          __host__ __device__ inline
  #define NGS_CPLX_SQRT          sqrt
  #define NGS_CPLX_ADDRSPACES(DEF)  DEF()
  // a complex atomic buffer is a plain pointer, added component-wise
  #define GLOBAL_ATOMIC_COMPLEX(T,name)   Complex<T> * name
  #define ATOMIC_ADD_COMPLEX(name,i,val)  atomicAdd(&name[i], val)

#elif defined(NGS_GPU_CPU)

  /* host reference backend: the kernel is ordinary c++, work-items of one
     group run as threads so that BARRIER and SHARED behave like on a gpu */

  /* GetProcAddress only sees exported symbols */
  #ifdef _WIN32
    #define NGS_KERNEL_EXPORT __declspec(dllexport)
  #else
    #define NGS_KERNEL_EXPORT
  #endif

  /* tinybla defines "thread" as an empty address-space macro, which would
     wipe out class thread in <thread>. Hide it across the includes, and keep
     an alias so the code below never spells std::thread. Works whichever way
     round the two preludes are concatenated. */
  #pragma push_macro("thread")
  #pragma push_macro("constant")
  #pragma push_macro("threadgroup")
  #undef thread
  #undef constant
  #undef threadgroup

  #include <atomic>
  #include <cmath>
  #include <cstddef>
  #include <deque>
  #include <thread>
  #include <vector>
  #include <mutex>
  #include <condition_variable>
  #include <utility>
  #include <type_traits>

  namespace ngs_cpu { using cpu_thread = ::std::thread; }

  #pragma pop_macro("threadgroup")
  #pragma pop_macro("constant")
  #pragma pop_macro("thread")

  namespace ngs_cpu
  {
    class Barrier
    {
      unsigned count, waiting = 0, generation = 0;
      std::mutex m;
      std::condition_variable cv;
    public:
      Barrier (unsigned n) : count(n) { }
      void Wait()
      {
        std::unique_lock<std::mutex> lk(m);
        unsigned g = generation;
        if (++waiting == count) { waiting = 0; generation++; cv.notify_all(); }
        else cv.wait (lk, [&]{ return generation != g; });
      }
    };

    struct Group
    {
      Barrier barrier;
      std::deque<Barrier> rowbarriers;   // one per (y,z) row of x-lanes
      std::vector<char> arena;
      std::size_t used = 0;
      void * last = nullptr;
      Group (unsigned sx, unsigned sy, unsigned sz, std::size_t sh)
        : barrier(sx*sy*sz), arena(sh)
      {
        for (unsigned i = 0; i < sy*sz; i++) rowbarriers.emplace_back (sx);
      }
    };

    thread_local unsigned lid[3] = {0,0,0}, gid[3] = {0,0,0};
    thread_local unsigned gsz[3] = {1,1,1}, ngr[3] = {1,1,1};
    thread_local Group * group = nullptr;

    // all work-items reach this in the same order, so one of them bumps
    inline void * SharedAlloc (std::size_t bytes)
    {
      Group & g = *group;
      g.barrier.Wait();
      if (lid[0]==0 && lid[1]==0 && lid[2]==0)
        {
          std::size_t off = (g.used + 15) & ~std::size_t(15);
          g.last = g.arena.data() + off;
          g.used = off + bytes;
        }
      g.barrier.Wait();
      return g.last;
    }

    // plain float has no fetch_add, atomic_ref gives the gpu's atomicAdd
    template <typename T> inline T AtomicAdd (T * p, T val)
    { return std::atomic_ref<T>(*p).fetch_add (val, std::memory_order_relaxed); }
    // for spin-waits: let the thread being waited for run
    template <typename T> inline T AtomicLoad (T * p)
    { std::this_thread::yield(); return std::atomic_ref<T>(*p).load (std::memory_order_relaxed); }

    // unpack void** using the static kernel signature - no libffi needed
    template <typename T> inline T Get (void * p)
    {
      if constexpr (std::is_pointer_v<T>) return (T)p;
      else return *(T*)p;
    }
    template <typename... A, std::size_t... I>
    inline void CallImpl (void (*f)(A...), void ** a, std::index_sequence<I...>)
    { f (Get<A>(a[I])...); }
    template <typename... A>
    inline void Call (void (*f)(A...), void ** a)
    { CallImpl (f, a, std::index_sequence_for<A...>{}); }
  }

  extern "C" NGS_KERNEL_EXPORT void _ngs_launch (void (*thunk)(void**), void ** args,
                               unsigned gx, unsigned gy, unsigned gz,
                               unsigned sx, unsigned sy, unsigned sz,
                               std::size_t dynshared)
  {
    using namespace ngs_cpu;
    unsigned nthreads = sx*sy*sz;
    Group g (sx, sy, sz, 65536 + dynshared);

    auto worker = [&] (unsigned t)
    {
      lid[0] = t % sx; lid[1] = (t/sx) % sy; lid[2] = t/(sx*sy);
      gsz[0] = sx; gsz[1] = sy; gsz[2] = sz;
      ngr[0] = gx; ngr[1] = gy; ngr[2] = gz;
      group = &g;
      for (unsigned bz = 0; bz < gz; bz++)
        for (unsigned by = 0; by < gy; by++)
          for (unsigned bx = 0; bx < gx; bx++)
            {
              gid[0] = bx; gid[1] = by; gid[2] = bz;
              if (t == 0) g.used = 0;
              g.barrier.Wait();
              thunk (args);
              g.barrier.Wait();
            }
    };

    std::vector<ngs_cpu::cpu_thread> threads;
    for (unsigned t = 1; t < nthreads; t++) threads.emplace_back (worker, t);
    worker (0);
    for (auto & th : threads) th.join();
  }

  #define KERNEL(name, ...)                                             \
    static void name(__VA_ARGS__);                                      \
    extern "C" NGS_KERNEL_EXPORT void name##_ngsthunk (void ** _a)      \
    { ngs_cpu::Call (&name, _a); }                                      \
    static void name(__VA_ARGS__)
  #define KERNEL_BOUNDS(name, maxthreads, mingroups, ...) KERNEL(name, __VA_ARGS__)

  #define GLOBAL(T,name)      T * name
  #define GLOBAL_IN(T,name)   const T * name
  #define VALUE(T,name)       T name

  #define GLOBAL_PTR(T)       T *
  #define LOCAL_PTR(T)        T *
  #define GPU_FUNC

  typedef unsigned int   uint;
  typedef unsigned short ushort;

  #define GLOBAL_ID_X   (ngs_cpu::gid[0]*ngs_cpu::gsz[0] + ngs_cpu::lid[0])
  #define GLOBAL_ID_Y   (ngs_cpu::gid[1]*ngs_cpu::gsz[1] + ngs_cpu::lid[1])
  #define GLOBAL_ID_Z   (ngs_cpu::gid[2]*ngs_cpu::gsz[2] + ngs_cpu::lid[2])
  #define LOCAL_ID_X    ngs_cpu::lid[0]
  #define LOCAL_ID_Y    ngs_cpu::lid[1]
  #define LOCAL_ID_Z    ngs_cpu::lid[2]
  #define GROUP_ID_X    ngs_cpu::gid[0]
  #define GROUP_ID_Y    ngs_cpu::gid[1]
  #define GROUP_ID_Z    ngs_cpu::gid[2]
  #define GROUP_SIZE_X  ngs_cpu::gsz[0]
  #define GROUP_SIZE_Y  ngs_cpu::gsz[1]
  #define GROUP_SIZE_Z  ngs_cpu::gsz[2]
  #define NUM_GROUPS_X  ngs_cpu::ngr[0]
  #define NUM_GROUPS_Y  ngs_cpu::ngr[1]
  #define NUM_GROUPS_Z  ngs_cpu::ngr[2]

  #define CONSTANT_ARRAY(T,name,N)  static const T name[N]

  #define SHARED(T,name,N)       T * name = (T*) ngs_cpu::SharedAlloc(sizeof(T)*(N))
  #define SHARED_2D(T,name,N,M)  T (*name)[M] = (T(*)[M]) ngs_cpu::SharedAlloc(sizeof(T)*(N)*(M))
  #define BARRIER()              ngs_cpu::group->barrier.Wait()

  #define GLOBAL_ATOMIC(T,name)  T * name
  #define ATOMIC_PTR(T)          T *
  #define ATOMIC_ADD(p,val)      ngs_cpu::AtomicAdd(p, val)
  #define ATOMIC_LOAD(p)         ngs_cpu::AtomicLoad(p)

  #define SIMD_BARRIER()         ngs_cpu::group->rowbarriers[ngs_cpu::lid[1] + ngs_cpu::gsz[1]*ngs_cpu::lid[2]].Wait()
  #define DEVICE_FENCE()         std::atomic_thread_fence(std::memory_order_seq_cst)

  #define NGS_CPLX_FUNC          inline
  #define NGS_CPLX_SQRT          std::sqrt
  #define NGS_CPLX_ADDRSPACES(DEF)  DEF()
  #define GLOBAL_ATOMIC_COMPLEX(T,name)   Complex<T> * name
  #define ATOMIC_ADD_COMPLEX(name,i,val)  ngs_cpu::AtomicAdd(&name[i], val)

#else   // metal

  #include <metal_stdlib>
  using namespace metal;

  #define NGS_INDICES \
    , uint3 _ngs_gtid  [[thread_position_in_grid]] \
    , uint3 _ngs_ltid  [[thread_position_in_threadgroup]] \
    , uint3 _ngs_grpid [[threadgroup_position_in_grid]] \
    , uint3 _ngs_grpsz [[threads_per_threadgroup]] \
    , uint3 _ngs_ngrps [[threadgroups_per_grid]]

  #define KERNEL(name, ...)   kernel void name(__VA_ARGS__ NGS_INDICES)
  #define KERNEL_BOUNDS(name, maxthreads, mingroups, ...) \
    [[max_total_threads_per_threadgroup(maxthreads)]] kernel void name(__VA_ARGS__ NGS_INDICES)

  #define GLOBAL(T,name)      device T * name
  #define GLOBAL_IN(T,name)   device const T * name
  #define VALUE(T,name)       constant T & name

  #define GLOBAL_PTR(T)       device T *
  #define LOCAL_PTR(T)        threadgroup T *
  #define GPU_FUNC

  #define GLOBAL_ID_X   _ngs_gtid.x
  #define GLOBAL_ID_Y   _ngs_gtid.y
  #define GLOBAL_ID_Z   _ngs_gtid.z
  #define LOCAL_ID_X    _ngs_ltid.x
  #define LOCAL_ID_Y    _ngs_ltid.y
  #define LOCAL_ID_Z    _ngs_ltid.z
  #define GROUP_ID_X    _ngs_grpid.x
  #define GROUP_ID_Y    _ngs_grpid.y
  #define GROUP_ID_Z    _ngs_grpid.z
  #define GROUP_SIZE_X  _ngs_grpsz.x
  #define GROUP_SIZE_Y  _ngs_grpsz.y
  #define GROUP_SIZE_Z  _ngs_grpsz.z
  #define NUM_GROUPS_X  _ngs_ngrps.x
  #define NUM_GROUPS_Y  _ngs_ngrps.y
  #define NUM_GROUPS_Z  _ngs_ngrps.z

  #define CONSTANT_ARRAY(T,name,N)  constant T name[N]

  #define SHARED(T,name,N)       alignas(16) threadgroup T name[N]
  #define SHARED_2D(T,name,N,M)  alignas(16) threadgroup T name[N][M]
  #define BARRIER()              threadgroup_barrier(mem_flags::mem_threadgroup)

  #define GLOBAL_ATOMIC(T,name)  device metal::atomic<T> * name
  #define ATOMIC_PTR(T)          device metal::atomic<T> *
  #define ATOMIC_ADD(p,val)      metal::atomic_fetch_add_explicit(p, val, metal::memory_order_relaxed)
  #define ATOMIC_LOAD(p)         metal::atomic_load_explicit(p, metal::memory_order_relaxed)

  #define SIMD_BARRIER()         simdgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup)
  #define DEVICE_FENCE()         atomic_thread_fence(mem_flags::mem_device, metal::memory_order_seq_cst, metal::thread_scope_device)

  #define NGS_CPLX_FUNC          inline
  #define NGS_CPLX_SQRT          sqrt
  // references to a struct carry an address space, one overload each
  #define NGS_CPLX_ADDRSPACES(DEF)  DEF(thread) DEF(device) DEF(threadgroup)
  // metal has no atomic<Complex>: the buffer is viewed as 2n reals
  #define GLOBAL_ATOMIC_COMPLEX(T,name)   device metal::atomic<T> * name
  #define ATOMIC_ADD_COMPLEX(name,i,val)  \
    { ATOMIC_ADD(&name[2*(i)], (val).re); ATOMIC_ADD(&name[2*(i)+1], (val).im); }

#endif


/*
  Complex numbers, the same struct on every backend so results agree
  bit for bit. Layout re,im as std::complex, so buffers transfer as is.
  Use Complex<float> / Complex<double> explicitly. Mixed arithmetic with
  a real converts the real to the value_type; write T(2) rather than 2
  where the compiler cannot deduce it.

    Complex<float> z(1, 2), w = z*conj(z) + T(0.5);   Norm(z) == |z|^2

  A buffer accumulated with atomics is declared GLOBAL_ATOMIC_COMPLEX(T,y)
  and written with ATOMIC_ADD_COMPLEX(y, i, val).
*/

template <typename T>
struct Complex
{
  typedef T value_type;
  T re, im;
  Complex () = default;
  NGS_CPLX_FUNC Complex (T are) : re(are), im(0) { }
  NGS_CPLX_FUNC Complex (T are, T aim) : re(are), im(aim) { }
};

// the scalar operand in a non-deduced position, so ints and doubles convert
#define NGS_CPLX_SCAL(T) typename Complex<T>::value_type

template <typename T> NGS_CPLX_FUNC T real (Complex<T> z) { return z.re; }
template <typename T> NGS_CPLX_FUNC T imag (Complex<T> z) { return z.im; }
template <typename T> NGS_CPLX_FUNC Complex<T> conj (Complex<T> z) { return Complex<T>(z.re, -z.im); }
template <typename T> NGS_CPLX_FUNC T Norm (Complex<T> z) { return z.re*z.re + z.im*z.im; }
template <typename T> NGS_CPLX_FUNC T abs (Complex<T> z) { return NGS_CPLX_SQRT(z.re*z.re + z.im*z.im); }

template <typename T> NGS_CPLX_FUNC Complex<T> operator- (Complex<T> a)
{ return Complex<T>(-a.re, -a.im); }

template <typename T> NGS_CPLX_FUNC Complex<T> operator+ (Complex<T> a, Complex<T> b)
{ return Complex<T>(a.re+b.re, a.im+b.im); }
template <typename T> NGS_CPLX_FUNC Complex<T> operator+ (Complex<T> a, NGS_CPLX_SCAL(T) b)
{ return Complex<T>(a.re+b, a.im); }
template <typename T> NGS_CPLX_FUNC Complex<T> operator+ (NGS_CPLX_SCAL(T) a, Complex<T> b)
{ return Complex<T>(a+b.re, b.im); }

template <typename T> NGS_CPLX_FUNC Complex<T> operator- (Complex<T> a, Complex<T> b)
{ return Complex<T>(a.re-b.re, a.im-b.im); }
template <typename T> NGS_CPLX_FUNC Complex<T> operator- (Complex<T> a, NGS_CPLX_SCAL(T) b)
{ return Complex<T>(a.re-b, a.im); }
template <typename T> NGS_CPLX_FUNC Complex<T> operator- (NGS_CPLX_SCAL(T) a, Complex<T> b)
{ return Complex<T>(a-b.re, -b.im); }

template <typename T> NGS_CPLX_FUNC Complex<T> operator* (Complex<T> a, Complex<T> b)
{ return Complex<T>(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re); }
template <typename T> NGS_CPLX_FUNC Complex<T> operator* (Complex<T> a, NGS_CPLX_SCAL(T) b)
{ return Complex<T>(a.re*b, a.im*b); }
template <typename T> NGS_CPLX_FUNC Complex<T> operator* (NGS_CPLX_SCAL(T) a, Complex<T> b)
{ return Complex<T>(a*b.re, a*b.im); }

// scaled to avoid overflow in the denominator
template <typename T> NGS_CPLX_FUNC Complex<T> operator/ (Complex<T> a, Complex<T> b)
{
  if ((b.re < 0 ? -b.re : b.re) >= (b.im < 0 ? -b.im : b.im))
    {
      T r = b.im / b.re, den = b.re + b.im*r;
      return Complex<T>((a.re + a.im*r) / den, (a.im - a.re*r) / den);
    }
  T r = b.re / b.im, den = b.re*r + b.im;
  return Complex<T>((a.re*r + a.im) / den, (a.im*r - a.re) / den);
}
template <typename T> NGS_CPLX_FUNC Complex<T> operator/ (Complex<T> a, NGS_CPLX_SCAL(T) b)
{ return Complex<T>(a.re/b, a.im/b); }
template <typename T> NGS_CPLX_FUNC Complex<T> operator/ (NGS_CPLX_SCAL(T) a, Complex<T> b)
{ return Complex<T>(a) / b; }

template <typename T> NGS_CPLX_FUNC bool operator== (Complex<T> a, Complex<T> b)
{ return a.re == b.re && a.im == b.im; }
template <typename T> NGS_CPLX_FUNC bool operator!= (Complex<T> a, Complex<T> b)
{ return !(a == b); }

// tinybla accumulates through fma(a,b,c)
template <typename T> NGS_CPLX_FUNC Complex<T> fma (Complex<T> a, Complex<T> b, Complex<T> c)
{ return a*b + c; }
template <typename T> NGS_CPLX_FUNC Complex<T> fma (NGS_CPLX_SCAL(T) a, Complex<T> b, Complex<T> c)
{ return a*b + c; }
template <typename T> NGS_CPLX_FUNC Complex<T> fma (Complex<T> a, NGS_CPLX_SCAL(T) b, Complex<T> c)
{ return a*b + c; }

#define NGS_CPLX_COMPOUND(AS)                                                        \
  template <typename T> NGS_CPLX_FUNC AS Complex<T> & operator+= (AS Complex<T> & a, Complex<T> b) \
  { a.re += b.re; a.im += b.im; return a; }                                          \
  template <typename T> NGS_CPLX_FUNC AS Complex<T> & operator-= (AS Complex<T> & a, Complex<T> b) \
  { a.re -= b.re; a.im -= b.im; return a; }                                          \
  template <typename T> NGS_CPLX_FUNC AS Complex<T> & operator*= (AS Complex<T> & a, Complex<T> b) \
  { Complex<T> t = Complex<T>(a.re, a.im) * b; a.re = t.re; a.im = t.im; return a; } \
  template <typename T> NGS_CPLX_FUNC AS Complex<T> & operator/= (AS Complex<T> & a, Complex<T> b) \
  { Complex<T> t = Complex<T>(a.re, a.im) / b; a.re = t.re; a.im = t.im; return a; } \
  template <typename T> NGS_CPLX_FUNC AS Complex<T> & operator+= (AS Complex<T> & a, NGS_CPLX_SCAL(T) b) \
  { a.re += b; return a; }                                                           \
  template <typename T> NGS_CPLX_FUNC AS Complex<T> & operator-= (AS Complex<T> & a, NGS_CPLX_SCAL(T) b) \
  { a.re -= b; return a; }                                                           \
  template <typename T> NGS_CPLX_FUNC AS Complex<T> & operator*= (AS Complex<T> & a, NGS_CPLX_SCAL(T) b) \
  { a.re *= b; a.im *= b; return a; }                                                \
  template <typename T> NGS_CPLX_FUNC AS Complex<T> & operator/= (AS Complex<T> & a, NGS_CPLX_SCAL(T) b) \
  { a.re /= b; a.im /= b; return a; }

NGS_CPLX_ADDRSPACES(NGS_CPLX_COMPOUND)

#ifdef __CUDACC__
template <typename T> __device__ inline void atomicAdd (Complex<T> * p, Complex<T> v)
{ atomicAdd (&p->re, v.re); atomicAdd (&p->im, v.im); }
#elif defined(NGS_GPU_CPU)
namespace ngs_cpu
{
  template <typename T> inline void AtomicAdd (Complex<T> * p, Complex<T> v)
  { AtomicAdd (&p->re, v.re); AtomicAdd (&p->im, v.im); }
}
#endif

