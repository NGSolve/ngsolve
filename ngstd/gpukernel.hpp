#ifndef FILE_GPUKERNEL
#define FILE_GPUKERNEL

/*********************************************************************/
/* File:   gpukernel.hpp                                             */
/* Date:   29. Aug. 2025                                             */
/*********************************************************************/

/*
  Common kernel syntax for cuda and metal.

  Prepend code_gpukernel to the source given to Device::CompileSource,
  the same way code_tinybla is prepended. A kernel is then written once:

      KERNEL(saxpy, GLOBAL_IN(float,x), GLOBAL(float,y),
                    VALUE(float,a), VALUE(int,n))
      {
        int i = GLOBAL_ID_X;
        if (i < n) y[i] = a*x[i] + y[i];
      }

  Arguments are positional: the i-th declared argument is the i-th entry
  of the args list passed to Queue::Launch. Metal assigns buffer indices
  in declaration order, cuda takes them as kernel parameters.
*/

namespace ngs_gpu
{
  inline const char * code_gpukernel = R"RAW(

#ifdef __CUDACC__

  #define NGS_INDICES
  #define KERNEL(name, ...)   extern "C" __global__ void name(__VA_ARGS__ NGS_INDICES)

  #define GLOBAL(T,name)      T * name
  #define GLOBAL_IN(T,name)   const T * name
  #define VALUE(T,name)       T name

  #define GLOBAL_PTR(T)       T *
  #define LOCAL_PTR(T)        T *
  #define GPU_FUNC            __device__

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

  #define SHARED(T,name,N)  __shared__ T name[N]
  #define BARRIER()         __syncthreads()

#elif defined(NGS_GPU_CPU)

  /* host reference backend: the kernel is ordinary c++, work-items of one
     group run as threads so that BARRIER and SHARED behave like on a gpu */

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

  #include <cstddef>
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
      std::vector<char> arena;
      std::size_t used = 0;
      void * last = nullptr;
      Group (unsigned n, std::size_t sh) : barrier(n), arena(sh) { }
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

  extern "C" void _ngs_launch (void (*thunk)(void**), void ** args,
                               unsigned gx, unsigned gy, unsigned gz,
                               unsigned sx, unsigned sy, unsigned sz,
                               std::size_t dynshared)
  {
    using namespace ngs_cpu;
    unsigned nthreads = sx*sy*sz;
    Group g (nthreads, 65536 + dynshared);

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
    extern "C" void name##_ngsthunk (void ** _a)                        \
    { ngs_cpu::Call (&name, _a); }                                      \
    static void name(__VA_ARGS__)

  #define GLOBAL(T,name)      T * name
  #define GLOBAL_IN(T,name)   const T * name
  #define VALUE(T,name)       T name

  #define GLOBAL_PTR(T)       T *
  #define LOCAL_PTR(T)        T *
  #define GPU_FUNC

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

  #define SHARED(T,name,N)  T * name = (T*) ngs_cpu::SharedAlloc(sizeof(T)*(N))
  #define BARRIER()         ngs_cpu::group->barrier.Wait()

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

  #define SHARED(T,name,N)  threadgroup T name[N]
  #define BARRIER()         threadgroup_barrier(mem_flags::mem_threadgroup)

#endif

)RAW";
}

#endif
