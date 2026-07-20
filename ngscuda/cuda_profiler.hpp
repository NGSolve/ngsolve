#ifndef NGS_CUDA_PROFILER_HPP
#define NGS_CUDA_PROFILER_HPP

#include <algorithm>
#include <core/profiler.hpp>
#include <initializer_list>
#include "cuda_ngstd.hpp"

namespace ngs_cuda
{
  using ngcore::Array;
  extern int gpu_clock;
// #define NGS_CUDA_DEVICE_TIMERS

#ifdef __CUDACC__
#ifdef NGS_CUDA_DEVICE_TIMERS
  constexpr int N_MAX_DEVICE_TIMERS = 16;
  struct DevTraceData {
    unsigned long long start;   // clock64 at region start (0 = entry unused)
    unsigned long long stop;    // clock64 at region stop
    union UserData {
      double d;
      float f[2];
      unsigned u[2];
    } data;
    unsigned blockNr;
    unsigned short timer_nr;
  };

  struct DevTimerData {
    unsigned long long time[N_MAX_DEVICE_TIMERS];
  };

  struct DevTraceBlockData {
    unsigned long long start;         // clock64 at block start (per-SM, high-res)
    unsigned long long stop;          // clock64 at block stop
    unsigned long long start_global;  // %globaltimer at block start (global anchor)
    unsigned long long stop_global;   // %globaltimer at block stop
    unsigned start_smid;
    unsigned stop_smid;
  };

  constexpr unsigned N_MAX_TRACER_OBJECTS = 4*1024*1024;
  constexpr int N_MAX_BLOCKS = 65*1024;

  // each thread block reserves trace entries in chunks of this size from one
  // global atomic counter; entry N_MAX_TRACER_OBJECTS of d_trace_data serves
  // as overflow sink once the capacity is exhausted
  constexpr unsigned TRACE_CHUNK_SIZE = 1024;

  struct DevTraceState {
    unsigned chunk_counter;   // next free chunk of d_trace_data
    unsigned nblocks;         // total number of blocks, written by block 0
  };

  extern __device__ DevTimerData d_timer_data[];
  extern __device__ DevTraceData d_trace_data[];
  extern __device__ DevTraceBlockData d_block_data[];
  extern __device__ DevTraceState d_trace_state;
  extern Array<DevTimerData> timer_data;
  extern Array<DevTraceData> trace_data;
  extern Array<DevTraceBlockData> block_data;


  __device__ inline unsigned GetSMID() {
    unsigned ret;
    asm("mov.u32 %0, %smid;" : "=r"(ret) );
    return ret;
  }

  // %globaltimer is a nanosecond timer that is synchronized across all SMs
  // (unlike %clock64, which is a per-SM cycle counter). volatile prevents the
  // compiler from reordering or merging the start/stop reads.
  __device__ inline unsigned long long GetGlobalTimer() {
    unsigned long long ret;
    asm volatile("mov.u64 %0, %%globaltimer;" : "=l"(ret) );
    return ret;
  }

  struct DeviceBlockRegionTracer {
    int blockNr, threadNr;
    unsigned trace_pos, trace_end;   // this block's current chunk of trace entries

    __device__ DeviceBlockRegionTracer(dim3 gridDim, dim3 blockDim, dim3 blockIdx, dim3 threadIdx)
      : blockNr(blockIdx.x+gridDim.x*blockIdx.y), threadNr(threadIdx.x+blockDim.x*threadIdx.y),
        trace_pos(0), trace_end(0)   // empty chunk -> first AllocTraceEntry reserves one
    {
      if(threadNr == 0)
      {
        // read both clocks back-to-back so the pair marks the same instant
        d_block_data[blockNr].start = clock64();
        d_block_data[blockNr].start_global = GetGlobalTimer();
        d_block_data[blockNr].start_smid = GetSMID();
        if(blockNr == 0)
          d_trace_state.nblocks = gridDim.x*gridDim.y;
      }
      __syncwarp();
    }

    // reserve the next trace entry in this block's current chunk, grabbing a
    // new chunk from the global counter when it is used up; only called from
    // the tracing thread (threadNr == 0), so no synchronization on the local
    // counters is needed
    __device__ unsigned AllocTraceEntry()
    {
      if(trace_pos == trace_end)
      {
        trace_pos = TRACE_CHUNK_SIZE * atomicAdd(&d_trace_state.chunk_counter, 1);
        trace_end = trace_pos + TRACE_CHUNK_SIZE;
      }
      if(trace_pos >= N_MAX_TRACER_OBJECTS)
      {
        trace_pos = N_MAX_TRACER_OBJECTS;       // park on the overflow sink and
        trace_end = N_MAX_TRACER_OBJECTS + 1;   // don't reserve further chunks
        return N_MAX_TRACER_OBJECTS;
      }
      return trace_pos++;
    }

    __device__ ~DeviceBlockRegionTracer()
    {
      if(threadNr == 0)
      {
        d_block_data[blockNr].stop = clock64();
        d_block_data[blockNr].stop_global = GetGlobalTimer();
        d_block_data[blockNr].stop_smid = GetSMID();
      }
    }
  };

  struct DeviceRegionTimer
  {
    int timer_nr, id;

    __device__ DeviceRegionTimer( int timer_nr_, int id_=-1)
      : timer_nr(timer_nr_), id(id_ == -1 ? blockIdx.x : id_)
    {
      if(threadIdx.x ==0)
        atomicAdd(&d_timer_data[id].time[timer_nr], -clock64());
    }

    __device__ ~DeviceRegionTimer()
    {
      if(threadIdx.x ==0)
        atomicAdd(&d_timer_data[id].time[timer_nr], clock64());
    }

  };

  struct DeviceRegionTracer
  {
    bool active;
    unsigned id;
    DevTraceData::UserData & data;

    __device__ DeviceRegionTracer( DeviceBlockRegionTracer & tr, int timer_nr )
      : active(tr.threadNr == 0),
        id(active ? tr.AllocTraceEntry() : N_MAX_TRACER_OBJECTS),
        data(d_trace_data[id].data)
    {
      if(active)
      {
        d_trace_data[id].blockNr = tr.blockNr;
        d_trace_data[id].timer_nr = timer_nr;
        d_trace_data[id].start = clock64();   // written last: marks the entry as used
      }
    }

    __device__ ~DeviceRegionTracer()
    {
      if(active)
        d_trace_data[id].stop = clock64();
    }

  };
#else
  struct DeviceBlockRegionTracer
  {
    __device__ DeviceBlockRegionTracer( dim3, dim3, dim3, dim3 ) {}
  };
  struct DeviceRegionTimer
  {
    __device__ DeviceRegionTimer( int , int ) {}
  };
  struct DeviceRegionTracer
  {
    __device__ DeviceRegionTracer( DeviceBlockRegionTracer & tr, int timer_nr ) { }
  };
#endif // NGS_CUDA_DEVICE_TIMERS

#endif // __CUDACC__


  namespace detail {
  struct GPUTimerUserData
  {
    int timernr;
    bool is_start;
  };

  // encode UserData directly in void* pointer (8 bytes available)
  inline void * EncodeCallbackData( bool start, int nr )
  {
    static_assert(sizeof(GPUTimerUserData) <= sizeof(void*), "cannot encode UserData in void*");
    GPUTimerUserData data{nr, start};
    return *reinterpret_cast<void**>(&data);
  }

  inline void CUDART_CB Callback(cudaStream_t stream, cudaError_t status, void *userData)
  {
    const GPUTimerUserData & ud = *reinterpret_cast<GPUTimerUserData*>(&userData);
    if(!ngcore::trace) return;
    if(ud.is_start)
      ngcore::trace->StartGPU(ud.timernr);
    else
      ngcore::trace->StopGPU(ud.timernr);
  }

  inline void CUDART_CB CallbackKernelStart(cudaStream_t stream, cudaError_t status, void *userData)
  {
    ngcore::TTimePoint & t = *reinterpret_cast<ngcore::TTimePoint*>(userData);
    t = ngcore::GetTimeCounter();
  }

  inline void StartGPUTimer(int nr) {
    cudaStreamAddCallback(0, Callback, EncodeCallbackData(true, nr), 0);
  }

  inline void StopGPUTimer(int nr) {
    cudaStreamAddCallback(0, Callback, EncodeCallbackData(false, nr), 0);
  }
  } // namspace detail

  class CudaRegionTimer
  {
    const ngcore::Timer<> & timer;
    bool is_already_stopped = false;
    ngcore::TTimePoint t_kernel_start;

    static bool is_cuda_timer_enabled;
    const std::initializer_list<const char *> * timer_names;

  public:
    CudaRegionTimer (ngcore::Timer<> & atimer, std::initializer_list<const char*> * atimer_names = nullptr) : timer(atimer), timer_names(atimer_names) {
      timer.Start();
      if(is_cuda_timer_enabled)
        detail::StartGPUTimer(timer);
#ifdef NGS_CUDA_DEVICE_TIMERS
      cudaStreamAddCallback(0, detail::CallbackKernelStart, &t_kernel_start, 0);
#endif // NGS_CUDA_DEVICE_TIMERS
    }
    ~CudaRegionTimer ()
    {
      if(!is_already_stopped)
        Stop();
    }
    CudaRegionTimer() = delete;
    CudaRegionTimer(const CudaRegionTimer &) = delete;
    CudaRegionTimer(CudaRegionTimer &&) = delete;
    void operator=(const CudaRegionTimer &) = delete;
    void operator=(CudaRegionTimer &&) = delete;

    void Stop() {
      if(is_cuda_timer_enabled)
        detail::StopGPUTimer(timer);
#ifdef NGS_CUDA_DEVICE_TIMERS
      ProcessTracingData(timer_names);
#endif // NGS_CUDA_DEVICE_TIMERS
      timer.Stop();
      is_already_stopped = true;
    }

    void ProcessTracingData(const std::initializer_list<const char *> * timer_names);

    ngcore::TTimePoint TimeD2H( long long ns )
    {
      double t_seconds = 1e-9*ns + 1e-6; // %globaltimer is in ns; add 1us for kernel startup
      return t_kernel_start + t_seconds / ngcore::seconds_per_tick;
    }

    static void SetCudaTimer(bool enabled){
      is_cuda_timer_enabled = enabled;
    }
  };

  void TimeProfiler();

#ifdef NGS_CUDA_DEVICE_TIMERS
  // device addresses of the tracing buffers, for jit-compiled kernels
  // (separate cuda modules cannot link against the __device__ symbols)
  void * GetDevTraceDataPtr();
  void * GetDevTraceBlockDataPtr();
  void * GetDevTraceStatePtr();
#endif
}

#endif // NGS_CUDA_PROFILER_HPP
