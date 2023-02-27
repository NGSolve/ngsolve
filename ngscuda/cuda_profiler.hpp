#ifndef NGS_CUDA_PROFILER_HPP
#define NGS_CUDA_PROFILER_HPP

#include <algorithm>
#include <core/profiler.hpp>
#include "cuda_ngstd.hpp"

namespace ngs_cuda
{
  using ngcore::Array;
  extern int gpu_clock;

#ifdef __CUDACC__
  // #define NGS_CUDA_DEVICE_TIMERS
#ifdef NGS_CUDA_DEVICE_TIMERS
  constexpr int N_MAX_DEVICE_TIMERS = 16;
  struct DevTraceData {
    unsigned long long start[N_MAX_DEVICE_TIMERS];
    unsigned long long stop[N_MAX_DEVICE_TIMERS];
    unsigned short start_smid[N_MAX_DEVICE_TIMERS];
    unsigned short stop_smid[N_MAX_DEVICE_TIMERS];
    union UserData {
      double d;
      float f[2];
      unsigned u[2];
    } data;
    unsigned blockNr;
  };

  struct DevTimerData {
    unsigned long long time[N_MAX_DEVICE_TIMERS];
  };

  struct DevTraceBlockData {
    unsigned long long start;
    unsigned long long stop;
    unsigned start_smid;
    unsigned stop_smid;
  };

  constexpr int N_MAX_TRACER_OBJECTS = 1024*1024;
  constexpr int N_MAX_BLOCKS = 65*1024;
  extern __device__ DevTimerData d_timer_data[];
  extern __device__ DevTraceData d_trace_data[];
  extern __device__ DevTraceBlockData d_block_data[];
  extern Array<DevTimerData> timer_data;
  extern Array<DevTraceData> trace_data;
  extern Array<DevTraceBlockData> block_data;


  __device__ inline unsigned GetSMID() {
    unsigned ret;
    asm("mov.u32 %0, %smid;" : "=r"(ret) );
    return ret;
  }

  struct DeviceBlockRegionTracer {
    int blockNr, threadNr;
    __device__ DeviceBlockRegionTracer(int nblocks, int blockNr_, int threadNr_)
      : blockNr(blockNr_), threadNr(threadNr_)
    {
      if(threadNr == 0)
      {
        d_block_data[blockNr].start = clock64();
        d_block_data[blockNr].start_smid = GetSMID();
        if(blockNr == 0)
        {
          d_trace_data[N_MAX_TRACER_OBJECTS].blockNr = nblocks;
        }
      }
      __syncwarp();
    }

    __device__ ~DeviceBlockRegionTracer()
    {
      d_block_data[blockNr].stop = clock64();
      d_block_data[blockNr].stop_smid = GetSMID();
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
    int timer_nr, id;
    bool active;
    DevTraceData::UserData & data;

    __device__ DeviceRegionTracer( const DeviceBlockRegionTracer & tr, int timer_nr_, int id_)
      : active(tr.threadNr == 0), timer_nr(timer_nr_), id(id_), data(d_trace_data[id_].data)
    {
      if(active)
      {
        d_trace_data[id].start[timer_nr] = clock64();
        d_trace_data[id].start_smid[timer_nr] = GetSMID();
        d_trace_data[id].blockNr = tr.blockNr;
      }
    }

    __device__ ~DeviceRegionTracer()
    {
      if(active)
      {
        d_trace_data[id].stop[timer_nr] = clock64();
        d_trace_data[id].stop_smid[timer_nr] = GetSMID();
      }
    }

  };
#else
  struct DeviceBlockRegionTracer
  {
    __device__ DeviceBlockRegionTracer( int , int, int ) {}    
  };
  struct DeviceRegionTimer
  {
    __device__ DeviceRegionTimer( int , int ) {}
  };
  struct DeviceRegionTracer
  {
    // __device__ DeviceRegionTracer( bool active_, int timer_nr_, int id_) {}
    __device__ DeviceRegionTracer( const DeviceBlockRegionTracer & tr, int timer_nr_, int id_) { }  
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

  public:
    CudaRegionTimer (ngcore::Timer<> & atimer) : timer(atimer) {
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
      ProcessTracingData();
#endif // NGS_CUDA_DEVICE_TIMERS
      timer.Stop();
      is_already_stopped = true;
    }

    void ProcessTracingData();

    ngcore::TTimePoint TimeD2H( long long ticks )
    {
      double t_seconds = 1.0*ticks/gpu_clock + 1e-6; // add 1us for kernel startup
      return t_kernel_start + t_seconds / ngcore::seconds_per_tick;
    }

    static void SetCudaTimer(bool enabled){
      is_cuda_timer_enabled = enabled;
    }
  };

  void TimeProfiler();
}

#endif // NGS_CUDA_PROFILER_HPP
