#include <ngstd.hpp>
#include <cuda_ngstd.hpp>
#include "cuda_profiler.hpp"

using namespace std;

namespace ngs_cuda
{

  int gpu_clock = 0;

#ifdef NGS_CUDA_DEVICE_TIMERS
  constexpr bool SKIP_BLOCKS_WITH_NO_TRACES = true;
  constexpr bool CHECK_ALL_TRACE_ENTRIES = true;

  __device__ DevTimerData d_timer_data[N_MAX_BLOCKS * N_MAX_DEVICE_TIMERS];
  __device__ DevTraceData d_trace_data[N_MAX_TRACER_OBJECTS+1];
  __device__ DevTraceBlockData d_block_data[N_MAX_BLOCKS];


  Array<DevTimerData> timer_data{N_MAX_BLOCKS * N_MAX_DEVICE_TIMERS};
  Array<DevTraceData> trace_data{N_MAX_TRACER_OBJECTS+1};
  Array<DevTraceBlockData> block_data{N_MAX_BLOCKS};
#endif // NGS_CUDA_DEVICE_TIMERS


  bool CudaRegionTimer :: is_cuda_timer_enabled = false;

  void CudaRegionTimer :: ProcessTracingData()
  {
#ifdef NGS_CUDA_DEVICE_TIMERS
      static ngcore::Timer<> t_overhead("CUDA Tracing overhead");
      static ngcore::Timer<> t_memcopy("memcopy");
      static ngcore::Timer<> t_cleanup("cleanup");
      cudaDeviceSynchronize();
      t_overhead.Start();
      detail::StartGPUTimer(t_overhead);
      t_memcopy.Start();
      cudaMemcpyFromSymbol (&trace_data.Last(), d_trace_data, sizeof(DevTraceData), sizeof(DevTraceData)*(trace_data.Size()-1));
      t_memcopy.Stop();
      auto nblocks = trace_data.Last().blockNr;
      if(nblocks>0)
      {
      t_memcopy.Start();
        cudaMemcpyFromSymbol(trace_data.Data(), d_trace_data, sizeof(DevTraceData)*(trace_data.Size()-1));
      t_memcopy.Stop();
      t_memcopy.Start();
        cudaMemcpyFromSymbol(timer_data.Data(), d_timer_data, sizeof(DevTimerData)*timer_data.Size());
      t_memcopy.Stop();
      t_memcopy.Start();
        cudaMemcpyFromSymbol(block_data.Data(), d_block_data, sizeof(DevTraceBlockData)*nblocks);
      t_memcopy.Stop();
        cudaDeviceSynchronize();

        auto blocks = block_data.Range(nblocks);

        unsigned max_smid = 0;
        for(auto b : blocks)
        {
          max_smid = std::max(max_smid, b.start_smid);
          max_smid = std::max(max_smid, b.stop_smid);
        }

        // cout << "max smid: "  << max_smid << endl;;
        Array<unsigned long long> smid_start_time(max_smid+1);
        smid_start_time = std::numeric_limits<unsigned long long>::max();
        for(auto b : blocks)
          smid_start_time[b.start_smid] = std::min(smid_start_time[b.start_smid], b.start);

        for(auto & b : blocks)
        {
          b.start -= smid_start_time[b.start_smid];
          b.stop -= smid_start_time[b.stop_smid];
        }

        if(ngcore::trace)
        {
          Array<int> container_sm;
          for(int i : Range(max_smid+1))
            container_sm.Append(ngcore::trace->AddUserContainer("SM " + ToString(i)));

          Array<bool> use_block(blocks.Size());
          use_block = true;
          if constexpr(SKIP_BLOCKS_WITH_NO_TRACES)
          {
            use_block = false;
            for(auto & tr : trace_data)
            {
              if(tr.start[0]==0)
              {
                if constexpr(CHECK_ALL_TRACE_ENTRIES)
                  continue;
                else
                  break;
              }
              use_block[tr.blockNr] = true;
            }
          }

          Array<int> container_block(nblocks);
          container_block = -1;
          for(int i : Range(blocks))
          {
            if(!use_block[i])
              continue;
            auto & b = blocks[i];
            string name = "Block " + ToString(i);
            if(b.start_smid == b.stop_smid)
              container_block[i] = ngcore::trace->AddUserContainer("Block " + ToString(i), container_sm[b.start_smid]);
            else
            {
              name += ", SMIDs: " + ToString(b.start_smid) + ',' + ToString(b.stop_smid);
              container_block[i] = ngcore::trace->AddUserContainer(name);
            }
            ngcore::trace->AddUserEvent({TimeD2H(b.start), TimeD2H(b.stop), "Active", container_block[i], i});
          }

          for(int id : Range(trace_data))
          {
            auto & tr = trace_data[id];
            if(tr.start[0]==0)
            {
              if constexpr(CHECK_ALL_TRACE_ENTRIES)
                continue;
              else
                break;
            }
            for(int ti : Range(N_MAX_DEVICE_TIMERS))
            {
              if(tr.start[ti]==0)
                break;
              auto t0 = smid_start_time[tr.start_smid[ti]];
              auto t1 = smid_start_time[tr.stop_smid[ti]];
              stringstream s;
              s << "Timer " << ti;
              if(tr.start_smid[ti] != tr.stop_smid[ti])
                s << "SMIDs: " << tr.start_smid[ti] << ',' << tr.stop_smid[ti];

              ngcore::trace->AddUserEvent({TimeD2H(tr.start[ti]-t0), TimeD2H(tr.stop[ti]-t1), s.str(), container_block[tr.blockNr], id});
            }
          }
        }

        // cleanup for next run
        t_cleanup.Start();
        timer_data = DevTimerData{{0}};
        cudaMemcpyToSymbol (d_timer_data, timer_data.Data(), sizeof(DevTimerData)*timer_data.Size());
        trace_data = DevTraceData{0,0,0};
        cudaMemcpyToSymbol (d_trace_data, trace_data.Data(), sizeof(DevTraceData)*trace_data.Size());
        block_data = DevTraceBlockData{0,0,0,0};
        cudaMemcpyToSymbol (d_block_data, block_data.Data(), sizeof(DevTraceBlockData)*block_data.Size());
        t_cleanup.Stop();
      }
      
      detail::StartGPUTimer(t_overhead);
      t_overhead.Stop();
#endif // NGS_CUDA_DEVICE_TIMERS
  }

  __global__ void SmallKernel (long long *clock)
  {
    if(clock)
      *clock = clock64();
  }

  void TimeProfiler() {
    static Timer t("cudaDeviceSynchronize");
    for(auto i : Range(10))
      SmallKernel<<<1,1>>>(nullptr);
    {
      RegionTimer rt(t);
      cudaDeviceSynchronize();
    }
    {
      static Timer t("Empty CudaRegionTimer"); RegionTimer rt(t);
      for(auto i : Range(10))
        CudaRegionTimer crt(t);
    }
    {
      RegionTimer rt(t);
      cudaDeviceSynchronize();
    }
    for(auto j : Range(5))
    {
      {
        static Timer t("CudaRegionTimer with small kernel"); RegionTimer rt(t);
        for(auto i : Range(100))
        {
          CudaRegionTimer crt(t);
          SmallKernel<<<1,1>>>(nullptr);
        }
      }
      {
        RegionTimer rt(t);
        cudaDeviceSynchronize();
      }
      {
        static Timer t("One CudaRegionTimer with 10 small kernel s"); RegionTimer rt(t);
        CudaRegionTimer crt(t);
        for(auto i : Range(100))
          SmallKernel<<<1,1>>>(nullptr);
      }
      {
        RegionTimer rt(t);
        cudaDeviceSynchronize();
      }
    }

  }
} // namespace ngs_cuda
