#include <ngstd.hpp>
#include <cuda_ngstd.hpp>
#include "cuda_profiler.hpp"

using namespace std;

namespace ngs_cuda
{

  int gpu_clock = 0;

#ifdef NGS_CUDA_DEVICE_TIMERS
  constexpr bool SKIP_BLOCKS_WITH_NO_TRACES = true;

  __device__ DevTimerData d_timer_data[N_MAX_BLOCKS * N_MAX_DEVICE_TIMERS];
  __device__ DevTraceData d_trace_data[N_MAX_TRACER_OBJECTS+1];
  __device__ DevTraceBlockData d_block_data[N_MAX_BLOCKS];
  __device__ DevTraceState d_trace_state;


  Array<DevTimerData> timer_data{N_MAX_BLOCKS * N_MAX_DEVICE_TIMERS};
  Array<DevTraceData> trace_data{N_MAX_TRACER_OBJECTS+1};
  Array<DevTraceBlockData> block_data{N_MAX_BLOCKS};

  void * GetDevTraceDataPtr()
  {
    void * p = nullptr;
    cudaGetSymbolAddress(&p, d_trace_data);
    return p;
  }
  void * GetDevTraceBlockDataPtr()
  {
    void * p = nullptr;
    cudaGetSymbolAddress(&p, d_block_data);
    return p;
  }
  void * GetDevTraceStatePtr()
  {
    void * p = nullptr;
    cudaGetSymbolAddress(&p, d_trace_state);
    return p;
  }
#endif // NGS_CUDA_DEVICE_TIMERS


  bool CudaRegionTimer :: is_cuda_timer_enabled = false;

  void CudaRegionTimer :: ProcessTracingData(const std::initializer_list<const char *> * timer_names)
  {
#ifdef NGS_CUDA_DEVICE_TIMERS
      static ngcore::Timer<> t_overhead("CUDA Tracing overhead");
      static ngcore::Timer<> t_memcopy("memcopy");
      static ngcore::Timer<> t_cleanup("cleanup");
      cudaDeviceSynchronize();
      t_overhead.Start();
      detail::StartGPUTimer(t_overhead);
      t_memcopy.Start();
      DevTraceState state;
      cudaMemcpyFromSymbol (&state, d_trace_state, sizeof(state));
      t_memcopy.Stop();
      auto nblocks = state.nblocks;
      size_t n_trace_entries = std::min(size_t(state.chunk_counter) * TRACE_CHUNK_SIZE,
                                        size_t(N_MAX_TRACER_OBJECTS));
      if(size_t(state.chunk_counter) * TRACE_CHUNK_SIZE > size_t(N_MAX_TRACER_OBJECTS))
        cerr << "CUDA tracing: trace buffer capacity exceeded, some events were dropped" << endl;
      if(nblocks>0)
      {
      t_memcopy.Start();
        if(n_trace_entries > 0)
          cudaMemcpyFromSymbol(trace_data.Data(), d_trace_data, sizeof(DevTraceData)*n_trace_entries);
      t_memcopy.Stop();
      t_memcopy.Start();
        cudaMemcpyFromSymbol(timer_data.Data(), d_timer_data, sizeof(DevTimerData)*timer_data.Size());
      t_memcopy.Stop();
      t_memcopy.Start();
        cudaMemcpyFromSymbol(block_data.Data(), d_block_data, sizeof(DevTraceBlockData)*nblocks);
      t_memcopy.Stop();
        cudaDeviceSynchronize();

        auto blocks = block_data.Range(0ul, nblocks);
        auto traces = trace_data.Range(0ul, n_trace_entries);

        unsigned max_smid = 0;
        for(auto b : blocks)
        {
          max_smid = std::max(max_smid, b.start_smid);
          max_smid = std::max(max_smid, b.stop_smid);
        }
        // cout << "max smid: "  << max_smid << endl;;

        // Per-SM linear calibration clock64 -> %globaltimer. clock64 is a cheap,
        // high-resolution per-SM cycle counter (used for the fine-grained timers)
        // but is not comparable across SMs; %globaltimer is global but coarse
        // (~tens of ns). Every block records both clocks at start and stop. Since
        // clock64 is one monotonic counter per SM, all blocks on an SM lie on a
        // single line; we fit it from the widest baseline on that SM (the min-
        // and max-clock64 anchors across all its blocks), so endpoint quantization
        // of %globaltimer is negligible.
        constexpr auto ull_max = std::numeric_limits<unsigned long long>::max();
        Array<unsigned long long> sm_c0(max_smid+1), sm_g0(max_smid+1); // min-clock anchor
        Array<unsigned long long> sm_c1(max_smid+1), sm_g1(max_smid+1); // max-clock anchor
        sm_c0 = ull_max; sm_g0 = 0; sm_c1 = 0; sm_g1 = 0;

        auto AddAnchor = [&](unsigned smid, unsigned long long c, unsigned long long g)
        {
          if(c < sm_c0[smid]) { sm_c0[smid] = c; sm_g0[smid] = g; }
          if(c > sm_c1[smid]) { sm_c1[smid] = c; sm_g1[smid] = g; }
        };
        for(auto b : blocks)
        {
          AddAnchor(b.start_smid, b.start, b.start_global);
          AddAnchor(b.stop_smid,  b.stop,  b.stop_global);
        }

        // map a clock64 value on the given SM to absolute %globaltimer nanoseconds
        auto ClockToGlobal = [&](unsigned smid, unsigned long long c) -> long long
        {
          auto c0 = sm_c0[smid], c1 = sm_c1[smid];
          if(c1 == c0)                 // degenerate: single anchor, no baseline
            return (long long) sm_g0[smid];
          double slope = double(sm_g1[smid] - sm_g0[smid]) / double(c1 - c0);
          return (long long) sm_g0[smid]
               + (long long) ((double(c) - double(c0)) * slope + 0.5);
        };

        // earliest global time overall (%globaltimer is comparable across SMs)
        long long global_start_time = std::numeric_limits<long long>::max();
        for(unsigned s = 0; s <= max_smid; s++)
          if(sm_c0[s] != ull_max)
            global_start_time = std::min(global_start_time, (long long) sm_g0[s]);

        if(ngcore::trace)
        {
          Array<int> container_sm;
          for(int i : Range(max_smid+1))
            container_sm.Append(ngcore::trace->AddUserContainer("SM " + ToString(i)));

          Array<bool> use_block(blocks.Size());
          use_block = true;
          if constexpr(SKIP_BLOCKS_WITH_NO_TRACES)
          {
            // chunks of different blocks interleave and the last chunk of each
            // block has an unused (zeroed) tail, so scan all reserved entries
            use_block = false;
            for(auto & tr : traces)
              if(tr.start != 0)
                use_block[tr.blockNr] = true;
          }

          // Instead of one container per block, pack blocks into a small number
          // of reusable lane-containers "SM i Block j" per SM. A lane can be
          // reused by a later block whose start time is after the stop time of
          // the last block already placed in that lane.
          Array<int> container_block(nblocks);
          container_block = -1;

          // process used blocks in order of their start time for tight packing
          Array<int> order;
          for(int i : Range(blocks))
            if(use_block[i])
              order.Append(i);
          QuickSort(order, [&](int a, int b) { return blocks[a].start < blocks[b].start; });

          // Pass 1: assign each block to a lane per SM, without creating any
          // containers yet - just record the lane and how many lanes each SM needs.
          Array<Array<unsigned long long>> lane_stop(max_smid+1);
          Array<int> lane_of_block(nblocks);
          lane_of_block = -1;
          for(int i : order)
          {
            auto & b = blocks[i];
            auto sm = b.start_smid;

            // find a lane whose last block already stopped by the time this one
            // starts (half-open intervals: touching blocks are not concurrent)
            int lane = -1;
            for(int l : Range(lane_stop[sm]))
              if(lane_stop[sm][l] <= b.start)
              {
                lane = l;
                break;
              }

            if(lane == -1)
            {
              lane = lane_stop[sm].Size();
              lane_stop[sm].Append(0ull);
            }

            lane_stop[sm][lane] = b.stop;
            lane_of_block[i] = lane;
          }

          // Pass 2: create the lane containers sorted by SM (and by lane within
          // each SM), so the generated container order is deterministic.
          Array<Array<int>> lane_container(max_smid+1);
          for(int sm : Range(max_smid+1))
            for(int lane : Range(lane_stop[sm].Size()))
              lane_container[sm].Append(ngcore::trace->AddUserContainer(
                    "SM " + ToString(sm) + " Block " + ToString(lane), container_sm[sm]));

          // Pass 3: emit the block "Active" events
          for(int i : order)
          {
            auto & b = blocks[i];
            container_block[i] = lane_container[b.start_smid][lane_of_block[i]];
            auto g_start = ClockToGlobal(b.start_smid, b.start) - global_start_time;
            auto g_stop  = ClockToGlobal(b.stop_smid,  b.stop)  - global_start_time;
            ngcore::trace->AddUserEvent({TimeD2H(g_start), TimeD2H(g_stop), "Active", container_block[i], i});
          }

          for(size_t id : Range(traces))
          {
            auto & tr = traces[id];
            if(tr.start == 0)
              continue;
            stringstream s;
            if(timer_names != nullptr && tr.timer_nr < timer_names->size())
              s << (*timer_names).begin()[tr.timer_nr] << ' ';
            else
              s << "Timer " << tr.timer_nr;

            // a block never migrates SMs, so its start_smid is valid for all
            // trace events recorded within it
            auto smid = blocks[tr.blockNr].start_smid;
            auto g_start = ClockToGlobal(smid, tr.start) - global_start_time;
            auto g_stop  = ClockToGlobal(smid, tr.stop)  - global_start_time;
            ngcore::trace->AddUserEvent({TimeD2H(g_start), TimeD2H(g_stop), s.str(), container_block[tr.blockNr], int(id)});
          }
        }

        // cleanup for next run (only the actually used entries need zeroing)
        t_cleanup.Start();
        timer_data = DevTimerData{{0}};
        cudaMemcpyToSymbol (d_timer_data, timer_data.Data(), sizeof(DevTimerData)*timer_data.Size());
        if(n_trace_entries > 0)
        {
          memset(trace_data.Data(), 0, sizeof(DevTraceData)*n_trace_entries);
          cudaMemcpyToSymbol (d_trace_data, trace_data.Data(), sizeof(DevTraceData)*n_trace_entries);
        }
        block_data = DevTraceBlockData{0,0,0,0,0,0};
        cudaMemcpyToSymbol (d_block_data, block_data.Data(), sizeof(DevTraceBlockData)*nblocks);
        DevTraceState zero_state{0, 0};
        cudaMemcpyToSymbol (d_trace_state, &zero_state, sizeof(zero_state));
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
