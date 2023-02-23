#ifdef CUDA

#include <ngstd.hpp>
#include <cuda_ngstd.hpp>

using namespace std;

namespace ngs_cuda
{

  // Print device properties
  void printDevProp(cudaDeviceProp devProp)
  {
    printf("Major revision number:         %d\n",  devProp.major);
    printf("Minor revision number:         %d\n",  devProp.minor);
    printf("Name:                          %s\n",  devProp.name);
    printf("Total global memory:           %lu\n",  devProp.totalGlobalMem);
    printf("Total shared memory per block: %lu\n",  devProp.sharedMemPerBlock);
    printf("Total registers per block:     %d\n",  devProp.regsPerBlock);
    printf("Warp size:                     %d\n",  devProp.warpSize);
    printf("Maximum memory pitch:          %lu\n",  devProp.memPitch);
    printf("Maximum threads per block:     %d\n",  devProp.maxThreadsPerBlock);
    for (int i = 0; i < 3; ++i)
      printf("Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
    for (int i = 0; i < 3; ++i)
      printf("Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
    printf("Clock rate:                    %d\n",  devProp.clockRate);
    printf("Total constant memory:         %lu\n",  devProp.totalConstMem);
    printf("Texture alignment:             %lu\n",  devProp.textureAlignment);
    printf("Concurrent copy and execution: %s\n",  (devProp.deviceOverlap ? "Yes" : "No"));
    printf("Number of multiprocessors:     %d\n",  devProp.multiProcessorCount);
    printf("Kernel execution timeout:      %s\n",  (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
    return;
  }

  int gpu_clock = 0;

  void InitCUDA (int verbose)
  {
    int devCount;
    printf("CUDA Device Query...\n"); 
    cudaGetDeviceCount(&devCount);
    if (devCount == 1)
      printf("There is %d CUDA device.\n", devCount);
    else
      printf("There are %d CUDA devices.\n", devCount);

    for (int i = 0; i < devCount; ++i)
      {
        // Get device properties
        cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, i);
        if(i==0)
          gpu_clock = devProp.clockRate * 1000;

        if (verbose == 1)
          {
            cout << "CUDA Device " << i << ": " << devProp.name 
                 << ", cap " << devProp.major << "." << devProp.minor << endl;
          }
        if (verbose >= 2)
          {
            cout << endl << "CUDA Device " << i << endl;
            printDevProp(devProp);
          }
      }

    int dev_id = 0;
    cudaGetDevice(&dev_id);
    if(verbose >= 1)
      cout << "Using device " << dev_id << endl;

    cudaDeviceSetSharedMemConfig ( cudaSharedMemBankSizeEightByte );
  }

  DevStackMemory stackmemory;

#ifdef NGS_CUDA_DEVICE_TIMERS
  __device__ DevTimerData d_timer_data[N_MAX_BLOCKS * N_MAX_DEVICE_TIMERS];
  __device__ DevTraceData d_trace_data[N_MAX_TRACER_OBJECTS+1];
  __device__ DevTraceBlockData d_block_data[N_MAX_BLOCKS];


  Array<DevTimerData> timer_data{N_MAX_BLOCKS * N_MAX_DEVICE_TIMERS};
  Array<DevTraceData> trace_data{N_MAX_TRACER_OBJECTS+1};
  Array<DevTraceBlockData> block_data{N_MAX_BLOCKS};
#endif // NGS_CUDA_DEVICE_TIMERS

  ngcore::Timer<> CudaRegionTimer :: t_tracing_overhead("CUDA Tracing overhead");

  void CudaRegionTimer :: ProcessTracingData()
  {
#ifdef NGS_CUDA_DEVICE_TIMERS
      cudaDeviceSynchronize();
      t_tracing_overhead.Start();
      StartGPUTimer(t_tracing_overhead);
      cudaMemcpyFromSymbol (&trace_data.Last(), d_trace_data, sizeof(DevTraceData), sizeof(DevTraceData)*(trace_data.Size()-1));
      auto nblocks = trace_data.Last().blockNr;
      if(nblocks>0)
      {
        cudaMemcpyFromSymbol (trace_data.Data(), d_trace_data, sizeof(DevTraceData)*(trace_data.Size()-1));
        cudaMemcpyFromSymbol (timer_data.Data(), d_timer_data, sizeof(DevTimerData)*timer_data.Size());
        cudaMemcpyFromSymbol (block_data.Data(), d_block_data, sizeof(DevTraceBlockData)*nblocks);

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
          for(auto i : Range(blocks))
          {
            auto & b = blocks[i];
            ngcore::trace->AddUserEvent(TimeD2H(b.start), TimeD2H(b.stop), "Block", i);
          }
          for(auto id : Range(trace_data))
          {
            auto & tr = trace_data[id];
            if(tr.start[0]==0)
              break;
            for(auto ti : Range(N_MAX_DEVICE_TIMERS))
            {
              if(tr.start[ti]==0)
                break;
              auto t0 = smid_start_time[tr.start_smid[ti]];
              auto t1 = smid_start_time[tr.stop_smid[ti]];
              stringstream s;
              s << "Timer " << ti;
              // s << "\tId " << id << ti;
              // s << "\tRange: " << tr.data.u[0] <<',' << tr.data.u[1] << ti;

              ngcore::trace->AddUserEvent(TimeD2H(tr.start[ti]-t0), TimeD2H(tr.stop[ti]-t1), s.str(), tr.blockNr);
            }
          }
        }

        // cleanup for next run
        timer_data = DevTimerData{{0}};
        cudaMemcpyToSymbol (d_timer_data, timer_data.Data(), sizeof(DevTimerData)*timer_data.Size());
        trace_data = DevTraceData{0,0,0};
        cudaMemcpyToSymbol (d_trace_data, trace_data.Data(), sizeof(DevTraceData)*trace_data.Size());
        block_data = DevTraceBlockData{0,0,0,0};
        cudaMemcpyToSymbol (d_block_data, block_data.Data(), sizeof(DevTraceBlockData)*block_data.Size());
      }
      // for (auto i : Range(100))
      //   std::cout << "block index " << trace_data[i].blockIdx[0] << ',' << trace_data[i].blockIdx[1] << std::endl;

      // cout << "block data " << block_data << endl;
      
      StartGPUTimer(t_tracing_overhead);
      t_tracing_overhead.Stop();
#endif // NGS_CUDA_DEVICE_TIMERS
  }


  DevBitArray :: DevBitArray (size_t asize)
  {
    SetSize (asize);
  }

  DevBitArray :: DevBitArray (const BitArray & ba)
  {
    (*this) = ba;
  }

  DevBitArray :: ~DevBitArray ()
  {
    if (size)
      cudaFree(dev_data);
  }

  DevBitArray & DevBitArray :: operator= (const BitArray &ba)
  {
    SetSize (ba.Size());

    if (!size)
      return *this;

    cudaMemcpy(dev_data, ba.Data(), Addr(size) + 1, cudaMemcpyHostToDevice);

    return *this;
  }

  void DevBitArray :: SetSize (size_t asize)
  {
    if (size == asize)
      return;

    if (size)
      cudaFree(dev_data);
    
    size = asize;
    cudaMalloc((void**) &dev_data, Addr(size) + 1);
  }

}


#endif
