#ifdef CUDA

#include <ngstd.hpp>
#include <cuda_ngstd.hpp>
#include "cuda_profiler.hpp"

using namespace std;

namespace ngs_cuda
{

  // Print device properties
  void printDevProp(cudaDeviceProp devProp, int device_id)
  {
    int clock_rate, kernel_exec_timeout, device_overlap;
    cudaDeviceGetAttribute(&clock_rate, cudaDevAttrClockRate, device_id);
    cudaDeviceGetAttribute(&kernel_exec_timeout, cudaDevAttrKernelExecTimeout, device_id);
    cudaDeviceGetAttribute(&device_overlap, cudaDevAttrConcurrentManagedAccess, device_id);

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
    printf("Clock rate:                    %d\n",  clock_rate);
    printf("Total constant memory:         %lu\n",  devProp.totalConstMem);
    printf("Texture alignment:             %lu\n",  devProp.textureAlignment);
    printf("Concurrent copy and execution: %s\n",  (device_overlap ? "Yes" : "No"));
    printf("Number of multiprocessors:     %d\n",  devProp.multiProcessorCount);
    printf("Kernel execution timeout:      %s\n",  (kernel_exec_timeout ? "Yes" : "No"));
    return;
  }


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
        int clock_rate;
        cudaDeviceGetAttribute(&clock_rate, cudaDevAttrClockRate, i);

        if(i==0)
          gpu_clock = clock_rate * 1000;

        if (verbose == 1)
          {
            cout << "CUDA Device " << i << ": " << devProp.name 
                 << ", cap " << devProp.major << "." << devProp.minor << endl;
          }
        if (verbose >= 2)
          {
            cout << endl << "CUDA Device " << i << endl;
            printDevProp(devProp, i);
          }
      }

    int dev_id = 0;
    cudaGetDevice(&dev_id);
    if(verbose >= 1)
      cout << "Using device " << dev_id << endl;

    cudaDeviceSetSharedMemConfig ( cudaSharedMemBankSizeEightByte );
  }

  DevStackMemory stackmemory;

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
