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


  void InitCUDA (int verbose)
  {
    int devCount;
    printf("CUDA Device Query...\n"); 
    cudaGetDeviceCount(&devCount);
    printf("There are %d CUDA devices.\n", devCount);

    for (int i = 0; i < devCount; ++i)
      {
	// Get device properties
	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, i);

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

    cudaDeviceSetSharedMemConfig ( cudaSharedMemBankSizeEightByte );
  }



}


#endif
