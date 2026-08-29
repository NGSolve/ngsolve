#ifndef FILE_CUDA_DEVICE_HPP
#define FILE_CUDA_DEVICE_HPP

/*
  CUDA implementation of the ngs_gpu interface (ngstd/gpuwrapper.hpp).

  Kernels are compiled at runtime with NVRTC and launched through the
  driver API, mirroring the metal backend. Entry points must be
  extern "C" __global__, otherwise the name is mangled and not findable.
*/

#include <gpuwrapper.hpp>

namespace ngs_cuda
{
  // installs the ngs_gpu device creator
  void InitCudaDevice();
}

#endif
