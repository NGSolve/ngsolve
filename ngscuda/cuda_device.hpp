#ifndef FILE_CUDA_DEVICE_HPP
#define FILE_CUDA_DEVICE_HPP

/*********************************************************************/
/* File:   cuda_device.hpp                                           */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   29. Aug. 2026                                             */
/*********************************************************************/

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

  // raw cuda device pointer of a buffer allocated by this backend,
  // for handing ngs_gpu storage to cuBLAS/cuSPARSE and typed kernels
  void * BufferDevPtr (ngs_gpu::Buffer & buf);
}

#endif
