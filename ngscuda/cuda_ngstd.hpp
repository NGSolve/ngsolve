#ifndef CUDA_NGSTD_HPP
#define CUDA_NGSTD_HPP

#include <cuda_runtime.h>
#include <ngstd.hpp>

#include "cuda_core.hpp"
#include "cuda_profiler.hpp"

namespace ngs_cuda
{
  using namespace ngstd;

  
  extern int gpu_clock;
  void InitCUDA (int verbose = 2);
  void WarmupCudaModule ();   // first launch of a kernel of this library

}

namespace ngcore {
  template <typename T>
  struct IsSafe<ngs_cuda::Dev<T>> {
    constexpr operator bool() const { return true; }
  };
}

namespace ngstd
{
  template <typename T>
  struct my_is_integral<ngs_cuda::Dev<T>> : my_is_integral<T>{};
}

#endif
