#ifndef CUDA_NGSTD_HPP
#define CUDA_NGSTD_HPP

#include <cuda_runtime.h>
#include <ngstd.hpp>

#include "cuda_core.hpp"

namespace ngs_cuda
{
  using namespace ngstd;

  
  void InitCUDA (int verbose = 2);

}

namespace ngcore {
  template <typename T>
  struct IsSafe<ngs_cuda::Dev<T>> {
    constexpr operator bool() const { return true; }
  };
}

namespace ngcore
{
  template <typename T>
  struct my_is_integral<ngs_cuda::Dev<T>> : my_is_integral<T>{};
}

#endif
