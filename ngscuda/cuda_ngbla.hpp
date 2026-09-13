#ifndef CUDA_NGBLA
#define CUDA_NGBLA

#include <cuda_runtime.h>

#include <vector.hpp>
#include <matrix.hpp>

#include "cuda_ngstd.hpp"


namespace ngbla
{
  using namespace ngs_cuda;

  // so that FlatVector<Dev<double>> & co. can be formed over device memory
  template<> struct is_scalar_type<Dev<double>> { static constexpr bool value = true; };
}

#endif
