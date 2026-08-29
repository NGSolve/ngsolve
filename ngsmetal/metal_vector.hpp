#ifndef FILE_METAL_VECTOR_HPP
#define FILE_METAL_VECTOR_HPP

/*
  MetalVector is the backend-independent ngla::DeviceVector on float,
  kept as an alias for the python bindings and existing code.
*/

#include <devicevector.hpp>

namespace ngsmetal
{
  using namespace ngla;

  typedef DeviceVector<float> MetalVector;
}

#endif
