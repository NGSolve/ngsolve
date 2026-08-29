#ifndef FILE_METAL_DEVICE_HPP
#define FILE_METAL_DEVICE_HPP

/*
  Metal implementation of the ngs_gpu interface (ngstd/gpuwrapper.hpp).
*/

#include <gpuwrapper.hpp>

namespace ngsmetal
{
  // installs the ngs_gpu device creator
  void InitMetalDevice();
}

#endif
