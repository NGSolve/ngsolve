#ifndef FILE_METAL_DEVICE_HPP
#define FILE_METAL_DEVICE_HPP

/*
  Metal implementation of the ngs_gpu interface (ngstd/gpuwrapper.hpp).
*/

#include <gpuwrapper.hpp>

namespace MTL { class Buffer; }

namespace ngsmetal
{
  // installs the ngs_gpu device creator
  void InitMetalDevice();

  // the MTL buffer behind an ngs_gpu buffer of this backend
  MTL::Buffer * GetMTLBuffer (const ngs_gpu::Buffer & buf);
}

#endif
