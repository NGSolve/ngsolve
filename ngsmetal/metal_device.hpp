#ifndef FILE_METAL_DEVICE_HPP
#define FILE_METAL_DEVICE_HPP

/*********************************************************************/
/* File:   metal_device.hpp                                          */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   29. Aug. 2026                                             */
/*********************************************************************/

/*
  Metal implementation of the ngs_gpu interface (ngstd/gpuwrapper.hpp).
*/

#include <gpuwrapper.hpp>

namespace MTL { class Buffer; class Device; class CommandQueue; }

namespace ngsmetal
{
  // installs the ngs_gpu device creator
  void InitMetalDevice();

  // raw Metal handles (lazily created singletons)
  MTL::Device * GetDevice();
  MTL::CommandQueue * GetCommandQueue();

  // the MTL buffer behind an ngs_gpu buffer of this backend
  MTL::Buffer * GetMTLBuffer (const ngs_gpu::Buffer & buf);
}

#endif
