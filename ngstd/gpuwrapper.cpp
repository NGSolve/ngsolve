/*********************************************************************/
/* File:   gpuwrapper.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   29. Aug. 2025                                             */
/*********************************************************************/

#include "gpuwrapper.hpp"

namespace ngs_gpu
{
  static DeviceCreator device_creator;
  static shared_ptr<Device> the_device;

  void SetDeviceCreator (DeviceCreator creator)
  {
    device_creator = std::move(creator);
    the_device.reset();
  }

  bool HasDevice() { return bool(device_creator); }

  shared_ptr<Device> GetDevice()
  {
    if (!the_device && device_creator)
      the_device = device_creator();
    return the_device;
  }
}
