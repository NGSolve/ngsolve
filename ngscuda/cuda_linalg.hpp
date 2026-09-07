#ifndef CUDA_LINALG_HPP
#define CUDA_LINALG_HPP

// partial override of overloaded function (MultAdd)
#pragma nv_diag_suppress 611
#pragma nv_diag_suppress 20013

#include <la.hpp>

#include <cuda_runtime.h>

#include "cuda_ngstd.hpp"

  
#include "cuda_ngbla.hpp"
#include "cuda_device.hpp"


namespace ngla
{
  using namespace ngs_cuda;

  void InitCuLinalg();


  // typed device access to any DeviceVector<double> (DeviceVectorWrapper, ...),
  // with the transfers the access implies

  inline Dev<double> * DevPtr (const DeviceVector<double> & v)
  {
    return (Dev<double>*)ngs_cuda::BufferDevPtr(*v.DevBufferRO()) + v.DevOffset();
  }

  // kernel reads and writes
  inline FlatVector<Dev<double>> FVDev (const DeviceVector<double> & v)
  {
    auto ptr = (Dev<double>*)ngs_cuda::BufferDevPtr(*v.DevBufferRW()) + v.DevOffset();
    return { v.Size(), ptr };
  }

  // kernel only reads
  inline FlatVector<Dev<double>> FVDevRO (const DeviceVector<double> & v)
  {
    return { v.Size(), DevPtr(v) };
  }


  class DevMatrix : public BaseMatrix
  {
  public:
    DevMatrix() { }

    AutoVector CreateRowVector() const override { return make_unique<DeviceVector<double>>(Width(), PreferredMemType()); }
    AutoVector CreateColVector() const override { return make_unique<DeviceVector<double>>(Height(), PreferredMemType()); }
  };

  shared_ptr<BaseMatrix> CreateDevMatrix (BaseMatrix &mat);

}


#endif
