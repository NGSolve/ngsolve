#ifndef FILE_METAL_BTDTB_HPP
#define FILE_METAL_BTDTB_HPP

#include <comp.hpp>
#include <gpuwrapper.hpp>
#include <devicevector.hpp>

using namespace ngcomp;



namespace ngsmetal
{
  template <typename REAL>
  class MetalBTDTBMatrix : public BaseMatrix
  {
    int h, w;
    int ne, warps;

    shared_ptr<ngs_gpu::Library> library;
    shared_ptr<ngs_gpu::Kernel> kernel;
    shared_ptr<ngs_gpu::Queue> queue;

    shared_ptr<ngs_gpu::Buffer> buffer_dofx, buffer_dofy;
    shared_ptr<ngs_gpu::Buffer> buffer_bmatx, buffer_bmaty;
    shared_ptr<ngs_gpu::Buffer> buffer_weights, buffer_Jacobi, buffer_JacobiDets;
    
  public:
    MetalBTDTBMatrix (const BaseMatrix& mat);


    
    
    AutoVector CreateRowVector() const override {
      return make_unique<DeviceVector<float>>(w, PreferredMemType());
    }
    AutoVector CreateColVector() const override {
      return make_unique<DeviceVector<float>>(h, PreferredMemType());
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    
  };
}


#endif
