#ifndef FILE_METAL_BTDTB_HPP
#define FILE_METAL_BTDTB_HPP

#include <comp.hpp>
#include "metal_vector.hpp"

using namespace ngcomp;



namespace ngsmetal
{
  class MetalBTDTBMatrix : public BaseMatrix
  {
    int h, w;
    MTL::Function* ApplyBTDTB_Func = nullptr;

    MTL::Buffer* buffer_dofx;
    MTL::Buffer* buffer_dofy;
    MTL::Buffer* buffer_bmatx;
    MTL::Buffer* buffer_bmaty;

    MTL::Buffer* buffer_weights;
    MTL::Buffer* buffer_Jacobi;
    MTL::Buffer* buffer_JacobiDets;
    MTL::Buffer* debug;


    MTL::ComputePipelineState* pipelineState;
 
    int ne, BS_els;
    
  public:
    MetalBTDTBMatrix (const BaseMatrix& mat);


    
    
    AutoVector CreateRowVector() const override {
      return make_unique<MetalVector>(w);
    }
    AutoVector CreateColVector() const override {
      return make_unique<MetalVector>(h);
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override;
    
  };
}


#endif
