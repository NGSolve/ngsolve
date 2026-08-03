#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>


#include "metal_vector.hpp"


namespace ngsmetal
{

  MTL::ComputePipelineState* saxpy_pipelineState = nullptr;
  MTL::ComputePipelineState* scale_pipelineState = nullptr;



  auto lam = []()
  {
    
    cout << "init vec kernels" << endl;
    
    const char* shaderSource = R"(
        #include <metal_stdlib>
        using namespace metal;

        kernel void saxpy(device const float* x [[buffer(0)]],
                          device float* y       [[buffer(1)]],
                          constant float& a     [[buffer(2)]],
                          uint id               [[thread_position_in_grid]]) {
            y[id] = a * x[id] + y[id];
        }

        kernel void scale(device float* x       [[buffer(0)]],
                          constant float& a     [[buffer(1)]],
                          uint id                [[thread_position_in_grid]]) {
            x[id] *= a;
        }
    )";


    
    NS::Error* error = nullptr;
    NS::String* mslCode = NS::String::string(shaderSource, NS::UTF8StringEncoding);
    MTL::Library* library = GetDevice()->newLibrary(mslCode, nullptr, &error);

    if (!library) {
      throw Exception("Shader compile error: " 
                      +string(error->localizedDescription()->utf8String()));
    }

    
    {
      NS::String* funcName = NS::String::string("saxpy", NS::UTF8StringEncoding);
      auto saxpyFunc = library->newFunction(funcName);
      saxpy_pipelineState = GetDevice()->newComputePipelineState(saxpyFunc, &error);            
    }
    {
      NS::String* funcName = NS::String::string("scale", NS::UTF8StringEncoding);
      auto scaleFunc = library->newFunction(funcName);
      scale_pipelineState = GetDevice()->newComputePipelineState(scaleFunc, &error);      
    }
    return 0;
  }();


  
  MetalVector :: MetalVector (size_t s)
  {
    this -> size = s;
    buffer = GetDevice()->newBuffer(s*sizeof(float), MTL::ResourceStorageModeShared);
  }


  MetalVector :: MetalVector (const BaseVector& v)
  {
    this -> size = v.Size();
    buffer = GetDevice()->newBuffer(this->size*sizeof(float), MTL::ResourceStorageModeShared);

    std::visit([&](auto vme) {
      
      FlatVector<float> me(this->size, (float*)buffer->contents());
      auto you = v.FV<decltype(vme)>();      

      if constexpr (requires { me(0) = you(0); })
        {
          me = you;
        }
      else 
        throw Exception("MetalVector::MetalVector - illegal vector type");
      
    }, v.GetScalarType());
  }
  

  BaseVector & MetalVector :: Scale (double scal)
  {
    MTL::CommandBuffer* commandBuffer = GetCommandQueue()->commandBuffer();
    MTL::ComputeCommandEncoder* encoder = commandBuffer->computeCommandEncoder();

    encoder->setComputePipelineState(scale_pipelineState);
    encoder->setBuffer(buffer, 0, 0);
    float fscal = scal;
    encoder->setBytes(&fscal, sizeof(float), 1);

    MTL::Size gridSize = MTL::Size(size, 1, 1);
    NS::UInteger maxThreads = scale_pipelineState->maxTotalThreadsPerThreadgroup();
    NS::UInteger threadsGroupDim = (size < maxThreads) ? size : maxThreads;
    MTL::Size threadgroupSize = MTL::Size(threadsGroupDim, 1, 1);
      
    encoder->dispatchThreads(gridSize, threadgroupSize);
    encoder->endEncoding();

    CommitAsync(commandBuffer);    
    return *this;
  }

  BaseVector & MetalVector :: Add (double scal, const BaseVector & v2)
  {
    const MetalVector &mv2 = dynamic_cast<const MetalVector&>(v2);
      
    NS::Error* error = nullptr;      
      
    MTL::CommandBuffer* commandBuffer = GetCommandQueue()->commandBuffer();
    MTL::ComputeCommandEncoder* encoder = commandBuffer->computeCommandEncoder();

    encoder->setComputePipelineState(saxpy_pipelineState);
    encoder->setBuffer(mv2.buffer, 0, 0);
    encoder->setBuffer(buffer, 0, 1);      
    float fscal = scal;
    encoder->setBytes(&fscal, sizeof(float), 2);

    MTL::Size gridSize = MTL::Size(size, 1, 1);
    NS::UInteger maxThreads = saxpy_pipelineState->maxTotalThreadsPerThreadgroup();
    NS::UInteger threadsGroupDim = (size < maxThreads) ? size : maxThreads;
    MTL::Size threadgroupSize = MTL::Size(threadsGroupDim, 1, 1);
      
    encoder->dispatchThreads(gridSize, threadgroupSize);
    encoder->endEncoding();

    CommitAsync(commandBuffer);
    return *this;
  }
    

  

  
}
