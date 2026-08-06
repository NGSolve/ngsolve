#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>


#include "metal_vector.hpp"


namespace ngsmetal
{

  MTL::ComputePipelineState* saxpy_pipelineState = nullptr;
  MTL::ComputePipelineState* scale_pipelineState = nullptr;
  MTL::ComputePipelineState* setscalar_pipelineState = nullptr;
  MTL::ComputePipelineState* set_pipelineState = nullptr;

  MTL::ComputePipelineState* setscalar4_pipelineState = nullptr;


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

        kernel void set(device const float* x [[buffer(0)]],
                        device float* y       [[buffer(1)]],
                        constant float& a     [[buffer(2)]],
                        uint id               [[thread_position_in_grid]]) {
            y[id] = a * x[id];
        }

        kernel void scale(device float* x       [[buffer(0)]],
                          constant float& a     [[buffer(1)]],
                          uint id                [[thread_position_in_grid]]) {
            x[id] *= a;
        }

        kernel void setscalar(device float* x       [[buffer(0)]],
                          constant float& a     [[buffer(1)]],
                          uint id                [[thread_position_in_grid]]) {
            x[id] = a;
        }

        kernel void setscalar4(device float4* x       [[buffer(0)]],
                          constant float& a     [[buffer(1)]],
                          uint id                [[thread_position_in_grid]]) {
            x[id] = float4(a);
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
    {
      NS::String* funcName = NS::String::string("set", NS::UTF8StringEncoding);
      auto func = library->newFunction(funcName);
      set_pipelineState = GetDevice()->newComputePipelineState(func, &error);      
    }
    {
      NS::String* funcName = NS::String::string("setscalar", NS::UTF8StringEncoding);
      auto func = library->newFunction(funcName);
      setscalar_pipelineState = GetDevice()->newComputePipelineState(func, &error);      
    }
    {
      NS::String* funcName = NS::String::string("setscalar4", NS::UTF8StringEncoding);
      auto func = library->newFunction(funcName);
      setscalar4_pipelineState = GetDevice()->newComputePipelineState(func, &error);      
    }
    return 0;
  }();


  
  MetalVector :: MetalVector (size_t s, bool ashared)
    : shared(ashared)
  {
    this -> size = s;
    size_t s4 = (s+3) & ~(size_t)3;
    if (shared)
      buffer = GetDevice()->newBuffer(s4*sizeof(float), MTL::ResourceStorageModeShared);
    else
      buffer = GetDevice()->newBuffer(s4*sizeof(float), MTL::ResourceStorageModePrivate);
  }


  MetalVector :: MetalVector (const BaseVector& v, bool ashared)
    : shared(ashared)    
  {
    this -> size = v.Size();
    size_t s4 = (this->size+3) & ~(size_t)3;    
    if (shared)    
      buffer = GetDevice()->newBuffer(s4*sizeof(float), MTL::ResourceStorageModeShared);
    else
      buffer = GetDevice()->newBuffer(s4*sizeof(float), MTL::ResourceStorageModePrivate);


    if (!shared)
      {
        auto tmp = make_shared<MetalVector>(v, true);
        this -> Set (1, *tmp);
        return;
      }
    
    std::visit([&](auto vme) {
      
      FlatVector<float> me(this->size, (float*)buffer->contents());
      auto you = v.FV<decltype(vme)>();      

      if constexpr (requires { me(0) = you(0); })
        {
          me = you;
        }
      else 
        throw Exception("MetalVector::MetalVector - illegal vector type");

      size_t s4 = (this->size+3) & ~(size_t)3;    
      FlatVector<float> pad(s4-this->size, (float*)buffer->contents()+this->size);
      pad = 0.0;
    }, v.GetScalarType());
  }
  

  BaseVector & MetalVector :: Scale (double scal)
  {
    if (scal == 1.0)
      return *this;
    
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

  BaseVector & MetalVector :: SetScalar (double scal)
  {
    MTL::CommandBuffer* commandBuffer = GetCommandQueue()->commandBuffer();
    MTL::ComputeCommandEncoder* encoder = commandBuffer->computeCommandEncoder();

    encoder->setComputePipelineState(setscalar4_pipelineState);
    encoder->setBuffer(buffer, 0, 0);
    float fscal = scal;
    encoder->setBytes(&fscal, sizeof(float), 1);

    int size4 = (size+3)/4;
    MTL::Size gridSize = MTL::Size(size4, 1, 1);
    NS::UInteger maxThreads = setscalar_pipelineState->maxTotalThreadsPerThreadgroup();
    NS::UInteger threadsGroupDim = (size4 < maxThreads) ? size4 : maxThreads;
    MTL::Size threadgroupSize = MTL::Size(threadsGroupDim, 1, 1);
      
    encoder->dispatchThreads(gridSize, threadgroupSize);
    encoder->endEncoding();

    CommitAsync(commandBuffer);    
    return *this;
  }

  BaseVector & MetalVector :: Add (double scal, const BaseVector & v2)
  {
    auto mvp = dynamic_cast<const MetalVector*> (&v2);
    if (!mvp)
      {
        auto tmp = make_shared<MetalVector>(v2, true);
        Add (scal, *tmp);
        return *this;
      }
    
    const MetalVector &mv2 = dynamic_cast<const MetalVector&>(v2);
      
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
    

  BaseVector & MetalVector :: Set (double scal, const BaseVector & v2)
  {
    auto mvp = dynamic_cast<const MetalVector*> (&v2);
    if (!mvp)
      {
        auto tmp = make_shared<MetalVector>(v2, true);
        Set (scal, *tmp);
        return *this;
      }

    
    const MetalVector &mv2 = dynamic_cast<const MetalVector&>(v2);
      
    MTL::CommandBuffer* commandBuffer = GetCommandQueue()->commandBuffer();
    MTL::ComputeCommandEncoder* encoder = commandBuffer->computeCommandEncoder();

    encoder->setComputePipelineState(set_pipelineState);
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
