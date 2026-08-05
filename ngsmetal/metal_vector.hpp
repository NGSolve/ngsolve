#ifndef FILE_METAL_VECTOR_HPP
#define FILE_METAL_VECTOR_HPP

#include "ngsmetal.hpp"
#include <basevector.hpp>

namespace ngsmetal
{
  using namespace ngla;

  
  class MetalVector : public S_BaseVector<float>
  {
    bool shared;
    MTL::Buffer* buffer;
    mutable MTL::CommandBuffer* lastCommandBuffer = nullptr; // Track the latest GPU job
    
  public:
    MetalVector (size_t s, bool shared=true);
    MetalVector (const BaseVector& v, bool shared=true);    
    MTL::Buffer* GetBuffer() const { return buffer; }

    virtual BaseVector & Scale (double scal) override;
    virtual BaseVector & SetScalar (double scal) override;
    virtual BaseVector & Set (double scal, const BaseVector & v) override;
    virtual BaseVector & Add (double scal, const BaseVector & v2) override;
    

    virtual size_t Size() const  { return size; }
    virtual void * Memory () const override
    {
      WaitUntilCompleted(); // Block CPU until GPU operations are done
      if (!shared) throw Exception("Buffer of private MetalVector is not host accessible");
      return buffer->contents();
    }
    
    virtual AutoVector CreateVector () const override
    {
      return make_unique<MetalVector>(Size(), false);
    }

    virtual ostream& Print (ostream &ost) const override
    {
      WaitUntilCompleted(); // Block CPU until GPU operations are done      
      ost << "MetalVector, size = " << size << endl;
      if (!shared) throw Exception("Buffer of private MetalVector is not host accessible");      
      float * values = (float*)buffer->contents();      
      for (size_t i = 0; i < size; i++)
        ost << values[i] << "\n";
      return ost;
    }



    // Helper to ensure GPU finishes work before CPU access
    void WaitUntilCompleted() const
    {
      if (lastCommandBuffer)
        {
          lastCommandBuffer->waitUntilCompleted();
          lastCommandBuffer->release(); // Release retained reference
          lastCommandBuffer = nullptr;
        }
    }

    // Helper to store and commit command buffers
    void CommitAsync(MTL::CommandBuffer* cb)
    {
      // If there was a previously tracked buffer that hasn't completed, clean it up
      if (lastCommandBuffer)
        {
          lastCommandBuffer->release();
        }

      lastCommandBuffer = cb;
      lastCommandBuffer->retain(); // Keep it alive so we can wait on it later
      lastCommandBuffer->commit(); // Dispatch to GPU asynchronously (non-blocking)
    }
    
  };
}

#endif
