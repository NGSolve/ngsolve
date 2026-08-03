#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>



#include <comp.hpp>

#include "ngsmetal.hpp"

namespace ngsmetal
{

  MTL::CommandQueue* commandQueue = nullptr;
  MTL::Device* device = nullptr;
  
  void InitNgsMetal()
  {
    if (device) return;
    
    cout << "InitNgsMetal" << endl;
    
        // Top-level autorelease pool handles memory lifetime management
    NS::AutoreleasePool* pool = NS::AutoreleasePool::alloc()->init();

    // 2. Obtain default Metal device & Command Queue
    device = MTL::CreateSystemDefaultDevice();
    if (!device) {
        std::cerr << "Metal is not supported on this system.\n";
        pool->release();
        return;
    }
    
    commandQueue = device->newCommandQueue();
  }

  MTL::CommandQueue* GetCommandQueue()
  {
    InitNgsMetal();
    return commandQueue;
  }
  
  MTL::Device* GetDevice()
  {
    InitNgsMetal();
    return device;
  }
  
}
