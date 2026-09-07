/*********************************************************************/
/* File:   metal_device.cpp                                          */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   29. Aug. 2026                                             */
/*********************************************************************/

#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION
#include <Foundation/Foundation.hpp>
#include <Metal/Metal.hpp>

#include <iostream>
#include <stdexcept>
#include <string>

#include "metal_device.hpp"
#include <IOKit/IOKitLib.h>

namespace ngsmetal
{
  using namespace ngs_gpu;

  static MTL::Device * device = nullptr;
  static MTL::CommandQueue * commandQueue = nullptr;

  MTL::Device * GetDevice()
  {
    if (!device)
      {
        device = MTL::CreateSystemDefaultDevice();
        if (device)
          commandQueue = device->newCommandQueue();
        else
          std::cerr << "Metal is not supported on this system.\n";
      }
    return device;
  }

  MTL::CommandQueue * GetCommandQueue()
  {
    GetDevice();
    return commandQueue;
  }

  static void Err (const std::string & msg)
  { throw std::runtime_error ("ngsmetal: " + msg); }

  static NS::String * Str (const std::string & s)
  { return NS::String::string (s.c_str(), NS::UTF8StringEncoding); }


  class MetalBuffer : public Buffer
  {
    MTL::Device * dev;
    MTL::CommandQueue * queue;
    MTL::Buffer * buf;

  public:
    MetalBuffer (MTL::Device * adev, MTL::CommandQueue * aqueue,
                 size_t bytes, MemType mt)
      : Buffer(bytes, mt), dev(adev), queue(aqueue)
    {
      auto mode = (mt == MemType::Shared)
        ? MTL::ResourceStorageModeShared : MTL::ResourceStorageModePrivate;
      buf = dev->newBuffer (bytes, mode);
      if (!buf) Err ("newBuffer failed");
    }

    ~MetalBuffer() { buf->release(); }

    MTL::Buffer * Get() const { return buf; }

  protected:
    void * DoHostPtr() const override
    { return (memtype == MemType::Shared) ? buf->contents() : nullptr; }

    void DoH2D (const void * src, size_t bytes, size_t offset) override
    {
      if (offset+bytes > size) Err ("H2D out of range");
      if (memtype == MemType::Shared)
        {
          std::memcpy ((char*)buf->contents()+offset, src, bytes);
          return;
        }
      auto stage = dev->newBuffer (src, bytes, MTL::ResourceStorageModeShared);
      Blit (stage, 0, buf, offset, bytes);
      stage->release();
    }

    void DoD2H (void * dst, size_t bytes, size_t offset) const override
    {
      if (offset+bytes > size) Err ("D2H out of range");
      if (memtype == MemType::Shared)
        {
          std::memcpy (dst, (const char*)buf->contents()+offset, bytes);
          return;
        }
      auto stage = dev->newBuffer (bytes, MTL::ResourceStorageModeShared);
      Blit (buf, offset, stage, 0, bytes);
      std::memcpy (dst, stage->contents(), bytes);
      stage->release();
    }

  private:
    // private storage needs a blit, runs synchronously
    void Blit (MTL::Buffer * from, size_t foff, MTL::Buffer * to, size_t toff,
               size_t bytes) const
    {
      auto cb = queue->commandBuffer();
      auto enc = cb->blitCommandEncoder();
      enc->copyFromBuffer (from, foff, to, toff, bytes);
      enc->endEncoding();
      cb->commit();
      cb->waitUntilCompleted();
    }
  };


  class MetalKernel : public Kernel
  {
    std::string name;
    MTL::ComputePipelineState * pso;
  public:
    MetalKernel (const std::string & aname, MTL::ComputePipelineState * apso)
      : name(aname), pso(apso) { }
    ~MetalKernel() { pso->release(); }

    string Name() const override { return name; }
    MTL::ComputePipelineState * Get() const { return pso; }

    KernelInfo Info (size_t groupsize) const override
    {
      KernelInfo ki;
      ki.shared_bytes = pso->staticThreadgroupMemoryLength();
      ki.max_threads_per_group = pso->maxTotalThreadsPerThreadgroup();
      return ki;
    }
  };


  class MetalLibrary : public Library
  {
    MTL::Device * dev;
    MTL::Library * lib;
  public:
    MetalLibrary (MTL::Device * adev, MTL::Library * alib) : dev(adev), lib(alib) { }
    ~MetalLibrary() { lib->release(); }

    shared_ptr<Kernel> DoGetKernel (const string & name) override
    {
      auto func = lib->newFunction (Str(name));
      if (!func) Err ("no kernel '" + name + "' in library");

      NS::Error * error = nullptr;
      auto pso = dev->newComputePipelineState (func, &error);
      func->release();
      if (!pso)
        Err ("pipeline state for '" + name + "': " +
             (error ? error->localizedDescription()->utf8String() : "unknown"));
      return std::make_shared<MetalKernel> (name, pso);
    }
  };


  class MetalQueue : public Queue
  {
    MTL::CommandQueue * queue;
    mutable MTL::CommandBuffer * pending = nullptr;

  public:
    MetalQueue (MTL::CommandQueue * aqueue) : queue(aqueue) { }
    ~MetalQueue() { if (pending) pending->release(); }

    void DoFinish() override
    {
      if (!pending) return;
      pending->waitUntilCompleted();

      // metal reports kernel faults asynchronously, only here
      std::string msg;
      if (auto err = pending->error())
        msg = err->localizedDescription()->utf8String();
      bool failed = (pending->status() == MTL::CommandBufferStatusError);

      pending->release();
      pending = nullptr;

      if (failed || !msg.empty())
        Err ("kernel execution failed: " + (msg.empty() ? "unknown" : msg));
    }

  private:
    void Encode (MTL::ComputeCommandEncoder * enc, Kernel & kernel, Dim3 groups, Dim3 groupsize,
                 const std::vector<KernelArg> & args, size_t dynamic_group_memory)
    {
      auto & mk = dynamic_cast<MetalKernel&> (kernel);
      enc->setComputePipelineState (mk.Get());

      for (size_t i = 0; i < args.size(); i++)
        {
          auto & a = args[i];
          if (a.GetKind() == KernelArg::Kind::Buffer)
            {
              auto mb = dynamic_cast<MetalBuffer*> (a.GetBuffer());
              if (!mb) Err ("kernel argument is not a metal buffer");
              enc->setBuffer (mb->Get(), a.Offset(), i);
            }
          else
            enc->setBytes (a.Data(), a.NBytes(), i);
        }

      if (dynamic_group_memory)
        enc->setThreadgroupMemoryLength (dynamic_group_memory, 0);

      enc->dispatchThreadgroups (MTL::Size(groups.x, groups.y, groups.z),
                                 MTL::Size(groupsize.x, groupsize.y, groupsize.z));
    }

    void Commit (MTL::CommandBuffer * cb)
    {
      if (pending) pending->release();
      pending = cb;
      pending->retain();
      pending->commit();
    }

  protected:
    void DoLaunch (Kernel & kernel, Dim3 groups, Dim3 groupsize,
                   const std::vector<KernelArg> & args,
                   size_t dynamic_group_memory) override
    {
      auto cb = queue->commandBuffer();
      auto enc = cb->computeCommandEncoder();
      Encode (enc, kernel, groups, groupsize, args, dynamic_group_memory);
      enc->endEncoding();
      Commit (cb);
    }

    // one command buffer, one serial encoder: dispatches run in order
    void DoReplay (const Program & prog) override
    {
      auto cb = queue->commandBuffer();
      auto enc = cb->computeCommandEncoder();
      for (auto & n : prog.Nodes())
        Encode (enc, *n.kernel, n.groups, n.groupsize, n.args, n.dynamic_group_memory);
      enc->endEncoding();
      Commit (cb);
    }
  };


  class MetalDevice : public Device
  {
    MTL::Device * dev;
    MTL::CommandQueue * queue;
    shared_ptr<Queue> defqueue;

  public:
    MetalDevice (MTL::Device * adev, MTL::CommandQueue * aqueue)
      : dev(adev), queue(aqueue),
        defqueue(std::make_shared<MetalQueue>(aqueue)) { }

    string Name() const override { return dev->name()->utf8String(); }
    bool HasFloat64() const override { return false; }       // no fp64 in Metal
    bool IsUnifiedMemory() const override { return dev->hasUnifiedMemory(); }
    size_t MaxThreadsPerGroup() const override
    { return dev->maxThreadsPerThreadgroup().width; }
    size_t SimdWidth() const override { return 32; }
    // gpu cores from the IOKit registry (Metal has no API for it)
    size_t ComputeUnits() const override
    {
      static size_t cores = [] () -> size_t
      {
        size_t n = 8;
        io_service_t srv = IOServiceGetMatchingService (kIOMainPortDefault, IOServiceMatching ("AGXAccelerator"));
        if (srv)
          {
            if (CFTypeRef val = IORegistryEntryCreateCFProperty (srv, CFSTR("gpu-core-count"), kCFAllocatorDefault, 0))
              {
                int v = 0;
                if (CFGetTypeID(val) == CFNumberGetTypeID() && CFNumberGetValue ((CFNumberRef)val, kCFNumberIntType, &v) && v > 0)
                  n = v;
                CFRelease (val);
              }
            IOObjectRelease (srv);
          }
        return n;
      } ();
      return cores;
    }

    shared_ptr<Buffer> DoNewBuffer (size_t bytes, MemType mt) override
    { return std::make_shared<MetalBuffer> (dev, queue, bytes, mt); }

    shared_ptr<Library> DoCompileSource (const string & source) override
    {
      NS::Error * error = nullptr;
      auto lib = dev->newLibrary (Str(source), nullptr, &error);
      if (!lib)
        Err (string("shader compile error: ") +
             (error ? error->localizedDescription()->utf8String() : "unknown"));
      return std::make_shared<MetalLibrary> (dev, lib);
    }

    shared_ptr<Queue> DefaultQueue() override { return defqueue; }
  };


  MTL::Buffer * GetMTLBuffer (const ngs_gpu::Buffer & buf)
  {
    auto mb = dynamic_cast<const MetalBuffer*> (&buf);
    if (!mb) Err ("buffer was not created by the metal device");
    return mb->Get();
  }


  void InitMetalDevice()
  {
    ngs_gpu::SetDeviceCreator ([]() -> shared_ptr<ngs_gpu::Device>
    {
      auto dev = GetDevice();
      if (!dev) return nullptr;
      return std::make_shared<MetalDevice> (dev, GetCommandQueue());
    });
  }
}
