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

#include <mach/mach_time.h>

#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <core/paje_trace.hpp>

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


  /*
    Buffers are sub-allocated from heaps. A fresh MTL::Buffer costs a page
    fault per 16 kB on first touch (~20 GB/s, and it does not parallelise),
    a heap keeps its pages wired across allocations, so a buffer released and
    another one of any size allocated later fills at memory bandwidth.
    Heaps grow on demand in chunks (or the request size, if larger), empty
    heaps beyond one spare per storage mode are released.
  */
  class MetalHeapPool
  {
    MTL::Device * dev;
    std::mutex mutex;
    std::vector<MTL::Heap*> heaps[2];   // [0] private, [1] shared
    static constexpr size_t CHUNK = size_t(256) << 20;

    static bool IsShared (MTL::ResourceOptions mode)
    { return ((mode >> 4) & 0xF) == MTL::StorageModeShared; }
    static int Index (MTL::ResourceOptions mode) { return IsShared(mode) ? 1 : 0; }

  public:
    MetalHeapPool (MTL::Device * adev) : dev(adev) { }
    ~MetalHeapPool()
    {
      for (auto & hs : heaps)
        for (auto h : hs) h->release();
    }

    MTL::Buffer * NewBuffer (size_t bytes, MTL::ResourceOptions mode)
    {
      std::lock_guard<std::mutex> lock(mutex);
      auto & hs = heaps[Index(mode)];
      auto sa = dev->heapBufferSizeAndAlign (bytes, mode);
      for (auto h : hs)
        if (h->maxAvailableSize (sa.align) >= sa.size)
          if (auto b = h->newBuffer (bytes, mode)) return b;

      auto hd = MTL::HeapDescriptor::alloc()->init();
      hd->setSize (std::max (CHUNK, sa.size));
      hd->setStorageMode (IsShared(mode) ? MTL::StorageModeShared : MTL::StorageModePrivate);
      hd->setHazardTrackingMode (MTL::HazardTrackingModeTracked);
      auto heap = dev->newHeap (hd);
      hd->release();
      if (!heap) return dev->newBuffer (bytes, mode);   // out of heap space, plain buffer
      hs.push_back (heap);
      return heap->newBuffer (bytes, mode);
    }

    // release the buffer, drop empty heaps beyond one spare
    void Release (MTL::Buffer * buf)
    {
      std::lock_guard<std::mutex> lock(mutex);
      auto heap = buf->heap();
      auto & hs = heaps[Index(buf->resourceOptions())];
      buf->release();
      if (!heap) return;
      bool spare = false;
      for (size_t i = 0; i < hs.size(); )
        {
          if (hs[i]->usedSize() == 0)
            {
              if (spare) { hs[i]->release(); hs.erase (hs.begin()+i); continue; }
              spare = true;
            }
          i++;
        }
    }

    size_t HeldBytes()
    {
      std::lock_guard<std::mutex> lock(mutex);
      size_t n = 0;
      for (auto & hs : heaps) for (auto h : hs) n += h->size();
      return n;
    }
  };

  // a heap buffer released and returned to the pool when it goes out of scope
  struct PoolBuffer
  {
    MetalHeapPool & pool;
    MTL::Buffer * buf;
    PoolBuffer (MetalHeapPool & apool, size_t bytes, MTL::ResourceOptions mode)
      : pool(apool), buf(apool.NewBuffer (bytes, mode))
    { if (!buf) Err ("newBuffer failed"); }
    ~PoolBuffer() { pool.Release (buf); }
    MTL::Buffer * operator->() const { return buf; }
    operator MTL::Buffer*() const { return buf; }
  };


  class MetalQueue;

  class MetalBuffer : public Buffer
  {
    shared_ptr<MetalHeapPool> pool;
    shared_ptr<MetalQueue> queue;   // transfers run (and are traced) on it
    MTL::Buffer * buf;

  public:
    MetalBuffer (shared_ptr<MetalHeapPool> apool, shared_ptr<MetalQueue> aqueue,
                 size_t bytes, MemType mt)
      : Buffer(bytes, mt), pool(std::move(apool)), queue(std::move(aqueue))
    {
      auto mode = (mt == MemType::Shared)
        ? MTL::ResourceStorageModeShared : MTL::ResourceStorageModePrivate;
      buf = pool->NewBuffer (bytes, mode);
      if (!buf) Err ("newBuffer failed");
    }

    ~MetalBuffer() { pool->Release (buf); }

    MTL::Buffer * Get() const { return buf; }

  protected:
    void * DoHostPtr() const override
    { return (memtype == MemType::Shared) ? buf->contents() : nullptr; }

    void DoH2D (const void * src, size_t bytes, size_t offset) override;
    void DoD2H (void * dst, size_t bytes, size_t offset) const override;
    void DoFill (size_t bytes, size_t offset, const FillFunc & fill) override;

  private:
    // private storage needs a blit, runs synchronously
    void Blit (MTL::Buffer * from, size_t foff, MTL::Buffer * to, size_t toff,
               size_t bytes, const std::string & label) const;
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

    static constexpr size_t TRACE_CAPACITY = 4096;
    ngcore::TraceContainer tracer{"GPU metal"};
    struct Traced { std::string label; MTL::CommandBuffer * cb; int value; };
    std::vector<Traced> traced;

    // mach_absolute_time units per second
    static double MachPerSec()
    {
      static const double f = []
      { mach_timebase_info_data_t t; mach_timebase_info (&t);
        return 1e9 * t.denom / t.numer; } ();
      return f;
    }

    void FlushTrace()
    {
      if (traced.empty()) return;

      // GPUStartTime/GPUEndTime are mach_absolute_time in seconds. ngcore ticks
      // run at a different rate (1 GHz on M4, 24 MHz mach timebase), and the
      // startup calibration of seconds_per_tick is coarse, so map through the
      // rate measured since the first flush, anchored at the current reading
      static const ngcore::TTimePoint tick0 = ngcore::GetTimeCounter();
      static const unsigned long long mach0 = mach_absolute_time();
      ngcore::TTimePoint tick = ngcore::GetTimeCounter();
      unsigned long long mach = mach_absolute_time();
      double rate = (mach - mach0 > MachPerSec()/20)          // > 50 ms baseline
        ? double(tick - tick0) / double(mach - mach0)
        : 1.0 / (ngcore::seconds_per_tick * MachPerSec());

      auto ToTicks = [&] (double sec)
        {
          return ngcore::TTimePoint ((long long)tick +
                                     (long long)((sec*MachPerSec() - double(mach)) * rate));
        };

      for (auto & [label, cb, value] : traced)
        {
          tracer.AddTicks (label, ToTicks(cb->GPUStartTime()), ToTicks(cb->GPUEndTime()), value);
          cb->release();
        }
      traced.clear();
    }

  public:
    MetalQueue (MTL::CommandQueue * aqueue) : queue(aqueue) { }
    MTL::CommandQueue * Get() const { return queue; }
    ~MetalQueue()
    {
      for (auto & t : traced) t.cb->release();
      if (pending) pending->release();
    }

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

      // the queue is serial, so every traced buffer completed with the last
      FlushTrace();

      if (failed || !msg.empty())
        Err ("kernel execution failed: " + (msg.empty() ? "unknown" : msg));
    }

    // a transfer command buffer: committed, traced like a launch, waited for
    void Transfer (MTL::CommandBuffer * cb, const std::string & label, size_t bytes)
    {
      Commit (cb, label, TransferValue(bytes));
      DoFinish();
    }

    // host-side transfer (shared storage), on the same trace row
    void TraceHost (const std::string & label, size_t bytes,
                    ngcore::TTimePoint t0, ngcore::TTimePoint t1)
    { tracer.AddTicks (label, t0, t1, TransferValue(bytes)); }

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

    void Commit (MTL::CommandBuffer * cb, const std::string & label, int value = 0)
    {
      if (tracer.Active())
        {
          if (traced.size() == TRACE_CAPACITY) DoFinish();
          cb->retain();
          traced.push_back ({label, cb, value});
        }
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
      Commit (cb, kernel.Name());
    }

    // one command buffer, one serial encoder: dispatches run in order
    void DoReplay (const Program & prog) override
    {
      auto cb = queue->commandBuffer();
      auto enc = cb->computeCommandEncoder();
      for (auto & n : prog.Nodes())
        Encode (enc, *n.kernel, n.groups, n.groupsize, n.args, n.dynamic_group_memory);
      enc->endEncoding();
      Commit (cb, "Program (" + std::to_string(prog.Size()) + " launches)");
    }
  };


  void MetalBuffer :: DoH2D (const void * src, size_t bytes, size_t offset)
  {
    if (offset+bytes > size) Err ("H2D out of range");
    auto label = TransferLabel ("H2D", bytes);
    if (memtype == MemType::Shared)
      {
        auto t0 = ngcore::GetTimeCounter();
        std::memcpy ((char*)buf->contents()+offset, src, bytes);
        queue->TraceHost (label, bytes, t0, ngcore::GetTimeCounter());
        return;
      }
    PoolBuffer stage (*pool, bytes, MTL::ResourceStorageModeShared);
    std::memcpy (stage->contents(), src, bytes);
    Blit (stage, 0, buf, offset, bytes, label);
  }

  // private storage only (shared is filled in place by Buffer::H2D)
  void MetalBuffer :: DoFill (size_t bytes, size_t offset, const FillFunc & fill)
  {
    PoolBuffer stage (*pool, bytes, MTL::ResourceStorageModeShared);
    fill (stage->contents(), 0, bytes);
    Blit (stage, 0, buf, offset, bytes, TransferLabel ("H2D", bytes));
  }

  void MetalBuffer :: DoD2H (void * dst, size_t bytes, size_t offset) const
  {
    if (offset+bytes > size) Err ("D2H out of range");
    auto label = TransferLabel ("D2H", bytes);
    if (memtype == MemType::Shared)
      {
        auto t0 = ngcore::GetTimeCounter();
        std::memcpy (dst, (const char*)buf->contents()+offset, bytes);
        queue->TraceHost (label, bytes, t0, ngcore::GetTimeCounter());
        return;
      }
    PoolBuffer stage (*pool, bytes, MTL::ResourceStorageModeShared);
    Blit (buf, offset, stage, 0, bytes, label);
    std::memcpy (dst, stage->contents(), bytes);
  }

  void MetalBuffer :: Blit (MTL::Buffer * from, size_t foff, MTL::Buffer * to, size_t toff,
                            size_t bytes, const std::string & label) const
  {
    auto cb = queue->Get()->commandBuffer();
    auto enc = cb->blitCommandEncoder();
    enc->copyFromBuffer (from, foff, to, toff, bytes);
    enc->endEncoding();
    queue->Transfer (cb, label, bytes);
  }


  class MetalDevice : public Device
  {
    MTL::Device * dev;
    shared_ptr<MetalHeapPool> pool;
    shared_ptr<MetalQueue> defqueue;

  public:
    MetalDevice (MTL::Device * adev, MTL::CommandQueue * aqueue)
      : dev(adev), pool(std::make_shared<MetalHeapPool>(adev)),
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
    { return std::make_shared<MetalBuffer> (pool, defqueue, bytes, mt); }

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
