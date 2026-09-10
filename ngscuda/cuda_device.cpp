/*********************************************************************/
/* File:   cuda_device.cpp                                           */
/* Author: Joachim Schoeberl                                         */
/*         (developed with AI assistance, Claude Fable 5.1)          */
/* Date:   29. Aug. 2026                                             */
/*********************************************************************/

#include <cuda.h>
#include <nvrtc.h>

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <core/paje_trace.hpp>
#include <core/taskmanager.hpp>

#include "cuda_device.hpp"

namespace ngs_cuda
{
  using namespace ngs_gpu;

  static void Check (CUresult res, const std::string & what)
  {
    if (res == CUDA_SUCCESS) return;
    const char * msg = nullptr;
    cuGetErrorString (res, &msg);
    throw std::runtime_error ("ngscuda: " + what + ": " + (msg ? msg : "unknown"));
  }

  static void Check (nvrtcResult res, const std::string & what)
  {
    if (res == NVRTC_SUCCESS) return;
    throw std::runtime_error ("ngscuda: " + what + ": " + nvrtcGetErrorString(res));
  }


  class CudaQueue;

  class CudaBuffer : public Buffer
  {
    CUdeviceptr ptr;
    shared_ptr<CudaQueue> queue;   // transfers run (and are traced) on it
  public:
    CudaBuffer (size_t bytes, MemType mt, shared_ptr<CudaQueue> aqueue)
      : Buffer(bytes, mt), queue(std::move(aqueue))
    {
      if (mt == MemType::Shared)
        Check (cuMemAllocManaged (&ptr, bytes, CU_MEM_ATTACH_GLOBAL), "cuMemAllocManaged");
      else
        Check (cuMemAlloc (&ptr, bytes), "cuMemAlloc");
    }

    ~CudaBuffer() { cuMemFree (ptr); }

    CUdeviceptr Get() const { return ptr; }

  protected:
    void * DoHostPtr() const override
    {
      // managed memory is addressable from the host
      return (memtype == MemType::Shared) ? (void*)ptr : nullptr;
    }

    uintptr_t DoDevicePtr() const override { return uintptr_t(ptr); }

    void DoH2D (const void * src, size_t bytes, size_t offset) override;
    void DoD2H (void * dst, size_t bytes, size_t offset) const override;
    void DoFill (size_t bytes, size_t offset, const FillFunc & fill) override;
  };


  class CudaKernel : public Kernel
  {
    std::string name;
    CUfunction func;
  public:
    CudaKernel (const std::string & aname, CUfunction afunc)
      : name(aname), func(afunc) { }

    string Name() const override { return name; }
    CUfunction Get() const { return func; }

    KernelInfo Info (size_t groupsize) const override
    {
      KernelInfo ki;
      auto attr = [&] (CUfunction_attribute a) { int v = 0; cuFuncGetAttribute (&v, a, func); return size_t(v); };
      ki.registers = attr (CU_FUNC_ATTRIBUTE_NUM_REGS);
      ki.local_bytes = attr (CU_FUNC_ATTRIBUTE_LOCAL_SIZE_BYTES);
      ki.shared_bytes = attr (CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES);
      ki.max_threads_per_group = attr (CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK);
      if (groupsize)
        {
          int blocks = 0;
          cuOccupancyMaxActiveBlocksPerMultiprocessor (&blocks, func, int(groupsize), 0);
          ki.max_groups_per_unit = blocks;
        }
      return ki;
    }
  };


  class CudaLibrary : public Library
  {
    CUmodule module;
  public:
    CudaLibrary (CUmodule amodule) : module(amodule) { }
    ~CudaLibrary() { cuModuleUnload (module); }

    shared_ptr<Kernel> DoGetKernel (const string & name) override
    {
      CUfunction func;
      auto res = cuModuleGetFunction (&func, module, name.c_str());
      if (res != CUDA_SUCCESS)
        throw std::runtime_error
          ("ngscuda: no kernel '" + name + "' in module"
           " (entry points must be declared extern \"C\" __global__)");
      return std::make_shared<CudaKernel> (name, func);
    }
  };


  // the stream the runtime-API side of ngscuda launches on; the graph
  // capture machinery (cuda_core.hpp) redirects it temporarily.
  // cudaStream_t and CUstream are the same underlying type.
  extern CUstream ngs_cuda_stream;

  class CudaQueue : public Queue
  {
    CUstream stream = nullptr;
    bool owned = false;
    bool tracking = false;    // follow ngs_cuda_stream at launch time
    CUstream capture_stream = nullptr;   // graphs are captured here, see DoReplay
    CUstream forced = nullptr;           // launches go here while set

    CUstream Current() const { return forced ? forced : (tracking ? ngs_cuda_stream : stream); }

    /*
      Transfers from pageable memory are staged by the driver through a
      pinned buffer with a single-threaded memcpy (~10 GB/s here). Own
      pinned ring: chunks filled with a parallel memcpy while the previous
      chunk is on the bus.
    */
    static constexpr size_t STAGE_CHUNK = size_t(32) << 20;
    static constexpr int STAGE_N = 2;
    static constexpr size_t STAGE_MIN = size_t(1) << 20;   // below: driver path
    void * stage[STAGE_N] = { };
    CUevent stage_event[STAGE_N] = { };

    void InitStaging()
    {
      if (stage[0]) return;
      for (int i = 0; i < STAGE_N; i++)
        {
          Check (cuMemHostAlloc (&stage[i], STAGE_CHUNK, CU_MEMHOSTALLOC_PORTABLE), "cuMemHostAlloc");
          Check (cuEventCreate (&stage_event[i], CU_EVENT_DISABLE_TIMING), "cuEventCreate");
          Check (cuEventRecord (stage_event[i], Current()), "cuEventRecord");
        }
    }

    static void ParallelMemcpy (void * dst, const void * src, size_t bytes)
    {
      ngcore::ParallelForRange (bytes, [&] (ngcore::IntRange r)
        { std::memcpy ((char*)dst+r.First(), (const char*)src+r.First(), r.Size()); });
    }

    // produce(chunk, off, n) writes bytes [off, off+n) into the pinned chunk
    void UploadStaged (CUdeviceptr dst, size_t bytes, const Buffer::FillFunc & produce)
    {
      InitStaging();
      int k = 0;
      for (size_t off = 0; off < bytes; off += STAGE_CHUNK, k = (k+1) % STAGE_N)
        {
          size_t n = std::min (STAGE_CHUNK, bytes-off);
          Check (cuEventSynchronize (stage_event[k]), "cuEventSynchronize");   // chunk free again
          produce (stage[k], off, n);
          Check (cuMemcpyHtoDAsync (dst+off, stage[k], n, Current()), "cuMemcpyHtoDAsync");
          Check (cuEventRecord (stage_event[k], Current()), "cuEventRecord");
        }
    }

    void DownloadStaged (void * dst, CUdeviceptr src, size_t bytes)
    {
      InitStaging();
      int k = 0;
      size_t prev_off = 0, prev_n = 0; int prev_k = -1;
      for (size_t off = 0; off < bytes; off += STAGE_CHUNK, k = (k+1) % STAGE_N)
        {
          size_t n = std::min (STAGE_CHUNK, bytes-off);
          Check (cuEventSynchronize (stage_event[k]), "cuEventSynchronize");
          Check (cuMemcpyDtoHAsync (stage[k], src+off, n, Current()), "cuMemcpyDtoHAsync");
          Check (cuEventRecord (stage_event[k], Current()), "cuEventRecord");
          if (prev_k >= 0)     // copy the previous chunk out while this one is on the bus
            {
              Check (cuEventSynchronize (stage_event[prev_k]), "cuEventSynchronize");
              ParallelMemcpy ((char*)dst+prev_off, stage[prev_k], prev_n);
            }
          prev_off = off; prev_n = n; prev_k = k;
        }
      Check (cuEventSynchronize (stage_event[prev_k]), "cuEventSynchronize");
      ParallelMemcpy ((char*)dst+prev_off, stage[prev_k], prev_n);
    }

    static constexpr size_t TRACE_CAPACITY = 4096;
    ngcore::TraceContainer tracer{"GPU cuda"};
    std::vector<CUevent> trace_events;     // start/stop pair per slot
    std::vector<std::string> trace_labels;
    std::vector<int> trace_values;
    size_t trace_slots = 0;
    CUevent trace_anchor = nullptr;

    // records a start event on the stream, returns the slot, -1 if not traced
    int BeginTrace (const std::string & label, int value = 0)
    {
      if (!tracer.Active()) return -1;
      // events recorded into a graph carry no readable time
      CUstreamCaptureStatus capturing;
      if (cuStreamIsCapturing (Current(), &capturing) != CUDA_SUCCESS ||
          capturing != CU_STREAM_CAPTURE_STATUS_NONE) return -1;

      if (trace_events.empty())
        {
          trace_events.resize (2*TRACE_CAPACITY);
          trace_labels.resize (TRACE_CAPACITY);
          trace_values.resize (TRACE_CAPACITY);
          for (auto & ev : trace_events)
            Check (cuEventCreate (&ev, CU_EVENT_DEFAULT), "cuEventCreate");
          Check (cuEventCreate (&trace_anchor, CU_EVENT_DEFAULT), "cuEventCreate");
        }
      if (trace_slots == TRACE_CAPACITY) FlushTrace();
      Check (cuEventRecord (trace_events[2*trace_slots], Current()), "cuEventRecord");
      trace_labels[trace_slots] = label;
      trace_values[trace_slots] = value;
      return int(trace_slots);
    }

    void EndTrace (int slot)
    {
      if (slot < 0) return;
      Check (cuEventRecord (trace_events[2*slot+1], Current()), "cuEventRecord");
      trace_slots++;
    }

    void FlushTrace()
    {
      if (!trace_slots) return;
      Check (cuEventRecord (trace_anchor, Current()), "cuEventRecord");
      Check (cuEventSynchronize (trace_anchor), "cuEventSynchronize");
      tracer.Anchor (0);
      for (size_t i = 0; i < trace_slots; i++)
        {
          float t0 = 0, t1 = 0;
          cuEventElapsedTime (&t0, trace_events[2*i], trace_anchor);
          cuEventElapsedTime (&t1, trace_events[2*i+1], trace_anchor);
          tracer.AddInterval (trace_labels[i], -1e-3*t0, -1e-3*t1, trace_values[i]);
        }
      trace_slots = 0;
    }

  public:
    CudaQueue() : owned(true)
    { Check (cuStreamCreate (&stream, CU_STREAM_NON_BLOCKING), "cuStreamCreate"); }

    // a synchronous transfer on the queue's stream, traced like a launch
    template <typename F>
    void Transfer (const std::string & label, size_t bytes, F copy)
    {
      int slot = BeginTrace (label, TransferValue(bytes));
      copy (Current());
      EndTrace (slot);
      DoFinish();   // sync + flush, the trace may end without a Finish
    }

    // host-side transfer (managed memory), on the same trace row
    void TraceHost (const std::string & label, size_t bytes,
                    ngcore::TTimePoint t0, ngcore::TTimePoint t1)
    { tracer.AddTicks (label, t0, t1, TransferValue(bytes)); }

    void H2D (CUdeviceptr dst, const void * src, size_t bytes)
    {
      Transfer (TransferLabel ("H2D", bytes), bytes, [&] (CUstream s)
        {
          if (bytes < STAGE_MIN)
            Check (cuMemcpyHtoDAsync (dst, src, bytes, s), "cuMemcpyHtoDAsync");
          else
            UploadStaged (dst, bytes, [&] (void * chunk, size_t off, size_t n)
                          { ParallelMemcpy (chunk, (const char*)src+off, n); });
        });
    }

    void Fill (CUdeviceptr dst, size_t bytes, const Buffer::FillFunc & fill)
    {
      Transfer (TransferLabel ("H2D", bytes), bytes, [&] (CUstream s)
        {
          if (bytes < STAGE_MIN)
            {
              std::vector<char> tmp (bytes);
              fill (tmp.data(), 0, bytes);
              Check (cuMemcpyHtoDAsync (dst, tmp.data(), bytes, s), "cuMemcpyHtoDAsync");
              Check (cuStreamSynchronize (s), "cuStreamSynchronize");   // tmp dies here
            }
          else
            UploadStaged (dst, bytes, fill);
        });
    }

    void D2H (void * dst, CUdeviceptr src, size_t bytes)
    {
      Transfer (TransferLabel ("D2H", bytes), bytes, [&] (CUstream s)
        {
          if (bytes < STAGE_MIN)
            Check (cuMemcpyDtoHAsync (dst, src, bytes, s), "cuMemcpyDtoHAsync");
          else
            DownloadStaged (dst, src, bytes);
        });
    }

    // launches follow the current ngs_cuda_stream, so they stay ordered
    // with cuda libraries and are recorded during graph capture
    struct TrackNgsStream { };
    CudaQueue (TrackNgsStream) : tracking(true) { }

    ~CudaQueue()
    {
      for (int i = 0; i < STAGE_N; i++)
        if (stage[i]) { cuMemFreeHost (stage[i]); cuEventDestroy (stage_event[i]); }
      for (auto ev : trace_events) cuEventDestroy (ev);
      if (trace_anchor) cuEventDestroy (trace_anchor);
      if (owned) cuStreamDestroy (stream);
      if (capture_stream) cuStreamDestroy (capture_stream);
    }

    void DoFinish() override
    {
      Check (cuStreamSynchronize (Current()), "cuStreamSynchronize");
      FlushTrace();
    }

  protected:
    void DoLaunch (Kernel & kernel, Dim3 groups, Dim3 groupsize,
                   const std::vector<KernelArg> & args,
                   size_t dynamic_group_memory) override
    {
      auto & ck = dynamic_cast<CudaKernel&> (kernel);

      // cuLaunchKernel wants pointers to the argument values,
      // for a buffer that is the device pointer itself
      std::vector<CUdeviceptr> devptrs (args.size());
      std::vector<void*> params (args.size());

      for (size_t i = 0; i < args.size(); i++)
        {
          auto & a = args[i];
          if (a.GetKind() == KernelArg::Kind::Buffer)
            {
              auto cb = dynamic_cast<CudaBuffer*> (a.GetBuffer());
              if (!cb) throw std::runtime_error ("ngscuda: kernel argument is not a cuda buffer");
              devptrs[i] = cb->Get() + a.Offset();
              params[i] = &devptrs[i];
            }
          else
            params[i] = const_cast<void*> (a.Data());
        }

      int slot = BeginTrace (ck.Name());
      Check (cuLaunchKernel (ck.Get(),
                             groups.x, groups.y, groups.z,
                             groupsize.x, groupsize.y, groupsize.z,
                             dynamic_group_memory, Current(),
                             params.data(), nullptr), "cuLaunchKernel");
      EndTrace (slot);
    }

    struct GraphCache
    {
      CUgraph graph = nullptr;
      CUgraphExec exec = nullptr;
      ~GraphCache()
      {
        if (exec) cuGraphExecDestroy (exec);
        if (graph) cuGraphDestroy (graph);
      }
    };

    // the launches are captured into a graph on the first replay. Capture
    // needs a non-default stream (the queue may follow the legacy default
    // stream), nothing runs during capture; the graph is then launched on
    // the queue's stream
    void DoReplay (const Program & prog) override
    {
      auto cache = std::static_pointer_cast<GraphCache> (prog.backend_cache);
      if (!cache)
        {
          if (!capture_stream)
            Check (cuStreamCreate (&capture_stream, CU_STREAM_NON_BLOCKING), "cuStreamCreate");
          cache = std::make_shared<GraphCache>();
          forced = capture_stream;
          try
            {
              Check (cuStreamBeginCapture (capture_stream, CU_STREAM_CAPTURE_MODE_THREAD_LOCAL), "cuStreamBeginCapture");
              for (auto & n : prog.Nodes())
                DoLaunch (*n.kernel, n.groups, n.groupsize, n.args, n.dynamic_group_memory);
              Check (cuStreamEndCapture (capture_stream, &cache->graph), "cuStreamEndCapture");
            }
          catch (...) { forced = nullptr; throw; }
          forced = nullptr;
          Check (cuGraphInstantiate (&cache->exec, cache->graph, 0), "cuGraphInstantiate");
          prog.backend_cache = cache;
        }
      int slot = BeginTrace ("Program (" + std::to_string(prog.Size()) + " launches)");
      Check (cuGraphLaunch (cache->exec, Current()), "cuGraphLaunch");
      EndTrace (slot);
    }
  };


  void CudaBuffer :: DoH2D (const void * src, size_t bytes, size_t offset)
  {
    if (offset+bytes > size)
      throw std::runtime_error ("ngscuda: H2D out of range");
    auto label = TransferLabel ("H2D", bytes);
    if (memtype == MemType::Shared)
      {
        auto t0 = ngcore::GetTimeCounter();
        std::memcpy ((char*)ptr+offset, src, bytes);
        queue->TraceHost (label, bytes, t0, ngcore::GetTimeCounter());
      }
    else
      queue->H2D (ptr+offset, src, bytes);
  }

  // device memory only (managed is filled in place by Buffer::H2D)
  void CudaBuffer :: DoFill (size_t bytes, size_t offset, const FillFunc & fill)
  {
    if (offset+bytes > size)
      throw std::runtime_error ("ngscuda: H2D out of range");
    queue->Fill (ptr+offset, bytes, fill);
  }

  void CudaBuffer :: DoD2H (void * dst, size_t bytes, size_t offset) const
  {
    if (offset+bytes > size)
      throw std::runtime_error ("ngscuda: D2H out of range");
    auto label = TransferLabel ("D2H", bytes);
    if (memtype == MemType::Shared)
      {
        auto t0 = ngcore::GetTimeCounter();
        std::memcpy (dst, (const char*)ptr+offset, bytes);
        queue->TraceHost (label, bytes, t0, ngcore::GetTimeCounter());
      }
    else
      queue->D2H (dst, ptr+offset, bytes);
  }


  // if NGS_CUDA_DUMP_PTX is set, the generated ptx and a summary of the
  // compiled kernels are written to that directory ("1" -> current dir)
  static const char * PtxDumpDir()
  {
    const char * dir = getenv ("NGS_CUDA_DUMP_PTX");
    if (!dir || !*dir) return nullptr;
    return (strcmp(dir, "1") == 0) ? "." : dir;
  }

  // names of the '.entry' points, in the order they appear in the ptx
  static std::vector<std::string> PtxEntryNames (const std::string & ptx)
  {
    std::vector<std::string> names;
    for (size_t pos = 0; (pos = ptx.find (".entry", pos)) != std::string::npos; )
      {
        pos += 6;
        size_t beg = ptx.find_first_not_of (" \t\n", pos);
        if (beg == std::string::npos) break;
        size_t end = ptx.find_first_of (" \t\n(", beg);
        names.push_back (ptx.substr (beg, end-beg));
        pos = beg;
      }
    return names;
  }

  static void DumpPtx (const std::string & dir, const std::string & source,
                       const std::string & ptx, const std::string & log,
                       const std::string & arch, size_t cubinsize, CUmodule module)
  {
    static int counter = 0;
    std::string base = dir + "/ngsgpu_" + std::to_string(counter++);

    std::ofstream (base+".ptx") << ptx;

    std::ofstream info (base+".info");
    info << "source:      ngsgpu.cu, " << source.size() << " bytes\n"
         << "options:     " << arch << " --std=c++17 --device-as-default-execution-space\n"
         << "ptx:         " << ptx.size() << " bytes\n"
         << "cubin:       " << cubinsize << " bytes\n";
    if (log.find_first_not_of (" \t\n") != std::string::npos)
      info << "compile log:\n" << log << "\n";

    for (auto & name : PtxEntryNames (ptx))
      {
        CUfunction func;
        if (cuModuleGetFunction (&func, module, name.c_str()) != CUDA_SUCCESS)
          continue;
        auto attr = [&] (CUfunction_attribute a)
        { int v = 0; cuFuncGetAttribute (&v, a, func); return v; };

        info << "\nkernel " << name << "\n"
             << "  registers:              " << attr (CU_FUNC_ATTRIBUTE_NUM_REGS) << "\n"
             << "  local bytes:            " << attr (CU_FUNC_ATTRIBUTE_LOCAL_SIZE_BYTES) << "\n"
             << "  const bytes:            " << attr (CU_FUNC_ATTRIBUTE_CONST_SIZE_BYTES) << "\n"
             << "  static shared bytes:    " << attr (CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES) << "\n"
             << "  max threads per block:  " << attr (CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK) << "\n"
             << "  ptx version:            " << attr (CU_FUNC_ATTRIBUTE_PTX_VERSION) << "\n"
             << "  binary version:         " << attr (CU_FUNC_ATTRIBUTE_BINARY_VERSION) << "\n";

        for (int bs : { 32, 64, 128, 256, 512, 1024 })
          {
            int blocks = 0;
            if (cuOccupancyMaxActiveBlocksPerMultiprocessor (&blocks, func, bs, 0) == CUDA_SUCCESS)
              info << "  blocks/SM at " << bs << " threads: " << blocks << "\n";
          }
      }
  }


  class CudaDevice : public Device
  {
    CUdevice dev;
    shared_ptr<CudaQueue> defqueue;
    int ccmajor, ccminor;

    int Attr (CUdevice_attribute a) const
    { int v = 0; cuDeviceGetAttribute (&v, a, dev); return v; }

  public:
    CudaDevice (CUdevice adev) : dev(adev)
    {
      ccmajor = Attr (CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR);
      ccminor = Attr (CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR);
      defqueue = std::make_shared<CudaQueue> (CudaQueue::TrackNgsStream{});
    }

    string Name() const override
    {
      char buf[256] = "";
      cuDeviceGetName (buf, sizeof(buf), dev);
      return buf;
    }

    bool HasFloat64() const override { return true; }
    bool IsUnifiedMemory() const override
    {
      return Attr (CU_DEVICE_ATTRIBUTE_INTEGRATED) != 0
        || Attr (CU_DEVICE_ATTRIBUTE_PAGEABLE_MEMORY_ACCESS_USES_HOST_PAGE_TABLES) != 0
        || Attr (CU_DEVICE_ATTRIBUTE_DIRECT_MANAGED_MEM_ACCESS_FROM_HOST) != 0;
    }
    size_t MaxThreadsPerGroup() const override
    { return Attr (CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK); }
    size_t SimdWidth() const override
    { return Attr (CU_DEVICE_ATTRIBUTE_WARP_SIZE); }
    size_t ComputeUnits() const override
    { return Attr (CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT); }

    shared_ptr<Buffer> DoNewBuffer (size_t bytes, MemType mt) override
    { return std::make_shared<CudaBuffer> (bytes, mt, defqueue); }

    shared_ptr<Queue> DefaultQueue() override { return defqueue; }

    shared_ptr<Library> DoCompileSource (const string & source) override
    {
      nvrtcProgram prog;
      Check (nvrtcCreateProgram (&prog, source.c_str(), "ngsgpu.cu",
                                 0, nullptr, nullptr), "nvrtcCreateProgram");

      // sm_XY: a cubin for the exact device, so loading needs no
      // PTX->SASS JIT (which lazy module loading would otherwise pay
      // at the first launch of every kernel)
      std::string arch = "--gpu-architecture=sm_"
        + std::to_string(ccmajor) + std::to_string(ccminor);
      const char * opts[] = { arch.c_str(), "--std=c++17",
                              "--device-as-default-execution-space" };

      auto res = nvrtcCompileProgram (prog, 3, opts);

      size_t logsize = 0;
      nvrtcGetProgramLogSize (prog, &logsize);
      std::string log (logsize, '\0');
      if (logsize) nvrtcGetProgramLog (prog, log.data());
      while (!log.empty() && log.back() == '\0') log.pop_back();

      if (res != NVRTC_SUCCESS)
        {
          nvrtcDestroyProgram (&prog);
          // the error lines first: nvrtc logs may start with pages of warnings
          std::string errors;
          for (size_t pos = 0; pos < log.size(); )
            {
              size_t nl = log.find ('\n', pos);
              std::string line = log.substr (pos, nl == std::string::npos ? std::string::npos : nl-pos);
              if (line.find (" error") != std::string::npos) errors += line + "\n";
              if (nl == std::string::npos) break;
              pos = nl+1;
            }
          throw std::runtime_error (std::string("ngscuda: kernel compile error (") + nvrtcGetErrorString(res) + ")\n"
                                    + errors + "full nvrtc log:\n" + log);
        }

      const char * dumpdir = PtxDumpDir();

      std::string ptx;
      if (dumpdir)
        {
          size_t ptxsize = 0;
          Check (nvrtcGetPTXSize (prog, &ptxsize), "nvrtcGetPTXSize");
          ptx.resize (ptxsize);
          Check (nvrtcGetPTX (prog, ptx.data()), "nvrtcGetPTX");
          while (!ptx.empty() && ptx.back() == '\0') ptx.pop_back();
        }

      std::string image;
      size_t cubinsize = 0;
      if (nvrtcGetCUBINSize (prog, &cubinsize) == NVRTC_SUCCESS && cubinsize > 0)
        {
          image.resize (cubinsize);
          Check (nvrtcGetCUBIN (prog, image.data()), "nvrtcGetCUBIN");
        }
      else
        {
          cubinsize = 0;
          if (ptx.empty())
            {
              size_t ptxsize = 0;
              Check (nvrtcGetPTXSize (prog, &ptxsize), "nvrtcGetPTXSize");
              ptx.resize (ptxsize);
              Check (nvrtcGetPTX (prog, ptx.data()), "nvrtcGetPTX");
            }
          image = ptx;
        }
      nvrtcDestroyProgram (&prog);

      CUmodule module;
      Check (cuModuleLoadData (&module, image.data()), "cuModuleLoadData");

      if (dumpdir)
        DumpPtx (dumpdir, source, ptx, log, arch, cubinsize, module);
      return std::make_shared<CudaLibrary> (module);
    }
  };


  void * BufferDevPtr (ngs_gpu::Buffer & buf)
  {
    auto p = buf.DevicePtr();
    if (!p) throw std::runtime_error ("ngscuda: buffer is not a cuda buffer");
    return (void*)p;
  }


  void InitCudaDevice()
  {
    ngs_gpu::SetDeviceCreator ([]() -> shared_ptr<ngs_gpu::Device>
    {
      Check (cuInit(0), "cuInit");

      CUdevice dev;
      Check (cuDeviceGet (&dev, 0), "cuDeviceGet");

      // the primary context is the one the runtime API uses,
      // so this stays interoperable with the rest of ngscuda
      CUcontext ctx;
      Check (cuDevicePrimaryCtxRetain (&ctx, dev), "cuDevicePrimaryCtxRetain");
      Check (cuCtxSetCurrent (ctx), "cuCtxSetCurrent");

      return std::make_shared<CudaDevice> (dev);
    });
  }
}
