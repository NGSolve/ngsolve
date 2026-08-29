/*********************************************************************/
/* File:   cuda_device.cpp                                           */
/* Date:   29. Aug. 2025                                             */
/*********************************************************************/

#include <cuda.h>
#include <nvrtc.h>

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

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


  class CudaBuffer : public Buffer
  {
    CUdeviceptr ptr;
  public:
    CudaBuffer (size_t bytes, MemType mt) : Buffer(bytes, mt)
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

    void DoH2D (const void * src, size_t bytes, size_t offset) override
    {
      if (offset+bytes > size)
        throw std::runtime_error ("ngscuda: H2D out of range");
      if (memtype == MemType::Shared)
        std::memcpy ((char*)ptr+offset, src, bytes);
      else
        Check (cuMemcpyHtoD (ptr+offset, src, bytes), "cuMemcpyHtoD");
    }

    void DoD2H (void * dst, size_t bytes, size_t offset) const override
    {
      if (offset+bytes > size)
        throw std::runtime_error ("ngscuda: D2H out of range");
      if (memtype == MemType::Shared)
        std::memcpy (dst, (const char*)ptr+offset, bytes);
      else
        Check (cuMemcpyDtoH (dst, ptr+offset, bytes), "cuMemcpyDtoH");
    }
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
  };


  class CudaLibrary : public Library
  {
    CUmodule module;
  public:
    CudaLibrary (CUmodule amodule) : module(amodule) { }
    ~CudaLibrary() { cuModuleUnload (module); }

    shared_ptr<Kernel> GetKernel (const string & name) override
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


  class CudaQueue : public Queue
  {
    CUstream stream;
  public:
    CudaQueue() { Check (cuStreamCreate (&stream, CU_STREAM_NON_BLOCKING), "cuStreamCreate"); }
    ~CudaQueue() { cuStreamDestroy (stream); }

    void Finish() override
    { Check (cuStreamSynchronize (stream), "cuStreamSynchronize"); }

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

      Check (cuLaunchKernel (ck.Get(),
                             groups.x, groups.y, groups.z,
                             groupsize.x, groupsize.y, groupsize.z,
                             dynamic_group_memory, stream,
                             params.data(), nullptr), "cuLaunchKernel");
    }
  };


  class CudaDevice : public Device
  {
    CUdevice dev;
    shared_ptr<Queue> defqueue;
    int ccmajor, ccminor;

    int Attr (CUdevice_attribute a) const
    { int v = 0; cuDeviceGetAttribute (&v, a, dev); return v; }

  public:
    CudaDevice (CUdevice adev) : dev(adev)
    {
      ccmajor = Attr (CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR);
      ccminor = Attr (CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR);
      defqueue = std::make_shared<CudaQueue>();
    }

    string Name() const override
    {
      char buf[256] = "";
      cuDeviceGetName (buf, sizeof(buf), dev);
      return buf;
    }

    bool HasFloat64() const override { return true; }
    bool IsUnifiedMemory() const override
    { return Attr (CU_DEVICE_ATTRIBUTE_MANAGED_MEMORY) != 0; }
    size_t MaxThreadsPerGroup() const override
    { return Attr (CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK); }
    size_t SimdWidth() const override
    { return Attr (CU_DEVICE_ATTRIBUTE_WARP_SIZE); }

    shared_ptr<Buffer> NewBuffer (size_t bytes, MemType mt) override
    { return std::make_shared<CudaBuffer> (bytes, mt); }

    shared_ptr<Queue> DefaultQueue() override { return defqueue; }

    shared_ptr<Library> CompileSource (const string & source) override
    {
      nvrtcProgram prog;
      Check (nvrtcCreateProgram (&prog, source.c_str(), "ngsgpu.cu",
                                 0, nullptr, nullptr), "nvrtcCreateProgram");

      std::string arch = "--gpu-architecture=compute_"
        + std::to_string(ccmajor) + std::to_string(ccminor);
      const char * opts[] = { arch.c_str(), "--std=c++17", "--extended-lambda" };

      auto res = nvrtcCompileProgram (prog, 3, opts);
      if (res != NVRTC_SUCCESS)
        {
          size_t logsize = 0;
          nvrtcGetProgramLogSize (prog, &logsize);
          std::string log (logsize, '\0');
          if (logsize) nvrtcGetProgramLog (prog, log.data());
          nvrtcDestroyProgram (&prog);
          throw std::runtime_error ("ngscuda: kernel compile error:\n" + log);
        }

      size_t ptxsize = 0;
      Check (nvrtcGetPTXSize (prog, &ptxsize), "nvrtcGetPTXSize");
      std::string ptx (ptxsize, '\0');
      Check (nvrtcGetPTX (prog, ptx.data()), "nvrtcGetPTX");
      nvrtcDestroyProgram (&prog);

      CUmodule module;
      Check (cuModuleLoadData (&module, ptx.c_str()), "cuModuleLoadData");
      return std::make_shared<CudaLibrary> (module);
    }
  };


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
