/*********************************************************************/
/* File:   gpu_cpu.cpp                                               */
/* Date:   29. Aug. 2025                                             */
/*********************************************************************/

/*
  Reference backend on the host. Kernels written in the common syntax
  (gpukernel.hpp) are compiled with the host c++ compiler and run with
  one thread per work-item, so results can be compared against cuda or
  metal without a gpu.
*/

#include "gpuwrapper.hpp"
#include "gpukernel.hpp"

#include <dlfcn.h>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <stdexcept>
#include <thread>

namespace ngs_gpu
{
  namespace
  {
    typedef void (*thunk_function) (void**);
    typedef void (*launch_function) (thunk_function, void**,
                                     unsigned, unsigned, unsigned,
                                     unsigned, unsigned, unsigned, size_t);

    void Err (const std::string & msg)
    { throw std::runtime_error ("ngs_gpu (cpu): " + msg); }


    class CpuBuffer : public Buffer
    {
      std::vector<char> mem;
    public:
      CpuBuffer (size_t bytes, MemType mt) : Buffer(bytes, mt), mem(bytes) { }
      void * Data() const { return (void*) mem.data(); }
    protected:
      void * DoHostPtr() const override { return (void*) mem.data(); }
      void DoH2D (const void * src, size_t bytes, size_t offset) override
      {
        if (offset+bytes > size) Err ("H2D out of range");
        std::memcpy ((char*)mem.data()+offset, src, bytes);
      }
      void DoD2H (void * dst, size_t bytes, size_t offset) const override
      {
        if (offset+bytes > size) Err ("D2H out of range");
        std::memcpy (dst, mem.data()+offset, bytes);
      }
    };


    class CpuKernel : public Kernel
    {
      std::string name;
      thunk_function thunk;
    public:
      CpuKernel (const std::string & aname, thunk_function athunk)
        : name(aname), thunk(athunk) { }
      string Name() const override { return name; }
      thunk_function Get() const { return thunk; }
    };


    class CpuLibrary : public Library
    {
      void * handle;
      std::filesystem::path dir;
    public:
      CpuLibrary (void * ahandle, std::filesystem::path adir)
        : handle(ahandle), dir(adir) { }
      ~CpuLibrary()
      {
        if (handle) dlclose (handle);
        if (std::getenv ("NGS_GPU_KEEP")) return;
        std::error_code ec;
        std::filesystem::remove_all (dir, ec);
      }

      launch_function Launcher()
      {
        auto f = (launch_function) dlsym (handle, "_ngs_launch");
        if (!f) Err ("_ngs_launch not found - was the prelude prepended?");
        return f;
      }

      shared_ptr<Kernel> DoGetKernel (const string & name) override
      {
        auto f = (thunk_function) dlsym (handle, (name+"_ngsthunk").c_str());
        if (!f) Err ("no kernel '" + name + "' in library"
                     " (entry points must be declared with KERNEL(...))");
        return std::make_shared<CpuKernel> (name, f);
      }
    };


    class CpuQueue : public Queue
    {
      launch_function launch;
    public:
      CpuQueue (launch_function alaunch) : launch(alaunch) { }
      void Finish() override { }     // launches run synchronously

    protected:
      void DoLaunch (Kernel & kernel, Dim3 groups, Dim3 groupsize,
                     const std::vector<KernelArg> & args,
                     size_t dynamic_group_memory) override
      {
        auto & ck = dynamic_cast<CpuKernel&> (kernel);

        // buffers are passed as the pointer itself, values as a pointer to them
        std::vector<void*> argv (args.size());
        for (size_t i = 0; i < args.size(); i++)
          {
            auto & a = args[i];
            if (a.GetKind() == KernelArg::Kind::Buffer)
              {
                auto cb = dynamic_cast<CpuBuffer*> (a.GetBuffer());
                if (!cb) Err ("kernel argument is not a cpu buffer");
                argv[i] = (char*)cb->Data() + a.Offset();
              }
            else
              argv[i] = const_cast<void*> (a.Data());
          }

        launch (ck.Get(), argv.data(),
                groups.x, groups.y, groups.z,
                groupsize.x, groupsize.y, groupsize.z,
                dynamic_group_memory);
      }
    };


    class CpuDevice : public Device
    {
      shared_ptr<Queue> defqueue;     // created with the first library
      std::vector<shared_ptr<Library>> libs;
    public:
      string Name() const override { return "CPU reference"; }
      bool HasFloat64() const override { return true; }
      bool IsUnifiedMemory() const override { return true; }
      size_t MaxThreadsPerGroup() const override { return 1024; }
      size_t SimdWidth() const override { return 1; }

      shared_ptr<Buffer> NewBuffer (size_t bytes, MemType mt) override
      { return std::make_shared<CpuBuffer> (bytes, mt); }

      shared_ptr<Queue> DefaultQueue() override
      {
        if (!defqueue) Err ("compile a source first, the cpu queue needs its runtime");
        return defqueue;
      }

      shared_ptr<Library> CompileSource (const string & source) override
      {
        namespace fs = std::filesystem;
        auto dir = fs::temp_directory_path() /
          ("ngsgpu_" + std::to_string(::getpid()) + "_" + std::to_string(libs.size()));
        fs::create_directories (dir);

        auto src = dir / "kernel.cpp";
        auto lib = dir / "kernel.so";
        std::ofstream(src) << source;

        const char * cxx = std::getenv ("NGS_GPU_CXX");
        const char * flags = std::getenv ("NGS_GPU_KEEP") ? " -O0 -g " : " -O2 ";
        std::string cmd = std::string(cxx ? cxx : "c++")
          + " -std=c++20 -fPIC -shared -DNGS_GPU_CPU -pthread " + flags
          + src.string() + " -o " + lib.string() + " 2>" + (dir/"log").string();

        if (std::system (cmd.c_str()) != 0)
          {
            std::ifstream log (dir/"log");
            std::string txt ((std::istreambuf_iterator<char>(log)),
                             std::istreambuf_iterator<char>());
            Err ("kernel compile error:\n" + txt);
          }

        if (std::getenv ("NGS_GPU_KEEP"))
          std::cout << "ngs_gpu (cpu): generated kernel in " << src.string() << std::endl;

        void * handle = dlopen (lib.string().c_str(), RTLD_NOW | RTLD_LOCAL);
        if (!handle) Err (std::string("dlopen failed: ") + dlerror());

        auto res = std::make_shared<CpuLibrary> (handle, dir);
        if (!defqueue) defqueue = std::make_shared<CpuQueue> (res->Launcher());
        libs.push_back (res);
        return res;
      }
    };
  }

  shared_ptr<Device> GetCpuDevice()
  {
    static shared_ptr<Device> dev = std::make_shared<CpuDevice>();
    return dev;
  }
}
