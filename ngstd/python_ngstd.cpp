#ifdef NGS_PYTHON

#include "python_ngstd.hpp"
#include "gpuwrapper.hpp"
#include "gpukernel.hpp"
#include <Python.h>
#include <pybind11/numpy.h>



PythonEnvironment pyenv;


using std::ostringstream;
namespace ngstd {
  bool have_numpy = false;
}

const char* docu_string(const char* str)
{
  if(getenv("NETGEN_DOCUMENTATION_RST_FORMAT"))
    return str;
  std::string replacement(str);
  bool replaced = false;
  while(true)
    {
      auto start_pos = replacement.find(":any:`");
      if(start_pos==std::string::npos)
        break;
      else
        replaced = true;
      auto rest = replacement.substr(start_pos+6); //first character after ":any:`"
      auto inner_end = rest.find("<");
      auto end = rest.find("`");
      if(inner_end==std::string::npos)
        inner_end = end;
      replacement.replace(start_pos,end+7,rest.substr(0,inner_end)); 
    }
  if(!replaced)
    return str;
  char * newchar = new char[replacement.size()+1];
  std::copy(replacement.begin(),replacement.end(),newchar);
  newchar[replacement.size()] = '\0';
  return newchar;
}


template<typename T, typename TSCAL, typename TCLASS>
void PyDefVecBuffer( TCLASS & c )
{
  c.def_buffer([](T &self) -> py::buffer_info {
      return py::buffer_info
        (
         &self[0],                                     /* Pointer to buffer */
         sizeof(TSCAL),                                /* Size of one scalar */
         py::format_descriptor<TSCAL>::format(),       /* Python struct-style format descriptor */
         1,                                            /* Number of dimensions */
         { self.Size() },                              /* Buffer dimensions */
         { sizeof(TSCAL)*(self.Addr(1)-self.Addr(0)) } /* Strides (in bytes) for each index */
         );
    });
}




static void ExportGPU (py::module & m)
{
  using namespace ngs_gpu;

  py::enum_<MemType> (m, "MemType")
    .value("Device", MemType::Device)
    .value("Shared", MemType::Shared)
    ;

  py::class_<Buffer, shared_ptr<Buffer>> (m, "GPUBuffer")
    .def_property_readonly("size", &Buffer::Size)
    .def_property_readonly("memtype", &Buffer::GetMemType)
    .def_property_readonly("host_visible",
                           [](Buffer & self) { return self.HostPtr() != nullptr; })
    .def("H2D", [](Buffer & self, py::array_t<float, py::array::c_style> a, size_t offset)
         {
           if (offset + a.nbytes() > self.Size())
             throw Exception("GPUBuffer.H2D: does not fit");
           self.H2D (static_cast<const void*>(a.data()), (size_t)a.nbytes(), offset);
         }, py::arg("array"), py::arg("offset")=0)
    .def("D2H", [](Buffer & self, size_t n, size_t offset)
         {
           if (offset + n*sizeof(float) > self.Size())
             throw Exception("GPUBuffer.D2H: does not fit");
           py::array_t<float> a(n);
           self.D2H (static_cast<void*>(a.mutable_data()), n*sizeof(float), offset);
           return a;
         }, py::arg("n"), py::arg("offset")=0)
    .def("__str__", [](Buffer & self)
         { return "GPUBuffer, size = " + ToString(self.Size()) + " bytes"; })
    ;

  py::class_<Kernel, shared_ptr<Kernel>> (m, "GPUKernel")
    .def_property_readonly("name", &Kernel::Name)
    .def("__str__", [](Kernel & self) { return "GPUKernel " + self.Name(); })
    ;

  py::class_<Library, shared_ptr<Library>> (m, "GPULibrary")
    .def("GetKernel", &Library::GetKernel, py::arg("name"))
    ;

  py::class_<Queue, shared_ptr<Queue>> (m, "GPUQueue")
    .def("Launch", [](Queue & self, shared_ptr<Kernel> kernel,
                      std::vector<unsigned> groups, std::vector<unsigned> groupsize,
                      py::list args, size_t group_memory)
         {
           auto dim = [](std::vector<unsigned> & v, const char * what)
           {
             if (v.empty() || v.size() > 3)
               throw Exception(string("GPUQueue.Launch: ")+what+" needs 1 to 3 entries");
             v.resize(3, 1);
             return Dim3(v[0], v[1], v[2]);
           };

           std::vector<KernelArg> kargs;
           for (auto arg : args)
             {
               if (py::isinstance<Buffer>(arg))
                 kargs.push_back (KernelArg(py::cast<shared_ptr<Buffer>>(arg)));
               else if (py::isinstance<py::float_>(arg))
                 kargs.push_back (KernelArg(py::cast<float>(arg)));   // Metal is fp32
               else if (py::isinstance<py::int_>(arg))
                 kargs.push_back (KernelArg(py::cast<int32_t>(arg)));
               else
                 throw Exception("GPUQueue.Launch: argument must be GPUBuffer, float or int");
             }

           self.Launch (*kernel, dim(groups, "groups"), dim(groupsize, "groupsize"),
                        kargs, group_memory);
         },
         py::arg("kernel"), py::arg("groups"), py::arg("groupsize"),
         py::arg("args")=py::list(), py::arg("group_memory")=0)
    .def("Finish", &Queue::Finish)
    ;

  py::class_<Device, shared_ptr<Device>> (m, "GPUDevice")
    .def_property_readonly("name", &Device::Name)
    .def_property_readonly("has_float64", &Device::HasFloat64)
    .def_property_readonly("unified_memory", &Device::IsUnifiedMemory)
    .def_property_readonly("max_threads_per_group", &Device::MaxThreadsPerGroup)
    .def_property_readonly("simd_width", &Device::SimdWidth)
    .def("NewBuffer", &Device::NewBuffer,
         py::arg("bytes"), py::arg("memtype")=MemType::Shared)
    .def("CompileSource", &Device::CompileSource, py::arg("source"))
    .def("DefaultQueue", &Device::DefaultQueue)
    .def("__str__", [](Device & self) { return "GPUDevice " + self.Name(); })
    ;

  m.def("GetGPUDevice", &GetDevice,
        "the GPU device, None if no backend (ngsmetal / ngscuda) is loaded");
  m.def("HasGPUDevice", &HasDevice);
  m.def("GetCPUDevice", &GetCpuDevice,
        "host reference backend, kernels compiled with the host c++ compiler");
  m.def("SetGPUDevice", [] (shared_ptr<Device> dev)
        { SetDeviceCreator ([dev]() { return dev; }); }, py::arg("device"),
        "run gpu vectors and matrices on this device, e.g. GetCPUDevice()");

  // common kernel syntax, prepend to the source given to CompileSource
  m.attr("GPUKernelPrelude") = code_gpukernel;
}


void NGS_DLL_HEADER  ExportNgstd(py::module & m) {
  try {
      auto numpy = py::module::import("numpy");
      have_numpy = !numpy.is_none();
  }
  catch(...) {}

  
  std::string nested_name = "ngstd";

  ExportGPU (m);

  py::class_<DummyArgument>(m, "DummyArgument")
    .def("__bool__", []( DummyArgument &self ) { return false; } )
    .def("__repr__", [] ( DummyArgument & self) { return "<ngsolve.ngstd.DummyArgument>"; })
    ;
  
  py::class_<ngstd::LocalHeap> (m, "LocalHeap", "A heap for fast memory allocation")
     .def(py::init<size_t,const char*>(), "size"_a=1000000, "name"_a="PyLocalHeap")
    ;

  py::class_<ngstd::HeapReset>
    (m, "HeapReset","stores heap-pointer on init, and resets it on exit")
    .def(py::init<LocalHeap&>(), py::arg("lh"))
    ;
  
  m.def("TestFlagsConversion", []( Flags flags) { cout << flags << endl; }, py::arg("flags") );

  py::class_<ngstd::IntRange> (m, "IntRange")
    .def( py::init<size_t,size_t>())
    .def("__str__", &ToString<IntRange>)
    .def("__iter__", [] (ngstd::IntRange & i)
         { return py::make_iterator(i.begin(), i.end()); },
         py::keep_alive<0,1>())
    .def("__contains__", [] (const IntRange & self, int i)
         { return i >= self.begin() && i < self.end(); })
    .def_property_readonly("start", [](IntRange& self) { return self.First();})
    .def_property_readonly("stop", [](IntRange& self) { return self.Next();})
    .def_property_readonly("step", [](IntRange& self) { return 1; })
    ;
  
  py::class_<Archive, shared_ptr<Archive>> (m, "Archive")
      /*
    .def("__init__", [](const string & filename, bool write,
                              bool binary) -> shared_ptr<Archive>
                           {
                             if(binary) {
                               if (write)
                                return make_shared<BinaryOutArchive> (filename);
                              else
                                return make_shared<BinaryInArchive> (filename);
                              }
                              else {
                                if (write)
                                  return make_shared<TextOutArchive> (filename);
                                else
                                  return make_shared<TextInArchive> (filename);
                              }
                           })
      */
      .def(py::init<> ([](const string & filename, bool write,
                          bool binary) -> shared_ptr<Archive>
                       {
                         if(binary) {
                           if (write)
                             return make_shared<BinaryOutArchive> (filename);
                           else
                             return make_shared<BinaryInArchive> (filename);
                         }
                         else {
                           if (write)
                             return make_shared<TextOutArchive> (filename);
                           else
                             return make_shared<TextInArchive> (filename);
                         }
                       }), py::arg("filename"), py::arg("write"), py::arg("binary"))
    .def("__and__" , [](shared_ptr<Archive> & self, Array<int> & a) 
                                         { cout << "output array" << endl;
                                           *self & a; return self; }, py::arg("array"))
  ;

  m.def("_PickleMemory", [](py::object pickler, MemoryView& view)
        {
          py::buffer_info bi((char*) view.Ptr(), view.Size());
          pickler.attr("write")(py::bytes("\xf0"));
          size_t size = view.Size();
          pickler.attr("write")(py::bytes((char*) & size, sizeof(size_t)));
          pickler.attr("write")(py::memoryview(bi));
        }, py::arg("pickler"), py::arg("view"));
  m.def("_UnpickleMemory", [](py::object unpickler)
        {
          auto size = *(size_t*) PyBytes_AsString(unpickler.attr("read")(sizeof(size_t)).ptr());
          char* mem = new char[size];
          constexpr int BUFFER_SIZE = 8 * 1024 * 1024; // read 8 MB
          size_t n = 0;
          while (n + BUFFER_SIZE < size)
            {
              auto buffer = unpickler.attr("read")(BUFFER_SIZE);
              memcpy(&mem[n], PyBytes_AsString(buffer.ptr()), BUFFER_SIZE);
              n += BUFFER_SIZE;
            }
          auto buffer = unpickler.attr("read")(size-n);
          memcpy(&mem[n], PyBytes_AsString(buffer.ptr()), size-n);
          unpickler.attr("append")(MemoryView(mem,size));
        }, py::arg("unpickler"));

  py::class_<MemoryView>(m, "_MemoryView")
    .def(py::pickle([](MemoryView& mv) -> py::tuple
                    {
                      if(have_numpy)
                          return py::make_tuple(true,
                                                py::array_t<char> (py::buffer_info((char*) mv.Ptr(),
                                                                                   mv.Size())));
                      else
                          return py::make_tuple(false,
                                                py::bytes((char*) mv.Ptr(),
                                                          mv.Size()));
                    },
                    [](py::tuple state)
                    {
                      if(py::cast<bool>(state[0]))
                        {
                          if(!have_numpy)
                            throw Exception("Data was pickled using numpy, need numpy to unpickle it!");
                          auto array = py::cast<py::array_t<char>>(state[1]);
                          auto size = array.size();
                          char* mem = new char[size];
                          memcpy(mem, array.data(0), size);
                          return MemoryView (mem, size);
                        }
                      else
                        {
                          auto bytes = py::cast<py::bytes>(state[1]);
                          char* buffer;
                          py::ssize_t size;
                          PYBIND11_BYTES_AS_STRING_AND_SIZE(bytes.ptr(), &buffer, &size);
                          char *mem = new char[size];
                          memcpy(mem, buffer, size);
                          return MemoryView(mem, (size_t) size);
                        }
                    }))
    ;

}

#endif // NGS_PYTHON
