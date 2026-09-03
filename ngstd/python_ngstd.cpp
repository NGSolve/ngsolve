#ifdef NGS_PYTHON

#include "python_ngstd.hpp"
#include "gpuwrapper.hpp"
#include "gpukernel.hpp"
#include "tinybla.hpp"
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

  // a device buffer with a numpy dtype: sizes and offsets count elements
  struct PyGPUBuffer
  {
    shared_ptr<Buffer> buf;
    py::dtype dtype;
    size_t ItemSize() const { return dtype.itemsize(); }
    size_t Size() const { return buf->Size() / ItemSize(); }
  };

  // numpy dtype -> kernel argument type
  static auto argtype = [] (const py::dtype & dt)
  {
    ArgType::Kind k = ArgType::Unknown;
    switch (dt.kind())
      {
      case 'f': k = ArgType::Float; break;
      case 'i': k = ArgType::Int; break;
      case 'u': k = ArgType::UInt; break;
      case 'b': k = ArgType::Bool; break;
      default: break;
      }
    return ArgType (k, unsigned(dt.itemsize()));
  };

  py::class_<PyGPUBuffer, shared_ptr<PyGPUBuffer>> (m, "GPUBuffer", py::buffer_protocol())
    .def_property_readonly("size", &PyGPUBuffer::Size, "number of elements")
    .def_property_readonly("nbytes", [](PyGPUBuffer & self) { return self.buf->Size(); })
    .def_property_readonly("dtype", [](PyGPUBuffer & self) { return self.dtype; })
    .def_property_readonly("memtype", [](PyGPUBuffer & self) { return self.buf->GetMemType(); })
    .def_property_readonly("host_visible",
                           [](PyGPUBuffer & self) { return self.buf->HostPtr() != nullptr; })
    .def("__len__", &PyGPUBuffer::Size)
    .def("H2D", [](PyGPUBuffer & self, py::object obj, size_t offset)
         {
           auto a = py::module::import("numpy").attr("ascontiguousarray")
             (obj, py::arg("dtype")=self.dtype).cast<py::array>();
           size_t n = a.size();
           if (offset + n > self.Size())
             throw Exception("GPUBuffer.H2D: does not fit");
           self.buf->H2D (a.data(), n*self.ItemSize(), offset*self.ItemSize());
         }, py::arg("array"), py::arg("offset")=0,
         "upload an array (converted to the buffer's dtype), offset in elements")
    .def("D2H", [](PyGPUBuffer & self, py::object nobj, size_t offset)
         {
           if (offset > self.Size())
             throw Exception("GPUBuffer.D2H: offset out of range");
           size_t n = nobj.is_none() ? self.Size()-offset : py::cast<size_t>(nobj);
           if (offset + n > self.Size())
             throw Exception("GPUBuffer.D2H: does not fit");
           py::array a(self.dtype, std::vector<py::ssize_t>{ py::ssize_t(n) });
           self.buf->D2H (a.mutable_data(), n*self.ItemSize(), offset*self.ItemSize());
           return a;
         }, py::arg("n")=py::none(), py::arg("offset")=0,
         "download n elements (default: all) as a numpy array, offset in elements")
    // numpy view without a copy, host visible buffers only
    .def_buffer([](PyGPUBuffer & self) -> py::buffer_info
         {
           void * ptr = self.buf->HostPtr();
           if (!ptr)
             throw Exception("GPUBuffer is not host visible, use D2H");
           return py::buffer_info (ptr, py::ssize_t(self.ItemSize()),
                                   py::cast<std::string>(self.dtype.attr("char")),
                                   1, { py::ssize_t(self.Size()) }, { py::ssize_t(self.ItemSize()) });
         })
    .def("numpy", [](py::object self)
         {
           if (!py::cast<PyGPUBuffer&>(self).buf->HostPtr())
             throw Exception("GPUBuffer.numpy: not host visible, use D2H");
           return py::module::import("numpy").attr("asarray")(self);
         }, "numpy view of a host visible buffer, no copy")
    .def("__str__", [](PyGPUBuffer & self)
         { return "GPUBuffer, size = " + ToString(self.Size()) + ", dtype = "
             + py::cast<std::string>(py::str(self.dtype)); })
    ;

  py::class_<Kernel, shared_ptr<Kernel>> (m, "GPUKernel")
    .def_property_readonly("name", &Kernel::Name)
    .def("Info", [](Kernel & self, size_t groupsize)
         {
           auto ki = self.Info(groupsize);
           py::dict d;
           d["registers"] = ki.registers; d["local_bytes"] = ki.local_bytes;
           d["shared_bytes"] = ki.shared_bytes;
           d["max_threads_per_group"] = ki.max_threads_per_group;
           d["max_groups_per_unit"] = ki.max_groups_per_unit;
           return d;
         }, py::arg("groupsize") = 0, "resource usage: registers, local, shared, occupancy")
    .def_property_readonly("signature", [](Kernel & self) -> py::object
         {
           auto sig = self.Signature();
           if (!sig) return py::none();
           py::list l;
           for (auto & a : sig->args)
             l.append (py::make_tuple (a.kind == KernelSignature::Buffer ? "buffer" :
                                       a.kind == KernelSignature::Value ? "value" : "?",
                                       a.type.ToString(), a.name));
           return l;
         }, "declared arguments as (kind, type, name), None if not declared with KERNEL(...)")
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

           auto npgeneric = py::module::import("numpy").attr("generic");
           std::vector<KernelArg> kargs;
           for (auto arg : args)
             {
               if (py::isinstance<PyGPUBuffer>(arg))
                 {
                   auto & b = py::cast<PyGPUBuffer&>(arg);
                   kargs.push_back (KernelArg(b.buf).WithType(argtype(b.dtype)));
                 }
               else if (py::isinstance(arg, npgeneric))
                 {
                   // numpy scalar: passed with its own dtype (np.float64 -> double)
                   auto bytes = py::cast<std::string>(arg.attr("tobytes")());
                   auto dt = py::dtype::from_args(arg.attr("dtype"));
                   kargs.push_back (KernelArg(bytes.data(), bytes.size(), argtype(dt)));
                 }
               else if (py::isinstance<py::float_>(arg))
                 kargs.push_back (KernelArg(py::cast<float>(arg)));   // Metal is fp32
               else if (py::isinstance<py::int_>(arg))
                 kargs.push_back (KernelArg(py::cast<int32_t>(arg)));
               else
                 throw Exception("GPUQueue.Launch: argument must be GPUBuffer, numpy scalar, float or int");
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
    .def("NewBuffer", [](Device & self, size_t size, py::object dtype, MemType mt)
         {
           auto dt = py::dtype::from_args(dtype);
           return make_shared<PyGPUBuffer> (PyGPUBuffer{ self.NewBuffer(size*dt.itemsize(), mt), dt });
         }, py::arg("size"), py::arg("dtype"), py::arg("memtype")=MemType::Shared,
         "buffer of size elements of dtype")
    .def("NewBuffer", [](Device & self, py::object arraylike, MemType mt)
         {
           // dtype, size and contents from the array, multi-dim arrays are flattened
           auto a = py::module::import("numpy").attr("ascontiguousarray")(arraylike).cast<py::array>();
           auto buf = self.NewBuffer (size_t(a.nbytes()), mt);
           buf->H2D (a.data(), size_t(a.nbytes()));
           return make_shared<PyGPUBuffer> (PyGPUBuffer{ buf, a.dtype() });
         }, py::arg("array"), py::arg("memtype")=MemType::Shared,
         "buffer with the dtype and contents of an array-like")
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

  // small linear algebra inside a kernel, prepend after GPUKernelPrelude
  m.attr("TinyBlaPrelude") = code_tinybla;
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
