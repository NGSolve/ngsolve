#ifdef NGS_PYTHON

#include "python_ngstd.hpp"
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



void NGS_DLL_HEADER  ExportNgstd(py::module & m) {
  try {
      auto numpy = py::module::import("numpy");
      have_numpy = !numpy.is_none();
  }
  catch(...) {}

  
  std::string nested_name = "ngstd";

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
    .def(py::pickle([](MemoryView& mv)
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
