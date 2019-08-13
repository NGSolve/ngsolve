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
  
  py::class_<PajeTrace >(m, "Tracer")
    .def("SetTraceThreads", &PajeTrace::SetTraceThreads)
    .def("SetTraceThreadCounter", &PajeTrace::SetTraceThreadCounter)
    .def("SetMaxTracefileSize", &PajeTrace::SetMaxTracefileSize)
    ;


  py::class_<ngstd::LocalHeap> (m, "LocalHeap", "A heap for fast memory allocation")
     .def(py::init<size_t,const char*>(), "size"_a=1000000, "name"_a="PyLocalHeap")
    ;

  py::class_<ngstd::HeapReset>
    (m, "HeapReset","stores heap-pointer on init, and resets it on exit")
    .def(py::init<LocalHeap&>(), py::arg("lh"))
    ;
  
  py::class_<ngstd::BitArray, shared_ptr<BitArray>> (m, "BitArray")
    // .def(py::init<size_t>()) // not doing the right thing ????? JS
    .def(py::init([] (size_t n) { return make_shared<BitArray>(n); }),py::arg("n"))
    .def(py::init([] (const BitArray& a) { return make_shared<BitArray>(a); } ), py::arg("ba"))
    .def(py::init([] (const vector<bool> & a)
                  {
                    auto ba = make_shared<BitArray>(a.size());
                    ba->Clear();
                    for (size_t i = 0; i < a.size(); i++)
                      if (a[i]) ba->Set(i);
                    return ba;
                  } ), py::arg("vec"))
    .def("__str__", &ToString<BitArray>)
    .def("__len__", &BitArray::Size)
    .def("__getitem__", [] (BitArray & self, int i) 
                                         {
                                           if (i < 0 || i >= self.Size())
                                             throw py::index_error();
                                           return self.Test(i); 
                                         }, py::arg("pos"), "Returns bit from given position")
    .def("__setitem__", [] (BitArray & self, int i, bool b) 
                                         {
                                           if (i < 0 || i >= self.Size())
                                             throw py::index_error();
                                           if (b) self.Set(i); else self.Clear(i); 
                                         }, py::arg("pos"), py::arg("value"), "Clear/Set bit at given position")

    .def("__setitem__", [] (BitArray & self, py::slice inds, bool b) 
                                         {
                                           size_t start, step, n;
                                           InitSlice( inds, self.Size(), start, step, n );
                                           if (start == 0 && n == self.Size() && step == 1)
                                             { // base branch
                                               if (b)
                                                 self.Set();
                                               else
                                                 self.Clear();
                                             }
                                           else
                                             {
                                               if (b)
                                                 for (size_t i=0; i<n; i++, start+=step)
                                                   self.Set(start);
                                               else
                                                 for (size_t i=0; i<n; i++, start+=step)
                                                   self.Clear(start);
                                             }
                                         }, py::arg("inds"), py::arg("value"), "Clear/Set bit at given positions")
    
    .def("__setitem__", [](BitArray & self,  IntRange range, bool b)
      {
        if (b)
          for (size_t i : range)
            self.Set(i);
        else
          for (size_t i : range)
            self.Clear(i);
      }, py::arg("range"), py::arg("value"), "Set value for range of indices" )
    
    .def("NumSet", [] (BitArray & self) { return self.NumSet(); })
    .def("Set", [] (BitArray & self, py::object in)
                                   {
                                     if (py::isinstance<DummyArgument>(in))
                                       self.Set();
                                     else if (py::extract<int>(in).check())
                                     {
                                       int i = py::extract<int>(in)();
                                       if (i < 0 || i >= self.Size()) 
                                         throw py::index_error();
                                       self.Set(i);
                                     }
                                     else
                                       throw py::value_error();
                                   }, py::arg("i") = DummyArgument(), "Set bit at given position")
    .def("Clear", [] (BitArray & self, py::object in)
                                   {
                                     if (py::isinstance<DummyArgument>(in))
                                       self.Clear();
                                     else if (py::extract<int>(in).check())
                                     {
                                       int i = py::extract<int>(in)();
                                       if (i < 0 || i >= self.Size()) 
                                         throw py::index_error();
                                       self.Clear(i);
                                     }
                                     else
                                       throw py::value_error();
                                   }, py::arg("i") = DummyArgument(), "Clear bit at given position")


    .def(py::self | py::self)
    .def(py::self & py::self)
    .def(py::self |= py::self)
    .def(py::self &= py::self)
    .def(~py::self)
    ;

  m.def("TestFlagsConversion", []( Flags flags) { cout << flags << endl; }, py::arg("flags") );

  py::class_<ngstd::IntRange> (m, "IntRange")
    .def( py::init<size_t,size_t>())
    .def("__str__", &ToString<IntRange>)
    .def("__iter__", [] (ngstd::IntRange & i)
         { return py::make_iterator(i.begin(), i.end()); },
         py::keep_alive<0,1>())
    .def_property_readonly("start", [](IntRange& self) { return self.First();})
    .def_property_readonly("stop", [](IntRange& self) { return self.Next();})
    .def_property_readonly("step", [](IntRange& self) { return 1; })
    ;

  py::class_<Timer> (m, "Timer")
    .def(py::init<const string&>())
    .def("Start", &Timer::Start, "start timer")
    .def("Stop", &Timer::Stop, "stop timer")
    .def_property_readonly("time", &Timer::GetTime, "returns time")
    .def("__enter__", &Timer::Start)
    .def("__exit__", [](Timer& t, py::object, py::object, py::object)
                     {
                       t.Stop();
                     })
    ;
  
  m.def("Timers",
	  []() 
	   {
	     py::list timers;
	     for (int i = 0; i < NgProfiler::SIZE; i++)
	       if (!NgProfiler::timers[i].name.empty())
               {
                 py::dict timer;
                 timer["name"] = py::str(NgProfiler::timers[i].name);
                 timer["time"] = py::float_(NgProfiler::GetTime(i));
                 timer["counts"] = py::int_(NgProfiler::GetCounts(i));
                 timer["flops"] = py::float_(NgProfiler::GetFlops(i));
                 timer["Gflop/s"] = py::float_(NgProfiler::GetFlops(i)/NgProfiler::GetTime(i)*1e-9);
                 timers.append(timer);
               }
	     return timers;
	   }, "Returns list of timers"
	   );

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

  m.def("RunWithTaskManager", 
          [](py::object lam)
                           {
                             cout << IM(3) << "running Python function with task-manager:" << endl;
                             RunWithTaskManager ([&] () { lam(); });
                           }, py::arg("lam"), docu_string(R"raw_string(
Parameters:

lam : object
  input function

)raw_string"))
          ;

  m.def("SetNumThreads", &TaskManager::SetNumThreads, py::arg("threads"), docu_string(R"raw_string(
Set number of threads

Parameters:

threads : int
  input number of threads

)raw_string") );

  // local TaskManager class to be used as context manager in Python
  class ParallelContextManager {
      int num_threads;
    public:
      ParallelContextManager() : num_threads(0) {};
      ParallelContextManager(size_t pajesize) : num_threads(0) {
        TaskManager::SetPajeTrace(pajesize > 0);
        PajeTrace::SetMaxTracefileSize(pajesize);
      }
      void Enter() {num_threads = EnterTaskManager(); }
      void Exit(py::object exc_type, py::object exc_value, py::object traceback) {
          ExitTaskManager(num_threads);
      }
    };

  py::class_<ParallelContextManager>(m, "TaskManager")
    .def(py::init<>())
    .def(py::init<size_t>(), "pajetrace"_a, "Run paje-tracer, specify buffersize in bytes")
    .def("__enter__", &ParallelContextManager::Enter)
    .def("__exit__", &ParallelContextManager::Exit)
    .def("__timing__", &TaskManager::Timing)
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


PYBIND11_MODULE(libngstd, m) {
  m.attr("__name__") = "ngstd";
  ExportNgstd(m);
}


#endif // NGS_PYTHON
