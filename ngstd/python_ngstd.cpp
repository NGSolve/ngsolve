#ifdef NGS_PYTHON

#include "python_ngstd.hpp"
#include <Python.h>
#include <pybind11/numpy.h>



PythonEnvironment pyenv;


using std::ostringstream;
namespace ngstd {
  bool have_numpy = false;
}

void SetFlag(Flags &flags, string s, py::object value) 
{
  if (py::isinstance<py::dict>(value))
    {             
      py::dict vdd(value);
      // call recursively to set dictionary
      for (auto item : vdd) {
        string name = item.first.cast<string>();
        py::object val = py::reinterpret_borrow<py::object>(item.second);
        SetFlag(flags, name, val);
      }
      return;
    }

  if (py::isinstance<py::bool_>(value))
    flags.SetFlag(s, value.cast<bool>());

  if (py::isinstance<py::float_>(value))
    flags.SetFlag(s, value.cast<double>());

  if (py::isinstance<py::int_>(value))
    flags.SetFlag(s, double(value.cast<int>()));

  if (py::isinstance<py::str>(value))
    flags.SetFlag(s, value.cast<string>());

  if (py::isinstance<py::list>(value))
    {             
      py::list vdl(value);
      if (py::len(vdl) > 0)
      {
        if(py::CheckCast<double>(vdl[0]))
          flags.SetFlag(s, makeCArray<double>(vdl));
        if(py::isinstance<py::str>(vdl[0]))
          flags.SetFlag(s, makeCArray<string>(vdl));
      }
      else
      {
        Array<string> dummystr;
        Array<double> dummydbl;
        flags.SetFlag(s,dummystr);
        flags.SetFlag(s,dummydbl);
      }
    }

  if (py::isinstance<py::tuple>(value))
    {
      py::tuple vdt(value);
      if (py::isinstance<py::float_>(value))
        flags.SetFlag(s, makeCArray<double>(vdt));
      if (py::isinstance<py::int_>(value))
        flags.SetFlag(s, makeCArray<double>(vdt));
      if (py::isinstance<py::str>(value))
        flags.SetFlag(s, makeCArray<string>(vdt));
    }
}

Flags NGS_DLL_HEADER CreateFlagsFromKwArgs(py::object pyclass, const py::kwargs& kwargs, py::list info)
{
  auto flags_doc = pyclass.attr("__flags_doc__")();
  py::dict flags_dict;

  if (kwargs.contains("flags"))
    {
      cout  << IM(2) << "WARNING: using flags as kwarg is deprecated in " << py::str(pyclass)
           << ", use the flag arguments as kwargs instead!" << endl;
      auto addflags = py::cast<py::dict>(kwargs["flags"]);
      for (auto item : addflags)
        flags_dict[item.first.cast<string>().c_str()] = item.second;
    }
  for (auto item : kwargs)
    if (!flags_doc.contains(item.first.cast<string>().c_str()) &&
        !(item.first.cast<string>() == "flags"))
        cout << IM(2) << "WARNING: kwarg '" << item.first.cast<string>()
             << "' is an undocumented flags option for class "
             << std::string(py::str(pyclass)) << ", maybe there is a typo?" << endl;

  py::dict special;
  if(py::hasattr(pyclass,"__special_treated_flags__"))
      special = pyclass.attr("__special_treated_flags__")();
  for (auto item : kwargs)
    {
      auto name = item.first.cast<string>();
      if (name != "flags")
        {
          if(!special.contains(name.c_str()))
            flags_dict[name.c_str()] = item.second;
        }
    }

  auto flags = py::extract<Flags>(flags_dict)();

  for (auto item : kwargs)
    {
      auto name = item.first.cast<string>();
      if (name != "flags")
        {
          if(special.contains(name.c_str()))
            special[name.c_str()](item.second, &flags, info);
        }
    }
  return flags;
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

  m.def("TestFlags", [] (py::dict const &d) { cout << py::extract<Flags>(d)() << endl; } );

  py::class_<DummyArgument>(m, "DummyArgument")
    .def("__bool__", []( DummyArgument &self ) { return false; } )
    .def("__repr__", [] ( DummyArgument & self) { return "<ngsolve.ngstd.DummyArgument>"; })
    ;
  
  py::class_<PajeTrace >(m, "Tracer")
    .def("SetTraceThreads", &PajeTrace::SetTraceThreads)
    .def("SetTraceThreadCounter", &PajeTrace::SetTraceThreadCounter)
    .def("SetMaxTracefileSize", &PajeTrace::SetMaxTracefileSize)
    ;




  
  py::class_<FlatArray<double> > class_flatarrayd (m, "FlatArrayD", py::buffer_protocol());
  class_flatarrayd.def(py::init<size_t, double *>());
  PyDefVector<FlatArray<double>, double>(m, class_flatarrayd);
  PyDefVecBuffer<FlatArray<double>,double>(class_flatarrayd);  
  PyDefToString<FlatArray<double>>(m, class_flatarrayd);
  
  py::class_<Array<double>, FlatArray<double> >(m, "ArrayD")
    .def(py::init([] (int n) { return new Array<double>(n); }), py::arg("n"), "Makes array of given length")
    .def(py::init([] (std::vector<double> const & x)
                  {
                    int s = x.size();
                    Array<double> a(s);
                    for (int i = 0; i < s; i++)
                      a[i] = x[i];
                    return a;
                  }), py::arg("vec"), "Makes array with given vector")
    .def("__rand__" ,  []( Array<double> & a, shared_ptr<Archive> & arch )
                                         { cout << "output d array" << endl;
                                           *arch & a; return arch; })
    .def("print", [](Array<double> &a) { cout << a << endl; }, "print the array" )
    ;

  py::class_<FlatArray<int> > class_flatarrayi (m, "FlatArrayI", py::buffer_protocol());
  PyDefVector<FlatArray<int>, int>(m, class_flatarrayi);
  PyDefToString<FlatArray<int> >(m, class_flatarrayi);
  PyDefVecBuffer<FlatArray<int>,int>(class_flatarrayi);
  class_flatarrayi.def(py::init<size_t, int *>());

  py::class_<Array<int>, FlatArray<int> >(m, "ArrayI")
    .def(py::init([] (int n) { return new Array<int>(n); }),py::arg("n"), "Makes array of given length")
    .def(py::init([] (std::vector<int> const & x)
                  {
                    int s = x.size();
                    Array<int> tmp(s);
                    for (int i = 0; i < s; i++)
                      tmp[i] = x[i]; 
                    return tmp;
                  }), py::arg("vec"), "Makes array with given vector")
    ;

  py::class_<FlatArray<size_t> > class_flatarrayst (m, "FlatArray_sizet", py::buffer_protocol());
  PyDefVector<FlatArray<size_t>, size_t>(m, class_flatarrayst);
  PyDefToString<FlatArray<size_t> >(m, class_flatarrayst);
  PyDefVecBuffer<FlatArray<size_t>, size_t>(class_flatarrayst);
  class_flatarrayst.def(py::init<size_t, size_t *>());
  
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

  

  py::class_<Flags>(m, "Flags")
    .def(py::init<>())
    .def("__str__", &ToString<Flags>)
    .def(py::init([](py::object & obj) {
          Flags flags;
          py::dict d(obj);          
          SetFlag (flags, "", d);
          return flags;
        }), py::arg("obj"), "Create Flags by given object")
    /*
    .def("__init__", [] (Flags &f, py::object & obj) {
         py::dict d(obj);
         new (&f) Flags();
         SetFlag(f, "", d);
     })
    */
    .def("__getstate__", [] (py::object self_object) {
        auto self = self_object.cast<Flags>();
        stringstream str;
        self.SaveFlags(str);
        return py::make_tuple(py::cast(str.str())); 
      }, "Return state of the flags")
    .def("__setstate__", [] (Flags & self, py::tuple state) {
        string s = state[0].cast<string>();
        stringstream str(s);
        new (&self) Flags();
        self.LoadFlags(str);
      }, py::arg("state"), "Set the state of the flags" )
    .def("Set",[](Flags & self,const py::dict & aflags)->Flags&
    {      
      cout << "we call Set(dict)" << endl;
      SetFlag(self, "", aflags);
      return self;
    }, py::arg("aflag"), "Set the flags by given dict")

    .def("Set",[](Flags & self, const char * akey, const py::object & value)->Flags&
    {             
      cout << "we call Set(key,obj)" << endl; 
        SetFlag(self, akey, value);
        return self;
    }, py::arg("akey"), py::arg("value"), "Set flag by given value.")

    .def("__getitem__", [](Flags & self, const string& name) -> py::object {

	  if(self.NumListFlagDefined(name))
	    return py::cast(self.GetNumListFlag(name));

	  if(self.StringListFlagDefined(name))
	    return py::cast(self.GetStringListFlag(name));
	 
	  if(self.NumFlagDefined(name))
	    return py::cast(*self.GetNumFlagPtr(name));
	  
	  if(self.StringFlagDefined(name))
	    return py::cast(self.GetStringFlag(name));

	  if(self.FlagsFlagDefined(name))
	    return py::cast(self.GetFlagsFlag(name));

	  return py::cast(self.GetDefineFlag(name));
      }, py::arg("name"), "Return flag by given name")
  ;

  m.def("TestFlagsConversion", []( Flags flags) { cout << flags << endl; }, py::arg("flags") );
  py::implicitly_convertible<py::dict, Flags>();

  py::class_<ngstd::IntRange> py_intrange (m, "IntRange");
  py_intrange.def( py::init<size_t,size_t>());
  py_intrange.def("__str__", &ToString<IntRange>);
  py_intrange.def("__iter__", [] (ngstd::IntRange & i)
      { return py::make_iterator(i.begin(), i.end()); },
      py::keep_alive<0,1>()
    );

  py::class_<Timer> (m, "Timer")
    .def(py::init<const string&>())
    .def("Start", &Timer::Start, "start timer")
    .def("Stop", &Timer::Stop, "stop timer")
    ;
  
  m.def("Timers",
	  []() 
	   {
	     py::list timers;
	     for (int i = 0; i < NgProfiler::SIZE; i++)
	       if (!NgProfiler::names[i].empty())
               {
                 py::dict timer;
                 timer["name"] = py::str(NgProfiler::names[i]);
                 timer["time"] = py::float_(NgProfiler::GetTime(i));
                 timer["counts"] = py::int_(NgProfiler::GetCounts(i));
                 timer["flops"] = py::int_(NgProfiler::GetFlops(i));
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
  py::class_<MemoryView>(m, "_MemoryView");

  
  py::class_<PyMPI_Comm> (m, "MPI_Comm")
    .def_property_readonly ("rank", &PyMPI_Comm::Rank)
    .def_property_readonly ("size", &PyMPI_Comm::Size)
    .def("Barrier", [](PyMPI_Comm c) { MyMPI_Barrier(c.comm); })
#ifdef PARALLEL
    .def("WTime", [](PyMPI_Comm c) { return MPI_Wtime(); })
#endif
    .def("Sum", [](PyMPI_Comm c, double x) { return MyMPI_AllReduce(x, MPI_SUM, c.comm); })
    .def("Min", [](PyMPI_Comm c, double x) { return MyMPI_AllReduce(x, MPI_MIN, c.comm); })
    .def("Max", [](PyMPI_Comm c, double x) { return MyMPI_AllReduce(x, MPI_MAX, c.comm); })
    .def("Sum", [](PyMPI_Comm c, int x) { return MyMPI_AllReduce(x, MPI_SUM, c.comm); })
    .def("Min", [](PyMPI_Comm c, int x) { return MyMPI_AllReduce(x, MPI_MIN, c.comm); })
    .def("Max", [](PyMPI_Comm c, int x) { return MyMPI_AllReduce(x, MPI_MAX, c.comm); })
    .def("Sum", [](PyMPI_Comm c, size_t x) { return MyMPI_AllReduce(x, MPI_SUM, c.comm); })
    .def("Min", [](PyMPI_Comm c, size_t x) { return MyMPI_AllReduce(x, MPI_MIN, c.comm); })
    .def("Max", [](PyMPI_Comm c, size_t x) { return MyMPI_AllReduce(x, MPI_MAX, c.comm); })
    ;

  
  
  m.def("MPI_Init", [&]()
        {
          const char * progname = "ngslib";
          typedef const char * pchar;
          pchar ptrs[2] = { progname, nullptr };
          pchar * pptr = &ptrs[0];
          
          static MyMPI mympi(1, (char**)pptr);
          return PyMPI_Comm(ngs_comm);
        });

}


PYBIND11_MODULE(libngstd, m) {
  m.attr("__name__") = "ngstd";
  ExportNgstd(m);
}


#endif // NGS_PYTHON
