#ifdef NGS_PYTHON

// #include <string>
// #include <ostream>
// #include <type_traits>

#include "python_ngstd.hpp"

#ifdef PARALLEL
bool MPIManager::initialized_by_me = false;
static MPIManager mpiman;
#endif

PythonEnvironment pyenv;


using std::ostringstream;

void SetFlag(Flags &flags, string s, py::object value) 
{
  py::dict vdd(value);
  if (vdd.check())
    {             
      // call recursively to set dictionary
      for (auto item : vdd) {
        string name = item.first.cast<string>();
        py::object val(item.second, true);
        SetFlag(flags, name, val);
      }
      return;
    }
  py::bool_ vb(value);
  if (vb.check() && bool(vb))
    flags.SetFlag(s);

  py::float_ vd(value);
  if (vd.check())
    flags.SetFlag(s, double(vd));         

  py::int_ vi(value);
  if (vi.check())
    flags.SetFlag(s, int(vi));         


  py::str vs(value);
  if (vs.check())
    flags.SetFlag(s, string(vs));

  py::list vdl(value);
  if (vdl.check())
    {             
      if (py::len(vdl) > 0)
      {
        py::float_ d0(vdl[0]);
        py::int_ i0(vdl[0]);
        py::str s0(vdl[0]);
        if(d0.check())
          flags.SetFlag(s, makeCArray<double>(vdl));
        if(i0.check())
          flags.SetFlag(s, makeCArray<double>(vdl));
        if (s0.check())
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

  py::tuple vdt(value);
  if (vdt.check())
    {
      py::float_ d0(vdt[0]);
      py::int_ i0(vdl[0]);
      py::str s0(vdt[0]);
      if (d0.check())
        flags.SetFlag(s, makeCArray<double>(vdt));
      if(i0.check())
        flags.SetFlag(s, makeCArray<double>(vdt));
      if (s0.check())
        flags.SetFlag(s, makeCArray<string>(vdt));
    }
}

void NGS_DLL_HEADER  ExportNgstd(py::module & m) {

  py::class_<MPIManager>(m, "MPIManager")
    .def("InitMPI", &MPIManager::InitMPIB)
    .def("Barrier", &MPIManager::Barrier)
    .def("GetWT", &MPIManager::GetWT)
    .def("GetRank", &MPIManager::GetRank)
    .def("GetNP", &MPIManager::GetNP)
    ;

  std::string nested_name = "ngstd";

  m.def("TestFlags", [] (py::dict const &d) { cout << py::extract<Flags>(d)() << endl; } );

  py::class_<DummyArgument>(m, "DummyArgument")
    .def("__bool__", []( DummyArgument &self ) { return false; } )
    ;
  
  py::class_<PajeTrace >(m, "Tracer")
    .def("SetTraceThreads", &PajeTrace::SetTraceThreads)
    .def("SetTraceThreadCounter", &PajeTrace::SetTraceThreadCounter)
    .def("SetMaxTracefileSize", &PajeTrace::SetMaxTracefileSize)
    ;

  py::class_<FlatArray<double> > class_flatarrayd (m, "FlatArrayD");
  class_flatarrayd.def(py::init<int, double *>());
  PyDefVector<FlatArray<double>, double>(m, class_flatarrayd);
  PyDefToString<FlatArray<double>>(m, class_flatarrayd);
  
  py::class_<Array<double>, FlatArray<double> >(m, "ArrayD")
    .def(py::init<int>())
    .def("__init__", [](Array<double> &a, std::vector<double> const & x)
                           {
                             int s = x.size();
                             new (&a) Array<double>(s);
                             for (int i = 0; i < s; i++)
                               a[i] = x[i]; 
                           })
    .def("__rand__" ,  []( Array<double> & a, shared_ptr<Archive> & arch )
                                         { cout << "output d array" << endl;
                                           *arch & a; return arch; })
    .def("print", [](Array<double> &a) { cout << a << endl; } )
    ;

  py::class_<FlatArray<int> > class_flatarrayi (m, "FlatArrayI");
  PyDefVector<FlatArray<int>, int>(m, class_flatarrayi);
  PyDefToString<FlatArray<int> >(m, class_flatarrayi);
  class_flatarrayi.def(py::init<int, int *>());

  py::class_<Array<int>, FlatArray<int> >(m, "ArrayI")
    .def(py::init<int>())
    .def("__init__", [](std::vector<int> const & x)
                           {
                             int s = x.size();
                             shared_ptr<Array<int>> tmp (new Array<int>(s));
                             for (int i = 0; i < s; i++)
                               (*tmp)[i] = x[i]; 
                             return tmp;
                           })
    ;

  py::class_<ngstd::LocalHeap> (m, "LocalHeap", "A heap for fast memory allocation")
     .def(py::init<size_t,const char*>(), "size"_a=1000000, "name"_a="PyLocalHeap")
    ;

  py::class_<ngstd::HeapReset>
    (m, "HeapReset","stores heap-pointer on init, and resets it on exit")
    .def(py::init<LocalHeap&>())
    ;
  
  py::class_<ngstd::BitArray, shared_ptr<BitArray>> (m, "BitArray")
    .def(py::init<int>())
    .def(py::init<const BitArray&>())
    .def("__str__", &ToString<BitArray>)
    .def("__len__", &BitArray::Size)
    .def("__getitem__", [] (BitArray & self, int i) 
                                         {
                                           if (i < 0 || i >= self.Size())
                                             throw py::index_error();
                                           return self.Test(i); 
                                         })
    .def("__setitem__", [] (BitArray & self, int i, bool b) 
                                         {
                                           if (i < 0 || i >= self.Size())
                                             throw py::index_error();
                                           if (b) self.Set(i); else self.Clear(i); 
                                         })

    .def("Set", [] (BitArray & self, int i)
                                 {
                                   if (i < 0 || i >= self.Size()) 
                                     throw py::index_error();
                                   self.Set(i); 
                                 })
    .def("Clear", [] (BitArray & self, int i)
                                   {
                                   if (i < 0 || i >= self.Size()) 
                                     throw py::index_error();
                                   self.Clear(i); 
                                   })
    .def("__ior__", [] (BitArray & self, BitArray & other)
                                 {
                                   self.Or(other); 
                                   return self;
                                 })
    .def("__iand__", [] (BitArray & self, BitArray & other)
                                 {
                                   self.And(other); 
                                   return self;
                                 })
    ;

  py::class_<Flags>(m, "Flags")
    .def(py::init<>())
    .def("__str__", &ToString<Flags>)
    .def("__init__", [] (Flags &f, py::object & obj) {
         py::dict d(obj);
         new (&f) Flags();
         SetFlag(f, "", d);
         // cout << f << endl;
     })
    .def("Set",[](Flags & self,const py::dict & aflags)->Flags&
    {      
      cout << "we call Set(dict)" << endl;
      SetFlag(self, "", aflags);
      return self;
    })

    .def("Set",[](Flags & self, const char * akey, const py::object & value)->Flags&
    {             
      cout << "we call Set(key,obj)" << endl; 
        SetFlag(self, akey, value);
        return self;
    })

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
	})
  ;

  m.def("TestFlagsConversion", []( Flags flags) { cout << flags << endl; } );
  py::implicitly_convertible<py::dict, Flags>();

  py::class_<ngstd::IntRange> py_intrange (m, "IntRange");
  py_intrange.def( py::init<int,int>());
  py_intrange.def("__str__", &ToString<IntRange>);
  PyDefIterable2<IntRange>(m, py_intrange);

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
	   }
	   );

  py::class_<Archive, shared_ptr<Archive>> (m, "Archive")
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

    .def("__and__" , [](shared_ptr<Archive> & self, Array<int> & a) 
                                         { cout << "output array" << endl;
                                           *self & a; return self; })
  ;

  m.def("RunWithTaskManager", 
          [](py::object lam)
                           {
                             cout << IM(3) << "running Python function with task-manager:" << endl;
                             RunWithTaskManager ([&] () { lam(); });
                           })
          ;

  m.def("SetNumThreads", &TaskManager::SetNumThreads );

  // local TaskManager class to be used as context manager in Python
  class ParallelContextManager {
      int num_threads;
    public:
      ParallelContextManager() : num_threads(0) {};
      void Enter() {num_threads = EnterTaskManager(); }
      void Exit(py::object exc_type, py::object exc_value, py::object traceback) {
          ExitTaskManager(num_threads);
      }
    };

  py::class_<ParallelContextManager>(m, "TaskManager")
    .def(py::init<>())
    .def("__enter__", &ParallelContextManager::Enter)
    .def("__exit__", &ParallelContextManager::Exit)
    ;
}


PYBIND11_PLUGIN(libngstd) {
  py::module m("ngstd", "pybind ngstd");
  ExportNgstd(m);
  return m.ptr();
}


#endif // NGS_PYTHON
