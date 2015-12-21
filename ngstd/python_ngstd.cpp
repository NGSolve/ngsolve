#ifdef NGS_PYTHON

// #include <string>
// #include <ostream>
// #include <type_traits>

#include "python_ngstd.hpp"

PythonEnvironment pyenv;


using std::ostringstream;


void SetFlag(Flags &flags, const char * s, bp::object value) 
{
  bp::extract<bool> vb(value);
  if (vb.check() && vb())
    flags.SetFlag(s);

  bp::extract<double> vd(value);
  if (vd.check())
    flags.SetFlag(s, vd());         

  bp::extract<char *> vs(value);
  if (vs.check())
    flags.SetFlag(s, vs());

  bp::extract<bp::list> vdl(value);
  if (vdl.check())
    {             
      bp::extract<double> d0(vdl()[0]);
      bp::extract<char *> s0(vdl()[0]);
      if(d0.check())
        flags.SetFlag(s, makeCArray<double>(value));
      if (s0.check())
        flags.SetFlag(s, makeCArray<string>(value));
    }

  bp::extract<bp::tuple> vdt(value);
  if (vdt.check())
    {
      bp::extract<double> d0(vdt()[0]);
      bp::extract<char *> s0(vdt()[0]);
      if (d0.check())
        flags.SetFlag(s, makeCArray<double>(value));
      if (s0.check())
        flags.SetFlag(s, makeCArray<string>(value));
    }
}

struct FlagsFromPythonDict 
{
  static void* convertible(PyObject* obj_ptr) {
    if (PyMapping_Check(obj_ptr)) {
      return obj_ptr;
    } else {
      return NULL;
    }
  }
  
  static void construct(PyObject* obj_ptr,
                        bp::converter::rvalue_from_python_stage1_data* data) {
    bp::dict aflags(bp::handle<>(bp::borrowed(obj_ptr)));
    Flags self;
    // cout << "py converter" << endl;
    for (int i=0; i<bp::len(aflags); i++) {
      char * s = bp::extract<char *>(aflags.keys()[i]);          
      SetFlag(self, s, aflags.values()[i]);
    }
    typedef bp::converter::rvalue_from_python_storage<Flags> storage_t;
    storage_t* the_storage = reinterpret_cast<storage_t*>(data);
    void* memory_chunk = the_storage->storage.bytes;
    /* Flags* output = */ new (memory_chunk) Flags(std::move(self) ); // Use the contents of obj to populate output, e.g. using extract<>
    data->convertible = memory_chunk;
  }
  
  FlagsFromPythonDict() {
    bp::converter::registry::push_back(&convertible, &construct, bp::type_id<Flags>());
  }
};


void TranslateException (const Exception & ex)
{
  string err = string("NGSolve exception: ")+ex.What();
  PyErr_SetString(PyExc_RuntimeError, err.c_str());
}


void NGS_DLL_HEADER  ExportNgstd() {
  std::string nested_name = "ngstd";
  if( bp::scope() )
    nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".ngstd");
  
  bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
  
  cout << IM(1) << "exporting ngstd as " << nested_name << endl;
  bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
  parent.attr("ngstd") = module ;
  
  bp::scope ngbla_scope(module);

  bp::def("TestFlags", FunctionPointer( [] (bp::dict const &d) { cout << bp::extract<Flags>(d)() << endl; } ) );
  

  bp::register_exception_translator<Exception>(&TranslateException);

  bp::class_<PajeTrace >("Tracer", bp::no_init)
    .def("SetTraceThreads", &PajeTrace::SetTraceThreads)
    .def("SetTraceThreadCounter", &PajeTrace::SetTraceThreadCounter)
    .def("SetMaxTracefileSize", &PajeTrace::SetMaxTracefileSize)
    ;
  bp::class_<FlatArray<double> >("FlatArrayD")
    .def(PyDefVector<FlatArray<double>, double>()) 
    .def(PyDefToString<FlatArray<double> >())
    .def(bp::init<int, double *>())
    ;
  
  bp::class_<Array<double>, bp::bases<FlatArray<double> > >("ArrayD")
    .def(bp::init<int>())
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](bp::list const & x)
                           {
                             int s = bp::len(x);
                             auto tmp = make_shared<Array<double>> (s);
                             for (int i = 0; i < s; i++)
                               (*tmp)[i] = bp::extract<double> (x[i]); 
                             return tmp;
                           })))
    .def("__rand__" , FunctionPointer( []( Array<double> & a, shared_ptr<Archive> & arch )
                                         { cout << "output d array" << endl;
                                           *arch & a; return arch; }));
    ;

  bp::class_<FlatArray<int> >("FlatArrayI")
    .def(PyDefVector<FlatArray<int>, int>()) 
    .def(PyDefToString<FlatArray<int> >())
    .def(bp::init<int, int *>())
    ;

  bp::class_<Array<int>, bp::bases<FlatArray<int> > >("ArrayI")
    .def(bp::init<int>())
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](bp::list const & x)
                           {
                             int s = bp::len(x);
                             shared_ptr<Array<int>> tmp (new Array<int>(s));
                             for (int i = 0; i < s; i++)
                               (*tmp)[i] = bp::extract<int> (x[i]); 
                             return tmp;
                           })))
    ;
  
  bp::class_<ngstd::LocalHeap, boost::noncopyable> 
    ("LocalHeap", 
     bp::init<size_t,const char*>
     ((bp::arg("self"), bp::arg("size")=1000000, bp::arg("name")="PyLocalHeap"),
      "A heap for fast memory allocation"))
    ;

  bp::class_<ngstd::HeapReset>
    ("HeapReset",bp::init<LocalHeap&>("stores heap-pointer on init, and resets it on exit"))
    // .def(bp::init<const HeapReset&>())
    // .def("__enter__", FunctionPointer([](HeapReset & lh) { cout << "enter" << endl; }))
    // .def("__exit__", FunctionPointer([](HeapReset & lh, bp::object x, bp::object y, bp::object z) { cout << "exit" << endl; }))    
    ;

  bp::class_<ngstd::BitArray> ("BitArray")
    .def(bp::init<int>())
    .def(bp::init<const BitArray&>())
    .def("__str__", &ToString<BitArray>)
    .def("__len__", &BitArray::Size)
    .def("__getitem__", FunctionPointer ([] (BitArray & self, int i) 
                                         {
                                           if (i < 0 || i >= self.Size())
                                             bp::exec("raise IndexError()\n");
                                           return self.Test(i); 
                                         }))
    .def("__setitem__", FunctionPointer ([] (BitArray & self, int i, bool b) 
                                         {
                                           if (i < 0 || i >= self.Size())
                                             bp::exec("raise IndexError()\n");
                                           if (b) self.Set(i); else self.Clear(i); 
                                         }))

    .def("Set", FunctionPointer ([] (BitArray & self, int i)
                                 {
                                   if (i < 0 || i >= self.Size()) 
                                     bp::exec("raise IndexError()\n");
                                   self.Set(i); 
                                 }))
    .def("Clear", FunctionPointer ([] (BitArray & self, int i)
                                   {
                                   if (i < 0 || i >= self.Size()) 
                                     bp::exec("raise IndexError()\n");
                                   self.Clear(i); 
                                   }))
    .def("__ior__", FunctionPointer ([] (BitArray & self, BitArray & other)
                                 {
                                   self.Or(other); 
                                   return self;
                                 }))
    .def("__iand__", FunctionPointer ([] (BitArray & self, BitArray & other)
                                 {
                                   self.And(other); 
                                   return self;
                                 }))
    ;

  bp::class_<ngstd::Flags, shared_ptr<Flags>, boost::noncopyable> ("Flags", bp::no_init)
    .def("__init__", bp::make_constructor (FunctionPointer ([](const bp::dict & aflags) 
                                                            {
      shared_ptr<Flags> self = make_shared<Flags>();
      for (int i = 0; i < bp::len(aflags); i++)
      {   
            char * s = bp::extract<char *>(aflags.keys()[i]);          
            SetFlag(*self, s, aflags.values()[i]);
      }
      return self;
                })))

    .def("Set", FunctionPointer([](Flags & self,const bp::dict & aflags)->Flags&
    {      
      cout << "we call Set(dict)" << endl;
      for (int i = 0; i < bp::len(aflags); i++)
      {   
          char * s = bp::extract<char *>(aflags.keys()[i]);          
          SetFlag(self, s, aflags.values()[i]);
      }
      return self;
    }), bp::return_value_policy<bp::reference_existing_object>())

    .def("Set", FunctionPointer([](Flags & self, const char * akey, const bp::object & value)->Flags&
    {             
      cout << "we call Set(key,obj)" << endl; 
        SetFlag(self, akey, value);
        return self;
    }), bp::return_value_policy<bp::reference_existing_object>()       
        )

    .def("__getitem__", FunctionPointer( [](Flags & self, const string& name) -> bp::object {

	  if(self.NumListFlagDefined(name))
	    return bp::object(self.GetNumListFlag(name));

	  if(self.StringListFlagDefined(name))
	    return bp::object(self.GetStringListFlag(name));
	 
	  if(self.NumFlagDefined(name))
	    return bp::object(*self.GetNumFlagPtr(name));
	  
	  if(self.StringFlagDefined(name))
	    return bp::object(self.GetStringFlag(name));

	  if(self.FlagsFlagDefined(name))
	    return bp::object(self.GetFlagsFlag(name));

	  return bp::object(self.GetDefineFlag(name));
	}))


    .def("__str__", &ToString<Flags>)   
   ;


  bp::class_<ngstd::IntRange>
    ("IntRange", bp::init<int,int>())
    // .def(PyDefIterable<IntRange,int>())
    .def(PyDefIterable2<IntRange>())
    .def("__str__", &ToString<IntRange>)
    ;
  
  FlagsFromPythonDict();

  bp::def("Timers", FunctionPointer
	  ([]() 
	   {
	     bp::dict timers;
	     for (int i = 0; i < NgProfiler::SIZE; i++)
	       if (!NgProfiler::names[i].empty())
		 timers[NgProfiler::names[i]] = NgProfiler::GetTime(i);
	     return timers;
	   }
	   ));
  
  
  FlagsFromPythonDict();

  bp::class_<Archive, shared_ptr<Archive>, boost::noncopyable> ("Archive", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](const string & filename, bool write,
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
                           })))

    .def("__and__" , FunctionPointer( [](shared_ptr<Archive> & self, Array<int> & a) 
                                         { cout << "output array" << endl;
                                           *self & a; return self; }))
  ;

  
  bp::def("RunWithTaskManager", 
          FunctionPointer ([](bp::object lam)
                           {
                             cout << IM(3) << "running Python function with task-manager:" << endl;
                             RunWithTaskManager ([&] () { lam(); });
                           }))
          ;
}


BOOST_PYTHON_MODULE(libngstd) {
  ExportNgstd();
}




#endif // NGS_PYTHON
