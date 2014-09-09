#ifdef NGS_PYTHON

// #include <string>
// #include <ostream>
// #include <type_traits>

#include "python_ngstd.hpp"

using std::ostringstream;

template<typename T>
Array<T> & makeCArray(const bp::object & obj)
{     
    Array<T> * C_vdL = new Array<T>(bp::len(obj));    
    for (int i = 0; i < bp::len(obj); i++)    
        (*C_vdL)[i] = bp::extract<T>(obj[i]);        
    return *C_vdL;
}


void SetFlag(Flags &flags, const char * s, bp::object value) {
    bp::extract<double> vd(value);
    if (vd.check())
      {
        // cout << "is double" << endl;
        flags.SetFlag(s, vd());         
      }

    bp::extract<char *> vs(value);
    if (vs.check())
      {
        // cout << "is char*" << endl;
        flags.SetFlag(s, vs());
      }

    bp::extract<bp::list> vdl(value);
    if (vdl.check())
    {             
      // cout << "is list" << endl;
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
      // cout << "is tuple" << endl;
        bp::extract<double> d0(vdt()[0]);
        bp::extract<char *> s0(vdt()[0]);
        if (d0.check())
            flags.SetFlag(s, makeCArray<double>(value));
        if (s0.check())
            flags.SetFlag(s, makeCArray<string>(value));
    }
}

struct FlagsFromPythonDict{
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


void ExportNgstd() {
    std::string nested_name = "ngstd";
    if( bp::scope() )
      nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".ngstd");
    
    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << "exporting ngstd as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("ngstd") = module ;

    bp::scope ngbla_scope(module);


    bp::def("TestFlags", FunctionPointer( [] (bp::dict const &d) { cout << bp::extract<Flags>(d)() << endl; } ) );

  bp::class_<FlatArray<double> >("FlatArrayD")
    .def(PyDefVector<FlatArray<double>, double>()) 
    .def(PyDefToString<FlatArray<double> >())
    .def(bp::init<int, double *>())
    ;
    
  bp::class_<Array<double>, bp::bases<FlatArray<double> > >("ArrayD")
    .def(bp::init<int>())
    .def("__init__", bp::make_constructor (FunctionPointer ([](bp::list const & x)
                                                            {
                  int s = bp::len(x);
                  shared_ptr<Array<double>> tmp (new Array<double>(s));
                  for (int i = 0; i < s; i++)
                    (*tmp)[i] = bp::extract<double> (x[i]); 
                  return tmp;
                })))
    ;

  bp::class_<FlatArray<int> >("FlatArrayI")
    .def(PyDefVector<FlatArray<int>, int>()) 
    .def(PyDefToString<FlatArray<int> >())
    .def(bp::init<int, int *>())
    ;

  bp::class_<Array<int>, bp::bases<FlatArray<int> > >("ArrayI")
    .def(bp::init<int>())
    .def("__init__", bp::make_constructor (FunctionPointer ([](bp::list const & x)
                {
                  int s = bp::len(x);
                  shared_ptr<Array<int>> tmp (new Array<int>(s));
                  for (int i = 0; i < s; i++)
                    (*tmp)[i] = bp::extract<int> (x[i]); 
                  return tmp;
                })))
    ;
    
  bp::class_<ngstd::LocalHeap>
    ("LocalHeap",bp::init<size_t,const char*>()
     //,(bp::arg("self"), bp::arg("size")=1000000, bp::arg("name")="PyLocalHeap")
     )
    ;

  bp::class_<ngstd::HeapReset>
    ("HeapReset",bp::init<LocalHeap&>())
    // .def(bp::init<const HeapReset&>())
    // .def("__enter__", FunctionPointer([](HeapReset & lh) { cout << "enter" << endl; }))
    // .def("__exit__", FunctionPointer([](HeapReset & lh, bp::object x, bp::object y, bp::object z) { cout << "exit" << endl; }))    
    ;

  bp::class_<ngstd::BitArray>
    ("BitArray")
    .def("__str__", &ToString<BitArray>)
  ;

  bp::class_<ngstd::Flags, shared_ptr<Flags>, boost::noncopyable> ("Flags", bp::no_init)
    .def("__init__", bp::make_constructor (FunctionPointer ([](bp::dict const & aflags) 
                                                            {
      cout << "Calling Flags constructor with dict input" << endl;
      shared_ptr<Flags> self = make_shared<Flags>();
      for (int i = 0; i < bp::len(aflags); i++)
      {   
            char * s = bp::extract<char *>(aflags.keys()[i]);          
            SetFlag(*self, s, aflags.values()[i]);
      }
      cout << "flags from dict: " << endl << *self << endl;
      return self;
                })))

    .def("Set", FunctionPointer([](Flags & self,bp::dict const & aflags)->Flags&
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

    .def("__str__", &ToString<Flags>)   
   ;


  bp::class_<ngstd::IntRange>
    ("IntRange", bp::init<int,int>())
    // .def(PyDefIterable<IntRange,int>())
    .def(PyDefIterable2<IntRange>())
    .def("__str__", &ToString<IntRange>)
    ;

    FlagsFromPythonDict();
}



BOOST_PYTHON_MODULE(libngstd) {
  ExportNgstd();
}




#endif // NGS_PYTHON
