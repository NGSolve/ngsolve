#ifdef NGS_PYTHON

// #include <string>
// #include <ostream>
// #include <type_traits>

#include "python_ngstd.hpp"

using std::ostringstream;


BOOST_PYTHON_MODULE(Ngstd) 
{
  cout << "init ngstd - py" << endl;

  bp::class_<FlatArray<double> >("FlatArrayD")
    .def(PyDefVector<FlatArray<double>, double>()) 
    .def(PyDefToString<FlatArray<double> >())
    .def(bp::init<int, double *>())
    ;
    
  bp::class_<Array<double>, bp::bases<FlatArray<double> > >("ArrayD")
    .def(bp::init<int>())
    ;

  bp::class_<FlatArray<int> >("FlatArrayI")
    .def(PyDefVector<FlatArray<int>, int>()) 
    .def(PyDefToString<FlatArray<int> >())
    .def(bp::init<int, int *>())
    ;

  bp::class_<Array<int>, bp::bases<FlatArray<int> > >("ArrayI")
    .def(bp::init<int>())
    ;
    
  bp::class_<ngstd::LocalHeap>
    ("LocalHeap",bp::init<size_t,const char*>())
    ;

  bp::class_<ngstd::HeapReset>
    ("HeapReset",bp::init<LocalHeap&>())
    ;

  bp::class_<ngstd::Flags>
    ("Flags")
    ;

  bp::class_<ngstd::IntRange>
    ("IntRange", bp::init<int,int>())
    .def(PyDefToString<IntRange>())
    .def(PyDefList<IntRange,int>())
    ;

    
}



BOOST_PYTHON_MODULE(libngstd)
{
  cout << "running pyinit_libngstd" << endl;
  PyImport_AppendInittab("Ngstd", PyInit_Ngstd);  
}

#endif // NGS_PYTHON
