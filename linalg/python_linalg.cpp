#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <la.hpp>
using namespace ngla;






BOOST_PYTHON_MODULE(libngla) {

  cout << "init ngla - py" << endl;
  
  bp::class_<BaseVector,boost::noncopyable>("BaseVector", bp::no_init)
    .def("__str__", &ToString<BaseVector>)
    .add_property("size", &BaseVector::Size)
    .def("CreateVector", &BaseVector::CreateVector,
         bp::return_value_policy<bp::manage_new_object>())

    .def("Assign", FunctionPointer([](BaseVector & self, BaseVector & v2, double s)->void { self.Set(s, v2); }))
    .def("Add", FunctionPointer([](BaseVector & self, BaseVector & v2, double s)->void { self.Add(s, v2); }))
    ;
  bp::def("InnerProduct", FunctionPointer([](BaseVector & v1, BaseVector & v2)->double { return InnerProduct(v1,v2); }))
    ;
  

  typedef BaseMatrix BM;
  typedef BaseVector BV;

  bp::class_<BaseMatrix,boost::noncopyable>("BaseMatrix", bp::no_init)
    .def("__str__", &ToString<BaseMatrix>)
    .add_property("height", &BaseMatrix::Height)
    .add_property("width", &BaseMatrix::Width)
    
    .def("Mult",        FunctionPointer( [](BM &m, BV &x, BV &y, double s) { m.Mult (x,y); y *= s; }) )
    .def("MultAdd",     FunctionPointer( [](BM &m, BV &x, BV &y, double s) { m.MultAdd (s, x, y); }))
    // .def("MultTrans",   FunctionPointer( [](BM &m, BV &x, BV &y, double s) { y  = s*Trans(m)*x; }) )
    // .def("MultTransAdd",FunctionPointer( [](BM &m, BV &x, BV &y, double s) { y += s*Trans(m)*x; }) )

    .def("Inverse", FunctionPointer( [](BM &m) { return m.InverseMatrix(); }),
         bp::return_value_policy<bp::manage_new_object>())
    ;
  
}




struct Init {
  Init() 
  { 
    cout << "adding module 'ngla' to py-inittab" << endl;
    PyImport_AppendInittab("ngla", PyInit_libngla); 
  }
};
static Init init;







#endif // NGS_PYTHON
