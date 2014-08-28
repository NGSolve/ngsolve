#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <la.hpp>
using namespace ngla;





void ExportNgla() {
    std::string nested_name = "ngla";
    if( bp::scope() )
         nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".ngla");
    
    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << "exporting ngla as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("ngla") = module ;

    bp::scope ngla_scope(module);

    bp::object expr_module = bp::import("expr");
    bp::object expr_namespace = expr_module.attr("__dict__");





  
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



int ExportNgbla();

BOOST_PYTHON_MODULE(libngbla) {
  ExportNgbla();
  ExportNgla();
}






#endif // NGS_PYTHON
