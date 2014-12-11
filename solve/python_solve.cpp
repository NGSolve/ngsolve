#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <solve.hpp>
using namespace ngsolve;



extern void ExportBVP();

void ExportNgsolve() {
    std::string nested_name = "solve";
    if( bp::scope() )
      nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".solve");
    
    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << "exporting solve as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("solve") = module ;

    bp::scope local_scope(module);

    bp::def ("Redraw", 
             FunctionPointer([](bool blocking) {Ng_Redraw(blocking);}),
             (bp::arg("blocking")=false)
             );

    ExportBVP();
}



BOOST_PYTHON_MODULE(libsolve) {
  ExportNgsolve();
}


#endif
