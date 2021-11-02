#include <comp.hpp>
using namespace ngcomp;
using namespace ngfem;
#include <python_comp.hpp>

#include "myFESpace.hpp"

/*
PYBIND11_MODULE(mymodule,m) {
  // import ngsolve such that python base classes are defined
  auto ngs = py::module::import("ngsolve");
  
  ExportMyFESpace(m);
}
*/

   
extern "C" void mymodule(py::object & res) {
  cout << "called mymodule" << endl;
  // import ngsolve such that python base classes are defined    
  auto ngs = py::module::import("ngsolve");    
  static py::module::module_def def;    
  py::module m = py::module::create_extension_module("", "", &def);    
    
  ExportMyFESpace(m);                                                                         
  res = m;    
}    

