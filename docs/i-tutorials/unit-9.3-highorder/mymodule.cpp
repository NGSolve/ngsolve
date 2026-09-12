#include <comp.hpp>
#include <python_comp.hpp>

// using namespace ngcomp;

#include "myHOElement.cpp"
#include "myHOFESpace.cpp"

   
// NGCORE_API_EXPORT: gcc gives pybind11 types hidden visibility, and passes it
// on to this function - without it the symbol is not exported and not found
extern "C" NGCORE_API_EXPORT void mymodule(py::object & res) {
  cout << "called mymodule" << endl;
  // import ngsolve such that python base classes are defined    
  auto ngs = py::module::import("ngsolve");    

  static py::module::module_def def;    
  py::module m = py::module::create_extension_module("", "", &def);    

  using namespace ngcomp;
  ExportFESpace<MyHighOrderFESpace>(m, "MyHighOrderFESpace");

  
  res = m;    
}    

