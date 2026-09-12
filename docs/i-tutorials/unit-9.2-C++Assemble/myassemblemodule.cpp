#include <comp.hpp>
#include <python_comp.hpp>

#include "myIntegrator.cpp"
#include "myAssembling.cpp"

   
// NGCORE_API_EXPORT: gcc gives pybind11 types hidden visibility, and passes it
// on to this function - without it the symbol is not exported and not found
extern "C" NGCORE_API_EXPORT void mymodule(py::object & res) {
  // import ngsolve such that python base classes are defined    
  auto ngs = py::module::import("ngsolve");    

  static py::module::module_def def;    
  py::module m = py::module::create_extension_module("", "", &def);    
    
  // ExportMyFESpace(m);

  using namespace ngfem;
  py::class_<MyLaplaceIntegrator, shared_ptr<MyLaplaceIntegrator>, BilinearFormIntegrator>
    (m, "MyLaplace")
    .def(py::init<shared_ptr<CoefficientFunction>>())
    ;
  py::class_<MySourceIntegrator, shared_ptr<MySourceIntegrator>, LinearFormIntegrator>
    (m, "MySource")
    .def(py::init<shared_ptr<CoefficientFunction>>())
    ;
  
  m.def ("MyAssembleMatrix",
         &ngcomp::MyAssembleMatrix,
         py::arg("fes"),py::arg("integrator"), py::arg("parallel")=false);
  
  res = m;    
}    

