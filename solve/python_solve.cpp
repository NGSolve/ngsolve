#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <solve.hpp>
using namespace ngcomp;

// extern AutoPtr<ngsolve::PDE> pde;


void PyExportNgSolve(BasePythonEnvironment & py_env) 
{
  PyExportSymbolTable< FESpace* > ();
  PyExportSymbolTable< GridFunction * > ();
  PyExportSymbolTable< BilinearForm * > ();
  PyExportSymbolTableStdTypes< double > ();
  PyExportSymbolTableStdTypes< double* > ();


  py_env["PDE"] = bp::class_<PDE> ("PDE")
    // .def("SetSelfPtr", FunctionPointer([](PDE & self) { pde.Reset (&self, false); return; }))
    .def("Load", static_cast<void(ngsolve::PDE::*)(const string &, const bool, const bool)> 
         (&ngsolve::PDE::LoadPDE),
         (boost::python::arg("filename"), 
          boost::python::arg("meshload")=0, 
          boost::python::arg("nogeometryload")=0))
    .def("Solve", &ngsolve::PDE::Solve)
    .def("Spaces", &PDE::GetSpaceTable, bp::return_value_policy<bp::reference_existing_object>())
    .def("Variables", &PDE::GetVariableTable, bp::return_value_policy<bp::reference_existing_object>())
    .def("Constants", &PDE::GetConstantTable, bp::return_value_policy<bp::reference_existing_object>())
    .add_property ("spaces", FunctionPointer([](PDE & self) { return PyRef<SymbolTable<FESpace*>>(self.GetSpaceTable()); }))
    .add_property ("gridfunctions", FunctionPointer([](PDE & self) { return PyRef<SymbolTable<GridFunction*>>(self.GetGridFunctionTable()); }))
    .add_property ("bilinearforms", FunctionPointer([](PDE & self) { return PyRef<SymbolTable<BilinearForm*>>(self.GetBilinearFormTable()); }))
    ;

}


#endif
