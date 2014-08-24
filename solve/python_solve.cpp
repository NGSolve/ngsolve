#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <solve.hpp>
using namespace ngcomp;

template <typename T> 
auto_ptr<T> make_auto_ptr(T * ptr)
{
  return auto_ptr<T> (ptr);
}

BOOST_PYTHON_MODULE(Ngsolve) {
  cout << "init py - ngsolve" << endl;
  PyExportSymbolTable<shared_ptr<FESpace>> ();
  PyExportSymbolTable<shared_ptr<GridFunction>> ();
  PyExportSymbolTable<shared_ptr<BilinearForm>> ();
  PyExportSymbolTableStdTypes< double > ();
  PyExportSymbolTableStdTypes< double* > ();


  bp::class_<PDE> ("PDE")
    .def("Load", static_cast<void(ngsolve::PDE::*)(const string &, const bool, const bool)> 
         (&ngsolve::PDE::LoadPDE),
         (boost::python::arg("filename"), 
          boost::python::arg("meshload")=0, 
          boost::python::arg("nogeometryload")=0))
    .def("Mesh",  static_cast<MeshAccess&(ngsolve::PDE::* const)(int)>(&PDE::GetMeshAccess),
         bp::return_value_policy<bp::reference_existing_object>(),
         (bp::arg("nr")=0))
    .def("Solve", &ngsolve::PDE::Solve)
    .def("Spaces", &PDE::GetSpaceTable, bp::return_value_policy<bp::reference_existing_object>())
    .def("Variables", &PDE::GetVariableTable, bp::return_value_policy<bp::reference_existing_object>())
    .def("Constants", &PDE::GetConstantTable, bp::return_value_policy<bp::reference_existing_object>())
    // .add_property ("spaces", FunctionPointer([](PDE & self) { return PyRef<SymbolTable<shared_ptr<FESpace>>>(self.GetSpaceTable()); }))
    .add_property ("spaces", FunctionPointer([](PDE & self) { return self.GetSpaceTable(); }))
    .add_property ("gridfunctions", FunctionPointer([](PDE & self) { return self.GetGridFunctionTable(); }))
    .add_property ("bilinearforms", FunctionPointer([](PDE & self) { return self.GetBilinearFormTable(); }))
    ;
}


#endif
