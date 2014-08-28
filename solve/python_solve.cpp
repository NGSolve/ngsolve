#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <solve.hpp>
using namespace ngsolve;





void ExportNgsolve() {
    std::string nested_name = "ngsolve";
    if( bp::scope() )
      nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".ngsolve");
    
    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << "exporting ngstd as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("ngsolve") = module ;

    bp::scope local_scope(module);



  PyExportSymbolTable<shared_ptr<FESpace>> ();
  PyExportSymbolTable<shared_ptr<GridFunction>> ();
  PyExportSymbolTable<shared_ptr<BilinearForm>> ();
  PyExportSymbolTable<shared_ptr<LinearForm>> ();
  PyExportSymbolTable<shared_ptr<NumProc>> ();
  PyExportSymbolTable<double> ();
  PyExportSymbolTable<shared_ptr<double>> ();


  bp::class_<NumProc, shared_ptr<NumProc>,boost::noncopyable> ("NumProc", bp::no_init)
    // .def("__str__", &ToString<NumProc>)
    .def ("__str__", 
          FunctionPointer([](NumProc & np) -> string
                          {
                            stringstream str;
                            np.PrintReport (str);
                            return str.str();
                          }));
  ;

  // bp::class_<NumProc>("NumProc", bp::no_init);

  bp::class_<PDE> ("PDE", bp::init<>())
    .def(bp::init<const string&>())
    .def("Load", static_cast<void(ngsolve::PDE::*)(const string &, const bool, const bool)> 
         (&ngsolve::PDE::LoadPDE),
         (boost::python::arg("filename"), 
          boost::python::arg("meshload")=0, 
          boost::python::arg("nogeometryload")=0))
    .def("Mesh",  static_cast<MeshAccess&(ngsolve::PDE::* const)(int)>(&PDE::GetMeshAccess),
         bp::return_value_policy<bp::reference_existing_object>(),
         (bp::arg("nr")=0))
    .def("Solve", &ngsolve::PDE::Solve)
    /*
    .def("__setitem__", FunctionPointer([](PDE & self, string name, shared_ptr<GridFunction> gf)
                                        {
                                          self.AddGridFunction (name, gf);
                                        }))
    */
    .def("Add", FunctionPointer([](PDE & self, shared_ptr<GridFunction> gf)
                                {
                                  self.AddGridFunction (gf->GetName(), gf);
                                }))

    .add_property ("constants", FunctionPointer([](PDE & self) { return self.GetConstantTable(); }))
    .add_property ("variables", FunctionPointer([](PDE & self) { return self.GetVariableTable(); }))
    .add_property ("spaces", FunctionPointer([](PDE & self) { return self.GetSpaceTable(); }))
    .add_property ("gridfunctions", FunctionPointer([](PDE & self) { return self.GetGridFunctionTable(); }))
    .add_property ("bilinearforms", FunctionPointer([](PDE & self) { return self.GetBilinearFormTable(); }))
    .add_property ("linearforms", FunctionPointer([](PDE & self) { return self.GetLinearFormTable(); }))
    .add_property ("numprocs", FunctionPointer([](PDE & self) { return self.GetNumProcTable(); }))
    ;
}




#endif
