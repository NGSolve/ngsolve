#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <solve.hpp>
using namespace ngsolve;

class NumProc1 : public NumProc
{
public:
  NumProc1 (PDE & pde, const Flags & flags) : NumProc (pde, flags)
  { ; }
  // virtual void Do (LocalHeap & lh) { ; }
};

class NumProcWrap : public NumProc1, public bp::wrapper<NumProc1> {
public:
  NumProcWrap (PDE & pde, const Flags & flags) : NumProc1(pde, flags) { ; }
  void Do(LocalHeap & lh)  {
    this->get_override("Do")(lh);
  }
};


// crashes 
class PyNumProc : public NumProc, public bp::wrapper<NumProc> {
public:
  PyNumProc (PDE & apde, const Flags & flags) : NumProc(pde,flags) { ; }

  void Do(LocalHeap & lh) {
    this->get_override("Do")(lh);
  }
};

static PyNumProc *python_np;

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

  // TestExport();

  bp::def("RegisterNumProc", FunctionPointer( [] (PyNumProc & np, string label, int dim) { 
      python_np = &np;
      GetNumProcs().AddNumProc (label, dim, FunctionPointer( [] (PDE & pde, const Flags & flags) { return shared_ptr<NumProc>(python_np);}  ));
//       RegisterNumProc<PyNumProc>(name, dim); return;
      } ));
      

  bp::class_<NumProc, shared_ptr<NumProc>, boost::noncopyable> ("NumProc", bp::no_init)
    ;

  // die geht
  bp::class_<NumProcWrap,bp::bases<NumProc>,boost::noncopyable>("NumProc1", bp::init<PDE&, const Flags&>())
    .def("Do", bp::pure_virtual(&NumProc1::Do)) 
    ;
  

  // die geht nicht
  bp::class_<PyNumProc,boost::noncopyable> ("NumProc2", bp::init<PDE&, const Flags&>())
    .def("Do", bp::pure_virtual(&NumProc::Do))
    .def ("__str__", 
          FunctionPointer([](NumProc & np) -> string
                          {
                            stringstream str;
                            np.PrintReport (str);
                            return str.str();
                          }));
  ;

  
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

    .def("Add", FunctionPointer([](PDE & self, PyNumProc &np)
                                {
                                  cout << "add pynumproc - ref" << endl;
                                  self.AddNumProc ("pynumproc", shared_ptr<NumProc> (&np, &NOOP_Deleter));
                                }))
    
    .def("Add", FunctionPointer([](PDE & self, NumProc1 & np)
                                {
                                  cout << "add numproc1 - at addr " << &np << endl;
                                  self.AddNumProc ("pynumproc", shared_ptr<NumProc> (&np, &NOOP_Deleter));
                                }))

    
    .def("Add", FunctionPointer([](PDE & self, shared_ptr<NumProc> np)
                                {
                                  cout << "add numproc - sp" << endl;
                                  self.AddNumProc ("pynumproc", np);
                                }))


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
