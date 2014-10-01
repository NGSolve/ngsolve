#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <solve.hpp>
using namespace ngsolve;


// static void NOOP_Deleter(void *) { ; }

class PyNumProc : public NumProc
{

public:
  PyNumProc (PDE & pde, const Flags & flags) : NumProc (pde, flags) { ; }
  shared_ptr<PDE> GetPDE() const { return shared_ptr<PDE> (&pde,&NOOP_Deleter); }
  // virtual void Do (LocalHeap & lh) { ; }
};

class NumProcWrap : public PyNumProc, public bp::wrapper<PyNumProc> {
public:
  NumProcWrap (PDE & pde, const Flags & flags) : PyNumProc(pde, flags) { ; }
  void Do(LocalHeap & lh)  {
    cout << "numproc wrap - do" << endl;
    this->get_override("Do")(lh);
  }
};


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
  PyExportSymbolTable<shared_ptr<CoefficientFunction>> ();
  PyExportSymbolTable<shared_ptr<GridFunction>> ();
  PyExportSymbolTable<shared_ptr<BilinearForm>> ();
  PyExportSymbolTable<shared_ptr<LinearForm>> ();
  PyExportSymbolTable<shared_ptr<Preconditioner>> ();
  PyExportSymbolTable<shared_ptr<NumProc>> ();
  PyExportSymbolTable<double> ();
  PyExportSymbolTable<shared_ptr<double>> ();


  bp::class_<NumProc, shared_ptr<NumProc>,bp::bases<NGS_Object>,boost::noncopyable> ("NumProc", bp::no_init)
    ;

  // die geht
  bp::class_<NumProcWrap,bp::bases<NumProc>,boost::noncopyable>("PyNumProc", bp::init<PDE&, const Flags&>())
    .def("Do", bp::pure_virtual(&PyNumProc::Do)) 
    .add_property("pde", &PyNumProc::GetPDE)
    ;
  

  
  bp::class_<PDE,shared_ptr<PDE>> ("PDE", bp::init<>())
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
    
//     .def("Add", FunctionPointer([](PDE & self, NumProc1 & np)
//                                 {
//                                   cout << "add numproc1 - at addr " << &np << endl;
//                                   self.AddNumProc ("pynumproc", shared_ptr<NumProc> (&np, &NOOP_Deleter));
//                                 }))

    
    .def("Add", FunctionPointer([](PDE & self, shared_ptr<NumProc> np)
                                {
                                  cout << "add numproc - sp" << endl;
                                  self.AddNumProc ("pynumproc", np);
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<FESpace> space)
                                {
                                  self.AddFESpace (space->GetName(), space);
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<GridFunction> gf)
                                {
                                  self.AddGridFunction (gf->GetName(), gf);
                                }))

    .def("Add", FunctionPointer([](PDE & self, shared_ptr<BilinearForm> bf)
                                {
                                  self.AddBilinearForm (bf->GetName(), bf);
                                }))

    .def("Add", FunctionPointer([](PDE & self, const bp::list &l)
                                {
                                    for (int i=0; i<bp::len(l); i++)
                                    {
                                        bp::extract<shared_ptr<PyNumProc>> np(l[i]);
                                        if(np.check())
                                        {
                                            self.AddNumProc (np()->GetName(), np());
                                            continue;
                                        }

                                        bp::extract<shared_ptr<NumProc>> pnp(l[i]);
                                        if(np.check())
                                        {
                                            self.AddNumProc (pnp()->GetName(), pnp());
                                            continue;
                                        }

                                        bp::extract<shared_ptr<GridFunction>> gf(l[i]);
                                        if(gf.check())
                                        {
                                            self.AddGridFunction (gf()->GetName(), gf());
                                            continue;
                                        }

                                        cout << "warning: unknown object at position " << i << endl;
                                    }
                                }))

    .add_property ("constants", FunctionPointer([](PDE & self) { return bp::object(self.GetConstantTable()); }))
    .add_property ("variables", FunctionPointer([](PDE & self) { return bp::object(self.GetVariableTable()); }))
    .add_property ("coefficients", FunctionPointer([](PDE & self) { return bp::object(self.GetCoefficientTable()); }))
    .add_property ("spaces", FunctionPointer([](PDE & self) { return bp::object(self.GetSpaceTable()); }))
    .add_property ("gridfunctions", FunctionPointer([](PDE & self) { return bp::object(self.GetGridFunctionTable()); }))
    .add_property ("bilinearforms", FunctionPointer([](PDE & self) { return bp::object(self.GetBilinearFormTable()); }))
    .add_property ("linearforms", FunctionPointer([](PDE & self) { return bp::object(self.GetLinearFormTable()); }))
    .add_property ("preconditioners", FunctionPointer([](PDE & self) { return bp::object(self.GetPreconditionerTable()); }))
    .add_property ("numprocs", FunctionPointer([](PDE & self) { return bp::object(self.GetNumProcTable()); }))
    ;
}

/*
BOOST_PYTHON_MODULE(libngsolve) {
  ExportNgsolve();
}
*/


#endif
