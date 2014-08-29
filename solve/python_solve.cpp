#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <solve.hpp>
using namespace ngsolve;

struct TestBase {
    virtual void f() const = 0;
    virtual void g() const { cout << "g" << endl; };
    virtual void h() const { cout << "h" << endl; };
};

struct TestBaseWrap : TestBase, bp::wrapper<TestBase> {
    void f() const {
        this->get_override("f")();
    }
};

void TestExport() {


  bp::class_<TestBaseWrap, boost::noncopyable>("TestBase")
    .def("f", bp::pure_virtual(&TestBase::f)) 
  ;

  bp::def("callf", FunctionPointer( [] (const TestBase &b) { b.f(); } ) );
}

class PyNumProc : public NumProc, public bp::wrapper<NumProc> {
    public:
        PyNumProc (PDE & apde, const Flags & flags) : NumProc(pde,flags) {
           cout << "pynumproc constructor" << endl; 
        }

        virtual void Do(LocalHeap & lh) {
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

  TestExport();

  bp::def("RegisterNumProc", FunctionPointer( [] (PyNumProc & np, string label, int dim) { 
      python_np = &np;
      GetNumProcs().AddNumProc (label, dim, FunctionPointer( [] (PDE & pde, const Flags & flags) { return (NumProc*) python_np;}  ));
//       RegisterNumProc<PyNumProc>(name, dim); return;
      } ));
      
      

  bp::class_<PyNumProc, shared_ptr<PyNumProc>,boost::noncopyable> ("NumProc", bp::no_init)
    // .def("__str__", &ToString<NumProc>)
    .def(bp::init<PDE &, const Flags &>())
    .def("Do", bp::pure_virtual(&PyNumProc::Do))
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
    /*
    .def("__setitem__", FunctionPointer([](PDE & self, string name, shared_ptr<GridFunction> gf)
                                        {
                                          self.AddGridFunction (name, gf);
                                        }))
    */
    .def("Add", FunctionPointer([](PDE & self, shared_ptr<PyNumProc> np)
                                {
                                cout << "add numproc" << endl;
                                  self.AddNumProc ("pynumproc", &*np);
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
