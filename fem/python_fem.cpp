#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <fem.hpp>
using namespace ngfem;



BOOST_PYTHON_MODULE(Ngfem) {

  cout << "init py - ngfem" << endl;

  bp::enum_<ELEMENT_TYPE>("ELEMENT_TYPE")
    .value("POINT", ET_POINT)     .value("SEGM", ET_SEGM)
    .value("TRIG", ET_TRIG)       .value("QUAD", ET_QUAD)
    .value("TET", ET_TET)         .value("PRISM", ET_PRISM)
    .value("PYRAMID", ET_PYRAMID) .value("HEX", ET_HEX)
    .export_values()
    ;

  bp::class_<FiniteElement, boost::noncopyable>("FiniteElement", bp::no_init)
    .add_property("ndof", &FiniteElement::GetNDof)    
    .add_property("order", &FiniteElement::Order)    
    .add_property("type", &FiniteElement::ElementType)    
    // .add_property("classname", &FiniteElement::ClassName)   // crashes ???
    ;
  
  bp::def("H1FE", FunctionPointer
          ([](ELEMENT_TYPE et, int order)
           {
             FiniteElement * fe = NULL;
             switch (et)
               {
               case ET_TRIG: fe = new H1HighOrderFE<ET_TRIG>(order); break;
               case ET_QUAD: fe = new H1HighOrderFE<ET_QUAD>(order); break;
               case ET_TET: fe = new H1HighOrderFE<ET_TET>(order); break;
               default: cerr << "cannot make fe " << et << endl;
               }
             return fe;
           }),
          bp::return_value_policy<bp::manage_new_object>()
          );

  bp::def("L2FE", FunctionPointer
          ([](ELEMENT_TYPE et, int order)
           {
             FiniteElement * fe = NULL;
             switch (et)
               {
               case ET_TRIG: fe = new L2HighOrderFE<ET_TRIG>(order); break;
               case ET_QUAD: fe = new L2HighOrderFE<ET_QUAD>(order); break;
               case ET_TET: fe = new L2HighOrderFE<ET_TET>(order); break;
               default: cerr << "cannot make fe " << et << endl;
               }
             return fe;
           }),
          bp::return_value_policy<bp::manage_new_object>()
          );
    
  bp::class_<ElementTransformation, boost::noncopyable>("ElementTransformation", bp::no_init)
    .def ("IsBoundary", &ElementTransformation::Boundary)
    .add_property("spacedim", &ElementTransformation::SpaceDim)
    ;
  
  bp::class_<BilinearFormIntegrator, shared_ptr<BilinearFormIntegrator>, boost::noncopyable>
    ("BilinearFormIntegrator", bp::no_init)
    .def("CalcElementMatrix", 
         static_cast<void(BilinearFormIntegrator::*) (const FiniteElement&, 
                                                      const ElementTransformation&,
                                                      FlatMatrix<double>&,LocalHeap&)const>
         (&BilinearFormIntegrator::CalcElementMatrix))
    ;
  bp::class_<LinearFormIntegrator, shared_ptr<LinearFormIntegrator>, boost::noncopyable>
    ("LinearFormIntegrator", bp::no_init)
    .def("CalcElementVector", 
         static_cast<void(LinearFormIntegrator::*)(const FiniteElement&, const ElementTransformation&, FlatVector<double>&,LocalHeap&)const>
         (&LinearFormIntegrator::CalcElementVector))
    ;

  bp::def("CreateBFI", FunctionPointer 
          ([](string name, int dim, CoefficientFunction & coef)
           {
             BilinearFormIntegrator * itor =
               GetIntegrators().CreateBFI (name, dim, &coef);
             
             if (!itor) cerr << "undefined integrator '" << name 
                             << "' in " << dim << " dimension having 1 coefficient"
                             << endl;
             
             return shared_ptr<BilinearFormIntegrator> (itor);
           }),
          (bp::arg("name")=NULL,bp::arg("dim")=2,bp::arg("coef")))
    ;

  bp::def("CreateLFI", FunctionPointer
          ([](string name, int dim, CoefficientFunction & coef)
           {
             LinearFormIntegrator * itor =
               GetIntegrators().CreateLFI (name, dim, &coef);
             
             if (!itor) cerr << "undefined integrator '" << name 
                             << "' in " << dim << " dimension having 1 coefficient"
                             << endl;
             
             return shared_ptr<LinearFormIntegrator> (itor);
           }),
          (bp::arg("name")=NULL,bp::arg("dim")=2,bp::arg("coef")))
    ;
  
  bp::class_<CoefficientFunction, boost::noncopyable>("CoefficientFunction", bp::no_init)
    ;

  bp::class_<ConstantCoefficientFunction,bp::bases<CoefficientFunction>>
    ("ConstantCF", bp::init<double>())
    ;



  bp::def("testlist", FunctionPointer
          ([](bp::object const & x) { cout << "any object"; }))
    ;

  bp::def("testlist", FunctionPointer
          ([](bp::list const & x) 
           { 
             cout << "got list of length "  << bp::len(x) << endl;
             for (int i = 0; i < bp::len(x); i++)
               {
                 cout << "list[" << i << "] = ";
                 CoefficientFunction & cf = bp::extract<CoefficientFunction&> (x[i]);
                 cf.PrintReport(cout);
               }
           }));
}

BOOST_PYTHON_MODULE(libngfem)
{
  cout << "execute 'import Ngstd' to import py-ngstd module" << endl;
  cout << "execute 'import Ngbla' to import py-ngbla module" << endl;
  cout << "execute 'import Ngfem' to import py-ngfem module" << endl;
}

struct Init {
  Init() 
  { 
    PyImport_AppendInittab("Ngfem", PyInit_Ngfem); 
  }
};
static Init init;



#endif
