#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <fem.hpp>
using namespace ngfem;
using ngfem::ELEMENT_TYPE;

void ExportNgfem() {
    std::string nested_name = "ngfem";
    if( bp::scope() )
      nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".ngfem");
    
    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << "exporting ngfem as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("ngfem") = module ;

    bp::scope ngbla_scope(module);


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
          ([](string name, int dim, shared_ptr<CoefficientFunction> coef)
           {
             // shared_ptr<CoefficientFunction> coef(&rcoef); // , NOOP_Deleter);
             auto bfi = GetIntegrators().CreateBFI (name, dim, coef);
             if (!bfi) cerr << "undefined integrator '" << name 
                            << "' in " << dim << " dimension having 1 coefficient"
                            << endl;
             
             return bfi;
           }),
          (bp::arg("name")=NULL,bp::arg("dim")=2,bp::arg("coef")))
    ;

  bp::def("CreateLFI", FunctionPointer
          ([](string name, int dim, shared_ptr<CoefficientFunction> coef)
           {
             auto lfi = GetIntegrators().CreateLFI (name, dim, coef);
             
             if (!lfi) cerr << "undefined integrator '" << name 
                             << "' in " << dim << " dimension having 1 coefficient"
                             << endl;
             
             return lfi;
           }),
          (bp::arg("name")=NULL,bp::arg("dim")=2,bp::arg("coef")))
    ;
  
  bp::class_<CoefficientFunction, shared_ptr<CoefficientFunction>, boost::noncopyable> 
    ("CoefficientFunction", bp::no_init)
    ;

  bp::class_<ConstantCoefficientFunction,bp::bases<CoefficientFunction>,
    shared_ptr<ConstantCoefficientFunction>, boost::noncopyable>
    ("ConstantCF", bp::init<double>())
    ;
  
  bp::implicitly_convertible 
    <shared_ptr<ConstantCoefficientFunction>, 
    shared_ptr<CoefficientFunction> >(); 

  
  // we better get rid of the template argument !
  bp::class_<DomainVariableCoefficientFunction<2>,bp::bases<CoefficientFunction>, 
    shared_ptr<DomainVariableCoefficientFunction<2>>, boost::noncopyable>
    ("VariableCF", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](string str)
                           {
                             // mem-leak -> smart pointer !!!
                             EvalFunction * ef = new EvalFunction (str);
                             return shared_ptr<DomainVariableCoefficientFunction<2>>
                               (new DomainVariableCoefficientFunction<2> (*ef));
                           })))
    ;

  bp::implicitly_convertible
    <shared_ptr<DomainVariableCoefficientFunction<2>>, 
    shared_ptr<CoefficientFunction> >(); 


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


void ExportNgstd();
void ExportNgbla();

BOOST_PYTHON_MODULE(libngfem) {
  ExportNgstd();
  ExportNgbla();
  ExportNgfem();
}



#endif
