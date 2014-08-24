#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <fem.hpp>
using namespace ngfem;

struct ElementTransformationWrap : ElementTransformation, bp::wrapper<ElementTransformation>
{
  void CalcJacobian (const IntegrationPoint & ip,
                     FlatMatrix<> dxdxi) const
  {
    ; // this -> get_override("CalcJacobian")(ip, dxdxi);
  }
  
};


class Base
{
public:
  Base () { cout << "constr Base" << endl; }
  virtual void f() { cout << "hi from base"; }
  virtual void g() = 0;
  virtual ~Base() { cout << "destr base"; }
  
};

class Derived : public Base
{
public:
  Derived () { cout << "constr Derived" << endl; }
  virtual void f() { cout << "hi from derived"; }
  virtual void g() { cout << "g is overloaded"; }
  virtual ~Derived() { cout << "destr derived"; }
};



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
  
  bp::class_<ElementTransformation, boost::noncopyable>("ElementTransformation", bp::no_init)
    .def ("IsBoundary", &ElementTransformation::Boundary)
    .add_property("spacedim", &ElementTransformation::SpaceDim)
    ;



  bp::class_<BilinearFormIntegrator, boost::noncopyable>("BilinearFormIntegrator", bp::no_init)
    .def("CalcElementMatrix", 
         static_cast<void(BilinearFormIntegrator::*) (const FiniteElement&, 
                                                      const ElementTransformation&,
                                                     FlatMatrix<double>&,LocalHeap&)const>
         (&BilinearFormIntegrator::CalcElementMatrix))
    ;

  bp::class_<CoefficientFunction, boost::noncopyable>("CoefficientFunction", bp::no_init)
    ;

  bp::class_<ConstantCoefficientFunction,bp::bases<CoefficientFunction>>
    ("ConstantCoefficientFunction", bp::init<double>())
    ;

  bp::class_<LaplaceIntegrator<2>, bp::bases<BilinearFormIntegrator>>
    ("LaplaceIntegrator_2d", bp::init<CoefficientFunction*>())  
    ;

  // testing
  bp::class_<Base, boost::noncopyable>("Base", bp::no_init)
    .def ("f", &Base::f)
    .def ("g", &Base::g)
    ;
  bp::class_<Derived, bp::bases<Base>>("Derived")
    ;


}


#endif
