#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <fem.hpp>
#include <mutex>
using namespace ngfem;
using ngfem::ELEMENT_TYPE;

struct PythonCoefficientFunction : public CoefficientFunction {
    PythonCoefficientFunction() { }

    virtual double EvaluateXYZ (double x, double y, double z) const = 0;

    bp::list GetCoordinates(const BaseMappedIntegrationPoint &bip ) {
        double x[3]{0};
        int dim = bip.GetTransformation().SpaceDim();
        const DimMappedIntegrationPoint<1,double> *ip1;
        const DimMappedIntegrationPoint<2,double> *ip2;
        const DimMappedIntegrationPoint<3,double> *ip3;
        switch(dim) {

            case 1:
                ip1 = static_cast<const DimMappedIntegrationPoint<1,double>*>(&bip);
                x[0] = ip1->GetPoint()[0];
                break;
            case 2:
                ip2 = static_cast<const DimMappedIntegrationPoint<2,double>*>(&bip);
                x[0] = ip2->GetPoint()[0];
                x[1] = ip2->GetPoint()[1];
                break;
            case 3:
                ip3 = static_cast<const DimMappedIntegrationPoint<3,double>*>(&bip);
                x[0] = ip3->GetPoint()[0];
                x[1] = ip3->GetPoint()[1];
                x[2] = ip3->GetPoint()[2];
                break;
            default:
                break;
        }
        bp::list list;
        int i;
        for(i=0; i<dim; i++)
            list.append(x[i]);
        for(i=0; i<3; i++)
            list.append(0.0);
        return list;
    }
};

class PythonCFWrap : public PythonCoefficientFunction , public bp::wrapper<PythonCoefficientFunction> {
    static std::mutex m;
    public:
        PythonCFWrap () : PythonCoefficientFunction() { ; }
        double EvaluateXYZ (double x, double y, double z) const {
            return this->get_override("EvaluateXYZ")(x,y,z);
        }

        double Evaluate (const BaseMappedIntegrationPoint & bip) const {
            double ret = 0;
            m.lock();
            try { 
                ret = this->get_override("Evaluate")(boost::ref(bip)); 
            }
            catch (bp::error_already_set const &) {
                PyErr_Print();
            }
            catch(...) {
                cout << "caught Exception in PythonCoefficientFunction::Evaluate" << endl;
            }
            m.unlock();
            return ret;
        }
};

std::mutex PythonCFWrap::m;




template <typename OP, typename OPC> 
class cl_UnaryOpCF : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  OP lam;
  OPC lamc;
public:
  cl_UnaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                OP alam, OPC alamc)
    : c1(ac1), lam(alam), lamc(alamc) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex(); }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return lam (c1->Evaluate(ip));
  }

  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return lamc (c1->EvaluateComplex(ip));
  }

  virtual double EvaluateConst () const
  {
    return lam (c1->EvaluateConst());
  }
};

template <typename OP, typename OPC> 
shared_ptr<CoefficientFunction> UnaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                          OP lam, OPC lamc)
{
  return shared_ptr<CoefficientFunction> (new cl_UnaryOpCF<OP,OPC> (c1, lam, lamc));
}


template <typename OP> 
class cl_BinaryOpCF : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1, c2;
  OP lam;
public:
  cl_BinaryOpCF (shared_ptr<CoefficientFunction> ac1, 
                 shared_ptr<CoefficientFunction> ac2, 
                 OP alam)
    : c1(ac1), c2(ac2), lam(alam) { ; }
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return lam (c1->Evaluate(ip), c2->Evaluate(ip));
  }
  virtual double EvaluateConst () const
  {
    return lam (c1->EvaluateConst(), c2->EvaluateConst());
  }
  virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                        FlatMatrix<> result) const
  {
#ifdef VLA
    double hmem[ir.Size()];
    FlatMatrix<> temp(ir.Size(), 1, hmem);
#else
    Matrix<> temp(ir.Size(), 1);
#endif

    c1->Evaluate (ir, result);
    c2->Evaluate (ir, temp);
    for (int i = 0; i < ir.Size(); i++)
      result(i,0) = lam (result(i,0), temp(i,0));
  }
};

template <typename OP> 
INLINE shared_ptr<CoefficientFunction> BinaryOpCF(shared_ptr<CoefficientFunction> c1, 
                                                  shared_ptr<CoefficientFunction> c2, 
                                                  OP lam)
{
  return shared_ptr<CoefficientFunction> (new cl_BinaryOpCF<OP> (c1, c2, lam));
}


class ScaleCoefficientFunction : public CoefficientFunction
{
  double scal;
  shared_ptr<CoefficientFunction> c1;
public:
  ScaleCoefficientFunction (double ascal, 
                            shared_ptr<CoefficientFunction> ac1)
    : scal(ascal), c1(ac1) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual int Dimension() const { return c1->Dimension(); }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->Evaluate(ip);
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->EvaluateComplex(ip);
  }
  virtual double EvaluateConst () const
  {
    return scal * c1->EvaluateConst();
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
};


class ScaleCoefficientFunctionC : public CoefficientFunction
{
  Complex scal;
  shared_ptr<CoefficientFunction> c1;
public:
  ScaleCoefficientFunctionC (Complex ascal, 
                            shared_ptr<CoefficientFunction> ac1)
    : scal(ascal), c1(ac1) { ; }
  
  virtual bool IsComplex() const { return true; }
  virtual int Dimension() const { return c1->Dimension(); }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    throw Exception ("real Evaluate called for complex CF");
  }
  virtual Complex EvaluateComplex (const BaseMappedIntegrationPoint & ip) const 
  {
    return scal * c1->EvaluateComplex(ip);    
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<Complex> result) const
  {
    c1->Evaluate (ip, result);
    result *= scal;
  }
    
};


class MultScalVecCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;  // scalar
  shared_ptr<CoefficientFunction> c2;  // vector
public:
  MultScalVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                  shared_ptr<CoefficientFunction> ac2)
    : c1(ac1), c2(ac2) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  virtual int Dimension() const { return c2->Dimension(); }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    throw Exception ("double MultScalVecCF::Evaluate called");
  }

  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    Vec<1> v1;
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, result);
    result *= v1(0);
  }
};


class MultVecVecCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  shared_ptr<CoefficientFunction> c2;
public:
  MultVecVecCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                 shared_ptr<CoefficientFunction> ac2)
    : c1(ac1), c2(ac2) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex() || c2->IsComplex(); }
  virtual int Dimension() const { return 1; }
  
  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    Vec<1> res;
    Evaluate (ip, res);
    return res(0);
  }
  virtual void Evaluate(const BaseMappedIntegrationPoint & ip,
                        FlatVector<> result) const
  {
    Vector<> v1(c1->Dimension()), v2(c2->Dimension());
    c1->Evaluate (ip, v1);
    c2->Evaluate (ip, v2);
    result(0) = InnerProduct (v1, v2);
  }
};


class ComponentCoefficientFunction : public CoefficientFunction
{
  shared_ptr<CoefficientFunction> c1;
  int comp;
public:
  ComponentCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                int acomp)
    : c1(ac1), comp(acomp) { ; }
  
  virtual bool IsComplex() const { return c1->IsComplex(); }
  virtual int Dimension() const { return 1; }

  virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
  {
    Vector<> v1(c1->Dimension());
    c1->Evaluate (ip, v1);
    return v1(comp);
  }
};





shared_ptr<CoefficientFunction> MakeCoefficient (bp::object py_coef)
{
  if (bp::extract<shared_ptr<CoefficientFunction>>(py_coef).check())
    return bp::extract<shared_ptr<CoefficientFunction>>(py_coef)();
  else if (bp::extract<double>(py_coef).check())
    return make_shared<ConstantCoefficientFunction> 
      (bp::extract<double>(py_coef)());
  else
    {
      bp::exec("raise KeyError()\n");
      return nullptr;
    }
}

Array<shared_ptr<CoefficientFunction>> MakeCoefficients (bp::object py_coef)
{
  Array<shared_ptr<CoefficientFunction>> tmp;
  if (bp::extract<bp::list>(py_coef).check())
    {
      auto l = bp::extract<bp::list>(py_coef)();
      for (int i = 0; i < bp::len(l); i++)
        tmp += MakeCoefficient(l[i]);
    }
  else if (bp::extract<bp::tuple>(py_coef).check())
    {
      auto l = bp::extract<bp::tuple>(py_coef)();
      for (int i = 0; i < bp::len(l); i++)
        tmp += MakeCoefficient(l[i]);
    }
  else
    tmp += MakeCoefficient(py_coef);
  return move(tmp);
}


void ExportCoefficientFunction()
{
   typedef CoefficientFunction CF;
  typedef shared_ptr<CF> SPCF;

  bp::class_<CoefficientFunction, shared_ptr<CoefficientFunction>, boost::noncopyable> 
    ("CoefficientFunction", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](double val) -> shared_ptr<CoefficientFunction>
                           {
                             return make_shared<ConstantCoefficientFunction> (val);
                           })))
    .def("Evaluate", static_cast<double (CoefficientFunction::*)(const BaseMappedIntegrationPoint &) const>(&CoefficientFunction::Evaluate))

    .add_property("dim", &CoefficientFunction::Dimension)    
    
    .def("__getitem__", FunctionPointer( [](SPCF & self, int comp) -> SPCF
                                         {
                                           return make_shared<ComponentCoefficientFunction> (self, comp); 
                                         }))

    // coefficient expressions
    .def ("__add__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           { return BinaryOpCF (c1, c2, [](double a, double b) { return a+b; });} ))
    .def ("__add__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a+b; }); }))
    .def ("__radd__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a+b; }); }))

    .def ("__sub__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           { return BinaryOpCF (c1, c2, [](double a, double b) { return a-b; });} ))
    .def ("__sub__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a-b; }); }))
    .def ("__rsub__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return b-a; }); }))

    .def ("__mul__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           { 
             if (c1->Dimension() > 1 && c2->Dimension() > 1)
               return make_shared<MultVecVecCoefficientFunction> (c1, c2);
             if (c1->Dimension() == 1 && c2->Dimension() > 1)
               return make_shared<MultScalVecCoefficientFunction> (c1, c2);
             if (c1->Dimension() > 1 && c2->Dimension() == 1)
               return make_shared<MultScalVecCoefficientFunction> (c2, c1);
             return BinaryOpCF (c1, c2, [](double a, double b) { return a*b; });
           } ))
    .def ("__mul__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return make_shared<ScaleCoefficientFunction> (val, coef); }))
    .def ("__rmul__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return make_shared<ScaleCoefficientFunction> (val, coef); }))
    .def ("__mul__", FunctionPointer 
          ([] (SPCF coef, Complex val) -> SPCF
           { return make_shared<ScaleCoefficientFunctionC> (val, coef); }))
    .def ("__rmul__", FunctionPointer 
          ([] (SPCF coef, Complex val) -> SPCF
           { 
             if (val.imag() == 0)
               return make_shared<ScaleCoefficientFunction> (val.real(), coef); 
             else
               return make_shared<ScaleCoefficientFunctionC> (val, coef); 
           }))

    // { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
    // [](double a, double b) { return a*b; }); }))

    .def ("__neg__", FunctionPointer 
          ([] (SPCF coef) -> SPCF
           { return make_shared<ScaleCoefficientFunction> (-1, coef); }))
    ;

  bp::def ("sin", FunctionPointer 
           ([] (SPCF coef) -> SPCF
            { return UnaryOpCF (coef, 
                                [](double a) { return sin(a); },
                                [](Complex a) { return sin(a); }); }));
  bp::def ("sin", FunctionPointer ([] (double d) -> double { return sin (d); }));
  bp::def ("sin", FunctionPointer ([] (Complex d) -> Complex { return sin (d); }));

  bp::def ("exp", FunctionPointer 
           ([] (SPCF coef) -> SPCF
            { return UnaryOpCF (coef, 
                                [](double a) { return exp(a); },
                                [](Complex a) { return exp(a); }); }));
  bp::def ("exp", FunctionPointer ([] (double d) -> double { return exp (d); }));

  
  bp::class_<ConstantCoefficientFunction,bp::bases<CoefficientFunction>,
    shared_ptr<ConstantCoefficientFunction>, boost::noncopyable>
    ("ConstantCF", bp::init<double>())
    ;
  bp::implicitly_convertible 
    <shared_ptr<ConstantCoefficientFunction>, shared_ptr<CoefficientFunction> >(); 

  class CoordCoefficientFunction : public CoefficientFunction
  {
    int dir;
  public:
    CoordCoefficientFunction (int adir) : dir(adir) { ; }
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
    {
      return ip.GetPoint()(dir);
    }
    virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                          FlatMatrix<> result) const
    {
      for (int i = 0; i < ir.Size(); i++)
        result(i,0) = ir[i].GetPoint()(dir);
    }
  };

  bp::class_<CoordCoefficientFunction,bp::bases<CoefficientFunction>,
    shared_ptr<CoordCoefficientFunction>, boost::noncopyable>
    ("CoordCF", bp::init<int>())
    ;
  bp::implicitly_convertible 
    <shared_ptr<CoordCoefficientFunction>, shared_ptr<CoefficientFunction> >(); 
}



void NGS_DLL_HEADER ExportNgfem() {
    std::string nested_name = "fem";
    if( bp::scope() )
      nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".fem");

    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << "exporting fem as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("fem") = module ;

    bp::scope ngbla_scope(module);


  bp::enum_<ELEMENT_TYPE>("ET")
    .value("POINT", ET_POINT)     .value("SEGM", ET_SEGM)
    .value("TRIG", ET_TRIG)       .value("QUAD", ET_QUAD)
    .value("TET", ET_TET)         .value("PRISM", ET_PRISM)
    .value("PYRAMID", ET_PYRAMID) .value("HEX", ET_HEX)
    .export_values()
    ;


  bp::class_<FiniteElement, shared_ptr<FiniteElement>, boost::noncopyable>("FiniteElement", bp::no_init)
    .add_property("ndof", &FiniteElement::GetNDof)    
    .add_property("order", &FiniteElement::Order)    
    .add_property("type", &FiniteElement::ElementType)    
    .add_property("dim", &FiniteElement::Dim)    
    .add_property("classname", &FiniteElement::ClassName)  
    .def("__str__", &ToString<FiniteElement>)
    ;

  bp::class_<BaseScalarFiniteElement, shared_ptr<BaseScalarFiniteElement>, 
    bp::bases<FiniteElement>, boost::noncopyable>("ScalarFE", bp::no_init)
    .def("CalcShape",
         FunctionPointer
         ([] (const BaseScalarFiniteElement & fe, double x, double y, double z)
          {
            IntegrationPoint ip(x,y,z);
            Vector<> v(fe.GetNDof());
            switch (fe.Dim())
              {
              case 1:
                dynamic_cast<const ScalarFiniteElement<1>&> (fe).CalcShape(ip, v); break;
              case 2:
                dynamic_cast<const ScalarFiniteElement<2>&> (fe).CalcShape(ip, v); break;
              case 3:
                dynamic_cast<const ScalarFiniteElement<3>&> (fe).CalcShape(ip, v); break;
              default:
                ;
              }
            return v;
          }),
         (bp::arg("self"),bp::arg("x"),bp::arg("y")=0.0,bp::arg("z")=0.0)
         )
    ;

  bp::implicitly_convertible 
    <shared_ptr<BaseScalarFiniteElement>, 
    shared_ptr<FiniteElement> >(); 


  bp::def("H1FE", FunctionPointer
          ([](ELEMENT_TYPE et, int order)
           {
             BaseScalarFiniteElement * fe = nullptr;
             switch (et)
               {
               case ET_TRIG: fe = new H1HighOrderFE<ET_TRIG>(order); break;
               case ET_QUAD: fe = new H1HighOrderFE<ET_QUAD>(order); break;
               case ET_TET: fe = new H1HighOrderFE<ET_TET>(order); break;
               default: cerr << "cannot make fe " << et << endl;
               }
             return shared_ptr<BaseScalarFiniteElement>(fe);
           })
          );

  bp::def("L2FE", FunctionPointer
          ([](ELEMENT_TYPE et, int order)
           {
             BaseScalarFiniteElement * fe = nullptr;
             switch (et)
               {
               case ET_TRIG: fe = new L2HighOrderFE<ET_TRIG>(order); break;
               case ET_QUAD: fe = new L2HighOrderFE<ET_QUAD>(order); break;
               case ET_TET: fe = new L2HighOrderFE<ET_TET>(order); break;
               default: cerr << "cannot make fe " << et << endl;
               }
             return shared_ptr<BaseScalarFiniteElement>(fe);
           })
          );
    
  bp::class_<ElementTransformation, boost::noncopyable>("ElementTransformation", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](ELEMENT_TYPE et, Matrix<> pmat) -> ElementTransformation*
                           {
                             return new FE_ElementTransformation<2,2> (et, pmat); 
                           })))
    .def ("IsBoundary", &ElementTransformation::Boundary)
    .add_property("spacedim", &ElementTransformation::SpaceDim)
    ;
  
  bp::class_<BilinearFormIntegrator, shared_ptr<BilinearFormIntegrator>, boost::noncopyable>
    ("BFI", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](string name, int dim, bp::object py_coef, bool imag)
                           {
                             Array<shared_ptr<CoefficientFunction>> coef = MakeCoefficients(py_coef);
                             auto bfi = GetIntegrators().CreateBFI (name, dim, coef);

                             if (!bfi) cerr << "undefined integrator '" << name 
                                            << "' in " << dim << " dimension" << endl;

                             if (imag)
                               bfi = make_shared<ComplexBilinearFormIntegrator> (bfi, Complex(0,1));

                             return bfi;
                           }),
          bp::default_call_policies(),        // need it to use named arguments
          (bp::arg("name")=NULL,bp::arg("dim")=-1,bp::arg("coef"),bp::arg("imag")=false)))
    
    .def("__str__", &ToString<BilinearFormIntegrator>)

    .def("CalcElementMatrix", 
         static_cast<void(BilinearFormIntegrator::*) (const FiniteElement&, 
                                                      const ElementTransformation&,
                                                      FlatMatrix<double>,LocalHeap&)const>
         (&BilinearFormIntegrator::CalcElementMatrix))

    .def("CalcElementMatrix",
         FunctionPointer([] (const BilinearFormIntegrator & self, 
                             const FiniteElement & fe, const ElementTransformation & trafo,
                             LocalHeap & lh)
                         {
                           Matrix<> mat(fe.GetNDof());
                           self.CalcElementMatrix (fe, trafo, mat, lh);
                           return mat;
                         }),
         bp::default_call_policies(),        // need it to use named arguments
         (bp::arg("bfi")=NULL, bp::arg("fel"),bp::arg("trafo"),bp::arg("localheap")))

    .def("CalcElementMatrix",
         FunctionPointer([] (const BilinearFormIntegrator & self, 
                             const FiniteElement & fe, const ElementTransformation & trafo,
                             int heapsize)
                         {
                           Matrix<> mat(fe.GetNDof());
                           while (true)
                             {
                               try
                                 {
                                   LocalHeap lh(heapsize);
                                   self.CalcElementMatrix (fe, trafo, mat, lh);
                                   return mat;
                                 }
                               catch (LocalHeapOverflow ex)
                                 {
                                   heapsize *= 10;
                                 }
                             };
                         }),
         bp::default_call_policies(),        // need it to use named arguments
         (bp::arg("bfi")=NULL, bp::arg("fel"),bp::arg("trafo"),bp::arg("heapsize")=10000))
    ;


  bp::def("CompoundBFI", 
          (FunctionPointer ([]( shared_ptr<BilinearFormIntegrator> bfi, int comp )
                            {
                                return make_shared<CompoundBilinearFormIntegrator>(bfi, comp);
                            })),
           bp::default_call_policies(),     // need it to use named arguments
           (bp::arg("bfi")=NULL, bp::arg("comp")=0)
      );

  bp::def("BlockBFI", 
          (FunctionPointer ([]( shared_ptr<BilinearFormIntegrator> bfi, int dim, int comp )
                            {
                                return make_shared<BlockBilinearFormIntegrator>(bfi, dim, comp);
                            })),
           bp::default_call_policies(),     // need it to use named arguments
           (bp::arg("bfi")=NULL, bp::arg("dim")=2, bp::arg("comp")=0)
      );


  bp::class_<CompoundBilinearFormIntegrator,bp::bases<BilinearFormIntegrator>,
             shared_ptr<CompoundBilinearFormIntegrator>, boost::noncopyable>
      ("CompoundBilinearFormIntegrator", bp::no_init);

  bp::class_<BlockBilinearFormIntegrator,bp::bases<BilinearFormIntegrator>,
             shared_ptr<BlockBilinearFormIntegrator>, boost::noncopyable>
      ("BlockBilinearFormIntegrator", bp::no_init);

  bp::implicitly_convertible 
    <shared_ptr<CompoundBilinearFormIntegrator>, 
    shared_ptr<BilinearFormIntegrator> >(); 

  bp::implicitly_convertible 
    <shared_ptr<BlockBilinearFormIntegrator>, 
    shared_ptr<BilinearFormIntegrator> >(); 


  bp::class_<LinearFormIntegrator, shared_ptr<LinearFormIntegrator>, boost::noncopyable>
    ("LFI", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](string name, int dim, 
                              bp::object py_coef,
                              bp::object definedon, bool imag, const Flags & flags)
                           {
                             Array<shared_ptr<CoefficientFunction>> coef = MakeCoefficients(py_coef);
                             auto lfi = GetIntegrators().CreateLFI (name, dim, coef);
                             
                             if (bp::extract<bp::list> (definedon).check())
                               lfi -> SetDefinedOn (makeCArray<int> (definedon));
 
                             if (imag)
                               lfi = make_shared<ComplexLinearFormIntegrator> (lfi, Complex(0,1));

                             // cout << "LFI: Flags = " << flags << endl;
                             if (!lfi) cerr << "undefined integrator '" << name 
                                            << "' in " << dim << " dimension having 1 coefficient"
                                            << endl;
                             return lfi;
                           }),
          bp::default_call_policies(),     // need it to use named arguments
          (bp::arg("name")=NULL,bp::arg("dim")=-1,
           bp::arg("coef"),bp::arg("definedon")=bp::object(), 
           bp::arg("imag")=false, bp::arg("flags")=bp::dict()))
         )

    .def("__init__", bp::make_constructor
         (FunctionPointer ([](string name, int dim, bp::list coefs_list,
                              bp::object definedon, bool imag, const Flags & flags)
                           {
                             Array<shared_ptr<CoefficientFunction> > coefs = makeCArray<shared_ptr<CoefficientFunction>> (coefs_list);
                             auto lfi = GetIntegrators().CreateLFI (name, dim, coefs);
                             
                             if (bp::extract<bp::list> (definedon).check())
                               lfi -> SetDefinedOn (makeCArray<int> (definedon));
 
                             if (imag)
                               lfi = make_shared<ComplexLinearFormIntegrator> (lfi, Complex(0,1));

                             // cout << "LFI: Flags = " << flags << endl;
                             if (!lfi) cerr << "undefined integrator '" << name 
                                            << "' in " << dim << " dimension having 1 coefficient"
                                            << endl;
                             return lfi;
                           }),
          bp::default_call_policies(),     // need it to use named arguments
          (bp::arg("name")=NULL,bp::arg("dim")=-1,
           bp::arg("coef"),bp::arg("definedon")=bp::object(), 
           bp::arg("imag")=false, bp::arg("flags")=bp::dict()))
        )

    .def("CalcElementVector", 
         static_cast<void(LinearFormIntegrator::*)(const FiniteElement&, const ElementTransformation&, FlatVector<double>,LocalHeap&)const>
         (&LinearFormIntegrator::CalcElementVector))
    ;



  bp::def("CompoundLFI", 
          (FunctionPointer ([]( shared_ptr<LinearFormIntegrator> lfi, int comp )
                            {
                                return make_shared<CompoundLinearFormIntegrator>(lfi, comp);
                            })),
           bp::default_call_policies(),     // need it to use named arguments
           (bp::arg("lfi")=NULL, bp::arg("comp")=0)
      );

  bp::def("BlockLFI", 
          (FunctionPointer ([]( shared_ptr<LinearFormIntegrator> lfi, int dim, int comp )
                            {
                                return make_shared<BlockLinearFormIntegrator>(lfi, dim, comp);
                            })),
           bp::default_call_policies(),     // need it to use named arguments
           (bp::arg("lfi")=NULL, bp::arg("dim")=2, bp::arg("comp")=0)
      );


  bp::class_<CompoundLinearFormIntegrator,bp::bases<LinearFormIntegrator>,
             shared_ptr<CompoundLinearFormIntegrator>, boost::noncopyable>
      ("CompoundLinearFormIntegrator", bp::no_init);

  bp::class_<BlockLinearFormIntegrator,bp::bases<LinearFormIntegrator>,
             shared_ptr<BlockLinearFormIntegrator>, boost::noncopyable>
      ("BlockLinearFormIntegrator", bp::no_init);

  bp::implicitly_convertible 
    <shared_ptr<CompoundLinearFormIntegrator>, 
    shared_ptr<LinearFormIntegrator> >(); 

  bp::implicitly_convertible 
    <shared_ptr<BlockLinearFormIntegrator>, 
    shared_ptr<LinearFormIntegrator> >(); 


  /*
  bp::def("CreateBFI", FunctionPointer 
          ([](string name, int dim, shared_ptr<CoefficientFunction> coef)
           {
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
  */
  
  bp::class_<BaseMappedIntegrationPoint, boost::noncopyable>( "BaseMappedIntegrationPoint", bp::no_init);  




  ExportCoefficientFunction ();
  /*
  typedef CoefficientFunction CF;
  typedef shared_ptr<CF> SPCF;

  bp::class_<CoefficientFunction, shared_ptr<CoefficientFunction>, boost::noncopyable> 
    ("CoefficientFunction", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](double val) -> shared_ptr<CoefficientFunction>
                           {
                             return make_shared<ConstantCoefficientFunction> (val);
                           })))
    .def("Evaluate", static_cast<double (CoefficientFunction::*)(const BaseMappedIntegrationPoint &) const>(&CoefficientFunction::Evaluate))



    // coefficient expressions

    .def ("__add__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           { return BinaryOpCF (c1, c2, [](double a, double b) { return a+b; });} ))
    .def ("__add__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a+b; }); }))
    .def ("__radd__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a+b; }); }))

    .def ("__sub__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           { return BinaryOpCF (c1, c2, [](double a, double b) { return a-b; });} ))
    .def ("__sub__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a-b; }); }))
    .def ("__rsub__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return b-a; }); }))

    .def ("__mul__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           { return BinaryOpCF (c1, c2, [](double a, double b) { return a*b; });} ))
    .def ("__mul__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a*b; }); }))
    .def ("__rmul__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a*b; }); }))

    .def ("__neg__", FunctionPointer 
          ([] (SPCF c1) -> SPCF
           { return UnaryOpCF (c1, [](double a) { return -a; });} ))
    ;
  bp::def ("sin", FunctionPointer 
           ([] (SPCF coef) -> SPCF
            { return UnaryOpCF (coef, [](double a) { return sin(a); }); }));
  bp::def ("sin", FunctionPointer ([] (double d) -> double { return sin (d); }));

  bp::def ("exp", FunctionPointer 
           ([] (SPCF coef) -> SPCF
            { return UnaryOpCF (coef, [](double a) { return exp(a); }); }));
  bp::def ("exp", FunctionPointer ([] (double d) -> double { return exp (d); }));

  
  bp::class_<ConstantCoefficientFunction,bp::bases<CoefficientFunction>,
    shared_ptr<ConstantCoefficientFunction>, boost::noncopyable>
    ("ConstantCF", bp::init<double>())
    ;
  bp::implicitly_convertible 
    <shared_ptr<ConstantCoefficientFunction>, shared_ptr<CoefficientFunction> >(); 


  class CoordCoefficientFunction : public CoefficientFunction
  {
    int dir;
  public:
    CoordCoefficientFunction (int adir) : dir(adir) { ; }
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
    {
      return ip.GetPoint()(dir);
    }
  };

  bp::class_<CoordCoefficientFunction,bp::bases<CoefficientFunction>,
    shared_ptr<CoordCoefficientFunction>, boost::noncopyable>
    ("CoordCF", bp::init<int>())
    ;
  bp::implicitly_convertible 
    <shared_ptr<CoordCoefficientFunction>, shared_ptr<CoefficientFunction> >(); 


  */



  bp::class_<PythonCFWrap ,bp::bases<CoefficientFunction>, shared_ptr<PythonCFWrap>, boost::noncopyable>("PythonCF", bp::init<>())
//     .def("MyEvaluate", bp::pure_virtual(static_cast<double (PythonCoefficientFunction::*)(double,double) const>(&PythonCoefficientFunction::MyEvaluate))) 
    .def("Evaluate", bp::pure_virtual(static_cast<double (CoefficientFunction::*)(const BaseMappedIntegrationPoint &) const>(&CoefficientFunction::Evaluate))) 
    .def("GetCoordinates", &PythonCoefficientFunction::GetCoordinates)
//     .def("MyEvaluate", bp::pure_virtual(&PythonCoefficientFunction::MyEvaluate)) 
    ;

  bp::implicitly_convertible 
    <shared_ptr<PythonCFWrap>, 
    shared_ptr<CoefficientFunction> >(); 

  
  bp::class_<DomainVariableCoefficientFunction,bp::bases<CoefficientFunction>, 
    shared_ptr<DomainVariableCoefficientFunction>, boost::noncopyable>
    ("VariableCF", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](string str)
                           {
                             auto ef = make_shared<EvalFunction> (str);
                             return make_shared<DomainVariableCoefficientFunction> 
                               (Array<shared_ptr<EvalFunction>> ({ ef }));
                           })))
    ;

  bp::implicitly_convertible
    <shared_ptr<DomainVariableCoefficientFunction>, 
    shared_ptr<CoefficientFunction> >(); 


  bp::class_<DomainConstantCoefficientFunction,bp::bases<CoefficientFunction>, 
    shared_ptr<DomainConstantCoefficientFunction>, boost::noncopyable>
    ("DomainConstantCF", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](bp::object coefs)->shared_ptr<DomainConstantCoefficientFunction>
                           {
                             Array<double> darray (makeCArray<double> (coefs));
                             return make_shared<DomainConstantCoefficientFunction> (darray);
                           })))
    ;

  bp::implicitly_convertible
    <shared_ptr<DomainConstantCoefficientFunction>, 
    shared_ptr<CoefficientFunction> >(); 







}


void ExportNgstd();
void ExportNgbla();

BOOST_PYTHON_MODULE(libngfem) {
  // ExportNgstd();
  // ExportNgbla();
  ExportNgfem();
}



#endif
