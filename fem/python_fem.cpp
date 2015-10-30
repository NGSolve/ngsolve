#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include "../ngstd/bspline.hpp"
#include <fem.hpp>
#include <mutex>
using namespace ngfem;
using ngfem::ELEMENT_TYPE;

#include "pml.hpp"

namespace ngfem
{
  extern SymbolTable<double> * constant_table_for_FEM;
  SymbolTable<double> pmlpar;
}


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






/*
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
*/

shared_ptr<CoefficientFunction> MakeCoefficient (bp::object val)
{
  bp::extract<shared_ptr<CoefficientFunction>> ecf(val);
  if (ecf.check()) return ecf();

  bp::extract<double> ed(val);
  if (ed.check()) 
    return make_shared<ConstantCoefficientFunction> (ed());

  bp::extract<bp::list> el(val);
  if (el.check())
    {
      Array<shared_ptr<CoefficientFunction>> cflist(bp::len(el()));
      for (int i : Range(cflist))
        cflist[i] = MakeCoefficient(el()[i]);
      return make_shared<DomainWiseCoefficientFunction>(move(cflist));
    }

  bp::extract<bp::tuple> et(val);
  if (et.check())
    {
      Array<shared_ptr<CoefficientFunction>> cflist(bp::len(et()));
      for (int i : Range(cflist))
        cflist[i] = MakeCoefficient(et()[i]);
      return make_shared<VectorialCoefficientFunction>(move(cflist));
    }


  throw Exception ("cannot make coefficient");
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

typedef CoefficientFunction CF;
typedef shared_ptr<CF> SPCF;


template <typename FUNC>
void ExportStdMathFunction(string name)
{
  bp::def (name.c_str(), FunctionPointer 
           ([] (bp::object x) -> bp::object
            {
              FUNC func;
              if (bp::extract<SPCF>(x).check())
                {
                  auto coef = bp::extract<SPCF>(x)();
                  return bp::object(UnaryOpCF(coef, func, func));
                }
              bp::extract<double> ed(x);
              if (ed.check()) return bp::object(func(ed()));
              if (bp::extract<Complex> (x).check())
                return bp::object(func(bp::extract<Complex> (x)()));
              throw Exception ("can't compute sin");
            }));
}


using std::exp;
template <int D, typename SCAL>
AutoDiff<D,SCAL> exp (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = std::exp(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * res.Value();
  return res;
}



struct GenericSin {
  template <typename T> T operator() (T x) const { return sin(x); }
  AutoDiff<1> operator() (AutoDiff<1> x) const { throw Exception ("sin(AD) not implemented"); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception ("sin(ADD) not implemented"); }
};
struct GenericCos {
  template <typename T> T operator() (T x) const { return cos(x); }
  AutoDiff<1> operator() (AutoDiff<1> x) const { throw Exception ("cos(AD) not implemented"); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception ("cos(ADD) not implemented"); }
};
struct GenericTan {
  template <typename T> T operator() (T x) const { return tan(x); }
  AutoDiff<1> operator() (AutoDiff<1> x) const { throw Exception ("tan(AD) not implemented"); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception ("tan(ADD) not implemented"); }
};
struct GenericExp {
  template <typename T> T operator() (T x) const { return exp(x); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception ("exp(ADD) not implemented"); }
};
struct GenericLog {
  template <typename T> T operator() (T x) const { return log(x); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception ("log(ADD) not implemented"); }
};
struct GenericSqrt {
  template <typename T> T operator() (T x) const { return sqrt(x); }
};
struct GenericConj {
  template <typename T> T operator() (T x) const { return Conj(x); } // from bla
  AutoDiff<1> operator() (AutoDiff<1> x) const { throw Exception ("Conj(AD) not implemented"); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception ("Conj(AD) not implemented"); }
};


  template <int D>
  class NormalVectorCF : public CoefficientFunction
  {
  public:
    NormalVectorCF () { ; }
    virtual int Dimension() const { return D; }

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
    {
      return 0;
    }
    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> res) const 
    {
      if (ip.Dim() != D)
        throw Exception("illegal dim of normal vector");
      res = static_cast<const DimMappedIntegrationPoint<D>&>(ip).GetNV();
    }
  };



void ExportCoefficientFunction()
{
  bp::class_<CoefficientFunction, shared_ptr<CoefficientFunction>, boost::noncopyable> 
    ("CoefficientFunction", bp::no_init)

    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](bp::object val, bp::object dims) 
                           {
                             auto coef = MakeCoefficient(val);
                             if (dims)
                               {
                                 Array<int> cdims = makeCArray<int> (dims);
                                 dynamic_pointer_cast<VectorialCoefficientFunction> (coef)->SetDimensions(cdims);
                               }
                             return coef;
                           }),
          bp::default_call_policies(),        // need it to use named arguments
          (bp::arg("coef"),bp::arg("dims")=bp::object())
          ))
    
    .def("__call__", FunctionPointer
	 ([] (CoefficientFunction & self, BaseMappedIntegrationPoint & mip) -> bp::object
	  {
	    if (!self.IsComplex())
	      {
                if (self.Dimension() == 1)
                  return bp::object(self.Evaluate(mip));
                Vector<> vec(self.Dimension());
                self.Evaluate (mip, vec);
                return bp::tuple(vec);
	      }
	    else
	      {
		if (self.Dimension() == 1)
		  return bp::object(self.EvaluateComplex(mip));
                Vector<Complex> vec(self.Dimension());
                self.Evaluate (mip, vec);
                return bp::tuple(vec);
	      }
	  }))
    .add_property("dim", &CoefficientFunction::Dimension)

    .add_property("dims", &CoefficientFunction::Dimensions)    
    
    .def("__getitem__", FunctionPointer( [](SPCF self, int comp) -> SPCF
                                         {
                                           if (comp < 0 || comp >= self->Dimension())
                                             bp::exec("raise IndexError()\n");
                                           return make_shared<ComponentCoefficientFunction> (self, comp); 
                                         }))
    .def("__getitem__", FunctionPointer( [](SPCF self, bp::tuple comps) -> SPCF
                                         {
                                           if (bp::len(comps) != 2)
                                             bp::exec("raise IndexError()\n");
                                           Array<int> dims = self->Dimensions();
                                           if (dims.Size() != 2)
                                             bp::exec("raise IndexError()\n");
                                           
                                           int c1 = bp::extract<int> (comps[0]);
                                           int c2 = bp::extract<int> (comps[1]);
                                           if (c1 < 0 || c2 < 0 || c1 >= dims[0] || c2 >= dims[1])
                                             bp::exec("raise IndexError()\n");

                                           int comp = c1 * dims[1] + c2;
                                           return make_shared<ComponentCoefficientFunction> (self, comp);
                                         }))

    // coefficient expressions
    .def ("__add__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF { return c1+c2; } ))
    .def ("__add__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           {
             return coef + make_shared<ConstantCoefficientFunction>(val);
           }))
    .def ("__radd__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return coef + make_shared<ConstantCoefficientFunction>(val); }))

    .def ("__sub__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           { return c1-c2; }))

    .def ("__sub__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return coef - make_shared<ConstantCoefficientFunction>(val); }))

    .def ("__rsub__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return make_shared<ConstantCoefficientFunction>(val) - coef; }))

    .def ("__mul__", FunctionPointer 
          ([] (SPCF c1, SPCF c2) -> SPCF
           {
             if (c1->Dimensions().Size() == 2 && c2->Dimensions().Size() == 2)
               return make_shared<MultMatMatCoefficientFunction> (c1, c2);
             if (c1->Dimension() > 1 && c2->Dimension() > 1)
               return make_shared<MultVecVecCoefficientFunction> (c1, c2);
             if (c1->Dimension() == 1 && c2->Dimension() > 1)
               return make_shared<MultScalVecCoefficientFunction> (c1, c2);
             if (c1->Dimension() > 1 && c2->Dimension() == 1)
               return make_shared<MultScalVecCoefficientFunction> (c2, c1);
             return BinaryOpCF (c1, c2, 
                                [](double a, double b) { return a*b; },
                                [](Complex a, Complex b) { return a*b; },
                                [](double a, double b, double & dda, double & ddb) { dda = b; ddb = a; },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = 1; ddbdb = 0; }
                                );
           } ))

    .def ("InnerProduct", FunctionPointer
          ([] (SPCF c1, SPCF c2) -> SPCF
           { 
             return make_shared<MultVecVecCoefficientFunction> (c1, c2);
           }))
          

    /*
      // it's using the complex functions anyway ...
    .def ("__mul__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { 
             return make_shared<ScaleCoefficientFunction> (val, coef); 
           }))
    .def ("__rmul__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return make_shared<ScaleCoefficientFunction> (val, coef); }))
    */
    .def ("__mul__", FunctionPointer 
          ([] (SPCF coef, Complex val) -> SPCF
           { 
             if (val.imag() == 0)
               return make_shared<ScaleCoefficientFunction> (val.real(), coef); 
             else
               return make_shared<ScaleCoefficientFunctionC> (val, coef); 
           }))
    .def ("__rmul__", FunctionPointer 
          ([] (SPCF coef, Complex val) -> SPCF
           { 
             if (val.imag() == 0)
               return make_shared<ScaleCoefficientFunction> (val.real(), coef); 
             else
               return make_shared<ScaleCoefficientFunctionC> (val, coef); 
           }))

    .def ("__truediv__", FunctionPointer 
          ([] (SPCF coef, SPCF coef2) -> SPCF
           { return BinaryOpCF (coef, coef2,
                                [](double a, double b) { return a/b; },
                                [](Complex a, Complex b) { return a/b; },
                                [](double a, double b, double & dda, double & ddb) { dda = 1.0/b; ddb = -a/(b*b); },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = -1.0/(b*b); ddbdb = 2*a/(b*b*b); }
                                ); }))

    .def ("__truediv__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return a/b; },
                                [](Complex a, Complex b) { return a/b; },
                                [](double a, double b, double & dda, double & ddb) { dda = 1.0/b; ddb = -a/(b*b); },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 0; ddadb = -1.0/(b*b); ddbdb = 2*a/(b*b*b); }
                                ); }))


    .def ("__rtruediv__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return BinaryOpCF (coef, make_shared<ConstantCoefficientFunction>(val), 
                                [](double a, double b) { return b/a; },
                                [](Complex a, Complex b) { return b/a; },
                                [](double a, double b, double & dda, double & ddb) { dda = -b/(a*a); ddb = 1.0/a; },
                                [](double a, double b, double & ddada, double & ddadb, double & ddbdb) 
                                { ddada = 2*b/(a*a*a); ddadb = -b/(a*a); ddbdb = 0; }
                                ); }))

    .def ("__neg__", FunctionPointer 
          ([] (SPCF coef) -> SPCF
           { return make_shared<ScaleCoefficientFunction> (-1, coef); }))

    .add_property ("trans", FunctionPointer
                   ([] (SPCF coef) -> SPCF
                    { return make_shared<TransposeCoefficientFunction> (coef); }))

    .def ("Compile", FunctionPointer
          ([] (SPCF coef) -> SPCF
           {
             return Compile (coef);
           }))
    ;

  ExportStdMathFunction<GenericSin>("sin");
  ExportStdMathFunction<GenericCos>("cos");
  ExportStdMathFunction<GenericTan>("tan");
  ExportStdMathFunction<GenericExp>("exp");
  ExportStdMathFunction<GenericLog>("log");
  ExportStdMathFunction<GenericSqrt>("sqrt");
  ExportStdMathFunction<GenericConj>("Conj");

  
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


  class MeshSizeCF : public CoefficientFunction
  {
  public:
    MeshSizeCF () { ; }
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
    {
      return pow(ip.GetMeasure(), 1.0/ip.Dim());
    }
  };


  class SpecialCoefficientFunctions
  {
  public:
    shared_ptr<CoefficientFunction> GetMeshSizeCF ()
    { return make_shared<MeshSizeCF>(); }

    shared_ptr<CoefficientFunction> GetNormalVectorCF (int dim)
    { 
      if (dim == 2)
        return make_shared<NormalVectorCF<2>>(); 
      else
        return make_shared<NormalVectorCF<3>>(); 
    }
  };

  bp::class_<SpecialCoefficientFunctions> ("SpecialCFCreator", bp::no_init)
    .add_property("mesh_size", 
                  &SpecialCoefficientFunctions::GetMeshSizeCF)
    .def("normal", &SpecialCoefficientFunctions::GetNormalVectorCF)
    ;
  static SpecialCoefficientFunctions specialcf;
  
  bp::scope().attr("specialcf") = bp::object(bp::ptr(&specialcf));

  bp::class_<BSpline> ("BSpline", bp::no_init)
   .def("__init__", bp::make_constructor 
        (FunctionPointer ([](int order, bp::list knots, bp::list vals)
                           {
                             return new BSpline (order, 
                                                 makeCArray<double> (knots),
                                                 makeCArray<double> (vals));
                           })))
    .def("__call__", &BSpline::Evaluate)
    .def("__call__", FunctionPointer
         ([](const BSpline & sp, SPCF coef) -> SPCF 
          {
            return UnaryOpCF (coef,
                              sp,
                              [](Complex a) { return 0; });                             
          }))
    .def("Integrate", 
         FunctionPointer([](const BSpline & sp) { return new BSpline(sp.Integrate()); }),
         bp::return_value_policy<bp::manage_new_object>()
         )
    .def("Differentiate", 
         FunctionPointer([](const BSpline & sp) { return new BSpline(sp.Differentiate()); }),
         bp::return_value_policy<bp::manage_new_object>()
         )
    ;
}





// *************************************** Export FEM ********************************


void NGS_DLL_HEADER ExportNgfem() {
    std::string nested_name = "fem";
    if( bp::scope() )
      nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".fem");

    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << IM(1) << "exporting fem as " << nested_name << endl;
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

  bp::enum_<NODE_TYPE>("NODE_TYPE")
    .value("VERTEX", NT_VERTEX)
    .value("EDGE", NT_EDGE)
    .value("FACE", NT_FACE)
    .value("CELL", NT_CELL)
    .value("ELEMENT", NT_ELEMENT)
    .value("FACET", NT_FACET)
    .export_values()
    ;


  bp::class_<ElementTopology> ("ElementTopology", bp::init<ELEMENT_TYPE>())
    .add_property("name", 
                  static_cast<const char*(ElementTopology::*)()> (&ElementTopology::GetElementName))
    .add_property("vertices", FunctionPointer([](ElementTopology & self)
                                              {
                                                bp::list verts;
                                                ELEMENT_TYPE et = ET_TRIG;
                                                const POINT3D * pts = self.GetVertices(et);
                                                int dim = self.GetSpaceDim(et);
                                                for (int i : Range(self.GetNVertices(et)))
                                                  {
                                                    switch (dim)
                                                      {
                                                      case 2: 
                                                        {
                                                          bp::list v; v.append(pts[i][0]); v.append(pts[i][1]);
                                                          verts.append(bp::tuple(v)); break;
                                                        }
                                                      default:
                                                        ;
                                                      }
                                                  }
                                                return verts;
                                              }));
    ;
    
  bp::class_<FiniteElement, shared_ptr<FiniteElement>, boost::noncopyable>
    ("FiniteElement", "any finite element", bp::no_init)
    .add_property("ndof", &FiniteElement::GetNDof, "number of degrees of freedom of element")    
    .add_property("order", &FiniteElement::Order, "maximal polynomial order of element")    
    .add_property("type", &FiniteElement::ElementType, "geometric type of element")    
    .add_property("dim", &FiniteElement::Dim, "spatial dimension of element")    
    .add_property("classname", &FiniteElement::ClassName, "name of element family")  
    .def("__str__", &ToString<FiniteElement>)
    ;

  bp::class_<BaseScalarFiniteElement, shared_ptr<BaseScalarFiniteElement>, 
    bp::bases<FiniteElement>, boost::noncopyable>
      ("ScalarFE", "a scalar-valued finite element", bp::no_init)
    .def("CalcShape",
         FunctionPointer
         ([] (const BaseScalarFiniteElement & fe, double x, double y, double z)
          {
            IntegrationPoint ip(x,y,z);
            Vector<> v(fe.GetNDof());
            /*
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
            */
            fe.CalcShape (ip, v);
            return v;
          }),
         (bp::arg("self"),bp::arg("x"),bp::arg("y")=0.0,bp::arg("z")=0.0)
         )
    .def("CalcDShape",
         FunctionPointer
         ([] (const BaseScalarFiniteElement & fe, const BaseMappedIntegrationPoint & mip)
          {
            Matrix<> mat(fe.GetNDof(), fe.Dim());
            switch (fe.Dim())
              {
              case 1:
                dynamic_cast<const ScalarFiniteElement<1>&> (fe).
                  CalcMappedDShape(static_cast<const MappedIntegrationPoint<1,1>&> (mip), mat); break;
              case 2:
                dynamic_cast<const ScalarFiniteElement<2>&> (fe).
                  CalcMappedDShape(static_cast<const MappedIntegrationPoint<2,2>&> (mip), mat); break;
              case 3:
                dynamic_cast<const ScalarFiniteElement<3>&> (fe).
                  CalcMappedDShape(static_cast<const MappedIntegrationPoint<3,3>&> (mip), mat); break;
              default:
                ;
              }
            return mat;
          }),
         (bp::arg("self"),bp::arg("mip"))
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
           }),
          "creates an H1 finite element of given geometric shape and polynomial order"
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


  bp::class_<BaseMappedIntegrationPoint, boost::noncopyable>( "BaseMappedIntegrationPoint", bp::no_init)
    .def("__str__", FunctionPointer
         ([] (const BaseMappedIntegrationPoint & bmip)
          {
            stringstream str;
            str << "p = " << bmip.GetPoint() << endl;
            switch (bmip.Dim())
              {
              case 1: 
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<1,1>&>(bmip);
                  str << "jac = " << mip.GetJacobian() << endl;
                  break;
                }
              case 2: 
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<2,2>&>(bmip);
                  str << "jac = " << mip.GetJacobian() << endl;
                  break;
                }
              case 3: 
                {
                  auto & mip = static_cast<const MappedIntegrationPoint<3,3>&>(bmip);
                  str << "jac = " << mip.GetJacobian() << endl;
                  break;
                }
              default:
                ;
              }
            str << "measure = " << bmip.GetMeasure() << endl;
            return str.str();
          }))
    .add_property("measure", &BaseMappedIntegrationPoint::GetMeasure)
    ;

    
  bp::class_<ElementTransformation, boost::noncopyable>("ElementTransformation", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](ELEMENT_TYPE et, bp::object vertices) -> ElementTransformation*
                           {
                             int nv = ElementTopology::GetNVertices(et);
                             int dim = bp::len(vertices[0]);
                             Matrix<> pmat(nv,dim);
                             for (int i : Range(nv))
                               for (int j : Range(dim))
                                 pmat(i,j) = bp::extract<double> (vertices[i][j])();
                             switch (Dim(et))
                               {
                               case 1:
                                 return new FE_ElementTransformation<1,1> (et, pmat); 
                               case 2:
                                 return new FE_ElementTransformation<2,2> (et, pmat); 
                               case 3:
                                 return new FE_ElementTransformation<3,3> (et, pmat); 
                               default:
                                 throw Exception ("cannot create ElementTransformation");
                               }
                           }),
          bp::default_call_policies(),        // need it to use named arguments
          (bp::arg("et")=ET_TRIG,bp::arg("vertices"))))
    .add_property("is_boundary", &ElementTransformation::Boundary)
    .add_property("spacedim", &ElementTransformation::SpaceDim)
    .def ("__call__", FunctionPointer
          ([] (ElementTransformation & self, double x, double y, double z)
           {
             
             return &self(IntegrationPoint(x,y,z), global_alloc);
           }),
          (bp::args("self"), bp::args("x"), bp::args("y")=0, bp::args("z")=0),
          bp::return_value_policy<bp::manage_new_object>())
    ;


  bp::class_<DifferentialOperator, shared_ptr<DifferentialOperator>, boost::noncopyable>
    ("DifferentialOperator", bp::no_init)
    ;

  
  bp::class_<BilinearFormIntegrator, shared_ptr<BilinearFormIntegrator>, boost::noncopyable>
    ("BFI", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](string name, int dim, bp::object py_coef, bp::object definedon, bool imag)
                           {
                             Array<shared_ptr<CoefficientFunction>> coef = MakeCoefficients(py_coef);
                             auto bfi = GetIntegrators().CreateBFI (name, dim, coef);

                             if (!bfi) cerr << "undefined integrator '" << name 
                                            << "' in " << dim << " dimension" << endl;

                             if (bp::extract<bp::list> (definedon).check())
                               bfi -> SetDefinedOn (makeCArray<int> (definedon));

                             if (imag)
                               bfi = make_shared<ComplexBilinearFormIntegrator> (bfi, Complex(0,1));

                             return bfi;
                           }),
          bp::default_call_policies(),        // need it to use named arguments
          (bp::arg("name")=NULL,bp::arg("dim")=-1,bp::arg("coef"),
           bp::arg("definedon")=bp::object(),bp::arg("imag")=false)))
    
    .def("__str__", &ToString<BilinearFormIntegrator>)

    .def("Evaluator", &BilinearFormIntegrator::GetEvaluator)
    .def("DefinedOn", &Integrator::DefinedOn)

    /*
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
    */

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

                             if (!lfi) throw Exception(string("undefined integrator '")+name+
                                                       "' in "+ToString(dim)+ " dimension having 1 coefficient");
                             
                             if (bp::extract<bp::list> (definedon).check())
                               lfi -> SetDefinedOn (makeCArray<int> (definedon));
 
                             if (imag)
                               lfi = make_shared<ComplexLinearFormIntegrator> (lfi, Complex(0,1));

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

    .def("DefinedOn", &Integrator::DefinedOn)
    .def("CalcElementVector", 
         static_cast<void(LinearFormIntegrator::*)(const FiniteElement&, const ElementTransformation&, FlatVector<double>,LocalHeap&)const>
         (&LinearFormIntegrator::CalcElementVector))
    .def("CalcElementVector",
         FunctionPointer([] (const LinearFormIntegrator & self, 
                             const FiniteElement & fe, const ElementTransformation & trafo,
                             int heapsize)
                         {
                           Vector<> vec(fe.GetNDof());
                           while (true)
                             {
                               try
                                 {
                                   LocalHeap lh(heapsize);
                                   self.CalcElementVector (fe, trafo, vec, lh);
                                   return vec;
                                 }
                               catch (LocalHeapOverflow ex)
                                 {
                                   heapsize *= 10;
                                 }
                             };
                         }),
         bp::default_call_policies(),        // need it to use named arguments
         (bp::arg("lfi")=NULL, bp::arg("fel"),bp::arg("trafo"),bp::arg("heapsize")=10000))
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



  ExportCoefficientFunction ();


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





  bp::def ("SetPMLParameters", 
           FunctionPointer([] (double rad, double alpha)
                           {
                             cout << "set pml parameters, r = " << rad << ", alpha = " << alpha << endl;
                             constant_table_for_FEM = &pmlpar;
                             pmlpar.Set("pml_r", rad);
                             pmlpar.Set("pml_alpha", alpha);
                             SetPMLParameters();
                           }),
           (bp::arg("rad")=1,bp::arg("alpha")=1))
    ;
    
                           
                           

}


void ExportNgstd();
void ExportNgbla();

BOOST_PYTHON_MODULE(libngfem) {
  // ExportNgstd();
  // ExportNgbla();
  ExportNgfem();
}



#endif
