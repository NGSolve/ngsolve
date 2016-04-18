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

  bp::extract<Complex> ec(val);
  if (ec.check()) 
    return make_shared<ConstantCoefficientFunctionC> (ec());

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

  // return move(tmp);  // clang recommends not to move it ...
  return tmp;
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
                  return bp::object(UnaryOpCF(coef, func, func, FUNC::Name()));
                }
              bp::extract<double> ed(x);
              if (ed.check()) return bp::object(func(ed()));
              if (bp::extract<Complex> (x).check())
                return bp::object(func(bp::extract<Complex> (x)()));
              throw Exception ("can't compute math-function");
            }));
}



template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> exp (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = std::exp(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * res.Value();
  return res;
}

template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> sin (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = std::sin(x.Value());
  SCAL c = std::cos(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * c;
  return res;
}

template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> cos (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = std::cos(x.Value());
  SCAL ms = -std::sin(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * ms;
  return res;
}


// df(u)/dx  = exp(x) * du/dx
// d^2 f(u) / dx^2 = exp(x) * (du/dx)^2 + exp(x) * d^2u /dx^2
template <int D>
INLINE AutoDiffDiff<D> exp (AutoDiffDiff<D> x)
{
  AutoDiffDiff<D> res;
  res.Value() = std::exp(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * res.Value();
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = (x.DValue(k) * x.DValue(l)+x.DDValue(k,l)) * res.Value();
  return res;
}

template <int D>
INLINE AutoDiffDiff<D> log (AutoDiffDiff<D> x)
{
  AutoDiffDiff<D> res;
  res.Value() = std::log(x.Value());
  double xinv = 1.0/x.Value();
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * xinv;
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = -xinv*xinv*x.DValue(k) * x.DValue(l) + xinv * x.DDValue(k,l);
  return res;
}



template <int D>
INLINE AutoDiffDiff<D> sin (AutoDiffDiff<D> x)
{
  AutoDiffDiff<D> res;
  double s = std::sin(x.Value());
  double c = std::cos(x.Value());
  
  res.Value() = s;
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * c;
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = -s * x.DValue(k) * x.DValue(l) + c * x.DDValue(k,l);
  return res;
}


template <int D>
INLINE AutoDiffDiff<D> cos (AutoDiffDiff<D> x)
{
  AutoDiffDiff<D> res;
  double s = std::sin(x.Value());
  double c = std::cos(x.Value());
  
  res.Value() = c;
  for (int k = 0; k < D; k++)
    res.DValue(k) = -s * x.DValue(k);
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = -c * x.DValue(k) * x.DValue(l) - s * x.DDValue(k,l);
  return res;
}


template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> tan (AutoDiff<D,SCAL> x)
{ return sin(x) / cos(x); }

template <int D>
INLINE AutoDiffDiff<D> tan (AutoDiffDiff<D> x)
{ return sin(x) / cos(x); }

template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> atan (AutoDiff<D,SCAL> x)
{
  AutoDiff<D> res;
  double a = std::atan(x.Value());
  res.Value() = a;
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k)/(1+x.Value()*x.Value()) ;
  return res;
}


template <int D>
INLINE AutoDiffDiff<D> atan (AutoDiffDiff<D> x)
{
  AutoDiffDiff<D> res;
  double a = std::atan(x.Value());
  res.Value() = a;
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k)/(1+x.Value()*x.Value()) ;
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = -2*x.Value()/((1+x.Value()*x.Value())*(1+x.Value()*x.Value())) * x.DValue(k) * x.DValue(l) + x.DDValue(k,l)/(1+x.Value()*x.Value());
  return res;
}



struct GenericBSpline {
  shared_ptr<BSpline> sp;
  GenericBSpline( const BSpline &asp ) : sp(make_shared<BSpline>(asp)) {;}
  GenericBSpline( shared_ptr<BSpline> asp ) : sp(asp) {;}
  template <typename T> T operator() (T x) const { return (*sp)(x); }
  Complex operator() (Complex x) const { return (*sp)(x.real()); }
};
struct GenericSin {
  template <typename T> T operator() (T x) const { return sin(x); }
  static string Name() { return "sin"; }
};
struct GenericCos {
  template <typename T> T operator() (T x) const { return cos(x); }
  static string Name() { return "cos"; }
};
struct GenericTan {
  template <typename T> T operator() (T x) const { return tan(x); }
  static string Name() { return "tan"; }
};
struct GenericExp {
  template <typename T> T operator() (T x) const { return exp(x); }
  static string Name() { return "exp"; }
};
struct GenericLog {
  template <typename T> T operator() (T x) const { return log(x); }
  static string Name() { return "log"; }
};
struct GenericATan {
  template <typename T> T operator() (T x) const { return atan(x); }
  static string Name() { return "atan"; }
};
struct GenericSqrt {
  template <typename T> T operator() (T x) const { return sqrt(x); }
  static string Name() { return "sqrt"; }
};
struct GenericConj {
  template <typename T> T operator() (T x) const { return Conj(x); } // from bla
  static string Name() { return "conj"; }
  AutoDiff<1> operator() (AutoDiff<1> x) const { throw Exception ("Conj(..) is not complex differentiable"); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception ("Conj(..) is not complex differentiable"); }
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

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<> res) const 
    {
      if (ir[0].Dim() != D)
        throw Exception("illegal dim of normal vector");
      FlatMatrixFixWidth<D> resD(res);
      for (int i = 0; i < ir.Size(); i++)
        resD.Row(i) = static_cast<const DimMappedIntegrationPoint<D>&>(ir[i]).GetNV();
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
        string miptype;
        if(code.is_simd)
          miptype = "SIMD<DimMappedIntegrationPoint<"+ToString(D)+">>*";
        else
          miptype = "DimMappedIntegrationPoint<"+ToString(D)+">*";
        auto nv_expr = CodeExpr("static_cast<const "+miptype+">(&ip)->GetNV()");
        auto nv = Var("tmp", index);
        code.body += nv.Assign(nv_expr);
        for( int i : Range(D))
          code.body += Var(index,i).Assign(nv(i));
    }
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
    {
      // static Timer t("NormalVec::EvalSIMD"); RegionTimer reg(t);            
      auto & ir22 = static_cast<const SIMD_MappedIntegrationRule<2,2>&> (ir);
      for (int i = 0; i < ir.Size(); i++)
        for (int j = 0; j < D; j++)
          values.Get(j,i) = ir22[i].GetNV()(j).Data();
    }
    
  };

  template <int D>
  class TangentialVectorCF : public CoefficientFunction
  {
  public:
    TangentialVectorCF () { ; }
    virtual int Dimension() const { return D; }

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
    {
      return 0;
    }
    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> res) const 
    {
      if (ip.Dim() != D)
        throw Exception("illegal dim of tangential vector");
      res = static_cast<const DimMappedIntegrationPoint<D>&>(ip).GetTV();
    }
  };




void ExportCoefficientFunction()
{
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<CoefficientFunction>);
  bp::class_<CoefficientFunction, shared_ptr<CoefficientFunction>, boost::noncopyable> 
    ("CoefficientFunction",
     "A CoefficientFunction (CF) is some function defined on a mesh.\n"
     "examples are coordinates x, y, z, domain-wise constants, solution-fields, ...\n"
     "CFs can be combined by mathematical operations (+,-,sin(), ...) to form new CFs",
     bp::no_init)

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
          (bp::arg("coef"),bp::arg("dims")=bp::object())),
         "Construct a CoefficientFunction from either one of\n"
         "  a scalar (float or complex)\n"
         "  a tuple of scalars and or CFs to define a vector-valued CF\n"
         "     use dims=(h,w) to define matrix-valued CF\n"
         "  a list of scalars and or CFs to define a domain-wise CF"
         )
    .def("__str__", &ToString<CoefficientFunction>)

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
                Vector<Complex> vec(self.Dimension());
                self.Evaluate (mip, vec);
                if (vec.Size()==1) return bp::object(vec(0));
                return bp::tuple(vec);
	      }
	  }),
         (bp::arg("coef"),bp::arg("mip")),
         "evaluate CF at a mapped integrationpoint mip. mip can be generated by calling mesh(x,y,z)")
    .add_property("dim", &CoefficientFunction::Dimension,
                  "number of components of CF")

    .add_property("dims", &CoefficientFunction::Dimensions,
                  "shape of CF:  (dim) for vector, (h,w) for matrix")    
    
    .def("__getitem__", FunctionPointer( [](SPCF self, int comp) -> SPCF
                                         {
                                           if (comp < 0 || comp >= self->Dimension())
                                             bp::exec("raise IndexError()\n");
                                           return make_shared<ComponentCoefficientFunction> (self, comp); 
                                         }),
         (bp::arg("coef"),bp::arg("comp")),         
         "returns component comp of vectorial CF")
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
             return c1*c2;
           } ))

    .def ("InnerProduct", FunctionPointer
          ([] (SPCF c1, SPCF c2) -> SPCF
           { 
             return InnerProduct (c1, c2);
           }))
          
    .def("Norm", FunctionPointer ( [](SPCF x) { return NormCF(x); }))

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
               return val.real() * coef;
             else
               return val * coef;
           }))
    .def ("__rmul__", FunctionPointer 
          ([] (SPCF coef, Complex val) -> SPCF
           { 
             if (val.imag() == 0)
               return val.real() * coef; 
             else
               return val * coef;
           }))

    .def ("__truediv__", FunctionPointer 
          ([] (SPCF coef, SPCF coef2) -> SPCF
           { return coef/coef2;
           }))

    .def ("__truediv__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return coef / make_shared<ConstantCoefficientFunction>(val); }))

    .def ("__rtruediv__", FunctionPointer 
          ([] (SPCF coef, double val) -> SPCF
           { return make_shared<ConstantCoefficientFunction>(val) / coef; }))

    .def ("__neg__", FunctionPointer 
          ([] (SPCF coef) -> SPCF
           { return -1.0 * coef; }))

    .add_property ("trans", FunctionPointer
                   ([] (SPCF coef) -> SPCF
                    {
                      return TransposeCF(coef);
                    }),
                   "transpose of matrix-valued CF")

    .def ("Compile", FunctionPointer
          ([] (SPCF coef, bool realcompile) -> SPCF
           { return Compile (coef, realcompile); }),
          (bp::args("self"), bp::args("realcompile")=false),
          "compile list of individual steps, experimental improvement for deep trees")
    ;

  ExportStdMathFunction<GenericSin>("sin");
  ExportStdMathFunction<GenericCos>("cos");
  ExportStdMathFunction<GenericTan>("tan");
  ExportStdMathFunction<GenericExp>("exp");
  ExportStdMathFunction<GenericLog>("log");
  ExportStdMathFunction<GenericATan>("atan");
  ExportStdMathFunction<GenericSqrt>("sqrt");
  ExportStdMathFunction<GenericConj>("Conj");
  
  bp::def ("IfPos", FunctionPointer 
           ([] (SPCF c1, bp::object then_obj, bp::object else_obj) -> SPCF
            {
              return IfPos(c1,
                           MakeCoefficient(then_obj),
                           MakeCoefficient(else_obj));
            } ))
    ;
  
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<ConstantCoefficientFunction>);
  bp::class_<ConstantCoefficientFunction,bp::bases<CoefficientFunction>,
    shared_ptr<ConstantCoefficientFunction>, boost::noncopyable>
    ("ConstantCF", "same as CoefficientFunction(c), obsolete", bp::init<double>())
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
      result.Col(0) = ir.GetPoints().Col(dir);
      /*
      for (int i = 0; i < ir.Size(); i++)
        result(i,0) = ir[i].GetPoint()(dir);
      */
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
        auto v = Var(index);
        if(dir==0) code.body += v.Assign(CodeExpr("ip.GetPoint()(0)"));
        if(dir==1) code.body += v.Assign(CodeExpr("ip.GetPoint()(1)"));
        if(dir==2) code.body += v.Assign(CodeExpr("ip.GetPoint()(2)"));
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
    {
      // static Timer t("CoordCF::EvalSIMD"); RegionTimer reg(t);      
      // auto & ir22 = static_cast<const SIMD_MappedIntegrationRule<2,2>&> (ir);
      auto points = ir.GetPoints();
      for (int i = 0; i < ir.Size(); i++)
        values.Get(i) = points.Get(i, dir);
        // values.Get(i) = ir22[i].Point()(dir).Data();
    }
  };

  bp::class_<CoordCoefficientFunction,bp::bases<CoefficientFunction>,
    shared_ptr<CoordCoefficientFunction>, boost::noncopyable>
    ("CoordCF", "CoefficientFunction for x, y, z", bp::init<int>())
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

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, AFlatMatrix<double> values) const
    {
      for(int i : Range(ir))
        values.Get(i) = pow(ir[i].GetMeasure(), 1.0/2.0).Data();//TODO FIX!!!!  (DIM=2.0)
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
      code.body += Var(index).Assign( CodeExpr("pow(ip.GetMeasure(), 1.0/ip.Dim())"));
    }
  };


  class SpecialCoefficientFunctions
  {
  public:
    shared_ptr<CoefficientFunction> GetMeshSizeCF ()
    { return make_shared<MeshSizeCF>(); }

    shared_ptr<CoefficientFunction> GetNormalVectorCF (int dim)
    { 
      switch(dim)
	{
	case 1:
	  return make_shared<NormalVectorCF<1>>(); 
	case 2:
	  return make_shared<NormalVectorCF<2>>(); 
	default:
	  return make_shared<NormalVectorCF<3>>(); 
	}
    }

    shared_ptr<CoefficientFunction> GetTangentialVectorCF (int dim)
    { 
      switch(dim)
	{
	case 1:
	  return make_shared<TangentialVectorCF<1>>(); 
	case 2:
	  return make_shared<TangentialVectorCF<2>>(); 
	default:
	  return make_shared<TangentialVectorCF<3>>(); 
	}
    }
  };

  bp::class_<SpecialCoefficientFunctions> ("SpecialCFCreator", bp::no_init)
    .add_property("mesh_size", 
                  &SpecialCoefficientFunctions::GetMeshSizeCF, "local mesh-size (approximate element diameter) as CF")
    .def("normal", &SpecialCoefficientFunctions::GetNormalVectorCF,
         "depending on contents: normal-vector to geometry or element\n"
         "space-dimension must be provided")
    .def("tangential", &SpecialCoefficientFunctions::GetTangentialVectorCF,
         "depending on contents: tangential-vector to element\n"
         "space-dimension must be provided")
    ;
  static SpecialCoefficientFunctions specialcf;
  
  bp::scope().attr("specialcf") = bp::object(bp::ptr(&specialcf));

  bp::class_<BSpline, shared_ptr<BSpline> > ("BSpline", bp::no_init)
   .def("__init__", bp::make_constructor 
        (FunctionPointer ([](int order, bp::list knots, bp::list vals)
                           {
                             return make_shared<BSpline> (order,
                                                 makeCArray<double> (knots),
                                                 makeCArray<double> (vals));
                           })),
        "B-Spline of a certain order, provide knot and value vectors")
    .def("__str__", &ToString<BSpline>)
    .def("__call__", &BSpline::Evaluate)
    .def("__call__", FunctionPointer
         ([](shared_ptr<BSpline> sp, SPCF coef) -> SPCF
          {
            return UnaryOpCF (coef, GenericBSpline(sp), GenericBSpline(sp));
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

  bp::docstring_options local_docstring_options(true, true, false);
  
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
                                                const POINT3D * pts = self.GetVertices();
                                                int dim = self.GetSpaceDim();
                                                for (int i : Range(self.GetNVertices()))
                                                  {
                                                    bp::list v;
                                                    for (int j = 0; j < dim; j++)
                                                      v.append(pts[i][j]);
                                                    verts.append (bp::tuple(v));
                                                  }
                                                return verts;
                                              }));
    ;
    
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<FiniteElement>);
  bp::class_<FiniteElement, shared_ptr<FiniteElement>, boost::noncopyable>
    ("FiniteElement", "any finite element", bp::no_init)
    .add_property("ndof", &FiniteElement::GetNDof, "number of degrees of freedom of element")    
    .add_property("order", &FiniteElement::Order, "maximal polynomial order of element")    
    .add_property("type", &FiniteElement::ElementType, "geometric type of element")    
    .add_property("dim", &FiniteElement::Dim, "spatial dimension of element")    
    .add_property("classname", &FiniteElement::ClassName, "name of element family")  
    .def("__str__", &ToString<FiniteElement>)
    ;

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<BaseScalarFiniteElement>);
  bp::class_<BaseScalarFiniteElement, shared_ptr<BaseScalarFiniteElement>, 
    bp::bases<FiniteElement>, boost::noncopyable>
      ("ScalarFE", "a scalar-valued finite element", bp::no_init)
    .def("CalcShape",
         FunctionPointer
         ([] (const BaseScalarFiniteElement & fe, double x, double y, double z)
          {
            IntegrationPoint ip(x,y,z);
            Vector<> v(fe.GetNDof());
            fe.CalcShape (ip, v);
            return v;
          }),
         (bp::arg("self"),bp::arg("x"),bp::arg("y")=0.0,bp::arg("z")=0.0)
         )
    .def("CalcShape",
         FunctionPointer
         ([] (const BaseScalarFiniteElement & fe, const BaseMappedIntegrationPoint & mip)
          {
            Vector<> v(fe.GetNDof());
            fe.CalcShape (mip.IP(), v);
            return v;
          }),
         (bp::arg("self"),bp::arg("mip"))
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


  bp::class_<IntegrationPoint>("IntegrationPoint", bp::no_init);

  bp::class_<IntegrationRule, boost::noncopyable>("IntegrationRule", bp::no_init)
    .def("__init__",  bp::make_constructor 
         (FunctionPointer ([](ELEMENT_TYPE et, int order) -> IntegrationRule *
                           {
                             return new IntegrationRule (et, order);
                           }),
          bp::default_call_policies(),
          (bp::arg("element type"), bp::arg("order"))))
    .def("__getitem__", FunctionPointer([](IntegrationRule & ir, int nr)->IntegrationPoint
                                        {
                                          if (nr < 0 || nr >= ir.Size())
                                            bp::exec("raise IndexError()\n");
                                          return ir[nr];
                                        }))
    .def("Integrate", FunctionPointer
         ([](IntegrationRule & ir, bp::object func) -> bp::object
          {
            bp::object sum;
            for (const IntegrationPoint & ip : ir)
              {
                bp::object val;
                switch (ir.Dim())
                  {
                  case 1:
                    val = func(ip(0)); break;
                  case 2:
                    val = func(ip(0), ip(1)); break;
                  case 3:
                    val = func(ip(0), ip(1), ip(2)); break;
                  default:
                    throw Exception("integration rule with illegal dimension");
                  }

                val = val * ip.Weight();
                if (sum == bp::object())
                  sum = val;
                else
                  sum = sum+val;
              }
            return sum;
          }))
    ;

  bp::class_<BaseMappedIntegrationPoint, boost::noncopyable>( "BaseMappedIntegrationPoint", bp::no_init)
    .def("__str__", FunctionPointer
         ([] (const BaseMappedIntegrationPoint & bmip)
          {
            stringstream str;
            str << "p = " << bmip.GetPoint() << endl;
            str << "jac = " << bmip.GetJacobian() << endl;
            /*
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
            */
            str << "measure = " << bmip.GetMeasure() << endl;
            return str.str();
          }))
    .add_property("measure", &BaseMappedIntegrationPoint::GetMeasure)
    .add_property("point", &BaseMappedIntegrationPoint::GetPoint)
    .add_property("jacobi", &BaseMappedIntegrationPoint::GetJacobian)
    // .add_property("trafo", &BaseMappedIntegrationPoint::GetTransformation)
    .add_property("trafo", bp::make_function( &BaseMappedIntegrationPoint::GetTransformation,
                                              bp::return_value_policy<bp::reference_existing_object>()))
    .add_property("elementid", FunctionPointer([](BaseMappedIntegrationPoint & mip)
                                               {
                                                 return mip.GetTransformation().GetElementId();
                                               }))
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
    .add_property("elementid", &ElementTransformation::GetElementId)
    .def ("__call__", FunctionPointer
          ([] (ElementTransformation & self, double x, double y, double z)
           {
             
             return &self(IntegrationPoint(x,y,z), global_alloc);
           }),
          (bp::args("self"), bp::args("x"), bp::args("y")=0, bp::args("z")=0),
          bp::return_value_policy<bp::manage_new_object>())
    .def ("__call__", FunctionPointer
          ([] (ElementTransformation & self, IntegrationPoint & ip)
           {
             return &self(ip, global_alloc);
           }),
          (bp::args("self"), bp::args("ip")),
          bp::return_value_policy<bp::manage_new_object>())
    ;


  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<DifferentialOperator>);
  bp::class_<DifferentialOperator, shared_ptr<DifferentialOperator>, boost::noncopyable>
    ("DifferentialOperator", bp::no_init)
    ;

  
  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<BilinearFormIntegrator>);
  bp::class_<BilinearFormIntegrator, shared_ptr<BilinearFormIntegrator>, boost::noncopyable>
    ("BFI", bp::no_init)
    .def("__init__", bp::make_constructor
         (FunctionPointer ([](string name, int dim, bp::object py_coef, bp::object definedon, bool imag,
                              string filename)
                           {
                             Array<shared_ptr<CoefficientFunction>> coef = MakeCoefficients(py_coef);
                             auto bfi = GetIntegrators().CreateBFI (name, dim, coef);

                             if (!bfi) cerr << "undefined integrator '" << name 
                                            << "' in " << dim << " dimension" << endl;

                             if (bp::extract<bp::list> (definedon).check())
                               bfi -> SetDefinedOn (makeCArray<int> (definedon));

                             if (bp::extract<BitArray> (definedon).check())
                               {
                                 // cout << "defon with bitarray: " << flush;
                                 // cout << bp::extract<BitArray> (definedon)() << endl;
                                 bfi -> SetDefinedOn (bp::extract<BitArray> (definedon)());
                               }

                             if (filename.length())
                               {
                                 cout << "set integrator filename: " << filename << endl;
                                 bfi -> SetFileName (filename);
                               }
                             if (imag)
                               bfi = make_shared<ComplexBilinearFormIntegrator> (bfi, Complex(0,1));

                             return bfi;
                           }),
          bp::default_call_policies(),        // need it to use named arguments
          (bp::arg("name")=NULL,bp::arg("dim")=-1,bp::arg("coef"),
           bp::arg("definedon")=bp::object(),bp::arg("imag")=false, bp::arg("filename")="")))
    
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


  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<CompoundBilinearFormIntegrator>);
  bp::class_<CompoundBilinearFormIntegrator,bp::bases<BilinearFormIntegrator>,
             shared_ptr<CompoundBilinearFormIntegrator>, boost::noncopyable>
      ("CompoundBilinearFormIntegrator", bp::no_init);

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<BlockBilinearFormIntegrator>);
  bp::class_<BlockBilinearFormIntegrator,bp::bases<BilinearFormIntegrator>,
             shared_ptr<BlockBilinearFormIntegrator>, boost::noncopyable>
      ("BlockBilinearFormIntegrator", bp::no_init);

  bp::implicitly_convertible 
    <shared_ptr<CompoundBilinearFormIntegrator>, 
    shared_ptr<BilinearFormIntegrator> >(); 

  bp::implicitly_convertible 
    <shared_ptr<BlockBilinearFormIntegrator>, 
    shared_ptr<BilinearFormIntegrator> >(); 


  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<LinearFormIntegrator>);
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


  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<CompoundLinearFormIntegrator>);
  bp::class_<CompoundLinearFormIntegrator,bp::bases<LinearFormIntegrator>,
             shared_ptr<CompoundLinearFormIntegrator>, boost::noncopyable>
      ("CompoundLinearFormIntegrator", bp::no_init);

  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<BlockLinearFormIntegrator>);
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


  REGISTER_PTR_TO_PYTHON_BOOST_1_60_FIX(shared_ptr<DomainConstantCoefficientFunction>);
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
