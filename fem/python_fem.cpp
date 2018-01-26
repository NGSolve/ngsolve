#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include "../ngstd/bspline.hpp"
#include <fem.hpp>
#include <mutex>
using namespace ngfem;
using ngfem::ELEMENT_TYPE;

#include "pml.hpp"

#include "tpintrule.hpp"
namespace ngfem
{
  extern SymbolTable<double> * constant_table_for_FEM;
  SymbolTable<double> pmlpar;
}


struct PythonCoefficientFunction : public CoefficientFunction {
  PythonCoefficientFunction() : CoefficientFunction(1,false) { ; }

    virtual double EvaluateXYZ (double x, double y, double z) const = 0;
  
    py::list GetCoordinates(const BaseMappedIntegrationPoint &bip ) {
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
        py::list list;
        int i;
        for(i=0; i<dim; i++)
            list.append(py::float_(x[i]));
        for(i=0; i<3; i++)
            list.append(py::float_(0.0));
        return list;
    }
};

typedef CoefficientFunction CF;

shared_ptr<CF> MakeCoefficient (py::object val)
{
  py::extract<shared_ptr<CF>> ecf(val);
  if (ecf.check()) return ecf();
  
  // a numpy.complex converts itself to a real, and prints a warning
  // thus we check for it first
  if (string(py::str(val.get_type())) == "<class 'numpy.complex128'>")
    return make_shared<ConstantCoefficientFunctionC> (val.cast<Complex>());

  if(py::CheckCast<double>(val))
    return make_shared<ConstantCoefficientFunction> (val.cast<double>());
  if(py::CheckCast<Complex>(val)) 
    return make_shared<ConstantCoefficientFunctionC> (val.cast<Complex>());

  if (py::isinstance<py::list>(val))
    {
      py::list el(val);
      Array<shared_ptr<CoefficientFunction>> cflist(py::len(el));
      for (int i : Range(cflist))
        cflist[i] = MakeCoefficient(el[i]);
      return MakeDomainWiseCoefficientFunction(move(cflist));
    }

  if (py::isinstance<py::tuple>(val))
    {
      py::tuple et(val);
      Array<shared_ptr<CoefficientFunction>> cflist(py::len(et));
      for (int i : Range(cflist))
        cflist[i] = MakeCoefficient(et[i]);
      return MakeVectorialCoefficientFunction(move(cflist));
    }


  throw Exception ("cannot make coefficient");
}

Array<shared_ptr<CoefficientFunction>> MakeCoefficients (py::object py_coef)
{
  Array<shared_ptr<CoefficientFunction>> tmp;
  if (py::isinstance<py::list>(py_coef))
    {
      auto l = py_coef.cast<py::list>();
      for (int i = 0; i < py::len(l); i++)
        tmp += MakeCoefficient(l[i]);
    }
  else if (py::isinstance<py::tuple>(py_coef))
    {
      auto l = py_coef.cast<py::tuple>();
      for (int i = 0; i < py::len(l); i++)
        tmp += MakeCoefficient(l[i]);
    }
  else
    tmp += MakeCoefficient(py_coef);

  // return move(tmp);  // clang recommends not to move it ...
  return tmp;
}



template <typename FUNC>
void ExportStdMathFunction(py::module &m, string name)
{
  m.def (name.c_str(), 
           [] (py::object x) -> py::object
            {
              FUNC func;
              if (py::extract<shared_ptr<CF>>(x).check())
                {
                  auto coef = py::extract<shared_ptr<CF>>(x)();
                  return py::cast(UnaryOpCF(coef, func, /* func, */ FUNC::Name()));
                }
              py::extract<double> ed(x);
              if (ed.check()) return py::cast(func(ed()));
              if (py::extract<Complex> (x).check())
                return py::cast(func(py::extract<Complex> (x)()));
              throw py::type_error (string("can't compute math-function, type = ")
                                    + typeid(FUNC).name());
            });
}


namespace ngfem
{
  void ExportUnaryFunction2 (py::module & m, string name,
                             std::function<shared_ptr<CoefficientFunction>(shared_ptr<CoefficientFunction>)> creator,
                             std::function<double(double)> func_real,
                             std::function<Complex(Complex)> func_complex)
  {
    m.def (name.c_str(),
           [creator, func_real, func_complex] (py::object x) -> py::object
           {
             if (py::extract<shared_ptr<CF>>(x).check())
               {
                 auto coef = py::extract<shared_ptr<CF>>(x)();
                 return py::cast(creator(coef));
             }
             
             py::extract<double> ed(x);
             if (ed.check()) return py::cast(func_real(ed()));
             if (py::extract<Complex> (x).check())
               return py::cast(func_complex(py::extract<Complex> (x)()));
             
             throw py::type_error ("can't compute unary math-function");
           });         
  }


  void ExportBinaryFunction2 (py::module & m, string name,
                              std::function<shared_ptr<CoefficientFunction>(shared_ptr<CoefficientFunction>,
                                                                            shared_ptr<CoefficientFunction>)> creator,
                              std::function<double(double,double)> func_real,
                              std::function<Complex(Complex,Complex)> func_complex)
  {
    m.def (name.c_str(),
           [creator, func_real, func_complex] (py::object x, py::object y) -> py::object
           {
             if (py::extract<shared_ptr<CF>>(x).check() && py::extract<shared_ptr<CF>>(y).check())
               {
                 auto coefx = py::extract<shared_ptr<CF>>(x)();
                 auto coefy = py::extract<shared_ptr<CF>>(y)();
                 return py::cast(creator(coefx, coefy));
             }
             
             py::extract<double> edx(x);
             py::extract<double> edy(y);
             if (edx.check() && edy.check()) return py::cast(func_real(edx(), edy()));
             if (py::extract<Complex> (x).check() && py::extract<Complex> (y).check())
               return py::cast(func_complex(py::extract<Complex> (x)(), py::extract<Complex> (y)()));
             
             throw py::type_error ("can't compute binary math-function");
           });         
  }


}
                          


template <typename FUNC>
void ExportStdMathFunction2(py::module &m, string name)
{
  m.def (name.c_str(), 
         [] (py::object x, py::object y) -> py::object
         {
           FUNC func;
           if (py::extract<shared_ptr<CF>>(x).check() || py::extract<shared_ptr<CF>>(y).check())
             {
               shared_ptr<CoefficientFunction> cx = MakeCoefficient(x);
               shared_ptr<CoefficientFunction> cy = MakeCoefficient(y);
               return py::cast(BinaryOpCF(cx, cy, func,
                                          [](bool a, bool b) { return a||b; }, 'X' /* FUNC::Name() */));
             }
           py::extract<double> dx(x), dy(y);
           if (dx.check() && dy.check()) return py::cast(func(dx(), dy()));
           py::extract<Complex> cx(x), cy(y);
           if (cx.check() && cy.check()) return py::cast(func(cx(), cy()));
           throw py::type_error (string("can't compute binary math-function")+typeid(FUNC).name());
         });
}





struct GenericBSpline {
  shared_ptr<BSpline> sp;
  GenericBSpline( const BSpline &asp ) : sp(make_shared<BSpline>(asp)) {;}
  GenericBSpline( shared_ptr<BSpline> asp ) : sp(asp) {;}
  template <typename T> T operator() (T x) const { return (*sp)(x); }
  Complex operator() (Complex x) const { return (*sp)(x.real()); }
  SIMD<double> operator() (SIMD<double> x) const
  { return SIMD<double>([&](int i)->double { return (*sp)(x[i]); } );}
  SIMD<Complex> operator() (SIMD<Complex> x) const
  { return SIMD<double>([&](int i)->double { return (*sp)(x.real()[i]); } );}
  AutoDiff<1,SIMD<double>> operator() (AutoDiff<1,SIMD<double>> x) const { throw ExceptionNOSIMD ("AutoDiff for bspline not supported"); }
  AutoDiffDiff<1,SIMD<double>> operator() (AutoDiffDiff<1,SIMD<double>> x) const { throw ExceptionNOSIMD ("AutoDiffDiff for bspline not supported"); }  
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
struct GenericFloor {
  template <typename T> T operator() (T x) const { return floor(x); }
  Complex operator() (Complex x) const { throw Exception("no floor for Complex"); }  
  // SIMD<double> operator() (SIMD<double> x) const { throw ExceptionNOSIMD("no floor for simd"); }
  SIMD<Complex> operator() (SIMD<Complex> x) const { throw ExceptionNOSIMD("no floor for simd"); }  
  // AutoDiff<1> operator() (AutoDiff<1> x) const { throw Exception("no floor for AD"); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception("no floor for ADD"); }
  static string Name() { return "floor"; }
};
struct GenericCeil {
  template <typename T> T operator() (T x) const { return ceil(x); }
  Complex operator() (Complex x) const { throw Exception("no ceil for Complex"); }  
  // SIMD<double> operator() (SIMD<double> x) const { throw ExceptionNOSIMD("no ceil for simd"); }
  SIMD<Complex> operator() (SIMD<Complex> x) const { throw ExceptionNOSIMD("no ceil for simd"); }  
  // AutoDiff<1> operator() (AutoDiff<1> x) const { throw Exception("no ceil for AD"); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception("no ceil for ADD"); }
  static string Name() { return "ceil"; }
};

struct GenericConj {
  template <typename T> T operator() (T x) const { return Conj(x); } // from bla
  static string Name() { return "conj"; }
  SIMD<double> operator() (SIMD<double> x) const { return x; }
  template<typename T>
  AutoDiff<1,T> operator() (AutoDiff<1,T> x) const { throw Exception ("Conj(..) is not complex differentiable"); }
  template<typename T>
  AutoDiffDiff<1,T> operator() (AutoDiffDiff<1,T> x) const { throw Exception ("Conj(..) is not complex differentiable"); }
};

struct GenericATan2 {
  double operator() (double x, double y) const { return atan2(x,y); }
  template <typename T1, typename T2> T1 operator() (T1 x, T2 y) const
  { throw Exception (string("atan2 not available for type ")+typeid(T1).name());  }
  static string Name() { return "atan2"; }
};

struct GenericPow {
  double operator() (double x, double y) const { return pow(x,y); }
  Complex operator() (Complex x, Complex y) const { return pow(x,y); }
  template <typename T1, typename T2> T1 operator() (T1 x, T2 y) const
  {
    throw Exception (string("pow not available for type ")+typeid(T1).name());
  }    
  static string Name() { return "pow"; }
};





  template <int D>
  class NormalVectorCF : public CoefficientFunctionNoDerivative
  {
  public:
    NormalVectorCF () : CoefficientFunctionNoDerivative(D,false) { ; }
    // virtual int Dimension() const { return D; }

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
      const TPMappedIntegrationRule * tpir = dynamic_cast<const TPMappedIntegrationRule *>(&ir);
       if(!tpir)
       {
         if (ir[0].Dim() != D)
           throw Exception("illegal dim of normal vector");
         FlatMatrixFixWidth<D> resD(res);
         for (int i = 0; i < ir.Size(); i++)
           resD.Row(i) = static_cast<const DimMappedIntegrationPoint<D>&>(ir[i]).GetNV();
       }
       else
       {
         int facet = tpir->GetFacet();
         auto & mir = *tpir->GetIRs()[facet];
         int dim = mir[0].Dim();
         int ii = 0;
         res = 0.0;
         if(facet == 0)
         {
           if(dim == 1)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)
                 res.Row(ii++).Range(0,dim) = static_cast<const DimMappedIntegrationPoint<1>&>(mir[i]).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 2)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(0,dim) = static_cast<const DimMappedIntegrationPoint<2>&>(mir[i]).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 3)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(0,dim) = static_cast<const DimMappedIntegrationPoint<3>&>(mir[i]).GetNV();//res1.Row(i).Range(0,dim);
         }
         else
         {
           if(dim == 1)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)
                 res.Row(ii++).Range(D-dim,D) = static_cast<const DimMappedIntegrationPoint<1>&>(mir[j]).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 2)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(D-dim,D) = static_cast<const DimMappedIntegrationPoint<2>&>(mir[j]).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 3)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(D-dim,D) = static_cast<const DimMappedIntegrationPoint<3>&>(mir[j]).GetNV();//res1.Row(i).Range(0,dim);
         }
      }
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, FlatMatrix<Complex> res) const 
    {
      if (ir[0].Dim() != D)
	throw Exception("illegal dim of normal vector");
      for (int i = 0; i < ir.Size(); i++)
	res.Row(i) = static_cast<const DimMappedIntegrationPoint<D>&>(ir[i]).GetNV();
    }
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
        string miptype;
        if(code.is_simd)
          miptype = "SIMD<DimMappedIntegrationPoint<"+ToLiteral(D)+">>*";
        else
          miptype = "DimMappedIntegrationPoint<"+ToLiteral(D)+">*";
        auto nv_expr = CodeExpr("static_cast<const "+miptype+">(&ip)->GetNV()");
        auto nv = Var("tmp", index);
        code.body += nv.Assign(nv_expr);
        for( int i : Range(D))
          code.body += Var(index,i).Assign(nv(i));
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
    {
      for (size_t i = 0; i < ir.Size(); i++)
        for (size_t j = 0; j < D; j++)
          values(j,i) = static_cast<const SIMD<DimMappedIntegrationPoint<D>>&>(ir[i]).GetNV()(j).Data();
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    {
      Evaluate (ir, values);
    }
    
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir, 
                                AFlatMatrix<double> values, AFlatMatrix<double> deriv) const
    {
      Evaluate (ir, values);
      deriv = 0.0;
    }
    
    virtual void EvaluateDeriv (const SIMD_BaseMappedIntegrationRule & ir,
                                FlatArray<AFlatMatrix<>*> input,
                                FlatArray<AFlatMatrix<>*> dinput,
                                AFlatMatrix<> result,
                                AFlatMatrix<> deriv) const
    {
      Evaluate (ir, result);
      deriv = 0.0;
    }
    
  };

  template <int D>
  class TangentialVectorCF : public CoefficientFunctionNoDerivative
  {
  public:
    TangentialVectorCF () : CoefficientFunctionNoDerivative(D,false) { ; }
    // virtual int Dimension() const { return D; }

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
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
        string miptype;
        if(code.is_simd)
          miptype = "SIMD<DimMappedIntegrationPoint<"+ToLiteral(D)+">>*";
        else
          miptype = "DimMappedIntegrationPoint<"+ToLiteral(D)+">*";
        auto tv_expr = CodeExpr("static_cast<const "+miptype+">(&ip)->GetTV()");
        auto tv = Var("tmp", index);
        code.body += tv.Assign(tv_expr);
        for( int i : Range(D))
          code.body += Var(index,i).Assign(tv(i));
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
    {
      for (size_t i = 0; i < ir.Size(); i++)
        for (size_t j = 0; j < D; j++)
          values(j,i) = static_cast<const SIMD<DimMappedIntegrationPoint<D>>&>(ir[i]).GetTV()(j).Data();
    }
  };




void ExportCoefficientFunction(py::module &m)
{

  py::class_<CoefficientFunction, shared_ptr<CoefficientFunction>>
    (m, "CoefficientFunction",
R"raw(A CoefficientFunction (CF) is some function defined on a mesh.
Examples are coordinates x, y, z, domain-wise constants, solution-fields, ...
CFs can be combined by mathematical operations (+,-,sin(), ...) to form new CFs
Parameters:

val : can be one of the following:

  scalar (float or complex):
    Creates a constant CoefficientFunction with value val

  tuple of scalars or CoefficientFunctions:
    Creates a vector or matrix valued CoefficientFunction, use dims=(h,w)
    for matrix valued CF
  list of scalars or CoefficientFunctions:
    Creates a domain-wise CF, use with generator expressions and mesh.GetMaterials()
    and mesh.GetBoundaries()
)raw")
    .def(py::init([] (py::object val, py::object dims)
        {
          shared_ptr<CoefficientFunction> coef;
          
          py::extract<shared_ptr<CF>> ecf(val);
          if (ecf.check())
            coef = UnaryOpCF (ecf(), [](auto x) { return x; }, " ");
          else
            coef = MakeCoefficient(val);
          if(dims)
            {
              try {
                Array<int> cdims = makeCArray<int> (dims);
                coef->SetDimensions(cdims);
              }
              catch (py::type_error){ }
            }
          return coef;
        }),
        py::arg("coef"),py::arg("dims")=DummyArgument(),
         "Construct a CoefficientFunction from either one of\n"
         "  a scalar (float or complex)\n"
         "  a tuple of scalars and or CFs to define a vector-valued CF\n"
         "     use dims=(h,w) to define matrix-valued CF\n"
         "  a list of scalars and or CFs to define a domain-wise CF"
        )
    .def("__str__",  [](CF& self) { return ToString<>(self);})

    .def("__call__", [] (CF& self, BaseMappedIntegrationPoint & mip) -> py::object
	  {
	    if (!self.IsComplex())
	      {
                if (self.Dimension() == 1)
                  return py::cast(self.Evaluate(mip));
                Vector<> vec(self.Dimension());
                self.Evaluate (mip, vec);
                py::tuple res(self.Dimension());
                for (auto i : Range(vec))
                  res[i] = py::cast(vec[i]);
                return res;
	      }
	    else
	      {
                Vector<Complex> vec(self.Dimension());
                self.Evaluate (mip, vec);
                if (vec.Size()==1) return py::cast(vec(0));
                py::tuple res(self.Dimension());
                for (auto i : Range(vec))
                  res[i] = py::cast(vec[i]);
                return res;
	      }
	  },
         py::arg("mip"),
         "evaluate CF at a mapped integrationpoint mip. mip can be generated by calling mesh(x,y,z)")
    .def_property_readonly("dim",
         [] (CF& self) { return self.Dimension(); } ,
                  "number of components of CF")

    /*
    .def_property_readonly("dims",
         [] (PyCF self) { return self->Dimensions(); } ,
                  "shape of CF:  (dim) for vector, (h,w) for matrix")    
    */
    .def_property("dims",
                  [] (shared_ptr<CF> self) { return Array<int>(self->Dimensions()); } ,
                  [] (shared_ptr<CF> self, py::tuple tup) { self->SetDimensions(makeCArray<int>(tup)); } ,
                  "shape of CF:  (dim) for vector, (h,w) for matrix")    
    
    .def("__getitem__",  [](shared_ptr<CF> self, int comp)
                                         {
                                           if (comp < 0 || comp >= self->Dimension())
                                             throw py::index_error();
                                           return MakeComponentCoefficientFunction (self, comp);
                                         },
         py::arg("comp"),         
         "returns component comp of vectorial CF")
    .def("__getitem__",  [](shared_ptr<CF> self, py::tuple comps)
                                         {
                                           if (py::len(comps) != 2)
                                             throw py::index_error();
                                           FlatArray<int> dims = self->Dimensions();
                                           if (dims.Size() != 2)
                                             throw py::index_error();
                                           
                                           int c1 = py::extract<int> (comps[0])();
                                           int c2 = py::extract<int> (comps[1])();
                                           if (c1 < 0 || c2 < 0 || c1 >= dims[0] || c2 >= dims[1])
                                             throw py::index_error();

                                           int comp = c1 * dims[1] + c2;
                                           return MakeComponentCoefficientFunction (self, comp);
                                         })

    // coefficient expressions
    .def ("__add__", [] (shared_ptr<CF> c1, shared_ptr<CF> c2) { return c1+c2; } )
    .def ("__add__", [] (shared_ptr<CF> coef, double val)
           {
             return coef + make_shared<ConstantCoefficientFunction>(val);
           })
    .def ("__radd__", [] (shared_ptr<CF> coef, double val)
           { return coef + make_shared<ConstantCoefficientFunction>(val); })

    .def ("__sub__", [] (shared_ptr<CF> c1, shared_ptr<CF> c2)
           { return c1-c2; })

    .def ("__sub__", [] (shared_ptr<CF> coef, double val)
           { return coef - make_shared<ConstantCoefficientFunction>(val); })

    .def ("__rsub__", [] (shared_ptr<CF> coef, double val)
           { return make_shared<ConstantCoefficientFunction>(val) - coef; })

    .def ("__mul__", [] (shared_ptr<CF> c1, shared_ptr<CF> c2)
           {
             return c1*c2;
           } )

    .def ("__pow__", [] (shared_ptr<CF> c1, int p)
           {
             shared_ptr<CF> one = make_shared<ConstantCoefficientFunction>(1.0);
             if(p==0) return one;

             unsigned n = abs(p);
             shared_ptr<CF> square = c1;
             shared_ptr<CF> res;

             // exponentiation by squaring
             while(n)
             {
               if(n%2)
               {
                 // check if res was not yet assigned any value
                 res = res ? res*square : square;
               }
               square = square*square;
               n /= 2;
             }

             if(p<0)
               return one/res;
             else
               return res;
           } )

    .def ("__pow__", [] (shared_ptr<CF> c1, shared_ptr<CF> c2)
           {
             GenericPow func;
             return BinaryOpCF(c1, c2, func,
                               [](bool a, bool b) { return a||b; }, 'X' /* FUNC::Name() */);
           } )

    .def ("__pow__", [] (shared_ptr<CF> c1, double val)
           {
             GenericPow func;
	     auto c2 = make_shared<ConstantCoefficientFunction>(val);
             return BinaryOpCF(c1, c2, func,
                               [](bool a, bool b) { return a||b; }, 'X' /* FUNC::Name() */);
           } )  

    .def ("InnerProduct", [] (shared_ptr<CF> c1, shared_ptr<CF> c2)
           { 
             return InnerProduct (c1, c2);
           })
    
    .def("Norm",  [](shared_ptr<CF> x) { return NormCF(x); })
    
    .def ("Other",
          [](shared_ptr<CF> x) { return MakeOtherCoefficientFunction(x); },
          "evaluate on other element, as needed for DG jumps")
    
    // it's using the complex functions anyway ...
    // it seems to take the double-version now
    .def ("__mul__", [] (shared_ptr<CF> coef, double val)
           {
             return val * coef; 
           })
    .def ("__rmul__", [] (shared_ptr<CF> coef, double val)
           { return val * coef; }
           )

    .def ("__mul__", [] (shared_ptr<CF> coef, Complex val)
           {
             if (val.imag() == 0)
               return val.real() * coef;
             else
               return val * coef;
           })
    .def ("__rmul__", [] (shared_ptr<CF> coef, Complex val)
           { 
             if (val.imag() == 0)
               return val.real() * coef;
             else
               return val * coef;
           })

    .def ("__truediv__", [] (shared_ptr<CF> coef, shared_ptr<CF> coef2)
           { return coef/coef2;
           })

    .def ("__truediv__", [] (shared_ptr<CF> coef, double val)
           // { return coef.Get() * make_shared<ConstantCoefficientFunction>(1/val); })
           { return (1/val) * coef; })

    .def ("__truediv__", [] (shared_ptr<CF> coef, Complex val)
           { return (1.0/val) * coef; })

    .def ("__rtruediv__", [] (shared_ptr<CF> coef, double val)
           { return make_shared<ConstantCoefficientFunction>(val) / coef; })
    .def ("__rtruediv__", [] (shared_ptr<CF> coef, Complex val)
           { return make_shared<ConstantCoefficientFunctionC>(val) / coef; })

    .def ("__neg__", [] (shared_ptr<CF> coef)
           { return -1.0 * coef; })

    .def_property_readonly ("trans", [] (shared_ptr<CF> coef)
                    {
                      return TransposeCF(coef);
                    },
                   "transpose of matrix-valued CF")
    .def_property_readonly ("real", [](shared_ptr<CF> coef) { return Real(coef); }, "real part of CF")
    .def_property_readonly ("imag", [](shared_ptr<CF> coef) { return Imag(coef); }, "imaginary part of CF")

    .def ("Compile", [] (shared_ptr<CF> coef, bool realcompile, int maxderiv, bool wait)
           { return Compile (coef, realcompile, maxderiv, wait); },
           py::arg("realcompile")=false,
           py::arg("maxderiv")=2,
           py::arg("wait")=false,
          "compile list of individual steps, experimental improvement for deep trees")
    ;

  ExportStdMathFunction<GenericSin>(m, "sin");
  ExportStdMathFunction<GenericCos>(m, "cos");
  ExportStdMathFunction<GenericTan>(m, "tan");
  ExportStdMathFunction<GenericExp>(m, "exp");
  ExportStdMathFunction<GenericLog>(m, "log");
  ExportStdMathFunction<GenericATan>(m, "atan");
  ExportStdMathFunction<GenericSqrt>(m, "sqrt");
  ExportStdMathFunction<GenericFloor>(m, "floor");
  ExportStdMathFunction<GenericCeil>(m, "ceil");
  ExportStdMathFunction<GenericConj>(m, "Conj");

  ExportStdMathFunction2<GenericATan2>(m, "atan2");
  ExportStdMathFunction2<GenericPow>(m, "pow");

  m.def ("IfPos", [] (shared_ptr<CF> c1, py::object then_obj, py::object else_obj)
            {
              return IfPos(c1,
                           MakeCoefficient(then_obj),
                           MakeCoefficient(else_obj));
            } ,docu_string(R"raw_string(Returns new CoefficientFunction with values then_obj if c1 is positive and else_obj else.

Parameters:

c1 : ngsolve.CoefficientFunction
  Indicator function
+
then_obj : object
  Values of new CF if c1 is positive, object must be implicitly convertible to
  ngsolve.CoefficientFunction. See help(:any:`CoefficientFunction` ) for information.

else_obj : object
  Values of new CF if c1 is not positive, object must be implicitly convertible to
  ngsolve.CoefficientFunction. See help(:any:`CoefficientFunction` ) for information.
)raw_string"))
    ;
  
         typedef shared_ptr<ParameterCoefficientFunction> spParameterCF;
         py::class_<ParameterCoefficientFunction, spParameterCF, CF>
    (m, "Parameter", docu_string(R"raw_string(CoefficientFunction with a modifiable value

Parameters:

val : float
  Parameter value
)raw_string"))
    .def ("__init__",
          [] (ParameterCoefficientFunction *instance, double val)
                            {
                              new (instance) ParameterCoefficientFunction(val);
                            })
         .def ("Set", [] (spParameterCF cf, double val)  { cf->SetValue (val); },
          "modify parameter value")
    .def ("Get", [] (spParameterCF cf)  { return cf->GetValue(); },
          "return parameter value")
    ;

  m.def("CoordCF", 
        [] (int direction)
        { return MakeCoordinateCoefficientFunction(direction); },
        "CoefficientFunction for x, y, z"
        );
  
  class MeshSizeCF : public CoefficientFunctionNoDerivative
  {
  public:
    MeshSizeCF () : CoefficientFunctionNoDerivative(1, false) { ; }
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
    {
      if (ip.IP().FacetNr() != -1) // on a boundary facet of the element
        {
          double det = 1;
          switch (ip.Dim())
            {
            case 1: det = fabs (static_cast<const MappedIntegrationPoint<1,1>&> (ip).GetJacobiDet()); break;
            case 2: det = fabs (static_cast<const MappedIntegrationPoint<2,2>&> (ip).GetJacobiDet()); break;
            case 3: det = fabs (static_cast<const MappedIntegrationPoint<3,3>&> (ip).GetJacobiDet()); break;
            default:
              throw Exception("Illegal dimension in MeshSizeCF");
            }
          return det/ip.GetMeasure();
        }
      
      switch (ip.Dim() - int(ip.VB()))
        {
        case 0: throw Exception ("don't have mesh-size on 0-D boundary");
        case 1: return fabs (static_cast<const ScalMappedIntegrationPoint<>&> (ip).GetJacobiDet());
        case 2: return pow (fabs (static_cast<const ScalMappedIntegrationPoint<>&> (ip).GetJacobiDet()), 1.0/2);
        case 3: default:
          return pow (fabs (static_cast<const ScalMappedIntegrationPoint<>&> (ip).GetJacobiDet()), 1.0/3);
        }
      // return pow(ip.GetMeasure(), 1.0/(ip.Dim());
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const
    {
      if (ir[0].IP().FacetNr() != -1)
        for(size_t i : Range(ir))
          values(i) =  fabs (ir[i].GetJacobiDet()) / ir[i].GetMeasure();
      else
        for(size_t i : Range(ir))
          values(i) =  pow(fabs (ir[i].GetJacobiDet()), 1.0/ir.DimElement()).Data();
    }

    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    {
      Evaluate (ir, values);
    }    

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
      if(code.is_simd)
      {
        string type = "SIMD<double>";
        code.body += Var(index).Declare(type);
        code.body += "if (mir[0].IP().FacetNr() != -1)\n{";
        code.body +=  Var(index).Assign( CodeExpr("fabs (ip.GetJacobiDet()) / ip.GetMeasure()"), false );
        code.body += "}else\n";
        code.body += Var(index).Assign( CodeExpr("pow(fabs(ip.GetJacobiDet()), 1.0/mir.DimElement())"), false);
      }
      else
      {
        code.body += Var(index).Declare( "double" );
        code.body += R"CODE_(
        {
          double tmp_res = 0.0;
          if (ip.IP().FacetNr() != -1)
          {
          double det = 1;
          switch (ip.Dim())
            {
            case 1: det = fabs (static_cast<const MappedIntegrationPoint<1,1>&> (ip).GetJacobiDet()); break;
            case 2: det = fabs (static_cast<const MappedIntegrationPoint<2,2>&> (ip).GetJacobiDet()); break;
            case 3: det = fabs (static_cast<const MappedIntegrationPoint<3,3>&> (ip).GetJacobiDet()); break;
            default:
              throw Exception("Illegal dimension in MeshSizeCF");
            }
          tmp_res = det/ip.GetMeasure();
          }
          else
          {
          switch (ip.Dim()) {
            case 1:  tmp_res =      fabs (static_cast<const MappedIntegrationPoint<1,1>&> (ip).GetJacobiDet()); break;
            case 2:  tmp_res = pow (fabs (static_cast<const MappedIntegrationPoint<2,2>&> (ip).GetJacobiDet()), 1.0/2); break;
            default: tmp_res = pow (fabs (static_cast<const MappedIntegrationPoint<3,3>&> (ip).GetJacobiDet()), 1.0/3);
            }
          }
        )CODE_" + Var(index).S() + " = tmp_res;\n}\n;";
      }
    }
  };


  class SpecialCoefficientFunctions
  {
  public:
    shared_ptr<CF> GetMeshSizeCF ()
    { return make_shared<MeshSizeCF>(); }

    shared_ptr<CF> GetNormalVectorCF (int dim)
    { 
      switch(dim)
	{ 
	case 1:
	  return make_shared<NormalVectorCF<1>>();
	case 2:
	  return make_shared<NormalVectorCF<2>>();
	case 3:
	  return make_shared<NormalVectorCF<3>>();
	case 4:
	  return make_shared<NormalVectorCF<4>>();
	case 5:
	  return make_shared<NormalVectorCF<5>>();
	case 6:
	  return make_shared<NormalVectorCF<6>>();
        default:
          throw Exception (string("Normal-vector not implemented for dimension")
                           +ToString(dim));
	}
    }

    shared_ptr<CF> GetTangentialVectorCF (int dim)
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

  py::class_<SpecialCoefficientFunctions> (m, "SpecialCFCreator")
    .def_property_readonly("mesh_size", 
                  &SpecialCoefficientFunctions::GetMeshSizeCF, "local mesh-size (approximate element diameter) as CF")
    .def("normal", &SpecialCoefficientFunctions::GetNormalVectorCF,
         "depending on contents: normal-vector to geometry or element\n"
         "space-dimension must be provided")
    .def("tangential", &SpecialCoefficientFunctions::GetTangentialVectorCF,
         "depending on contents: tangential-vector to element\n"
         "space-dimension must be provided")
    ;
  static SpecialCoefficientFunctions specialcf;
  
  m.attr("specialcf") = py::cast(&specialcf);

  py::class_<BSpline, shared_ptr<BSpline> > (m, "BSpline",R"raw(BSpline of arbitrary order

Parameters:

order : int
  order of the BSpline

knots : list of float

vals : list of float

)raw")
   .def("__init__",
        [](BSpline *instance, int order, py::list knots, py::list vals)
                           {
                             new (instance) BSpline (order,
                                                 makeCArray<double> (knots),
                                                 makeCArray<double> (vals));
                           },
        "B-Spline of a certain order, provide knot and value vectors")
    .def("__str__", &ToString<BSpline>)
    .def("__call__", &BSpline::Evaluate)
    .def("__call__", [](shared_ptr<BSpline> sp, shared_ptr<CF> coef)
          {
            return UnaryOpCF (coef, GenericBSpline(sp) /* , GenericBSpline(sp) */);
          })
    .def("Integrate", 
         [](const BSpline & sp) { return make_shared<BSpline>(sp.Integrate()); })
    .def("Differentiate", 
         [](const BSpline & sp) { return make_shared<BSpline>(sp.Differentiate()); })
    ;
}




// *************************************** Export FEM ********************************


void NGS_DLL_HEADER ExportNgfem(py::module &m) {

  py::enum_<ELEMENT_TYPE>(m, "ET")
    .value("POINT", ET_POINT)     .value("SEGM", ET_SEGM)
    .value("TRIG", ET_TRIG)       .value("QUAD", ET_QUAD)
    .value("TET", ET_TET)         .value("PRISM", ET_PRISM)
    .value("PYRAMID", ET_PYRAMID) .value("HEX", ET_HEX)
    .export_values()
    ;

  py::enum_<NODE_TYPE>(m, "NODE_TYPE")
    .value("VERTEX", NT_VERTEX)
    .value("EDGE", NT_EDGE)
    .value("FACE", NT_FACE)
    .value("CELL", NT_CELL)
    .value("ELEMENT", NT_ELEMENT)
    .value("FACET", NT_FACET)
    .export_values()
    ;


  py::class_<ElementTopology> (m, "ElementTopology")
    .def(py::init<ELEMENT_TYPE>())
    .def_property_readonly("name", 
                  static_cast<const char*(ElementTopology::*)()> (&ElementTopology::GetElementName))
    .def_property_readonly("vertices", [](ElementTopology & self)
                                              {
                                                py::list verts;
                                                const POINT3D * pts = self.GetVertices();
                                                int dim = self.GetSpaceDim();
                                                for (int i : Range(self.GetNVertices()))
                                                  {
                                                    py::list v;
                                                    for (int j = 0; j < dim; j++)
                                                      v.append(py::cast(pts[i][j]));
                                                    verts.append (py::tuple(v));
                                                  }
                                                return verts;
                                              });
    ;
    
  py::class_<FiniteElement, shared_ptr<FiniteElement>>
    (m, "FiniteElement", "any finite element")
    .def_property_readonly("ndof", &FiniteElement::GetNDof, "number of degrees of freedom of element")    
    .def_property_readonly("order", &FiniteElement::Order, "maximal polynomial order of element")    
    .def_property_readonly("type", &FiniteElement::ElementType, "geometric type of element")    
    .def_property_readonly("dim", &FiniteElement::Dim, "spatial dimension of element")    
    .def_property_readonly("classname", &FiniteElement::ClassName, "name of element family")  
    .def("__str__", &ToString<FiniteElement>)
    // .def("__timing__", [] (FiniteElement & fel) { return py::cast(fel.Timing()); })
    .def("__timing__", &FiniteElement::Timing)
    ;

  py::class_<BaseScalarFiniteElement, shared_ptr<BaseScalarFiniteElement>, 
    FiniteElement>
      (m, "ScalarFE", "a scalar-valued finite element")
    .def("CalcShape",
         [] (const BaseScalarFiniteElement & fe, double x, double y, double z)
          {
            IntegrationPoint ip(x,y,z);
            Vector<> v(fe.GetNDof());
            fe.CalcShape (ip, v);
            return v;
          },
         py::arg("x"),py::arg("y")=0.0,py::arg("z")=0.0)
    .def("CalcShape",
         [] (const BaseScalarFiniteElement & fe, const BaseMappedIntegrationPoint & mip)
          {
            Vector<> v(fe.GetNDof());
            fe.CalcShape (mip.IP(), v);
            return v;
          },
         py::arg("mip"))
    .def("CalcDShape",
         [] (const BaseScalarFiniteElement & fe, const BaseMappedIntegrationPoint & mip)
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
          },
         py::arg("mip"))
    ;



  py::implicitly_convertible 
    <BaseScalarFiniteElement, 
    FiniteElement >(); 


  m.def("H1FE", [](ELEMENT_TYPE et, int order)
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
           },
          "creates an H1 finite element of given geometric shape and polynomial order"
          );

  m.def("L2FE", [](ELEMENT_TYPE et, int order)
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
           }
          );


  py::class_<IntegrationPoint>(m, "IntegrationPoint");

  py::class_<IntegrationRule>(m, "IntegrationRule")
    .def("__init__",
         [](IntegrationRule *instance, ELEMENT_TYPE et, int order)
                           {
                             new (instance) IntegrationRule (et, order);
                           },
          py::arg("element type"), py::arg("order"))
    
    .def("__init__",
         [](IntegrationRule *instance, py::list points, py::list weights)
         {
           IntegrationRule * ir = new (instance) IntegrationRule ();
           for (size_t i = 0; i < len(points); i++)
             {
               py::object pnt = points[i];
               IntegrationPoint ip;
               for (int j = 0; j < len(pnt); j++)
                 ip(j) = py::extract<double> (py::tuple(pnt)[j])();
               ip.SetWeight(py::extract<double> (weights[i])());
               ir -> Append (ip);
             }
         },
         py::arg("points"), py::arg("weights"))
    .def("__str__", &ToString<IntegrationRule>)
    .def("__getitem__", [](IntegrationRule & ir, int nr)
                                        {
                                          if (nr < 0 || nr >= ir.Size())
                                            throw py::index_error();
                                          return ir[nr];
                                        })
    .def("Integrate", [](IntegrationRule & ir, py::object func) -> py::object
          {
            py::object sum;
            bool first = true;
            for (const IntegrationPoint & ip : ir)
              {
                py::object val;
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

                val = val.attr("__mul__")(py::cast((double)ip.Weight()));
                if (first)
                  sum = val;
                else
                  sum = sum.attr("__add__")(val);
                first = false;
              }
            return sum;
          })
    ;

  py::class_<BaseMappedIntegrationPoint>(m, "BaseMappedIntegrationPoint")
    .def("__str__",
         [] (const BaseMappedIntegrationPoint & bmip)
          {
            stringstream str;
            if (bmip.IsComplex())
            {
              str << "p = " << bmip.GetPointComplex() << endl;
              str << "jac = " << bmip.GetJacobianComplex() << endl;
            }
            else 
            {
              str << "p = " << bmip.GetPoint() << endl;
              str << "jac = " << bmip.GetJacobian() << endl;
            }
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
          })
    .def_property_readonly("measure", &BaseMappedIntegrationPoint::GetMeasure)
    .def_property_readonly("point", &BaseMappedIntegrationPoint::GetPoint)
    .def_property_readonly("jacobi", &BaseMappedIntegrationPoint::GetJacobian)
    // .def_property_readonly("trafo", &BaseMappedIntegrationPoint::GetTransformation)
    .def_property_readonly("trafo", &BaseMappedIntegrationPoint::GetTransformation)
    .def_property_readonly("elementid", [](BaseMappedIntegrationPoint & mip)
                                               {
                                                 return mip.GetTransformation().GetElementId();
                                               })
    ;

  py::class_<ElementTransformation, shared_ptr<ElementTransformation>>(m, "ElementTransformation")
    .def(py::init([] (ELEMENT_TYPE et, py::list vertices)
        -> shared_ptr<ElementTransformation>
        {
          int nv = ElementTopology::GetNVertices(et);
          int dim = py::len(vertices[0]);
          Matrix<> pmat(nv,dim);
          for (int i : Range(nv))
            for (int j : Range(dim))
              pmat(i,j) = py::extract<double> (vertices[py::int_(i)][py::int_(j)])();
          switch (Dim(et))
            {
            case 1:
              return make_shared<FE_ElementTransformation<1,1>> (et, pmat);
            case 2:
              return make_shared<FE_ElementTransformation<2,2>> (et, pmat);
            case 3:
              return make_shared<FE_ElementTransformation<3,3>> (et, pmat);
            default:
              throw Exception ("cannot create ElementTransformation");
            }
        }),
        py::arg("et")=ET_TRIG,py::arg("vertices"))
    .def_property_readonly("VB", &ElementTransformation::VB)
    .def_property_readonly("spacedim", &ElementTransformation::SpaceDim)
    .def_property_readonly("elementid", &ElementTransformation::GetElementId)
    .def ("__call__", [] (shared_ptr<ElementTransformation> self, double x, double y, double z)
           {
             
             return &(*self)(IntegrationPoint(x,y,z), global_alloc);
           },
          py::arg("x"), py::arg("y")=0, py::arg("z")=0,
          py::return_value_policy::reference)
    .def ("__call__", [] (shared_ptr<ElementTransformation> self, IntegrationPoint & ip)
           {
             return &(*self)(ip, global_alloc);
           },
          py::arg("ip"),
          py::return_value_policy::reference)
    ;


  py::class_<DifferentialOperator, shared_ptr<DifferentialOperator>>
    (m, "DifferentialOperator")
    ;

  typedef BilinearFormIntegrator BFI;
  auto bfi_class = py::class_<BFI, shared_ptr<BFI>> (m, "BFI");
  bfi_class
    .def(py::init([bfi_class] (const string name, py::object py_coef, int dim, bool imag,
                      string filename, py::kwargs kwargs)
        -> shared_ptr<BilinearFormIntegrator>
        {
          auto flags = CreateFlagsFromKwArgs(bfi_class,kwargs);
          Array<shared_ptr<CoefficientFunction>> coef = MakeCoefficients(py_coef);
          auto bfi = GetIntegrators().CreateBFI (name, dim, coef);

          if (!bfi) cerr << "undefined integrator '" << name
                         << "' in " << dim << " dimension" << endl;

          if (filename.length())
            {
              cout << "set integrator filename: " << filename << endl;
              bfi -> SetFileName (filename);
            }
          bfi -> SetFlags (flags);
          if (imag)
            bfi = make_shared<ComplexBilinearFormIntegrator> (bfi, Complex(0,1));
          bfi_class.attr("__initialize__")(bfi,**kwargs);
          return bfi;
        }),
        py::arg("name")="",
        py::arg("coef"),py::arg("dim")=-1,
        py::arg("imag")=false, py::arg("filename")=""
        )
    .def_static("__flags_doc__", [] ()
         {
           return py::dict
             (
              py::arg("dim") = "int = -1\n"
              "Dimension of integrator. If -1 then dim is set when integrator is\n"
              "added to BilinearForm",
              py::arg("definedon") = "ngsolve.Region\n"
              "Region the integrator is defined on. Regions can be obtained by i.e.\n"
              "mesh.Materials('regexp') or mesh.Boundaries('regexp'). If not set\n"
              "integration is done on all volume elements",
              py::arg("definedonelem") = "ngsolve.BitArray\n"
              "Element wise integrator definition."
              );
         })
    .def_static("__special_treated_flags__", [] ()
         {
           return py::dict
             (
              py::arg("definedonelem") = py::cpp_function([](py::object,Flags*,py::list) { ; }),
              py::arg("definedon") = py::cpp_function ([] (py::object, Flags*, py::list) { ; })
              );
         })
    .def("__initialize__", [] (shared_ptr<BFI> self, py::kwargs kwargs)
         {
           if(kwargs.contains("definedon"))
             {
               auto definedon = kwargs["definedon"];
               auto definedon_list = py::extract<py::list>(definedon);
               if (definedon_list.check())
                 {
                   Array<int> defon = makeCArray<int> (definedon_list());
                   for (int & d : defon) d--;
                   self->SetDefinedOn (defon);
                 }
               else if (py::extract<BitArray> (definedon).check())
                 self->SetDefinedOn (py::extract<BitArray> (definedon)());
               else if (!py::extract<DummyArgument>(definedon).check())
                 throw Exception (string ("cannot handle definedon of type <todo>"));
             }
           if(kwargs.contains("definedonelem"))
               self->SetDefinedOnElements(py::cast<shared_ptr<BitArray>>(kwargs["definedonelem"]));
         })
    .def("__str__",  [](shared_ptr<BFI> self) { return ToString<BilinearFormIntegrator>(*self); } )

    .def("Evaluator",  [](shared_ptr<BFI> self, string name ) { return self->GetEvaluator(name); } )
    // .def("DefinedOn", &Integrator::DefinedOn)
    .def("GetDefinedOn",  [] (shared_ptr<BFI> self) -> const BitArray &{ return self->GetDefinedOn(); } ,
         py::return_value_policy::reference)

    .def("SetDefinedOnElements",  [](shared_ptr<BFI> self, shared_ptr<BitArray> ba )
                                                  { self->SetDefinedOnElements (ba); } )
    .def("SetIntegrationRule", [] (shared_ptr<BFI> self, ELEMENT_TYPE et, IntegrationRule ir)
         {
           self -> SetIntegrationRule(et,ir);
           return self;
         })
    .def("CalcElementMatrix",
         [] (shared_ptr<BFI> self,
             const FiniteElement & fe, const ElementTransformation &trafo,
             int heapsize, bool complex)
                         {
                           while (true)
                             {
                               try
                                 {
                                   LocalHeap lh(heapsize);
                                   if (complex)
                                     {
                                       Matrix<Complex> mat(fe.GetNDof() * self->GetDimension());
                                       self->CalcElementMatrix(fe,trafo,mat,lh);
                                       return py::cast(mat);
                                     }
                                   else
                                     {
                                       Matrix<> mat(fe.GetNDof() * self->GetDimension());
                                       self->CalcElementMatrix (fe, trafo, mat, lh);
                                       return py::cast(mat);
                                     }
                                 }
                               catch (LocalHeapOverflow ex)
                                 {
                                   heapsize *= 10;
                                 }
                             }
                         },
         py::arg("fel"),py::arg("trafo"),py::arg("heapsize")=10000, py::arg("complex") = false)
    ;


  m.def("CompoundBFI", 
          []( shared_ptr<BFI> bfi, int comp )
                            {
                                return make_shared<CompoundBilinearFormIntegrator>(bfi, comp);
                            },
           py::arg("bfi")=NULL, py::arg("comp")=0
      );

  m.def("BlockBFI", 
          []( shared_ptr<BFI> bfi, int dim, int comp )
                            {
                                return make_shared<BlockBilinearFormIntegrator>(bfi, dim, comp);
                            },
           py::arg("bfi")=NULL, py::arg("dim")=2, py::arg("comp")=0
      );

  typedef LinearFormIntegrator LFI;
  py::class_<LFI, shared_ptr<LFI>>
    (m, "LFI")
    .def(py::init([] (string name, int dim,
                      py::object py_coef,
                      py::object definedon, bool imag, const Flags & flags,
                      py::object definedonelem)
                  {
                    Array<shared_ptr<CoefficientFunction>> coef = MakeCoefficients(py_coef);
                    auto lfi = GetIntegrators().CreateLFI (name, dim, coef);

                    if (!lfi) throw Exception(string("undefined integrator '")+name+
                                              "' in "+ToString(dim)+ " dimension having 1 coefficient");

                    if(hasattr(definedon,"Mask"))
                      {
                        auto vb = py::cast<VorB>(definedon.attr("VB")());
                        if(vb != lfi->VB())
                          throw Exception(string("LinearFormIntegrator ") + name + " not defined for " +
                                          (vb==VOL ? "VOL" : (vb==BND ? "BND" : "BBND")));
                        lfi->SetDefinedOn(py::cast<BitArray>(definedon.attr("Mask")()));
                      }
                    if (py::extract<py::list> (definedon).check())
                      {
                        Array<int> defon = makeCArray<int> (definedon);
                        for (int & d : defon) d--;
                        lfi -> SetDefinedOn (defon);
                      }
                    if (! py::extract<DummyArgument> (definedonelem).check())
                      lfi -> SetDefinedOnElements (py::extract<shared_ptr<BitArray>>(definedonelem)());

                    if (imag)
                      lfi = make_shared<ComplexLinearFormIntegrator> (lfi, Complex(0,1));
                    return lfi;
                  }),
         py::arg("name")=NULL,py::arg("dim")=-1,
         py::arg("coef"),py::arg("definedon")=DummyArgument(),
         py::arg("imag")=false, py::arg("flags")=py::dict(),
         py::arg("definedonelements")=DummyArgument())

    .def("__str__",  [](shared_ptr<LFI> self) { return ToString<LinearFormIntegrator>(*self); } )
    
    // .def("GetDefinedOn", &Integrator::GetDefinedOn)
    .def("GetDefinedOn",  [] (shared_ptr<LFI> self) -> const BitArray &{ return self->GetDefinedOn(); } ,
         py::return_value_policy::reference)
    .def("SetDefinedOnElements",  [](shared_ptr<LFI> self, shared_ptr<BitArray> ba )
                                                  { self->SetDefinedOnElements (ba); } )
    .def("SetIntegrationRule", [](shared_ptr<LFI> self, ELEMENT_TYPE et, IntegrationRule ir)
         {
           self->SetIntegrationRule(et,ir);
           return self;
         })

    .def("CalcElementVector", 
        static_cast<void(LinearFormIntegrator::*)(const FiniteElement&, const ElementTransformation&, FlatVector<double>, LocalHeap&)const>(&LinearFormIntegrator::CalcElementVector))
    .def("CalcElementVector",
         [] (shared_ptr<LFI>  self, const FiniteElement & fe, const ElementTransformation& trafo,
             int heapsize, bool complex)
         {
           while (true)
             {
               try
                 {
                   LocalHeap lh(heapsize);
                   if (complex)
                     {
                       Vector<Complex> vec(fe.GetNDof() * self->GetDimension());
                       self->CalcElementVector(fe,trafo,vec,lh);
                       return py::cast(vec);
                     }
                   else
                     {
                       Vector<> vec(fe.GetNDof() * self->GetDimension());
                       self->CalcElementVector (fe, trafo, vec, lh);
                       return py::cast(vec);
                     }
                 }
               catch (LocalHeapOverflow ex)
                 {
                   heapsize *= 10;
                 }
             };
         },
         py::arg("fel"),py::arg("trafo"),py::arg("heapsize")=10000, py::arg("complex")=false)
    ;



  m.def("CompoundLFI", 
          []( shared_ptr<LFI> lfi, int comp )
                            {
                                return shared_ptr<LFI>(make_shared<CompoundLinearFormIntegrator>(lfi, comp));
                            },
           "lfi"_a=NULL, py::arg("comp")=0);

  m.def("BlockLFI", 
          []( shared_ptr<LFI> lfi, int dim, int comp )
                            {
                                return shared_ptr<LFI>(make_shared<BlockLinearFormIntegrator>(lfi, dim, comp));
                            },
           "lfi"_a=NULL, py::arg("dim")=2, py::arg("comp")=0);


  ExportCoefficientFunction (m);

  m.def ("SetPMLParameters", 
           [] (double rad, double alpha)
                           {
                             cout << "set pml parameters, r = " << rad << ", alpha = " << alpha << endl;
                             constant_table_for_FEM = &pmlpar;
                             pmlpar.Set("pml_r", rad);
                             pmlpar.Set("pml_alpha", alpha);
                             SetPMLParameters();
                           },
           py::arg("rad")=1,py::arg("alpha")=1);
    
                           
  m.def("GenerateL2ElementCode", &GenerateL2ElementCode);

}


PYBIND11_MODULE(libngfem, m) {
  m.attr("__name__") = "fem";
  ExportNgfem(m);
}



#endif
