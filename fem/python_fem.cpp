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

// TODO:
// class PythonCFWrap : public PythonCoefficientFunction , public py::wrapper<PythonCoefficientFunction> {
//     static std::mutex m;
//     public:
//         PythonCFWrap () : PythonCoefficientFunction() { ; }
//         double EvaluateXYZ (double x, double y, double z) const {
//             return this->get_override("EvaluateXYZ")(x,y,z);
//         }
// 
//         double Evaluate (const BaseMappedIntegrationPoint & bip) const {
//             double ret = 0;
//             m.lock();
//             try { 
//                 ret = this->get_override("Evaluate")(boost::ref(bip)); 
//             }
//             catch (py::error_already_set const &) {
//                 PyErr_Print();
//             }
//             catch(...) {
//                 cout << "caught Exception in PythonCoefficientFunction::Evaluate" << endl;
//             }
//             m.unlock();
//             return ret;
//         }
// };
// 
// std::mutex PythonCFWrap::m;






/*
shared_ptr<CoefficientFunction> MakeCoefficient (py::object py_coef)
{
  if (py::extract<shared_ptr<CoefficientFunction>>(py_coef).check())
    return py::extract<shared_ptr<CoefficientFunction>>(py_coef)();
  else if (py::extract<double>(py_coef).check())
    return make_shared<ConstantCoefficientFunction> 
      (py::extract<double>(py_coef)());
  else
    {
      py::exec("raise KeyError()\n");
      return nullptr;
    }
}
*/
typedef CoefficientFunction CF;
typedef PyWrapper<CoefficientFunction> PyCF;

PyCF MakeCoefficient (py::object val)
{
  py::extract<PyCF> ecf(val);
  if (ecf.check()) return ecf();

  py::extract<double> ed(val);
  if (ed.check()) 
    return PyCF(make_shared<ConstantCoefficientFunction> (ed()));

  py::extract<Complex> ec(val);
  if (ec.check()) 
    return PyCF(make_shared<ConstantCoefficientFunctionC> (ec()));

  py::extract<py::list> el(val);
  if (el.check())
    {
      Array<shared_ptr<CoefficientFunction>> cflist(py::len(el()));
      for (int i : Range(cflist))
        cflist[i] = MakeCoefficient(el()[i]).Get();
      return PyCF(MakeDomainWiseCoefficientFunction(move(cflist)));
    }

  py::extract<py::tuple> et(val);
  if (et.check())
    {
      Array<shared_ptr<CoefficientFunction>> cflist(py::len(et()));
      for (int i : Range(cflist))
        cflist[i] = MakeCoefficient(et()[i]).Get();
      return PyCF(MakeVectorialCoefficientFunction(move(cflist)));
    }


  throw Exception ("cannot make coefficient");
}

Array<shared_ptr<CoefficientFunction>> MakeCoefficients (py::object py_coef)
{
  Array<shared_ptr<CoefficientFunction>> tmp;
  if (py::extract<py::list>(py_coef).check())
    {
      auto l = py::extract<py::list>(py_coef)();
      for (int i = 0; i < py::len(l); i++)
        tmp += MakeCoefficient(l[i]).Get();
    }
  else if (py::extract<py::tuple>(py_coef).check())
    {
      auto l = py::extract<py::tuple>(py_coef)();
      for (int i = 0; i < py::len(l); i++)
        tmp += MakeCoefficient(l[i]).Get();
    }
  else
    tmp += MakeCoefficient(py_coef).Get();

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
              if (py::extract<PyCF>(x).check())
                {
                  auto coef = py::extract<PyCF>(x)();
                  return py::cast(PyCF(UnaryOpCF(coef.Get(), func, func, FUNC::Name())));
                }
              py::extract<double> ed(x);
              if (ed.check()) return py::cast(func(ed()));
              if (py::extract<Complex> (x).check())
                return py::cast(func(py::extract<Complex> (x)()));
              throw py::type_error ("can't compute math-function");
            });
}

template <typename FUNC>
void ExportStdMathFunction2(py::module &m, string name)
{
  m.def (name.c_str(), 
         [] (py::object x, py::object y) -> py::object
         {
           FUNC func;
           if (py::extract<PyCF>(x).check() || py::extract<PyCF>(y).check())
             {
               auto cx = MakeCoefficient(x);
               auto cy = MakeCoefficient(y);
               return py::cast(PyCF(BinaryOpCF(cx.Get(), cy.Get(), func, func, func, func,
                                               [](bool a, bool b) { return a||b; }, 'X' /* FUNC::Name() */)));
             }
           py::extract<double> dx(x), dy(y);
           if (dx.check() && dy.check()) return py::cast(func(dx(), dy()));
           // py::extract<Complex> cx(x), cy(y);
           // if (cx.check() && cy.check()) return py::cast(func(cx(), cy()));
           throw py::type_error ("can't compute math-function");
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
  SIMD<double> operator() (SIMD<double> x) const { return x; }
  AutoDiff<1> operator() (AutoDiff<1> x) const { throw Exception ("Conj(..) is not complex differentiable"); }
  AutoDiffDiff<1> operator() (AutoDiffDiff<1> x) const { throw Exception ("Conj(..) is not complex differentiable"); }
};

struct GenericATan2 {
  double operator() (double x, double y) const { return atan2(x,y); }
  double operator() (double x, double y, double & dx, double & dy) const { throw Exception ("atan2 deriv not available");  }
  double operator() (double x, double y, double & ddx, double & dxdy, double & ddy ) const
  { throw Exception ("atan2 dderiv not available");  }
  template <typename T1, typename T2> T1 operator() (T1 x, T2 y) const { throw Exception ("atan2 not available");  }
  static string Name() { return "atan2"; }
};

struct GenericPow {
  double operator() (double x, double y) const { return pow(x,y); }
  double operator() (double x, double y, double & dx, double & dy) const { throw Exception ("pow deriv not available");  }
  double operator() (double x, double y, double & ddx, double & dxdy, double & ddy ) const
  { throw Exception ("pow dderiv not available");  }
  template <typename T1, typename T2> T1 operator() (T1 x, T2 y) const { throw Exception ("pow not available");  }
  static string Name() { return "pow"; }
};





  template <int D>
  class NormalVectorCF : public CoefficientFunction
  {
  public:
    NormalVectorCF () : CoefficientFunction(D,false) { ; }
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
         int dim = tpir->GetIRs()[facet]->operator[](0).Dim();
         int ii = 0;
         res = 0.0;
         if(facet == 0)
           if(dim == 1)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)
                 res.Row(ii++).Range(0,dim) = static_cast<const DimMappedIntegrationPoint<1>&>(tpir->GetIRs()[facet]->operator[](i)).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 2)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(0,dim) = static_cast<const DimMappedIntegrationPoint<2>&>(tpir->GetIRs()[facet]->operator[](i)).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 3)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(0,dim) = static_cast<const DimMappedIntegrationPoint<3>&>(tpir->GetIRs()[facet]->operator[](i)).GetNV();//res1.Row(i).Range(0,dim);
         else
           if(dim == 1)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)
                 res.Row(ii++).Range(D-dim,D) = static_cast<const DimMappedIntegrationPoint<1>&>(tpir->GetIRs()[facet]->operator[](j)).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 2)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(D-dim,D) = static_cast<const DimMappedIntegrationPoint<2>&>(tpir->GetIRs()[facet]->operator[](j)).GetNV();//res1.Row(i).Range(0,dim);
           if(dim == 3)
             for(int i=0;i<tpir->GetIRs()[0]->Size();i++)
               for(int j=0;j<tpir->GetIRs()[1]->Size();j++)          
                 res.Row(ii++).Range(D-dim,D) = static_cast<const DimMappedIntegrationPoint<3>&>(tpir->GetIRs()[facet]->operator[](j)).GetNV();//res1.Row(i).Range(0,dim);
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
          miptype = "SIMD<DimMappedIntegrationPoint<"+ToString(D)+">>*";
        else
          miptype = "DimMappedIntegrationPoint<"+ToString(D)+">*";
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
    
  };

  template <int D>
  class TangentialVectorCF : public CoefficientFunction
  {
  public:
    TangentialVectorCF () : CoefficientFunction(D,false) { ; }
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
  };


class CoordCoefficientFunction : public T_CoefficientFunction<CoordCoefficientFunction>
  {
    int dir;
    typedef T_CoefficientFunction<CoordCoefficientFunction> BASE;
  public:
    CoordCoefficientFunction (int adir) : BASE(1, false), dir(adir) { ; }
    using BASE::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
    {
      return ip.GetPoint()(dir);
    }
    virtual void Evaluate(const BaseMappedIntegrationRule & ir,
                          FlatMatrix<> result) const
    {
       const TPMappedIntegrationRule * tpmir = dynamic_cast<const TPMappedIntegrationRule *>(&ir);
       if(!tpmir)
       {
           result.Col(0) = ir.GetPoints().Col(dir);
           return;
       }
       if(dir<=2)
       {
         int ii = 0;
         for(int i=0;i<tpmir->GetIRs()[0]->Size();i++)
           for(int j=0;j<tpmir->GetIRs()[1]->Size();j++)
             result(ii++,0) = tpmir->GetIRs()[0]->GetPoints().Col(dir)(i);
       }
       else
       {
         // int ii = 0;
         // for(int i=0;i<tpmir->GetIRs()[0]->Size();i++)
           // for(int j=0;j<tpmir->GetIRs()[1]->Size();j++)
             // result(ii++,0) = tpmir->GetIRs()[1]->GetPoints().Col(dir-3)(j);
 
             //int ii = 0;
         for(int i=0;i<tpmir->GetIRs()[0]->Size();i++)
           //for(int j=0;j<tpmir->GetIRs()[1]->Size();j++)
             result.Col(0).Rows(i*tpmir->GetIRs()[1]->Size(),(i+1)*tpmir->GetIRs()[1]->Size()) = tpmir->GetIRs()[1]->GetPoints().Col(dir-3);
      }
    }
    virtual void Evaluate(const BaseMappedIntegrationRule & ir,
			  FlatMatrix<Complex> result) const
    {
      result.Col(0) = ir.GetPoints().Col(dir);
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
        auto v = Var(index);
        if(dir==0) code.body += v.Assign(CodeExpr("ip.GetPoint()(0)"));
        if(dir==1) code.body += v.Assign(CodeExpr("ip.GetPoint()(1)"));
        if(dir==2) code.body += v.Assign(CodeExpr("ip.GetPoint()(2)"));
    }

    template <typename T>
    void T_Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<T>> values) const
    {
      auto points = ir.GetPoints();
      size_t nv = ir.Size();
      __assume (nv > 0);
      for (size_t i = 0; i < nv; i++)
        values(i) = SIMD<double> (points.Get(i, dir));
    }
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    {
      Evaluate (ir, values);
    }
    
  };

namespace ngstd {
  template <>
  struct PyWrapperTraits<CoordCoefficientFunction> {
    typedef PyWrapperDerived<CoordCoefficientFunction, ngfem::CoefficientFunction> type;
  };
}


void ExportCoefficientFunction(py::module &m)
{
  py::class_<PyWrapper<CoefficientFunction>>
    (m, "CoefficientFunction",
     "A CoefficientFunction (CF) is some function defined on a mesh.\n"
     "examples are coordinates x, y, z, domain-wise constants, solution-fields, ...\n"
     "CFs can be combined by mathematical operations (+,-,sin(), ...) to form new CFs"
    )

    .def("__init__",
         [](PyCF *instance, py::object val, py::object dims)
                           {
                             auto coef = new (instance) PyCF(MakeCoefficient(val));
                             if (dims)
                               {
                                 try {
                                   Array<int> cdims = makeCArray<int> (dims);
                                   coef->Get()->SetDimensions(cdims);
                                 }
                                 catch (py::type_error){ }
                               }
                           },
          py::arg("coef"),py::arg("dims")=DummyArgument(),
         "Construct a CoefficientFunction from either one of\n"
         "  a scalar (float or complex)\n"
         "  a tuple of scalars and or CFs to define a vector-valued CF\n"
         "     use dims=(h,w) to define matrix-valued CF\n"
         "  a list of scalars and or CFs to define a domain-wise CF"
         )
    .def("__str__", FunctionPointer( [](PyCF self) { return ToString<CoefficientFunction>(*self.Get());}))

    .def("__call__", FunctionPointer
	 ([] (PyCF self_wrapper, BaseMappedIntegrationPoint & mip) -> py::object
	  {
            CF & self = *self_wrapper.Get();
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
	  }),
         py::arg("mip"),
         "evaluate CF at a mapped integrationpoint mip. mip can be generated by calling mesh(x,y,z)")
    .def_property_readonly("dim",
         [] (PyCF self) { return self->Dimension(); } ,
                  "number of components of CF")

    /*
    .def_property_readonly("dims",
         [] (PyCF self) { return self->Dimensions(); } ,
                  "shape of CF:  (dim) for vector, (h,w) for matrix")    
    */
    .def_property("dims",
                  [] (PyCF self) { return self->Dimensions(); } ,
                  [] (PyCF self, py::tuple tup) { self->SetDimensions(makeCArray<int>(tup)); } ,
                  "shape of CF:  (dim) for vector, (h,w) for matrix")    
    
    .def("__getitem__", FunctionPointer( [](PyCF self, int comp) -> PyCF
                                         {
                                           if (comp < 0 || comp >= self->Dimension())
                                             throw py::index_error();
                                           return PyCF(MakeComponentCoefficientFunction (self.Get(), comp));
                                         }),
         py::arg("comp"),         
         "returns component comp of vectorial CF")
    .def("__getitem__", FunctionPointer( [](PyCF self, py::tuple comps) -> PyCF
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
                                           return PyCF(MakeComponentCoefficientFunction (self.Get(), comp));
                                         }))

    // coefficient expressions
    .def ("__add__", FunctionPointer 
          ([] (PyCF c1, PyCF c2) -> PyCF { return c1.Get()+c2.Get(); } ))
    .def ("__add__", FunctionPointer 
          ([] (PyCF coef, double val) -> PyCF
           {
             return coef.Get() + make_shared<ConstantCoefficientFunction>(val);
           }))
    .def ("__radd__", FunctionPointer 
          ([] (PyCF coef, double val) -> PyCF
           { return coef.Get() + make_shared<ConstantCoefficientFunction>(val); }))

    .def ("__sub__", FunctionPointer 
          ([] (PyCF c1, PyCF c2) -> PyCF
           { return c1.Get()-c2.Get(); }))

    .def ("__sub__", FunctionPointer 
          ([] (PyCF coef, double val) -> PyCF
           { return coef.Get() - make_shared<ConstantCoefficientFunction>(val); }))

    .def ("__rsub__", FunctionPointer 
          ([] (PyCF coef, double val) -> PyCF
           { return make_shared<ConstantCoefficientFunction>(val) - coef.Get(); }))

    .def ("__mul__", FunctionPointer 
          ([] (PyCF c1, PyCF c2) -> PyCF
           {
             return c1.Get()*c2.Get();
           } ))

    .def ("InnerProduct", FunctionPointer
          ([] (PyCF c1, PyCF c2) -> PyCF
           { 
             return InnerProduct (c1.Get(), c2.Get());
           }))
          
    .def("Norm", FunctionPointer ( [](PyCF x) -> PyCF { return NormCF(x.Get()); }))

    /*
      // it's using the complex functions anyway ...
    .def ("__mul__", FunctionPointer 
          ([] (PyCF coef, double val) -> PyCF
           { 
             return make_shared<ScaleCoefficientFunction> (val, coef); 
           }))
    .def ("__rmul__", FunctionPointer 
          ([] (PyCF coef, double val) -> PyCF
           { return make_shared<ScaleCoefficientFunction> (val, coef); }))
    */
    .def ("__mul__", FunctionPointer 
          ([] (PyCF coef, Complex val) -> PyCF
           { 
             if (val.imag() == 0)
               return val.real() * coef.Get();
             else
               return val * coef.Get();
           }))
    .def ("__rmul__", FunctionPointer 
          ([] (PyCF coef, Complex val) -> PyCF
           { 
             if (val.imag() == 0)
               return val.real() * coef.Get();
             else
               return val * coef.Get();
           }))

    .def ("__truediv__", FunctionPointer 
          ([] (PyCF coef, PyCF coef2) -> PyCF
           { return coef.Get()/coef2.Get();
           }))

    .def ("__truediv__", FunctionPointer 
          ([] (PyCF coef, double val) -> PyCF
           { return coef.Get() / make_shared<ConstantCoefficientFunction>(val); }))

    .def ("__rtruediv__", FunctionPointer 
          ([] (PyCF coef, double val) -> PyCF
           { return make_shared<ConstantCoefficientFunction>(val) / coef.Get(); }))

    .def ("__neg__", FunctionPointer 
          ([] (PyCF coef) -> PyCF
           { return -1.0 * coef.Get(); }))

    .def_property_readonly ("trans", FunctionPointer
                   ([] (PyCF coef) -> PyCF
                    {
                      return TransposeCF(coef.Get());
                    }),
                   "transpose of matrix-valued CF")

    .def ("Compile", FunctionPointer
          ([] (PyCF coef, bool realcompile) -> PyCF
           { return Compile (coef.Get(), realcompile); }),
           py::arg("realcompile")=false,
          "compile list of individual steps, experimental improvement for deep trees")
    ;

  ExportStdMathFunction<GenericSin>(m, "sin");
  ExportStdMathFunction<GenericCos>(m, "cos");
  ExportStdMathFunction<GenericTan>(m, "tan");
  ExportStdMathFunction<GenericExp>(m, "exp");
  ExportStdMathFunction<GenericLog>(m, "log");
  ExportStdMathFunction<GenericATan>(m, "atan");
  ExportStdMathFunction<GenericSqrt>(m, "sqrt");
  ExportStdMathFunction<GenericConj>(m, "Conj");

  ExportStdMathFunction2<GenericATan2>(m, "atan2");
  ExportStdMathFunction2<GenericPow>(m, "pow");
  
  m.def ("IfPos", FunctionPointer 
           ([] (PyCF c1, py::object then_obj, py::object else_obj) -> PyCF
            {
              return IfPos(c1.Get(),
                           MakeCoefficient(then_obj).Get(),
                           MakeCoefficient(else_obj).Get());
            } ))
    ;
  
  typedef PyWrapper<ConstantCoefficientFunction> PyConstCF;
  py::class_<PyConstCF,PyCF>
    (m, "ConstantCF", "same as CoefficientFunction(c), obsolete")
     .def("__init__", 
          [](PyConstCF *instance, double value)
                             {
                               new (instance) PyConstCF(make_shared<ConstantCoefficientFunction>(value));
                             })
    ;

  // py::implicitly_convertible 
  // <shared_ptr<ConstantCoefficientFunction>, shared_ptr<CoefficientFunction> >(); 

  typedef PyWrapper<ParameterCoefficientFunction> PyParameterCF;
  py::class_<PyParameterCF,PyCF>
    (m, "Parameter", "CoefficientFunction with a modifiable value")
    .def ("__init__",
          [] (PyParameterCF *instance, double val)
                            {
                              new (instance) PyParameterCF(make_shared<ParameterCoefficientFunction>(val));
                            })
    .def ("Set",
          FunctionPointer ([] (PyParameterCF cf, double val)
                           {
                             cf->SetValue (val);
                           }),
          "modify parameter value")
    ;

  // py::implicitly_convertible 
  // <shared_ptr<ParameterCoefficientFunction>, shared_ptr<CoefficientFunction> >(); 

  

  typedef PyWrapper<CoordCoefficientFunction> PyCoordCF;
  py::class_<PyCoordCF,PyCF>
    (m, "CoordCF", "CoefficientFunction for x, y, z")
     .def("__init__",
          [](PyCoordCF *instance, int direction)
                             {
                               new (instance) PyCoordCF(make_shared<CoordCoefficientFunction>(direction));
                             })
    ;

  class MeshSizeCF : public CoefficientFunction
  {
  public:
    MeshSizeCF () : CoefficientFunction(1, false) { ; }
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

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const {
      if(code.is_simd)
        code.body += Var(index).Assign( CodeExpr("pow(ip.GetJacobiDet(), 1.0/mir.DimElement())"));
      else
      {
        code.body += Var(index).Declare( "double" );
        code.body += R"CODE_(
        {
          double tmp_res = 0.0;
          switch (ip.Dim()) {
            case 1:  tmp_res =      fabs (static_cast<const MappedIntegrationPoint<1,1>&> (ip).GetJacobiDet()); break;
            case 2:  tmp_res = pow (fabs (static_cast<const MappedIntegrationPoint<2,2>&> (ip).GetJacobiDet()), 1.0/2); break;
            default: tmp_res = pow (fabs (static_cast<const MappedIntegrationPoint<3,3>&> (ip).GetJacobiDet()), 1.0/3);
        }
        )CODE_" + Var(index).S() + " = tmp_res;\n}\n;";
      }
    }
  };


  class SpecialCoefficientFunctions
  {
  public:
    PyCF GetMeshSizeCF ()
    { return PyCF(make_shared<MeshSizeCF>()); }

    PyCF GetNormalVectorCF (int dim)
    { 
      switch(dim)
	{
	case 1:
	  return PyCF(make_shared<NormalVectorCF<1>>());
	case 2:
	  return PyCF(make_shared<NormalVectorCF<2>>());
	default:
	  return PyCF(make_shared<NormalVectorCF<3>>());
	}
    }

    PyCF GetTangentialVectorCF (int dim)
    { 
      switch(dim)
	{
	case 1:
	  return PyCF(make_shared<TangentialVectorCF<1>>());
	case 2:
	  return PyCF(make_shared<TangentialVectorCF<2>>());
	default:
	  return PyCF(make_shared<TangentialVectorCF<3>>());
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

  py::class_<BSpline, shared_ptr<BSpline> > (m, "BSpline")
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
    .def("__call__", FunctionPointer
         ([](shared_ptr<BSpline> sp, PyCF coef) -> PyCF
          {
            return UnaryOpCF (coef.Get(), GenericBSpline(sp), GenericBSpline(sp));
          }))
    .def("Integrate", 
         FunctionPointer([](const BSpline & sp) { return make_shared<BSpline>(sp.Integrate()); }))
    .def("Differentiate", 
         FunctionPointer([](const BSpline & sp) { return make_shared<BSpline>(sp.Differentiate()); }))
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


  m.def("H1FE", FunctionPointer
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

  m.def("L2FE", FunctionPointer
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


  py::class_<IntegrationPoint>(m, "IntegrationPoint");

  py::class_<IntegrationRule>(m, "IntegrationRule")
    .def("__init__",
         [](IntegrationRule *instance, ELEMENT_TYPE et, int order)
                           {
                             new (instance) IntegrationRule (et, order);
                           },
          py::arg("element type"), py::arg("order"))
    .def("__getitem__", [](IntegrationRule & ir, int nr)
                                        {
                                          if (nr < 0 || nr >= ir.Size())
                                            throw py::index_error();
                                          return ir[nr];
                                        })
    .def("Integrate", FunctionPointer
         ([](IntegrationRule & ir, py::object func) -> py::object
          {
            py::object sum;
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

                  // TODO: fix!!!
//                 val = val * py::cast((double)ip.Weight());
//                 if (sum == DummyArgument())
//                   sum = val;
//                 else
//                   sum = sum+py::cast(val);
              }
            return sum;
          }))
    ;

  py::class_<BaseMappedIntegrationPoint>(m, "BaseMappedIntegrationPoint")
    .def("__str__",
         [] (const BaseMappedIntegrationPoint & bmip)
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

    
  typedef PyWrapper<ElementTransformation> PyElementTransformation;

  py::class_<PyElementTransformation>(m, "ElementTransformation")
    .def("__init__",
         [](PyElementTransformation *instance, ELEMENT_TYPE et, py::list vertices)
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
                                 new (instance) PyElementTransformation( new FE_ElementTransformation<1,1> (et, pmat));
                               case 2:
                                 new (instance) PyElementTransformation( new FE_ElementTransformation<2,2> (et, pmat));
                               case 3:
                                 new (instance) PyElementTransformation( new FE_ElementTransformation<3,3> (et, pmat));
                               default:
                                 throw Exception ("cannot create ElementTransformation");
                               }
                           },
          py::arg("et")=ET_TRIG,py::arg("vertices"))
    .def_property_readonly("VB", &ElementTransformation::VB)
    .def_property_readonly("spacedim", &ElementTransformation::SpaceDim)
    .def_property_readonly("elementid", &ElementTransformation::GetElementId)
    .def ("__call__", FunctionPointer
          ([] (PyElementTransformation & self, double x, double y, double z)
           {
             
             return &(*self)(IntegrationPoint(x,y,z), global_alloc);
           }),
          py::arg("x"), py::arg("y")=0, py::arg("z")=0,
          py::return_value_policy::reference)
    .def ("__call__", FunctionPointer
          ([] (PyElementTransformation & self, IntegrationPoint & ip)
           {
             return &(*self)(ip, global_alloc);
           }),
          py::arg("ip"),
          py::return_value_policy::reference)
    ;


  py::class_<DifferentialOperator, shared_ptr<DifferentialOperator>>
    (m, "DifferentialOperator")
    ;

  
  typedef PyWrapper<BilinearFormIntegrator> PyBFI;
  py::class_<PyBFI>
    (m, "BFI")
    .def("__init__",
         [](PyBFI *instance, string name, int dim, py::object py_coef, py::object definedon, bool imag,
                              string filename, Flags flags)
                           {
                             Array<shared_ptr<CoefficientFunction>> coef = MakeCoefficients(py_coef);
                             auto bfi = GetIntegrators().CreateBFI (name, dim, coef);

                             if (!bfi) cerr << "undefined integrator '" << name 
                                            << "' in " << dim << " dimension" << endl;

                             auto definedon_list = py::extract<py::list>(definedon);
                             if (definedon_list.check())
                               {
                                 Array<int> defon = makeCArray<int> (definedon_list());
                                 for (int & d : defon) d--;
                                 bfi -> SetDefinedOn (defon); 
                               }
                             else if (py::extract<BitArray> (definedon).check())
                               bfi -> SetDefinedOn (py::extract<BitArray> (definedon)());
                             else if (!py::extract<DummyArgument>(definedon).check())
                               throw Exception (string ("cannot handle definedon of type <todo>"));

                             if (filename.length())
                               {
                                 cout << "set integrator filename: " << filename << endl;
                                 bfi -> SetFileName (filename);
                               }
                             bfi -> SetFlags (flags);
                             if (imag)
                               bfi = make_shared<ComplexBilinearFormIntegrator> (bfi, Complex(0,1));
                             new (instance) PyBFI(bfi);
                           },
           py::arg("name")="",py::arg("dim")=-1,py::arg("coef"),
           py::arg("definedon")=DummyArgument(),py::arg("imag")=false, py::arg("filename")="", py::arg("flags") = py::dict())
    
    .def("__str__", FunctionPointer( [](PyBFI self) { return ToString<BilinearFormIntegrator>(*self.Get()); } ))

    .def("Evaluator", FunctionPointer( [](PyBFI self, string name ) { return self->GetEvaluator(name); } ))
    // .def("DefinedOn", &Integrator::DefinedOn)
    .def("GetDefinedOn", FunctionPointer
         ( [] (PyBFI self) -> const BitArray &{ return self->GetDefinedOn(); } ),
         py::return_value_policy::reference)
    

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
         (py::arg("bfi")=NULL, py::arg("fel"),py::arg("trafo"),py::arg("localheap")))
    */

    .def("CalcElementMatrix",
         FunctionPointer([] (PyBFI self,
                             const FiniteElement & fe, const PyElementTransformation & trafo,
                             int heapsize)
                         {
                           Matrix<> mat(fe.GetNDof());
                           while (true)
                             {
                               try
                                 {
                                   LocalHeap lh(heapsize);
                                   self->CalcElementMatrix (fe, *trafo, mat, lh);
                                   return mat;
                                 }
                               catch (LocalHeapOverflow ex)
                                 {
                                   heapsize *= 10;
                                 }
                             };
                         }),
         py::arg("fel"),py::arg("trafo"),py::arg("heapsize")=10000)
    ;


  m.def("CompoundBFI", 
          []( PyBFI bfi, int comp ) -> PyBFI
                            {
                                return make_shared<CompoundBilinearFormIntegrator>(bfi.Get(), comp);
                            },
           py::arg("bfi")=NULL, py::arg("comp")=0
      );

  m.def("BlockBFI", 
          []( PyBFI bfi, int dim, int comp ) -> PyBFI
                            {
                                return make_shared<BlockBilinearFormIntegrator>(bfi.Get(), dim, comp);
                            },
           py::arg("bfi")=NULL, py::arg("dim")=2, py::arg("comp")=0
      );


  py::class_<PyWrapper<CompoundBilinearFormIntegrator>,PyBFI>
      (m, "CompoundBilinearFormIntegrator");

  py::class_<PyWrapper<BlockBilinearFormIntegrator>,PyBFI>
      (m, "BlockBilinearFormIntegrator");


  typedef PyWrapper<LinearFormIntegrator> PyLFI;
  py::class_<PyLFI>
    (m, "LFI")
    .def("__init__", 
         [](PyLFI *instance, string name, int dim, 
                              py::object py_coef,
                              py::object definedon, bool imag, const Flags & flags)
                           {
                             Array<shared_ptr<CoefficientFunction>> coef = MakeCoefficients(py_coef);
                             auto lfi = GetIntegrators().CreateLFI (name, dim, coef);

                             if (!lfi) throw Exception(string("undefined integrator '")+name+
                                                       "' in "+ToString(dim)+ " dimension having 1 coefficient");

                             
                             auto definedon_list = py::extract<py::list>(definedon);
                             if (definedon_list.check())
                               lfi -> SetDefinedOn (makeCArray<int> (definedon_list()));
 
                             if (imag)
                               lfi = make_shared<ComplexLinearFormIntegrator> (lfi, Complex(0,1));

                             new (instance) PyLFI(lfi);
                           },
           py::arg("name")=NULL,py::arg("dim")=-1,
           py::arg("coef"),py::arg("definedon")=DummyArgument(), 
           py::arg("imag")=false, py::arg("flags")=py::dict())

    /*
    .def("__init__", py::make_constructor
         (FunctionPointer ([](string name, int dim, py::list coefs_list,
                              py::object definedon, bool imag, const Flags & flags)
                           {
                             Array<shared_ptr<CoefficientFunction> > coefs = makeCArray<shared_ptr<CoefficientFunction>> (coefs_list);
                             auto lfi = GetIntegrators().CreateLFI (name, dim, coefs);
                             
                             if (py::extract<py::list> (definedon).check())
                               lfi -> SetDefinedOn (makeCArray<int> (definedon));
 
                             if (imag)
                               lfi = make_shared<ComplexLinearFormIntegrator> (lfi, Complex(0,1));

                             // cout << "LFI: Flags = " << flags << endl;
                             if (!lfi) cerr << "undefined integrator '" << name 
                                            << "' in " << dim << " dimension having 1 coefficient"
                                            << endl;
                             return lfi;
                           }),
          (py::arg("name")=NULL,py::arg("dim")=-1,
           py::arg("coef"),py::arg("definedon")=DummyArgument(), 
           py::arg("imag")=false, py::arg("flags")=py::dict()))
        )
    */

    .def("__str__", FunctionPointer( [](PyLFI self) { return ToString<LinearFormIntegrator>(*self.Get()); } ))
    
    // .def("GetDefinedOn", &Integrator::GetDefinedOn)
    .def("GetDefinedOn", FunctionPointer
         ( [] (PyLFI self) -> const BitArray &{ return self->GetDefinedOn(); } ),
         py::return_value_policy::reference)

    .def("CalcElementVector", 
         [] (PyLFI  self, const FiniteElement & fe, const PyElementTransformation & trafo, FlatVector<double> v, LocalHeap &lh) { self->CalcElementVector(fe, *trafo, v, lh); } )
    .def("CalcElementVector",
         [] (PyLFI  self,
                             const FiniteElement & fe, const PyElementTransformation & trafo,
                             int heapsize)
                         {
                           Vector<> vec(fe.GetNDof());
                           while (true)
                             {
                               try
                                 {
                                   LocalHeap lh(heapsize);
                                   self->CalcElementVector (fe, *trafo, vec, lh);
                                   return vec;
                                 }
                               catch (LocalHeapOverflow ex)
                                 {
                                   heapsize *= 10;
                                 }
                             };
                         },
         py::arg("fel"),py::arg("trafo"),py::arg("heapsize")=10000)
     ;



  m.def("CompoundLFI", 
          []( PyLFI lfi, int comp )
                            {
                                return PyLFI(make_shared<CompoundLinearFormIntegrator>(lfi.Get(), comp));
                            },
           "lfi"_a=NULL, py::arg("comp")=0);

  m.def("BlockLFI", 
          []( PyLFI lfi, int dim, int comp )
                            {
                                return PyLFI(make_shared<BlockLinearFormIntegrator>(lfi.Get(), dim, comp));
                            },
           "lfi"_a=NULL, py::arg("dim")=2, py::arg("comp")=0);


  ExportCoefficientFunction (m);


//   typedef PyWrapperDerived<PythonCFWrap, CoefficientFunction> PythonCFWrapWrap;
//   py::class_<PythonCFWrapWrap ,PyCF>(m, "PythonCF")
//     .def("Evaluate",
//         py::pure_virtual( FunctionPointer( [] ( PythonCFWrapWrap self, const BaseMappedIntegrationPoint & ip )
//         {
//           return static_cast<CoefficientFunction &>(*self.Get()).Evaluate(ip);
//         } )))
//     .def("GetCoordinates", FunctionPointer( [] ( PythonCFWrapWrap self, const BaseMappedIntegrationPoint & ip ) 
//         {
//           return self->GetCoordinates(ip);
//         }))
//     ;

  typedef PyWrapper<DomainVariableCoefficientFunction> PyDomVarCF;
  py::class_<PyDomVarCF,PyCF>
    (m, "VariableCF")
    .def("__init__",
         [](PyDomVarCF *instance, string str)
                           {
                             auto ef = make_shared<EvalFunction> (str);
                             new (instance) PyDomVarCF(make_shared<DomainVariableCoefficientFunction>
                               (Array<shared_ptr<EvalFunction>> ({ ef })));
                           })
    ;

  typedef PyWrapper<DomainConstantCoefficientFunction> PyDomConstCF;
  py::class_<PyDomConstCF,PyCF>
    (m, "DomainConstantCF")
    .def("__init__",
         [](PyDomConstCF *instance, py::object coefs)
                           {
                             Array<double> darray (makeCArray<double> (coefs));
                             new (instance) PyDomConstCF(make_shared<DomainConstantCoefficientFunction> (darray));
                           })
    ;

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
    
                           
                           

}


PYBIND11_PLUGIN(libngfem) {
  py::module m("fem", "pybind fem");
  ExportNgfem(m);
  return m.ptr();
}



#endif
