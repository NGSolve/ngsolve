#ifdef NGS_PYTHON
#include <core/python_ngcore.hpp>
#include "../ngstd/python_ngstd.hpp"
#include "python_fem.hpp"
#include "../ngstd/bspline.hpp"

#include "hdivdivfe.hpp"
#include <fem.hpp>
#include <comp.hpp>
#include <mutex>
using namespace ngfem;
using ngfem::ELEMENT_TYPE;

#include "pml.hpp"

#include "voxelcoefficientfunction.hpp"

#include "tpintrule.hpp"
#include "coefficient_stdmath.hpp"

namespace ngfem
{
  extern SymbolTable<double> * constant_table_for_FEM;
  SymbolTable<double> pmlpar;

  shared_ptr<CoefficientFunction> MakeCoefficient (py::object val)
  {
    py::extract<shared_ptr<CoefficientFunction>> ecf(val);
    if (ecf.check()) return ecf();

    // a numpy.complex converts itself to a real, and prints a warning
    // thus we check for it first
    if (string(py::str(val.get_type())) == "<class 'numpy.complex128'>")
      return make_shared<ConstantCoefficientFunctionC> (val.cast<Complex>());

    try
      {
        double value = val.cast<double>();
        if (value)
          return make_shared<ConstantCoefficientFunction> (value);
        else
          return ZeroCF(Array<int>());
      }
    catch(py::cast_error) {}
    try { return make_shared<ConstantCoefficientFunctionC> (val.cast<Complex>()); }
    catch(py::cast_error) {}

    if (py::isinstance<py::list>(val))
      {
        auto el = py::cast<py::list>(val);
        Array<shared_ptr<CoefficientFunction>> cflist(el.size());
        for (int i : Range(cflist.Size()))
          cflist[i] = MakeCoefficient(el[i]);
        return MakeDomainWiseCoefficientFunction(move(cflist));
      }

    if (py::isinstance<py::tuple>(val))
      {
        auto et = py::cast<py::tuple>(val);
        Array<shared_ptr<CoefficientFunction>> cflist(et.size());
        for (int i : Range(cflist.Size()))
          cflist[i] = MakeCoefficient(et[i]);
        return MakeVectorialCoefficientFunction(move(cflist));
      }
    throw std::invalid_argument(string("Cannot make CoefficientFunction from ") + string(py::str(val)) + " of type " + string(py::str(val.get_type())));
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


#include "integratorcf.hpp"


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
  void DoArchive(Archive& ar) { ar & sp; }
};

template <> void
cl_UnaryOpCF<GenericBSpline>::GenerateCode(Code &code, FlatArray<int> inputs, int index) const
{
  // bspline.hpp is not automatically included. so we should include it:
  code.top+= "#include <bspline.hpp>\n";
  
  stringstream s;
  s << "reinterpret_cast<BSpline*>(" << code.AddPointer(lam.sp.get()) << ")";
  code.body += Var(index,0,1).Assign(s.str());
  
  Var(index).Assign( Var(inputs[0],0,c1->Dimensions()).Func(Var(index,0,1).S()+"->operator()"));
}

template <> shared_ptr<CoefficientFunction>
cl_UnaryOpCF<GenericBSpline>::Diff(const CoefficientFunction * var,
                                   shared_ptr<CoefficientFunction> dir) const
{
  if (this == var) return dir;
  return UnaryOpCF(c1, GenericBSpline(lam.sp->Differentiate())) * c1->Diff(var, dir);
}






/*
struct GenericConj {
  template <typename T> T operator() (T x) const { return Conj(x); } // from bla
  static string Name() { return "conj"; }
  SIMD<double> operator() (SIMD<double> x) const { return x; }
  template<typename T>
  AutoDiff<1,T> operator() (AutoDiff<1,T> x) const { throw Exception ("Conj(..) is not complex differentiable"); }
  template<typename T>
  AutoDiffDiff<1,T> operator() (AutoDiffDiff<1,T> x) const { throw Exception ("Conj(..) is not complex differentiable"); }
  void DoArchive(Archive& ar) {}
};
*/

struct GenericATan2 {
  template <typename T1, typename T2> T1 operator() (T1 y, T2 x) const { return atan2(y,x);  }
  SIMD<Complex> operator() (SIMD<Complex> y,SIMD<Complex> x) const { throw Exception("atan not available for SIMD<complex>"); }
  Complex operator() (Complex y,Complex x) const { throw Exception("atan not available for complex"); }
  static string Name() { return "atan2"; }
  void DoArchive(Archive& ar) {}
};

template <> 
shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericATan2>::Diff(const CoefficientFunction * var,
                                  shared_ptr<CoefficientFunction> dir) const
{
  if (var == this) return dir;    
  return (c1->Diff(var,dir)*c2 - c2->Diff(var,dir)*c1) / (c1*c1+c2*c2);
}


struct GenericPow {
  double operator() (double x, double y) const { return pow(x,y); }
  Complex operator() (Complex x, Complex y) const { return pow(x,y); }
  template <typename T1, typename T2> T1 operator() (T1 x, T2 y) const
  {
      return exp (log(x)*y);
  }    
  static string Name() { return "pow"; }
  void DoArchive(Archive& ar) {}
};




template <> shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericPow>::Diff(const CoefficientFunction * var,
                                 shared_ptr<CoefficientFunction> dir) const
{
  if (this == var) return dir;
  // return UnaryOpCF(c1,GenericLog(),"log")*c2->Diff(var, dir)*BinaryOpCF(c1,c2,GenericPow(), "pow") + c2*c1->Diff(var,dir)/c1*BinaryOpCF(c1,c2,GenericPow(), "pow");
  return log(c1)*c2->Diff(var, dir)*BinaryOpCF(c1,c2,GenericPow(), "pow") + c2*c1->Diff(var,dir)/c1*BinaryOpCF(c1,c2,GenericPow(), "pow");
}

template <> shared_ptr<CoefficientFunction>
cl_BinaryOpCF<GenericPow>::DiffJacobi(const CoefficientFunction * var) const
{
  if (this == var) return make_shared<ConstantCoefficientFunction>(1);
  /// auto loga = UnaryOpCF( c1, GenericLog(), "log");
  // auto exp_b_loga = UnaryOpCF(c2*loga, MakeSGenericExp(), "exp");
  auto exp_b_loga = exp(c2*log(c1));
  return exp_b_loga->DiffJacobi(var);
}




  template <int D>
  class ReferenceCoordinateCF : public CoefficientFunctionNoDerivative
  {
  public:
    ReferenceCoordinateCF () : CoefficientFunctionNoDerivative(D,false) { ; }

    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const 
    {
      return 0;
    }

    using CoefficientFunctionNoDerivative::Evaluate;
    virtual void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> res) const 
    {
      if (ip.DimElement() != D)
        throw Exception("illegal dim of tangential vector");
      res.Range(0,D) = ip.IP().Point().Range(0,D);
    }
    /*
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
    */
  };





void ExportCoefficientFunction(py::module &m)
{
  m.def ("IfPos", [] (shared_ptr<CF> c1, py::object then_obj, py::object else_obj)
            {
              return IfPos(c1,
                           MakeCoefficient(then_obj),
                           MakeCoefficient(else_obj));
            }, py::arg("c1"), py::arg("then_obj"), py::arg("else_obj") ,docu_string(R"raw_string(
Returns new CoefficientFunction with values then_obj if c1 is positive and else_obj else.

Parameters:

c1 : ngsolve.CoefficientFunction
  Indicator function

then_obj : object
  Values of new CF if c1 is positive, object must be implicitly convertible to
  ngsolve.CoefficientFunction. See help(:any:`CoefficientFunction` ) for information.

else_obj : object
  Values of new CF if c1 is not positive, object must be implicitly convertible to
  ngsolve.CoefficientFunction. See help(:any:`CoefficientFunction` ) for information.

)raw_string"))
    ;
  
  m.def("CoordCF", 
        [] (int direction)
        { return MakeCoordinateCoefficientFunction(direction); }, py::arg("direction"),
        docu_string(R"raw_string(
CoefficientFunction for x, y, z.

Parameters:

direction : int
  input direction

)raw_string"));
  
  class MeshSizeCF : public CoefficientFunctionNoDerivative
  {
  public:
    MeshSizeCF () : CoefficientFunctionNoDerivative(1, false) { ; }
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      if (ip.IP().FacetNr() != -1) // on a boundary facet of the element
        {
          double det = 1;
          switch (ip.DimSpace())
            {
            case 1: det = fabs (static_cast<const MappedIntegrationPoint<1,1>&> (ip).GetJacobiDet()); break;
            case 2: det = fabs (static_cast<const MappedIntegrationPoint<2,2>&> (ip).GetJacobiDet()); break;
            case 3: det = fabs (static_cast<const MappedIntegrationPoint<3,3>&> (ip).GetJacobiDet()); break;
            default:
              throw Exception("Illegal dimension in MeshSizeCF");
            }
          return det/ip.GetMeasure();
        }
      
      // switch (ip.DimSpace() - int(ip.VB()))
      switch (ip.DimElement())  
        {
        case 0: throw Exception ("don't have mesh-size on 0-D boundary");
        case 1: return fabs (static_cast<const ScalMappedIntegrationPoint<>&> (ip).GetJacobiDet());
        case 2: return pow (fabs (static_cast<const ScalMappedIntegrationPoint<>&> (ip).GetJacobiDet()), 1.0/2);
        case 3: default:
          return pow (fabs (static_cast<const ScalMappedIntegrationPoint<>&> (ip).GetJacobiDet()), 1.0/3);
        }
      // return pow(ip.GetMeasure(), 1.0/(ip.DimSpace());
    }

      using CoefficientFunctionNoDerivative::Evaluate;
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values) const override
    {
      if (ir[0].IP().FacetNr() != -1)
        for(size_t i : Range(ir))
          values(i) =  fabs (ir[i].GetJacobiDet()) / ir[i].GetMeasure();
      else
        for(size_t i : Range(ir))
          values(i) =  pow(fabs (ir[i].GetJacobiDet()), 1.0/ir.DimElement()).Data();
    }

    /*
    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, FlatArray<AFlatMatrix<double>*> input,
                           AFlatMatrix<double> values) const
    {
      Evaluate (ir, values);
    }    
    */
    
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
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
          switch (ip.DimSpace())
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
          switch (ip.DimSpace()) {
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

    
    shared_ptr<CF> GetReferenceCoordinateCF (int dim)
    {
      switch(dim)
	{
	case 1:
	  return make_shared<ReferenceCoordinateCF<1>>();
	case 2:
	  return make_shared<ReferenceCoordinateCF<2>>();          
	default:
	  return make_shared<ReferenceCoordinateCF<3>>();          
	}
    }
    
    shared_ptr<CF> GetNormalVectorCF (int dim)
    {
      return NormalVectorCF(dim);
    }

    shared_ptr<CF> GetTangentialVectorCF (int dim, bool consistent)
    {
      return TangentialVectorCF(dim, consistent);
    }

    shared_ptr<CF> GetJacobianMatrixCF (int dims, int dimr)
    {
      return JacobianMatrixCF(dims,dimr);
    }

    shared_ptr<CF> GetWeingartenCF (int dim)
    {
      return WeingartenCF(dim);
    }

    shared_ptr<CF> GetVertexTangentialVectorsCF (int dim)
    {
      return VertexTangentialVectorsCF(dim);
    }

    shared_ptr<CF> GetEdgeFaceTangentialVectorsCF (int dim)
    {
      return EdgeFaceTangentialVectorsCF(dim);
    }
    
    shared_ptr<CF> GetEdgeCurvatureCF (int dim)
    {
      return EdgeCurvatureCF(dim);
    }

  };

  
  ExportStdMathFunction_<GenericSin>(m, "sin", "Sine of argument in radians");
  ExportStdMathFunction_<GenericCos>(m, "cos", "Cosine of argument in radians");
  ExportStdMathFunction_<GenericTan>(m, "tan", "Tangent of argument in radians");
  ExportStdMathFunction_<GenericSinh>(m, "sinh", "Hyperbolic sine of argument in radians");
  ExportStdMathFunction_<GenericCosh>(m, "cosh", "Hyperbolic cosine of argument in radians");
  ExportStdMathFunction_<GenericExp>(m, "exp", "Exponential function");
  ExportStdMathFunction_<GenericLog>(m, "log", "Logarithm function");
  ExportStdMathFunction_<GenericATan>(m, "atan", "Inverse tangent in radians");
  ExportStdMathFunction_<GenericACos>(m, "acos", "Inverse cosine in radians");
  ExportStdMathFunction_<GenericASin>(m, "asin", "Inverse sine in radians");
  ExportStdMathFunction_<GenericSqrt>(m, "sqrt", "Square root function");
  ExportStdMathFunction_<GenericErf>(m, "erf", "Error function");
  ExportStdMathFunction_<GenericFloor>(m, "floor", "Round to next lower integer");
  ExportStdMathFunction_<GenericCeil>(m, "ceil", "Round to next greater integer");
  // ExportStdMathFunction<GenericConj>(m, "Conj", "Conjugate imaginary part of complex number");
  // ExportStdMathFunction<GenericIdentity>(m, " ", "Passes value through");

  ExportStdMathFunction2<GenericATan2>(m, "atan2", "Four quadrant inverse tangent in radians", "y", "x");
  ExportStdMathFunction2<GenericPow>(m, "pow", "Power function");

  py::class_<SpecialCoefficientFunctions> (m, "SpecialCFCreator")
    .def_property_readonly("mesh_size", 
                  &SpecialCoefficientFunctions::GetMeshSizeCF, "local mesh-size (approximate element diameter) as CF")
    .def("xref", &SpecialCoefficientFunctions::GetReferenceCoordinateCF, py::arg("dim"),
         "element reference-coordinates")
    .def("normal", &SpecialCoefficientFunctions::GetNormalVectorCF, py::arg("dim"),
         "depending on contents: normal-vector to geometry or element\n"
         "space-dimension must be provided")
    .def("tangential", &SpecialCoefficientFunctions::GetTangentialVectorCF, py::arg("dim"), py::arg("consistent")=false,
         "depending on contents: tangential-vector to element\n"
         "space-dimension must be provided")
    .def("JacobianMatrix", [] (SpecialCoefficientFunctions& self, int dim)
	  {
            return self.GetJacobianMatrixCF(dim,dim);
	  },
         py::arg("dim"),
         "Jacobian matrix of transformation to physical element\n"
         "space-dimension must be provided")
    .def("JacobianMatrix", [] (SpecialCoefficientFunctions& self, int dimr, int dims)
	  {
            if (dimr < dims)
              throw Exception("In SpecialCoefficientFunctions::GetJacobianMatrixCF: dimr > dims");
            return self.GetJacobianMatrixCF(dims,dimr);
	  },
         py::arg("dimr"), py::arg("dims"),
         "Jacobian matrix of transformation to physical element\n"
         "space-dimensions dimr >= dims must be provided")
    .def("Weingarten", &SpecialCoefficientFunctions::GetWeingartenCF, py::arg("dim"),
         "Weingarten tensor \n"
         "space-dimension must be provided")
    .def("VertexTangentialVectors", &SpecialCoefficientFunctions::GetVertexTangentialVectorsCF, py::arg("dim"),
         "VertexTangentialVectors \n"
         "space-dimension must be provided")
    .def("EdgeFaceTangentialVectors", &SpecialCoefficientFunctions::GetEdgeFaceTangentialVectorsCF, py::arg("dim"),
         "EdgeFaceTangentialVectors \n"
         "space-dimension must be provided")
    .def("EdgeCurvature", &SpecialCoefficientFunctions::GetEdgeCurvatureCF, py::arg("dim"),
         "EdgeCurvature \n"
         "space-dimension must be provided")
    ;
  static SpecialCoefficientFunctions specialcf;
  
  m.attr("specialcf") = py::cast(&specialcf);

  auto cf_class = py::class_<CoefficientFunction, shared_ptr<CoefficientFunction>>
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
)raw", py::dynamic_attr())
    .def(py::init([] (py::dict data)
                  {
                    CoefficientFunction* cf;
                    py::list lst = py::cast<py::list>(data["childs"]);
                    lst.append(data["data"]);
                    lst.append(data["version_stored"]);
                    lst.append(data["version_needed"]);
                    PyArchive<BinaryInArchive> ar(lst);
                    ar & cf;
                    return cf;
                  }))
    .def(py::init([] (py::object val, optional<py::tuple> dims)
        {
          shared_ptr<CoefficientFunction> coef;
          
          py::extract<shared_ptr<CF>> ecf(val);
          if (ecf.check())
            // coef = UnaryOpCF (ecf(), GenericIdentity{}, " ");
            coef = CreateWrapperCF(ecf()); // to have a different object from input
          else
            coef = MakeCoefficient(val);
          if(dims.has_value())  
            { // can now reshape without additional wrapper
              auto cdims = makeCArray<int> (*dims);
              int dimension = 1;
              for (int d : cdims) dimension *= d;
              if (coef->Dimension() != dimension)
                throw Exception("dims does not fit to dimension of CoefficientFunction");

              coef->SetDimensions(cdims);
            }
          return coef;
        }),
      py::arg("coef"),py::arg("dims")=nullptr,
         "Construct a CoefficientFunction from either one of\n"
         "  a scalar (float or complex)\n"
         "  a tuple of scalars and or CFs to define a vector-valued CF\n"
         "     use dims=(h,w) to define matrix-valued CF\n"
         "  a list of scalars and or CFs to define a domain-wise CF"
        )
    .def("__str__",  [](CF& self) { return ToString<>(self);})

    .def("__call__", [] (CF& self, BaseMappedIntegrationPoint & mip) -> py::object
	  {
            auto Eval = [&] (auto vec) -> py::object
              {
                self.Evaluate (mip, vec);
                if (self.Dimensions().Size() == 0)
                  return py::cast(vec(0));
                py::tuple res(self.Dimension());
                for (auto i : Range(vec))
                  res[i] = py::cast(vec[i]);
                return std::move(res);
              };

	    if (!self.IsComplex())
              return Eval(Vector<> (self.Dimension()));
            else
              return Eval(Vector<Complex> (self.Dimension()));
            
            /*
	    if (!self.IsComplex())
	      {
                Vector<> vec(self.Dimension());
                self.Evaluate (mip, vec);
                if (self.Dimensions().Size() == 0)
                  return py::cast(vec(0));
                py::tuple res(self.Dimension());
                for (auto i : Range(vec))
                  res[i] = py::cast(vec[i]);
                return std::move(res);
	      }
	    else
	      {
                Vector<Complex> vec(self.Dimension());
                self.Evaluate (mip, vec);
                if (vec.Dimensions().Size()==0)
                  return py::cast(vec(0));
                py::tuple res(self.Dimension());
                for (auto i : Range(vec))
                  res[i] = py::cast(vec[i]);
                return std::move(res);
	      }
            */
	  },
         py::arg("mip"),
         "evaluate CF at a mapped integrationpoint mip. mip can be generated by calling mesh(x,y,z)")
    .def("__call__", [](shared_ptr<CF> self, double x, optional<double> y, optional<double> z)
         {
           Vector<> pnt;
           if (y && z)
             pnt = Vector<>{ x, y.value_or(0), z.value_or(0) };
           else if (y)
             pnt = Vector<>{ x, y.value_or(0) };
           else
             pnt = Vector<>{ x };
           return make_shared<ngcomp::PointEvaluationFunctional>(self, pnt);
         }, py::arg("x"), py::arg("y")=nullopt, py::arg("z")=nullopt)
    .def_property_readonly("dim",
         [] (CF& self) { return self.Dimension(); } ,
                  "number of components of CF")
    .def_property("dims",
                  [] (shared_ptr<CF> self) { return Array<int>(self->Dimensions()); } ,
                  [] (shared_ptr<CF> self, py::tuple tup) { self->SetDimensions(makeCArray<int>(tup)); } ,
                  "shape of CF:  (dim) for vector, (h,w) for matrix")
    .def("Reshape", [] (shared_ptr<CF> self, py::tuple tup) { return self->Reshape(makeCArray<int>(tup)); } ,
         "reshape CF:  (dim) for vector, (h,w) for matrix")
    .def("ExtendDimension", [] (shared_ptr<CF> self, py::tuple dims, optional<py::tuple> pos, optional<py::tuple> stride) {
        Array<int> cpos, cstride;
        if (stride) cstride = makeCArray<int>(*stride);
        if (pos) cpos = makeCArray<int>(*pos);
        return MakeExtendDimensionCoefficientFunction(self, makeCArray<int>(dims),
                                                  move(cpos), move(cstride)); } ,
      "Extend shape by 0-padding", py::arg("dims"), py::arg("pos")=nullopt, py::arg("stride")=nullopt)
    .def_property_readonly("is_complex",
                           [] (CF &  self) { return self.IsComplex(); },
                           "is CoefficientFunction complex-valued ?")
    
    .def_property("spacedim",
                  [] (CF & self) { return self.SpaceDim(); },
                  [] (CF & self, int dim) { self.SetSpaceDim(dim); })
    
    .def("__getitem__",  [](shared_ptr<CF> self, int comp)
                                         {
                                           if (comp < 0 || comp >= self->Dimension())
                                             throw py::index_error();
                                           return MakeComponentCoefficientFunction (self, comp);
                                         },
         py::arg("comp"),
         "returns component comp of vectorial CF")
    .def("__getitem__",  [](shared_ptr<CF> self, py::slice inds)
         {
           FlatArray<int> dims = self->Dimensions();
           if (dims.Size() != 1)
             throw py::index_error();

           size_t start, step, n;
           InitSlice( inds, dims[0], start, step, n );
           int first = start;
           Array<int> num = { int(n) };
           Array<int> dist = { int(step) };
           /*
             if (c1 < 0 || c2 < 0 || c1 >= dims[0] || c2 >= dims[1])
             throw py::index_error();
           */
           return MakeSubTensorCoefficientFunction (self, first, move(num), move(dist));
         }, py::arg("components"))

    /*
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
                                         }, py::arg("components"))
    */
//    .def("__getitem__",  [](shared_ptr<CF> self, tuple<int,int> comps)
//         {
//           FlatArray<int> dims = self->Dimensions();
//           if (dims.Size() != 2)
//             throw py::index_error();
//
//           auto [c1,c2] = comps;
//           if (c1 < 0 || c2 < 0 || c1 >= dims[0] || c2 >= dims[1])
//             throw py::index_error();
//
//           int comp = c1 * dims[1] + c2;
//           return MakeComponentCoefficientFunction (self, comp);
//         }, py::arg("components"))
//
//    .def("__getitem__",  [](shared_ptr<CF> self, tuple<py::slice,int> comps)
//         {
//           FlatArray<int> dims = self->Dimensions();
//           if (dims.Size() != 2)
//             throw py::index_error();
//
//           auto [inds,c2] = comps;
//           size_t start, step, n;
//           InitSlice( inds, dims[0], start, step, n );
//           int first = start*dims[1]+c2;
//           Array<int> num = { int(n) };
//           Array<int> dist = { int(step)*dims[1] };
//           /*
//             if (c1 < 0 || c2 < 0 || c1 >= dims[0] || c2 >= dims[1])
//             throw py::index_error();
//           */
//           return MakeSubTensorCoefficientFunction (self, first, move(num), move(dist));
//         }, py::arg("components"))
    
//    .def("__getitem__",  [](shared_ptr<CF> self, tuple<int,py::slice> comps)
//         {
//           FlatArray<int> dims = self->Dimensions();
//           if (dims.Size() != 2)
//             throw py::index_error();
//
//           auto [c1,inds] = comps;
//           size_t start, step, n;
//           InitSlice( inds, dims[1], start, step, n );
//           // cout << "get row " << c1 << ", start=" << start << ", step = " << step << ", n = " << n << endl;
//           int first = start+c1*dims[1];
//           Array<int> num = { int(n) };
//           Array<int> dist = { int(step) };
//           /*
//             if (c1 < 0 || c2 < 0 || c1 >= dims[0] || c2 >= dims[1])
//             throw py::index_error();
//           */
//           return MakeSubTensorCoefficientFunction (self, first, move(num), move(dist));
//         }, py::arg("components"))
    
//    .def("__getitem__",  [](shared_ptr<CF> self, tuple<py::slice,py::slice> comps)
//         {
//           FlatArray<int> dims = self->Dimensions();
//           if (dims.Size() != 2)
//             throw py::index_error();
//
//           auto [inds1,inds2] = comps;
//           size_t start1, step1, n1;
//           InitSlice( inds1, dims[0], start1, step1, n1 );
//           size_t start2, step2, n2;
//           InitSlice( inds2, dims[1], start2, step2, n2 );
//
//           int first = start1*dims[1]+start2;
//           Array<int> num = { int(n1), int(n2) };
//           Array<int> dist = { int(step1)*dims[1], int(step2) };
//           /*
//             if (c1 < 0 || c2 < 0 || c1 >= dims[0] || c2 >= dims[1])
//             throw py::index_error();
//           */
//           return MakeSubTensorCoefficientFunction (self, first, move(num), move(dist));
//         }, py::arg("components"))
    

//        .def("__getitem__",  [](shared_ptr<CF> self, tuple<py::object,py::object,py::object> comps)
//         {
//           FlatArray<int> dims = self->Dimensions();
//           if (dims.Size() != 3)
//             throw py::index_error();
//
//           auto [comp0,comp1,comp2] = comps;
//           Array<py::object> components = {comp0, comp1, comp2};
//           int numslice = 0;
//           int first = 0;
//           Array<int> num;
//           Array<int> dist;
//           int numcontracted = 0;
//           for (auto i : Range(3))
//             {
//               if (py::extract<int> (components[i]).check())
//               {
//                 int c = components[i].cast<int>();
//                 if (c < 0 || c >= dims[i]) throw py::index_error();
//                 for (auto j : Range(dist.Size()))
//                     dist[j] *= dims[i];
//
//                 first *= dims[i];
//                 first += c;
//               }
//             else if (py::extract<shared_ptr<CoefficientFunction>> (components[i]).check())
//               {
//                 shared_ptr<CoefficientFunction> vec = components[i].cast<shared_ptr<CoefficientFunction>>();
//                 if ( vec->Dimensions().Size() != 1 || vec->Dimensions()[0] != self->Dimensions()[i-numcontracted] )
//                   throw py::index_error();
//
//                 self = MakeSingleContractionCoefficientFunction (self, vec, i-numcontracted);
//                 numcontracted++;
//               }
//             else if (py::extract<py::slice> (components[i]).check())
//               {
//                 py::slice slice = components[i].cast<py::slice>();
//                 numslice++;
//                 size_t start, step, n;
//                 InitSlice( slice, dims[i], start, step, n );
//                 num.Append( int(n) );
//                 for (auto j : Range(dist.Size()))
//                   dist[j] *= dims[i];
//                 dist.Append( int(step) );
//
//                 first *= dims[i];
//                 first += start;
//               }
//             else
//               throw Exception("Invalid object. Only integers and slices are allowed");
//             }
//
//           if (numslice == 0)
//             return MakeComponentCoefficientFunction (self, first);
//           else
//             return MakeSubTensorCoefficientFunction (self, first, move(num), move(dist));
//
//         }, py::arg("components"))

    
//    .def("__getitem__",  [](shared_ptr<CF> self, tuple<int,int,int,int> comps)
//         {
//           FlatArray<int> dims = self->Dimensions();
//           if (dims.Size() != 4)
//             throw py::index_error();
//
//           auto [c1,c2,c3,c4] = comps;
//           if (c1 < 0 || c2 < 0 || c3 < 0 || c4 < 0 ||
//               c1 >= dims[0] || c2 >= dims[1] || c3 >= dims[2] || c4 >= dims[3])
//             throw py::index_error();
//
//           int comp = ((c1 * dims[1] + c2) * dims[2] + c3) * dims[3] + c4;
//           return MakeComponentCoefficientFunction (self, comp);
//         }, py::arg("components"))

    

//    .def("__getitem__",  [](shared_ptr<CF> self, tuple<py::slice,py::slice,py::slice,py::slice> comps)
//         {
//           FlatArray<int> dims = self->Dimensions();
//           if (dims.Size() != 4)
//             throw py::index_error();
//
//           auto [inds1,inds2,inds3,inds4] = comps;
//           size_t start1, step1, n1;
//           InitSlice( inds1, dims[0], start1, step1, n1 );
//           size_t start2, step2, n2;
//           InitSlice( inds2, dims[1], start2, step2, n2 );
//           size_t start3, step3, n3;
//           InitSlice( inds3, dims[2], start3, step3, n3 );
//           size_t start4, step4, n4;
//           InitSlice( inds4, dims[3], start4, step4, n4 );
//
//           int first = ((start1*dims[1]+start2)*dims[2]+start3)*dims[3] + start4;
//           Array<int> num = { int(n1), int(n2), int(n3), int(n4) };
//           Array<int> dist = { int(step1)*dims[1]*dims[2]*dims[3], int(step2)*dims[2]*dims[3], int(step3)*dims[3], int(step4) };
//
//           return MakeSubTensorCoefficientFunction (self, first, move(num), move(dist));
//         }, py::arg("components"))

    
//    .def("__getitem__",  [](shared_ptr<CF> self, tuple<shared_ptr<CF>, shared_ptr<CF>> vectors)
//         {
//           FlatArray<int> dims = self->Dimensions();
//           if (dims.Size() != 2)
//             throw py::index_error();
//
//           auto [v0,v1] = vectors;
//           return MakeVectorContractionCoefficientFunction (self, Array{v0,v1});
//         })
//
//    .def("__getitem__",  [](shared_ptr<CF> self, tuple<shared_ptr<CF>, shared_ptr<CF>, shared_ptr<CF>> vectors)
//         {
//           FlatArray<int> dims = self->Dimensions();
//           if (dims.Size() != 3)
//             throw py::index_error();
//
//           auto [v0,v1,v2] = vectors;
//           return MakeVectorContractionCoefficientFunction (self, Array{v0,v1,v2});
//         })

    .def("__getitem__",  [](shared_ptr<CoefficientFunction> self, py::tuple comps)
         {
           FlatArray<int> dims = self->Dimensions();

           // process ellipses
           size_t ellipse_count = 0;
           size_t ellipse_pos = 0;
           for (auto i : Range(comps.size()))
             if (py::extract<py::ellipsis>(comps[i]).check())
               {
                 ellipse_count++;
                 ellipse_pos = i;
               }

           if (ellipse_count > 1)
              throw Exception(ToString(ellipse_count) + " ellipses detected, but only one is allowed.");
           else if (ellipse_count == 1)
             {
               py::list new_comps{};
               for (auto i : Range(comps.size()))
                 {
                   if (i == ellipse_pos)
                       for (auto j : Range(dims.Size() - comps.size() + 1))
                           new_comps.append(py::slice(0, dims[ellipse_pos + j], 1));
                   else
                       new_comps.append(comps[i]);
                 }
               comps = py::tuple(new_comps);
             }
//           cout << "comps: " << comps << endl;

           if (comps.size() != dims.Size())
             throw Exception("Too few indices or slices. Maybe use an ellipse '...'");

           // detect full contractions beforehand
           Array<shared_ptr<CoefficientFunction>> vecs;
           vecs.SetAllocSize(comps.size());
           for (auto i : Range(comps.size()))
               if (!py::extract<int>(comps[i]).check() && // prevent conversion from int to Constant CF
                   py::extract<shared_ptr<CoefficientFunction>>(comps[i]).check())
                 {
                   shared_ptr<CoefficientFunction> vec = comps[i].cast<shared_ptr<CoefficientFunction>>();
                   if (vec->Dimensions().Size() != 1 ||
                       vec->Dimensions()[0] != self->Dimensions()[i])
                       throw py::index_error();
                   vecs.Append(vec);
                 }
           if (vecs.Size() == dims.Size())
               return MakeVectorContractionCoefficientFunction (self, move(vecs));

           int numslice = 0;
           int first = 0;
           Array<int> num;
           Array<int> dist;
           int numcontracted = 0;
           for (auto i : Range(comps.size()))
             {
               if (py::extract<int>(comps[i]).check())
                 {
                   int c = comps[i].cast<int>();
                   if (c < 0 || c >= dims[i])
                       throw py::index_error();
                   for (auto j : Range(dist.Size()))
                       dist[j] *= dims[i];

                   first *= dims[i];
                   first += c;
                 }
               else if (py::extract<shared_ptr<CoefficientFunction>>(comps[i]).check())
                 {
                   shared_ptr<CoefficientFunction> vec =
                           comps[i].cast<shared_ptr<CoefficientFunction>>();
                   if (vec->Dimensions().Size() != 1 ||
                       vec->Dimensions()[0] != self->Dimensions()[i - numcontracted])
                       throw py::index_error();

                   self = MakeSingleContractionCoefficientFunction(self, vec, i - numcontracted);
                   numcontracted++;
                 }
               else if (py::extract<py::slice>(comps[i]).check())
                 {
                   py::slice slice = comps[i].cast<py::slice>();
                   numslice++;
                   size_t start, step, n;
                   InitSlice(slice, dims[i], start, step, n);
                   num.Append(int(n));
                   for (auto j : Range(dist.Size()))
                       dist[j] *= dims[i];
                   dist.Append(int(step));

                   first *= dims[i];
                   first += start;
                 }
               else
                 throw Exception("Invalid object. Only integers, slices and (single) ellipses are allowed");
             }

           if (numslice == 0)
             return MakeComponentCoefficientFunction(self, first);
           else
             return MakeSubTensorCoefficientFunction(self, first, move(num), move(dist));
         })


    // coefficient expressions
    .def ("__add__", [] (shared_ptr<CF> c1, shared_ptr<CF> c2) { return c1+c2; }, py::arg("cf") )
    .def ("__add__", [] (shared_ptr<CF> coef, double val)
           {return coef + make_shared<ConstantCoefficientFunction>(val); }, py::arg("value"))
    .def ("__add__", [] (shared_ptr<CF> coef, Complex val)
           {return coef + make_shared<ConstantCoefficientFunctionC> (val); }, py::arg("value"))
    .def ("__radd__", [] (shared_ptr<CF> coef, double val)
          { return coef + make_shared<ConstantCoefficientFunction>(val); }, py::arg("value"))
    .def ("__radd__", [] (shared_ptr<CF> coef, Complex val)
          { return coef + make_shared<ConstantCoefficientFunctionC> (val); }, py::arg("value"))

    .def ("__sub__", [] (shared_ptr<CF> c1, shared_ptr<CF> c2)
          { return c1-c2; }, py::arg("cf"))

    .def ("__sub__", [] (shared_ptr<CF> coef, double val)
          { return coef - make_shared<ConstantCoefficientFunction>(val); }, py::arg("value"))

    .def ("__rsub__", [] (shared_ptr<CF> coef, double val)
          { return make_shared<ConstantCoefficientFunction>(val) - coef; }, py::arg("value"))

    .def ("__mul__", [] (shared_ptr<CF> c1, shared_ptr<CF> c2)
           {
             return c1*c2;
           }, py::arg("cf") )

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
           }, py::arg("exponent") )
    .def ("__pow__", [m] (shared_ptr<CF> c1, double c2) { return m.attr("pow")(c1, c2); })
    .def ("__rpow__", [m] (shared_ptr<CF> c1, double c2) { return m.attr("pow")(c2, c1); })
    .def ("__pow__", [m] (shared_ptr<CF> c1, shared_ptr<CF> c2) { return m.attr("pow")(c1, c2); })

    .def ("InnerProduct", [] (shared_ptr<CF> c1, shared_ptr<CF> c2)
           { 
             return InnerProduct (c1, c2);
           }, py::arg("cf"), docu_string(R"raw_string( 
Returns InnerProduct with another CoefficientFunction.

Parameters:

cf : ngsolve.CoefficientFunction
  input CoefficientFunction

 )raw_string"))
    
    .def("Norm",  NormCF, "Returns Norm of the CF")

    .def("Eig", EigCF, "Returns eigenvectors and eigenvalues of matrix-valued CF")
    
    .def ("Other", MakeOtherCoefficientFunction,
          "Evaluate on other element, as needed for DG jumps")
    .def ("Operator", [](shared_ptr<CF> coef, string name) {
        return coef->Operator(name);
      })
    .def ("Derive",
          // &CoefficientFunction::Diff,
          [] (shared_ptr<CF> coef, shared_ptr<CF> var, shared_ptr<CF> dir)
          {
            cout << "warning: Derive is deprecated, use Diff instead" << endl;
            return coef->Diff(var.get(), dir);
          },
          "depricated: use 'Diff' instead", 
          py::arg("variable"), py::arg("direction")=1.0)
    .def ("Diff", // &CoefficientFunction::Diff,
          [] (shared_ptr<CF> coef, shared_ptr<CF> var, shared_ptr<CF> dir)
          {
            if (dir)
              return coef->Diff(var.get(), dir);
            else
              return coef->DiffJacobi(var.get());
            /*
            if (var->Dimension() == 1)
              return coef->Diff(var.get(), make_shared<ConstantCoefficientFunction>(1));
            else
              {
                if (coef->Dimension() != 1)
                  throw Exception("cannot differentiate vectorial CFs by vectrial CFs");
                int dim = var->Dimension();
                Array<shared_ptr<CoefficientFunction>> ddi(dim), ei(dim);
                auto zero = ZeroCF(Array<int>());//make_shared<ConstantCoefficientFunction>(0);
                auto one =  make_shared<ConstantCoefficientFunction>(1);
                for (int i = 0; i < dim; i++)
                  {
                    ei = zero;
                    ei[i] = one;
                    auto vec = MakeVectorialCoefficientFunction (Array<shared_ptr<CoefficientFunction>>(ei));
                    vec->SetDimensions(var->Dimensions());
                    ddi[i] = coef->Diff(var.get(), vec);
                  }
                auto dvec = MakeVectorialCoefficientFunction (move(ddi));
                dvec->SetDimensions(var->Dimensions());
                return dvec;
              }
            */
          },
          "Compute directional derivative with respect to variable",
          py::arg("variable"), py::arg("direction")=nullptr)

    .def ("DiffShape", [] (shared_ptr<CF> coef, shared_ptr<CF> dir, std::vector<shared_ptr<CoefficientFunction>> Eulerian)
          {
            DiffShapeCF shape;
            for (auto gf : Eulerian)
              shape.Eulerian_gridfunctions.Append(gf.get());
            return coef->Diff (&shape, dir);
            
            // return coef->Diff (shape.get(), dir);
          },
          "Compute shape derivative in direction", 
          py::arg("direction")=1.0,  py::arg("Eulerian")=std::vector<shared_ptr<CoefficientFunction>> ())
    
    // it's using the complex functions anyway ...
    // it seems to take the double-version now
    .def ("__mul__", [] (shared_ptr<CF> coef, double val)
           {
             return val * coef; 
           }, py::arg("value"))
    .def ("__rmul__", [] (shared_ptr<CF> coef, double val)
          { return val * coef; }, py::arg("value")
           )

    .def ("__mul__", [] (shared_ptr<CF> coef, Complex val)
           {
             if (val.imag() == 0)
               return val.real() * coef;
             else
               return val * coef;
           }, py::arg("value"))

    .def("__mul__", [](shared_ptr<CoefficientFunction> cf, DifferentialSymbol & dx)
         {
           // return make_shared<SumOfIntegrals>(make_shared<Integral> (cf, dx));
           return make_shared<SumOfIntegrals> (dx.MakeIntegral (cf));  // overloaded in ngsxfem
         })
    
    .def ("__rmul__", [] (shared_ptr<CF> coef, Complex val)
           { 
             if (val.imag() == 0)
               return val.real() * coef;
             else
               return val * coef;
           }, py::arg("value"))

    .def ("__truediv__", [] (shared_ptr<CF> coef, shared_ptr<CF> coef2)
           { return coef/coef2;
           }, py::arg("cf"))

    .def ("__truediv__", [] (shared_ptr<CF> coef, double val)
           // { return coef.Get() * make_shared<ConstantCoefficientFunction>(1/val); })
          { return (1/val) * coef; }, py::arg("value"))

    .def ("__truediv__", [] (shared_ptr<CF> coef, Complex val)
          { return (1.0/val) * coef; }, py::arg("value"))

    .def ("__rtruediv__", [] (shared_ptr<CF> coef, double val)
          { return make_shared<ConstantCoefficientFunction>(val) / coef; }, py::arg("value"))
    .def ("__rtruediv__", [] (shared_ptr<CF> coef, Complex val)
          { return make_shared<ConstantCoefficientFunctionC>(val) / coef; }, py::arg("value"))

    .def ("__neg__", [] (shared_ptr<CF> coef)
           { return -1.0 * coef; })

    .def_property_readonly ("trans", [] (shared_ptr<CF> coef)
                    {
                      return TransposeCF(coef);
                    },
                   "transpose of matrix-valued CF")
    .def ("TensorTranspose", [](shared_ptr<CF> coef, py::tuple ordering)
          {
            return MakeTensorTransposeCoefficientFunction (coef, makeCArray<int>(ordering));
          })
    .def_property_readonly ("real", [](shared_ptr<CF> coef) { return Real(coef); }, "real part of CF")
    .def_property_readonly ("imag", [](shared_ptr<CF> coef) { return Imag(coef); }, "imaginary part of CF")

    .def ("Freeze", [] (shared_ptr<CF> coef)
          { return Freeze(coef); },
          "don't differentiate this expression")

    .def ("Compile", [] (shared_ptr<CF> coef, bool realcompile, int maxderiv, bool wait, bool keep_files)
           { return Compile (coef, realcompile, maxderiv, wait, keep_files); },
           py::arg("realcompile")=false,
           py::arg("maxderiv")=2,
          py::arg("wait")=false,
          py::arg("keep_files")=false,
          py::call_guard<py::gil_scoped_release>(), docu_string(R"raw_string(
Compile list of individual steps, experimental improvement for deep trees

Parameters:

realcompile : bool
  True -> Compile to C++ code

maxderiv : int
  input maximal derivative

wait : bool
  True -> Waits until the previous Compile call is finished before start compiling

keep_files : bool
  True -> Keep temporary files

)raw_string"))


    .def_property_readonly("data",
                           [] (shared_ptr<CF> cf)
                           {
                             PyArchive<BinaryOutArchive> ar;
                             ar << cf.get();
                             auto list = ar.WriteOut();
                             py::dict ret;
                             ret["version_needed"] = list.attr("pop")();
                             ret["version_stored"] = list.attr("pop")();
                             ret["data"] = list.attr("pop")();
                             ret["childs"] = list;
                             return ret;
                           })
    
    .def (NGSPickle<CoefficientFunction>())
    ;

  m.def("Cross", [] (shared_ptr<CF> cf1, shared_ptr<CF> cf2) { return CrossProduct(cf1, cf2); });
  m.def("Sym", [] (shared_ptr<CF> cf) { return SymmetricCF(cf); });
  m.def("Skew", [] (shared_ptr<CF> cf) { return SkewCF(cf); });
  m.def("Trace", [] (shared_ptr<CF> cf) { return TraceCF(cf); });
  m.def("Id", [] (int dim, int order) { return IdentityCF(dim, order); }, py::arg("dim"), py::arg("tensor_order") = 2,
        "Identity matrix of given dimension");
  m.def("Id", [] (const Array<int>& dims) { return IdentityCF(dims); }, py::arg("dims"),
        "Identity tensor for a space with dimensions 'dims', ie. the result is of 'dims + dims'");
  m.def("Zero", [] (const Array<int>& dims) { return ZeroCF(move(dims)); });
  m.def("Inv", [] (shared_ptr<CF> cf) { return InverseCF(cf); });
  m.def("Cof", [] (shared_ptr<CF> cf) { return CofactorCF(cf); });
  m.def("Det", [] (shared_ptr<CF> cf) { return DeterminantCF(cf); });
  m.def("Conj", [] (shared_ptr<CF> cf) { return ConjCF(cf); }, "complex-conjugate");

  m.def("MinimizationCF", [](shared_ptr<CF> expr, py::object startingpoint, optional<double> tol,
                       optional<double> rtol, optional<int> maxiter){

          // First check for a GridFunction. This case is important for GFs on
          // (abstract) compound spaces.
          py::extract<shared_ptr<ngcomp::GridFunction>> egf(startingpoint);
          if (egf.check()) {
            const auto gf = egf();
            if (gf->GetFESpace()->GetEvaluator())
              return CreateMinimizationCF(expr, gf, tol, rtol, maxiter);
            else {
              // Probably a GF on a generic compound space
              Array<shared_ptr<CoefficientFunction>> stps(gf->GetNComponents());
              for (int comp : Range(stps))
                stps[comp] = gf->GetComponent(comp);
              return CreateMinimizationCF(expr, stps, tol, rtol, maxiter);
            }
          }

          // If the given value is a list or a tuple... This must be handled
          // explicitly to prevent implicit conversion to a CoefficientFunction,
          // which is problematic when a list of GridFunctions is given.
          if (py::isinstance<py::list>(startingpoint) || py::isinstance<py::tuple>(startingpoint))
          {
            Array<shared_ptr<CoefficientFunction>> stps{};
            for (auto val : startingpoint)
            {
              py::extract<shared_ptr<CoefficientFunction>> vcf(val);
              if (vcf.check())
                stps.Append(vcf());
              else
                throw std::invalid_argument(
                    string("Cannot make CoefficientFunction from ")
                    + string(py::str(val)) + " of type "
                    + string(py::str(val.get_type())));
            }
            return CreateMinimizationCF(expr, stps, tol, rtol, maxiter);
          }

          // Attemp std. conversion to a CoefficientFunction
          py::extract<shared_ptr<CoefficientFunction>> ecf(startingpoint);
          if (ecf.check())
            return CreateMinimizationCF(expr, ecf(), tol, rtol, maxiter);

          throw std::invalid_argument(
              string("Failed to convert startingpoint ")
              + string(py::str(startingpoint))
              + " to a CoefficientFunction");
        },
        py::arg("expression"), py::arg("startingpoint"), py::arg("tol") = 1e-6,
        py::arg("rtol") = 0.0, py::arg("maxiter") = 20, docu_string(R"raw_string(
Creates a CoefficientFunction that returns the solution to a minimization problem.
Convergence failure is indicated by returning NaN(s).

Parameters:

expression : CoefficientFunction
  the objective function to be minimized

startingpoint: CoefficientFunction, list/tuple of CoefficientFunctions
  the initial guess for the iterative solution of the nonlinear problem

tol: double
  absolute tolerance

rtol: double
  relative tolerance

maxiter: int
  maximum iterations

)raw_string"));

  
  m.def("Einsum", [](string index_signature, py::args args, const py::kwargs &kwargs) {
            Array<shared_ptr<CoefficientFunction>> cfs(args.size());
            for (auto i : Range(cfs.Size()))
                cfs[i] = py::extract<shared_ptr<CoefficientFunction>>(args[i])();
            map<string, bool> options{};
            for (const auto kv : kwargs)
                options[py::extract<string>(kv.first)()] = py::extract<bool>(kv.second)();
            return EinsumCF(move(index_signature), move(cfs), options);
        }, py::arg("einsum_signature"), docu_string(R"raw_string(
Generic tensor product in the spirit of numpy's \"einsum\" feature.

Parameters:

einsum_signature: str
  specification of the tensor product in numpy's "einsum" notation
  
args: 
  CoefficientFunctions
  
kwargs:
  "expand_einsum" (true) -- expand nested "einsums" for later optimization
  "optimize_path" (false) -- try to reorder product for greater efficiency
  "optimize_identities" (false) -- try to eliminate identity tensors
  "use_legacy_ops" (false) -- fall back to existing CFs implementing certain blas operations where possible

)raw_string"));
  
  m.def("LeviCivitaSymbol", &LeviCivitaCF);
  
  m.def("NewtonCF", [](shared_ptr<CF> expr, py::object startingpoint, optional<double> tol,
                       optional<double> rtol, optional<int> maxiter){

          // First check for a GridFunction. This case is important for GFs on
          // generic compound spaces.
          py::extract<shared_ptr<ngcomp::GridFunction>> egf(startingpoint);
          if (egf.check()) {
            const auto gf = egf();
            if (gf->GetFESpace()->GetEvaluator())
              return CreateNewtonCF(expr, gf, tol, rtol, maxiter);
            else {
              // Probably a GF on a generic compound space
              Array<shared_ptr<CoefficientFunction>> stps(gf->GetNComponents());
              for (int comp : Range(stps))
                stps[comp] = gf->GetComponent(comp);
              return CreateNewtonCF(expr, stps, tol, rtol, maxiter);
            }
          }

          // If the given value is a list or a tuple... This must be handled
          // explicitly to prevent implicit conversion to a CoefficientFunction,
          // which is problematic when a list of GridFunctions is given.
          if (py::isinstance<py::list>(startingpoint) || py::isinstance<py::tuple>(startingpoint))
          {
            Array<shared_ptr<CoefficientFunction>> stps{};
            for (auto val : startingpoint)
            {
              py::extract<shared_ptr<CoefficientFunction>> vcf(val);
              if (vcf.check())
                stps.Append(vcf());
              else
                throw std::invalid_argument(
                    string("Cannot make CoefficientFunction from ")
                    + string(py::str(val)) + " of type "
                    + string(py::str(val.get_type())));
            }
            return CreateNewtonCF(expr, stps, tol, rtol, maxiter);
          }

          // Attemp std. conversion to a CoefficientFunction
          py::extract<shared_ptr<CoefficientFunction>> ecf(startingpoint);
          if (ecf.check())
            return CreateNewtonCF(expr, ecf(), tol, rtol, maxiter);

          throw std::invalid_argument(
              string("Failed to convert startingpoint ")
              + string(py::str(startingpoint))
              + " to a CoefficientFunction");
      },
      py::arg("expression"), py::arg("startingpoint"), py::arg("tol") = 1e-6,
      py::arg("rtol") = 0.0, py::arg("maxiter") = 10, docu_string(R"raw_string(
Creates a CoefficientFunction that returns the solution to a nonlinear problem.
Convergence failure is indicated by returning NaN(s).

Parameters:

expression : CoefficientFunction
  the residual of the nonlinear equation

startingpoint: CoefficientFunction, list/tuple of CoefficientFunctions
  the initial guess for the iterative solution of the nonlinear problem

tol: double
  absolute tolerance

rtol: double
  relative tolerance

maxiter: int
  maximum iterations

)raw_string"));

  py::implicitly_convertible<double, CoefficientFunction>();
  py::implicitly_convertible<Complex, CoefficientFunction>();
  py::implicitly_convertible<int, CoefficientFunction>();
  py::implicitly_convertible<py::tuple, CoefficientFunction>();
  py::implicitly_convertible<py::list, CoefficientFunction>();

  
  py::class_<ngcomp::PointEvaluationFunctional, shared_ptr<ngcomp::PointEvaluationFunctional>> (m, "PointEvaluationFunctional")
    .def("Assemble", &ngcomp::PointEvaluationFunctional::Assemble)
    ;

  
  if(have_numpy)
    {
      cf_class.def("__call__", [](shared_ptr<CF> self, py::array_t<MeshPoint> points) -> py::array
         {
           auto pts = points.unchecked<1>(); // pts has array access without bounds checks
           size_t npoints = pts.shape(0);
           py::array np_array;
           if (!self->IsComplex())
             {
               Array<double> vals(npoints * self->Dimension());
               constexpr size_t maxp = 16;
               ParallelForRange(Range(npoints), [&](IntRange r)
                           {
                             LocalHeapMem<50000> lh("CF evaluate");
                             Matrix<SIMD<double>> simdvals(self->Dimension(), maxp / SIMD<double>::Size());
                             IntegrationRule ir;

                             try
                               {
                                 for (auto i = r.begin(); i < r.end(); )
                                   {
                                     HeapReset hr(lh);
                                     auto& mp = pts(i);
                                     if(mp.nr == -1)
                                       throw Exception("Meshpoint " + to_string(i) + " not in mesh!");
                                     auto ei = ElementId(mp.vb, mp.nr);
                                     auto& trafo = mp.mesh->GetTrafo(ei, lh);
                                     ir.SetSize(0);
                                     
                                     auto first = i;
                                     ir.Append (IntegrationPoint(mp.x, mp.y, mp.z));
                                     i++;
                                     while (i < r.end() && ElementId(pts(i).vb, pts(i).nr) == ei && i < first+maxp)
                                       {
                                         ir.Append(IntegrationPoint(pts(i).x, pts(i).y, pts(i).z));
                                         i++;
                                       }
                                     SIMD_IntegrationRule simd_ir(ir, lh);
                                     auto& mir = trafo(simd_ir, lh);                                 
                                     self->Evaluate(mir, simdvals.Cols(0, simd_ir.Size()));
                                     
                                     FlatMatrix<double> fm(ir.Size(), self->Dimension(), &vals[first*self->Dimension()]);
                                     SliceMatrix<> simdfm(self->Dimension(), ir.Size(), simdvals.Width()*SIMD<double>::Size(),
                                                          (double*)&simdvals(0,0));
                                     fm = Trans(simdfm);
                                   }
                               }
                             catch (ExceptionNOSIMD e)
                               {
                                 for (auto i = r.begin(); i < r.end(); )
                                   {
                                     HeapReset hr(lh);
                                     auto& mp = pts(i);
                                     if(mp.nr == -1)
                                       throw Exception("Meshpoint " + to_string(i) + " not in mesh!");
                                     auto ei = ElementId(mp.vb, mp.nr);
                                     auto& trafo = mp.mesh->GetTrafo(ei, lh);
                                     ir.SetSize(0);
                                     
                                     auto first = i;
                                     ir.Append (IntegrationPoint(mp.x, mp.y, mp.z));
                                     i++;
                                     while (i < r.end() && ElementId(pts(i).vb, pts(i).nr) == ei && i < first+maxp)
                                       {
                                         ir.Append(IntegrationPoint(pts(i).x, pts(i).y, pts(i).z));
                                         i++;
                                       }
                                     auto& mir = trafo(ir, lh);
                                     FlatMatrix<double> fm(ir.Size(), self->Dimension(), &vals[first*self->Dimension()]);
                                     self->Evaluate(mir, fm);
                                   }
                               }
                           });
               np_array = MoveToNumpyArray(vals);
             }
           else
             {
               Array<Complex> vals(npoints * self->Dimension());
               ParallelFor(Range(npoints), [&](size_t i)
                           {
                             LocalHeapMem<1000> lh("CF evaluate");
                             auto& mp = pts(i);
                             auto& trafo = mp.mesh->GetTrafo(ElementId(mp.vb, mp.nr), lh);
                             auto& mip = trafo(IntegrationPoint(mp.x,mp.y,mp.z),lh);
                             FlatVector<Complex> fv(self->Dimension(), &vals[i*self->Dimension()]);
                             self->Evaluate(mip, fv);
                           });
               np_array = MoveToNumpyArray(vals);
             }
           return np_array.attr("reshape")(npoints, self->Dimension());
         });
    }

  typedef shared_ptr<ParameterCoefficientFunction<double>> spParameterCF;
  py::class_<ParameterCoefficientFunction<double>, spParameterCF, CF>
    (m, "Parameter", docu_string(R"raw_string(
CoefficientFunction with a modifiable value

Parameters:

value : float
  Parameter value

)raw_string"))
    .def (py::init<double>(),
          py::arg("value"), "Construct a ParameterCF from a scalar")
    .def(NGSPickle<ParameterCoefficientFunction<double>>())
    .def ("Set", [] (spParameterCF cf, double val)  { cf->SetValue (val); }, py::arg("value"),
          docu_string(R"raw_string(
Modify parameter value.

Parameters:

value : double
  input scalar  

)raw_string"))
    .def ("Get", [] (spParameterCF cf)  { return cf->GetValue(); },
          "return parameter value")
    ;

  using spParCFC = shared_ptr<ParameterCoefficientFunction<Complex>>;
  py::class_<ParameterCoefficientFunction<Complex>, spParCFC, CF>
    (m, "ParameterC", docu_string(R"raw_string(
CoefficientFunction with a modifiable complex value

Parameters:

value : complex
  Parameter value

)raw_string"))
    .def (py::init<Complex>(),
          py::arg("value"), "Construct a ParameterCF from a scalar")
    .def(NGSPickle<ParameterCoefficientFunction<Complex>>())
    .def("Set", [] (spParCFC cf, Complex val)
      { cf->SetValue (val); },
      py::arg("value"),
      docu_string(R"raw_string(
Modify parameter value.

Parameters:

value : complex
  input scalar

)raw_string"))
    .def ("Get", [] (spParameterCF cf)  { return cf->GetValue(); },
          "return parameter value")
    ;

  py::class_<PlaceholderCoefficientFunction,
             shared_ptr<PlaceholderCoefficientFunction>,
             CoefficientFunction>
    (m, "PlaceholderCF")
    .def(py::init<shared_ptr<CF>>())
    .def("Set", &PlaceholderCoefficientFunction::Set)
    ;

  py::class_<BSpline, shared_ptr<BSpline> > (m, "BSpline",R"raw(
BSpline of arbitrary order

Parameters:

order : int
  order of the BSpline

knots : list
  list of float

vals : list
  list of float

)raw")
    .def(py::init
         ([](int order, py::list knots, py::list vals)
          {
            return make_shared<BSpline> (order,
                                         makeCArray<double> (knots),
                                         makeCArray<double> (vals));
          }), py::arg("order"), py::arg("knots"), py::arg("vals"),
        "B-Spline of a certain order, provide knot and value vectors")
    .def("__str__", &ToString<BSpline>)
    .def("__call__", [] (shared_ptr<BSpline> self, double pt)
          {
            return self->Evaluate(pt);
          },
      py::arg("pt"),
      py::return_value_policy::move)
    .def("__call__", [](shared_ptr<BSpline> sp, shared_ptr<CF> coef)
          {
            return UnaryOpCF (coef, GenericBSpline(sp) /* , GenericBSpline(sp) */);
          }, py::arg("cf"))
    .def("Integrate", 
         [](const BSpline & sp) { return make_shared<BSpline>(sp.Integrate()); }, "Integrate the BSpline")
    .def("Differentiate", 
         [](const BSpline & sp) { return make_shared<BSpline>(sp.Differentiate()); }, "Differentiate the BSpline")
    ;

  m.def ("LoggingCF", LoggingCF, py::arg("cf"), py::arg("logfile")="stdout");
  m.def ("CacheCF", CacheCF, py::arg("cf"));
}


// *************************************** Export FEM ********************************


void NGS_DLL_HEADER ExportNgfem(py::module &m) {

  py::enum_<ELEMENT_TYPE>(m, "ET", "Enumeration of all supported element types.")
    .value("POINT", ET_POINT)     .value("SEGM", ET_SEGM)
    .value("TRIG", ET_TRIG)       .value("QUAD", ET_QUAD)
    .value("TET", ET_TET)         .value("PRISM", ET_PRISM)
    .value("PYRAMID", ET_PYRAMID) .value("HEX", ET_HEX)
    .export_values()
    ;

  py::enum_<NODE_TYPE>(m, "NODE_TYPE", "Enumeration of all supported node types.")
    .value("VERTEX", NT_VERTEX)
    .value("EDGE", NT_EDGE)
    .value("FACE", NT_FACE)
    .value("CELL", NT_CELL)
    .value("ELEMENT", NT_ELEMENT)
    .value("FACET", NT_FACET)
    .export_values()
    ;


  py::class_<ElementTopology> (m, "ElementTopology", docu_string(R"raw_string(
Element Topology

Parameters:

et : ngsolve.fem.ET
  input element type

)raw_string"))
    .def(py::init<ELEMENT_TYPE>(), py::arg("et"))
    .def_property_readonly("name", 
                           static_cast<const char*(ElementTopology::*)()> (&ElementTopology::GetElementName), "Name of the element topology")
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
                                              }, "Vertices of the element topology");
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

  py::class_<MixedFiniteElement, shared_ptr<MixedFiniteElement>, FiniteElement>
    (m, "MixedFE", "pair of finite elements for trial and test-functions")
    .def(py::init<const FiniteElement&, const FiniteElement&>())
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
         py::arg("x"),py::arg("y")=0.0,py::arg("z")=0.0,docu_string(R"raw_string(
Parameters:

x : double
  input x value

y : double
  input y value

z : double
  input z value

)raw_string"))
    .def("CalcShape",
         [] (const BaseScalarFiniteElement & fe, const BaseMappedIntegrationPoint & mip)
          {
            Vector<> v(fe.GetNDof());
            fe.CalcShape (mip.IP(), v);
            return v;
          },
         py::arg("mip"),docu_string(R"raw_string(
Parameters:

mip : ngsolve.BaseMappedIntegrationPoint
  input mapped integration point

)raw_string"))
    .def("CalcDShape",
         [] (const BaseScalarFiniteElement & fe, const BaseMappedIntegrationPoint & mip)
          {
            Matrix<> mat(fe.GetNDof(), fe.Dim());
            fe.CalcMappedDShape(mip, mat);
            /*
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
            */
            return mat;
          },
         py::arg("mip"),docu_string(R"raw_string(
Computes derivative of the shape in an integration point.

Parameters:

mip : ngsolve.BaseMappedIntegrationPoint
  input mapped integration point

)raw_string"))
        .def("CalcDShape",
         [] (const BaseScalarFiniteElement & fe, double x, double y, double z)
         {
           IntegrationPoint ip(x,y,z);
           Matrix<> mat(fe.GetNDof(), fe.Dim());
           fe.CalcDShape(ip, mat);
           return mat;
          },
         py::arg("x"),py::arg("y")=0.0,py::arg("z")=0.0,docu_string(R"raw_string(
Computes derivative of the shape in an integration point.

Parameters:

x : double
  input x value

y : double
  input y value

z : double
  input z value

)raw_string"))
    ;


  py::class_<BaseHCurlFiniteElement, shared_ptr<BaseHCurlFiniteElement>, 
             FiniteElement>
    (m, "HCurlFE", "an H(curl) finite element")
    .def("CalcShape",
         [] (const BaseHCurlFiniteElement & fe, double x, double y, double z)
         {
           IntegrationPoint ip(x,y,z);
           Matrix<> mat(fe.GetNDof(), fe.Dim());
           fe.CalcShape (ip, mat);
           return mat;
         },
         py::arg("x"),py::arg("y")=0.0,py::arg("z")=0.0)
    .def("CalcShape",
         [] (const BaseHCurlFiniteElement & fe, const BaseMappedIntegrationPoint & mip)
          {
            Matrix<> mat(fe.GetNDof(), fe.Dim());
            fe.CalcMappedShape (mip, mat);
            return mat;
          },
         py::arg("mip"))
    .def("CalcCurlShape",
         [] (const BaseHCurlFiniteElement & fe, const BaseMappedIntegrationPoint & mip)
          {
            Matrix<> mat(fe.GetNDof(), fe.Dim());
            fe.CalcMappedCurlShape(mip, mat);
            return mat;
          },
         py::arg("mip"))
    ;

    py::class_<BaseHDivFiniteElement, shared_ptr<BaseHDivFiniteElement>, 
             FiniteElement>
    (m, "HDivFE", "an H(div) finite element")
    .def("CalcShape",
         [] (const BaseHDivFiniteElement & fe, double x, double y, double z)
         {
           IntegrationPoint ip(x,y,z);
           Matrix<> mat(fe.GetNDof(), fe.Dim());
           fe.CalcShape (ip, mat);
           return mat;
         },
         py::arg("x"),py::arg("y")=0.0,py::arg("z")=0.0)
     .def("CalcDivShape",
         [] (const BaseHDivFiniteElement & fe, double x, double y, double z)
         {
           IntegrationPoint ip(x,y,z);
           Vector<> v(fe.GetNDof());
           fe.CalcDivShape(ip,v);
           return v;
         },
         py::arg("x"),py::arg("y")=0.0,py::arg("z")=0.0)
    ;

    py::class_<BaseHDivDivFiniteElement, shared_ptr<BaseHDivDivFiniteElement>, 
               FiniteElement>
      (m, "HDivDivFE", "an H(div div) finite element")
      .def("CalcShape",
         [] (const BaseHDivDivFiniteElement & fe, double x, double y, double z)
         {
           IntegrationPoint ip(x,y,z);
           Matrix<> mat(fe.GetNDof(), fe.Dim()*(1+fe.Dim())/2);
           fe.CalcShape (ip, mat);
           return mat;
         },
         py::arg("x"),py::arg("y")=0.0,py::arg("z")=0.0)
     .def("CalcDivShape",
         [] (const BaseHDivDivFiniteElement & fe, double x, double y, double z)
         {
           IntegrationPoint ip(x,y,z);
           Matrix<> mat(fe.GetNDof(), fe.Dim());
           fe.CalcDivShape (ip, mat);
           return mat;
         },
         py::arg("x"),py::arg("y")=0.0,py::arg("z")=0.0)
    ;


  py::implicitly_convertible <BaseScalarFiniteElement, FiniteElement >(); 

  m.def("H1FE", [](ELEMENT_TYPE et, int order)
        {
          SwitchET (et, [order] (auto et2) -> shared_ptr<BaseScalarFiniteElement>
                    {
                      constexpr ELEMENT_TYPE ET = et2.ElementType();
                      return make_shared<H1HighOrderFE<ET>> (order);
                    });
        },
        py::arg("et"), py::arg("order"),
        docu_string(R"raw_string(Creates an H1 finite element of given geometric shape and polynomial order.

Parameters:

et : ngsolve.fem.ET
  input element type

order : int
  input polynomial order

)raw_string")
          );

  m.def("L2FE", [](ELEMENT_TYPE et, int order)
        {
          SwitchET (et, [order] (auto et2) -> shared_ptr<BaseScalarFiniteElement>
                    {
                      constexpr ELEMENT_TYPE ET = et2.ElementType();                      
                      return make_shared<L2HighOrderFE<ET>> (order);
                    });
        },
        py::arg("et"), py::arg("order"),
        docu_string(R"raw_string(Creates an L2 finite element of given geometric shape and polynomial order.

Parameters:

et : ngsolve.fem.ET
  input element type

order : int
  input polynomial order

)raw_string")
          );


  py::class_<IntegrationPoint>(m, "IntegrationPoint")
    .def_property_readonly("point", [](IntegrationPoint& self)
                           {
                             return py::make_tuple(self.Point()[0], self.Point()[1], self.Point()[2]);
                           }, "Integration point coordinates as tuple, has always x,y and z component, which do not have meaning in lesser dimensions")
    .def_property_readonly("weight", &IntegrationPoint::Weight, "Weight of the integration point")
    ;

  py::class_<IntegrationRule>(m, "IntegrationRule", docu_string(R"raw_string(
Integration rule

2 __init__ overloads


1)

Parameters:

element type : ngsolve.fem.ET
  input element type

order : int
  input order of integration rule


2)

Parameters:

points : list
  input list of integration points

weights : list
  input list of integration weights

)raw_string"))
    .def(py::init
         ([](ELEMENT_TYPE et, int order)
          { return new IntegrationRule (et, order); }),
         py::arg("element type"), py::arg("order"))
    
    .def(py::init
         ([](py::list points, py::list weights)
          {
            IntegrationRule * ir = new IntegrationRule ();
            for (size_t i = 0; i < len(points); i++)
              {
                py::object pnt = points[i];
                IntegrationPoint ip;
                ip.SetNr(i);
                for (int j = 0; j < len(pnt); j++)
                  ip(j) = py::extract<double> (py::tuple(pnt)[j])();
                if(weights.size() == 0)
                  ip.SetWeight(0.);
                else
                  ip.SetWeight(py::extract<double> (weights[i])());
                ir -> Append (ip);
              }
            return ir;
          }),
         py::arg("points"), py::arg("weights") = py::list())
    .def("__str__", &ToString<IntegrationRule>)
    .def("__getitem__", [](IntegrationRule & ir, int nr)
                                        {
                                          if (nr < 0 || nr >= ir.Size())
                                            throw py::index_error();
                                          return ir[nr];
                                        }, py::arg("nr"), "Return integration point at given position")
    .def("__len__", &IntegrationRule::Size)
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
          }, py::arg("func"), "Integrates a given function")
    .def_property_readonly("weights", [] (IntegrationRule& self)
                           {
                             py::list weights;
                             for (auto ip : self)
                               weights.append(ip.Weight());
                             return weights;
                           }, "Weights of IntegrationRule")
    .def_property_readonly("points", [] (IntegrationRule& self)
                           {
                             py::list points;
                             for(auto ip : self)
                               switch(self.Dim())
                                 {
                                 case 1:
                                   points.append(py::make_tuple(ip.Point()[0]));
                                   break;
                                 case 2:
                                   points.append(py::make_tuple(ip.Point()[0],
                                                                ip.Point()[1]));
                                   break;
                                 default:
                                   points.append(py::make_tuple(ip.Point()[0],
                                                                ip.Point()[1],
                                                                ip.Point()[2]));
                                 }
                             return points;
                           }, "Points of IntegrationRule as tuple")
    ;


  py::class_<MeshPoint>(m, "MeshPoint")
    .def_property_readonly("pnt", [](MeshPoint& p) { return py::make_tuple(p.x,p.y,p.z); }, "Gives coordinates of point on reference triangle. One can create a MappedIntegrationPoint using the ngsolve.fem.BaseMappedIntegrationPoint constructor. For physical coordinates the coordinate CoefficientFunctions x,y,z can be evaluated in the MeshPoint")
    .def_property_readonly("mesh", [](MeshPoint& p) { return p.mesh; })
    .def_property_readonly("vb", [](MeshPoint& p) { return p.vb; })
    .def_property_readonly("nr", [](MeshPoint& p) { return p.nr; })
    ;

  if (have_numpy)
  {
    py::detail::npy_format_descriptor<MeshPoint>::register_dtype({
        py::detail::field_descriptor { "x", offsetof(MeshPoint, x), sizeof(double),
            py::format_descriptor<double>::format(), py::detail::npy_format_descriptor<double>::dtype() },
          py::detail::field_descriptor { "y", offsetof(MeshPoint, y), sizeof(double),
              py::format_descriptor<double>::format(), py::detail::npy_format_descriptor<double>::dtype() },
            py::detail::field_descriptor { "z", offsetof(MeshPoint, z), sizeof(double),
                py::format_descriptor<double>::format(), py::detail::npy_format_descriptor<double>::dtype() },
            py::detail::field_descriptor { "meshptr", offsetof(MeshPoint, mesh), sizeof(double),
                py::format_descriptor<double>::format(), py::detail::npy_format_descriptor<double>::dtype() },
            py::detail::field_descriptor { "VorB", offsetof(MeshPoint, vb), sizeof(int),
                py::format_descriptor<int>::format(), py::detail::npy_format_descriptor<int>::dtype() },
              py::detail::field_descriptor {"nr", offsetof(MeshPoint, nr), sizeof(int),
                  py::format_descriptor<int>::format(), py::detail::npy_format_descriptor<int>::dtype()}});
  }

  py::class_<BaseMappedIntegrationPoint>(m, "BaseMappedIntegrationPoint")
    .def(py::init([](MeshPoint& pnt)
                  {
                    if(pnt.nr == -1){
                      cout << "WARNING: MeshPoint not in mesh, can't convert to BaseMappedIntegrationPoint!" << endl;
                      throw Exception("Meshpoint at (" + to_string(pnt.x) + ", " + to_string(pnt.y) + ", " +
                                      to_string(pnt.z) + ") not in mesh!");
                    }
                    auto& trafo = pnt.mesh->GetTrafo(ElementId(pnt.vb, pnt.nr), global_alloc);
                    auto& mip = trafo(IntegrationPoint(pnt.x,pnt.y,pnt.z),global_alloc);
                    mip.SetOwnsTrafo(true);
                    return &mip;
                  }))
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
    .def_property_readonly("measure", &BaseMappedIntegrationPoint::GetMeasure, "Measure of the mapped integration point ")
    .def_property_readonly("point", &BaseMappedIntegrationPoint::GetPoint, "Point of the mapped integration point")
    .def_property_readonly("jacobi", &BaseMappedIntegrationPoint::GetJacobian, "jacobian of the mapped integration point")
    // .def_property_readonly("trafo", &BaseMappedIntegrationPoint::GetTransformation)
    .def_property_readonly("trafo", &BaseMappedIntegrationPoint::GetTransformation, "Transformation of the mapped integration point")
    .def_property_readonly("elementid", [](BaseMappedIntegrationPoint & mip)
                                               {
                                                 return mip.GetTransformation().GetElementId();
                                               }, "Element ID of the mapped integration point")
    ;

  py::implicitly_convertible<MeshPoint, BaseMappedIntegrationPoint>();

  auto eltrans_class = py::class_<ElementTransformation, shared_ptr<ElementTransformation>>(m, "ElementTransformation")
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
    .def_property_readonly("VB", &ElementTransformation::VB, "VorB (VOL, BND, BBND, BBBND)")
    .def_property_readonly("spacedim", &ElementTransformation::SpaceDim, "Space dimension of the element transformation")
    .def_property_readonly("elementid", &ElementTransformation::GetElementId, "Element ID of the element transformation")
    .def_property_readonly("curved", &ElementTransformation::IsCurvedElement, "Is mapping non-affine ?")    
    .def ("__call__", [] (shared_ptr<ElementTransformation> self, double x, double y, double z)
           {
             return &(*self)(IntegrationPoint(x,y,z), global_alloc);
           },
          py::arg("x"), py::arg("y")=0, py::arg("z")=0,
          py::return_value_policy::reference)
    .def ("__call__", [] (shared_ptr<ElementTransformation> self, IntegrationPoint & ip)
           {
             return &(*self)(ip, global_alloc);
           }, py::arg("ip"), py::return_value_policy::reference)
    ;

  if(have_numpy)
    {
      eltrans_class.def("__call__", [](ElementTransformation& self, IntegrationRule& ir)
                                    {
                                      Array<MeshPoint> pts;
                                      pts.SetAllocSize(ir.Size());
                                      for(auto& ip : ir)
                                        pts.Append(MeshPoint{ip(0), ip(1),
                                                             ip(2),
                                                             (ngcomp::MeshAccess*) self.GetMesh(),
                                                             self.VB(),
                                                             self.GetElementNr()});
                                      return MoveToNumpyArray(pts);
                                    });
    }


  py::class_<DifferentialOperator, shared_ptr<DifferentialOperator>>
    (m, "DifferentialOperator")
    .def("__timing__", [&] (const DifferentialOperator & diffop,
                            const FiniteElement & fel, const ElementTransformation & trafo,
                            const IntegrationRule & ir)
         {
           LocalHeap lh(1000000);
           auto & mir = trafo(ir, lh);
           return diffop.Timing (fel, mir);
         })
    ;

  typedef BilinearFormIntegrator BFI;
  auto bfi_class = py::class_<BFI, shared_ptr<BFI>> (m, "BFI", docu_string(R"raw_string(
Bilinear Form Integrator

Parameters:

name : string
  Name of the bilinear form integrator.

py_coef : object
  CoefficientFunction of the bilinear form.

dim : int
  dimension of the bilinear form integrator

imag : bool
  Multiplies BFI with 1J

filename : string
  filename 

kwargs : kwargs
  For a description of the possible kwargs have a look a bit further down.

)raw_string"), py::dynamic_attr());
  bfi_class
    .def(py::init([bfi_class] (const string name, py::object py_coef, int dim, bool imag,
                      string filename, py::kwargs kwargs)
        -> shared_ptr<BilinearFormIntegrator>
        {
          auto flags = CreateFlagsFromKwArgs(kwargs, bfi_class);
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
    .def_property("simd_evaluate",
                  [](shared_ptr<BFI> self) { return self->SimdEvaluate(); },
                  [](shared_ptr<BFI> self, bool b) { return self->SetSimdEvaluate(b); },                  
                  "SIMD evaluate ?"
      )

    .def("__str__",  [](shared_ptr<BFI> self) { return ToString<BilinearFormIntegrator>(*self); } )

    .def("Evaluator",  [](shared_ptr<BFI> self, string name ) { return self->GetEvaluator(name); }, py::arg("name"), docu_string(R"raw_string(
Returns requested evaluator

Parameters:

name : string
  input name of requested evaluator

)raw_string") )
    // .def("DefinedOn", &Integrator::DefinedOn)
    .def("GetDefinedOn",  [] (shared_ptr<BFI> self) -> const BitArray &{ return self->GetDefinedOn(); } ,
         py::return_value_policy::reference, "Returns a BitArray where the bilinear form is defined on")

    .def("SetDefinedOnElements",  [](shared_ptr<BFI> self, shared_ptr<BitArray> ba )
         { self->SetDefinedOnElements (ba); }, py::arg("bitarray"), docu_string(R"raw_string( 
Set the elements on which the bilinear form is defined on.

Parameters:

bitarray : ngsolve.ngstd.BitArray
  input bitarray

)raw_string") )
    .def("SetIntegrationRule", [] (shared_ptr<BFI> self, ELEMENT_TYPE et, IntegrationRule ir)
         {
           self -> SetIntegrationRule(et,ir);
           return self;
         }, py::arg("et"), py::arg("intrule"), docu_string(R"raw_string( 
Set integration rule of the bilinear form.

Parameters:

et : ngsolve.fem.Element_Type
  input element type

intrule : ngsolve.fem.Integrationrule
  input integration rule

)raw_string"))
    .def("CalcElementMatrix",
         [] (shared_ptr<BFI> self,
             const FiniteElement & fe, const ElementTransformation &trafo,
             size_t heapsize, bool complex)
                         {
                           while (true)
                             {
                               LocalHeap lh(heapsize);
                               try
                                 {
                                   if (complex)
                                     {                                       
                                       Matrix<Complex> mat(fe.GetNDof() * self->GetDimension());
                                       self->CalcElementMatrix(fe,trafo,mat,lh);
                                       return py::cast(mat);
                                     }
                                   else
                                     {
                                       const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fe);
                                       const FiniteElement & fe_trial = mixedfe ? mixedfe->FETrial() : fe;
                                       const FiniteElement & fe_test = mixedfe ? mixedfe->FETest() : fe;
                                       
                                       Matrix<> mat(fe_test.GetNDof() * self->GetDimension(),
                                                    fe_trial.GetNDof() * self->GetDimension());
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
         py::arg("fel"),py::arg("trafo"),py::arg("heapsize")=10000, py::arg("complex") = false, docu_string(R"raw_string( 
Calculate element matrix of a specific element.

Parameters:

fel : ngsolve.fem.FiniteElement
  input finite element

trafo : ngsolve.fem.ElementTransformation
  input element transformation

heapsize : int
  input heapsize

complex : bool
  input complex

)raw_string"))
    .def("CalcLinearizedElementMatrix",
         [] (shared_ptr<BFI> self,
             const FiniteElement & fe, FlatVector<double> vec, 
             const ElementTransformation &trafo,
             size_t heapsize)
                         {
                           while (true)
                             {
                               LocalHeap lh(heapsize);
                               try
                                 {
                                  const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fe);
                                  const FiniteElement & fe_trial = mixedfe ? mixedfe->FETrial() : fe;
                                  const FiniteElement & fe_test = mixedfe ? mixedfe->FETest() : fe;
                                  
                                  Matrix<> mat(fe_test.GetNDof() * self->GetDimension(),
                                              fe_trial.GetNDof() * self->GetDimension());
                                  self->CalcLinearizedElementMatrix (fe, trafo, vec, mat, lh);
                                  return py::cast(mat);
                                 }
                               catch (LocalHeapOverflow ex)
                                 {
                                   heapsize *= 10;
                                 }
                             }
                         },
         py::arg("fel"),py::arg("vec"),py::arg("trafo"),py::arg("heapsize")=10000, docu_string(R"raw_string( 
Calculate (linearized) element matrix of a specific element.

Parameters:

fel : ngsolve.fem.FiniteElement
  input finite element

vec : Vector
  linearization argument

trafo : ngsolve.fem.ElementTransformation
  input element transformation

heapsize : int
  input heapsize
)raw_string"))
    .def("ApplyElementMatrix",
         [] (shared_ptr<BFI> self,
             const FiniteElement & fe, FlatVector<double> vec, 
             const ElementTransformation &trafo,
             size_t heapsize)
                         {
                           while (true)
                             {
                               LocalHeap lh(heapsize);
                               try
                                 {
                                  const MixedFiniteElement * mixedfe = dynamic_cast<const MixedFiniteElement*> (&fe);
                                  // const FiniteElement & fe_trial = mixedfe ? mixedfe->FETrial() : fe;
                                  const FiniteElement & fe_test = mixedfe ? mixedfe->FETest() : fe;
                                  Vector<> vecy(fe_test.GetNDof() * self->GetDimension());
                                  self->ApplyElementMatrix (fe, trafo, vec, vecy, 0, lh);
                                  return py::cast(vecy);
                                 }
                               catch (LocalHeapOverflow ex)
                                 {
                                   heapsize *= 10;
                                 }
                             }
                         },
         py::arg("fel"),py::arg("vec"),py::arg("trafo"),py::arg("heapsize")=10000, docu_string(R"raw_string( 
Apply element matrix of a specific element.

Parameters:

fel : ngsolve.fem.FiniteElement
  input finite element

vec : Vector
  evaluation argument

trafo : ngsolve.fem.ElementTransformation
  input element transformation

heapsize : int
  input heapsize

)raw_string"))
    ;

  m.def("CompoundBFI", 
          []( shared_ptr<BFI> bfi, int comp )
                            {
                                return make_shared<CompoundBilinearFormIntegrator>(bfi, comp);
                            },
           py::arg("bfi")=NULL, py::arg("comp")=0, docu_string(R"raw_string(
Compound Bilinear Form Integrator

Parameters:

bfi : ngsolve.fem.BFI
  input bilinear form integrator

comp : int
  input component

)raw_string"));

  m.def("BlockBFI", 
          []( shared_ptr<BFI> bfi, int dim, int comp )
                            {
                                return make_shared<BlockBilinearFormIntegrator>(bfi, dim, comp);
                            },
           py::arg("bfi")=NULL, py::arg("dim")=2, py::arg("comp")=0
      , docu_string(R"raw_string(
Block Bilinear Form Integrator

Parameters:

bfi : ngsolve.fem.BFI
  input bilinear form integrator

dim : int
  input dimension of block bilinear form integrator

comp : int
  input comp

)raw_string"));

  typedef LinearFormIntegrator LFI;
  py::class_<LFI, shared_ptr<LFI>>
    (m, "LFI", docu_string(R"raw_string(
Linear Form Integrator

Parameters:

name : string
  Name of the linear form integrator.

dim : int
  dimension of the linear form integrator

coef : object
  CoefficientFunction of the bilinear form.

definedon : object
  input region where the linear form is defined on

imag : bool
  Multiplies LFI with 1J

flags : ngsolve.ngstd.Flags
  input flags

definedonelem : object
  input definedonelem

)raw_string"), py::dynamic_attr())
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
    .def_property("simd_evaluate",
                  [](shared_ptr<LFI> self) { return self->SimdEvaluate(); },
                  [](shared_ptr<LFI> self, bool b) { return self->SetSimdEvaluate(b); },                  
                  "SIMD evaluate ?")
    // .def("GetDefinedOn", &Integrator::GetDefinedOn)
    .def("GetDefinedOn",  [] (shared_ptr<LFI> self) -> const BitArray &{ return self->GetDefinedOn(); } ,
         py::return_value_policy::reference, "Reterns regions where the lienar form integrator is defined on.")
    .def("SetDefinedOnElements",  [](shared_ptr<LFI> self, shared_ptr<BitArray> ba )
         { self->SetDefinedOnElements (ba); }, py::arg("ba"), docu_string(R"raw_string(
Set the elements on which the linear form integrator is defined on

Parameters:

ba : ngsolve.ngstd.BitArray
  input bit array ( 1-> defined on, 0 -> not defoned on)

)raw_string"))
    .def("SetIntegrationRule", [](shared_ptr<LFI> self, ELEMENT_TYPE et, IntegrationRule ir)
         {
           self->SetIntegrationRule(et,ir);
           return self;
         }, py::arg("et"), py::arg("ir"), docu_string(R"raw_string(
Set a different integration rule for elements of type et

Parameters:

et : ngsolve.fem.ET
  input element type

ir : ngsolve.fem.IntegrationRule
  input integration rule

)raw_string"))

    .def("CalcElementVector", 
         static_cast<void(LinearFormIntegrator::*)(const FiniteElement&, const ElementTransformation&, FlatVector<double>, LocalHeap&)const>(&LinearFormIntegrator::CalcElementVector), py::arg("fel"), py::arg("trafo"), py::arg("vec"), py::arg("lh"))
    .def("CalcElementVector",
         [] (shared_ptr<LFI>  self, const FiniteElement & fe, const ElementTransformation& trafo,
             size_t heapsize, bool complex)
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
           "lfi"_a=NULL, py::arg("comp")=0, docu_string(R"raw_string(
Compound Linear Form Integrator

Parameters:

lfi : ngsolve.fem.LFI
  input linear form integrator

comp : int
  input component

)raw_string"));

  m.def("BlockLFI", 
          []( shared_ptr<LFI> lfi, int dim, int comp )
                            {
                                return shared_ptr<LFI>(make_shared<BlockLinearFormIntegrator>(lfi, dim, comp));
                            },
           "lfi"_a=NULL, py::arg("dim")=2, py::arg("comp")=0, docu_string(R"raw_string(
Block Linear Form Integrator

Parameters:

lfi : ngsolve.fem.LFI
  input bilinear form integrator

dim : int
  input dimension of block linear form integrator

comp : int
  input comp

)raw_string"));


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
         py::arg("rad")=1,py::arg("alpha")=1, docu_string(R"raw_string(
Parameters:

rad : double
  input radius of PML

alpha : double
  input damping factor of PML

)raw_string"));
    
                           
  m.def("GenerateL2ElementCode", &GenerateL2ElementCode);

  m.def("VoxelCoefficient",
        [](py::tuple pystart, py::tuple pyend, py::array values,
           bool linear, py::object trafocf)
        -> shared_ptr<CoefficientFunction>
        {
          shared_ptr<CoefficientFunction> trafo;
          try { trafo = MakeCoefficient(trafocf); }
          catch(...) { trafo=nullptr; }
          Array<string> allowed_types = { "float64", "complex128" };
          if(!allowed_types.Contains(py::cast<string>(values.dtype().attr("name"))))
            throw Exception("Only float64 and complex128 dtype arrays allowed!");
          Array<double> start, end;
          Array<size_t> dim_vals;
          for(auto val : pystart)
            start.Append(py::cast<double>(val));
          for(auto val : pyend)
            end.Append(py::cast<double>(val));
          for(auto dim : Range(values.ndim()))
            dim_vals.Insert(0,values.shape(dim));

          if(values.dtype().kind() == 'c')
            {
              auto c_array = py::cast<py::array_t<Complex>>(values.attr("ravel")());
              Array<Complex> vals(c_array.size());
              for(auto i : Range(vals))
                vals[i] = c_array.at(i);
              return make_shared<VoxelCoefficientFunction<Complex>>
                (start, end, dim_vals, move(vals), linear, trafo);
            }
          auto d_array = py::cast<py::array_t<double>>(values.attr("ravel")());
          Array<double> vals(values.size());
          for(auto i : Range(vals))
            vals[i] = d_array.at(i);
          return make_shared<VoxelCoefficientFunction<double>>
              (start, end, dim_vals, move(vals), linear, trafo);
        }, py::arg("start"), py::arg("end"), py::arg("values"),
        py::arg("linear")=true, py::arg("trafocf")=DummyArgument(), R"delimiter(CoefficientFunction defined on a grid.

Start and end mark the cartesian boundary of domain. The function will be continued by a constant function outside of this box. Inside a cartesian grid will be created by the dimensions of the numpy input array 'values'. This array must have the dimensions of the mesh and the values stored as:
x1y1z1, x2y1z1, ..., xNy1z1, x1y2z1, ...

If linear is True the function will be interpolated linearly between the values. Otherwise the nearest voxel value is taken.

)delimiter");

      const string header = R"CODE(
#include <comp.hpp>
#include <python_ngstd.hpp>

using namespace ngcomp;

extern "C" {

  NGCORE_API_EXPORT void init(py::object & res)
  {
    static py::module::module_def def;
    py::module m = py::module::create_extension_module("", "", &def);

    // BEGIN USER DEFINED CODE

)CODE";
      const string footer = R"CODE(
    // END USER DEFINED CODE
    res = m;
  }
}
)CODE";
      const string docu = R"raw_string(
Utility function to compile c++ code with python bindings at run-time.

Parameters:

code: c++ code snippet ( add_header=True ) or a complete .cpp file ( add_header=False )

init_function_name (default = "init"): Function, which is called after the compiled code is loaded. The prototype must match:
extern "C" void init_function_name(pybind11::object & res);

add_header (default = True): wrap the code snippet with the template
)raw_string" + header + footer;

  m.def("CompilePythonModule",
       [header, footer](string code, string init_function_name, bool add_header)
       {
           py::object result;
           typedef void (*init_function_type)(py::object & res);

           if(add_header)
             code = header + code + footer;

           auto library = CompileCode( {code}, {""} );
           auto func = library->GetFunction<init_function_type>(init_function_name);
           func(result);
           library.release(); // TODO: bind lifetime of "library" to python object "result"
           return result;
       },
       py::arg("code"),
       py::arg("init_function_name")="init",
       py::arg("add_header")=true,
       docu.c_str()
    );
  m.def("CompilePythonModule",
       [](filesystem::path file, string init_function_name)
       {
           py::object result;
           typedef void (*init_function_type)(py::object & res);

           auto library = CompileCode( {file}, {""} );
           auto func = library->GetFunction<init_function_type>(init_function_name);
           func(result);
           library.release(); // TODO: bind lifetime of "library" to python object "result"
           return result;
       },
       py::arg("file"),
       py::arg("init_function_name")="init",
R"raw_string(
Utility function to compile a c++ file with python bindings at run-time.

Parameters:

file: c++ code file (type: pathlib.Path)

init_function_name (default = "init"): Function, which is called after the compiled code is loaded. The prototype must match:
extern "C" void init_function_name(pybind11::object & res);
)raw_string"
    );



}

#endif
