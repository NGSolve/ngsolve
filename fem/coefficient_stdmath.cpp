#include <core/register_archive.hpp>
// #include <fem.hpp>
#include <coefficient.hpp>
#include "scalarfe.hpp"
#include <../ngstd/evalfunc.hpp>
#include <algorithm>
#ifdef NGS_PYTHON
#include <core/python_ngcore.hpp> // for shallow archive
#include "../ngstd/python_ngstd.hpp"
#include "python_fem.hpp"
#endif // NGS_PYTHON

#include "coefficient_stdmath.hpp"


namespace ngfem
{

#ifdef NGS_PYTHON
  NGS_DLL_HEADER void ExportStdMathFunctions(py::module &m);
  
  void ExportStdMathFunctions(py::module &m)
  {
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
  }
#endif
  

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericSqrt>::Diff(const CoefficientFunction * var,
                                  shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    return CWMult (0.5/sqrt(c1), c1->Diff(var, dir));
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericSqrt>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return 0.5/sqrt(c1) * c1->DiffJacobi(var, cache);
  }

  
  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericSin>::Diff(const CoefficientFunction * var,
                                 shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    // return UnaryOpCF(c1, GenericCos(), "cos") * c1->Diff(var, dir);
    // return CWMult (UnaryOpCF(c1, GenericCos(), "cos"), c1->Diff(var, dir));
    return CWMult (cos(c1), c1->Diff(var, dir));
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericSin>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return cos(c1) * c1->DiffJacobi(var, cache);
  }
  
  
  
  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericCos>::Diff(const CoefficientFunction * var,
                                 shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    // return -1 * UnaryOpCF(c1, GenericSin(), "sin") * c1->Diff(var, dir);
    // return CWMult (-1 * UnaryOpCF(c1, GenericSin(), "sin"),  c1->Diff(var, dir));
    return CWMult (-sin(c1),  c1->Diff(var, dir));
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericCos>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return -sin(c1) * c1->DiffJacobi(var, cache);
  }

  
  
  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericTan>::Diff(const CoefficientFunction * var,
                                 shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    return 1.0 / (UnaryOpCF(c1, GenericCos(), "cos")*UnaryOpCF(c1, GenericCos(), "cos")) * c1->Diff(var, dir);
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericTan>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return 1.0 / (UnaryOpCF(c1, GenericCos(), "cos")*UnaryOpCF(c1, GenericCos(), "cos")) * c1->DiffJacobi(var, cache);
  }
  

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericASin>::Diff(const CoefficientFunction * var,
                                  shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    return make_shared<ConstantCoefficientFunction>(1)/UnaryOpCF(make_shared<ConstantCoefficientFunction>(1) - c1 * c1, GenericSqrt(), "sqrt") * c1->Diff(var, dir);
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericASin>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return make_shared<ConstantCoefficientFunction>(1)/UnaryOpCF(make_shared<ConstantCoefficientFunction>(1) - c1 * c1, GenericSqrt(), "sqrt") * c1->DiffJacobi(var, cache);
  }


  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericACos>::Diff(const CoefficientFunction * var,
                                  shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    return make_shared<ConstantCoefficientFunction>(-1)/UnaryOpCF(make_shared<ConstantCoefficientFunction>(1) - c1*c1, GenericSqrt(), "sqrt") * c1->Diff(var, dir);
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericACos>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return make_shared<ConstantCoefficientFunction>(-1)/UnaryOpCF(make_shared<ConstantCoefficientFunction>(1) - c1*c1, GenericSqrt(), "sqrt") * c1->DiffJacobi(var, cache);
  }


  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericATan>::Diff(const CoefficientFunction * var,
                                  shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    return make_shared<ConstantCoefficientFunction>(1) / (c1*c1 + make_shared<ConstantCoefficientFunction>(1)) * c1->Diff(var, dir);
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericATan>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return make_shared<ConstantCoefficientFunction>(1) / (c1*c1 + make_shared<ConstantCoefficientFunction>(1)) * c1->DiffJacobi(var,cache);
  }

  


  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericSinh>::Diff(const CoefficientFunction * var,
                                  shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    // return UnaryOpCF(c1, GenericCosh(), "cosh") * c1->Diff(var, dir);
    // return CWMult (UnaryOpCF(c1, GenericCosh(), "cosh"), c1->Diff(var, dir));
    return CWMult (cosh(c1), c1->Diff(var, dir));
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericSinh>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return cosh(c1) * c1->DiffJacobi(var, cache);
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericCosh>::Diff(const CoefficientFunction * var,
                                  shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    // return UnaryOpCF(c1, GenericSinh(), "sinh") * c1->Diff(var, dir);
    // return CWMult (UnaryOpCF(c1, GenericSinh(), "sinh"), c1->Diff(var, dir));
    return CWMult (sinh(c1), c1->Diff(var, dir));
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericCosh>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return sinh(c1) * c1->DiffJacobi(var, cache);
  }



  
  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericExp>::Diff(const CoefficientFunction * var,
                                 shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    // return UnaryOpCF(c1, GenericExp(), "exp") * c1->Diff(var, dir);
    return CWMult (UnaryOpCF(c1, GenericExp(), "exp"), c1->Diff(var, dir));
  }
  
  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericExp>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction> (1);
    return const_cast<cl_UnaryOpCF<GenericExp>*>(this)->shared_from_this() * c1->DiffJacobi(var, cache);
  }



  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericLog>::Diff(const CoefficientFunction * var,
                                 shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    return c1->Diff(var, dir) / c1;
  }
  
  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericLog>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return make_shared<ConstantCoefficientFunction>(1.0)/c1 * c1->DiffJacobi(var,cache);
  }


  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericErf>::Diff(const CoefficientFunction * var,
                                 shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    return CWMult (2. / sqrt(M_PI) * exp(- c1 * c1), c1->Diff(var, dir));
  }

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericErf>::DiffJacobi(const CoefficientFunction * var, T_DJC & cache) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return 2. / sqrt(M_PI) * exp(- c1 * c1) * c1->DiffJacobi(var,cache);
  }
  


  template <typename FUNC>
  shared_ptr<CoefficientFunction> MakeStdMathFunction(shared_ptr<CoefficientFunction> x)
  {
    // static RegisterClassForArchive<cl_UnaryOpCF<FUNC>, CoefficientFunction> reguopcf;
    
    FUNC func;
    return UnaryOpCF(x, func, FUNC::Name());
  }

  
  shared_ptr<CoefficientFunction> sqrt(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericSqrt>(x);
  }

  shared_ptr<CoefficientFunction> sin(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericSin>(x);
  }
  shared_ptr<CoefficientFunction> cos(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericCos>(x);
  }
  shared_ptr<CoefficientFunction> tan(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericTan>(x);
  }

  shared_ptr<CoefficientFunction> asin(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericASin>(x);
  }
  shared_ptr<CoefficientFunction> acos(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericACos>(x);
  }
  shared_ptr<CoefficientFunction> atan(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericATan>(x);
  }



  shared_ptr<CoefficientFunction> sinh(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericSinh>(x);
  }
  shared_ptr<CoefficientFunction> cosh(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericCosh>(x);
  }
  
  shared_ptr<CoefficientFunction> exp(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericExp>(x);
  }
  shared_ptr<CoefficientFunction> log(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericLog>(x);
  }

   shared_ptr<CoefficientFunction> erf(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericErf>(x);
  } 

  shared_ptr<CoefficientFunction> floor(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericFloor>(x);
  }
  shared_ptr<CoefficientFunction> ceil(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericCeil>(x);
  }

  
}
