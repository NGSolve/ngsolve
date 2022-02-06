#include <core/register_archive.hpp>
#include <fem.hpp>
#include <../ngstd/evalfunc.hpp>
#include <algorithm>
#ifdef NGS_PYTHON
#include <core/python_ngcore.hpp> // for shallow archive
#endif // NGS_PYTHON

#include "coefficient_stdmath.hpp"


namespace ngfem
{

  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericSin>::Diff(const CoefficientFunction * var,
                                 shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    // return UnaryOpCF(c1, GenericCos(), "cos") * c1->Diff(var, dir);
    // return CWMult (UnaryOpCF(c1, GenericCos(), "cos"), c1->Diff(var, dir));
    return CWMult (CosCF(c1), c1->Diff(var, dir));
  }
  


  
  
  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericCos>::Diff(const CoefficientFunction * var,
                                 shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    // return -1 * UnaryOpCF(c1, GenericSin(), "sin") * c1->Diff(var, dir);
    // return CWMult (-1 * UnaryOpCF(c1, GenericSin(), "sin"),  c1->Diff(var, dir));
    return CWMult (-SinCF(c1),  c1->Diff(var, dir));
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
  cl_UnaryOpCF<GenericExp>::DiffJacobi(const CoefficientFunction * var) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction> (1);
    return const_cast<cl_UnaryOpCF<GenericExp>*>(this)->shared_from_this() * c1->DiffJacobi(var);
  }



  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericLog>::Diff(const CoefficientFunction * var,
                                 shared_ptr<CoefficientFunction> dir) const
  {
    if (this == var) return dir;
    return c1->Diff(var, dir) / c1;
  }
  
  template <> shared_ptr<CoefficientFunction>
  cl_UnaryOpCF<GenericLog>::DiffJacobi(const CoefficientFunction * var) const
  {
    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    return make_shared<ConstantCoefficientFunction>(1.0)/c1 * c1->DiffJacobi(var);
  }





  


  template <typename FUNC>
  shared_ptr<CoefficientFunction> MakeStdMathFunction(shared_ptr<CoefficientFunction> x)
  {
    static RegisterClassForArchive<cl_UnaryOpCF<FUNC>, CoefficientFunction> reguopcf;
    
    FUNC func;
    return UnaryOpCF(x, func, FUNC::Name());
  }

  
  
  shared_ptr<CoefficientFunction> ExpCF(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericExp>(x);
  }
  shared_ptr<CoefficientFunction> SinCF(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericSin>(x);
  }
  shared_ptr<CoefficientFunction> CosCF(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericCos>(x);
  }
  shared_ptr<CoefficientFunction> LogCF(shared_ptr<CoefficientFunction> x)
  {
    return MakeStdMathFunction<GenericLog>(x);
  }

  
}
