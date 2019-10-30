#include <fem.hpp>
#include <diffop_impl.hpp>
  
namespace ngfem
{


  template <int D, typename FEL>
  shared_ptr<CoefficientFunction> DiffOpGradient<D,FEL> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir)
  {
    return -TransposeCF(dir->Operator("Grad")) * proxy;
    // return -dir->Operator("grad") * proxy;
  }

  template <int D, typename FEL>
  shared_ptr<CoefficientFunction> DiffOpGradientBoundary<D,FEL> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir)
  {
    return -TransposeCF(dir->Operator("Gradboundary")) * proxy;
    // return -dir->Operator("grad") * proxy;
  }



  
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradient<1> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradient<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradient<3> >;

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientBoundary<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradientBoundary<3> >;

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradVectorH1<1> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradVectorH1<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradVectorH1<3> >;

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradBoundaryVectorH1<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpGradBoundaryVectorH1<3> >;
}
