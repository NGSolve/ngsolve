#include <fem.hpp>
#include <diffop_impl.hpp>
  
namespace ngfem
{


  template <int D, typename FEL>
  shared_ptr<CoefficientFunction> DiffOpGradient<D,FEL> ::
  DiffShape (shared_ptr<CoefficientFunction> proxy,
             shared_ptr<CoefficientFunction> dir)
  {
    auto dir_proxy = dynamic_pointer_cast<ProxyFunction> (dir);
    // return (-1)*TransposeCF(dir_proxy->Deriv()) * proxy;
    return (-1)*dir_proxy->Deriv() * proxy;
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
