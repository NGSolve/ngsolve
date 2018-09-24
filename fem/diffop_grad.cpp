#include <fem.hpp>
#include <diffop_impl.hpp>
  
namespace ngfem
{ 
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
