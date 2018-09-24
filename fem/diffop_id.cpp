#include <fem.hpp>
#include <diffop_impl.hpp>
  
namespace ngfem
{ 
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpId<1>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpId<2>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpId<3>>;

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdBoundary<1>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdBoundary<2>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdBoundary<3>>;

  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<1,VOL>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<2,VOL>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<3,VOL>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<1,BND>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<2,BND>>;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpIdVectorH1<3,BND>>;
}
