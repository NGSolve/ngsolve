// #include <fem.hpp>

#include "bdbequations.hpp"
#include "diffop_impl.hpp"
  
namespace ngfem
{ 
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<1> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesse<3> >;

  // template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<1> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<2> >;
  template class NGS_DLL_HEADER T_DifferentialOperator<DiffOpHesseBoundary<3> >;
}
