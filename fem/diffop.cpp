/*********************************************************************/
/* File:   diffop.cpp                                                */
/* Author: Start                                                     */
/* Date:   24. Nov. 2009                                             */
/*********************************************************************/

 
#include <fem.hpp>
#include "diffop_impl.hpp"



namespace ngfem
{
  template class T_DifferentialOperator<DiffOpId<1> >;
  template class T_DifferentialOperator<DiffOpId<2> >;
  template class T_DifferentialOperator<DiffOpId<3> >;

  template class T_DifferentialOperator<DiffOpIdBoundary<1> >;
  template class T_DifferentialOperator<DiffOpIdBoundary<2> >;
  template class T_DifferentialOperator<DiffOpIdBoundary<3> >;


  template class T_DifferentialOperator<DiffOpGradient<1> >;
  template class T_DifferentialOperator<DiffOpGradient<2> >;
  template class T_DifferentialOperator<DiffOpGradient<3> >;

  template class T_DifferentialOperator<DiffOpIdEdge<2> >;
  template class T_DifferentialOperator<DiffOpIdEdge<3> >;

  template class T_DifferentialOperator<DiffOpCurlEdge<2> >;
  template class T_DifferentialOperator<DiffOpCurlEdge<3> >;

}
