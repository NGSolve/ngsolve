/*********************************************************************/
/* File:   l2hofe.cpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include <fem.hpp>
#include "l2hofe.hpp"

#include <l2hofe_impl.hpp>
#include <tscalarfe_impl.hpp>


namespace ngfem
{



  template <> inline void L2HighOrderFE<ET_POINT> :: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    mass(0) = 1;
  }

  template <> inline void L2HighOrderFE<ET_SEGM> :: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    for (int ix = 0; ix <= order; ix++)
      mass(ix) = 1.0 / (2*ix+1);
  }

  template <> inline void L2HighOrderFE<ET_TRIG> :: 
  GetDiagMassMatrix (FlatVector<> mass) const
  {
    for (int ix = 0, ii = 0; ix <= order; ix++)
      for (int iy = 0; iy <= order-ix; iy++, ii++)
        mass(ii) = 1.0 / ( (2*iy+1) * (2*ix+2*iy+2));
  }





  template NGS_DLL_HEADER class L2HighOrderFE<ET_POINT>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_SEGM>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_TRIG>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_QUAD>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_TET>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_PRISM>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_PYRAMID>;
  template NGS_DLL_HEADER class L2HighOrderFE<ET_HEX>;
 
  /*
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_POINT>, ET_POINT, DGFiniteElement<0> >;
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_SEGM>, ET_SEGM, DGFiniteElement<1> >;
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_TRIG>, ET_TRIG, DGFiniteElement<2> >;
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_QUAD>, ET_QUAD, DGFiniteElement<2> >;
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_TET>, ET_TET, DGFiniteElement<3> >;
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_PRISM>, ET_PRISM, DGFiniteElement<3> >;
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_HEX>, ET_HEX, DGFiniteElement<3> >;
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_PYRAMID>, ET_PYRAMID, DGFiniteElement<3> >;
  */

  
} // namespace







