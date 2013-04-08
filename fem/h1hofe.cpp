/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

 
#include <fem.hpp>

#include <h1hofe_impl.hpp>
#include <tscalarfe_impl.hpp>


namespace ngfem
{

  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_POINT>, ET_POINT>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_SEGM>, ET_SEGM>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_TRIG>, ET_TRIG>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_QUAD>, ET_QUAD>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_TET>, ET_TET>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_PRISM>, ET_PRISM>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_HEX>, ET_HEX>;
  template class T_ScalarFiniteElement2<H1HighOrderFE_Shape<ET_PYRAMID>, ET_PYRAMID>;

}
