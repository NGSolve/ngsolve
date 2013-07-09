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
  template NGS_DLL_HEADER class H1HighOrderFE<ET_POINT>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_SEGM>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_TRIG>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_QUAD>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_TET>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_PRISM>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_PYRAMID>;
  template NGS_DLL_HEADER class H1HighOrderFE<ET_HEX>;

  /*
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_POINT>, ET_POINT>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_SEGM>, ET_SEGM>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_TRIG>, ET_TRIG>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_QUAD>, ET_QUAD>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_TET>, ET_TET>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_PRISM>, ET_PRISM>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_HEX>, ET_HEX>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_PYRAMID>, ET_PYRAMID>;
  */
}
