/*********************************************************************/
/* File:   h1hofe_segm.cpp                                           */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include <h1hofe_impl.hpp>
#include <tscalarfe_impl.hpp>

namespace ngfem
{
  template class H1HighOrderFE<ET_TRIG>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_TRIG>, ET_TRIG>;
}
