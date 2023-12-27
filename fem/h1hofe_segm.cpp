/*********************************************************************/
/* File:   h1hofe_segm.cpp                                           */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

 
// #include <fem.hpp>
#include <h1hofe.hpp>

#include <h1hofe_impl.hpp>
#include <tscalarfe_impl.hpp>


namespace ngfem
{
  template class H1HighOrderFE<ET_SEGM>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_SEGM>, ET_SEGM>;
}
