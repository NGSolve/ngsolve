/*********************************************************************/
/* File:   hcurlhofe_prism.cpp                                       */
/* Author: Sabine Zaglmayr                                           */
/* Date:   20. Maerz 2003                                            */
/*                                                                   */
/* AutoDiff - revision: J. Schoeberl, March 2009                     */
/*********************************************************************/

// #include <fem.hpp>
#include <thcurlfe_impl.hpp>
#include <hcurlhofe_impl.hpp>

namespace ngfem
{
  template class HCurlHighOrderFE<ET_PRISM>;
  template class
  T_HCurlHighOrderFiniteElement<ET_PRISM, HCurlHighOrderFE_Shape<ET_PRISM>>;
}
