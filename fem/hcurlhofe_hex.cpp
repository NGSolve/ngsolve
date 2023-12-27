/*********************************************************************/
/* File:   hcurlhofe_hex.cpp                                         */
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
  template class HCurlHighOrderFE<ET_HEX>;
  template class
  T_HCurlHighOrderFiniteElement<ET_HEX, HCurlHighOrderFE_Shape<ET_HEX>>;
}
