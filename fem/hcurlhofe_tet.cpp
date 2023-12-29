/*********************************************************************/
/* File:   hcurlhofe_tet.cpp                                         */
/* Author: Sabine Zaglmayr                                           */
/* Date:   20. Maerz 2003                                            */
/*                                                                   */
/* AutoDiff - revision: J. Schoeberl, March 2009                     */
/*********************************************************************/


#include <thcurlfe_impl.hpp>
#include <hcurlhofe_impl.hpp>

namespace ngfem
{
  template class HCurlHighOrderFE<ET_TET>;
  template class
  T_HCurlHighOrderFiniteElement<ET_TET, HCurlHighOrderFE_Shape<ET_TET>>;
}
