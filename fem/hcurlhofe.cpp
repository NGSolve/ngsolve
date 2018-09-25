/*********************************************************************/
/* File:   hcurlhofe.cpp                                             */
/* Author: Sabine Zaglmayr                                           */
/* Date:   20. Maerz 2003                                            */
/*                                                                   */
/* AutoDiff - revision: J. Schoeberl, March 2009                     */
/*********************************************************************/

#define FILE_HCURLHOFE_CPP
#include <fem.hpp>    

namespace ngfem
{
  // template class HCurlHighOrderFE<ET_POINT>;
  // template class HCurlHighOrderFE<ET_SEGM>;
  template class HCurlHighOrderFE<ET_TRIG>;
  // template class HCurlHighOrderFE<ET_QUAD>;

  // template class
  // T_HCurlHighOrderFiniteElement<ET_POINT, HCurlHighOrderFE_Shape<ET_POINT>>;
  // template class
  // T_HCurlHighOrderFiniteElement<ET_SEGM, HCurlHighOrderFE_Shape<ET_SEGM>>;
  template class
  T_HCurlHighOrderFiniteElement<ET_TRIG, HCurlHighOrderFE_Shape<ET_TRIG>>;
  // template class
  // T_HCurlHighOrderFiniteElement<ET_QUAD, HCurlHighOrderFE_Shape<ET_QUAD>>;

}
