/*********************************************************************/
/* File:   h1hofe.cpp                                                */
/* Author: Start                                                      */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/
 
#include <h1hofe_impl.hpp>
#include <tscalarfe_impl.hpp>

namespace ngfem
{


  template <>
  bool H1HighOrderFE_Shape<ET_TET> :: GetDiagDualityMassInverse2 (FlatVector<> diag) const 
  {
    diag.Range(0,4) = 1.0;
    int ii = 4;
    for (int i = 0; i < N_EDGE; i++)
      for (int j = 2; j <= order_edge[i]; j++)
        diag(ii++) = (2*j-1)*(2*j)*(2*j-2);
    for (int f = 0; f < N_FACE; f++)
      if (int p = order_face[f][0]; p >= 3)
	ii += DubinerBasisOrthoBub::CalcNormInv (p-3, diag+ii);
    if (int p = order_cell[0][0]; p >= 4)
      DubinerBasis3DOrthoBub::CalcNormInv(p-4, diag+ii);
    return true;
  }


  
  template class H1HighOrderFE<ET_TET>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_TET>, ET_TET>;  
}
