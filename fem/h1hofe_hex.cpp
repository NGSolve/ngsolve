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
  bool H1HighOrderFE_Shape<ET_HEX> :: GetDiagDualityMassInverse2 (FlatVector<> diag) const 
  {
    diag.Range(0,8) = 1.0;
    int ii = 8;
    for (int i = 0; i < N_EDGE; i++)
      for (int j = 2; j <= order_edge[i]; j++)
        diag(ii++) = (2*j-1)*(2*j)*(2*j-2);
    for (int f = 0; f < N_FACE; f++)
      {
        IVec<2> p = order_face[f];
        for (int i = 2; i <= p[0]; i++)
          for (int j = 2; j <= p[1]; j++)
            diag(ii++) = 1.0*(2*j-1)*(2*j)*(2*j-2) * (2*i-1)*(2*i)*(2*i-2);
      }
    IVec<3> p = order_cell[0];    
    for (int i = 2; i <= p[0]; i++)
      for (int j = 2; j <= p[1]; j++)
        for (int k = 2; k <= p[2]; k++)
          diag(ii++) = 1.0*(2*j-1)*(2*j)*(2*j-2) * (2*i-1)*(2*i)*(2*i-2) * (2*k-1)*(2*k)*(2*k-2);
    
    return true;
  }

  
  template class H1HighOrderFE<ET_HEX>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_HEX>, ET_HEX>;  
}
