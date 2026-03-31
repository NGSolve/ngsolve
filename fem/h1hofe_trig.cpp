/*********************************************************************/
/* File:   h1hofe_segm.cpp                                           */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include <h1hofe_impl.hpp>
#include <tscalarfe_impl.hpp>

namespace ngfem
{

  template <>
  bool H1HighOrderFE_Shape<ET_TRIG> :: GetDiagDualityMassInverse2 (FlatVector<> diag) const 
  {
    diag.Range(0,3) = 1.0;
    int ii = 3;
    for (int i = 0; i < N_EDGE; i++)
      for (int j = 2; j <= order_edge[i]; j++)
        diag(ii++) = (2*j-1)*(2*j)*(2*j-2);
    IVec<2> p = order_face[0];
    for (int i = 0; i <= p[0]-3; i++)
      for (int j = 0; j <= p[0]-i-3; j++)
        diag(ii++) = 0.5*(5+2*i+2*j)*(4+2*i+j)*(j+1) * (2*i+3)*(2*i+4) / (i+1);
    // cout << "trig duality diag = " << diag << endl;
    return true;
  }


  
  template class H1HighOrderFE<ET_TRIG>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_TRIG>, ET_TRIG>;
}
