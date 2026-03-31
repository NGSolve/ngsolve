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
  bool H1HighOrderFE_Shape<ET_QUAD> :: GetDiagDualityMassInverse2 (FlatVector<> diag) const 
  {
    diag.Range(0,4) = 1.0;
    int ii = 4;
    for (int i = 0; i < N_EDGE; i++)
      for (int j = 2; j <= order_edge[i]; j++)
        diag(ii++) = (2*j-1)*(2*j)*(2*j-2);
    IVec<2> p = order_face[0];
    for (int i = 2; i <= p[0]; i++)
      for (int j = 2; j <= p[1]; j++)
        diag(ii++) = 1.0*(2*j-1)*(2*j)*(2*j-2) * (2*i-1)*(2*i)*(2*i-2);

    // cout << "quad duality diag = " << diag << endl;
    return true;
  }


  
  template class H1HighOrderFE<ET_QUAD>;
  template class T_ScalarFiniteElement<H1HighOrderFE_Shape<ET_QUAD>, ET_QUAD>;
}
