/*********************************************************************/
/* File:   l2hofe_trig.cpp                                           */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

// #define FILE_L2HOFE_TRIG_CPP
 
#include <tscalarfe_impl.hpp>
#include <l2hofe_impl.hpp>
#include "l2hofefo.hpp"

namespace ngfem
{

  /*
  template <> inline void L2HighOrderFE<ET_TRIG> ::
  GetDiagMassMatrix(FlatVector<> mass) const
  {
    for (int ix = 0, ii = 0; ix <= order; ix++)
      for (int iy = 0; iy <= order - ix; iy++, ii++)
	mass(ii) = 1.0 / ((2 * ix + 1) * (2 * ix + 2 * iy + 2));
  }
  */

  
  template class L2HighOrderFE<ET_TRIG>;  
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_TRIG>, ET_TRIG, DGFiniteElement<ET_TRIG> >;

  //constexpr int MAX_FO_TRIG = 2;
  
  template<>
  ScalarFiniteElement<2> * CreateL2HighOrderFE<ET_TRIG> (int order, FlatArray<int> vnums, Allocator & lh)
  {
    DGFiniteElement<ET_TRIG> * hofe = nullptr;

    // fixed order only for the orientations of an oriented mesh: first vertex is the lowest
    if (vnums[0] < vnums[1] && vnums[0] < vnums[2])
      switch (order)
        {
        case 0: hofe = new (lh)  L2HighOrderFEFO<ET_TRIG,0, MeshOrientation> (); break;
        case 1: hofe = new (lh)  L2HighOrderFEFO<ET_TRIG,1, MeshOrientation> (); break;
        case 2: hofe = new (lh)  L2HighOrderFEFO<ET_TRIG,2, MeshOrientation> (); break;
        default: ;
        }

    if (!hofe)
      hofe = new (lh) L2HighOrderFE<ET_TRIG> (order);


    for (int j = 0; j < 3; j++)
      hofe->SetVertexNumber (j, vnums[j]);
    
    return hofe;
  }


  
}
