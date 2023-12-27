/*********************************************************************/
/* File:   l2hofe_trig.cpp                                           */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

// #define FILE_L2HOFE_TET_CPP

// #include <fem.hpp>
#include "l2hofe.hpp"
#include <tscalarfe_impl.hpp>
#include <l2hofe_impl.hpp>
#include "l2hofefo.hpp"


namespace ngfem
{
  /*
  template <> inline void L2HighOrderFE<ET_TET> ::
  GetDiagMassMatrix(FlatVector<> mass) const
  {
    for (int ix = 0, ii = 0; ix <= order; ix++)
      for (int iy = 0; iy <= order - ix; iy++)
        for (int iz = 0; iz <= order - ix-iy; iz++, ii++)
          mass(ii) = 1.0 / ((2 * ix + 1) * (2 * ix + 2 * iy + 2) * (2 * ix + 2 * iy + 2 * iz + 3));
  }
  */
  
  template class L2HighOrderFE<ET_TET>;  
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_TET>, ET_TET, DGFiniteElement<ET_TET> >;



  template<>
  ScalarFiniteElement<3> * CreateL2HighOrderFE<ET_TET> (int order, FlatArray<int> vnums, Allocator & lh)
  {
    DGFiniteElement<ET_TET> * hofe = nullptr;    
    if (vnums[0] < vnums[1] && vnums[1] < vnums[2] && vnums[1] < vnums[3])
    // if (false)
      { // new standard orientation
        if (vnums[2] < vnums[3])
          {
            switch (order)
              {
              case 0: hofe = new (lh)  L2HighOrderFEFO<ET_TET,0, FixedOrientation<0,1,2,3>> (); break;
              case 1: hofe = new (lh)  L2HighOrderFEFO<ET_TET,1, FixedOrientation<0,1,2,3>> (); break;
              case 2: hofe = new (lh)  L2HighOrderFEFO<ET_TET,2, FixedOrientation<0,1,2,3>> (); break;
                // case 3: return new (lh)  L2HighOrderFEFO<ET_TET,3, FixedOrientation<0,1,2,3>> (); break;
              default: ; 
              }
          }
        else
          {
            switch (order)
              {
              case 0: hofe = new (lh)  L2HighOrderFEFO<ET_TET,0, FixedOrientation<0,1,3,2>> (); break;
              case 1: hofe = new (lh)  L2HighOrderFEFO<ET_TET,1, FixedOrientation<0,1,3,2>> (); break;
              case 2: hofe = new (lh)  L2HighOrderFEFO<ET_TET,2, FixedOrientation<0,1,3,2>> (); break;
                // case 3: return new (lh)  L2HighOrderFEFO<ET_TET,3, FixedOrientation<0,1,3,2>> (); break;                
              default: ; 
              }
          }
      }

    if (!hofe)
      hofe = new (lh) L2HighOrderFE<ET_TET> (order); 
    for (int j = 0; j < 4; j++)
      hofe->SetVertexNumber (j, vnums[j]);
    return hofe;
  }

  
}
