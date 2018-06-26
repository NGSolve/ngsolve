/*********************************************************************/
/* File:   l2hofe_trig.cpp                                           */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

#define FILE_L2HOFE_CPP

#include <fem.hpp>
#include <tscalarfe_impl.hpp>
#include <l2hofe_impl.hpp>
#include "l2hofefo.hpp"


namespace ngfem
{
  template class L2HighOrderFE<ET_TET>;  
  template class T_ScalarFiniteElement<L2HighOrderFE_Shape<ET_TET>, ET_TET, DGFiniteElement<3> >;



  template<>
  ScalarFiniteElement<3> * CreateL2HighOrderFE<ET_TET> (int order, FlatArray<int> vnums, Allocator & lh)
  {
    DGFiniteElement<3> * hofe = nullptr;    
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
