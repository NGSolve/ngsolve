/*********************************************************************/
/* File:   l2hofe.cpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

#define FILE_L2HOFE_CPP

#include <fem.hpp>
#include "l2hofefo.hpp"



namespace ngfem
{


  template<>
  ScalarFiniteElement<1> * CreateL2HighOrderFE<ET_SEGM> (int order, FlatArray<int> vnums, LocalHeap & lh)
  {
    DGFiniteElement<1> * hofe = 0;
    switch (order)
      {
      case 0: hofe = new (lh)  L2HighOrderFEFO<ET_SEGM,0> (); break;
      case 1: hofe = new (lh)  L2HighOrderFEFO<ET_SEGM,1> (); break;
      case 2: hofe = new (lh)  L2HighOrderFEFO<ET_SEGM,2> (); break;
      case 3: hofe = new (lh)  L2HighOrderFEFO<ET_SEGM,3> (); break;
      case 4: hofe = new (lh)  L2HighOrderFEFO<ET_SEGM,4> (); break;
      case 5: hofe = new (lh)  L2HighOrderFEFO<ET_SEGM,5> (); break;
      case 6: hofe = new (lh)  L2HighOrderFEFO<ET_SEGM,6> (); break;
      case 7: hofe = new (lh)  L2HighOrderFEFO<ET_SEGM,7> (); break;      
      case 8: hofe = new (lh)  L2HighOrderFEFO<ET_SEGM,8> (); break;
      default: hofe = new (lh) L2HighOrderFE<ET_SEGM> (order); break;
      }

    for (int j = 0; j < 2; j++)
      hofe->SetVertexNumber (j, vnums[j]);
    
    return hofe;
  }
  


  template<>
  ScalarFiniteElement<2> * CreateL2HighOrderFE<ET_TRIG> (int order, FlatArray<int> vnums, LocalHeap & lh)
  {
    DGFiniteElement<2> * hofe = 0;
    switch (order)
      {
      case 0: hofe = new (lh)  L2HighOrderFEFO<ET_TRIG,0> (); break;
      case 1: hofe = new (lh)  L2HighOrderFEFO<ET_TRIG,1> (); break;
      case 2: hofe = new (lh)  L2HighOrderFEFO<ET_TRIG,2> (); break;
      case 3: hofe = new (lh)  L2HighOrderFEFO<ET_TRIG,3> (); break;
      default: hofe = new (lh) L2HighOrderFE<ET_TRIG> (order); break;
      }

    for (int j = 0; j < 3; j++)
      hofe->SetVertexNumber (j, vnums[j]);
    
    return hofe;
  }
  
  
}

