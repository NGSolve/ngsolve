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
  ScalarFiniteElement<2> * CreateL2HighOrderFE<ET_QUAD> (int order, FlatArray<int> vnums, Allocator & lh)
  {
    DGFiniteElement<ET_QUAD> * hofe = new (lh) L2HighOrderFE<ET_QUAD> (order);
    for (int j = 0; j < 4; j++)
      hofe->SetVertexNumber (j, vnums[j]);
    return hofe;
  }

  
}

