/*********************************************************************/
/* File:   finiteelement.cpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


/* 
   Finite Element Definitions
*/

#include <bla.hpp>


#ifdef NETGEN_ELTRANS
#include <nginterface.h>
#include <nginterface_v2.hpp>
#endif


#include <fem.hpp>


#ifdef NETGEN_ELTRANS
namespace ngfem
{
  
  /*
  void ElementTransformation :: SetElement (bool aboundary, int aelnr, int aelindex)
  {
    boundary = aboundary;
    elnr = aelnr;
    elindex = aelindex;
    dim = Ng_GetDimension();
    if (boundary)
      iscurved = Ng_IsSurfaceElementCurved (elnr+1);
    else
      iscurved = Ng_IsElementCurved (elnr+1);
  }


  
  void ElementTransformation :: GetSort (FlatArray<int> sort) const
  {
  }
  */

}
#endif
