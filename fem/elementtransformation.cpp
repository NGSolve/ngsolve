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
    int vnums[12];
    if (boundary)
      Ng_GetSurfaceElement (elnr+1, vnums);
    else
      Ng_GetElement (elnr+1, vnums);
    
    switch (eltype)
      {
      case ET_TRIG:
	for (int i = 0; i < 3; i++) sort[i] = i;
	if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
	if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	// vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]]
	break; 

      case ET_TET:
	for (int i = 0; i < 4; i++) sort[i] = i;
	if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	if (vnums[sort[2]] > vnums[sort[3]]) Swap (sort[2], sort[3]);
	if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
	if (vnums[sort[1]] > vnums[sort[3]]) Swap (sort[1], sort[3]);
	if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);

	// vnums[sort[0]] < vnums[sort[1]] < vnums[sort[2]] < vnums[sort[3]]
	break; 

      case ET_PRISM:
	for (int i = 0; i < 6; i++) sort[i] = i;

	if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
	if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
	if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
        
	if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);
	if (vnums[sort[4]] > vnums[sort[5]]) Swap (sort[4], sort[5]);
	if (vnums[sort[3]] > vnums[sort[4]]) Swap (sort[3], sort[4]);
	break;

      default:
	throw Exception ("undefined eltype in ElementTransformation::GetSort()\n");
      }
  }




}
#endif
