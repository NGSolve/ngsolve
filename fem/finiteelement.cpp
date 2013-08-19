/*********************************************************************/
/* File:   finiteelement.cpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


/* 
   Finite Element Definitions
*/


#include <fem.hpp>

namespace ngfem
{
  
  // FiniteElement :: ~FiniteElement () { ; }


  string FiniteElement :: ClassName() const
  {
    return "FiniteElement"; 
  }

  void FiniteElement :: 
  PrecomputeShapes (const IntegrationRule & ir) 
  {
    ;
  }


  CompoundFiniteElement ::  CompoundFiniteElement (FlatArray<const FiniteElement*> afea)
    : FiniteElement (), fea(afea)
  {
    if (fea.Size() && fea[0])
      {
	// eltype = fea[0]->ElementType();
	ndof = 0;
	order = 0;
	for (int i = 0; i < fea.Size(); i++)
	  if (fea[i])
	    {
	      ndof += fea[i]->GetNDof();
	      order = max2 (order, fea[i]->Order());
	    }
	  else
	    {
	      cout << "WARNING: CompoundFE, undefined component" << i << endl;
	    }
      }
    else
      {
	throw Exception("WARNING: CompoundFE, undefined components");
	cout << "WARNING: CompoundFE, undefined components" << endl;
	ndof = 0;
	order = 0;
	// eltype = ET_TRIG;
      }
  }


  template class DummyFE<ET_SEGM>;
  template class DummyFE<ET_TRIG>;
  template class DummyFE<ET_QUAD>;
  template class DummyFE<ET_TET>;
  template class DummyFE<ET_PRISM>;
  template class DummyFE<ET_PYRAMID>;
  template class DummyFE<ET_HEX>;
}
