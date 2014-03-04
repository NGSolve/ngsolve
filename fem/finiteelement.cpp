/*********************************************************************/
/* File:   finiteelement.cpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

#define FILE_FINITEELEMENT_CPP

/* 
   Finite Element Definitions
*/


#include <fem.hpp>

namespace ngfem
{
  
  string FiniteElement :: ClassName() const
  {
    return "FiniteElement"; 
  }

  void FiniteElement :: Print (ostream & ost) const
  {
    ost << "Finite Element, type = " << ElementType() << ", ndof = " << ndof << endl;
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
      }
  }

  void CompoundFiniteElement :: Print (ostream & ost) const
  {
    ost << "CompoundFiniteElement" << endl;
    for (int i = 0; i < GetNComponents(); i++)
      (*this)[i].Print (ost);
  }
}

