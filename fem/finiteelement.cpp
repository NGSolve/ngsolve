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
  
  using namespace ngfem;


  /*
  NeedsUpdateException :: NeedsUpdateException ()
    : Exception ("FE needs update")
  { ; }
  */

  void FiniteElement :: GetInternalDofs (Array<int> & idofs) const
  {
    idofs.SetSize(0);
  }
  
  void FiniteElement :: GetDofs (Array<Dof> & dofs) const
  {
    throw Exception(string ("GetDofs not implemented for element ") + typeid(*this).name());
  }


  CompoundFiniteElement ::  CompoundFiniteElement (Array<const FiniteElement*> & afea)
    : FiniteElement (), fea(afea)
  {
    if (fea.Size() && fea[0])
      {
	int ii = 0;
	while ( ii < fea.Size() )
	  {
	    // do not use dummy elements
	    if ( dynamic_cast<const FE_SegmDummy *> (fea[ii]) )
	      { ii++; continue; }
	    dimspace = fea[ii]->SpatialDim();
	    eltype = fea[ii]->ElementType();
	    break;
	  }
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
	dimspace = 0;
	eltype = ET_TRIG;
      }
  }

  void CompoundFiniteElement :: GetInternalDofs (Array<int> & idofs) const
  {
    *testout << "compound, getinternal dofs" << endl;
    idofs.SetSize (0);
    ArrayMem<int,20> bidofs;
    int base = 0;
    for (int i = 0; i < fea.Size(); i++)
      {
	(*this)[i].GetInternalDofs (bidofs);
	for (int j = 0; j < bidofs.Size(); j++)
	  idofs.Append (base+bidofs[j]);
	base += (*this)[i].GetNDof();
      }
    *testout << idofs << endl;
  }






}
