
/*********************************************************************/
/* File:   discontinuous.cpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   Jun 2019                                                  */
/*********************************************************************/


#include <comp.hpp>

namespace ngcomp {
  
  DiscontinuousFESpace :: DiscontinuousFESpace (shared_ptr<FESpace> aspace, const Flags & flags)
    : FESpace(aspace->GetMeshAccess(), flags), space(aspace)
  {
    type = "Discontinuous" + space->type;
    evaluator[VOL] = space->GetEvaluator(VOL);
    flux_evaluator[VOL] = space->GetFluxEvaluator(VOL);
    integrator[VOL] = space->GetIntegrator(VOL);
    
    iscomplex = space->IsComplex();
    /*
      // not yet implemented ...
      if (space->LowOrderFESpacePtr() && false)
        {
          auto lo_flags = flags;
          lo_flags.SetFlag("order",1);
          low_order_space = make_shared<DiscontinuousFESpace>(space->LowOrderFESpacePtr(),lo_flags,used_idnrs);
        }
    */
    }
    
  void DiscontinuousFESpace :: Update (LocalHeap & lh)
  {      
    space->Update (lh);
    FESpace::Update(lh);
    
    first_element_dof.SetSize(ma->GetNE()+1);
    for (ElementId ei : ma->Elements(VOL))
      {
        HeapReset hr(lh);
        first_element_dof[ei.Nr()] = space->GetFE(ei, lh).GetNDof();
      }

    size_t ndof = 0;
    for (size_t i : Range(ma->GetNE()))
      {
        size_t nd = first_element_dof[i];
        first_element_dof[i] = ndof;
        ndof += nd;
      }
    first_element_dof[ma->GetNE()] = ndof;
    SetNDof(ndof);
  }
           
  FiniteElement& DiscontinuousFESpace :: GetFE (ElementId ei, Allocator & alloc) const
    {
      if (ei.IsVolume())
        return space->GetFE(ei,alloc);

      ELEMENT_TYPE et = ma->GetElType(ei);
      return SwitchET(et, [&alloc] (auto type) -> FiniteElement&
                      { return * new (alloc) DummyFE<type.ElementType()>(); });
    }

  
  void DiscontinuousFESpace :: GetDofNrs(ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize0();
    if (ei.IsVolume())
      dnums += IntRange(first_element_dof[ei.Nr()], first_element_dof[ei.Nr()+1]);
  }

  void DiscontinuousFESpace :: GetDofNrs (NodeId ni, Array<DofId> & dnums) const
  {
    throw Exception ("DiscontinuousFESpace :: GetDofNrs(NodeId) not implemented");
  }
  
  void DiscontinuousFESpace :: GetVertexDofNrs (int vnr,  Array<DofId> & dnums) const
  {
    throw Exception ("DiscontinuousFESpace :: GetVertexDofNrs not implemented");    
  }
  
  void DiscontinuousFESpace :: GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const
  {
    throw Exception ("DiscontinuousFESpace :: GetEdgeDofNrs not implemented");    
  }
    
  void DiscontinuousFESpace :: GetFaceDofNrs (int fanr, Array<DofId> & dnums) const
  {
    throw Exception ("DiscontinuousFESpace :: GetFacetDofNrs not implemented");        
  }
}
