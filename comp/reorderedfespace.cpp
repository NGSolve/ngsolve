/*********************************************************************/
/* File:   reordered.cpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   Jun 2019                                                  */
/*********************************************************************/


#include <comp.hpp>

namespace ngcomp {
  
  ReorderedFESpace :: ReorderedFESpace (shared_ptr<FESpace> aspace, const Flags & flags)
    : FESpace(aspace->GetMeshAccess(), flags), space(aspace)
  {
    type = "Reordered" + space->type;
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
          low_order_space = make_shared<ReorderedFESpace>(space->LowOrderFESpacePtr(),lo_flags,used_idnrs);
        }
    */
    }
    
  void ReorderedFESpace :: Update()
  {      
    space->Update();
    FESpace::Update();

    SetNDof(space->GetNDof());
    size_t ndof = space->GetNDof();
    dofmap.SetSize(ndof);
    dofmap = UNUSED_DOF;

    Array<DofId> dofs;
    size_t cnt = 0;
    for (auto nt : { NT_VERTEX, NT_EDGE, NT_FACE, NT_CELL })
      for (auto nr : Range(ma->GetNNodes(nt)))
        {
          space->GetDofNrs (NodeId(nt, nr), dofs);
          for (auto d : dofs)
            dofmap[d] = cnt++;
        }
    
    ctofdof.SetSize(ndof);
    for (auto i : Range(ndof))
      ctofdof[dofmap[i]] = space->GetDofCouplingType(i);
  }
           
  FiniteElement& ReorderedFESpace :: GetFE (ElementId ei, Allocator & alloc) const
    {
      return space->GetFE(ei,alloc);
    }

  
  void ReorderedFESpace :: GetDofNrs(ElementId ei, Array<DofId> & dnums) const
  {
    space->GetDofNrs (ei, dnums);
    for (auto & d : dnums)
      d = dofmap[d];
  }

  void ReorderedFESpace :: GetDofNrs (NodeId ni, Array<DofId> & dnums) const
  {
    throw Exception ("ReorderedFESpace :: GetDofNrs(NodeId) not implemented");
  }
  
  void ReorderedFESpace :: GetVertexDofNrs (int vnr,  Array<DofId> & dnums) const
  {
    throw Exception ("ReorderedFESpace :: GetVertexDofNrs not implemented");    
  }
  
  void ReorderedFESpace :: GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const
  {
    throw Exception ("ReorderedFESpace :: GetEdgeDofNrs not implemented");    
  }
    
  void ReorderedFESpace :: GetFaceDofNrs (int fanr, Array<DofId> & dnums) const
  {
    throw Exception ("ReorderedFESpace :: GetFacetDofNrs not implemented");        
  }
}
