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
           

  void ReorderedFESpace :: FinalizeUpdate ()
  {
    space->FinalizeUpdate();
    FESpace::FinalizeUpdate();

    /**
       If the underlying space is a CompoundFESpace, the CompoundFESpace's flags do
       usually not contain information about dirichlet-boundaries, so we have to
       set free_dofs manually.
     **/
    if (auto comp_fes = dynamic_pointer_cast<CompoundFESpace>(space)) {
      auto space_free_dofs = comp_fes->GetFreeDofs();
      free_dofs->Clear();
      auto external_space_free_dofs = comp_fes->GetFreeDofs(true);
      external_free_dofs->Clear();

      for (auto k : Range(GetNDof())) {
	if (space_free_dofs->Test(k))
	  { free_dofs->Set(dofmap[k]); }
	if (external_space_free_dofs->Test(k))
	  { external_free_dofs->Set(dofmap[k]); }
      }
    }
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
    space->GetDofNrs (ni, dnums);
    for (auto & d : dnums)
      d = dofmap[d];
  }
  
  void ReorderedFESpace :: GetVertexDofNrs (int vnr,  Array<DofId> & dnums) const
  {
    space->GetVertexDofNrs (vnr, dnums);
    for (auto & d : dnums)
      d = dofmap[d];
  }
  
  void ReorderedFESpace :: GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const
  {
    space->GetEdgeDofNrs (ednr, dnums);
    for (auto & d : dnums)
      d = dofmap[d];
  }
    
  void ReorderedFESpace :: GetFaceDofNrs (int fanr, Array<DofId> & dnums) const
  {
    space->GetFaceDofNrs (fanr, dnums);
    for (auto & d : dnums)
      d = dofmap[d];
  }
}
