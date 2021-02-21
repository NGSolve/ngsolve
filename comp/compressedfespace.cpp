#include <comp.hpp>
#include "compressedfespace.hpp"

namespace ngcomp
{

  CompressedFESpace::CompressedFESpace (shared_ptr<FESpace> bfes)
      : FESpace (bfes->GetMeshAccess(), bfes->GetFlags()), space(bfes)
  {
    type = "wrapped-" + space->type;
    for (auto vb : { VOL,BND, BBND, BBBND })
    {
      evaluator[vb] = space->GetEvaluator(vb);
      flux_evaluator[vb] = space->GetFluxEvaluator(vb);
      integrator[vb] = space->GetIntegrator(vb);
    }
    iscomplex = space->IsComplex();
  }

  void CompressedFESpace::Update()
  { 
    // space->Update(); // removed as it may override changed doftypes
    // we need it to update the space after refinement.
    // space should set dofs only in first update on level

    const int ndofall = space->GetNDof();
    all2comp.SetSize(ndofall);
    comp2all.SetSize(ndofall);

    if (active_dofs && active_dofs->Size() != ndofall)
      throw Exception("active_dofs size doesn't match FESpace (anymore?).\n[active_dofs->Size() = "+to_string(active_dofs->Size())+", ndofall = "+to_string(ndofall)+"]");

    int ndof = 0;
    for (int i : Range(ndofall))
    {
      if (   (   active_dofs  && (active_dofs->Test(i))) 
          || ( (!active_dofs) && ((space->GetDofCouplingType(i) & VISIBLE_DOF))))
      {
        comp2all[ndof] = i;
        all2comp[i] = ndof++;
      }
      else
      {
        // all2comp[i] = -1;
        if (space->GetDofCouplingType(i) == HIDDEN_DOF)
          all2comp[i] = NO_DOF_NR_CONDENSE;
        else
          all2comp[i] = NO_DOF_NR;
      }
    }
    comp2all.SetSize(ndof);

    ctofdof.SetSize(ndof);
    for (int i : Range(ndof))
      ctofdof[i] = space->GetDofCouplingType(comp2all[i]);


    (*testout) << "dof mapping of the wrapper space:" << endl;
    for (int i : Range(ndof))
      (*testout) << i << " -> " << comp2all[i] << endl;

    SetNDof(ndof);
    FESpace::FinalizeUpdate();
  }

  FiniteElement & CompressedFESpace::GetFE (ElementId ei, Allocator & lh) const
  {
    return space->GetFE(ei,lh);
  }

  void CompressedFESpace::GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    space->GetDofNrs(ei,dnums);
    WrapDofs(dnums);
  }

  void CompressedFESpace::GetElementDofsOfType (ElementId ei, Array<DofId> & dnums, COUPLING_TYPE ctype) const
  {
    space->GetElementDofsOfType(ei,dnums,ctype);
  }

  void CompressedFESpace::GetDofNrs (NodeId ni, Array<DofId> & dnums) const
  {
    space->GetDofNrs(ni,dnums);
    WrapDofs(dnums);
  }

  void CompressedFESpace::GetVertexDofNrs (int vnr, Array<DofId> & dnums) const
  {
    space->GetVertexDofNrs(vnr,dnums);
    WrapDofs(dnums);
  }

  /// get dofs on edge enr
  void CompressedFESpace::GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const
  {
    space->GetEdgeDofNrs(ednr,dnums);
    WrapDofs(dnums);
  }

  /// get dofs on face fnr
  void CompressedFESpace::GetFaceDofNrs (int fanr, Array<DofId> & dnums) const
  {
    space->GetFaceDofNrs(fanr,dnums);
    WrapDofs(dnums);
  }

  /// get dofs on element (=cell) elnr
  void CompressedFESpace::GetInnerDofNrs (int elnr, Array<DofId> & dnums) const
  {
    space->GetInnerDofNrs(elnr,dnums);
    WrapDofs(dnums);
  }



}
