
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

    DefineDefineFlag("BND");
    vb = flags.GetDefineFlag("BND") ? BND : VOL; 
    
    type = "Discontinuous" + space->type;

    for (VorB vorb : {VOL,BND,BBND})
    {
      evaluator[vorb] = space->GetEvaluator(vorb);
      flux_evaluator[vorb] = space->GetFluxEvaluator(vorb);
      integrator[vorb] = space->GetIntegrator(vorb);
    }    
    iscomplex = space->IsComplex();
    for(auto vb : Range(4))
      definedon[vb] = space->DefinedOn(VorB(vb));
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
    
  void DiscontinuousFESpace :: Update()
  {      
    space->Update();
    FESpace::Update();
    
    first_element_dof.SetSize(ma->GetNE(vb)+1);
    LocalHeap lh(10000, "discontinuous lh");
    // fes->Elements only iterates over definedon elements, here we need all!
    for (auto i : Range(ma->GetNE(vb))) {
      HeapReset hr(lh);
      first_element_dof[i] = space->GetFE(ElementId(vb, i), lh).GetNDof();
    }

    size_t ndof = 0;
    for (size_t i : Range(ma->GetNE(vb)))
      {
        size_t nd = first_element_dof[i];
        first_element_dof[i] = ndof;
        ndof += nd;
      }
    first_element_dof[ma->GetNE(vb)] = ndof;
    SetNDof(ndof);

    ctofdof.SetSize(ndof);
    ctofdof = LOCAL_DOF;
  }
           
  FiniteElement& DiscontinuousFESpace :: GetFE (ElementId ei, Allocator & alloc) const
    {
      if (ei.VB() == vb)
        return space->GetFE(ei,alloc);

      ELEMENT_TYPE et = ma->GetElType(ei);
      return SwitchET(et, [&alloc] (auto type) -> FiniteElement&
                      { return * new (alloc) DummyFE<type.ElementType()>(); });
    }

  
  void DiscontinuousFESpace :: GetDofNrs(ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize0();
    if (ei.VB() == vb)
      dnums += IntRange(first_element_dof[ei.Nr()], first_element_dof[ei.Nr()+1]);
  }

  void DiscontinuousFESpace :: GetDofNrs (NodeId ni, Array<DofId> & dnums) const
  {
    dnums.SetSize0();
    if (CoDimension(ni.GetType(), ma->GetDimension()) == int(vb))
      dnums += IntRange(first_element_dof[ni.GetNr()], first_element_dof[ni.GetNr()+1]);
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
