
/*********************************************************************/
/* File:   hidden.cpp                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   Feb 2021                                                  */
/*********************************************************************/


#include <comp.hpp>

namespace ngcomp {
  
  HiddenFESpace :: HiddenFESpace (shared_ptr<FESpace> aspace, const Flags & flags)
    : FESpace(aspace->GetMeshAccess(), flags), space(aspace)
  {
    type = "Hidden" + space->type;

    for (VorB vorb : {VOL,BND,BBND})
      {
        evaluator[vorb] = space->GetEvaluator(vorb);
        flux_evaluator[vorb] = space->GetFluxEvaluator(vorb);
        integrator[vorb] = space->GetIntegrator(vorb);
      }    
    iscomplex = space->IsComplex();
    /*
      // not yet implemented ...
      if (space->LowOrderFESpacePtr() && false)
        {
        auto lo_flags = flags;
          lo_flags.SetFlag("order",1);
          low_order_space = make_shared<HiddenFESpace>(space->LowOrderFESpacePtr(),lo_flags,used_idnrs);
        }
    */
  }
    
  void HiddenFESpace :: Update()
  {      
    space->Update();
    FESpace::Update();
    
    SetNDof (0);
  }
           
  FiniteElement& HiddenFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    return space->GetFE(ei,alloc);
  }

  
  void HiddenFESpace :: GetDofNrs(ElementId ei, Array<DofId> & dnums) const
  {
    space->GetDofNrs(ei, dnums);
    dnums = -2; 
  }

  void HiddenFESpace :: GetDofNrs (NodeId ni, Array<DofId> & dnums) const
  {
    space->GetDofNrs(ni, dnums);
    dnums = -2; 
  }
  
  void HiddenFESpace :: GetVertexDofNrs (int vnr,  Array<DofId> & dnums) const
  {
    throw Exception ("HiddenFESpace :: GetVertexDofNrs not implemented");    
  }
  
  void HiddenFESpace :: GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const
  {
    throw Exception ("HiddenFESpace :: GetEdgeDofNrs not implemented");    
  }
    
  void HiddenFESpace :: GetFaceDofNrs (int fanr, Array<DofId> & dnums) const
  {
    throw Exception ("HiddenFESpace :: GetFacetDofNrs not implemented");        
  }
}
