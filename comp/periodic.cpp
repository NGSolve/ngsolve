
/*********************************************************************/
/* File:   periodic.cpp                                              */
/* Author: Christopher Lackner                                       */
/* Date:   Feb. 2017                                                 */
/*********************************************************************/


#include <comp.hpp>

namespace ngcomp {

  PeriodicFESpace :: PeriodicFESpace (shared_ptr<FESpace> aspace, const Flags & flags)
    : FESpace(aspace->GetMeshAccess(), flags), space(aspace)
    {
      type = "Periodic" + space->type;
      for(auto vb : {VOL,BND,BBND})
        {
          evaluator[vb] = space->GetEvaluator(vb);
          flux_evaluator[vb] = space->GetFluxEvaluator(vb);
          integrator[vb] = space->GetIntegrator(vb);
        }
      iscomplex = space->IsComplex();
      // not yet working...
      if (space->LowOrderFESpacePtr() && false)
        {
          auto lo_flags = flags;
          lo_flags.SetFlag("order",1);
          low_order_space = make_shared<PeriodicFESpace>(space->LowOrderFESpacePtr(),lo_flags);
        }
    }
    
  void PeriodicFESpace :: Update (LocalHeap & lh)
    {      
      space->Update (lh);
      FESpace::Update(lh);
      dofmap.SetSize (space->GetNDof());
      for (int i = 0; i < dofmap.Size(); i++)
	dofmap[i] = i;

      int nid = ma->GetNetgenMesh()->GetIdentifications().GetMaxNr();

      // identifications are 1 based
      for (int idnr = 1; idnr<nid+1; idnr++)
        {
          Array<int> slave_dofnrs;
          Array<int> master_dofnrs;
          for (auto node_type : {NT_VERTEX, NT_EDGE, NT_FACE})
            {
              auto & periodic_nodes = ma->GetPeriodicNodes(node_type, idnr);
              for(auto node_pair : periodic_nodes)
                {
                  space->GetDofNrs(NodeId(node_type,node_pair[0]),master_dofnrs);
                  space->GetDofNrs(NodeId(node_type,node_pair[1]),slave_dofnrs);
                  for(auto i : Range(master_dofnrs.Size()))
                    {
                    dofmap[slave_dofnrs[i]] = dofmap[master_dofnrs[i]];
                    }
                }
            }
        }
      ctofdof.SetSize(dofmap.Size());
      for (auto i : Range(ctofdof.Size()))
        ctofdof[i] = space->GetDofCouplingType(i);
      for (int i = 0; i < dofmap.Size(); i++)
	if (dofmap[i] != i){
          ctofdof[i] = UNUSED_DOF;
	}
    }
    
  FiniteElement& PeriodicFESpace :: GetFE (ElementId ei, Allocator & alloc) const
    {
      auto & fe = space->GetFE(ei,alloc);
      const auto & ngel = ma->GetElement(ei);
      switch (ngel.GetType())
	{
	case ET_TRIG:
	  {
            auto hofe = dynamic_cast<VertexOrientedFE<ET_TRIG>*>(&fe);
            if(hofe)
              hofe->SetVertexNumbers(dofmap[ngel.Vertices()]);
            break;
	  }
	case ET_TET:
	  {
            auto hofe = dynamic_cast<VertexOrientedFE<ET_TET>*>(&fe);
            if(hofe)
              hofe->SetVertexNumbers(dofmap[ngel.Vertices()]);
            break;
	  }
        default:
          throw Exception("ElementType not implemented for PeriodicFESpace::GetFE");
	}
      return fe;
    }


  void PeriodicFESpace :: GetDofNrs(ElementId ei, Array<int> & dnums) const
    {
      space->GetDofNrs(ei,dnums);
      for (int i = 0; i< dnums.Size(); i++)
	dnums[i] = dofmap[dnums[i]];
    }

}
