
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
          Array<std::tuple<NodeId,NodeId>> periodic_nodes;
          GetPeriodicNodeIds(periodic_nodes,idnr);

          Array<int> slave_dofnrs;
          Array<int> master_dofnrs;
          for(auto node_pair : periodic_nodes)
            {
              NodeId id1, id2;
              std::tie(id1,id2) = node_pair;
              slave_dofnrs.SetSize(0);
              master_dofnrs.SetSize(0);
              space->GetDofNrs(id1,master_dofnrs);
              space->GetDofNrs(id2,slave_dofnrs);
              for(auto i : Range(master_dofnrs.Size())){
                dofmap[slave_dofnrs[i]] = dofmap[master_dofnrs[i]];
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


  void PeriodicFESpace :: GetPeriodicNodeIds(Array<std::tuple<NodeId,NodeId>> & node_ids,int idnr) const
    {
      Array<INT<2>> vertex_pairs;
      ma->GetPeriodicVertices(idnr,vertex_pairs);
      for (auto pair : vertex_pairs)
        node_ids.Append(std::make_tuple(NodeId(NT_VERTEX,pair[0]),NodeId(NT_VERTEX,pair[1])));

      Array<int> vertex_map(GetNDof());
      for (auto i : Range(vertex_map.Size()))
        vertex_map[i] = i;
      for (auto pair : vertex_pairs)
        vertex_map[pair[1]] = pair[0];
      
      // build vertex-pair to edge hashtable:
      HashTable<INT<2>, int> vp2e(ma->GetNEdges());
      
      for (int enr = 0; enr < ma->GetNEdges(); enr++)
        {
          int v1, v2;
          ma->GetEdgePNums (enr, v1, v2);
          if (v1 > v2) Swap (v1, v2);
          vp2e[INT<2>(v1,v2)] = enr;
        }

      for (int enr = 0; enr < ma->GetNEdges(); enr++)
        {
          int v1, v2;
          ma->GetEdgePNums (enr, v1, v2);
          int mv1 = vertex_map[v1];
          int mv2 = vertex_map[v2];
          if(mv1 != v1 && mv2 != v2)
            {
              if (mv1 > mv2) Swap(mv1,mv2);
              int menr = vp2e.Get(INT<2>(mv1,mv2));
              node_ids.Append(std::make_tuple(NodeId(NT_EDGE,menr),NodeId(NT_EDGE,enr)));
            }
        }

      // build vertex-triple to face hashtable
      HashTable<INT<3>, int> v2f(ma->GetNFaces());
      Array<int> pnums;
      for (auto fnr : Range(ma->GetNFaces()))
        {
          ma->GetFacePNums (fnr, pnums);
          INT<3> i3(pnums[0], pnums[1], pnums[2]);
          i3.Sort();
          v2f[i3] = fnr;
        }

      for (auto fnr : Range(ma->GetNFaces()))
        {
          ma->GetFacePNums(fnr,pnums);
          INT<3> mv(vertex_map[pnums[0]],vertex_map[pnums[1]],vertex_map[pnums[2]]);
          if(mv[0] != pnums[0] && mv[1] != pnums[1] && mv[2] != pnums[2])
            {
              mv.Sort();
              int mfnr = v2f[mv];
              node_ids.Append(std::make_tuple(NodeId(NT_FACE,mfnr),NodeId(NT_FACE,fnr)));
            }
        }
    }
    
  
}
