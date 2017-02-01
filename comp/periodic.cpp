#include <comp.hpp>
#include <python_ngstd.hpp>

namespace ngcomp {

  template<typename BASE>
  class PeriodicFESpace : public BASE
  {
    Array<int> dofmap; // mapping of dofs
    using BASE::GetNDof;
    using BASE::ma;
    using BASE::ctofdof;
    

  public:
    PeriodicFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
      : BASE(ama,flags)
    {
      ;
    }
    
    virtual ~PeriodicFESpace () { ; }
    virtual void Update (LocalHeap & lh) override
    { 
      BASE::Update (lh);
            
      dofmap.SetSize (GetNDof());
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
              FESpace::GetDofNrs(id1,master_dofnrs);
              FESpace::GetDofNrs(id2,slave_dofnrs);
              for(auto i : Range(master_dofnrs.Size())){
                dofmap[slave_dofnrs[i]] = dofmap[master_dofnrs[i]];
              }
            }
        }
      for (int i = 0; i < dofmap.Size(); i++)
	if (dofmap[i] != i){
	  ctofdof[i] = UNUSED_DOF;
	}
    }
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override
    {
      auto & fe = BASE::GetFE(ei,alloc);
      const auto & ngel = ma->GetElement(ei);
      switch (ngel.GetType())
	{
	case ET_TRIG:
	  {
            try
              {
                auto & hofe = dynamic_cast<VertexOrientedFE<ET_TRIG>&> (fe);
                hofe.SetVertexNumbers (dofmap[ngel.Vertices()]);
              }
            catch(exception e) { ; }
            break;
            
	  }
	case ET_TET:
	  {
            try
              {
                auto & hofe = dynamic_cast<VertexOrientedFE<ET_TET>&> (fe);
                hofe.SetVertexNumbers (dofmap[ngel.Vertices()]);
              }
            catch(exception e) { ; }
            break;
	  }
        default:
          throw Exception("ElementType not implemented for PeriodicFESpace::GetFE");
	}
      return fe;
    }


    virtual void GetDofNrs(ElementId ei, Array<int> & dnums) const override
    {
      BASE::GetDofNrs(ei,dnums);
      for (int i = 0; i< dnums.Size(); i++)
	dnums[i] = dofmap[dnums[i]];
    }

  private:
    void GetPeriodicNodeIds(Array<std::tuple<NodeId,NodeId>> & node_ids,int idnr)
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
    
  };

  static RegisterFESpace<PeriodicFESpace<H1HighOrderFESpace>> initperh1 ("perH1");
  static RegisterFESpace<PeriodicFESpace<HCurlHighOrderFESpace>> initperhcurl ("perHCurl");
  static RegisterFESpace<PeriodicFESpace<HDivHighOrderFESpace>> initperhdiv ("perHDiv");
  
}
