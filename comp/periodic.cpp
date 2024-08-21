
/*********************************************************************/
/* File:   periodic.cpp                                              */
/* Author: Christopher Lackner                                       */
/* Date:   Feb. 2017                                                 */
/*********************************************************************/


#include "periodic.hpp"

namespace ngcomp {

  PeriodicFESpace :: PeriodicFESpace (shared_ptr<FESpace> aspace, const Flags & flags, shared_ptr<Array<int>> aused_idnrs)
    : FESpace(aspace->GetMeshAccess(), flags), space(aspace), used_idnrs(aused_idnrs)
    {
      type = "Periodic" + space->type;
      for (auto vb : { VOL,BND, BBND, BBBND })
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
          low_order_space = make_shared<PeriodicFESpace>(space->LowOrderFESpacePtr(),lo_flags,used_idnrs);
        }
    }
    
  void PeriodicFESpace :: Update()
    {      
      space->Update();
      FESpace::Update();
      dofmap.SetSize (space->GetNDof());
      vertex_map.SetSize(ma->GetNP());
      for (int i = 0; i < dofmap.Size(); i++)
	dofmap[i] = i;
      for (int i : Range(vertex_map.Size()))
        vertex_map[i] = i;

      for (auto idnr : Range(ma->GetNPeriodicIdentifications()))
        {
	  if (used_idnrs->Size() && !used_idnrs->Contains(idnr)) continue;
          Array<int> minion_dofnrs;
          Array<int> master_dofnrs;
          for (auto node_type : {NT_VERTEX, NT_EDGE, NT_FACE})
            {
              const auto & periodic_nodes = ma->GetPeriodicNodes(node_type, idnr);
              if(node_type==NT_VERTEX)
                {
                  for (const auto& per_verts : periodic_nodes)
                    vertex_map[per_verts[1]] = vertex_map[per_verts[0]];
                }
              for(const auto& node_pair : periodic_nodes)
                {
                  space->GetDofNrs(NodeId(node_type,node_pair[0]),master_dofnrs);
                  space->GetDofNrs(NodeId(node_type,node_pair[1]),minion_dofnrs);
                  for(auto i : Range(master_dofnrs.Size()))
                    {
                    dofmap[minion_dofnrs[i]] = dofmap[master_dofnrs[i]];
		    DofMapped(minion_dofnrs[i],master_dofnrs[i],idnr);
                    }
                }
            }
        }
      bool changed = true;
      while(changed)
        {
          changed = false;
          for (auto i : Range(dofmap))
            {
              if(dofmap[dofmap[i]] != dofmap[i])
                {
                  dofmap[i] = dofmap[dofmap[i]];
                  changed = true;
                }
            }
        }
      changed = true;
      while(changed)
        {
          changed = false;
          for (auto i : Range(vertex_map))
            {
              if(vertex_map[vertex_map[i]] != vertex_map[i])
                {
                  vertex_map[i] = vertex_map[vertex_map[i]];
                  changed = true;
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
    
  void PeriodicFESpace :: DoArchive(Archive & archive)
    {
      FESpace::DoArchive(archive);
      archive.Shallow(space);
      archive & dofmap & vertex_map & used_idnrs;
    }

  FiniteElement& PeriodicFESpace :: GetFE (ElementId ei, Allocator & alloc) const
    {
      auto & fe = space->GetFE(ei,alloc);
      const auto & ngel = ma->GetElement(ei);

      /*
      SwitchET (ngel.GetType(), [&](auto et)
                {
                  if (auto hofe = dynamic_cast<VertexOrientedFE<et.ElementType()>*>(&fe))
                    hofe->SetVertexNumbers(vertex_map[ngel.Vertices()]);
                });
      */

      fe.SetVertexNumbers( ArrayMem<int,8> (vertex_map[ngel.Vertices()]) );

      /*
      switch (ngel.GetType())
	{
	case ET_TRIG:
	  {
            auto hofe = dynamic_cast<VertexOrientedFE<ET_TRIG>*>(&fe);
            if(hofe)
              hofe->SetVertexNumbers(vertex_map[ngel.Vertices()]);
            break;
	  }
	case ET_QUAD:
	  {
            auto hofe = dynamic_cast<VertexOrientedFE<ET_QUAD>*>(&fe);
            if(hofe)
              hofe->SetVertexNumbers(vertex_map[ngel.Vertices()]);
            break;
	  }
	case ET_TET:
	  {
            auto hofe = dynamic_cast<VertexOrientedFE<ET_TET>*>(&fe);
            if(hofe)
              hofe->SetVertexNumbers(vertex_map[ngel.Vertices()]);
            break;
	  }
	case ET_PRISM:
	  {
            auto hofe = dynamic_cast<VertexOrientedFE<ET_PRISM>*>(&fe);
            if(hofe)
              hofe->SetVertexNumbers(vertex_map[ngel.Vertices()]);
            break;
	  }
	case ET_HEX:
	  {
            auto hofe = dynamic_cast<VertexOrientedFE<ET_HEX>*>(&fe);
            if(hofe)
              hofe->SetVertexNumbers(vertex_map[ngel.Vertices()]);
            break;
	  }
        case ET_SEGM:
          {
            auto hofe = dynamic_cast<VertexOrientedFE<ET_SEGM>*>(&fe);
            if (hofe)
              hofe->SetVertexNumbers(vertex_map[ngel.Vertices()]);
            break;
          }
        default:
          throw Exception("ElementType not implemented for PeriodicFESpace::GetFE");
	}
      */
      return fe;
    }

  
  void PeriodicFESpace :: GetDofNrs(ElementId ei, Array<DofId> & dnums) const
    {
      space->GetDofNrs(ei,dnums);
      for (auto & d : dnums)
        if (IsRegularDof(d)) d = dofmap[d];
    }

    void PeriodicFESpace :: GetDofNrs (NodeId ni, Array<DofId> & dnums) const
    {
      space->GetDofNrs(ni, dnums);
      for (auto & d : dnums)
        if (IsRegularDof(d)) d = dofmap[d];
    }

    void PeriodicFESpace :: GetVertexDofNrs (int vnr,  Array<DofId> & dnums) const
    { 
      space->GetVertexDofNrs(vnr, dnums); 
      for (auto & d : dnums)
        if (IsRegularDof(d)) d = dofmap[d];
    }

    void PeriodicFESpace :: GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const
    { 
      space->GetEdgeDofNrs (ednr, dnums); 
      for (auto & d : dnums)
        if (IsRegularDof(d)) d = dofmap[d];
    }
    
    void PeriodicFESpace :: GetFaceDofNrs (int fanr, Array<DofId> & dnums) const
    { 
      space->GetFaceDofNrs(fanr, dnums); 
      for (auto & d : dnums)
        if (IsRegularDof(d)) d = dofmap[d];
    }

  ProxyNode PeriodicFESpace ::
  MakeProxyFunction (bool testfunction,
                     const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const 
  {
    /*
    auto proxy = GetBaseSpace()->MakeProxyFunction (testfunction, addblock);
    shared_ptr<FESpace> fes = dynamic_pointer_cast<FESpace> (const_cast<PeriodicFESpace*>(this)->shared_from_this());
    proxy.SetFESpace(fes);
    return proxy;
    */
    return GetBaseSpace()->MakeProxyFunction (testfunction,
                                              [&] (shared_ptr<ProxyFunction> proxy)
                                                {
                                                  shared_ptr<FESpace> fes = dynamic_pointer_cast<FESpace> (const_cast<PeriodicFESpace*>(this)->shared_from_this());
                                                  proxy->SetFESpace(fes);
                                                  return addblock (proxy);
                                                });
  }
  
  
  template<typename TSCAL>
  QuasiPeriodicFESpace<TSCAL> :: QuasiPeriodicFESpace(shared_ptr<FESpace> fespace, const Flags & flags, shared_ptr<Array<int>> aused_idnrs, shared_ptr<Array<TSCAL>> afactors) :
    PeriodicFESpace(fespace, flags, aused_idnrs), factors(afactors)
  {  }

  template<typename TSCAL>
  void QuasiPeriodicFESpace<TSCAL> :: Update()
  {
    space->Update();
    dof_factors.SetSize(space->GetNDof());

    dof_factors = 1.;
    master_dofs.SetSize(space->GetNDof());
    for(auto& md : master_dofs)
      md = {};
    PeriodicFESpace::Update();
  }

  template<typename TSCAL>
  void QuasiPeriodicFESpace<TSCAL> :: VTransformMR (ElementId ei, SliceMatrix<double> mat, TRANSFORM_TYPE tt) const
  {
    if constexpr(!is_same_v<TSCAL, double>)
        throw Exception("Shouldn't get here: QuasiPeriodicFESpace::TransformMR, space with complex factors should be complex");
    else
      {
        PeriodicFESpace::VTransformMR(ei, mat, tt);
        Array<int> dofnrs;
        space->GetDofNrs(ei, dofnrs);
        for (int i : Range(dofnrs.Size()))
          {
            if (dofnrs[i] != dofmap[dofnrs[i]])
              {
                if (tt & TRANSFORM_MAT_LEFT)
                  mat.Row(i) *= dof_factors[dofnrs[i]];
                if (tt & TRANSFORM_MAT_RIGHT)
                  mat.Col(i) *= dof_factors[dofnrs[i]];
              }
          }
      }
  }

  template<typename TSCAL>
  void QuasiPeriodicFESpace<TSCAL> :: VTransformMC (ElementId ei, SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const
  {
    PeriodicFESpace::VTransformMC(ei, mat, tt);
    Array<int> dofnrs;
    space->GetDofNrs(ei, dofnrs);
    for (int i : Range(dofnrs.Size()))
      {
	if (dofnrs[i] != dofmap[dofnrs[i]])
	  {
	    if (tt & TRANSFORM_MAT_LEFT)
	      mat.Row(i) *= conj(dof_factors[dofnrs[i]]);
	    if (tt & TRANSFORM_MAT_RIGHT)
	      mat.Col(i) *= dof_factors[dofnrs[i]];
	  }
      }
  }

  template<typename TSCAL>
  void QuasiPeriodicFESpace<TSCAL> :: VTransformVR (ElementId ei, SliceVector<double> vec, TRANSFORM_TYPE tt) const
  {
    if constexpr(!is_same_v<TSCAL, double>)
        throw Exception("Shouldn't get here: QuasiPeriodicFESpace::TransformMR, space with complex factors should be complex");
    else
      {
        PeriodicFESpace::VTransformVR(ei, vec, tt);
        Array<int> dofnrs;
        space->GetDofNrs(ei, dofnrs);
        for (auto i : Range(dofnrs))
          {
            if (dofnrs[i] != dofmap[dofnrs[i]])
              {
                if (tt == TRANSFORM_RHS)
                  vec[i] *= dof_factors[dofnrs[i]];
                else if (tt == TRANSFORM_SOL)
                  vec[i] *= dof_factors[dofnrs[i]];
                else // TRANSFORM_SOL_INVERSE
                  vec[i] /= dof_factors[dofnrs[i]];
              }
          }
      }
  }

  template<typename TSCAL>
  void QuasiPeriodicFESpace<TSCAL> :: VTransformVC (ElementId ei, SliceVector<Complex> vec, TRANSFORM_TYPE tt) const 
  {
    PeriodicFESpace::VTransformVC(ei, vec, tt);
    Array<int> dofnrs;
    space->GetDofNrs(ei, dofnrs);
    for (auto i : Range(dofnrs))
      {
	if (dofnrs[i] != dofmap[dofnrs[i]])
	  {
	    if (tt == TRANSFORM_RHS)
	      vec[i] *= conj(dof_factors[dofnrs[i]]);
	    else if (tt == TRANSFORM_SOL)
	      vec[i] *= dof_factors[dofnrs[i]];
            else // TRANSFORM_SOL_INVERSE
              vec[i] /= dof_factors[dofnrs[i]];
	  }
      }
  }

  template<typename TSCAL>
  void QuasiPeriodicFESpace<TSCAL> :: DofMapped(size_t from, size_t to, size_t idnr)
  {
    // if the same dofs are mapped twice by different identification numbers only multiply once with
    // the factor!
    auto& md = master_dofs[from];
    if(md.find(to) == md.end())
      {
        dof_factors[from] *= (*factors)[idnr];
        md.insert(to);
      }
  }

  template class QuasiPeriodicFESpace<double>;
  template class QuasiPeriodicFESpace<Complex>;
  
  static RegisterClassForArchive<PeriodicFESpace, FESpace> reg_periodic;
  static RegisterClassForArchive<QuasiPeriodicFESpace<double>, PeriodicFESpace> reg_qperiodic_d;
  static RegisterClassForArchive<QuasiPeriodicFESpace<Complex>, PeriodicFESpace> reg_qperiodic_c;
}
