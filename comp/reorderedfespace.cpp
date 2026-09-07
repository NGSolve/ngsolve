/*********************************************************************/
/* File:   reordered.cpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   Jun 2019                                                  */
/*********************************************************************/


#include "reorderedfespace.hpp"

namespace ngcomp {
  
  ReorderedFESpace :: ReorderedFESpace (shared_ptr<FESpace> aspace, const Flags & flags)
    : FESpace(aspace->GetMeshAccess(), flags), space(aspace)
  {
    type = "Reordered" + space->type;
    evaluator[VOL] = space->GetEvaluator(VOL);
    evaluator[BND] = space->GetEvaluator(BND);
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
    
  // reverse Cuthill-McKee on the graph with nodes 0..ne-1, two nodes adjacent
  // iff they share a key in el2key; negative keys are skipped
  static Array<int> RCMOrderFromTable (FlatTable<int> el2key, size_t nkeys)
  {
    size_t ne = el2key.Size();

    // key -> element table
    TableCreator<int> kecreator(nkeys);
    for ( ; !kecreator.Done(); kecreator++)
      for (size_t i = 0; i < ne; i++)
        for (auto k : el2key[i])
          if (k >= 0)
            kecreator.Add (k, i);
    Table<int> k2el = kecreator.MoveTable();

    // call func for every element sharing a key with elnr (each one once,
    // deduplicated via stamps; elnr itself is stamped out)
    Array<int> stamp(ne);
    stamp = -1;
    int stampcnt = 0;
    auto IterateNeighbours = [&] (int elnr, auto func)
    {
      stampcnt++;
      stamp[elnr] = stampcnt;
      for (auto k : el2key[elnr])
        if (k >= 0)
          for (auto nb : k2el[k])
            if (stamp[nb] != stampcnt)
              {
                stamp[nb] = stampcnt;
                func(nb);
              }
    };

    Array<int> degree(ne);
    for (size_t i = 0; i < ne; i++)
      {
        int cnt = 0;
        IterateNeighbours (i, [&cnt] (int nb) { cnt++; });
        degree[i] = cnt;
      }

    // Cuthill-McKee: BFS, children by ascending degree; reversed at the end
    Array<bool> done(ne);
    done = false;
    Array<int> mark(ne);
    mark = -1;
    int epoch = 0;
    Array<int> workq;
    workq.SetAllocSize(ne);

    auto BFS = [&] (int start)
    {
      epoch++;
      workq.SetSize0();
      workq.Append(start);
      mark[start] = epoch;
      for (size_t qi = 0; qi < workq.Size(); qi++)
        {
          size_t first = workq.Size();
          IterateNeighbours (workq[qi], [&] (int nb)
                             {
                               if (!done[nb] && mark[nb] != epoch)
                                 {
                                   mark[nb] = epoch;
                                   workq.Append(nb);
                                 }
                             });
          // children in ascending degree
          QuickSort (workq.Range(first, workq.Size()),
                     [&degree] (int a, int b) { return degree[a] < degree[b]; });
        }
    };

    Array<int> order;
    order.SetAllocSize(ne);
    for (size_t seed = 0; seed < ne; seed++)
      if (!done[seed])
        {
          // pseudo-peripheral start: end point of repeated BFS
          int start = seed;
          for (int it = 0; it < 2; it++)
            {
              BFS(start);
              start = workq.Last();
            }
          BFS(start);
          for (auto el : workq)
            {
              done[el] = true;
              order.Append(el);
            }
        }

    for (size_t i = 0; i < ne/2; i++)
      Swap (order[i], order[ne-1-i]);
    return order;
  }


  // element ordering for locality: RCM on the element graph of the space.
  // Two elements are adjacent if they share a mesh vertex or a dof: the dof
  // part sees identifications the mesh does not (periodic spaces), the vertex
  // part keeps the graph connected for discontinuous spaces
  static Array<int> RCMElementOrder (const FESpace & fes)
  {
    static Timer t("RCMElementOrder"); RegionTimer reg(t);

    auto ma = fes.GetMeshAccess();
    size_t ne = ma->GetNE(VOL);
    size_t nv = ma->GetNV();
    Array<DofId> dofs;

    TableCreator<int> creator(ne);
    for ( ; !creator.Done(); creator++)
      for (size_t i = 0; i < ne; i++)
        {
          ElementId ei(VOL, i);
          for (auto v : ma->GetElement(ei).Vertices())
            creator.Add (i, v);
          fes.GetDofNrs (ei, dofs);
          for (auto d : dofs)
            if (IsRegularDof(d))
              creator.Add (i, nv+d);
        }
    Table<int> el2key = creator.MoveTable();

    return RCMOrderFromTable (el2key, nv + fes.GetNDof());
  }


  void RCMReorderSubset (FlatArray<size_t> els, FlatTable<int> el2dof, size_t ndof)
  {
    if (els.Size() <= 2) return;

    // induced sub-table: row i = dofs of element els[i]
    TableCreator<int> creator(els.Size());
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(els))
        for (auto d : el2dof[els[i]])
          creator.Add (i, d);
    Table<int> sub = creator.MoveTable();

    // disconnected components are seeded in the current order of els,
    // so a locality-aware input order is preserved across components
    Array<int> perm = RCMOrderFromTable (sub, ndof);

    Array<size_t> tmp(els.Size());
    for (auto i : Range(els))
      tmp[i] = els[perm[i]];
    for (auto i : Range(els))
      els[i] = tmp[i];
  }


  void ReorderedFESpace :: Update()
  {
    space->Update();
    FESpace::Update();

    SetNDof(space->GetNDof());
    size_t ndof = space->GetNDof();
    size_t ne = ma->GetNE(VOL);
    Array<DofId> dofs;

    elorder = RCMElementOrder (*space);

    // first-touch dof numbering along the element order;
    // record cluster boundaries every 'step' elements (contiguous dof ranges)
    constexpr size_t step = 20;
    dofmap.SetSize(ndof);
    dofmap = -1;
    Array<size_t> bounds;
    size_t cnt = 0;
    for (auto i : Range(ne))
      {
        if (i % step == 0) bounds.Append(cnt);
        space->GetDofNrs (ElementId(VOL, elorder[i]), dofs);
        for (auto d : dofs)
          if (IsRegularDof(d) && dofmap[d] == -1)
            dofmap[d] = cnt++;
      }
    // dofs not touched by volume elements (surface-only, unused, ...)
    for (auto & d : dofmap)
      if (d == -1)
        d = cnt++;
    bounds.Append(cnt);

    ctofdof.SetSize(ndof);
    for (auto i : Range(ndof))
      ctofdof[dofmap[i]] = space->GetDofCouplingType(i);

    {
      // clusters are contiguous ranges of the first-touch numbering
      size_t ngroups = bounds.Size()-1;
      Array<int> sizes(ngroups);
      for (auto i : Range(ngroups))
        sizes[i] = bounds[i+1] - bounds[i];
      clusters = make_shared<Table<int>> (sizes);
      for (auto i : Range(ngroups))
        for (auto j : Range(sizes[i]))
          (*clusters)[i][j] = bounds[i] + j;
    }
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
	  { free_dofs->SetBit(dofmap[k]); }
	if (external_space_free_dofs->Test(k))
	  { external_free_dofs->SetBit(dofmap[k]); }
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
      if (IsRegularDof(d))
        d = dofmap[d];
  }

  void ReorderedFESpace :: GetDofNrs (NodeId ni, Array<DofId> & dnums) const
  {
    space->GetDofNrs (ni, dnums);
    for (auto & d : dnums)
      if (IsRegularDof(d))
        d = dofmap[d];
  }
  
  void ReorderedFESpace :: GetVertexDofNrs (int vnr,  Array<DofId> & dnums) const
  {
    space->GetVertexDofNrs (vnr, dnums);
    for (auto & d : dnums)
      if (IsRegularDof(d))
        d = dofmap[d];
  }
  
  void ReorderedFESpace :: GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const
  {
    space->GetEdgeDofNrs (ednr, dnums);
    for (auto & d : dnums)
      if (IsRegularDof(d))
        d = dofmap[d];
  }
    
  void ReorderedFESpace :: GetFaceDofNrs (int fanr, Array<DofId> & dnums) const
  {
    space->GetFaceDofNrs (fanr, dnums);
    for (auto & d : dnums)
      if (IsRegularDof(d))
        d = dofmap[d];
  }
}
