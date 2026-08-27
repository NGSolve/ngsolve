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
    
  Array<int> MortonElementOrder (const shared_ptr<MeshAccess> & ma)
  {
    size_t ne = ma->GetNE(VOL);

    // element order along Morton (Z-)curve of element centroids
    Array<uint64_t> code(ne);
    Vec<3> pmin(0.0), pmax(0.0);
    Array<Vec<3>> cent(ne);
    for (size_t elnr = 0; elnr < ne; elnr++)
      {
        auto vs = ma->GetElement(ElementId(VOL, elnr)).Vertices();
        Vec<3> c(0.0);
        for (auto v : vs)
          c += ma->GetPoint<3>(v);
        c *= 1.0 / vs.Size();
        cent[elnr] = c;
        if (elnr == 0) { pmin = c; pmax = c; }
        for (int j = 0; j < 3; j++)
          {
            pmin(j) = min2(pmin(j), c(j));
            pmax(j) = max2(pmax(j), c(j));
          }
      }

    // spread lower 21 bits with 2-bit gaps
    auto spreadbits = [](uint64_t x)
      {
        x &= 0x1fffff;
        x = (x | x << 32) & 0x001f00000000ffffULL;
        x = (x | x << 16) & 0x001f0000ff0000ffULL;
        x = (x | x << 8)  & 0x100f00f00f00f00fULL;
        x = (x | x << 4)  & 0x10c30c30c30c30c3ULL;
        x = (x | x << 2)  & 0x1249249249249249ULL;
        return x;
      };

    for (size_t elnr = 0; elnr < ne; elnr++)
      {
        uint64_t q[3];
        for (int j = 0; j < 3; j++)
          {
            double h = pmax(j) - pmin(j);
            double rel = (h > 0) ? (cent[elnr](j) - pmin(j)) / h : 0.0;
            q[j] = uint64_t(min2(rel, 1.0) * 2097151.0);   // 21 bits
          }
        code[elnr] = spreadbits(q[0]) | (spreadbits(q[1]) << 1) | (spreadbits(q[2]) << 2);
      }

    Array<int> elorder(ne);
    for (auto i : Range(ne))
      elorder[i] = i;
    QuickSortI (code, elorder);
    return elorder;
  }


  Array<int> RCMElementOrder (const shared_ptr<MeshAccess> & ma)
  {
    static Timer t("RCMElementOrder"); RegionTimer reg(t);
    size_t ne = ma->GetNE(VOL);
    size_t nv = ma->GetNV();

    // vertex -> element table
    TableCreator<int> vecreator(nv);
    for ( ; !vecreator.Done(); vecreator++)
      for (size_t i = 0; i < ne; i++)
        for (auto v : ma->GetElement(ElementId(VOL,i)).Vertices())
          vecreator.Add (v, i);
    Table<int> v2el = vecreator.MoveTable();

    // neighbours = elements sharing a vertex (deduplicated via stamps)
    Array<int> stamp(ne);
    stamp = -1;
    int stampcnt = 0;
    Array<int> nbs;
    auto GetNeighbours = [&] (int elnr, Array<int> & nbs)
    {
      nbs.SetSize0();
      stampcnt++;
      for (auto v : ma->GetElement(ElementId(VOL,elnr)).Vertices())
        for (auto nb : v2el[v])
          if (nb != elnr && stamp[nb] != stampcnt)
            {
              stamp[nb] = stampcnt;
              nbs.Append(nb);
            }
    };

    Array<int> degree(ne);
    for (size_t i = 0; i < ne; i++)
      {
        GetNeighbours(i, nbs);
        degree[i] = nbs.Size();
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
          GetNeighbours(workq[qi], nbs);
          QuickSort (nbs, [&] (int a, int b) { return degree[a] < degree[b]; });
          for (auto nb : nbs)
            if (!done[nb] && mark[nb] != epoch)
              {
                mark[nb] = epoch;
                workq.Append(nb);
              }
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


  Array<int> LocalityElementOrder (const shared_ptr<MeshAccess> & ma)
  {
    // topological (RCM) default; MortonElementOrder is the geometric alternative
    return RCMElementOrder (ma);
  }


  void ReorderedFESpace :: Update()
  {
    space->Update();
    FESpace::Update();

    SetNDof(space->GetNDof());
    size_t ndof = space->GetNDof();
    size_t ne = ma->GetNE(VOL);
    Array<DofId> dofs;

    elorder = LocalityElementOrder (ma);

    // first-touch dof numbering along the Morton element order;
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
