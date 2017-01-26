#include <comp.hpp>
#include "tpfes.hpp"
using namespace ngmg;

namespace ngcomp
{
  using ngfem::TPHighOrderFE;

  TPHighOrderFESpace::TPHighOrderFESpace (FlatArray<shared_ptr<FESpace>> spaces, const Flags & flags, bool parseflags, Array<int> * el_counts)
    : FESpace (spaces[0]->GetMeshAccess(), flags)
  {
    nmeshes = spaces.Size();
    fespaces.SetSize(nmeshes);
    fespaces=spaces;
    space_x = spaces[0];
    spaces_y.SetSize(1);
    spaces_y[0] = spaces[1];
    meshes.SetSize(nmeshes);
    ndofs.SetSize(nmeshes);
    nels.SetSize(nmeshes);
    nfacets.SetSize(nmeshes);
    ndof = 1;
    nel = 1;
    for (int i : Range(nmeshes)) 
    {
      ndofs[i] = spaces[i]->GetNDof();
      meshes[i] = spaces[i]->GetMeshAccess();      
      nels[i] = meshes[i]->GetNE();
      nfacets[i] = meshes[i]->GetNFacets();
      ndof *= ndofs[i];
      nel *= nels[i];
    }
    nelsyinverse = 1.0/nels[1];
    first_element_dof.SetSize(nel+1);
    
    LocalHeap lh(100000,"Setup TP Space"); 

    int ii=0;
    first_element_dof[0] = 0;
    for(int elx : Range(nels[0]) )
    {
      int ndofx = space_x->GetFE(ElementId(VOL,elx),lh).GetNDof();
      for(int ely : Range(nels[1]) )
      {
        HeapReset hr(lh);
        int ndofy = Space(elx)->GetFE(ElementId(VOL,ely),lh).GetNDof();
        first_element_dof[ii+1] = first_element_dof[ii] + ndofx*ndofy;
        ii++;
      }
    }
    Array<shared_ptr<DifferentialOperator>> evaluators(nmeshes);
    for (int i : Range(nmeshes))
      evaluators[i] = spaces[i]->GetEvaluator();
    int dim = 0;
    for (auto eval : evaluators)
      dim = max2(dim, eval->Dim());
    int difforder = evaluators[0]->DiffOrder();
    for (auto eval : evaluators)
      difforder = min2(difforder, eval->DiffOrder());
    int blockdim = 1;
    evaluator[VOL] = shared_ptr<DifferentialOperator>( new TPDifferentialOperator(evaluators, dim, blockdim, VOL, difforder) );
  }

  TPHighOrderFESpace::TPHighOrderFESpace (shared_ptr<FESpace> aspace_x,FlatArray<shared_ptr<FESpace>> aspaces_y, const Flags & flags, bool parseflags)
    : FESpace (aspace_x->GetMeshAccess(), flags)
  {
    nmeshes = 2;
    fespaces.SetSize(nmeshes);
    space_x = aspace_x;
    fespaces[0]=space_x;
    spaces_y.SetSize(aspaces_y.Size());
    spaces_y = aspaces_y;
    fespaces[1] = spaces_y[0];
    meshes.SetSize(nmeshes);
    ndofs.SetSize(nmeshes);
    nels.SetSize(nmeshes);
    nfacets.SetSize(nmeshes);
    Array<int> temp(2);
    ndof = 1;
    nel = 1;
    meshes[0] = fespaces[0]->GetMeshAccess();
    meshes[1] = fespaces[1]->GetMeshAccess();
    for (int i : Range(nmeshes)) 
    {
      nels[i] = meshes[i]->GetNE();
      nfacets[i] = meshes[i]->GetNFacets();
      nel *= nels[i];
    }
    nelsyinverse = 1.0/nels[1];
    ndof = 0;
    LocalHeap lh(100000,"Setup TP Space");
    for(int i=0;i<space_x->GetMeshAccess()->GetNE();i++)
      ndof += fespaces[0]->GetFE(ElementId(VOL,i),lh).GetNDof()*spaces_y[i]->GetNDof();
    
    first_element_dof.SetSize(nel+1);
    int ii=0;
    first_element_dof[0] = 0;
    for(int elx : Range(nels[0]) )
    {
      int ndofx = space_x->GetFE(ElementId(VOL,elx),lh).GetNDof();
      for(int ely : Range(nels[1]) )
      {
        HeapReset hr(lh);
        int ndofy = Space(elx)->GetFE(ElementId(VOL,ely),lh).GetNDof();
        first_element_dof[ii+1] = first_element_dof[ii] + ndofx*ndofy;
        ii++;
      }
    }
    
    Array<shared_ptr<DifferentialOperator>> evaluators(nmeshes);
    
    evaluators[0] = space_x->GetEvaluator();
    evaluators[1] = spaces_y[0]->GetEvaluator();
    int dim = 0;
    for (auto eval : evaluators)
      dim = max2(dim, eval->Dim());
    int difforder = evaluators[0]->DiffOrder();
    for (auto eval : evaluators)
      difforder = min2(difforder, eval->DiffOrder());
    int blockdim = 1;
    evaluator[VOL] = shared_ptr<DifferentialOperator>( new TPDifferentialOperator(evaluators, dim, blockdim, VOL, difforder) );
  }

  TPHighOrderFESpace::~TPHighOrderFESpace () { ; }

  SymbolTable<shared_ptr<DifferentialOperator>> TPHighOrderFESpace::GetAdditionalEvaluators () const
  {
    SymbolTable<shared_ptr<DifferentialOperator>> ops;
    ArrayMem<shared_ptr<DifferentialOperator>,2> gradx(2);
    ArrayMem<shared_ptr<DifferentialOperator>,2> grady(2);
    gradx[0] = space_x->GetFluxEvaluator();
    gradx[1] = spaces_y[0]->GetEvaluator();
    grady[0] = space_x->GetEvaluator();
    grady[1] = spaces_y[0]->GetFluxEvaluator();
    int dim = 0;
    for (auto eval : gradx)
      dim = max2(dim, eval->Dim());
    int difforder = gradx[0]->DiffOrder();
    for (auto eval : gradx)
      difforder = min2(difforder, eval->DiffOrder());
    int blockdim = 1;
    ops.Set("gradx", make_shared<TPDifferentialOperator>( gradx, dim, blockdim, VOL, difforder ));
    dim = 0;
    for (auto eval : grady)
      dim = max2(dim, eval->Dim());
    difforder = grady[0]->DiffOrder();
    for (auto eval : grady)
      difforder = min2(difforder, eval->DiffOrder());
    ops.Set("grady", make_shared<TPDifferentialOperator>( grady,dim,blockdim,VOL,difforder ));
    return ops;
  } 

  void TPHighOrderFESpace :: ReduceToXSpace(shared_ptr<GridFunction> gf_in, shared_ptr<GridFunction> gf_out,LocalHeap & clh,
                                              const function<void(shared_ptr<FESpace>,const FiniteElement &,const ElementTransformation &,FlatVector<>,FlatVector<>,LocalHeap&)> & func)
  {
    BaseVector & vec_in = gf_in->GetVector();
    BaseVector & vec_out = gf_out->GetVector();
    Array<int> dnums,dnumsx;
    for(int i=0;i<nels[0];i++)
    {
      FlatVector<> elvec_out(space_x->GetFE(ElementId(VOL,i),clh).GetNDof(),clh);
      elvec_out = 0.0;
      for(int j=0;j<nels[1];j++)
      {
        HeapReset hr(clh);
        int index = GetIndex(i,j);
        GetDofNrs(index,dnums);
        FlatVector<> elvec(dnums.Size(),clh);
        vec_in.GetIndirect(dnums,elvec);
        const TPHighOrderFE & tpfel = dynamic_cast<const TPHighOrderFE &>(GetFE(ElementId(index),clh));
        ElementTransformation & trafo = fespaces[1]->GetMeshAccess()->GetTrafo(j,clh);
        func(gf_in->GetFESpace(),tpfel,trafo,elvec,elvec_out,clh);
      }
      space_x->GetDofNrs(ElementId(VOL,i),dnumsx);
      gf_out->GetVector().SetIndirect(dnumsx,elvec_out);
    }
  }

  void TPHighOrderFESpace :: ProlongateFromXSpace(shared_ptr<GridFunction> gf_in, shared_ptr<GridFunction> gf_out,LocalHeap & lh)
  {
    BaseVector & vec_in = gf_in->GetVector();
    BaseVector & vec_out = gf_out->GetVector();
    Array<int> dnums,dnumsx;
    for(int i=0;i<nels[0];i++)
    {
      Vector<> elvec_in(space_x->GetFE(ElementId(VOL,i),lh).GetNDof());
      fespaces[0]->GetDofNrs(ElementId(VOL,i),dnumsx);
      vec_in.GetIndirect(dnumsx,elvec_in);
      for(int j=0;j<nels[1];j++)
      {
        HeapReset hr(lh);
        int index = GetIndex(i,j);
        GetDofNrs(index,dnums);
        FlatVector<> elvec_out(dnums.Size(),lh);
        FlatMatrix<> elmat_out(dnumsx.Size(),dnums.Size()/dnumsx.Size(),&elvec_out(0));
        elmat_out = 0.0;
        elmat_out.Col(0) = elvec_in;
        vec_out.SetIndirect(dnums,elvec_out);
      }
    }
  }
  
  void TPHighOrderFESpace::FinalizeUpdate(LocalHeap & lh) 
  {
    space_x->FinalizeUpdate(lh);
    for(auto fes : spaces_y)
      fes->FinalizeUpdate(lh);
    FESpace::FinalizeUpdate(lh);
    element_coloring[VOL] = Table<int>(nel,1);
    for (int i : Range(nel))
      element_coloring[VOL][i][0] = i;
  }

  void TPHighOrderFESpace::Update(LocalHeap & lh)
  {
    space_x->Update(lh);
    for (auto fes : spaces_y)
      fes->Update(lh);
    FESpace::Update(lh);
  }
  
  void TPHighOrderFESpace::UpdateDofTables()
  {
    throw Exception("TPHighOrderFESpace::UpdateDofTables() not implemented");
  }
  
  void TPHighOrderFESpace::UpdateCouplingDofArray() 
  {
    throw Exception("TPHighOrderFESpace::UpdateCouplingDofArray() not implemented");
  }
  
  size_t TPHighOrderFESpace::GetNDof () const throw() 
  {
    return ndof;
  }
  
  size_t TPHighOrderFESpace::GetNDofLevel (int level) const 
  {
    throw Exception("TPHighOrderFESpace::GetNDofLevel() not implemented");
  }
  
  FiniteElement & TPHighOrderFESpace::GetFE (ElementId ei, Allocator & lh) const
  {
    int elnr = ei.Nr();
    ArrayMem<int,2> elnums(2);
    GetIndices(elnr, elnums);
    ArrayMem<const FiniteElement *,2> els(2);
    const FiniteElement & ref_elx = space_x->GetFE( ElementId(VOL,elnums[0]), lh );
    els[0] = &ref_elx;
    const FiniteElement & ref_ely = Space(elnums[0])->GetFE( ElementId(VOL,elnums[1]), lh );
    els[1] = &ref_ely;
    TPHighOrderFE * fe = new (lh) TPHighOrderFE (els);
    return *fe;
  }
  
  const FiniteElement & TPHighOrderFESpace::GetFE (int elnr, LocalHeap & lh) const
  {
    ArrayMem<int,2> elnums(2);
    GetIndices(elnr, elnums);
    ArrayMem<const FiniteElement *,2> els(2);
    const FiniteElement & ref_elx = space_x->GetFE( ElementId(VOL,elnums[0]), lh );
    els[0] = &ref_elx;
    const FiniteElement & ref_ely = Space(elnums[0])->GetFE( ElementId(VOL,elnums[1]), lh );
    els[1] = &ref_ely;
    TPHighOrderFE * fe =  new (lh) TPHighOrderFE (els);
    return *fe;
  }
    
  ngfem::ElementTransformation & TPHighOrderFESpace::GetTrafo (ElementId ei, Allocator & lh) const
  {
     TPElementTransformation *trafo = new (lh) TPElementTransformation ( ei );
     ArrayMem<int,2> indices(2);
     GetIndices(ei.Nr(),indices);
     ArrayMem<ElementTransformation *,2> trafos(2);
     ElementId eiix(VOL,indices[0]);
     trafos[0] = &(space_x->GetMeshAccess()->GetTrafo(eiix,lh));
     ElementId eiiy(VOL,indices[1]);
     trafos[1] = &(Space(eiix.Nr())->GetMeshAccess()->GetTrafo(eiiy,lh));
     trafo->SetTrafos(trafos);
     return *trafo;
  }

  const FiniteElement & TPHighOrderFESpace::GetSFE (int elnr, LocalHeap & lh) const
  {
    throw Exception("TPHighOrderFESpace::GetSFE() not implemented");
  }

  const FiniteElement & TPHighOrderFESpace::GetFacetFE (int fnr, LocalHeap & lh) const
  {
    throw Exception("TPHighOrderFESpace::GetFacetFE() not implemented");
  }

  void TPHighOrderFESpace::GetDofRanges (ElementId ei, Array<IntRange> & dranges) const
  {
    throw Exception("TPHighOrderFESpace::GetDofRanges() not implemented");
  }

  void TPHighOrderFESpace::GetDofNrs(ngfem::ElementId ei, ngstd::Array<int>& dnums) const
  {
    dnums.SetSize(0);
    if(ei.VB() != VOL)
      return;
    ArrayMem<int,2> indices;
    ArrayMem<int,100> dnumsx, dnumsy;
    GetIndices(ei.Nr(),indices);
    space_x->GetDofNrs(indices[0],dnumsx);
    Space(indices[0])->GetDofNrs(indices[1],dnumsy);
    dnums.SetSize(dnumsx.Size()*dnumsy.Size());
    int ii=0;
    // new -> working!!!!
    for(int i=0;i<dnumsx.Size();i++)
      for(int j=0;j<dnumsy.Size();j++,ii++)
        dnums[ii] = Space(indices[0])->GetNDof()*dnumsx[i]+dnumsy[j];
  }

  void TPHighOrderFESpace::GetSliceDofNrs(ngfem::ElementId ei, int direction, ngstd::Array<int>& dnums) const
  {
    dnums.SetSize(0);
    if(ei.VB() != VOL)
      return;
    ArrayMem<int,100> dnumsx,dnumsy;
    int totsize = 0;
    if(direction == 1)
    {
      Array<int> alldnums(Space(ei.Nr())->GetNDof());
      space_x->GetDofNrs(ei,dnumsx);
      for(int el=0;el<nels[direction];el++)
      {
        Space(ei.Nr())->GetDofNrs(el,dnumsy);
        alldnums.Range(totsize,totsize+dnumsy.Size()) = dnumsy;
        totsize+=dnumsy.Size();
      }
      dnums.SetSize(dnumsx.Size()*alldnums.Size());
      int ii=0;
      for(int i=0;i<dnumsx.Size();i++)
        for(int j=0;j<alldnums.Size();j++,ii++)
          dnums[ii] = spaces_y[0]->GetNDof()*dnumsx[i]+alldnums[j];
    }
    else
    {
      Array<int> alldnums(space_x->GetNDof());
      spaces_y[0]->GetDofNrs(ei,dnumsy);
      for(int el=0;el<nels[direction];el++)
      {
        space_x->GetDofNrs(el,dnumsx);
        alldnums.Range(totsize,totsize+dnumsx.Size()) = dnumsx;
        totsize+=dnumsx.Size();
      }
      dnums.SetSize(dnumsy.Size()*alldnums.Size());
      int ii=0;
      for(int i=0;i<dnumsy.Size();i++)
        for(int j=0;j<alldnums.Size();j++,ii++)
          dnums[ii] = spaces_y[0]->GetNDof()*alldnums[j] + dnumsy[i];
    }
  }

  
  shared_ptr<Table<int>> TPHighOrderFESpace::CreateSmoothingBlocks (const Flags & precflags) const
  {
    throw Exception("TPHighOrderFESpace::CreateSmoothingBlocks() not implemented");
  }

  void TPHighOrderFESpace::GetVertexDofNrs (int vnr, Array<int> & dnums) const
  {
    throw Exception("TPHighOrderFESpace::GetVertexDofNumbers() not implemented");
  }

  void TPHighOrderFESpace::GetEdgeDofNrs (int ednr, Array<int> & dnums) const
  {
    throw Exception("TPHighOrderFESpace::GetEdgeDofNrs() not implemented");
  }

  void TPHighOrderFESpace::GetFaceDofNrs (int fanr, Array<int> & dnums) const
  {
    throw Exception("TPHighOrderFESpace::GetFaceDofNrs() not implemented");
  }

  void TPHighOrderFESpace::GetInnerDofNrs (int elnr, Array<int> & dnums) const
  {
    throw Exception("TPHighOrderFESpace::GetInnerDofNumbers() not implemented");
  }

  void IterateElementsTP (const FESpace & fes, VorB vb, LocalHeap & clh, 
            const function<void(ElementId,ElementId,LocalHeap&)> & func)
  {
    static mutex copyex_mutex;
    const TPHighOrderFESpace & festp = dynamic_cast<const TPHighOrderFESpace &>(fes);
    shared_ptr<FESpace> space_x = festp.Space(-1);
    shared_ptr<FESpace> space_y = festp.Space(0);
    auto & nels = festp.GetNels();
    const Table<int> & element_coloring0 = space_x->ElementColoring(vb);
    if (task_manager)
    {
      for (FlatArray<int> els_of_col : element_coloring0)
      {
        SharedLoop sl(els_of_col.Range());
        task_manager -> CreateJob
        ( [&] (const TaskInfo & ti) 
        {
          LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
          for (int mynr : sl)
          {
            for(int j=0;j<nels[1];j++)
            {
              HeapReset hr(lh);
              func (ElementId(vb,els_of_col[mynr]),ElementId(vb,j), lh);
            }
          }
        }
        );
      }
      return;
    }
    if(vb == VOL)
      for (int i = 0; i< nels[0];i++)
      {
        for(int j=0;j<nels[1];j++)
        {
          HeapReset hr(clh);
          func (ElementId (vb, i),ElementId (vb, j), clh);
        }
      }
  }
  
  void TPHighOrderFESpace::SolveM (CoefficientFunction & rho, BaseVector & vec,
                         LocalHeap & clh) const
  {
    static Timer tall("TPHighOrderFESpace::SolveM"); RegionTimer rall(tall);
    const Array<shared_ptr<FESpace> > & spaces = Spaces(0);
    int ndofyspace = spaces[1]->GetNDof();
    auto meshx = spaces[0]->GetMeshAccess();
    auto meshy = spaces[1]->GetMeshAccess();
    int nely = meshy->GetNE();
    auto & element_coloring0 = spaces[0]->ElementColoring(VOL);
    for (FlatArray<int> els_of_col : element_coloring0)
    {
      SharedLoop2 sl(els_of_col.Range());
      task_manager -> CreateJob
      ( [&] (const TaskInfo & ti) 
      {
        LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
        for (int mynr : sl)
        {
          HeapReset hr(lh);
          int elnrx = els_of_col[mynr];
          auto & felx = spaces[0]->GetFE(ElementId(elnrx),lh);
          int ndofx = felx.GetNDof();
          const ElementTransformation & xtrafo = meshx->GetTrafo(ElementId(elnrx), lh);
          const IntegrationRule & irx = SelectIntegrationRule(felx.ElementType(),2*felx.Order());
          BaseMappedIntegrationRule & mirx = xtrafo(irx, lh);
          FlatMatrix<> elmat_yslice(ndofx,ndofyspace,lh);
          Array<int> dnums_yslice(ndofx*ndofyspace, lh);
          GetSliceDofNrs(ElementId(elnrx), 1, dnums_yslice);
          vec.GetIndirect (dnums_yslice, elmat_yslice.AsVector());
          FlatMatrix<> shapesx(felx.GetNDof(),irx.Size(),lh);
          shapesx = 0.0;
          dynamic_cast<const BaseScalarFiniteElement &>(felx).CalcShape(irx,shapesx);
          for(int ip=0;ip<irx.Size();ip++)
            shapesx.Col(ip)*=sqrt(mirx[ip].GetWeight());
          FlatMatrix<> elmatx(felx.GetNDof(),lh);
          elmatx = shapesx*Trans(shapesx);
          CalcInverse(elmatx);
          FlatMatrix<> res1(ndofx,ndofyspace,lh);
          FlatMatrix<> out(ndofx,ndofyspace,lh);
          res1 = elmatx*elmat_yslice;
          int firstydof = 0;
          for(int j=0;j<nely;j++)
          {
            auto & fely = spaces[1]->GetFE(ElementId(j),lh);
            int ndofy = fely.GetNDof();
            const ElementTransformation & ytrafo = meshy->GetTrafo(ElementId(j), lh);
            const IntegrationRule & iry = SelectIntegrationRule(fely.ElementType(),2*fely.Order());
            BaseMappedIntegrationRule & miry = ytrafo(iry, lh);
            FlatMatrix<> shapesy(fely.GetNDof(),iry.Size(),lh);
            shapesy = 0.0;
            dynamic_cast<const BaseScalarFiniteElement &>(fely).CalcShape(iry,shapesy); 
            for(int ip=0;ip<iry.Size();ip++)
              shapesy.Col(ip)*=sqrt(miry[ip].GetWeight());
            FlatMatrix<> elmaty(fely.GetNDof(),lh);
            elmaty = shapesy*Trans(shapesy);
            CalcInverse(elmaty);
            IntRange dnumsy(firstydof, firstydof+ndofy);
            firstydof+=ndofy;
            out.Cols(dnumsy) = res1.Cols(dnumsy)*Trans(elmaty);
          }
          vec.SetIndirect(dnums_yslice,out.AsVector());
        }
      });
    }
  }

  void Transfer2StdMesh(const GridFunction * gfutp, GridFunction* gfustd)
  {
    static Timer tall("TPHighOrderFESpace::Transfer2StdMesh"); RegionTimer rall(tall);
    const shared_ptr<FESpace> fes = gfustd->GetFESpace();
    const shared_ptr<TPHighOrderFESpace> tpfes = dynamic_pointer_cast<TPHighOrderFESpace>(gfutp->GetFESpace());
    shared_ptr<FESpace> fesx = tpfes->Space(-1);
    shared_ptr<FESpace> fesy = tpfes->Space(0);
    LocalHeap lh(100000000,"heap");
    auto & els = dynamic_cast<TPHighOrderFE &> (tpfes->GetFE(ElementId(0),lh));
    IterateElementsTP(*tpfes,VOL,lh,
    [&] (ElementId ei0,ElementId ei1,LocalHeap & lh)
    {
      HeapReset hr(lh);
      ArrayMem<int,2> ind(2);
      ind[0] = ei0.Nr(); ind[1] = ei1.Nr();
      int elnr = tpfes->GetIndex(ind);
      TPHighOrderFE & tpfel = dynamic_cast<TPHighOrderFE&>(tpfes->GetFE(ElementId(elnr),lh));
      int elnrstd = elnr;
      const FiniteElement & fel = fes->GetFE(ElementId(elnrstd),lh);
      Array<const IntegrationRule *> irs(tpfel.elements.Size());
      for(int s=0;s<irs.Size();s++)
        irs[s] = &SelectIntegrationRule(tpfel.elements[s]->ElementType(),2*fel.Order());
      TPIntegrationRule ir(irs);
      const ElementTransformation & tptrafo = tpfes->GetTrafo(ElementId(elnr),lh);
      BaseMappedIntegrationRule & tpmir = tptrafo(ir, lh);
      IntegrationRule irstd;
      for(int s=0;s<ir(0).Size();s++)
        for(int t=0;t<ir(1).Size();t++)
        {
          if(fesx->GetMeshAccess()->GetDimension() == 1 )
              irstd.AddIntegrationPoint(IntegrationPoint(ir(0)[s](0),ir(1)[t](0),ir(1)[t](1), ir(0)[s].Weight()*ir(1)[t].Weight()));
          if(fesx->GetMeshAccess()->GetDimension() == 2 )
              irstd.AddIntegrationPoint(IntegrationPoint(ir(0)[s](0),ir(0)[s](1),ir(1)[t](0), ir(0)[s].Weight()*ir(1)[t].Weight()));
         }
      const ElementTransformation & trafo = fes->GetMeshAccess()->GetTrafo(ElementId(VOL,elnrstd),lh);
      BaseMappedIntegrationRule & mirstd = trafo(irstd,lh);
      shared_ptr<TPDifferentialOperator> diffop = dynamic_pointer_cast<TPDifferentialOperator>(tpfes->GetEvaluator());
      int niptp = ir(0).Size()*ir(1).Size();
      // Evalute \int u_tp * v_std :
      FlatMatrix<> flux(niptp,diffop->Dim(),lh);
      FlatVector<> coef(tpfel.GetNDof(),lh);
      const BaseVector & base = gfutp->GetVector();
      Array<int> dnums;
      tpfes->GetDofNrs(elnr,dnums);
      base.GetIndirect(dnums,coef);
      diffop->Apply(tpfel, tpmir, coef,flux,lh);
      // Build mass matrix for \int u_std * v_std :
      FlatMatrix<> elmat(fel.GetNDof(),lh);
      elmat = 0.0;
      FlatMatrix<> shapes(fel.GetNDof(),mirstd.Size(),lh);
      shapes = 0.0;
      dynamic_cast<const BaseScalarFiniteElement &>(fel).CalcShape(irstd,shapes);
      FlatMatrix<> shapes1(fel.GetNDof(),mirstd.Size(),lh);
        
      FlatVector<> elvec(fel.GetNDof(),lh),y(fel.GetNDof(),lh);
      shapes1 = shapes;
      for(int s=0;s<shapes.Width();s++)
        shapes.Col(s) *= mirstd[s].GetWeight();
      elmat = shapes1*Trans(shapes);
      elvec = shapes * flux;
      CalcInverse(elmat);
      y = elmat*elvec;
      BaseVector & baseout = gfustd->GetVector();
      fes->GetDofNrs(ElementId(VOL,elnrstd),dnums);
      baseout.SetIndirect(dnums,y);    
    });
  }

  void Transfer2TPMesh(const CoefficientFunction * cfstd, GridFunction* gfutp)
  {
    static Timer tall("TPHighOrderFESpace::Transfer2TPMesh"); RegionTimer rall(tall);
    const shared_ptr<TPHighOrderFESpace> tpfes = dynamic_pointer_cast<TPHighOrderFESpace>(gfutp->GetFESpace());
    LocalHeap lh(100000000,"heap");
    IterateElementsTP(*tpfes,VOL,lh,
    [&] (ElementId ei0,ElementId ei1,LocalHeap & lh)
    {
        HeapReset hr(lh);
        Array<int> ind(2);
        
        ind[0] = ei0.Nr(); ind[1] = ei1.Nr();
        int elnr = tpfes->GetIndex(ind);
        
        TPHighOrderFE & tpfel = dynamic_cast<TPHighOrderFE&>(tpfes->GetFE(ElementId(elnr),lh));
        ArrayMem<const IntegrationRule * , 2> irs(tpfel.elements.Size());
        for(int s=0;s<irs.Size();s++)
          irs[s] = &SelectIntegrationRule(tpfel.elements[s]->ElementType(),2*tpfel.elements[s]->Order());
        TPIntegrationRule ir(irs);
        const ElementTransformation & tptrafo = tpfes->GetTrafo(ElementId(elnr),lh);
        TPMappedIntegrationRule & tpmir = dynamic_cast<TPMappedIntegrationRule & >(tptrafo(ir, lh));
        int tpnip = irs[0]->Size()*irs[1]->Size();
        FlatMatrix<> result(tpnip,1,lh);
        cfstd->Evaluate(tpmir,result);
        FlatMatrix<> shapes(tpfel.GetNDof(),tpnip,lh);
        shapes = 0.0;
        tpfel.CalcShape(ir,shapes);
        FlatVector<> elvecy(tpfel.GetNDof(),lh), elvecx(tpfel.GetNDof(),lh);
        FlatMatrix<> elmat(tpfel.GetNDof(),lh);
        elmat = 0.0;
        FlatMatrix<> shapes1(tpfel.GetNDof(),tpnip,lh);
        shapes1 = shapes;
        TPMappedIntegrationRule & ttpmir = dynamic_cast<TPMappedIntegrationRule &>(tpmir);        
        for(int i=0,ii=0;i<ttpmir.GetIRs()[0]->Size();i++)
          for(int j=0;j<ttpmir.GetIRs()[1]->Size();j++,ii++)
            shapes.Col(ii) *= (*ttpmir.GetIRs()[0])[i].GetWeight()*(*ttpmir.GetIRs()[1])[j].GetWeight();
        elmat = shapes1*Trans(shapes);
        elvecx = shapes*result;
        CalcInverse(elmat);
        elvecy = elmat * elvecx;
        Array<int> dnums;
        tpfes->GetDofNrs(elnr,dnums);
        gfutp->GetVector().SetIndirect(dnums,elvecy);    
    });
  }
  
  }
