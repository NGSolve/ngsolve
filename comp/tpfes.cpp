#include <comp.hpp>
#include "tpfes.hpp"
#include "./../fem/hermiteelementtransformation.hpp"
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
      if( el_counts!=nullptr )
      {
        nels[i] = (*el_counts)[i];
        nfacets[i] = (*el_counts)[i+2];
      }
      else
      {
        nels[i] = meshes[i]->GetNE();
        nfacets[i] = meshes[i]->GetNFacets();
      }
      ndof *= ndofs[i];
      nel *= nels[i];
    }
    
    first_element_dof.SetSize(nel+1);
    
    LocalHeap lh(100000,"Setup TP Space");    
    int ii=0;
    first_element_dof[0] = 0;
    for(int elx : Range(nels[0]) )
    {
      int ndofx = space_x->GetFE(elx,lh).GetNDof();
      for(int ely : Range(nels[1]) )
      {
        HeapReset hr(lh);
        int ndofy = Space(elx)->GetFE(ely,lh).GetNDof();
        first_element_dof[ii+1] = first_element_dof[ii] + ndofx*ndofy;
        ii++;
      }
    }
    Array<shared_ptr<DifferentialOperator>> evaluators(nmeshes);
    for (int i : Range(nmeshes))
      evaluators[i] = spaces[i]->GetEvaluator();
    evaluator = shared_ptr<DifferentialOperator>( new TPDifferentialOperator(evaluators) );
  }

  TPHighOrderFESpace::TPHighOrderFESpace (shared_ptr<FESpace> aspace_x,FlatArray<shared_ptr<FESpace>> aspaces_y, const Flags & flags, bool parseflags, Array<int> * el_counts)
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
    for (int i : Range(nmeshes)) {
      if( el_counts )
      {
        nels[i] = (*el_counts)[i];
        nfacets[i] = (*el_counts)[i+2];        
      }
      else
      {
        nels[i] = meshes[i]->GetNE();
        nfacets[i] = meshes[i]->GetNFacets();
      }
      nel *= nels[i];
    }
    ndof = 0;
    LocalHeap lh(100000,"Setup TP Space");
    for(int i=0;i<space_x->GetMeshAccess()->GetNE();i++)
      ndof += fespaces[0]->GetFE(i,lh).GetNDof()*spaces_y[i]->GetNDof();
    
    first_element_dof.SetSize(nel+1);
    int ii=0;
    first_element_dof[0] = 0;
    for(int elx : Range(nels[0]) )
    {
      int ndofx = space_x->GetFE(elx,lh).GetNDof();
      for(int ely : Range(nels[1]) )
      {
        HeapReset hr(lh);
        int ndofy = Space(elx)->GetFE(ely,lh).GetNDof();
        first_element_dof[ii+1] = first_element_dof[ii] + ndofx*ndofy;
        ii++;
      }
    }
    
    Array<shared_ptr<DifferentialOperator>> evaluators(nmeshes);
    
    evaluators[0] = space_x->GetEvaluator();
    evaluators[1] = spaces_y[0]->GetEvaluator();
    evaluator = shared_ptr<DifferentialOperator>( new TPDifferentialOperator(evaluators) );
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
    ops.Set("gradx", make_shared<TPDifferentialOperator>( gradx ));
    ops.Set("grady", make_shared<TPDifferentialOperator>( grady ));
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
      FlatVector<> elvec_out(space_x->GetFE(i,clh).GetNDof(),clh);
      elvec_out = 0.0;
      for(int j=0;j<nels[1];j++)
      {
        HeapReset hr(clh);
        int index = GetIndex(i,j);
        GetDofNrs(index,dnums);
        FlatVector<> elvec(dnums.Size(),clh);
        vec_in.GetIndirect(dnums,elvec);
        const TPHighOrderFE & tpfel = dynamic_cast<const TPHighOrderFE &>(GetFE(index,clh));
        if(tpfel.elements[1]->ElementType() == ET_HERMITE)
        {
          ElementTransformation & trafo = dynamic_cast<TPElementTransformation &>(GetTrafo(i,clh)).GetTrafo(1);
          func(gf_in->GetFESpace(),tpfel,trafo,elvec,elvec_out,clh);
        }
        else
        {
          ElementTransformation & trafo = fespaces[1]->GetMeshAccess()->GetTrafo(j,clh);
          func(gf_in->GetFESpace(),tpfel,trafo,elvec,elvec_out,clh);
        }
      }
      space_x->GetDofNrs(i,dnumsx);
      gf_out->GetVector().SetIndirect(dnumsx,elvec_out);
    }
    
    // shared_ptr<FESpace> space_x = Space(-1);
    // shared_ptr<FESpace> space_y = Space(0);
    // const Table<int> & element_coloring0 = space_x->ElementColoring(VOL);
    // Array<int> dnums,dnumsx;
    // if (task_manager)
    // {
      // for (FlatArray<int> els_of_col : element_coloring0)
      // {
        // SharedLoop sl(els_of_col.Range());
        // task_manager -> CreateJob
        // ( [&] (const TaskInfo & ti) 
        // {
          // LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
          // for (int mynr : sl)
          // {
            // FlatVector<> elvec_out(space_x->GetFE(els_of_col[mynr],lh).GetNDof(),lh);
            // elvec_out = 0.0;
            // for(int j=0;j<nels[1];j++)
            // {
              // HeapReset hr(lh);
              // int index = GetIndex(els_of_col[mynr],j);
              // GetDofNrs(index,dnums);
              // FlatVector<> elvec(dnums.Size(),lh);
              // vec_in.GetIndirect(dnums,elvec);
              // const TPHighOrderFE & tpfel = dynamic_cast<const TPHighOrderFE &>(GetFE(index,lh));
              // ElementTransformation & trafo = fespaces[1]->GetMeshAccess()->GetTrafo(j,lh);
              // func(gf_in->GetFESpace(),tpfel,trafo,elvec,elvec_out,lh);
            // }
            // space_x->GetDofNrs(els_of_col[mynr],dnumsx);
            // gf_out->GetVector().SetIndirect(dnumsx,elvec_out);
          // }
        // });
      // }
      // return;
    // }
    
    
    
    
  }

  void TPHighOrderFESpace :: ProlongateFromXSpace(shared_ptr<GridFunction> gf_in, shared_ptr<GridFunction> gf_out,LocalHeap & lh)
  {
    BaseVector & vec_in = gf_in->GetVector();
    BaseVector & vec_out = gf_out->GetVector();
    Array<int> dnums,dnumsx;
    for(int i=0;i<nels[0];i++)
    {
      Vector<> elvec_in(space_x->GetFE(i,lh).GetNDof());
      fespaces[0]->GetDofNrs(i,dnumsx);
      vec_in.GetIndirect(dnumsx,elvec_in);      
      //indices[0] = i;
      for(int j=0;j<nels[1];j++)
      {
        HeapReset hr(lh);
        //indices[1] = j;
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
    element_coloring = Table<int>(nel,1);
    for (int i : Range(nel))
      element_coloring[i][0] = i;
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
  
  int TPHighOrderFESpace::GetNDof () const throw() 
  {
    return ndof;
  }
  
  int TPHighOrderFESpace::GetNDofLevel (int level) const 
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
    const FiniteElement & ref_elx = space_x->GetFE( elnums[0], lh );
    els[0] = &ref_elx;
    const FiniteElement & ref_ely = Space(elnums[0])->GetFE( elnums[1], lh );
    els[1] = &ref_ely;
    TPHighOrderFE * fe =  new (lh) TPHighOrderFE (els);
    return *fe;
  }
    
  ngfem::ElementTransformation & TPHighOrderFESpace::GetTrafo (ElementId ei, Allocator & lh) const
  {
     TPElementTransformation *trafo = new(lh) TPElementTransformation ( ei );
     ArrayMem<int,2> indices(2);
     GetIndices(ei.Nr(),indices);
     ArrayMem<ElementTransformation *,2> trafos(2);
     ElementId eiix(VOL,indices[0]);
     trafos[0] = &(space_x->GetMeshAccess()->GetTrafo(eiix,lh));
     ElementId eiiy(VOL,indices[1]);
     if( Space(eiix.Nr())->GetFE(0,lh).ElementType() == ET_HERMITE )
     {
//       if(Space(eiix.Nr())->GetSpacialDimension() == 1)
         trafos[1] = new (lh) HM_ElementTransformation(dynamic_pointer_cast<GlobalHermiteFESpace>(spaces_y[0])->AnsatzTemp(eiix.Nr()),&dynamic_pointer_cast<GlobalHermiteFESpace>(spaces_y[0])->AnsatzV(eiix.Nr()));
//       if(Space(eiix.Nr())->GetSpacialDimension() == 2)
//         trafos[1] = new (lh) HM_ElementTransformation(1.0,Vec<2>(0.0,0.0));
//       if(Space(eiix.Nr())->GetSpacialDimension() == 3)
//         trafos[1] = new (lh) HM_ElementTransformation(1.0,Vec<3>(0.0,0.0,0.0));
     }
     else
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

  void TPHighOrderFESpace::GetDofNrs (int elnr, Array<int> & dnums) const 
  {
    dnums.SetSize(first_element_dof[elnr+1]-first_element_dof[elnr]);
    for(int i=0;i<first_element_dof[elnr+1]-first_element_dof[elnr];i++)
      dnums[i] = first_element_dof[elnr]+i;
  }

  void TPHighOrderFESpace::GetSDofNrs (int selnr, Array<int> & dnums) const
  {
    dnums.SetSize(0);
  }

  Table<int> * TPHighOrderFESpace::CreateSmoothingBlocks (const Flags & precflags) const
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

//  void IterateElementsTP (const FESpace & fes, VorB vb, LocalHeap & clh, 
//			const function<void(FESpace::Element,FESpace::Element,LocalHeap&)> & func)
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
    // else if(vb == BND)
    // {
      // for (int i = 0; i< ma0->GetNSE();i++)
      // {
        // Array<int> temp_dnums0;
        // FESpace::Element el0(*spaces[0], ElementId (BND, i), temp_dnums0, clh);
        // for(int j=0;j<ma1->GetNE();j++)
        // {
          // Array<int> temp_dnums1;
          // HeapReset hr(clh);
          // FESpace::Element el1(*spaces[1], ElementId (VOL, j), temp_dnums1, clh);
          // func (move(el0),move(el1), clh);
        // }
      // }
      // for (int i = 0; i< ma0->GetNE();i++)
      // {
        // Array<int> temp_dnums0;
        // FESpace::Element el0(*spaces[0], ElementId (VOL, i), temp_dnums0, clh);
        // for(int j=0;j<ma1->GetNSE();j++)
        // {
          // Array<int> temp_dnums1;
          // HeapReset hr(clh);
          // FESpace::Element el1(*spaces[1], ElementId (BND, j), temp_dnums1, clh);
          // func (move(el0),move(el1), clh);
        // }
      // }
    // }
  }
  
  void TPHighOrderFESpace::SolveM (CoefficientFunction & rho, BaseVector & vec,
                         LocalHeap & lh) const
  {
    IterateElementsTP(*this, VOL, lh, 
      [&] (ElementId ei0,ElementId ei1, LocalHeap & lh)
      {
        Array<int> indices(2);
        indices[0] = ei0.Nr(); indices[1] = ei1.Nr();
        int elnr = GetIndex(indices);
        const Array<shared_ptr<FESpace> > & spaces = Spaces(ei0.Nr());
        TPElementTransformation & eltrans = dynamic_cast<TPElementTransformation &>(GetTrafo (elnr, lh));
        TPHighOrderFE & tpfel = dynamic_cast<TPHighOrderFE&>(GetFE(ElementId(elnr),lh));
        IntegrationRule irx(tpfel.elements[0]->ElementType(),2*tpfel.elements[0]->Order());
        IntegrationRule iry(tpfel.elements[1]->ElementType(),2*tpfel.elements[1]->Order());
        BaseMappedIntegrationRule & mirx = eltrans.GetTrafo(0)(irx,lh);
        BaseMappedIntegrationRule & miry = eltrans.GetTrafo(1)(iry,lh);
        FlatMatrix<> shapesx(tpfel.elements[0]->GetNDof(),irx.Size(),lh);
        FlatMatrix<> shapesy(tpfel.elements[1]->GetNDof(),iry.Size(),lh);
        shapesx = 0.0;
        shapesy = 0.0;
        dynamic_cast<const BaseScalarFiniteElement *>(tpfel.elements[0])->CalcShape(irx,shapesx);
        dynamic_cast<const BaseScalarFiniteElement *>(tpfel.elements[1])->CalcShape(iry,shapesy);
        FlatMatrix<> elmatx(tpfel.elements[0]->GetNDof(),lh);
        FlatMatrix<> elmaty(tpfel.elements[1]->GetNDof(),lh);
        for(int i=0;i<shapesx.Width();i++)
          shapesx.Col(i)*=sqrt(mirx[i].GetWeight());
        for(int i=0;i<shapesy.Width();i++)
          shapesy.Col(i)*=sqrt(miry[i].GetWeight());
        elmatx = shapesx*Trans(shapesx);
        elmaty = shapesy*Trans(shapesy);
        CalcInverse(elmatx);
        CalcInverse(elmaty);
        FlatVector<> elvecx(tpfel.GetNDof(),lh),elvecy(tpfel.GetNDof(),lh);
        Array<int> dnums;
        GetDofNrs (elnr, dnums);
        vec.GetIndirect(dnums,elvecx);
        FlatMatrix<> umat(tpfel.elements[0]->GetNDof(),tpfel.elements[1]->GetNDof(),&elvecx(0));
        FlatMatrix<> helper(tpfel.elements[0]->GetNDof(),tpfel.elements[1]->GetNDof(),lh);
        FlatMatrix<> helper2(tpfel.elements[0]->GetNDof(),tpfel.elements[1]->GetNDof(),&elvecy(0));
        helper = elmatx*umat;
        helper2 = helper*Trans(elmaty);
        vec.SetIndirect(dnums,elvecy);
      });
      
  }

  void Transfer2StdMesh(const GridFunction * gfutp, GridFunction* gfustd)
  {
    const shared_ptr<FESpace> fes = gfustd->GetFESpace();
    const shared_ptr<TPHighOrderFESpace> tpfes = dynamic_pointer_cast<TPHighOrderFESpace>(gfutp->GetFESpace());
    shared_ptr<FESpace> fesx = tpfes->Space(-1);
    shared_ptr<FESpace> fesy = tpfes->Space(0);
    LocalHeap lh(100000000,"heap");
    auto & els = dynamic_cast<TPHighOrderFE &> (tpfes->GetFE(ElementId(0),lh));
    if(els.elements[1]->ElementType() == ET_HERMITE)
    {
      IterateElementsTP(*tpfes,VOL,lh,
      [&] (ElementId ei0,ElementId ei1,LocalHeap & lh)
      {
        HeapReset hr(lh);
        Array<int> ind(2);
        ind[0] = ei0.Nr(); ind[1] = ei1.Nr();
        int elnr = tpfes->GetIndex(ind);
        TPHighOrderFE & tpfel = dynamic_cast<TPHighOrderFE&>(tpfes->GetFE(ElementId(elnr),lh));
        const shared_ptr<MeshAccess> helperma = dynamic_pointer_cast<GlobalHermiteFESpace>(fesy)->HelperSpace()->GetMeshAccess();
        for(int i = 0;i < helperma->GetNE() ; i++)
        {
          int elnrstd = elnr*helperma->GetNE() + i;
          const FiniteElement & fel = fes->GetFE(elnrstd,lh);
          Array<const IntegrationRule *> irs(tpfel.elements.Size());
          irs[0] = &SelectIntegrationRule(tpfel.elements[0]->ElementType(),2*fel.Order());
          const FiniteElement & felyhelper = dynamic_pointer_cast<GlobalHermiteFESpace> (fesy)->HelperSpace()->GetFE(ElementId(VOL,i),lh);
          irs[1] = &SelectIntegrationRule(felyhelper.ElementType(),2*max(tpfel.elements[1]->Order(),fel.Order()) );
          const ElementTransformation & tptrafo = tpfes->GetTrafo(ElementId(VOL,elnr),lh);
          shared_ptr<TPDifferentialOperator> diffop = dynamic_pointer_cast<TPDifferentialOperator>(tpfes->GetEvaluator());
          int niptp = irs[0]->Size()*irs[1]->Size();
          FlatMatrix<> elmat(fel.GetNDof(), lh);
          elmat = 0.0;
          FlatMatrix<> shapes(fel.GetNDof(),niptp,lh);
          FlatVector<> elvec(fel.GetNDof(),lh),y(fel.GetNDof(),lh);
          IntegrationRule irstd;
          IntegrationRule helperirv;
          for(int s=0;s<irs[0]->Size();s++)
            for(int t=0;t<irs[1]->Size();t++)
            {
              if(fesx->GetMeshAccess()->GetDimension() == 1 )
                  irstd.AddIntegrationPoint(IntegrationPoint((*irs[0])[s](0),(*irs[1])[t](0),(*irs[1])[t](1), (*irs[0])[s].Weight()*(*irs[1])[t].Weight()));
              if(fesx->GetMeshAccess()->GetDimension() == 2 )
                  irstd.AddIntegrationPoint(IntegrationPoint((*irs[0])[s](0),(*irs[0])[s](1),(*irs[1])[t](0), (*irs[0])[s].Weight()*(*irs[1])[t].Weight()));
            }
          
          for(int t=0;t<irs[1]->Size();t++)
            helperirv.AddIntegrationPoint( (helperma->GetTrafo(ElementId(VOL,i),lh)(*irs[1],lh))[t].GetPoint()(0) );
          irs[1] = &helperirv;
          TPIntegrationRule iir(irs);
          TPMappedIntegrationRule & ttpmir = dynamic_cast<TPMappedIntegrationRule &>(tptrafo(iir,lh));
          FlatVector<> flux(niptp,lh);
          
          FlatVector<> coef(tpfel.GetNDof(),lh);
          Array<int> dnums;
          tpfes->GetDofNrs(elnr,dnums);
          gfutp->GetVector().GetIndirect(dnums,coef);
          FlatMatrix<> shape0(tpfel.elements[0]->GetNDof(), irs[0]->Size(),lh);
          FlatMatrix<> coefmat(tpfel.elements[0]->GetNDof(),tpfel.elements[1]->GetNDof(),&coef(0));
          dynamic_cast<const BaseScalarFiniteElement *>(tpfel.elements[0])->CalcShape(*irs[0],shape0);
          FlatMatrix<> coefsipx(irs[0]->Size(),tpfel.elements[1]->GetNDof(), lh);
          coefsipx = Trans(shape0)*coefmat;
          if(tpfes->Space(1)->GetSpacialDimension() == 1)
          {
            double Tref = dynamic_cast<HM_ElementTransformation &>(dynamic_cast<const TPElementTransformation &>(tptrafo).GetTrafo(1)).GetAnsatzT();
            Vec<1> Vref = dynamic_cast<HM_ElementTransformation &>(dynamic_cast<const TPElementTransformation &>(tptrafo).GetTrafo(1)).GetAnsatzV();
            for(int u=0;u<irs[0]->Size();u++)
              dynamic_cast<const Distribution<1> *>(tpfel.elements[1])->Evaluate<NODAL>(coefsipx.Row(u),helperirv,flux.Range(u*irs[1]->Size(),(u+1)*irs[1]->Size()) ,Tref,Vref);
          }
          if(tpfes->Space(1)->GetSpacialDimension() == 2)
          {
            double Tref = dynamic_cast<HM_ElementTransformation &>(dynamic_cast<const TPElementTransformation &>(tptrafo).GetTrafo(1)).GetAnsatzT();
            Vec<2> Vref = dynamic_cast<HM_ElementTransformation &>(dynamic_cast<const TPElementTransformation &>(tptrafo).GetTrafo(1)).GetAnsatzV();
            for(int u=0;u<irs[0]->Size();u++)
              dynamic_cast<const Distribution<2> *>(tpfel.elements[1])->Evaluate<NODAL>(coefsipx.Row(u),helperirv,flux.Range(u*irs[1]->Size(),(u+1)*irs[1]->Size()) ,Tref,Vref);
          }          
          const ElementTransformation & trafo = fes->GetMeshAccess()->GetTrafo(ElementId(VOL,elnrstd),lh);
          BaseMappedIntegrationRule & mirstd = trafo(irstd,lh);
          shapes = 0.0;
          dynamic_cast<const BaseScalarFiniteElement &>(fel).CalcShape(irstd,shapes);
          FlatMatrix<> shapes1(fel.GetNDof(),mirstd.Size(),lh);
          shapes1 = shapes;
          for(int s=0;s<shapes.Width();s++)
            shapes.Col(s)*= mirstd[s].GetWeight();
          elmat = shapes1*Trans(shapes);
          elvec = shapes * flux;
          fes->GetDofNrs(elnrstd,dnums);
          CalcInverse(elmat);
          y = elmat*elvec;
          BaseVector & baseout = gfustd->GetVector();
          baseout.SetIndirect(dnums,y);
        }
      });
    }
    else
      IterateElementsTP(*tpfes,VOL,lh,
      [&] (ElementId ei0,ElementId ei1,LocalHeap & lh)
      {
        HeapReset hr(lh);
        Array<int> ind(2);
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
        fes->GetDofNrs(elnrstd,dnums);
        baseout.SetIndirect(dnums,y);    
      });
  }

  void Transfer2TPMesh(const CoefficientFunction * cfstd, GridFunction* gfutp)
  {
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
        if(tpfel.elements[1]->ElementType() == ET_HERMITE)
          irs[1] = &SelectIntegrationRule(tpfel.elements[1]->ElementType(),tpfel.elements[1]->Order()+1);
        TPIntegrationRule ir(irs);
        const ElementTransformation & tptrafo = tpfes->GetTrafo(ElementId(elnr),lh);
        TPMappedIntegrationRule & tpmir = dynamic_cast<TPMappedIntegrationRule & >(tptrafo(ir, lh));
        int tpnip = irs[0]->Size()*irs[1]->Size();
        FlatMatrix<> result(tpnip,1,lh);
        cfstd->Evaluate(tpmir,result);
        if(tpfel.elements[1]->ElementType() == ET_HERMITE)
        {
          bool hm;
          if(dynamic_pointer_cast<GlobalHermiteFESpace>(tpfes->Space(0))->GetSpacialDimension() == 1)
            hm = dynamic_cast<const Distribution<1> *>(tpfel.elements[1])->UseHermiteFunctions();
          if(dynamic_pointer_cast<GlobalHermiteFESpace>(tpfes->Space(0))->GetSpacialDimension() == 2)
            hm = dynamic_cast<const Distribution<2> *>(tpfel.elements[1])->UseHermiteFunctions();            
          if(!hm)
            for(int i=0,ii=0;i<tpmir.GetIRs()[0]->Size();i++)
              for(int j=0;j<tpmir.GetIRs()[1]->Size();j++,ii++)
                result(ii) *=exp( 0.5*L2Norm2( (*irs[1])[j].Point() ) );
        }
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