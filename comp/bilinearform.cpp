#include <comp.hpp>
#include <multigrid.hpp>

#include <parallelngs.hpp>
#include "../fem/h1lofe.hpp"

namespace ngcomp
{
  
  // dummy function header 
  void CalcEigenSystem (FlatMatrix<Complex> & elmat, 
                        FlatVector<Complex> & lami, 
                        FlatMatrix<Complex> & evecs)
  { ; }

  void MinusMultAB (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
  {
    c = -a * b | Lapack;
  }
  void MinusMultABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
  {
    c = -a * Trans(b) | Lapack;
  }
  void AddAB (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
  {
    c += a*b | Lapack;
  }

  template <typename T>
  inline void AInvBt (ngbla::FlatMatrix<T> a, ngbla::FlatMatrix<T> b)
  {
    // static Timer t("LapackAinvBt");
    // RegionTimer reg(t);
    LapackAInvBt (a, b, 'N');
  }

  // b <-- b A^{-1} 
  inline void AInvBt (ngbla::FlatMatrix<double> a, ngbla::FlatMatrix<double> b)
  {
    /*
    Matrix ha = a;
    Matrix hb = b;
    LapackAInvBt(ha, hb, 'N');
    */
    // static Timer t("AinvBt - new LU");
    // RegionTimer reg(t);    

    ArrayMem<int, 100> p(a.Height());
    CalcLU (a, p);
    SolveTransFromLU (a, p, Trans(b));
    
    // cout << "diff = " << L2Norm(hb-b) << endl;
  }
  
  template <typename SCAL>
  void CalcSchur (FlatMatrix<SCAL> a, 
                  FlatMatrix<SCAL> schur,
                  FlatArray<int> keep, FlatArray<int> eliminate)
  {
    // via Cholesky
    // A = L D L^t
    int nk = keep.Size();
    int ne = eliminate.Size();
    int n = nk+ne;
    Matrix<SCAL> temp(n);

    temp.Rows(0,ne).Cols(0,ne) = a.Rows(eliminate).Cols(eliminate);
    temp.Rows(0,ne).Cols(ne,n) = a.Rows(eliminate).Cols(keep);
    temp.Rows(ne,n).Cols(0,ne) = a.Rows(keep).Cols(eliminate);
    temp.Rows(ne,n).Cols(ne,n) = a.Rows(keep).Cols(keep);

    for (int i = 0; i < ne; i++)
      {
        SCAL d = temp(i,i);
        SCAL id = 1.0/d;

        for (int j = i+1; j < n; j++)
          for (int k = i+1; k < n; k++)
            temp(j,k) -= temp(j,i) * id * temp(i,k);
        for (int j = i+1; j < n; j++)
          temp(j,i) *= id;
      }
    schur = temp.Rows(ne,n).Cols(ne,n);
  }






  BilinearForm :: 
  BilinearForm (shared_ptr<FESpace> afespace,
                const string & aname,
                const Flags & flags)
    : NGS_Object(afespace->GetMeshAccess(), flags, aname), fespace(afespace)
  {
    fespace2 = NULL;

    multilevel = true;
    symmetric = flags.GetDefineFlag ("symmetric");

    // low_order_bilinear_form = NULL;
    linearform = NULL;


    SetGalerkin( flags.GetDefineFlag( "project" ) );
    SetNonAssemble (flags.GetDefineFlag ("nonassemble"));
    SetDiagonal (flags.GetDefineFlag ("diagonal"));
    if (flags.GetDefineFlag ("nonsym"))  SetSymmetric (0);
    if (flags.GetDefineFlag ("nonmultilevel")) SetMultiLevel (0);
    SetHermitean (flags.GetDefineFlag ("hermitean"));
    SetHermitean (flags.GetDefineFlag ("hermitian"));
    SetUnusedDiag (flags.GetNumFlag ("unuseddiag",0.0));
    SetEpsRegularization (flags.GetNumFlag ("regularization",0));
  
    SetPrint (flags.GetDefineFlag ("print"));
    SetPrintElmat (flags.GetDefineFlag ("printelmat"));
    SetElmatEigenValues (flags.GetDefineFlag ("elmatev")); 
    SetTiming (flags.GetDefineFlag ("timing"));
    SetEliminateInternal (flags.GetDefineFlag ("eliminate_internal") || flags.GetDefineFlag("condense"));
    SetEliminateHidden (flags.GetDefineFlag ("eliminate_hidden"));
    SetKeepInternal (eliminate_internal &&
                     !flags.GetDefineFlagX ("keep_internal").IsFalse() &&
                     !flags.GetDefineFlag ("nokeep_internal"));
    SetStoreInner (flags.GetDefineFlag ("store_inner"));
    precompute = flags.GetDefineFlag ("precompute");
    checksum = flags.GetDefineFlag ("checksum");
    spd = flags.GetDefineFlag ("spd");
    geom_free = flags.GetDefineFlag("geom_free");
    matrix_free_bdb = flags.GetDefineFlag("matrix_free_bdb");    
    nonlinear_matrix_free_bdb = flags.GetDefineFlag("nonlinear_matrix_free_bdb");    
    if (spd) symmetric = true;
    SetCheckUnused (!flags.GetDefineFlagX("check_unused").IsFalse());
  }


  BilinearForm :: 
  BilinearForm (shared_ptr<FESpace> afespace,
                shared_ptr<FESpace> afespace2, const string & aname,
                const Flags & flags)
    : NGS_Object(afespace->GetMeshAccess(), flags, aname), fespace(afespace), fespace2(afespace2)
  {
    if (fespace->GetMeshAccess() != fespace2->GetMeshAccess())
      throw Exception("Trial- and test-spaces must be defined on the same mesh");
    
    multilevel = true;
    galerkin = false;
    symmetric = false;
    spd = false;
    hermitean = false;

    SetEpsRegularization (0);

    // low_order_bilinear_form = NULL;
    linearform = NULL;

    timing = false;
    print = false;
    printelmat = false;
    elmat_ev = false;
    eliminate_internal = false;
    eliminate_hidden = false;
    // keep_internal = false;


    SetGalerkin( flags.GetDefineFlag( "project" ) );
    SetNonAssemble (flags.GetDefineFlag ("nonassemble"));
    SetDiagonal (flags.GetDefineFlag ("diagonal"));
    if (flags.GetDefineFlag ("nonsym"))  SetSymmetric (0);
    if (flags.GetDefineFlag ("nonmultilevel")) SetMultiLevel (0);
    SetHermitean (flags.GetDefineFlag ("hermitean"));
    SetHermitean (flags.GetDefineFlag ("hermitian"));
    SetUnusedDiag (flags.GetNumFlag ("unuseddiag",0.0));
  
    SetPrint (flags.GetDefineFlag ("print"));
    SetPrintElmat (flags.GetDefineFlag ("printelmat"));
    SetElmatEigenValues (flags.GetDefineFlag ("elmatev"));

    if (flags.GetDefineFlag ("timing")) SetTiming (1);
    if (flags.GetDefineFlag ("eliminate_internal") || flags.GetDefineFlag("condense")) SetEliminateInternal (1);
    if (flags.GetDefineFlag ("eliminate_hidden")) SetEliminateHidden (1);
    SetKeepInternal (eliminate_internal && 
                     !flags.GetDefineFlag ("nokeep_internal"));
    if (flags.GetDefineFlag ("store_inner")) SetStoreInner (1);
    geom_free = flags.GetDefineFlag("geom_free");
    matrix_free_bdb = flags.GetDefineFlag("matrix_free_bdb");
    
    precompute = flags.GetDefineFlag ("precompute");
    checksum = flags.GetDefineFlag ("checksum");
    SetCheckUnused (!flags.GetDefineFlagX("check_unused").IsFalse());    
  }


  
  shared_ptr<BilinearFormIntegrator> FixDimension (shared_ptr<BilinearFormIntegrator> bfi, int dim)
  {
    auto anydim = dynamic_pointer_cast<BilinearFormIntegratorAnyDim> (bfi);
    if (anydim) return anydim->GetBFI (dim);

    auto blockbfi = dynamic_pointer_cast<BlockBilinearFormIntegrator> (bfi);
    if (blockbfi)
      {
        auto newblockbfi =
          make_shared<BlockBilinearFormIntegrator> (FixDimension(blockbfi->BlockPtr(), dim), 
                                                    blockbfi->GetDim(), blockbfi->GetComp());
        newblockbfi -> SetDefinedOn (blockbfi->GetDefinedOn());
        return newblockbfi;
      }

    auto compbfi = dynamic_pointer_cast<CompoundBilinearFormIntegrator> (bfi);
    if (compbfi)
      {
        auto newbfi = make_shared<CompoundBilinearFormIntegrator> (FixDimension(compbfi->GetBFI(), dim),
                                                                   compbfi->GetComponent());
        newbfi -> SetDefinedOn (compbfi->GetDefinedOn());
        return newbfi;
      }

    return bfi;
  }
  
  
  BilinearForm & BilinearForm :: AddIntegrator (shared_ptr<BilinearFormIntegrator> bfi)
  {
    /*
    auto anydim = dynamic_pointer_cast<BilinearFormIntegratorAnyDim> (bfi);
    if (anydim) bfi = anydim->GetBFI (ma->GetDimension());
    */
    bfi = FixDimension(bfi, fespace->GetSpatialDimension());
    
    if (symmetric && bfi->IsSymmetric().IsFalse())
      throw Exception (string ("Adding non-symmetric integrator to symmetric bilinear-form\n")+
                       string ("bfi is ")+bfi->Name());

    parts.Append (bfi);

    if ((bfi->geom_free && nonassemble) || geom_free)
      {
        geom_free_parts += bfi;
        return *this;
      }
    
    if(bfi->SkeletonForm())
      {
        auto dgform = bfi -> GetDGFormulation();
        if (dgform.element_boundary) {
	  auto fbfi = dynamic_pointer_cast<FacetBilinearFormIntegrator> (bfi);
	  if (!fbfi)  throw Exception ("not a FacetBFI");
          elementwise_skeleton_parts.Append(fbfi);
          // #ifdef PARALLEL
	  mpi_facet_parts.Append(fbfi);
          // #endif
	}
	else
          {
            if (bfi->VB() > 1) throw Exception ("skeletonform makes sense only for VOL or BND");
            // facetwise_skeleton_parts[bfi->VB()].Append(bfi);
            auto fbfi = dynamic_pointer_cast<FacetBilinearFormIntegrator> (bfi);
            if (!fbfi)  throw Exception ("not a FacetBFI");
            facetwise_skeleton_parts[bfi->VB()] += fbfi;
            // #ifdef PARALLEL
	    if(bfi->VB()==VOL) { //BND-integrators are per definition not mpi!!
	      mpi_facet_parts.Append(fbfi);
	    }
            // #endif
          }
      }
    else
      {
        // the first parts are the symmetric ones ...
        if (bfi->IsSymmetric().IsTrue())
          {
            size_t pos = 0;
            while (pos < VB_parts[bfi->VB()].Size() && VB_parts[bfi->VB()][pos]->IsSymmetric().IsTrue())
              pos++;
            VB_parts[bfi->VB()].Insert(pos, bfi);            
          }
        else
          VB_parts[bfi->VB()].Append(bfi);
      }
    
    if (low_order_bilinear_form)
      low_order_bilinear_form -> AddIntegrator (parts.Last());
    return *this;
  }


  /*
    void BilinearForm :: AddIndependentIntegrator (BilinearFormIntegrator * bfi,
    const int master_surface,
    const int other,
    const bool deletable)
    {
    independent_parts.Append (bfi);
    independent_parts_deletable.Append(deletable);
    Vec<2,int> indices(master_surface,other); 
    independent_meshindex.Append(indices);
    if (low_order_bilinear_form)
    low_order_bilinear_form -> AddIndependentIntegrator (independent_parts.Last(),
    master_surface,other,
    deletable);
    }
  */

  BilinearForm :: ~BilinearForm ()
  { ; }

  void BilinearForm :: SetPrint (bool ap)
  { 
    print = ap; 
    if (low_order_bilinear_form)
      low_order_bilinear_form -> SetPrint (ap);
  }

  void BilinearForm :: SetPrintElmat (bool ap)
  { 
    printelmat = ap; 
    if (low_order_bilinear_form)
      low_order_bilinear_form -> SetPrintElmat (ap);
  }

  void BilinearForm :: SetElmatEigenValues (bool ee)
  { 
    elmat_ev = ee;
    if (low_order_bilinear_form)
      low_order_bilinear_form -> SetElmatEigenValues (ee);
  }

  void BilinearForm :: SetCheckUnused (bool b)
  {
    check_unused = b;
    if (low_order_bilinear_form)
      low_order_bilinear_form -> SetCheckUnused (b);
  }

  void BilinearForm :: AddSpecialElement (unique_ptr<SpecialElement> spel)
  {
    specialelements.Append (std::move(spel));
    special_element_coloring = nullptr;
    specialelements_timestamp = GetNextTimeStamp();
  }

  void BilinearForm :: DeleteSpecialElement(size_t index)
  {
    specialelements.DeleteElement(index);
    special_element_coloring = nullptr;    
    specialelements_timestamp = GetNextTimeStamp();
  }

  void BilinearForm :: DeleteSpecialElements()
  {
    // for(auto el : specialelements)
    // delete el;
    specialelements.DeleteAll();
    special_element_coloring = nullptr;    
    specialelements_timestamp = GetNextTimeStamp();
  }

  Table<int> & BilinearForm :: SpecialElementColoring() const
  {
    if (!special_element_coloring)
      {
        cout << "building special element coloring" << endl;
        NETGEN_TIMER_FROM_HERE("SpecialElementColoring");

        size_t ndof = GetTestSpace()->GetNDof();
        Array<MyMutex> locks(ndof);

        Array<int> col(specialelements.Size());
        col = -1;

        int maxcolor = 0;
        
        int basecol = 0;
        Array<unsigned int> mask(ndof);

        atomic<int> found(0);
        size_t cnt = specialelements.Size();

        while (found < cnt)
          {
            ParallelForRange
              (mask.Size(),
               [&] (IntRange myrange) { mask[myrange] = 0; });
            
            // size_t ne = ma->GetNE(vb);
            
            ParallelForRange
              (specialelements.Size(), [&] (IntRange myrange)
               {
                 Array<DofId> dofs;
                 size_t myfound = 0;
                 
                 for (size_t nr : myrange)
                   {
                     // ElementId el = { vb, nr };
                     // if (!DefinedOn(el)) continue;
                     if (col[nr] >= 0) continue;
                     
                     unsigned check = 0;
                     specialelements[nr]->GetDofNrs(dofs);
                     
                     for (int i = dofs.Size()-1; i >= 0; i--)
                       if (!IsRegularDof(dofs[i])) dofs.DeleteElement(i);
                     QuickSort (dofs);   // sort to avoid dead-locks
                     
                     for (auto d : dofs) 
                       locks[d].lock();
                     
                     for (auto d : dofs) 
                       check |= mask[d];
                     
                     if (check != UINT_MAX) // 0xFFFFFFFF)
                       {
                         myfound++;
                         unsigned checkbit = 1;
                         int color = basecol;
                         while (check & checkbit)
                           {
                             color++;
                             checkbit *= 2;
                           }
                         
                         col[nr] = color;
                         if (color > maxcolor) maxcolor = color;
                         
                         for (auto d : dofs) 
                           mask[d] |= checkbit;
                       }
                     
                     for (auto d : dofs) 
                       locks[d].unlock();
                   }
                 found += myfound;
               });
                 
            basecol += 8*sizeof(unsigned int); // 32;
          }

        Array<int> cntcol(maxcolor+1);
        cntcol = 0;

        // for (ElementId el : Elements(vb))
        for (int nr : Range(specialelements))
          cntcol[col[nr]]++;
        
        special_element_coloring = make_unique<Table<int>> (cntcol);
        Table<int> & coloring = *special_element_coloring;

	cntcol = 0;
        for (int nr : Range(specialelements))          
          coloring[col[nr]][cntcol[col[nr]]++] = nr;
        
        cout << "needed " << maxcolor+1 << " colors for special elements" << endl;
        cout << "coloring = " << coloring << endl;
      }        

    return *special_element_coloring;
  }

  
  void BilinearForm :: SetPreconditioner (Preconditioner * pre)
  {
    // cout << "SetPreconditioner, type fes = " << typeid(*fespace).name() << ", type pre = " << typeid(*pre).name() << endl;
    if (preconditioners.Contains(pre))
      throw Exception (string("preconditioner ")+typeid(*pre).name()+ " already registered in bfa");
    if (pre->GetFlags().GetDefineFlag("not_register_for_auto_update"))
      throw Exception (string("'not_register_for_auto_update' set, but preconditioner") + typeid(*pre).name() +" registers anyway");
    preconditioners.Append (pre);
  }

  void BilinearForm :: UnsetPreconditioner (Preconditioner* pre)
  {
    if(auto pos = preconditioners.Pos(pre); pos != preconditioners.ILLEGAL_POSITION)
      preconditioners.DeleteElement(pos);
  }


  MatrixGraph BilinearForm :: GetGraph (int level, bool symmetric)
  {
    static Timer timer ("BilinearForm::GetGraph");
    RegionTimer reg (timer);

    size_t ndof = fespace->GetNDof();
    size_t nf = ma->GetNFacets();
    size_t neV = ma->GetNE(VOL);
    size_t neB = ma->GetNE(BND);
    size_t neBB = ma->GetNE(BBND);
    // const Array<SpecialElement*> & specialelements = fespace->GetSpecialElements();
    size_t nspe = specialelements.Size();

    Array<DofId> dnums;
    Array<int> fnums; //facets of one element
    Array<int> elnums; //elements neighbouring one facet
    Array<int> nbelems; //neighbour elements


    int maxind = neV + neB + neBB + specialelements.Size();
    if (fespace->UsesDGCoupling()) maxind += nf;

    TableCreator<int> creator(maxind);
    for ( ; !creator.Done(); creator++)
      {

        /*
        task_manager->CreateJob
          ([&] (const TaskInfo & ti)

           {
             Array<int> dnums;

             auto myr = Range(ne).Split(ti.task_nr, ti.ntasks);
             for (auto i : myr)
               {
                 if (!fespace->DefinedOn (ma->GetElIndex(i))) continue;
                 
                 if (eliminate_internal)
                   fespace->GetDofNrs (i, dnums, EXTERNAL_DOF);
                 else
                   fespace->GetDofNrs (i, dnums);
                 
                 for (int d : dnums)
                   if (d != -1) creator.Add (i, d);
               }

             
             auto mysr = Range(nse).Split(ti.task_nr, ti.ntasks);
             for (auto i : mysr)
               {
                 if (!fespace->DefinedOnBoundary (ma->GetSElIndex(i))) continue;
                 
                 fespace->GetSDofNrs (i, dnums);
                 for (int d : dnums)
                   if (d != -1) creator.Add (ne+i, d);
               }
           });
        */
	for(VorB vb : {VOL, BND, BBND})
	  {
            size_t shift = (vb==VOL) ? 0 : ((vb==BND) ? neV : neV+neB);
            bool condensation_allowed = (vb == VOL) || ((neV==0) && (vb == BND));
	    ParallelForRange
              (ma->GetNE(vb), [&](IntRange r)
               {
                 Array<DofId> dnums;
                 for (auto i : r)
                   {
                     auto eid = ElementId(vb,i);
                     if (!fespace->DefinedOn (vb, ma->GetElIndex(eid))) continue;
                     
                     if (condensation_allowed && eliminate_internal)
                       fespace->GetDofNrs (eid, dnums, EXTERNAL_DOF);
                     else if (condensation_allowed && eliminate_hidden)
                       fespace->GetDofNrs (eid, dnums, VISIBLE_DOF);
                     else
                       fespace->GetDofNrs (eid, dnums);
                         
                     
                     for (DofId d : dnums)
                       if (IsRegularDof(d)) creator.Add (shift+i, d);
                   }
               });
	  }

        
        for (int i = 0; i < specialelements.Size(); i++)
          {
            specialelements[i]->GetDofNrs (dnums);
            QuickSort(dnums);
            int last = -1;
            for (int d : dnums)
              {
                if (d!=last && IsRegularDof(d))
                    creator.Add (neV+neB+neBB+i, d);
                last = d;
              }
          }

        if (fespace->UsesDGCoupling())
        {
          //add dofs of neighbour elements as well
          Array<DofId> dnums_dg;
          Array<int> elnums_per;

          for (int i = 0; i < nf; i++)
            {
              nbelems.SetSize(0);
              ma->GetFacetElements(i,elnums);

              // timerDG1.Stop();
              if(elnums.Size() < 2)
              {
                int facet2 = ma->GetPeriodicFacet(i);
                if(facet2 > i)
                {
                  ma->GetFacetElements (facet2, elnums_per);
		  // if the facet is identified across subdomain
		  // boundary, we only have the surface element
		  // and not the other volume element!
		  if (elnums_per.Size())
		    elnums.Append(elnums_per[0]);
                }
              }
              
              for (int k=0; k<elnums.Size(); k++)
                nbelems.Append(elnums[k]);

              dnums_dg.SetSize(0);
              for (int k=0;k<nbelems.Size();k++){
                int elnr=nbelems[k];
                if (!fespace->DefinedOn (VOL,ma->GetElIndex(ElementId(VOL,elnr)))) continue;
                fespace->GetDofNrs (ElementId(VOL,elnr), dnums);
                dnums_dg.Append(dnums);
              }
              QuickSort (dnums_dg);
              for (int j = 0; j < dnums_dg.Size(); j++)
                if (IsRegularDof(dnums_dg[j]) && (j==0 || (dnums_dg[j] != dnums_dg[j-1]) ))
                  creator.Add (neV+neB+neBB+nspe+i, dnums_dg[j]);
            }
        }

      }
    
    MatrixGraph * graph;

    if (!fespace2)
      {
        /*
        graph = new MatrixGraph (ndof, *creator.GetTable(), *creator.GetTable(), symmetric);
        delete creator.GetTable();
        */
        auto table = creator.MoveTable();
        graph = new MatrixGraph (ndof, ndof, table, table, symmetric);        
      }
    else
      {
        TableCreator<int> creator2(maxind);
        for ( ; !creator2.Done(); creator2++)
          {
            int offset=0;
	    for(VorB vb  : {VOL, BND, BBND})
	      {
		int nre = ma->GetNE(vb);
                offset += nre;
		for (int i = 0; i < nre; i++)
		  {
		    auto eid = ElementId(vb,i);
		    if (!fespace2->DefinedOn (vb,ma->GetElIndex(eid))) continue;
		    
		    if (vb == VOL && eliminate_internal)
		      fespace2->GetDofNrs (eid, dnums, EXTERNAL_DOF);
		    else if (vb == VOL && eliminate_hidden)
		      fespace2->GetDofNrs (eid, dnums, VISIBLE_DOF);
		    else
		      fespace2->GetDofNrs (eid, dnums);
		    
		    int shift = (vb==VOL) ? 0 : ((vb==BND) ? neV : neV + neB);
		    for (int d : dnums)
		      if (IsRegularDof(d)) creator2.Add (shift+i, d);
		  }
              }
		
            // just not tested ...
            for (int i = 0; i < specialelements.Size(); i++)
            {
              specialelements[i]->GetDofNrs (dnums);
		
              for (int j = 0; j < dnums.Size(); j++)
                if (dnums[j] != -1)
                  creator2.Add (offset+i, dnums[j]);
            }
            offset += specialelements.Size();	
            if (fespace2->UsesDGCoupling())
              //add dofs of neighbour elements as well
              for (int i = 0; i < nf; i++)
              {
                nbelems.SetSize(0);
                ma->GetFacetElements(i,elnums);
                for (int k=0; k<elnums.Size(); k++)
                  nbelems.Append(elnums[k]);
		
                for (int k=0;k<nbelems.Size();k++){
                  int elnr=nbelems[k];
                  if (!fespace2->DefinedOn (ElementId(VOL,elnr))) continue;
                  fespace2->GetDofNrs (ElementId(VOL,elnr), dnums);
                  for (int j = 0; j < dnums.Size(); j++)
                    if (dnums[j] != -1)
                      creator2.Add (offset+i, dnums[j]);
                }
              }
	  }

        /*
        graph = new MatrixGraph (fespace2->GetNDof(), *creator2.GetTable(), *creator.GetTable(), symmetric);        
        delete creator.GetTable();
        delete creator2.GetTable();
        */
        auto table = creator.MoveTable();
        auto table2 = creator2.MoveTable();
        graph = new MatrixGraph (fespace2->GetNDof(), fespace->GetNDof(),
                                 table2, table, symmetric);
      }
    
    graph -> FindSameNZE();
    return move(*graph);
  }






  void BilinearForm :: Assemble (LocalHeap & lh)
  {
    if (mats.Size() == ma->GetNLevels())
      return;


    if (nonassemble)
      {
        // mats.Append (make_shared<BilinearFormApplication> (shared_ptr<BilinearForm>(this, NOOP_Deleter)), lh);

        mats.SetSize(ma->GetNLevels());
        shared_ptr<BaseMatrix> app = 
          make_shared<BilinearFormApplication> (dynamic_pointer_cast<BilinearForm>(this->shared_from_this()), lh);
        cout << "craete bilinearformapplication" << endl;
        if (fespace->IsParallel())
          app = make_shared<ParallelMatrix> (app, GetTrialSpace()->GetParallelDofs(), GetTestSpace()->GetParallelDofs(), C2D);
        
        mats.Last() = app;
      
        if (precompute)
          {
            Array<DofId> dnums;
            precomputed_data.SetSize( max2(max2(ma->GetNE(VOL), ma->GetNE(BND)),ma->GetNE(BBND))*NumIntegrators() );
            precomputed_data = nullptr;

            LocalHeap lh (20000000, "biform - assemble");
            
            // int hasbound = 0;
            // int hasinner = 0;
            for (VorB vb : {VOL, BND, BBND})
	      {
		int ne = ma->GetNE(vb);
		if(VB_parts[vb].Size())
		  {
		    for (int i=0; i < ne; i++)
		      {
                        ElementId ei(vb, i);
			// auto & fel = fespace->GetFE(ei,lh);
			// auto & eltrans = ma->GetTrafo(ei,lh);
			fespace->GetDofNrs(ei,dnums);
                        /*
			for(int j=0; j<NumIntegrators(); j++)
			  {
			    const BilinearFormIntegrator & bfi = *parts[j];
                        */
                        /*
                        for (auto & bfi : VB_parts[vb])
                          {
			    // if(bfi.VB()!=vb) continue;
			    if(!bfi->DefinedOn (ma->GetElIndex(i))) continue;
                            
			    precomputed_data[i*NumIntegrators()+j] =
			      bfi->PrecomputeData(fel,eltrans,lh);
			  }
                        */
                        throw Exception ("precomputed currently not supported");
		      }
		  }
	      }
	  }
            
        
        if (timing)
          {
            Timer timer("bftimer");

            auto vecf = mats.Last()->CreateColVector();
            auto vecu = mats.Last()->CreateColVector();
          
            *vecu = 1;
            do
              {
                timer.Start();
                vecf = (*mats.Last()) * vecu;
                timer.Stop();
              }
            while (timer.GetTime() < 2.0);
          
            cout << " 1 application takes " << timer.GetTime() / timer.GetCounts()
                 << " seconds" << endl;
          }


        return;
      }
    

    if (geom_free)
      {
        AssembleGF(lh);
        return;
      }

    if (matrix_free_bdb || nonlinear_matrix_free_bdb)
      {
        AssembleBDB(lh, matrix_free_bdb);
        return;
      }
    
    try
      {

        if (low_order_bilinear_form)
          low_order_bilinear_form->Assemble(lh);

      }
    catch (Exception & e)
      {
        e.Append (string ("\nthrown by Do loworder ") +
                  string (GetName()));
        throw e;
      }

    
    try
      {
        AllocateMatrix ();
        graph_timestamp = GetNextTimeStamp();
      }
    catch (Exception & e)
      {
        e.Append (string ("\nthrown by allocate matrix ") +
                  string (GetName()));
        throw;
      }
    catch (exception & e)
      {
        throw Exception (e.what() + 
                         string ("\nthrown by allocate matrix ") +
                         string (GetName()));
      }


    DoAssemble(lh);


    if (timing)
      {
        double time;
        double starttime = WallTime();
      
        auto vecf = mats.Last()->CreateColVector();
        auto vecu = mats.Last()->CreateRowVector();
      
        vecu = 1;
        int steps = 0;
        do
          {
            // vecf = (*mats.Last()) * vecu;
            vecf = Transpose (*mats.Last()) * vecu;
            steps++;
            time = WallTime()-starttime;
          }
        while (time < 2.0);
          
        cout << " 1 application takes "
             << time / steps
             << " seconds" << endl;

        auto & lastmat = *mats.Last();
        int nze = lastmat.NZE();
        cout << "NZE = " << nze << ", MFLOP = " << double (nze * steps) / time * 1e-6 << endl;
        cout << "type = " << typeid(lastmat).name() << endl;
      }

    if (galerkin)
      GalerkinProjection();
  }


  void BilinearForm :: AssembleGF (LocalHeap & lh)
  {
    auto fesx = GetTrialSpace();
    auto fesy = GetTestSpace();
    auto ma = GetMeshAccess();


    shared_ptr<BaseMatrix> sum;
    
    for (VorB vb : {VOL, BND, BBND})
      {
        bool has_vb = false;
        for (auto bfi : geom_free_parts)
          if (bfi->VB() == vb)
            has_vb = true;
        if (!has_vb)
          continue;

        Array<short> classnr(ma->GetNE(vb));
        ma->IterateElements
          (vb, lh, [&] (auto el, LocalHeap & llh)
           {
             classnr[el.Nr()] = 
               SwitchET<ET_SEGM, ET_TRIG,ET_TET>
               (el.GetType(),
                [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
           });
        
        TableCreator<size_t> creator;
        for ( ; !creator.Done(); creator++)
          for (auto i : Range(classnr))
            creator.Add (classnr[i], i);
        Table<size_t> table = creator.MoveTable();
        
    
        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
            
            // size_t nr = classnr[elclass_inds[0]];
            ElementId ei(vb,elclass_inds[0]);
            auto & felx = GetTrialSpace()->GetFE (ei, lh);
            auto & fely = GetTestSpace()->GetFE (ei, lh);
            auto & trafo = GetTrialSpace()->GetMeshAccess()->GetTrafo(ei, lh);
            MixedFiniteElement fel(felx, fely);
        
            Matrix<> elmat(fely.GetNDof(), felx.GetNDof());

            
            bool done = false;
            while (!done)
              {
                try
                  {
                    elmat = 0.0;
                    bool symmetric_so_far = true;
                    for (auto bfi : geom_free_parts)
                      if (bfi->VB() == vb)
                        bfi->CalcElementMatrixAdd(fel, trafo, elmat, symmetric_so_far, lh);
                    done = true;
                  }
                catch (const ExceptionNOSIMD& e)
                  {
                    ;
                  }
              }

            /*
              elmat = 0.0;
              bool symmetric_so_far = true;
              for (auto bfi : geom_free_parts)
              if (bfi->VB() == vb)
              bfi->CalcElementMatrixAdd(fel, trafo, elmat, symmetric_so_far, lh);
            */
            
            Table<DofId> xdofs(elclass_inds.Size(), felx.GetNDof()),
              ydofs(elclass_inds.Size(), fely.GetNDof());
            
            Array<DofId> dnumsx, dnumsy;
            for (auto i : Range(elclass_inds))
              {
                ElementId ei(vb, elclass_inds[i]);
                fesx->GetDofNrs(ei, dnumsx);
                fesy->GetDofNrs(ei, dnumsy);
                xdofs[i] = dnumsx;
                ydofs[i] = dnumsy;
              }
            
            // cout << "elmat = " << elmat << endl;
            // cout << "xdofs = " << xdofs << endl;
            // cout << "ydofs = " << ydofs << endl;
            auto mat = make_shared<ConstantElementByElementMatrix>
              (fesy->GetNDof(), fesx->GetNDof(),
               elmat, std::move(ydofs), std::move(xdofs));
            
            if (sum)
              sum = make_shared<SumMatrix>(sum, mat);
            else
              sum = mat;
          }
      }
    
        // mats.Append(sum);
    mats.SetSize (ma->GetNLevels());
    mats.Last() = sum;
  }



  /*
    Stores product of B2 D B1 matrices

    B matrices ... sum of EBE const 
    D bock-diagonal
   */

  class ApplyIntegrationPoints : public BaseMatrix
  {
    Array<shared_ptr<CoefficientFunction>> coefs;
    Array<ProxyFunction*> trialproxies;
    
    size_t dimx, dimy;
    size_t nip;
    
  public:
    ApplyIntegrationPoints (Array<shared_ptr<CoefficientFunction>> acoefs,
                            const Array<ProxyFunction*> & atrialproxies,
                            size_t adimx, size_t adimy, size_t anip)
      : coefs(acoefs), trialproxies{atrialproxies}, dimx(adimx), dimy(adimy), nip(anip)

    { 
      // make my own code
      for (auto cf : coefs)
        {
          auto compiledcf = Compile (cf, false);
          Code code = compiledcf->GenerateProgram(0, false);
          /*
          cout << "Generated code:" << endl;
          cout << code.header << endl;
          cout << code.body << endl;
          */
        }
    }
      
    AutoVector CreateColVector() const override
    { return make_unique<VVector<double>> (nip*dimy); }
    AutoVector CreateRowVector() const override
    { return make_unique<VVector<double>> (nip*dimx); }

    virtual int VHeight() const override { return nip*dimy; }
    virtual int VWidth() const override { return nip*dimx; }
    
    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      static Timer t("ApplyIntegrationPoints"); RegionTimer reg(t);
      static Timer teval("ApplyIntegrationPoints eval");
      static Timer tmir("ApplyIntegrationPoints mir");      
      static Timer ttransx("ApplyIntegrationPoints transx");
      static Timer ttransy("ApplyIntegrationPoints transy");

      ParallelForRange
        (nip, [&](IntRange r)
         {
           constexpr size_t BS = 256;
           LocalHeap lh(1000000);

           for (size_t ii = r.First(); ii < r.Next(); ii+=BS)
             {
               IntRange r2(ii, min(ii+BS, r.Next()));               
               HeapReset hr(lh);

               SIMD_IntegrationRule simdir(r2.Size(), lh);               
               FE_ElementTransformation<2,2> trafo(ET_TRIG);
               SIMD_MappedIntegrationRule<2,2> simdmir(simdir, trafo, 0, lh); // don't actually compute

               // tmir.Stop();
           
               ProxyUserData ud(trialproxies.Size(), 0, lh);
               trafo.userdata = &ud;
               ScalarFE<ET_TRIG,1> dummyfe;
               ud.fel = &dummyfe;

               int starti = 0;
               for (auto proxy : trialproxies)
                 {
                   int nexti = starti + proxy->Evaluator()->Dim();
                   ud.AssignMemory (proxy, r2.Size(), proxy->Evaluator()->Dim(), lh);

                   // ttransx.Start();
                   SliceMatrix<double> (dimx, r2.Size(), SIMD<double>::Size()*simdmir.Size(),
                                        (double*)(ud.GetAMemory (proxy)).Data()) = 
                     x.FV<double>().AsMatrix(dimx, nip).Cols(r2).Rows(starti, nexti);
                   // ttransx.Stop();
                   
                   starti = nexti;
                 }
               
           
               // teval.Start();
               FlatMatrix<SIMD<double>> simdres(dimy, simdmir.Size(), lh);
               starti = 0;
               for (auto coef : coefs)
                 {
                   int nexti = starti + coef->Dimension();
                   coef -> Evaluate (simdmir, simdres.Rows(starti, nexti));
                   starti = nexti;
                 }
               // teval.Stop();
               
               // ttransy.Start();
               SliceMatrix<> res = y.FV<double>().AsMatrix(dimy, nip).Cols(r2);
               res = SliceMatrix<double> (dimy, r2.Size(), SIMD<double>::Size()*simdmir.Size(),
                                          (double*)simdres.Data());
               // ttransy.Stop();
             }
         }, TasksPerThread(3));
      
    }
  };

  
  void BilinearForm :: AssembleBDB (LocalHeap & lh, bool linear)
  {
    static Timer t("assemble-BDB"); RegionTimer reg(t);
    
    auto fesx = GetTrialSpace();
    auto fesy = GetTestSpace();
    auto ma = GetMeshAccess();


    Array<short> classnr(ma->GetNE(VOL));
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_SEGM, ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
       });
        
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();
    

    shared_ptr<BaseMatrix> sum;

    for (auto part : parts)
      {
        auto bfi = dynamic_pointer_cast<SymbolicBilinearFormIntegrator> (part);

        auto & trialproxies = bfi->TrialProxies();
        auto & testproxies = bfi->TestProxies();
        
        int dimx = 0, dimy = 0;
        for (auto proxy : trialproxies)
          dimx += proxy->Evaluator()->Dim();
        for (auto proxy : testproxies)
          dimy += proxy->Evaluator()->Dim();
        
        int dimxref = 0, dimyref = 0;
        for (auto proxy : trialproxies)
          dimxref += proxy->Evaluator()->DimRef();
        for (auto proxy : testproxies)
          dimyref += proxy->Evaluator()->DimRef();
        


        for (auto elclass_inds : table)
          {
            if (elclass_inds.Size() == 0) continue;
        
            ElementId ei(VOL,elclass_inds[0]);
            auto & felx = GetTrialSpace()->GetFE (ei, lh);
            auto & fely = GetTestSpace()->GetFE (ei, lh);
        
            MixedFiniteElement fel(felx, fely);
            // const IntegrationRule & ir = bfi->GetIntegrationRule(felx.ElementType(), felx.Order()+fely.Order());
            IntegrationRule ir;
            if (bfi->ElementVB() == VOL)
              {
                const IntegrationRule & volir = bfi->GetIntegrationRule(felx.ElementType(), felx.Order()+fely.Order());
                for (auto ip : volir)
                  ir += ip;
              }
            else
              {
                auto eltype = felx.ElementType();
                
                Facet2ElementTrafo transform(eltype, bfi->ElementVB()); 
                int nfacet = transform.GetNFacets();
                
                for (int k = 0; k < nfacet; k++)
                  {
                    HeapReset hr(lh);
                    ngfem::ELEMENT_TYPE etfacet = transform.FacetType (k);
                    IntegrationRule ir_facet(etfacet, felx.Order()+fely.Order() /* +bonus_intorder */);
                    IntegrationRule & ir_facet_vol = transform(k, ir_facet, lh);
                    for (auto ip : ir_facet_vol)
                      ir += ip;
                  }
              }
            
            Matrix<double,ColMajor> bmatx_(ir.Size()*dimxref, felx.GetNDof());
            Matrix<double,ColMajor> bmaty_(ir.Size()*dimyref, fely.GetNDof());
        
            for (int i : Range(ir.Size()))
              {
                int starti = i*dimxref;
                for (auto proxy : trialproxies)
                  {
                    auto diffopx = proxy->Evaluator();
                    int nexti = starti+diffopx->DimRef();
                    diffopx->CalcMatrix(felx, ir[i], bmatx_.Rows(starti, nexti), lh);
                    starti = nexti;
                  }
                
                starti = i*dimyref;
                for (auto proxy : testproxies)
                  {
                    auto diffopy = proxy->Evaluator();
                    int nexti = starti+diffopy->DimRef();
                    diffopy->CalcMatrix(fely, ir[i], bmaty_.Rows(starti, nexti), lh);
                    starti = nexti;
                  }
              }

            Matrix bmatx = bmatx_;
            Matrix bmaty = bmaty_;
            
            Table<DofId> xdofsin(elclass_inds.Size(), felx.GetNDof());
            Table<DofId> xdofsout(elclass_inds.Size(), bmatx.Height());

            Table<DofId> ydofsin(elclass_inds.Size(), fely.GetNDof());
            Table<DofId> ydofsout(elclass_inds.Size(), bmaty.Height());

            Array<DofId> dnumsx, dnumsy;
            for (auto i : Range(elclass_inds))
              {
                ElementId ei(VOL, elclass_inds[i]);
                fesx->GetDofNrs(ei, dnumsx);
                fesy->GetDofNrs(ei, dnumsy);
                xdofsin[i] = dnumsx;
                ydofsin[i] = dnumsy;
              }

            int nel = elclass_inds.Size();
            size_t nip = ir.Size()*nel;
            
            auto xa = xdofsout.AsArray();
            auto ya = ydofsout.AsArray();
            
            if (linear)
              {
                for (size_t i = 0; i < xa.Size(); i++)
                  xa[i] = i;
                for (size_t i = 0; i < ya.Size(); i++)
                  ya[i] = i;
              }
            else
              {
                for (size_t i = 0; i < nip; i++)
                  for (size_t j = 0; j < dimxref; j++)
                    xa[i*dimxref+j] = j*nip+i;
                for (size_t i = 0; i < nip; i++)
                  for (size_t j = 0; j < dimyref; j++)
                    ya[i*dimyref+j] = j*nip+i;
              }
            
            auto bx = make_shared<ConstantElementByElementMatrix>
              (nip*dimxref, fesx->GetNDof(),
               bmatx, std::move(xdofsout), std::move(xdofsin));
            
            auto by = make_shared<ConstantElementByElementMatrix>
              (nip*dimyref, fesy->GetNDof(),
               bmaty, std::move(ydofsout), std::move(ydofsin));


            shared_ptr<BaseMatrix> mat;
            
            if (linear)
              {
                Tensor<3> diag(elclass_inds.Size()*ir.Size(), dimyref, dimxref);
                for (auto i : Range(elclass_inds))
                  {
                    HeapReset hr(lh);
                    ElementId ei(VOL, elclass_inds[i]);
                    auto & trafo = ma->GetTrafo(ei, lh);
                    auto & mir = trafo(ir, lh);
                    if (bfi->ElementVB() != VOL) 
                      mir.ComputeNormalsAndMeasure (fel.ElementType());
                    FlatMatrix transx(dimx, dimxref, lh);
                    FlatMatrix transy(dimy, dimyref, lh);
                    Matrix prod(dimxref, dimyref);
                    
                    for (int j = 0; j < ir.Size(); j++)
                      {
                        auto diffopx = trialproxies[0]->Evaluator();
                        auto diffopy = testproxies[0]->Evaluator();
                        
                        diffopx->CalcTransformationMatrix(mir[j], transx, lh);
                        diffopy->CalcTransformationMatrix(mir[j], transy, lh);
                        
                        prod = Trans(transy) * transx;
                        prod *= mir[j].GetWeight();
                        diag(i*ir.Size()+j,STAR,STAR) = prod;
                      }
                  }
                
                auto diagmat = make_shared<BlockDiagonalMatrix> (std::move(diag));
                mat = TransposeOperator(by) * diagmat * bx;
              }
            else
              {
                Tensor<3> diagx(dimx, dimxref, nip);
                Tensor<3> diagy(dimy, dimyref, nip);
                
                for (auto i : Range(elclass_inds))
                  {
                    HeapReset hr(lh);
                    ElementId ei(VOL, elclass_inds[i]);
                    auto & trafo = ma->GetTrafo(ei, lh);
                    auto & mir = trafo(ir, lh);
                    if (bfi->ElementVB() != VOL) 
                      mir.ComputeNormalsAndMeasure (fel.ElementType());
                    
                    FlatMatrix transx(dimx, dimxref, lh);
                    FlatMatrix transy(dimy, dimyref, lh);
                    
                    transx = 0.0;
                    transy = 0.0;
                    for (int j = 0; j < ir.Size(); j++)
                      {
                        int starti = 0, startiref = 0;
                        for (auto proxy : trialproxies)
                          {
                            auto diffop = proxy->Evaluator();
                            int nexti = starti+diffop->Dim();
                            int nextiref = startiref+diffop->DimRef();
                            diffop->CalcTransformationMatrix(mir[j], transx.Rows(starti,nexti).Cols(startiref,nextiref), lh);
                            starti = nexti;
                            startiref = nextiref;
                          }
                        starti = 0; startiref = 0;
                        for (auto proxy : testproxies)
                          {
                            auto diffop = proxy->Evaluator();
                            int nexti = starti+diffop->Dim();
                            int nextiref = startiref+diffop->DimRef();                        
                            diffop->CalcTransformationMatrix(mir[j], transy.Rows(starti,nexti).Cols(startiref,nextiref), lh);
                            starti = nexti;
                            startiref = nextiref;
                          }

                        transy *= mir[j].GetWeight();
                        diagx(STAR,STAR,i*ir.Size()+j) = transx;
                        diagy(STAR,STAR,i*ir.Size()+j) = transy;
                      }
                  }
                
                shared_ptr<CoefficientFunction> coef = bfi -> GetCoefficientFunction();
                Array<shared_ptr<CoefficientFunction>> diffcfs;
                for (auto proxy : testproxies)
                  {
                    CoefficientFunction::T_DJC cache;
                    diffcfs += coef -> DiffJacobi(proxy, cache);                  
                  }
                
                auto ipop = make_shared<ApplyIntegrationPoints> (std::move(diffcfs), trialproxies, dimx, dimy, nip);
                
                auto diagmatx = make_shared<BlockDiagonalMatrixSoA> (std::move(diagx));
                auto diagmaty = make_shared<BlockDiagonalMatrixSoA> (std::move(diagy));
                
                mat = TransposeOperator(diagmaty * by) * ipop * (diagmatx * bx);
              }
            
            if (sum)
              sum = sum + mat;
            else
              sum = mat;
          }
      }

    mats.SetSize (ma->GetNLevels());
    mats.Last() = sum;
  }


  

  void BilinearForm :: ReAssemble (LocalHeap & lh, bool reallocate)
  {
    if (nonassemble)
      {
        Assemble(lh);
        return;
      }

    if (low_order_bilinear_form)
      low_order_bilinear_form->ReAssemble(lh);

    if (mats.Size() < ma->GetNLevels())
      {
        Assemble(lh);
        return;
      }

    if (specialelements_timestamp > graph_timestamp)
      {
        reallocate = true;
        cout << IM(3) << "reallocate due to changed special elements" << endl;
      }
    
    if (reallocate)
      {
        // delete mats.Last();
        mats.DeleteLast();

        Assemble(lh);
        return;
      }

    GetMatrix() = 0.0;
    DoAssemble(lh);

    if (galerkin)
      GalerkinProjection();
  }

  shared_ptr<BaseMatrix> BilinearForm :: GetMatrixPtr () const
  {
    if (!mats.Size())
      return nullptr;
    return mats.Last(); 
  }
  

  void BilinearForm :: PrintReport (ostream & ost) const
  {
    ost << "on space " << GetFESpace()->GetName() << endl
        << "symmetric   = " << symmetric << endl
        << "multilevel  = " << multilevel << endl
        << "nonassemble = " << nonassemble << endl
        << "printelmat = " << printelmat << endl
        << "elmatev    = " << elmat_ev << endl
        << "eliminate_internal = " << eliminate_internal << endl
        << "eliminate_hidden = " << eliminate_hidden << endl
        << "keep_internal = " << keep_internal << endl
        << "store_inner = " << store_inner << endl
        << "integrators: " << endl;
  
    for (int i = 0; i < parts.Size(); i++)
      ost << "  " << parts[i]->Name() << endl;
  }


  Array<MemoryUsage> BilinearForm :: GetMemoryUsage () const
  {
    Array<MemoryUsage> mu;
    if (low_order_bilinear_form)
      mu = low_order_bilinear_form -> GetMemoryUsage ();
  
    int olds = mu.Size();

    for (int i = 0; i < mats.Size(); i++)
      if (mats[i]) mu += mats[i]->GetMemoryUsage ();

    for (int i = olds; i < mu.Size(); i++)
      mu[i].AddName (string(" bf ")+GetName());
    return mu;
  }


  template <class SCAL>
  S_BilinearForm<SCAL> :: ~S_BilinearForm () { ; }

  
  template <class SCAL>
  void S_BilinearForm<SCAL> :: AllocateInternalMatrices ()
  {
    if (eliminate_internal && keep_internal)
      {
        VorB vb = VOL;
        if (!VB_parts[VOL].Size()) vb = BND;
        
        size_t ne = ma->GetNE(vb);
        size_t ndof = fespace->GetNDof();
        size_t dim = fespace->GetDimension();
        // ndof *= dim;

        
        Array<int> nidofs(ne), nodofs(ne);
        nidofs = 0;
        nodofs = 0;
        
        ParallelForRange
          (ne, [&] (IntRange r)
           {
             Array<DofId> dnums;
             for (size_t i : r)
               {
                 int ni = 0, no = 0;
                 ElementId ei(vb,i);
                 
                 if (fespace->DefinedOn(ei))
                   {
                     fespace->GetDofNrs (ei, dnums );
                     for (auto d : dnums)
                       {
                         auto ct = fespace->GetDofCouplingType(d);
                         if (ct & EXTERNAL_DOF)
                           no++;
                         if (ct & LOCAL_DOF)
                           ni++;
                       }
                   }
                 nidofs[i] = dim*ni;
                 nodofs[i] = dim*no;
               }
           });

        /*
        harmonicext = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
        if (!symmetric)
          harmonicexttrans = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
        else
          harmonicexttrans = make_shared<Transpose>(*harmonicext);
        innersolve = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
        if (store_inner)
          innermatrix = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
        */
        auto ext = make_shared<ElementByElementMatrix<SCAL>>(ndof, ndof, nidofs, nodofs, false, false, false);
        harmonicext = ext;
        harmonicext_ptr = ext.get();
        
        if (!symmetric)
          {
            auto extt = make_shared<ElementByElementMatrix<SCAL>>(ndof, ndof, nodofs, nidofs, false, false, false);
            harmonicexttrans = extt;
            harmonicexttrans_ptr = extt.get();
          }
        else
          {
            harmonicexttrans = make_shared<Transpose>(*harmonicext);
            harmonicexttrans_ptr = nullptr;
          }
        auto isol = make_shared<ElementByElementMatrix<SCAL>>(ndof, ndof, nidofs, nidofs, false, true, true);
        innersolve = isol; 
        innersolve_ptr = isol.get();
        
        if (store_inner)
          {
            auto imat = make_shared<ElementByElementMatrix<SCAL>>(ndof, ndof, nidofs, nidofs, false, true, true);
            innermatrix = imat;
            innermatrix_ptr = imat.get();
          }
        else
          innermatrix_ptr = nullptr;

        if (this->GetFESpace()->IsParallel())
          {
            harmonicext = make_shared<ParallelMatrix> (harmonicext,
                                                       this->GetTrialSpace()->GetParallelDofs(),
                                                       this->GetTrialSpace()->GetParallelDofs(), C2C);
            harmonicexttrans = make_shared<ParallelMatrix> (harmonicexttrans,
                                                            this->GetTestSpace()->GetParallelDofs(),
                                                            this->GetTestSpace()->GetParallelDofs(), D2D);
            innersolve = make_shared<ParallelMatrix> (innersolve,
                                                      this->GetTestSpace()->GetParallelDofs(),
                                                      this->GetTrialSpace()->GetParallelDofs(), D2C);
            if (innermatrix)
              innermatrix = make_shared<ParallelMatrix> (innermatrix,
                                                         this->GetTrialSpace()->GetParallelDofs(),
                                                         this->GetTestSpace()->GetParallelDofs(), C2D);
          }
        if(harmonicext)
          GetMemoryTracer().Track(*harmonicext, "HarmonicExt",
                                  *harmonicexttrans, "HarmonicExtTrans");
        if(innermatrix)
          GetMemoryTracer().Track(*innermatrix, "InnerMatrix");
      }
  }


  template <class SCAL>
  void S_BilinearForm<SCAL> :: Assemble_facetwise_skeleton_parts_VOL (Array<bool>& useddof, size_t & gcnt, LocalHeap & clh, const BaseVector * lin)
  {
    size_t ne = ma->GetNE(VOL);
    size_t nf = ma->GetNFacets();
    BitArray fine_facet(nf);
    fine_facet.Clear();
    Array<int> elfacets;
                
    for (int i = 0; i < ne; ++i)
    {
      auto elfacets = ma->GetElFacets(ElementId(VOL,i));
      for (auto f : elfacets) fine_facet.SetBit(f);
    }
                
    ProgressOutput progress(ma,string("assemble inner facet"), nf);
    for (auto colfacets : fespace->FacetColoring())
    {
      SharedLoop2 sl(colfacets.Size());
      ParallelJob
        ( [&] (const TaskInfo & ti) 
        {
          LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
          Array<int> elnums_per(2, lh);
          Array<int> dnums, dnums1, dnums2, elnums, fnums, vnums1, vnums2;
          for (int il : sl)
          {
            int i = colfacets[il];
            progress.Update();
            if (!fine_facet.Test(i)) continue;
            HeapReset hr(lh);
                                  
            int el1 = -1, el2 = -1;
                          
            ma->GetFacetElements(i,elnums);

            el1 = elnums[0];

            int fac2 = i;
            // timerDG1.Stop();
            if(elnums.Size() < 2)
            {
              int facet2 = ma->GetPeriodicFacet(i);
              if(facet2 > i)
              {
                ma->GetFacetElements (facet2, elnums_per);
                elnums.Append(elnums_per[0]);
                fac2 = facet2;
              }
              else
                continue;
            }
                          
            if(elnums.Size()<2) continue;
                          
            el2 = elnums[1];
                                  
            ElementId ei1(VOL, el1);
            ElementId ei2(VOL, el2);
                          
            fnums = ma->GetElFacets(ei1);
            int facnr1 = fnums.Pos(i);
                          
            fnums = ma->GetElFacets(ei2);
            int facnr2 = fnums.Pos(fac2);
                          
            const FiniteElement & fel1 = fespace->GetFE (ei1, lh);
            const FiniteElement & fel2 = fespace->GetFE (ei2, lh);
                          
            ElementTransformation & eltrans1 = ma->GetTrafo (ei1, lh);
            ElementTransformation & eltrans2 = ma->GetTrafo (ei2, lh);
                          
                          
            fespace->GetDofNrs (ei1, dnums1);
            dnums=dnums1;
            fespace->GetDofNrs (ei2, dnums2);
            dnums.Append(dnums2);
                          
            vnums1 = ma->GetElVertices (ei1);
            vnums2 = ma->GetElVertices (ei2);
            if(fel1.GetNDof() != dnums1.Size() || ((elnums.Size()>1) && (fel2.GetNDof() != dnums2.Size() )))
            {
              cout << "facet, neighbouring fel(1): GetNDof() = " << fel1.GetNDof() << endl;
              cout << "facet, neighbouring fel(2): GetNDof() = " << fel2.GetNDof() << endl;
              cout << "facet, neighbouring fel(1): dnums.Size() = " << fel1.GetNDof() << endl;
              cout << "facet, neighbouring fel(2): dnums.Size() = " << fel2.GetNDof() << endl;
              throw Exception ( "Inconsistent number of degrees of freedom " );
            }

            FlatVector<SCAL> elveclin(dnums.Size() * fespace->GetDimension(),lh);
            if (lin)
            {
              lin->GetIndirect(dnums,elveclin);
              int n1 = dnums1.Size()*fespace->GetDimension();
              FlatVector<SCAL> elveclin1(n1,&elveclin(0));
              FlatVector<SCAL> elveclin2(n1,&elveclin(n1));
              fespace->TransformVec(ei1,elveclin1,TRANSFORM_SOL);
              fespace->TransformVec(ei2,elveclin2,TRANSFORM_SOL);
            }
            
            for (auto & bfi : facetwise_skeleton_parts[VOL])
            {
              // if (!bfi->SkeletonForm()) continue;
              // if (bfi->VB() == BND) continue;
              if (!bfi->DefinedOn (ma->GetElIndex (ei1))) continue;
              if (!bfi->DefinedOn (ma->GetElIndex (ei2))) continue;
              if (!bfi->DefinedOnElement(i)) continue;                        
                              
              if (check_unused)
                for (auto d : dnums)
                  if (IsRegularDof(d)) useddof[d] = true;
                              
              int elmat_size = (dnums1.Size()+dnums2.Size())*fespace->GetDimension();
              FlatMatrix<SCAL> elmat(elmat_size, lh);

              auto & mapped_trafo1 = eltrans1.AddDeformation(bfi->GetDeformation().get(), lh);
              auto & mapped_trafo2 = eltrans2.AddDeformation(bfi->GetDeformation().get(), lh);

              if (lin)
              {
                // bfi->CalcFacetMatrix (fel1,facnr1,eltrans1,vnums1,
                //                       fel2,facnr2,eltrans2,vnums2, elmat, lh);
                bfi->CalcLinearizedFacetMatrix (fel1,facnr1,mapped_trafo1,vnums1, 
                                                fel2,facnr2,mapped_trafo2,vnums2, elveclin, elmat, lh);
              }
              else                
              {
                bfi->CalcFacetMatrix (fel1,facnr1,mapped_trafo1,vnums1,
                                      fel2,facnr2,mapped_trafo2,vnums2, elmat, lh);
              }
                              
              fespace->TransformMat (ei1, elmat.Rows(0,dnums1.Size()*fespace->GetDimension()), TRANSFORM_MAT_LEFT);
              fespace->TransformMat (ei2, elmat.Rows(dnums1.Size()*fespace->GetDimension(),elmat_size), TRANSFORM_MAT_LEFT);
              fespace->TransformMat (ei1, elmat.Cols(0,dnums1.Size()*fespace->GetDimension()), TRANSFORM_MAT_RIGHT);
              fespace->TransformMat (ei2, elmat.Cols(dnums1.Size()*fespace->GetDimension(),elmat_size), TRANSFORM_MAT_RIGHT);
                              
              if (printelmat)
              {
                testout->precision(8);
                                  
                (*testout) << "facet-elnum= " << i << endl;
                (*testout) << "integrator " << bfi->Name() << endl;
                (*testout) << "dnums1 = " << endl << dnums1 << endl;
                (*testout) << "dnums2 = " << endl << dnums2 << endl;
                (*testout) << "element1-index = " << eltrans1.GetElementIndex() << endl;
                (*testout) << "element2-index = " << eltrans2.GetElementIndex() << endl;
                (*testout) << "elmat = " << endl << elmat << endl;
              }
                              
              ArrayMem<int, 50> dnums;
              dnums.SetSize(0);
              dnums.Append(dnums1);
              dnums.Append(dnums2);
                              
              ArrayMem<int, 50> map(dnums.Size());
              for (int i = 0; i < map.Size(); i++) map[i] = i;
              QuickSortI (dnums, map);
                              
              ArrayMem<int,50> compressed_dnums;
              compressed_dnums.SetSize(0);
                              
              ArrayMem<int,50> dnums_to_compressed(dnums.Size());
              int compressed_dofs = 0;
              for (int i = 0; i < dnums.Size(); ++i)
              {
                if (i==0 || (dnums[map[i]] != dnums[map[i-1]]))
                {
                  compressed_dnums.Append(dnums[map[i]]);
                  dnums_to_compressed[map[i]] = compressed_dofs++;
                }
                else
                {
                  dnums_to_compressed[map[i]] = dnums_to_compressed[map[i-1]];
                }
              }
                              
              FlatMatrix<SCAL> compressed_elmat(compressed_dofs * fespace->GetDimension(), lh);
              compressed_elmat = 0.0;
              for (int i = 0; i < dnums.Size(); ++i)
                for (int j = 0; j < dnums.Size(); ++j)
                  compressed_elmat(dnums_to_compressed[i],dnums_to_compressed[j]) += elmat(i,j);
                              
              if (elmat_ev)
              {
                testout->precision(8);
                                  
                (*testout) << "elind1 = " << eltrans1.GetElementIndex() << endl;
                (*testout) << "elind2 = " << eltrans2.GetElementIndex() << endl;
#ifdef LAPACK
                LapackEigenSystem(compressed_elmat, lh);
#else
                Vector<SCAL> lami(compressed_elmat.Height());
                Matrix<SCAL> evecs(compressed_elmat.Height());
                                  
                CalcEigenSystem (compressed_elmat, lami, evecs);
                (*testout) << "lami = " << endl << lami << endl;
#endif
                // << "evecs = " << endl << evecs << endl;
              }
                              
              //                    for(int k=0; k<elmat.Height(); k++)
              //                      if(fabs(elmat(k,k)) < 1e-7 && dnums[k] != -1)
              //                        cout << "dnums " << dnums << " elmat " << elmat << endl; 
                              
              AddElementMatrix (compressed_dnums, compressed_dnums, compressed_elmat, ElementId(BND,i), false, lh);
            }
          }
        });
    }
    progress.Done();
    gcnt += nf;
  }

  

  template <class SCAL>
  void S_BilinearForm<SCAL> :: DoAssemble (LocalHeap & clh)
  {
    static Timer mattimer("Matrix assembling");

    static Timer mattimer_checkintegrators("Matrix assembling check integrators");
    static Timer mattimer1a("Matrix assembling initializing");
    static Timer mattimer_finalize("Matrix assembling finalize matrix");
    static Timer<> mattimer_VB[] = { Timer("Matrix assembling vol"),
                                   Timer("Matrix assembling bound"),
                                   Timer("Matrix assembling co dim 2") };
    
    static mutex addelemfacbnd_mutex;
    static mutex addelemfacin_mutex;
    static mutex printelmat_mutex;
    static mutex printmatasstatus2_mutex;
    static mutex printmatspecel_mutex;
    static mutex printmatspecel2_mutex;

    RegionTimer reg (mattimer);

    timestamp = ++global_timestamp;
    
    mattimer_checkintegrators.Start();
    // check if integrators fit to space
    for(VorB vb : {VOL,BND,BBND})
      {
	for(auto bfi : VB_parts[vb])
	  {
	    IterateElements
	      (*fespace,vb,clh, [&] (FESpace::Element el, LocalHeap & lh)
	       {
		 if(bfi->DefinedOn(el.GetIndex()) && bfi->DefinedOnElement(el.Nr()))
		   bfi->CheckElement(el.GetFE());
	       });
	    if(bfi->VB()==VOL && bfi->SkeletonForm())
	      if (!fespace->UsesDGCoupling())
		throw Exception("FESpace is not suitable for those integrators (try -dgjumps)");
	  }
      }
    mattimer_checkintegrators.Stop();    
    
    try
      {
        if (!MixedSpaces())
	  
          {
	    mattimer1a.Start();

            ma->PushStatus ("Assemble Matrix");
 
            size_t ndof = fespace->GetNDof();
            Array<bool> useddof(ndof);
            if (check_unused)
              useddof = false;

            size_t nf = ma->GetNFacets();

            GetMatrix().SetZero();
	    
            if (print)
              {
                *testout << " BILINEARFORM TEST:" << endl;
                *testout << " hasinner = " << (VB_parts[VOL].Size()>0) << endl;
                *testout << " hasouter = " << (VB_parts[BND].Size()>0) << endl;
		*testout << " hasbbnd = " << (VB_parts[BBND].Size()>0) << endl;
                *testout << " hasskeletoninner = " << (facetwise_skeleton_parts[VOL].Size()>0) << endl;
                *testout << " hasskeletonouter = " << (facetwise_skeleton_parts[BND].Size()>0) << endl;
              }
            
            size_t loopsteps = 0;
	    for (VorB vb : {VOL,BND,BBND})
	      if (VB_parts[vb].Size())
		loopsteps += ma->GetNE(vb);
            
	    if (facetwise_skeleton_parts[VOL].Size())
	      loopsteps += ma->GetNFacets();
            
	    if (facetwise_skeleton_parts[BND].Size())
	      loopsteps += ma->GetNSE();
            
	    if (specialelements.Size())
	      loopsteps += specialelements.Size(); 
	    
            size_t gcnt = 0; //global count (for all cases)
            
            for (auto pre : preconditioners)
              pre -> InitLevel(fespace->GetFreeDofs(eliminate_internal));

	    mattimer1a.Stop();

	    for (VorB vb : { VOL, BND, BBND })
	      {
		if (!VB_parts[vb].Size()) continue;
                
                RegionTimer reg(mattimer_VB[vb]);
		size_t ne = ma->GetNE(vb);
                
                if (vb==VOL && diagonal)
                  {
                    static Timer tdiag("AssembleDiag"); RegionTimer reg(tdiag);
                    
                    // double prevtime = WallTime();
                    // Array<DofId> dnums;

                    // for (int i = 0; i < ne; i++)
                    IterateElements
                      (*fespace, vb, clh,  [&] (FESpace::Element el, LocalHeap & lh)
                      {
                        // HeapReset hr(clh);
                        // ElementId id(vb, i);
                        // gcnt++;
			/*
                        if (clock()-prevtime > 0.1 * CLOCKS_PER_SEC)
                          {
                            cout << IM(3) << "\rassemble element " << i << "/" << ne << flush;
                            ma->SetThreadPercentage ( 100.0*gcnt / (loopsteps) );
                            prevtime = clock();
                          }
                        */
                        // clh.CleanUp(heapp);
			
                        // if (!fespace->DefinedOn (vb,ma->GetElIndex (id))) continue;
			
                        const FiniteElement & fel = fespace->GetFE (el, lh);
                        ElementTransformation & eltrans = ma->GetTrafo (el, lh);
                        // auto ei = ElementId(vb,i);
                        // fespace->GetDofNrs (ei, dnums);
                        FlatArray<int> dnums = el.GetDofs();
                        
                        FlatVector<SCAL> sum_diag(dnums.Size()*fespace->GetDimension(), lh);
                        sum_diag = 0;
                        
                        for (auto & bfip : VB_parts[vb])
                          {
                            const BilinearFormIntegrator & bfi = *bfip;
                            if (!bfi.DefinedOn (el.GetIndex())) continue;
                            if (!bfi.DefinedOnElement (el.Nr())) continue;
                            
                            FlatVector<double> diag(sum_diag.Size(), lh);
                            try
                              {
                                bfi.CalcElementMatrixDiag (fel, eltrans, diag, lh);
				
                                if (printelmat) 
                                  {
                                    lock_guard<mutex> guard(printelmat_mutex);                                    
                                    testout->precision(8);
                                    (*testout) << "elnum= " << el << endl;
                                    (*testout) << "integrator " << bfi.Name() << endl;
                                    (*testout) << "dnums = " << endl << dnums << endl;
                                    (*testout) << "diag-elmat = " << endl << diag << endl;
                                  }
                              }
                            catch (Exception & e)
                              {
                                e.Append (string(" in Assemble Element Matrix, bfi = ") + 
                                          bfi.Name() + string("\n"));
                                throw;
                              }
                            catch (exception & e)
                              {
                                throw (Exception (string(e.what()) +
                                                  string(" in Assemble Element Matrix, bfi = ") + 
                                                  bfi.Name() + string("\n")));
                              }
                            
                            sum_diag += diag;
                          }
			
                        AddDiagElementMatrix (dnums, sum_diag, 1, el.Nr(), lh);

                        if (check_unused)
                          for (int j = 0; j < dnums.Size(); j++)
                            if (IsRegularDof(dnums[j]))
                              useddof[dnums[j]] = true;
                      });
                    // cout << IM(3) << "\rassemble element " << ne << "/" << ne << endl;
                  }
                else // not diagonal
                  {
                    ProgressOutput progress(ma,string("assemble ") + ToString(vb) + string(" element"), ma->GetNE(vb));
                    /*
                    if ( (vb == VOL || (!VB_parts[VOL].Size() && vb==BND) ) && eliminate_internal && keep_internal)
                      {
                        harmonicext = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
                        if (!symmetric)
                          harmonicexttrans = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
                        else
                          harmonicexttrans = make_shared<Transpose>(*harmonicext);
                        innersolve = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
                        if (store_inner)
                          innermatrix = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
                      }
                    */
                    IterateElements
                      (*fespace, vb, clh,  [&] (FESpace::Element el, LocalHeap & lh)
                       {
                         if (elmat_ev && vb == VOL) 
                           *testout << " Assemble Element " << el.Nr() << endl;  
                         
                         progress.Update ();
			 
                         const FiniteElement & fel = fespace->GetFE (el, lh);
                         const ElementTransformation & eltrans = ma->GetTrafo (el, lh);
                         FlatArray<int> dnums = el.GetDofs();
                         
                         if (fel.GetNDof() != dnums.Size())
                           {
                             *testout << "Info from finite element: " << endl;
                             fel.Print (*testout);
                             (*testout) << "fel::GetNDof() = " << fel.GetNDof() << endl;
                             (*testout) << "dnums.Size() = " << dnums.Size() << endl;
                             (*testout) << "dnums = " << dnums << endl;
                             throw Exception ( string("Inconsistent number of degrees of freedom, vb="+ToString(vb)+" fel::GetNDof() = ") + ToString(fel.GetNDof()) + string(" != dnums.Size() = ") + ToString(dnums.Size()) + string("!") );
                           }
                         
                         int elmat_size = dnums.Size()*fespace->GetDimension();
                         FlatMatrix<SCAL> sum_elmat(elmat_size, lh);
			 bool elem_has_integrator = false;

                         {
                         static Timer elmattimer("calc elmats", NoTracing);
                         RegionTimer reg (elmattimer);
                         
                         if (printelmat || elmat_ev)
                           {
                             // need every part of the element matrix
                             sum_elmat = 0;
                             for (auto & bfip : VB_parts[vb])
                               {
                                 const BilinearFormIntegrator & bfi = *bfip;
                                 if (!bfi.DefinedOn (el.GetIndex())) continue;                        
                                 if (!bfi.DefinedOnElement (el.Nr())) continue;                        
                                 
                                 elem_has_integrator = true;
                                 
                                 HeapReset hr(lh);
                                 FlatMatrix<SCAL> elmat(elmat_size, lh);
                                 
                                 try
                                   {
                                     bfi.CalcElementMatrix (fel, eltrans, elmat, lh);
                                     
                                     if (printelmat)
                                       {
                                         lock_guard<mutex> guard(printelmat_mutex);
                                         testout->precision(8);
                                         *testout << "elnum = " << el << endl;
                                         *testout << "eltype = " << fel.ElementType() << endl;
                                         *testout << "integrator = " << bfi.Name() << endl;
                                         *testout << "dnums = " << endl << dnums << endl;
                                         *testout << "ct = ";
                                         for (auto d : dnums)
                                           if (!IsRegularDof(d)) *testout << "0 ";
                                           else *testout << fespace->GetDofCouplingType (d) << " ";
                                         *testout << endl;
                                         *testout << "element-index = " << eltrans.GetElementIndex() << endl;
                                         *testout << "elmat = " << endl << elmat << endl;
                                       }
                                     
                                     if (elmat_ev)
                                       LapackEigenSystem(elmat, lh);
                                   }
                                 catch (exception & e)
                                   {
                                     throw (Exception (string(e.what()) +
                                                       string("in Assemble Element Matrix, bfi = ") + 
                                                       bfi.Name() + string("\n")));
                                   }
                                 
                                 sum_elmat += elmat;
                               }
                           }
                         else
                           {
                             /*
                             for (auto & bfip : VB_parts[vb])
                               {
                                 const BilinearFormIntegrator & bfi = *bfip;
                                 if (!bfi.DefinedOn (el.GetIndex())) continue;                        
                                 if (!bfi.DefinedOnElement (el.Nr())) continue;                        
                                 
                                 elem_has_integrator = true;
                                 
                                 try
                                   {
                                     bfi.CalcElementMatrixAdd (fel, eltrans, sum_elmat, lh);
                                   }
                                 catch (exception & e)
                                   {
                                     throw (Exception (string(e.what()) +
                                                       string("in Assemble Element Matrix, bfi = ") + 
                                                       bfi.Name() + string("\n")));
                                   }
                               }
                             */
                             bool done = false;
                             while (!done)
                               {
                                 done = true;
                                 sum_elmat = 0;
                                 bool symmetric_so_far = true;
                                 for (auto & bfip : VB_parts[vb])
                                   {
                                     const BilinearFormIntegrator & bfi = *bfip;
                                     if (!bfi.DefinedOn (el.GetIndex())) continue;                        
                                     if (!bfi.DefinedOnElement (el.Nr())) continue;                        

                                     elem_has_integrator = true;
                                     
                                     try
                                       {
                                         // should we give an optional derformation to the integrators ? 
                                         auto & mapped_trafo = eltrans.AddDeformation(bfi.GetDeformation().get(), lh);
                                         bfi.CalcElementMatrixAdd (fel, mapped_trafo, sum_elmat, symmetric_so_far, lh);
                                       }
                                     catch (ExceptionNOSIMD & e)
                                       {
                                         done = false;
                                       }
                                   }
                               }
                           }
                         } 
                         
                         if (!elem_has_integrator) return;
                         
                         fespace->TransformMat (el, sum_elmat, TRANSFORM_MAT_LEFT_RIGHT);
			 
                         if (elmat_ev)
                           {
                             (*testout) << "sum matrix:" << endl;
                             LapackEigenSystem(sum_elmat, lh);
                           }

                         /*
                         Array<int> lhdofs(dnums.Size(), lh);
                         fespace->GetElementDofsOfType(el, lhdofs, HIDDEN_DOF);
                         bool has_hidden = lhdofs.Size() > 0;
                         */

                         bool has_hidden = false;
                         if (eliminate_hidden || eliminate_internal)
                           {
                             for (auto d : dnums)
                               if (fespace->GetDofCouplingType(d) & HIDDEN_DOF)
                                 has_hidden = true;
                           }

                         bool elim_only_hidden =
                           (!eliminate_internal) && eliminate_hidden && /* (lhdofs.Size() > 0)*/ has_hidden;
                         if ((vb == VOL || (!VB_parts[VOL].Size() && vb==BND) ) && (elim_only_hidden || eliminate_internal))
                           {
                             // static Timer t("static condensation", NoTracing);
                             // RegionTracer reg(TaskManager::GetThreadId(), t);    

                             // if (!fespace->CouplingTypeArrayAvailable())
                             // throw Exception ("need coupling types for static condensation");
                             static Timer statcondtimer("static condensation", NoTracing);
                             RegionTimer regstat (statcondtimer);
                             static Timer statcondtimer2("static condensation 2", NoTracing);

                             static Timer statcondtimer_mult("static condensation mult", NoTracing);
                             static Timer statcondtimer_inv("static condensation inv", NoTracing);
                             
                             // Array<int> idofs1(dnums.Size(), lh);
                             // fespace->GetElementDofsOfType (el, idofs1, elim_only_hidden ? HIDDEN_DOF : CONDENSABLE_DOF);
                             
                             Array<int> idofs1(dnums.Size(), lh), odofs1(dnums.Size(), lh);
                             idofs1.SetSize0(); odofs1.SetSize0();
                             
                             auto ctype = elim_only_hidden ? HIDDEN_DOF : CONDENSABLE_DOF;
                             for (auto i : Range(dnums))
                               {
                                 auto ct = fespace->GetDofCouplingType(dnums[i]);
                                 if (ct & ctype)
                                   idofs1.AppendHaveMem(i);
                                 else
                                   if (ct != UNUSED_DOF)
                                     odofs1.AppendHaveMem(i);
                               }
                             
                             if (printelmat) 
                               {
                                 lock_guard<mutex> guard(printelmat_mutex);
                                 *testout << "eliminate internal";
                                 if (elim_only_hidden)
                                  *testout << " (only hidden)" << endl;
                                 *testout << "idofs1 = " << idofs1 << endl;
                               }
                             
                             if (idofs1.Size())
                               {
                                 HeapReset hr (lh);
				 
                                 int size = sum_elmat.Height();
                                 // int dim = size / dnums.Size();
                                 int dim = fespace->GetDimension();
				 
                                 int sizei = dim * idofs1.Size();
                                 int sizeo = dim * odofs1.Size();
				 
                                 FlatArray<int> idofs (sizei, lh);
                                 FlatArray<int> odofs (sizeo, lh);

                                 for (int j = 0, k = 0; j < idofs1.Size(); j++)
                                   for (int jj = 0; jj < dim; jj++)
                                     idofs[k++] = dim*idofs1[j]+jj;
                                 
                                 for (int j = 0, k = 0; j < odofs1.Size(); j++)
                                   for (int jj = 0; jj < dim; jj++)
                                     odofs[k++] = dim*odofs1[j]+jj;

                                 if (printelmat)
                                   {
                                     lock_guard<mutex> guard(printelmat_mutex);
                                     (*testout) << "idofs = " << endl << idofs << endl;
                                     (*testout) << "odofs = " << endl << odofs << endl;
                                   }
                                 
                                 FlatMatrix<SCAL> 
                                   a = sum_elmat.Rows(odofs).Cols(odofs) | lh,
                                   b = sum_elmat.Rows(odofs).Cols(idofs) | lh,
                                   c = Trans(sum_elmat.Rows(idofs).Cols(odofs)) | lh,
                                   d = sum_elmat.Rows(idofs).Cols(idofs) | lh;
                                 
				 
                                 /*
                                   statcondtimer.AddFlops (double(sizei)*sizei*sizei/3);  // LU fact
                                   statcondtimer.AddFlops (double(sizei)*sizei*sizeo);  
                                   statcondtimer.AddFlops (double(sizei)*sizeo*sizeo);  
                                 */
				 
                                 // A := A - B D^{-1} C^T
                                 // new Versions, July 07
                                 if (elim_only_hidden || !keep_internal) 
                                   {
                                     AInvBt (d, b);    // b <--- b d^-1
                                     LapackMultAddABt (b, c, -1, a);                                 
                                   }
                                 else
                                   {
                                     /*
                                     Array<int> idnums1(dnums.Size(), lh), 
                                       ednums1(dnums.Size(), lh),
                                       hdnums1(dnums.Size(), lh);
                                     fespace->GetElementDofsOfType(el,idnums1,CONDENSABLE_DOF);
                                     fespace->GetElementDofsOfType(el,ednums1,EXTERNAL_DOF);
                                     fespace->GetElementDofsOfType(el,hdnums1,HIDDEN_DOF);

                                     if (! (ednums1 == odofs1) )
                                       cout << "they are different:" << endl
                                            << "ednums1 = " << endl << ednums1 << endl
                                            << "odofs1 = " << endl << odofs1 << endl;
                                     if (! (idnums1 == idofs1) )
                                       cout << "they are different:" << endl
                                            << "idnums1 = " << endl << ednums1 << endl
                                            << "idofs1 = " << endl << idofs1 << endl;
                                     
                                     for (auto d : Range(idnums1.Size()))
                                       idnums1[d] = dnums[idnums1[d]];
                                     for (auto d : Range(ednums1.Size()))
                                       ednums1[d] = dnums[ednums1[d]];
                                     for (auto ldof : hdnums1)
                                       idnums1[ldof] = NO_DOF_NR;
                                     */

                                     Array<DofId> idnums1(idofs1.Size(), lh), ednums1(odofs1.Size(), lh);
                                     for (int i : Range(idofs1))
                                       {
                                         DofId d = dnums[idofs1[i]];
                                         if (fespace->GetDofCouplingType(d) == HIDDEN_DOF)
                                           d = NO_DOF_NR_CONDENSE;
                                         idnums1[i] = d;
                                       }
                                     for (int i : Range(odofs1))
                                       ednums1[i] = dnums[odofs1[i]];
                                     
                                     
                                     RegionTimer regstat2 (statcondtimer2);

                                     Array<int> idnums(dim*idnums1.Size(), lh);
                                     Array<int> ednums(dim*ednums1.Size(), lh);
                                     idnums.SetSize0(); 
                                     ednums.SetSize0();

                                     for (DofId idnum1 : idnums1)
                                       if (IsRegularDof(idnum1))
                                         idnums += dim*IntRange(idnum1, idnum1+1);
                                       else
                                         for (size_t k = 0; k < dim; k++)
                                           idnums.AppendHaveMem(idnum1);

                                     for (DofId ednum1 : ednums1)
                                       ednums += dim * IntRange(ednum1, ednum1+1);
                                     
                                     if (store_inner)
                                     {
                                       if (has_hidden)
                                       {
                                         HeapReset hr(lh);
                                         Array<int> ldofs(idnums.Size(), lh); //LOCAL DOFs
                                         Array<int> hdofs(idnums.Size(), lh); //HIDDEN DOFs
                                         ldofs.SetSize0(); 
                                         hdofs.SetSize0();

                                         for (int i : Range(idnums))
                                           if (idnums[i] == NO_DOF_NR_CONDENSE)
                                             hdofs.AppendHaveMem(i);
                                           else
                                             ldofs.AppendHaveMem(i);

                                         FlatMatrix<SCAL> 
                                           da = d.Rows(ldofs).Cols(ldofs) | lh,
                                           db = d.Rows(ldofs).Cols(hdofs) | lh,
                                           dc = Trans(d.Rows(hdofs).Cols(ldofs)) | lh,
                                           dd = d.Rows(hdofs).Cols(hdofs) | lh;

                                         Array<int> cidnums(idnums.Size(), lh); //compressed idofs (no hidden)
                                         cidnums.SetSize0();
                                         for (DofId dof : idnums)
                                           if (IsRegularDof(dof))
                                             cidnums.AppendHaveMem(dof);

                                         AInvBt (dd, db);    // b <--- b d^-1
                                         LapackMultAddABt (db, dc, -1, da);   
                             
                                         innermatrix_ptr->AddElementMatrix(el.Nr(),cidnums,cidnums,da);
                                       }
                                       else
                                         innermatrix_ptr->AddElementMatrix(el.Nr(),idnums,idnums,d);
                                     }
                                     
                                     /*
                                       Matrix<SCAL> hd = d;
                                       Vector<SCAL> diag(d.Height());
                                       for (int i = 0; i < d.Height(); i++)
                                       diag(i) = sqrt(fabs(hd(i,i)));
                                       for (int i = 0; i < d.Height(); i++)
                                       {
                                       hd.Row(i) *= 1.0/diag(i);
                                       hd.Col(i) *= 1.0/diag(i);
                                       }
                                       Vector<SCAL> lam(d.Height());
                                       Matrix<SCAL> evecs(d.Height());
                                       CalcEigenSystem (hd, lam, evecs);
                                       cout << "lam = " << lam << endl;
                                     */

                                     {
                                       // RegionTimer reg (statcondtimer_inv);
                                       // RegionTracer rtr(TaskManager::GetThreadId(), statcondtimer_inv);    
                              
                                       // LapackInverse (d);
                                       CalcInverse (d);
                                     }
                                     FlatMatrix<SCAL> he (sizei, sizeo, lh);

                                     {
                                       RegionTimer reg (statcondtimer_mult);
                                       NgProfiler::AddThreadFlops (statcondtimer_mult, TaskManager::GetThreadId(),
                                                                   d.Height()*d.Width()*c.Width());
                                       
                                       // V1:
                                       // he = 0.0;
                                       // he -= d * Trans(c) | Lapack;
                                       // V2:
                                       // he = -d * Trans(c) | Lapack;
                                       // V3:
                                       // MinusMultABt (d, c, he);
                                       he = -d * Trans(c);
                                     }
                                     
                                     harmonicext_ptr->AddElementMatrix(el.Nr(),idnums,ednums,he);
                                     if (!symmetric)
                                       {
                                         FlatMatrix<SCAL> het (sizeo, sizei, lh);
                                         // het = -b*d | Lapack;
                                         // MinusMultAB (b, d, het);
                                         het = -b * d;
                                         harmonicexttrans_ptr->AddElementMatrix(el.Nr(),ednums,idnums,het);
                                       }
                                     
                                     innersolve_ptr->AddElementMatrix(el.Nr(),idnums,idnums,d);
                                     {
                                       RegionTimer reg (statcondtimer_mult);
                                       NgProfiler::AddThreadFlops (statcondtimer_mult, TaskManager::GetThreadId(),
                                                                   b.Height()*b.Width()*he.Width());
                                       // a += b * he | Lapack;
                                       // AddAB (b, he, a);
                                       a += b * he;
                                     }
                                     
                                     if (spd)
                                       { // more stable ? 
                                         FlatMatrix<SCAL> schur(odofs.Size(), lh);
                                         CalcSchur (sum_elmat, schur, odofs, idofs);
                                         a = schur;
                                       }
                                   }                             
                                 
                                 if (printelmat) 
                                   {
                                     testout->precision(8);
                                     (*testout) << "Schur elmat = " << endl << a << endl;
                                   }
                                 
                                 if (elmat_ev)
                                   {
                                     testout->precision(8);
                                     
                                     (*testout) << "EV of Schur complement:" << endl;
                                     LapackEigenSystem(a, lh);
                                   }
                                 
                                 sum_elmat.Rows(odofs).Cols(odofs) = a;
                                 if (linearform && (!keep_internal))
                                   {
                                     FlatVector<SCAL> elvec (size, lh);
                                     linearform -> GetVector().GetIndirect (dnums, elvec);
                                     FlatVector<SCAL> hfi(sizei, lh);
                                     FlatVector<SCAL> hfo(sizeo, lh);
                                     
                                     hfi = elvec(idofs);
                                     hfo = b * hfi;
                                     elvec(odofs) -= hfo;
                                     
                                     linearform->GetVector().SetIndirect (dnums, elvec);
                                   }
                                 
                                 for (int k = 0; k < idofs1.Size(); k++)
                                   dnums[idofs1[k]] = NO_DOF_NR;
                               }
                           }
                         if (printelmat)
                           {
                             lock_guard<mutex> guard(printelmat_mutex);
                             *testout<< "elem " << el << ", elmat = " << endl << sum_elmat << endl;
                           }
                         
                         AddElementMatrix (dnums, dnums, sum_elmat, el, false, lh);
			 
                         for (auto pre : preconditioners)
                           pre -> AddElementMatrix (dnums, sum_elmat, el, lh);

                         if (check_unused)
                           {
                             if (printelmat)
                               *testout << "set these as useddof: " << dnums << endl;
                             for (auto d : dnums)
                               if (IsRegularDof(d)) useddof[d] = true;
                           }
                         // timer3_VB[vb].Stop();
                       });
                    progress.Done();
                    
                    /*
			  if (linearform && keep_internal)
			  {
			  cout << IM(3) << "\rmodifying condensated rhs";
			  
			  linearform -> GetVector() += 
			  GetHarmonicExtensionTrans() * linearform -> GetVector();
			  
			  cout << IM(3) << "\t done" << endl;
			  }
			*/
                    
                    gcnt += ne;
                  }
              }


	    //simplify
	    // for (VorB vb : { VOL, BND })
            // {
            // if (!VB_skeleton_parts[vb].Size()) continue;

            // size_t ne = ma->GetNE(vb);
            /*
                if (vb == VOL)
                  {
                    bool has_element_wise = false;
                    bool has_facet_wise = false;
                    for (int j = 0; j < NumIntegrators(); j++)
                      if (parts[j] -> SkeletonForm())
                        {
                          auto dgform = parts[j] -> GetDGFormulation();
                          if (dgform.element_boundary)
                            has_element_wise = true;
                          if (!dgform.element_boundary && !parts[j]->BoundaryForm())
                            has_facet_wise = true;                        
                        }
            */
            if (facetwise_skeleton_parts[VOL].Size())
              Assemble_facetwise_skeleton_parts_VOL(useddof,gcnt,clh);
            
            // if (has_element_wise)
            if (elementwise_skeleton_parts.Size())
              {
                int cnt = 0;
                int ne = ma->GetNE(VOL);
                ParallelForRange
                  (IntRange(ne), [&] ( IntRange r )
                   {
                     LocalHeap lh = clh.Split();
                     
                     Array<int> dnums, dnums1, dnums2, elnums, elnums_per, fnums1, fnums2, vnums1, vnums2;
                     for (int el1 : r)
                       {
                         ElementId ei1(VOL, el1);
                         fnums1 = ma->GetElFacets(ei1);
                         for (int facnr1 : Range(fnums1))
                           {
                             HeapReset hr(lh);
                             
                             ma->GetFacetElements(fnums1[facnr1],elnums);

                             int facet2 = fnums1[facnr1];
                             // timerDG1.Stop();
                             if(elnums.Size() < 2)
                             {
                               facet2 = ma->GetPeriodicFacet(fnums1[facnr1]);
                               if(facet2 != fnums1[facnr1])
                               {
                                 ma->GetFacetElements (facet2, elnums_per);
                                 elnums.Append(elnums_per[0]);
                               }
                             }

                             
                             if (elnums.Size()<2)
                               {
                                 ma->GetFacetSurfaceElements (fnums1[facnr1], elnums);
                                 int sel = elnums[0];
                                 ElementId sei(BND, sel);
                                 
                                 const FiniteElement & fel = fespace->GetFE (ei1, lh);
                                 vnums1 = ma->GetElVertices (ei1);
                                 vnums2 = ma->GetElVertices (sei);     
                                 
                                 ElementTransformation & eltrans = ma->GetTrafo (ei1, lh);
                                 ElementTransformation & seltrans = ma->GetTrafo (sei, lh);
                                 
                                 fespace->GetDofNrs (ei1, dnums);
                                 if(fel.GetNDof() != dnums.Size())
                                   {
                                     (*testout) << "Surface fel::GetNDof() = " << fel.GetNDof() << endl;
                                     (*testout) << "dnums.Size() = " << dnums.Size() << endl;
                                     (*testout) << "dnums = " << dnums << endl;
                                     throw Exception ( string("Inconsistent number of degrees of freedom Surface fel::GetNDof() = ") + ToString(fel.GetNDof()) + string(" != dnums.Size() = ") + ToString(dnums.Size()) + string("!") );
                                   }

                                 /*
                                 for (int j = 0; j < NumIntegrators(); j++)
                                   {
                                     const BilinearFormIntegrator & bfi = *parts[j];
                                 */
                                 for (auto & bfi : elementwise_skeleton_parts)
                                   {
                                     // if (bfi.VB()!=VOL) continue;
                                     // if (!bfi.SkeletonForm()) continue;
                                     // if (!bfi.GetDGFormulation().element_boundary) continue;
				     if (!bfi->DefinedOn (ma->GetElIndex (ei1))) continue;
                                     if (!bfi->DefinedOnElement (el1)) continue;
                                     
                                     if (check_unused)
                                       for (auto d : dnums)
                                         if (IsRegularDof(d)) useddof[d] = true;
                                     
                                     int elmat_size = dnums.Size()*fespace->GetDimension();
                                     FlatMatrix<SCAL> elmat(elmat_size, lh);
                                     
                                     // dynamic_cast<const FacetBilinearFormIntegrator&>(*bfi).  
                                     auto & mapped_trafo = eltrans.AddDeformation(bfi->GetDeformation().get(), lh);
                                     auto & mapped_strafo = seltrans.AddDeformation(bfi->GetDeformation().get(), lh);

                                     bfi->CalcFacetMatrix (fel,facnr1,mapped_trafo,vnums1, mapped_strafo, vnums2, elmat, lh);
                                     
                                     fespace->TransformMat (ei1, elmat, TRANSFORM_MAT_LEFT_RIGHT);
                                     
                                     if (printelmat)
                                       {
                                         testout->precision(8);
                                         
                                         (*testout) << "surface-elnum= " << sel << endl;
                                         (*testout) << "integrator " << bfi->Name() << endl;
                                         (*testout) << "dnums = " << endl << dnums << endl;
                                         (*testout) << "element-index = " << eltrans.GetElementIndex() << endl;
                                         (*testout) << "elmat = " << endl << elmat << endl;
                                       }
                                     
                                     {
                                       // lock_guard<mutex> guard(addelemfacbnd_mutex);
                                       AddElementMatrix (dnums, dnums, elmat, ElementId(VOL,el1), true, lh);
                                     }
                                   } //end for (numintegrators)
                                 continue;
                               } // end if boundary facet
                             
                             
                             
                             int el2 = elnums[0] + elnums[1] - el1;
                             ElementId ei2(VOL, el2);
                             fnums2 = ma->GetElFacets(ei2);
                             int facnr2 = fnums2.Pos(facet2);
                             
                             {
                               lock_guard<mutex> guard(printmatasstatus2_mutex);
                               cnt++;
                               gcnt++;
                               if (cnt % 10 == 0)
                                 cout << IM(3) << "\rassemble inner facet element " << cnt << "/" << nf << flush;
                               ma->SetThreadPercentage ( 100.0*(gcnt) / (loopsteps) );
                             }
                             
                             const FiniteElement & fel1 = fespace->GetFE (ei1, lh);
                             const FiniteElement & fel2 = fespace->GetFE (ei2, lh);
                             
                             ElementTransformation & eltrans1 = ma->GetTrafo (ei1, lh);
                             ElementTransformation & eltrans2 = ma->GetTrafo (ei2, lh);
                             
                             
                             fespace->GetDofNrs (ei1, dnums1);
                             dnums=dnums1;
                             fespace->GetDofNrs (ei2, dnums2);
                             dnums.Append(dnums2);
                             
                             vnums1 = ma->GetElVertices (ei1);
                             vnums2 = ma->GetElVertices (ei2);
                             if(fel1.GetNDof() != dnums1.Size() || ((elnums.Size()>1) && (fel2.GetNDof() != dnums2.Size() )))
                               {
                                 cout << "facet, neighbouring fel(1): GetNDof() = " << fel1.GetNDof() << endl;
                                 cout << "facet, neighbouring fel(2): GetNDof() = " << fel2.GetNDof() << endl;
                                 cout << "facet, neighbouring fel(1): dnums.Size() = " << fel1.GetNDof() << endl;
                                 cout << "facet, neighbouring fel(2): dnums.Size() = " << fel2.GetNDof() << endl;
                                 throw Exception ( "Inconsistent number of degrees of freedom " );
                               }
                             /*
                             for (int j = 0; j < NumIntegrators(); j++)
                               {
                                 shared_ptr<BilinearFormIntegrator> bfi = parts[j];
                             */
                             for (auto & bfi : elementwise_skeleton_parts)
                               {
                                 // if (!bfi->SkeletonForm()) continue;
                                 // if (bfi->VB() != VOL) continue;
                                 if (!bfi->DefinedOn (ma->GetElIndex (ei1))) continue;
                                 if (!bfi->DefinedOn (ma->GetElIndex (ei2))) continue;
                                 if (!bfi->DefinedOnElement (el1)) continue;

                                 if (check_unused)
                                   for (auto d : dnums)
                                     if (IsRegularDof(d)) useddof[d] = true;
                                 
                                 int elmat_size = (dnums1.Size()+dnums2.Size())*fespace->GetDimension();
                                 FlatMatrix<SCAL> elmat(elmat_size, lh);
                                 
                                 // shared_ptr<FacetBilinearFormIntegrator> fbfi = 
                                 // dynamic_pointer_cast<FacetBilinearFormIntegrator>(bfi);
                                 
                                 auto & mapped_trafo1 = eltrans1.AddDeformation(bfi->GetDeformation().get(), lh);
                                 auto & mapped_trafo2 = eltrans2.AddDeformation(bfi->GetDeformation().get(), lh);
                                 // if (fbfi)
                                   {
                                     bfi->CalcFacetMatrix (fel1,facnr1,mapped_trafo1,vnums1,
                                                           fel2,facnr2,mapped_trafo2,vnums2, elmat, lh);
                                   }
                                   /*
                                     // currently CompoundBFI(FactBFI) not possible .... do we need it ? 
                                 else
                                   {
                                     shared_ptr<CompoundBilinearFormIntegrator> cbfi = 
                                       dynamic_pointer_cast<CompoundBilinearFormIntegrator>(bfi);
                                     
                                     if (!cbfi)
                                       throw Exception("neither compound nor facetbilinearformintegrator!");
                                     
                                     fbfi = dynamic_pointer_cast<FacetBilinearFormIntegrator>(cbfi->GetBFI());
                                     
                                     if (!fbfi)
                                       throw Exception("not a FacetBFI inside CompoundBFI!");
                                     
                                     int comp = cbfi->GetComponent();
                                     const CompoundFiniteElement & cfel1 =
                                       dynamic_cast<const CompoundFiniteElement&> (fel1);
                                     const CompoundFiniteElement & cfel2 =
                                       dynamic_cast<const CompoundFiniteElement&> (fel2);
                                     
                                     FlatMatrix<double> mat1(cfel1[comp].GetNDof()+cfel2[comp].GetNDof(), lh);
                                     fbfi->CalcFacetMatrix (
                                                            cfel1[comp], facnr1, eltrans1, vnums1,
                                                            cfel2[comp], facnr2, eltrans2, vnums2, mat1, lh);
                                     
                                     IntRange range1 = cfel1.GetRange (comp);
                                     IntRange range2 = cfel2.GetRange (comp) + cfel1.GetNDof();
                                     
                                     int nd1 = cfel1[comp].GetNDof();
                                     int nd2 = nd1 + cfel2[comp].GetNDof();
                                     
                                     elmat = 0.0;
                                     elmat.Rows (range1).Cols(range1) = mat1.Rows(0,nd1).Cols(0,nd1);
                                     elmat.Rows (range1).Cols(range2) = mat1.Rows(0,nd1).Cols(nd1,nd2);
                                     elmat.Rows (range2).Cols(range1) = mat1.Rows(nd1,nd2).Cols(0,nd1);
                                     elmat.Rows (range2).Cols(range2) = mat1.Rows(nd1,nd2).Cols(nd1,nd2);
                                   }
                                   */
                                 // *testout << "elmat : \n" << elmat << endl;
                                 
                                 fespace->TransformMat (ei1, elmat.Rows(0,dnums1.Size()*fespace->GetDimension()), TRANSFORM_MAT_LEFT);
                                 fespace->TransformMat (ei2, elmat.Rows(dnums1.Size()*fespace->GetDimension(),elmat_size), TRANSFORM_MAT_LEFT);
                                 fespace->TransformMat (ei1, elmat.Cols(0,dnums1.Size()*fespace->GetDimension()), TRANSFORM_MAT_RIGHT);
                                 fespace->TransformMat (ei2, elmat.Cols(dnums1.Size()*fespace->GetDimension(),elmat_size), TRANSFORM_MAT_RIGHT);
                                 
                                 if (printelmat)
                                   {
                                     testout->precision(8);
                                             
                                     (*testout) << "facet-elnum= " << fnums1[facnr1] << endl;
                                     (*testout) << "integrator " << bfi->Name() << endl;
                                     (*testout) << "dnums1 = " << endl << dnums1 << endl;
                                     (*testout) << "dnums2 = " << endl << dnums2 << endl;
                                     (*testout) << "element1-index = " << eltrans1.GetElementIndex() << endl;
                                     (*testout) << "element2-index = " << eltrans2.GetElementIndex() << endl;
                                     (*testout) << "elmat = " << endl << elmat << endl;
                                   }
                                 
                                 Array<int> dnums;
                                 dnums.SetSize(0);
                                 dnums.Append(dnums1);
                                 dnums.Append(dnums2);
                                 
                                 ArrayMem<int, 50> map(dnums.Size());
                                 for (int i = 0; i < map.Size(); i++) map[i] = i;
                                 QuickSortI (dnums, map);
                                 
                                 Array<int> compressed_dnums;
                                 compressed_dnums.SetSize(0);
                                 
                                 Array<int> dnums_to_compressed(dnums.Size());
                                 int compressed_dofs = 0;
                                 for (int i = 0; i < dnums.Size(); ++i)
                                   {
                                     if (i==0 || (dnums[map[i]] != dnums[map[i-1]]))
                                       {
                                         compressed_dnums.Append(dnums[map[i]]);
                                         dnums_to_compressed[map[i]] = compressed_dofs++;
                                       }
                                     else
                                       {
                                         dnums_to_compressed[map[i]] = dnums_to_compressed[map[i-1]];
                                       }
                                   }
                                 
                                 FlatMatrix<SCAL> compressed_elmat(compressed_dofs * fespace->GetDimension(), lh);
                                 compressed_elmat = 0.0;

                                 int dim = fespace->GetDimension();

                                 for (int i = 0; i < dnums.Size(); ++i)
                                   for (int j = 0; j < dnums.Size(); ++j)
                                     for (int di = 0; di < dim; ++di)
                                       for (int dj = 0; dj < dim; ++dj)
                                         compressed_elmat(dim*dnums_to_compressed[i]+di,dim*dnums_to_compressed[j]+dj) += elmat(i*dim+di,j*dim+dj);

                                 {
                                   // lock_guard<mutex> guard(addelemfacin_mutex);
                                   AddElementMatrix (compressed_dnums, compressed_dnums, compressed_elmat,
                                                     ElementId(BND,el1), true, lh);
                                 }
                               }
                           }
                       }                             
                   });
                cout << IM(3) << "\rassemble inner facet element " << nf << "/" << nf << endl;
              } // if (elementwise_skeleton_parts.Size())
                
            if (facetwise_skeleton_parts[BND].Size())
              {
                // cout << "check bnd" << endl;
                int cnt = 0;
                int ne = ma->GetNE(BND);
                ParallelForRange
                  ( IntRange(ne), [&] ( IntRange r )
                    {
                      LocalHeap lh = clh.Split();
                      Array<int> fnums, elnums, vnums, svnums, dnums;
                      
                      for (int i : r)
                        {
                          {
                            lock_guard<mutex> guard(printmatasstatus2_mutex);
                            cnt++;
                            gcnt++;
                            // if (cnt % 10 == 0)
                              // cout << "\rassemble facet surface element " << cnt << "/" << ne << flush;
                            ma->SetThreadPercentage ( 100.0*(gcnt) / (loopsteps) );
                          }
                          
                          HeapReset hr(lh);
                          ElementId sei(BND, i);
                              
                          // if (!fespace->DefinedOn (BND,ma->GetElIndex (sei))) continue;
                          fnums = ma->GetElFacets(sei);
                          int fac = fnums[0];
                          ma->GetFacetElements(fac,elnums);
                          int el = elnums[0];
                          ElementId ei(VOL, el);
                          fnums = ma->GetElFacets(ei);
                          const FiniteElement & fel = fespace->GetFE (ei, lh);
                          int facnr = 0;
                          for (int k=0; k<fnums.Size(); k++)
                            if(fac==fnums[k]) facnr = k;
                          vnums = ma->GetElVertices (ei);
                          svnums = ma->GetElVertices (sei);     
                              
                          ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
                          ElementTransformation & seltrans = ma->GetTrafo (sei, lh);
                              
                          fespace->GetDofNrs (ei, dnums);
                          if(fel.GetNDof() != dnums.Size())
                            {
                              (*testout) << "Surface fel::GetNDof() = " << fel.GetNDof() << endl;
                              (*testout) << "dnums.Size() = " << dnums.Size() << endl;
                              (*testout) << "dnums = " << dnums << endl;
                              throw Exception ( string("Inconsistent number of degrees of freedom Surface fel::GetNDof() = ") + ToString(fel.GetNDof()) + string(" != dnums.Size() = ") + ToString(dnums.Size()) + string("!") );
                            }

                          /*
                            for (int j = 0; j < NumIntegrators(); j++)
                            {
                            const BilinearFormIntegrator & bfi = *parts[j];
                          */
                          for (auto & bfi : facetwise_skeleton_parts[BND])
                            {
                              // if (bfi.VB() != BND) continue;
                              // if (!bfi.SkeletonForm()) continue;
				  
                              if (!bfi->DefinedOn (ma->GetElIndex(sei) )) continue;                
                              if (!bfi->DefinedOnElement (i)) continue;

                              if (check_unused)
                                for (int k = 0; k < dnums.Size(); k++)
                                  if (IsRegularDof(dnums[k]))
                                    useddof[dnums[k]] = true;
                                  
                              int elmat_size = dnums.Size()*fespace->GetDimension();
                              FlatMatrix<SCAL> elmat(elmat_size, lh);
				  
                              // original version did not compile on MacOS V
                              const FacetBilinearFormIntegrator & fbfi = 
                                dynamic_cast<const FacetBilinearFormIntegrator&>(*bfi);  
                              auto & mapped_trafo = eltrans.AddDeformation(fbfi.GetDeformation().get(), lh);
                              auto & mapped_strafo = seltrans.AddDeformation(fbfi.GetDeformation().get(), lh);
                              fbfi.CalcFacetMatrix (fel,facnr,mapped_trafo,vnums, mapped_strafo, svnums, elmat, lh);
				  
                              fespace->TransformMat (ei, elmat, TRANSFORM_MAT_LEFT_RIGHT);
				  
                              if (printelmat)
                                {
                                  testout->precision(8);
                                      
                                  (*testout) << "surface-elnum= " << i << endl;
                                  (*testout) << "integrator " << bfi->Name() << endl;
                                  (*testout) << "dnums = " << endl << dnums << endl;
                                  (*testout) << "element-index = " << eltrans.GetElementIndex() << endl;
                                  (*testout) << "elmat = " << endl << elmat << endl;
                                }
                                  
				  
                              if (elmat_ev)
                                {
                                  testout->precision(8);
                                      
                                  (*testout) << "elind = " << eltrans.GetElementIndex() << endl;
#ifdef LAPACK
                                  LapackEigenSystem(elmat, lh);
#else
                                  Vector<SCAL> lami(elmat.Height());
                                  Matrix<SCAL> evecs(elmat.Height());
                                      
                                  CalcEigenSystem (elmat, lami, evecs);
                                  (*testout) << "lami = " << endl << lami << endl;
#endif
                                  // << "evecs = " << endl << evecs << endl;
                                } 
                                  
                              //                    for(int k=0; k<elmat.Height(); k++)
                              //                      if(fabs(elmat(k,k)) < 1e-7 && dnums[k] != -1)
                              //                        cout << "dnums " << dnums << " elmat " << elmat << endl; 
                              {
                                // lock_guard<mutex> guard(addelemfacbnd_mutex);
                                AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), true, lh);
                              }
                            }//end for (numintegrators)
                        }//end for nse                  
                    });//end of parallel
                // cout << "\rassemble facet surface element " << ne << "/" << ne << endl;  
              } // if facetwise_skeleton_parts[BND].size
            
            ma->SetThreadPercentage ( 100.0 );
            
	    RegionTimer reg(mattimer_finalize);
            
            /*
              if(NumIndependentIntegrators() > 0)
              {
              DoAssembleIndependent(useddof,clh);
              }
            */


            clh.CleanUp();

            bool assembledspecialelements = false;
            

            static Timer tspecial("Special Elements");
            tspecial.Start();
            int nspecel = 0;
            ParallelForRange( IntRange(specialelements.Size()), [&] ( IntRange r )
                              {
              LocalHeap lh = clh.Split();
              Array<int> dnums;
              
              for (int i : r)
                {
                  {
                    lock_guard<mutex> guard(printmatspecel_mutex);
                    gcnt++;
                    nspecel++;
                    if (i % 10 == 0)
                      cout << IM(3) << "\rassemble special element " << nspecel << "/" << specialelements.Size() << flush;
                    ma->SetThreadPercentage ( 100.0*(gcnt) / (loopsteps) );
                  }
                  
                  const SpecialElement & el = *specialelements[i];
                
                  el.GetDofNrs (dnums);
                
                  FlatMatrix<SCAL> elmat(dnums.Size(), lh);
                  el.CalcElementMatrix (elmat, lh);
                  
                  {
                    // lock_guard<mutex> guard(printmatspecel2_mutex);
                    if (check_unused)
                      for (int j = 0; j < dnums.Size(); j++)
                        if (IsRegularDof(dnums[j]))
                          useddof[dnums[j]] = true;
                    
                    AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), true, lh);
                  }
                  
                  assembledspecialelements = true;
                  lh.CleanUp();
                }
            });
            if(assembledspecialelements) cout << IM(3) << "\rassemble special element " 
                                              << specialelements.Size() << "/" << specialelements.Size() << endl;
            tspecial.Stop();

            
            
            // add eps to avoid empty lines
            FlatMatrix<SCAL> elmat (fespace->GetDimension(), clh);
            elmat = 0;
            Array<int> dnums;
            dnums.SetSize(1);
            
            if (eps_regularization != 0)
              {
                for (int i = 0; i < elmat.Height(); i++)
                  elmat(i, i) = eps_regularization;
                for (int i = 0; i < ndof; i++)
                  {
                    dnums[0] = i; 
                    AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), false, clh);
                  }
              }
            if (unuseddiag != 0 && check_unused)
              {
                for (int i = 0; i < elmat.Height(); i++)
                  elmat(i, i) = unuseddiag;

                for (int i = 0; i < ndof; i++)
                  if (!useddof[i])
                    {
                      dnums[0] = i;
                      AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), false, clh);
                    }
              }
            
            for (auto pre : preconditioners)
              pre -> FinalizeLevel(&GetMatrix());
            
            if (print)
              (*testout) << "mat = " << endl << GetMatrix() << endl;

            if (check_unused)
              {
                size_t cntused = 0;
                for (bool u : useddof)
                  if (u) cntused++;
                
                if (cntused < useddof.Size())
                  cout << IM(4) << "used " << cntused
                       << ", unused = " << useddof.Size()-cntused
                       << ", total = " << useddof.Size() << endl;
              }
            
            int MASK = eliminate_internal ? EXTERNAL_DOF : (eliminate_hidden ? VISIBLE_DOF : ANY_DOF);
            bool first_time = true;

            auto comm = ma->GetCommunicator();
            if (comm.Size() == 1 && check_unused) 
              for (int i = 0; i < useddof.Size(); i++)
                if (useddof[i] != 
                    ((fespace->GetDofCouplingType(i) & MASK) != 0) )
                  {
                    *testout << "used dof inconsistency: " 
                             << " dof " << i << " used = " << int(useddof[i])
                             << " ct = " << fespace->GetDofCouplingType(i) << endl;
                  
                    if (first_time)
                      {
                        cerr << "used dof inconsistency" << endl;
                        cerr << "(silence this warning by setting BilinearForm(...check_unused=False) )" << endl;
                      }
                    first_time = false;
                  }
            

            ma->PopStatus ();
          }

        else // MixedSpaces()

          {
            cout << IM(3) << "assemble mixed bilinearform" << endl;
            ma->PushStatus ("Assemble Matrix");
      
            BaseMatrix & mat = GetMatrix();
            mat = 0.0;


            // bool hasbound = false;

	    for (VorB vb : { VOL, BND, BBND, BBBND })
              if (VB_parts[vb].Size())
                IterateElements 
                  (*fespace, vb, clh,          // coloring for 1 space is enough
                   [&] (ElementId ei, LocalHeap & lh)
                   {
                     const FiniteElement & fel1 = fespace->GetFE (ei, lh);
                     const FiniteElement & fel2 = fespace2->GetFE (ei, lh);
                     
                     Array<int> dnums1(fel1.GetNDof(), lh);
                     Array<int> dnums2(fel2.GetNDof(), lh);
                     const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
                     fespace->GetDofNrs (ei, dnums1);
                     fespace2->GetDofNrs (ei, dnums2);
          
                     FlatMatrix<SCAL> elmat(dnums2.Size()*fespace2->GetDimension(),
                                            dnums1.Size()*fespace->GetDimension(), lh);
                     for (auto & bfi : VB_parts[vb])
                       {
                         if (!bfi->DefinedOn (eltrans.GetElementIndex())) continue;
                         if (!bfi->DefinedOnElement (ei.Nr())) continue;
                         
                         auto & mapped_trafo = eltrans.AddDeformation(bfi->GetDeformation().get(), lh);
                         MixedFiniteElement fel(fel1, fel2);
                         bfi->CalcElementMatrix (fel, mapped_trafo, elmat, lh);
                         fespace->TransformMat(ei, elmat, TRANSFORM_MAT_RIGHT);
                         fespace2->TransformMat(ei, elmat, TRANSFORM_MAT_LEFT);
                         AddElementMatrix (dnums2, dnums1, elmat, ei, false, lh);
                       }
                   });
            
            /*
            if (hasbound)
              IterateElements 
                (*fespace, BND, clh,          // coloring for 1 space is enough
                 [&] (ElementId ei, LocalHeap & lh)
                 {
                   const FiniteElement & fel1 = fespace->GetFE (ei, lh);
                   const FiniteElement & fel2 = fespace2->GetFE (ei, lh);
                   
                   Array<int> dnums1(fel1.GetNDof(), lh);
                   Array<int> dnums2(fel2.GetNDof(), lh);
                   const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
                   fespace->GetDofNrs (ei, dnums1);
                   fespace2->GetDofNrs (ei, dnums2);
          
                   FlatMatrix<SCAL> elmat(dnums2.Size(), dnums1.Size(), lh);
                   for (int j = 0; j < NumIntegrators(); j++)
                     {
                       const BilinearFormIntegrator & bfi = *GetIntegrator(j);
                       if (!bfi.BoundaryForm()) continue;

                       // ArrayMem<const FiniteElement*,2> fea = { &fel1, &fel2 };
                       ArrayMem<const FiniteElement*,2> fea(2);
                       fea[0] = &fel1;
                       fea[1] = &fel2;
                       CompoundFiniteElement cfel(fea);

                       bfi.CalcElementMatrix (cfel, eltrans, elmat, lh);

  
                       // fespace->Transform (i, true, elmat, TRANSFORM_MAT_RIGHT);
                       // fespace2->Transform (i, true, elmat, TRANSFORM_MAT_LEFT);
                       
                       AddElementMatrix (dnums2, dnums1, elmat, ei, lh);
                       }
                       });
            */

            /*
            if (VB_parts[BND].Size()) throw Exception ("mixed biforms don't support boundary terms");
            if (VB_parts[BBND].Size()) throw Exception ("mixed biforms don't support bboundary terms");
            if (VB_parts[BBBND].Size()) throw Exception ("mixed biforms don't support bbboundary terms");
            */
            
            //if (facetwise_skeleton_parts[VOL].Size()) throw Exception ("mixed biforms don't support skeleton terms");

            if (facetwise_skeleton_parts[VOL].Size())
            { // inner facets
              // auto vb = VOL;
              size_t nf = ma->GetNFacets();
              size_t ne = ma->GetNE(VOL);
              BitArray fine_facet(nf);
              fine_facet.Clear();
              Array<int> elfacets;
                
              for (int i = 0; i < ne; ++i)
              {
                auto elfacets = ma->GetElFacets(ElementId(VOL,i));
                for (auto f : elfacets) fine_facet.SetBit(f);
              }
                
              ProgressOutput progress(ma,string("assemble inner facet"), nf);
              for (auto colfacets : fespace->FacetColoring())
              {
                SharedLoop2 sl(colfacets.Size());
                ParallelJob
                  ( [&] (const TaskInfo & ti) 
                    {
                      LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
                      Array<int> elnums_per(2, lh);
                      Array<int> dnums_trial, dnums_test, dnums1_trial, dnums1_test, 
                        dnums2_trial, dnums2_test, elnums, fnums, vnums1, vnums2;
                      for (int il : sl)
                      {
                        int i = colfacets[il];
                        progress.Update();
                        if (!fine_facet.Test(i)) continue;
                        HeapReset hr(lh);
                                  
                        int el1 = -1, el2 = -1;
                          
                        ma->GetFacetElements(i,elnums);

                        //if (i == asdf)

                        el1 = elnums[0];

                        int fac2 = i;
                        // timerDG1.Stop();
                        if(elnums.Size() < 2)
                        {
                          int facet2 = ma->GetPeriodicFacet(i);
                          if(facet2 > i)
                          {
                            ma->GetFacetElements (facet2, elnums_per);
                            elnums.Append(elnums_per[0]);
                            fac2 = facet2;
                          }
                          else
                            continue;
                        }
                          
                        if(elnums.Size()<2) continue;
                          
                        el2 = elnums[1];
                                  
                        ElementId ei1(VOL, el1);
                        ElementId ei2(VOL, el2);
                          
                        fnums = ma->GetElFacets(ei1);
                        int facnr1 = fnums.Pos(i);
                          
                        fnums = ma->GetElFacets(ei2);
                        int facnr2 = fnums.Pos(fac2);

                        const FiniteElement & fel1_trial = fespace->GetFE (ei1, lh);
                        const FiniteElement & fel1_test = fespace2->GetFE (ei1, lh);
                        MixedFiniteElement fel1(fel1_trial, fel1_test);
                        const FiniteElement & fel2_trial = fespace->GetFE (ei2, lh);
                        const FiniteElement & fel2_test = fespace2->GetFE (ei2, lh);
                        MixedFiniteElement fel2(fel2_trial, fel2_test);
                        // const FiniteElement & fel1 = fespace->GetFE (ei1, lh);
                        // const FiniteElement & fel2 = fespace->GetFE (ei2, lh);
                          
                        ElementTransformation & eltrans1 = ma->GetTrafo (ei1, lh);
                        ElementTransformation & eltrans2 = ma->GetTrafo (ei2, lh);
                          

                        fespace->GetDofNrs (ei1, dnums1_trial);
                        fespace->GetDofNrs (ei2, dnums2_trial);
                        
                        dnums_trial=dnums1_trial;
                        dnums_trial.Append(dnums2_trial);

                        fespace2->GetDofNrs (ei1, dnums1_test);
                        fespace2->GetDofNrs (ei2, dnums2_test);
                        
                        dnums_test=dnums1_test;
                        dnums_test.Append(dnums2_test);
                        
                        // fespace->GetDofNrs (ei1, dnums1);
                        // dnums=dnums1;
                        // fespace->GetDofNrs (ei2, dnums2);
                        // dnums.Append(dnums2);
                          
                        vnums1 = ma->GetElVertices (ei1);
                        vnums2 = ma->GetElVertices (ei2);

                        bool inconsistent_ndof = false;
                        inconsistent_ndof |= fel1_test.GetNDof() != dnums1_test.Size();
                        inconsistent_ndof |= fel2_test.GetNDof() != dnums2_test.Size();
                        inconsistent_ndof |= fel1_trial.GetNDof() != dnums1_trial.Size();
                        inconsistent_ndof |= fel2_trial.GetNDof() != dnums2_trial.Size();
                        
                        if(inconsistent_ndof)
                        {
                          cout << "facet, neighbouring fel(1): Trial NDof() = " << fel1_trial.GetNDof() << endl;
                          cout << "facet, neighbouring fel(1): Test NDof() = " << fel1_test.GetNDof() << endl;
                          cout << "facet, neighbouring fel(2): Trial NDof() = " << fel2_trial.GetNDof() << endl;
                          cout << "facet, neighbouring fel(2): Test NDof() = " << fel2_test.GetNDof() << endl;
                          cout << "facet, neighbouring fel(1): dnums1_trial.Size() = " << dnums1_trial.Size() << endl;
                          cout << "facet, neighbouring fel(1): dnums1_test.Size() = " << dnums1_test.Size() << endl;
                          cout << "facet, neighbouring fel(2): dnums2_trial.Size() = " << dnums2_trial.Size() << endl;
                          cout << "facet, neighbouring fel(2): dnums2_test.Size() = " << dnums2_test.Size() << endl;
                          throw Exception ( "Inconsistent number of degrees of freedom " );
                        }

                        for (auto & bfi : facetwise_skeleton_parts[VOL])
                        {
                          if (!bfi->DefinedOn (ma->GetElIndex (ei1))) continue;
                          if (!bfi->DefinedOn (ma->GetElIndex (ei2))) continue;
                          if (!bfi->DefinedOnElement(i)) continue;                        
                              
                          // if (check_unused)
                          //   for (auto d : dnums)
                          //     if (IsRegularDof(d)) useddof[d] = true;
                              
                          int elmat_width = (dnums1_trial.Size() + dnums2_trial.Size()) * fespace->GetDimension();
                          int elmat_height = (dnums1_test.Size() + dnums2_test.Size()) * fespace->GetDimension();
                          FlatMatrix<SCAL> elmat(elmat_height, elmat_width, lh);
                              
                          {
                            auto & mapped_trafo1 = eltrans1.AddDeformation(bfi->GetDeformation().get(), lh);
                            auto & mapped_trafo2 = eltrans2.AddDeformation(bfi->GetDeformation().get(), lh);
                            bfi->CalcFacetMatrix (fel1,facnr1,mapped_trafo1,vnums1,
                                                  fel2,facnr2,mapped_trafo2,vnums2, elmat, lh);
                          }
                              
                          fespace->TransformMat (ei1, elmat.Rows(0,dnums1_test.Size()*fespace->GetDimension()), TRANSFORM_MAT_LEFT);
                          fespace->TransformMat (ei2, elmat.Rows(dnums1_test.Size()*fespace->GetDimension(),elmat_height), TRANSFORM_MAT_LEFT);
                          fespace->TransformMat (ei1, elmat.Cols(0,dnums1_trial.Size()*fespace->GetDimension()), TRANSFORM_MAT_RIGHT);
                          fespace->TransformMat (ei2, elmat.Cols(dnums1_trial.Size()*fespace->GetDimension(),elmat_width), TRANSFORM_MAT_RIGHT);
                              
                          if (printelmat)
                          {
                            testout->precision(8);
                                  
                            (*testout) << "facet-elnum= " << i << endl;
                            (*testout) << "integrator " << bfi->Name() << endl;
                            (*testout) << "dnums1_trial = " << endl << dnums1_trial << endl;
                            (*testout) << "dnums1_test = " << endl << dnums1_test << endl;
                            (*testout) << "dnums2_trial = " << endl << dnums2_trial << endl;
                            (*testout) << "dnums2_test = " << endl << dnums2_test << endl;
                            (*testout) << "element1-index = " << eltrans1.GetElementIndex() << endl;
                            (*testout) << "element2-index = " << eltrans2.GetElementIndex() << endl;
                            (*testout) << "elmat = " << endl << elmat << endl;
                          }
                              
                          // ArrayMem<int, 50> dnums;
                          // dnums.SetSize(0);
                          // dnums.Append(dnums1);
                          // dnums.Append(dnums2);
                              
                          ArrayMem<int, 50> map_trial(dnums_trial.Size());
                          for (int i = 0; i < map_trial.Size(); i++) map_trial[i] = i;
                          QuickSortI (dnums_trial, map_trial);
                              
                          ArrayMem<int,50> compressed_dnums_trial;
                          compressed_dnums_trial.SetSize(0);
                              
                          ArrayMem<int,50> dnums_to_compressed_trial(dnums_trial.Size());
                          int compressed_dofs_trial = 0;
                          for (int i = 0; i < dnums_trial.Size(); ++i)
                          {
                            if (i==0 || (dnums_trial[map_trial[i]] != dnums_trial[map_trial[i-1]]))
                            {
                              compressed_dnums_trial.Append(dnums_trial[map_trial[i]]);
                              dnums_to_compressed_trial[map_trial[i]] = compressed_dofs_trial++;
                            }
                            else
                            {
                              dnums_to_compressed_trial[map_trial[i]] = dnums_to_compressed_trial[map_trial[i-1]];
                            }
                          }

                          ArrayMem<int, 50> map_test(dnums_test.Size());
                          for (int i = 0; i < map_test.Size(); i++) map_test[i] = i;
                          QuickSortI (dnums_test, map_test);
                              
                          ArrayMem<int,50> compressed_dnums_test;
                          compressed_dnums_test.SetSize(0);
                              
                          ArrayMem<int,50> dnums_to_compressed_test(dnums_test.Size());
                          int compressed_dofs_test = 0;
                          for (int i = 0; i < dnums_test.Size(); ++i)
                          {
                            if (i==0 || (dnums_test[map_test[i]] != dnums_test[map_test[i-1]]))
                            {
                              compressed_dnums_test.Append(dnums_test[map_test[i]]);
                              dnums_to_compressed_test[map_test[i]] = compressed_dofs_test++;
                            }
                            else
                            {
                              dnums_to_compressed_test[map_test[i]] = dnums_to_compressed_test[map_test[i-1]];
                            }
                          }
                              
                          FlatMatrix<SCAL> compressed_elmat(compressed_dofs_test * fespace->GetDimension(),
                                                            compressed_dofs_trial * fespace->GetDimension(), lh);
                          compressed_elmat = 0.0;
                          for (int i = 0; i < dnums_test.Size(); ++i)
                            for (int j = 0; j < dnums_trial.Size(); ++j)
                              compressed_elmat(dnums_to_compressed_test[i],dnums_to_compressed_trial[j]) += elmat(i,j);
                              
                          if (elmat_ev)
                          {
                            testout->precision(8);
                                  
                            (*testout) << "elind1 = " << eltrans1.GetElementIndex() << endl;
                            (*testout) << "elind2 = " << eltrans2.GetElementIndex() << endl;
#ifdef LAPACK
                            LapackEigenSystem(compressed_elmat, lh);
#else
                            Vector<SCAL> lami(compressed_elmat.Height());
                            Matrix<SCAL> evecs(compressed_elmat.Height());
                                  
                            CalcEigenSystem (compressed_elmat, lami, evecs);
                            (*testout) << "lami = " << endl << lami << endl;
#endif
                            // << "evecs = " << endl << evecs << endl;
                          }
                              
                          //                    for(int k=0; k<elmat.Height(); k++)
                          //                      if(fabs(elmat(k,k)) < 1e-7 && dnums[k] != -1)
                          //                        cout << "dnums " << dnums << " elmat " << elmat << endl; 

                          // cout << "compressed_elmat:" << compressed_elmat << endl;
                          // cout << "compressed_dnums_test:" << compressed_dnums_test << endl;
                          // cout << "compressed_dnums_trial:" << compressed_dnums_trial << endl;
                          
                          AddElementMatrix (compressed_dnums_test, compressed_dnums_trial, compressed_elmat, ElementId(BND,i), false, lh);
                        }
                      }
                    });
              }
              progress.Done();
            }
            
            
            //if (facetwise_skeleton_parts[BND].Size()) throw Exception ("mixed biforms don't support skeleton terms");
            if (facetwise_skeleton_parts[BND].Size())
              {
                // cout << "check bnd" << endl;
                // int cnt = 0;
                int ne = ma->GetNE(BND);
                ParallelForRange
                  ( IntRange(ne), [&] ( IntRange r )
                    {
                      LocalHeap lh = clh.Split();
                      Array<int> fnums, elnums, vnums, svnums, dnums_trial, dnums_test;

                      for (int i : r)
                        {
                          //{
                            //lock_guard<mutex> guard(printmatasstatus2_mutex);
                            //cnt++;
                            //gcnt++;
                            //// if (cnt % 10 == 0)
                              //// cout << "\rassemble facet surface element " << cnt << "/" << ne << flush;
                            //ma->SetThreadPercentage ( 100.0*(gcnt) / (loopsteps) );
                          //}

                          HeapReset hr(lh);
                          ElementId sei(BND, i);

                          // if (!fespace->DefinedOn (BND,ma->GetElIndex (sei))) continue;
                          fnums = ma->GetElFacets(sei);
                          int fac = fnums[0];
                          ma->GetFacetElements(fac,elnums);
                          int el = elnums[0];
                          ElementId ei(VOL, el);
                          fnums = ma->GetElFacets(ei);
                          //const FiniteElement & fel = fespace->GetFE (ei, lh);

                        const FiniteElement & fel_trial = fespace->GetFE (ei, lh);
                        const FiniteElement & fel_test = fespace2->GetFE (ei, lh);
                        MixedFiniteElement fel(fel_trial, fel_test);

                          int facnr = 0;
                          for (int k=0; k<fnums.Size(); k++)
                            if(fac==fnums[k]) facnr = k;
                          vnums = ma->GetElVertices (ei);
                          svnums = ma->GetElVertices (sei);

                          ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
                          ElementTransformation & seltrans = ma->GetTrafo (sei, lh);

                          //fespace->GetDofNrs (ei, dnums);

                        fespace->GetDofNrs (ei, dnums_trial);
                        fespace2->GetDofNrs (ei, dnums_test);

                        bool inconsistent_ndof = false;
                        inconsistent_ndof |= fel_test.GetNDof() != dnums_test.Size();
                        inconsistent_ndof |= fel_trial.GetNDof() != dnums_trial.Size();


                          //if(fel.GetNDof() != dnums.Size())
                          if(inconsistent_ndof)
                            {
                              (*testout) << "Surface fel::GetNDof() = " << fel.GetNDof() << endl;
                              (*testout) << "dnums_test.Size() = " << dnums_test.Size() << endl;
                              (*testout) << "dnums test = " << endl << dnums_test << endl;
                              (*testout) << "dnums_trial.Size() = " << dnums_trial.Size() << endl;
                              (*testout) << "dnums trial = " << endl << dnums_trial << endl;
                              throw Exception ( string("Inconsistent number of degrees of freedom Surface!") );
                            }

                          /*
                            for (int j = 0; j < NumIntegrators(); j++)
                            {
                            const BilinearFormIntegrator & bfi = *parts[j];
                          */
                          for (auto & bfi : facetwise_skeleton_parts[BND])
                            {
                              // if (bfi.VB() != BND) continue;
                              // if (!bfi.SkeletonForm()) continue;

                              if (!bfi->DefinedOn (ma->GetElIndex(sei) )) continue;
                              if (!bfi->DefinedOnElement (i)) continue;

                              //if (check_unused)
                                //for (int k = 0; k < dnums.Size(); k++)
                                  //if (IsRegularDof(dnums[k]))
                                    //useddof[dnums[k]] = true;

                              //int elmat_size = dnums.Size()*fespace->GetDimension();
                              //FlatMatrix<SCAL> elmat(elmat_size, lh);
                          FlatMatrix<SCAL> elmat(dnums_test.Size()*fespace->GetDimension(),
                                                 dnums_trial.Size()*fespace->GetDimension(), lh);

                              // original version did not compile on MacOS V
                              const FacetBilinearFormIntegrator & fbfi =
                                dynamic_cast<const FacetBilinearFormIntegrator&>(*bfi);
                              auto & mapped_trafo = eltrans.AddDeformation(fbfi.GetDeformation().get(), lh);
                              auto & mapped_strafo = seltrans.AddDeformation(fbfi.GetDeformation().get(), lh);
                              fbfi.CalcFacetMatrix (fel,facnr,mapped_trafo,vnums, mapped_strafo, svnums, elmat, lh);

                              fespace->TransformMat (ei, elmat, TRANSFORM_MAT_LEFT_RIGHT);

                              if (printelmat)
                                {
                                  testout->precision(8);

                                  (*testout) << "surface-elnum= " << i << endl;
                                  (*testout) << "integrator " << bfi->Name() << endl;
                                  (*testout) << "dnums test = " << endl << dnums_test << endl;
                                  (*testout) << "dnums trial = " << endl << dnums_trial << endl;
                                  (*testout) << "element-index = " << eltrans.GetElementIndex() << endl;
                                  (*testout) << "elmat = " << endl << elmat << endl;
                                }


                              if (elmat_ev)
                                {
                                  testout->precision(8);

                                  (*testout) << "elind = " << eltrans.GetElementIndex() << endl;
#ifdef LAPACK
                                  LapackEigenSystem(elmat, lh);
#else
                                  Vector<SCAL> lami(elmat.Height());
                                  Matrix<SCAL> evecs(elmat.Height());

                                  CalcEigenSystem (elmat, lami, evecs);
                                  (*testout) << "lami = " << endl << lami << endl;
#endif
                                  // << "evecs = " << endl << evecs << endl;
                                }

                              //                    for(int k=0; k<elmat.Height(); k++)
                              //                      if(fabs(elmat(k,k)) < 1e-7 && dnums[k] != -1)
                              //                        cout << "dnums " << dnums << " elmat " << elmat << endl;
                              {
                                // lock_guard<mutex> guard(addelemfacbnd_mutex);
                                AddElementMatrix (dnums_test, dnums_trial, elmat, ElementId(BND,i), true, lh);
                              }
                            }//end for (numintegrators)
                        }//end for nse
                    });//end of parallel
                // cout << "\rassemble facet surface element " << ne << "/" << ne << endl;
              } // if facetwise_skeleton_parts[BND].size

            if (elementwise_skeleton_parts.Size()) throw Exception ("mixed biforms don't support elementwise skeleton terms");
            
            if (print) *testout << "mat = " << mat << endl;

            
            
            //cout << endl;
          }

        //  WriteMatrix (*testout);

        if (checksum)
          cout << "|matrix| = " 
               << setprecision(16) << L2Norm (GetMatrix().AsVector()) << endl;
      }
    catch (Exception & e)
      {
        cout << "catch in AssembleBilinearform 2: " << e.What() << endl;
        e.Append (string ("in Assemble BilinearForm '") + 
                  GetName() + string("'\n"));
        throw;
      }
    
    catch (exception & e)
      {
        throw (Exception (string(e.what()) +
                          string("\n in Assemble BilinearForm '" + GetName() + "'\n")));
      }
    
  }

  
  template <class SCAL>
  void S_BilinearForm<SCAL> :: 
  ModifyRHS (BaseVector & f) const
  {
    if (keep_internal)
      f += *GetHarmonicExtensionTrans() * f;
  }

  template <class SCAL>
  void S_BilinearForm<SCAL> :: 
  ComputeInternal (BaseVector & u, const BaseVector & f, LocalHeap & clh) const
  {
    if (!eliminate_internal) return;

    static Timer timer ("Compute Internal");
    RegionTimer reg (timer);


    try
      {
        ma->PushStatus ("Compute Internal");

        int ne = ma->GetNE();

        if (VB_parts[VOL].Size())
          {
            if (keep_internal)
              {
                cout << IM(1) << "compute internal element ... ";
                
                //Set u_inner to zero
                for (int i = 0; i < ne; i++)
                  {
                    HeapReset hr(clh);
                    Array<int> dnums;
                    fespace->GetDofNrs (ElementId(VOL,i), dnums, LOCAL_DOF);            
                    FlatVector<SCAL> elu (dnums.Size(), clh);
                    elu = 0.0;
                    u.SetIndirect (dnums, elu);
                  }
                
                if (linearform)
                  u += *GetInnerSolve() * linearform->GetVector();
                else
                  u += *GetInnerSolve() * f;
                
                u += *GetHarmonicExtension() * u;
                cout << IM(1) << endl;
              }
            else
              {  
                ProgressOutput progress (ma, "compute internal element", ma->GetNE());

                IterateElements 
                  (*fespace, VOL, clh, 
                   [&] (ElementId ei, LocalHeap & lh)
                   {
                     progress.Update ();

                     const FiniteElement & fel = fespace->GetFE (ei, lh);
                     ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
                     // FlatArray<int> dnums = fespace->GetDofNrs (ei, lh);
                     Array<int> dnums (fel.GetNDof(), lh);
                     fespace->GetDofNrs (ei, dnums);

                     Array<int> idofs(dnums.Size(), lh);
                     fespace->GetDofNrs (ei, idofs, LOCAL_DOF);
                     if (!idofs.Size()) return;
                     for (int j = 0; j < idofs.Size(); j++)
                       idofs[j] = dnums.Pos(idofs[j]);
                     
                     int elmat_size = dnums.Size()*fespace->GetDimension();
                     FlatMatrix<SCAL> sum_elmat(elmat_size, lh);
                     sum_elmat = 0;
                     for (int j = 0; j < NumIntegrators(); j++)
                       {
                         HeapReset hr(lh);
                         const BilinearFormIntegrator & bfi = *parts[j];
                         
                         if (!bfi.VolumeForm()) continue;
                         if (!bfi.DefinedOn (eltrans.GetElementIndex())) continue;
                         if (!bfi.DefinedOnElement (ei.Nr())) continue;
                         
                         FlatMatrix<SCAL> elmat(elmat_size, lh);
                         bfi.CalcElementMatrix (fel, eltrans, elmat, lh);
                         
                         sum_elmat += elmat;
                       }
                
                     fespace->TransformMat (ei, sum_elmat, TRANSFORM_MAT_LEFT_RIGHT);
                     
                     if (idofs.Size())
                       {
                         int dim = sum_elmat.Height() / dnums.Size();
                         
                         FlatVector<SCAL> elf (sum_elmat.Height(), lh);
                         FlatVector<SCAL> elu (sum_elmat.Height(), lh);
                         FlatVector<SCAL> resi(dim*idofs.Size(), lh);
                         FlatVector<SCAL> wi(dim*idofs.Size(), lh);
                         
                         FlatMatrix<SCAL> ai(dim*idofs.Size(), lh);
                         
                         if (linearform)
                           linearform->GetVector().GetIndirect (dnums, elf);
                         else
                           f.GetIndirect (dnums, elf);
                         u.GetIndirect (dnums, elu);
                          
                         // compute residuum
                         elf -= sum_elmat * elu;
                          
                         for (int jj = 0; jj < idofs.Size(); jj++)
                           for (int j = 0; j < dim; j++)
                             {
                               for (int kk = 0; kk < idofs.Size(); kk++)
                                 for (int k = 0; k < dim; k++)
                                   {
                                     ai(jj*dim+j, kk*dim+k) =
                                       sum_elmat(idofs[jj]*dim+j, idofs[kk]*dim+k);
                                   }
                               resi(jj*dim+j) = elf(idofs[jj]*dim+j);
                             }
                         
                         // *testout << "residuum = " << resi << endl;
                         
                         LapackInverse (ai);
                         wi = ai * resi;
                         
                         // *testout << "inv_ai = " << endl << inv_ai << endl;
                         
                         for (int jj = 0; jj < idofs.Size(); jj++)
                           for (int j = 0; j < dim; j++)
                             elu(idofs[jj]*dim+j) += wi(jj*dim+j);
                         
                         //  *testout << "elu = " << endl << elu << endl;
                         
                         u.SetIndirect (dnums, elu);
                       }
                   });
                
                progress.Done();
                
              }//end of keep_internal-if
          }
        ma->PopStatus ();
      }

    catch (Exception & e)
      {
        e.Append (string("in ComputeInternal\n"));
        throw;
      }
    catch (exception & e)
      {
        throw (Exception (string(e.what()) +
                          string("\n in ComputeInternal\n")));
      }
  }










  template <class SCAL>
  void S_BilinearForm<SCAL> :: AssembleLinearization (const BaseVector & lin,
                                                      LocalHeap & clh, 
                                                      bool reallocate)
  {
    static Timer timer ("Assemble Linearization");
    static Timer timervol ("Assemble Linearization - volume");
    static Timer timerbound ("Assemble Linearization - boundary");
    static Timer timerbbound ("Assemble Linearization - co dim 2");
    static Timer timerspecial ("Assemble Linearization - SpecialElements");
    // static Timer timerVB[] = { timervol, timerbound, timerbbound };

    static mutex addelmatboundary1_mutex;

    lin.Cumulate();

    if (nonassemble) {
      shared_ptr<BaseMatrix> app =
        make_shared<LinearizedBilinearFormApplication> (dynamic_pointer_cast<BilinearForm>(this->shared_from_this()), &lin, clh);
      if (fespace->IsParallel())
        app = make_shared<ParallelMatrix> (app, GetTrialSpace()->GetParallelDofs(), GetTestSpace()->GetParallelDofs(), C2D);
      if (this->mats.Size() < this->ma->GetNLevels())
	mats.Append (app);
      else // we might be using a different vector now
	mats.Last() = app;
      return;
    }

    RegionTimer reg (timer);
    ma->PushStatus ("Assemble Linearization");

    if (specialelements_timestamp > graph_timestamp)
      {
        reallocate = true;
        cout << IM(3) << "reallocate due to changed special elements" << endl;
      }

    if(reallocate && this->mats.Size())
      this->mats.DeleteLast();

    if (this->mats.Size() < this->ma->GetNLevels())
      {
        AllocateMatrix();
        graph_timestamp = GetNextTimeStamp();
      }

    // timestamp = ++global_timestamp;
    timestamp = GetNextTimeStamp();

    try
      {
        if(MixedSpaces())
          {
            cout << IM(3) << "assemble linearization mixed bilinearform" << endl;
            auto& mat = GetMatrix();
            mat = 0.0;
            for(auto vb : {VOL,BND,BBND,BBBND})
              if(VB_parts[vb].Size())
                IterateElements
                  (*fespace, vb, clh,          // coloring for 1 space is enough
                   [&] (FESpace::Element ei, LocalHeap & lh)
                   {
                     const FiniteElement & fel1 = fespace->GetFE (ei, lh);
                     const FiniteElement & fel2 = fespace2->GetFE (ei, lh);

                     Array<int> dnums1(fel1.GetNDof(), lh);
                     Array<int> dnums2(fel2.GetNDof(), lh);
                     const ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
                     fespace->GetDofNrs (ei, dnums1);
                     fespace2->GetDofNrs (ei, dnums2);
                     FlatVector<SCAL> elveclin(dnums1.Size() * fespace->GetDimension(),lh);
                     lin.GetIndirect(dnums1,elveclin);
                     fespace->TransformVec(ei,elveclin,TRANSFORM_SOL);
                     FlatMatrix<SCAL> elmat(dnums2.Size(), dnums1.Size(), lh);
                     for (auto & bfi : VB_parts[vb])
                       {
                         if (!bfi->DefinedOn(fespace->GetMeshAccess()->GetElIndex(ei))) continue;
                         if (!bfi->DefinedOnElement(ei.Nr())) continue;
                         MixedFiniteElement fel(fel1, fel2);
                         bfi->CalcLinearizedElementMatrix (fel, eltrans, elveclin, elmat, lh);
                         /*
                           fespace->Transform (i, true, elmat, TRANSFORM_MAT_RIGHT);
                           fespace2->Transform (i, true, elmat, TRANSFORM_MAT_LEFT);
                         */
                         AddElementMatrix (dnums2, dnums1, elmat, ei, false, lh);
                       }
                   });
            return;
          }
        int ndof = fespace->GetNDof();
        Array<bool> useddof(ndof);
        useddof = false;
      
        BaseMatrix & mat = GetMatrix();
        mat = 0.0;
      
        cout << IM(3) << "Assemble linearization" << endl;
      
        Array<int> dnums;

        for (auto pre : preconditioners)
          pre -> InitLevel(fespace->GetFreeDofs(eliminate_internal));
	
        for (VorB vb : { VOL, BND, BBND, BBBND })
          if (VB_parts[vb].Size() || ((vb == VOL) && facetwise_skeleton_parts[BND].Size())) 
          {
            RegionTimer reg(timervol);
            ProgressOutput progress(ma,string("assemble ") + ToString(vb) + string(" element"), ma->GetNE(vb));

            /*
            if ( (vb == VOL || (!VB_parts[VOL].Size() && vb==BND) ) && eliminate_internal && keep_internal)
              {
                size_t ndof = fespace->GetNDof();
                size_t ne = ma->GetNE(vb);
                harmonicext = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
                if (!symmetric)
                  harmonicexttrans = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
                else
                  harmonicexttrans = make_shared<Transpose>(*harmonicext);
                innersolve = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
                if (store_inner)
                  innermatrix = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
              }
            */
            
            IterateElements 
              (*fespace, vb, clh,  [&] (FESpace::Element el, LocalHeap & lh)
               {
                 bool elem_has_integrator = false;
                 for (auto & bfi : VB_parts[vb])
                   if ((bfi->DefinedOn (el.GetIndex()))&&(bfi->DefinedOnElement (el.Nr())))
                     elem_has_integrator = true;

		 if ((vb == VOL) && facetwise_skeleton_parts[BND].Size())
		   elem_has_integrator = true;
				 
                 if (!elem_has_integrator)
                   return;
                 
                 const FiniteElement & fel = fespace->GetFE (el, lh);
                 ElementTransformation & eltrans = ma->GetTrafo (el, lh);
                 
                 Array<int> dnums(fel.GetNDof(), lh);
                 fespace->GetDofNrs (el, dnums);
                 
                 if(fel.GetNDof() != dnums.Size())
                   {
                     (*testout) << "fel::GetNDof() = " << fel.GetNDof() << endl;
                     (*testout) << "dnums.Size() = " << dnums.Size() << endl;
                     (*testout) << "dnums = " << dnums << endl;
                     throw Exception ( string("Inconsistent number of degrees of freedom fel::GetNDof() = ") + ToString(fel.GetNDof()) + string(" != dnums.Size() = ") + ToString(dnums.Size()) + string("!") );
                   }
                 for (auto d : dnums) 
                   if (IsRegularDof(d)) useddof[d] = true;
                 
                 FlatMatrix<SCAL> sum_elmat(dnums.Size()*fespace->GetDimension(), lh);
                 FlatMatrix<SCAL> elmat(dnums.Size()*fespace->GetDimension(), lh);
                 sum_elmat = 0;
                 
                 FlatVector<SCAL> elveclin (dnums.Size()*fespace->GetDimension(), lh);
                 lin.GetIndirect (dnums, elveclin);
                 fespace->TransformVec (el, elveclin, TRANSFORM_SOL);

                 for (auto & bfi : VB_parts[vb])
                   {
                     HeapReset hr(lh);
                     if (!bfi->DefinedOn (el.GetIndex())) continue;
                     if (!bfi->DefinedOnElement (el.Nr())) continue;
                     
                     try
                       {
                         auto & mapped_trafo = eltrans.AddDeformation(bfi->GetDeformation().get(), lh);
                         bfi->CalcLinearizedElementMatrix (fel, mapped_trafo, elveclin, elmat, lh);
                         
                         if (printelmat) 
                           {
                             lock_guard<mutex> guard(addelmatboundary1_mutex);
                             testout->precision(8);
                             if (vb != VOL)
                               {
                                 (*testout) << "surface-elnum = " << el.Nr() << endl;
                                 (*testout) << "boundary = " << ma->GetMaterial (el) << endl;
                               }
                             else
                               (*testout) << "elnum = " << el.Nr() << endl;
                             (*testout) << "eltype " << fel.ElementType() << endl;
                             (*testout) << "integrator " << bfi->Name() << endl;
                             (*testout) << "dnums = " << endl << dnums << endl;
                             (*testout) << "elveclin = " << endl << elveclin << endl;
                             (*testout) << "elmat = " << endl << elmat << endl;
                           }

                         if (elmat_ev)
                            {
                              testout->precision(8);
                              
                              (*testout) << "elind = " << eltrans.GetElementIndex() << endl;
#ifdef LAPACK
                              LapackEigenSystem(elmat, lh);
#else
                              Vector<SCAL> lami(elmat.Height());
                              Matrix<SCAL> evecs(elmat.Height());
                              
                              CalcEigenSystem (elmat, lami, evecs);
                              (*testout) << "lami = " << endl << lami << endl;
#endif
                              // << "evecs = " << endl << evecs << endl;
                            } 
                       }
                     catch (Exception & e)
                       {
                         e.Append (string("in Assemble Element Mat, bfi = ") + 
                                   bfi->Name() + string("\n"));
                         throw;
                       }
                     catch (exception & e)
                       {
                         throw (Exception (string(e.what()) +
                                           string("in Assemble Element Mat, bfi = ") + 
                                           bfi->Name() + string("\n")));
                       }
                     
                     sum_elmat += elmat;
                   }	 
		 
		 if((vb == VOL) && facetwise_skeleton_parts[BND].Size()) // skeleton = True on BND with element dofs
		   {
		     Array<int> fnums;
		     fnums = ma->GetElFacets(el);
		     
		     for(int facnr = 0; facnr < fnums.Size(); facnr++)//for (auto f : ma->GetElFacets(el))
		       {
			 Array<int> selnr;
			 ma->GetFacetSurfaceElements(fnums[facnr], selnr);
			 for(auto se : selnr) // se is sel nr, facnr is facet number index
			   {
			     HeapReset hr(lh);
			     ElementId sei(BND, se);
			     if (!fespace->DefinedOn (BND,ma->GetElIndex (sei))) continue;

			     auto ei = el;
			     Array<int>  vnums, svnums;
			     
			     vnums = ma->GetElVertices (ei);
			     svnums = ma->GetElVertices (sei);
			     
			     ElementTransformation & seltrans = ma->GetTrafo (sei, lh);
			     
			     if(fel.GetNDof() != dnums.Size())
			       {
				 (*testout) << "Surface fel::GetNDof() = " << fel.GetNDof() << endl;
				 (*testout) << "dnums.Size() = " << dnums.Size() << endl;
				 (*testout) << "dnums = " << dnums << endl;
				 throw Exception ( string("Inconsistent number of degrees of freedom Surface fel::GetNDof() = ") + ToString(fel.GetNDof()) + string(" != dnums.Size() = ") + ToString(dnums.Size()) + string("!") );
			       }

			     for (auto & bfi : facetwise_skeleton_parts[BND])
			       {
				 
				 if (!bfi->DefinedOn (ma->GetElIndex(sei) )) continue;                
                          
				 for (int k = 0; k < dnums.Size(); k++)
				   if (IsRegularDof(dnums[k]))
				     useddof[dnums[k]] = true;
                          
				 // int elmat_size = dnums.Size()*fespace->GetDimension();
				 
				 elmat = 0.0;
				 // original version did not compile on MacOS V
				 const FacetBilinearFormIntegrator & fbfi = 
				   dynamic_cast<const FacetBilinearFormIntegrator&>(*bfi);				 
				 
         auto & mapped_trafo = eltrans.AddDeformation(fbfi.GetDeformation().get(), lh);
         auto & mapped_strafo = seltrans.AddDeformation(fbfi.GetDeformation().get(), lh);
				 fbfi.CalcLinearizedFacetMatrix (fel,facnr,mapped_trafo,vnums, mapped_strafo, svnums, elveclin, elmat, lh);
				 fespace->TransformMat (ei, elmat, TRANSFORM_MAT_LEFT_RIGHT);
				 if (printelmat)
				   {
				     testout->precision(8);
                              
				     (*testout) << "surface-elnum= " << se << endl;
				     (*testout) << "integrator " << bfi->Name() << endl;
				     (*testout) << "dnums = " << endl << dnums << endl;
				     (*testout) << "element-index = " << eltrans.GetElementIndex() << endl;
				     (*testout) << "elmat = " << endl << elmat << endl;
				   }
                          
                          
				 if (elmat_ev)
				   {
				     testout->precision(8);
                              
				     (*testout) << "elind = " << eltrans.GetElementIndex() << endl;
#ifdef LAPACK
				     LapackEigenSystem(elmat, lh);
#else
				     Vector<SCAL> lami(elmat.Height());
				     Matrix<SCAL> evecs(elmat.Height());
                              
				     CalcEigenSystem (elmat, lami, evecs);
				     (*testout) << "lami = " << endl << lami << endl;
#endif				     
				   }
				 
				  sum_elmat += elmat;				  
			       } // for bf integrators			 
			   } //for elements in getsurfacetelement of facet		       
		       } // for facets of element
		   } //if VOL and skeleton bfi
                 
                 
                 fespace->TransformMat (el, sum_elmat, TRANSFORM_MAT_LEFT_RIGHT);
                 
                 if (printelmat)
                   *testout << "summat = " << sum_elmat << endl;
                 

                 // Array<int> lhdofs(dnums.Size(), lh);
                 // fespace->GetElementDofsOfType (el, lhdofs, HIDDEN_DOF);

                 bool has_hidden = false;
                 if (eliminate_hidden)
                   {
                     for (auto d : dnums)
                       if (fespace->GetDofCouplingType(d) & HIDDEN_DOF)
                         has_hidden = true;
                   }
                 
                 bool elim_only_hidden = (!eliminate_internal) && eliminate_hidden && /* (lhdofs.Size() > 0) */ has_hidden;

                 if ((vb == VOL || (!VB_parts[VOL].Size() && vb==BND) ) && (elim_only_hidden || eliminate_internal))
                   {
                     static Timer statcondtimer("static condensation", NoTracing);
                     RegionTimer regstat (statcondtimer);
                     
                     // ArrayMem<int,100> idofs, idofs1, odofs;
                     int i = el.Nr();

                     // fespace->GetElementDofsOfType (el, idofs1, elim_only_hidden ? HIDDEN_DOF : CONDENSABLE_DOF);
                     Array<int> idofs1(dnums.Size(), lh), odofs1(dnums.Size(), lh);
                     idofs1.SetSize0(); odofs1.SetSize0();
                     
                     auto ctype = elim_only_hidden ? HIDDEN_DOF : CONDENSABLE_DOF;
                     for (auto i : Range(dnums))
                       {
                         auto ct = fespace->GetDofCouplingType(dnums[i]);
                         if (ct & ctype)
                           idofs1.AppendHaveMem(i);
                         else
                           if (ct != UNUSED_DOF)
                             odofs1.AppendHaveMem(i);
                       }
                     
                     
                          
                     if (printelmat) 
                       {
                         *testout << "eliminate internal" << endl;
                         *testout << "idofs1 = " << idofs1 << endl;
                       }
                     
                    
                     if (idofs1.Size())
                       {
                         HeapReset hr (lh);
                         
                         // int size = sum_elmat.Height();
                         int dim = fespace->GetDimension();                         
                         // int dim = size / dnums.Size();

                         /*
                         idofs.SetSize (0);
                         for (int j = 0; j < idofs1.Size(); j++)
                           for (int jj = 0; jj < dim; jj++)
                             idofs.Append (dim*idofs1[j]+jj);
                         */
                         int sizei = dim * idofs1.Size();
                         int sizeo = dim * odofs1.Size();
                         
                         FlatArray<int> idofs (sizei, lh);
                         FlatArray<int> odofs (sizeo, lh);
                         
                         for (int j = 0, k = 0; j < idofs1.Size(); j++)
                           for (int jj = 0; jj < dim; jj++)
                             idofs[k++] = dim*idofs1[j]+jj;
                         
                         for (int j = 0, k = 0; j < odofs1.Size(); j++)
                           for (int jj = 0; jj < dim; jj++)
                             odofs[k++] = dim*odofs1[j]+jj;
                         
                         /*
                         odofs.SetSize (0);
                         for (int j = 0; j < size; j++)
                           if (!idofs.Contains(j))
                             odofs.Append(j);
                         */
                         
                         if (printelmat)
                           {
                             (*testout) << "idofs = " << endl << idofs << endl;
                             (*testout) << "odofs = " << endl << odofs << endl;
                           }

                         FlatMatrix<SCAL> a(sizeo, sizeo, lh);
                         FlatMatrix<SCAL> b(sizeo, sizei, lh);
                         FlatMatrix<SCAL> c(sizeo, sizei, lh);
                         FlatMatrix<SCAL> d(sizei, sizei, lh);
                         
                         a = sum_elmat.Rows(odofs).Cols(odofs);
                         b = sum_elmat.Rows(odofs).Cols(idofs);
                         c = Trans(sum_elmat.Rows(idofs).Cols(odofs));
                         d = sum_elmat.Rows(idofs).Cols(idofs);
                         
                              
                         // A := A - B D^{-1} C^T
                         // new Versions, July 07
                         if (elim_only_hidden || !keep_internal)
                           {
                             AInvBt (d, b);
                             LapackMultAddABt (b, c, -1, a);
                           }
                         else
                           {
                             /*
                             ArrayMem<int,50> idnums1, idnums;
                             ArrayMem<int,50> hdnums1;
                             ArrayMem<int,50> ednums1, ednums;
                             
                             fespace->GetElementDofsOfType(el,idnums1,CONDENSABLE_DOF);
                             fespace->GetElementDofsOfType(el,ednums1,EXTERNAL_DOF);
                             fespace->GetElementDofsOfType(el,hdnums1,HIDDEN_DOF);
                             for (auto d : Range(idnums1.Size()))
                               idnums1[d] = dnums[idnums1[d]];
                             for (auto d : Range(ednums1.Size()))
                               ednums1[d] = dnums[ednums1[d]];
                             for (auto ldof : hdnums1)
                               idnums1[ldof] = -1;
                             
                             for (int j = 0; j < idnums1.Size(); j++)
                               idnums += dim*IntRange(idnums1[j], idnums1[j]+1);
                             for (int j = 0; j < ednums1.Size(); j++)
                               ednums += dim * IntRange(ednums1[j], ednums1[j]+1);
                             */

                             Array<DofId> idnums1(idofs1.Size(), lh), ednums1(odofs1.Size(), lh);
                             for (int i : Range(idofs1))
                               {
                                 DofId d = dnums[idofs1[i]];
                                 if (fespace->GetDofCouplingType(d) == HIDDEN_DOF)
                                   d = NO_DOF_NR; // or _CONDENSE; ???
                                 idnums1[i] = d;
                               }
                             
                             for (int i : Range(odofs1))
                               ednums1[i] = dnums[odofs1[i]];

                             Array<int> idnums(dim*idnums1.Size(), lh);
                             Array<int> ednums(dim*ednums1.Size(), lh);
                             idnums.SetSize0(); 
                             ednums.SetSize0();
                             
                             for (DofId idnum1 : idnums1)
                               if (IsRegularDof(idnum1))
                                 idnums += dim*IntRange(idnum1, idnum1+1);
                               else
                                 for (size_t k = 0; k < dim; k++)
                                   idnums.AppendHaveMem(idnum1);
                             
                             for (DofId ednum1 : ednums1)
                               ednums += dim * IntRange(ednum1, ednum1+1);
                             
                             
                             if (store_inner)
                               static_cast<ElementByElementMatrix<SCAL>*>(innermatrix.get())
                                      ->AddElementMatrix(i,idnums,idnums,d);
                             
                             // LapackInverse (d);
                             CalcInverse (d);
                             FlatMatrix<SCAL> he (sizei, sizeo, lh);
                             // he=0.0;
                             // LapackMultAddABt (d, c, -1, he);
                             he = -d * Trans(c);
                             harmonicext_ptr->AddElementMatrix(i,idnums,ednums,he);
                                  
                             if (!symmetric)
                               {
                                 FlatMatrix<SCAL> het (sizeo, sizei, lh);
                                 // het = 0.0;
                                 // LapackMultAddAB (b, d, -1, het);
                                 het = -b * d;
                                 harmonicexttrans_ptr->AddElementMatrix(i,ednums,idnums,het);
                               }
                             innersolve_ptr->AddElementMatrix(i,idnums,idnums,d);
                             
                             // LapackMultAddAB (b, he, 1.0, a);
                             a += b*he;
                           }                                 
                         
                         
                         if (printelmat)
                           *testout << "schur matrix = " << a << endl;
                         
                         sum_elmat.Rows(odofs).Cols(odofs) = a;
                         
                         for (int k = 0; k < idofs1.Size(); k++)
                           dnums[idofs1[k]] = -1;
                       }
                   }
                 

                 AddElementMatrix (dnums, dnums, sum_elmat, el, false, lh);

                 for (auto pre : preconditioners)
                   pre -> AddElementMatrix (dnums, sum_elmat, el, lh);
               });

            progress.Done();
          }
    

        if (facetwise_skeleton_parts[VOL].Size()){
          size_t dummy;
          Assemble_facetwise_skeleton_parts_VOL(useddof, dummy, clh, &lin);
        }
        // if (facetwise_skeleton_parts[VOL].Size())
        //   throw Exception ("CalcLinearization for facet-wise skeleton-VOL not implemented");
        if (elementwise_skeleton_parts.Size())
          throw Exception ("CalcLinearization for element-wise skeleton-VOL not implemented");
                
        if (specialelements.Size())
          cout << IM(3) << "special elements: " << specialelements.Size() << endl;

        timerspecial.Start();
        static mutex specel_mutex;

        // auto& secoloring = SpecialElementColoring();   // needs more testing 
        ParallelForRange(IntRange(specialelements.Size()), [&](IntRange r)
        {
          LocalHeap lh = clh.Split();
          Array<DofId> dnums;
          for(auto i : r)
            {
              HeapReset hr(lh);
              const SpecialElement & el = *specialelements[i];
              el.GetDofNrs(dnums);
              FlatVector<SCAL> elvec(dnums.Size()*fespace->GetDimension(), lh);
              lin.GetIndirect (dnums, elvec);
              FlatMatrix<SCAL> elmat(dnums.Size() * fespace->GetDimension(), lh);
              el.CalcLinearizedElementMatrix(elvec, elmat, lh);

              {
                // lock_guard<mutex> guard(specel_mutex);
                for (int j = 0; j < dnums.Size(); j++)
                  if (IsRegularDof(dnums[j]))
                    useddof[dnums[j]] = true;
                AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), true, lh);
                for(auto pre : preconditioners)
                  pre->AddElementMatrix(dnums, elmat, ElementId(BND, i), lh);
              }
            }
        });
        timerspecial.Stop();
      
      
        // add eps to avoid empty lines
	HeapReset hr(clh);
        FlatMatrix<SCAL> elmat (fespace->GetDimension(), clh);
        elmat = 0;
        dnums.SetSize(1);
      
        if (eps_regularization != 0)
          {
            for (int i = 0; i < elmat.Height(); i++)
              elmat(i, i) = eps_regularization;
            for (int i = 0; i < ndof; i++)
              {
                dnums[0] = i; 
                AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), false, clh);
              }
          }
        if (unuseddiag != 0)
          {
            for (int i = 0; i < elmat.Height(); i++)
              elmat(i, i) = unuseddiag;
            for (int i = 0; i < ndof; i++)
              if (!useddof[i])
                {
                  dnums[0] = i;
                  AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), false, clh);
                }
          }
      
        if (print)
          {
            (*testout) << "mat = " << endl << GetMatrix() << endl;
            if (harmonicext) (*testout) << "harm ext = " << endl << *harmonicext << endl;
            if (harmonicexttrans) (*testout) << "harm ext = " << endl << *harmonicexttrans << endl;
            if (innermatrix) (*testout) << "harm ext = " << endl << *innermatrix << endl;
          }
        /*
          if (mat.Height() < 100)
          mat.Print (cout);
        */
        int cntused = 0;
        for (int i = 0; i < useddof.Size(); i++)
          if (useddof[i])
            cntused++;
        cout << IM(5) << "used " << cntused
             << ", unused = " << useddof.Size()-cntused
             << ", total = " << useddof.Size() << endl;

        for (int j = 0; j < preconditioners.Size(); j++)
          preconditioners[j] -> FinalizeLevel(&GetMatrix());
      }
    catch (Exception & e)
      {
        stringstream ost;
        ost << "in AssembleLinearization\n" << endl;
        e.Append (ost.str());
        throw;
      }
    catch (exception & e)
      {
        throw (Exception (string(e.what()) +
                          string("\n in AssembleLinearization\n")));
      }
    ma->PopStatus();
  }


 template<class SCAL>
 void S_BilinearForm<SCAL> :: AddMatrixTP(SCAL val,
                                          const BaseVector & x,
                                          BaseVector & y, LocalHeap & clh) const
 {
    static Timer timerall ("Apply Matrix1 (TP) - all");
    static Timer timervol ("Apply Matrix1 (TP) - volume");
    static Timer timerfac1 ("Apply Matrix1 (TP) - facets 1");
    static Timer timerfac2 ("Apply Matrix1 (TP) - facets 2");
    RegionTimer rall(timerall);
    // bool hasbound = false;
    bool hasinner = false;
    bool hasskeletonbound = false;
    bool hasskeletoninner = false;
    int volumeintegrals = -1;
    for(int j=0;j<parts.Size();j++)
    {
      const BilinearFormIntegrator & bfi = *GetIntegrator(j);
      if (bfi.BoundaryForm())
        {
          if (bfi.SkeletonForm())
            hasskeletonbound = true;
          // else
          // hasbound = true;
        }
      else
        if (bfi.SkeletonForm())
          hasskeletoninner = true; 
        else
        {
          hasinner = true;
          volumeintegrals = j;
        }
    }
    LocalHeap chelperheap(1000000000,"my x heap");
    const shared_ptr<TPHighOrderFESpace> & tpfes = dynamic_pointer_cast<TPHighOrderFESpace > (fespace);
    const Array<shared_ptr<FESpace> > & spaces = tpfes->Spaces(0);
    int dimspace = tpfes->GetDimension();
    const Table<int> & element_coloring0 = spaces[0]->ElementColoring(VOL);
    auto meshx = spaces[0]->GetMeshAccess();
    auto meshy = spaces[1]->GetMeshAccess();
    int nelx = meshx->GetNE();
    int nely = meshy->GetNE();
    int ndofxspace = spaces[0]->GetNDof();
    int ndofyspace = spaces[1]->GetNDof();
    if(hasinner)
    {
      timervol.Start();
      // if (task_manager)
      // {
        //HeapReset hr(chelperheap);
        for (FlatArray<int> els_of_col : element_coloring0)
        {
          SharedLoop2 sl(els_of_col.Range());
          // task_manager -> CreateJob
          ParallelJob
          ( [&] (const TaskInfo & ti) 
          {
            LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
            LocalHeap xheap = chelperheap.Split(ti.thread_nr, ti.nthreads);
            for (int mynr : sl)
            {
              HeapReset hr(lh);
              HeapReset hrx(xheap);
              int elnrx = els_of_col[mynr];
              auto & felx = spaces[0]->GetFE(ElementId(elnrx),lh);
              int ndofx = felx.GetNDof();
              double ndofxinv = 1.0/ndofx;
              const ElementTransformation & xtrafo = meshx->GetTrafo(ElementId(elnrx), lh);
              const IntegrationRule & ir = SelectIntegrationRule(felx.ElementType(),2*felx.Order());
              BaseMappedIntegrationRule & mir = xtrafo(ir, lh);
              FlatMatrix<> elvec_yslicemat(ndofx,ndofyspace*dimspace,xheap);
              Array<int> dnums_yslice(ndofx*ndofyspace, xheap);
              tpfes->GetSliceDofNrs(ElementId(elnrx), 1, dnums_yslice,xheap);
              x.GetIndirect (dnums_yslice, elvec_yslicemat.AsVector());
              static_cast<TensorProductBilinearFormIntegrator &>(*parts[volumeintegrals]).ApplyXElementMatrix(felx, xtrafo, elvec_yslicemat, &xheap,&mir, lh);
              int firstydof = 0;
              for(int j=0;j<nely;j++)
              {
                HeapReset hr(lh);
                ElementId elid(j+elnrx*nely);
                auto & tpfel = tpfes->GetFE(elid,lh);
                // int ndofy = spaces[1]->GetFE(ElementId(j),lh).GetNDof();
                int ndofy = ndofxinv*tpfel.GetNDof();
                IntRange dnumsy(firstydof, firstydof+dimspace*ndofy);
                firstydof+=dimspace*ndofy;
                const ElementTransformation & tptrafo = tpfes->GetTrafo(elid,lh);
                static_cast<TensorProductBilinearFormIntegrator &>(*parts[volumeintegrals]).ApplyYElementMatrix(tpfel,tptrafo,dnumsy,xtrafo.userdata,&mir,lh);
              }
              FlatMatrix<> elvecy_mat(ndofx,ndofyspace*dimspace,lh);
              static_cast<TensorProductBilinearFormIntegrator &>(*parts[volumeintegrals]).ApplyXElementMatrixTrans(felx,xtrafo,elvecy_mat,xtrafo.userdata,&mir,lh);
              //elvecy_mat *= (double)val;
              y.AddIndirect(dnums_yslice, elvecy_mat.AsVector());
            }
          }
          );
        }
      // }
      timervol.Stop();
    }
    // bool needs_facet_loop = false;
    // bool needs_element_boundary_loop = false;
    // bool neighbor_testfunction = false;
    int facetvolumeintegrals = -1;
    int facetboundaryintegrals = -1;
    if (hasskeletonbound||hasskeletoninner)
    {
      for (int j = 0; j < NumIntegrators(); j++)
      {
        if (parts[j] -> SkeletonForm())
        {
          auto dgform = parts[j] -> GetDGFormulation();
          if (!dgform.element_boundary && !parts[j]->BoundaryForm())
          {
            // needs_facet_loop = true;
            facetvolumeintegrals = j;
          }
          if (!dgform.element_boundary && parts[j]->BoundaryForm())
          {
            // needs_facet_loop = true;
            facetboundaryintegrals = j;
          }
          if (dgform.element_boundary)
          {
            throw Exception("Element boundary formulation is not implemented for tensor product spaces, please reformulate as skeleton integrals");
            // needs_element_boundary_loop = true;
          }
        }
      }
      // do we need locks for neighbor - testfunctions ?
      /*
      for (int j = 0; j < NumIntegrators(); j++)
        if (parts[j] -> SkeletonForm())
        {
          auto dgform = parts[j] -> GetDGFormulation();
          // if (dgform.neighbor_testfunction)
          // neighbor_testfunction = true;
        }
      */
    }
    else
      return;
      
    if(facetvolumeintegrals == -1 && facetboundaryintegrals == -1)
      return;
    // auto & nels = tpfes->GetNels();
    // auto & nfacets = tpfes->GetNFacets();
    timerfac1.Start();
    for (FlatArray<int> colfacets : spaces[0]->FacetColoring())
    {
	  //HeapReset hr(chelperheap);
      SharedLoop2 sl(colfacets.Range());
      // task_manager -> CreateJob
      ParallelJob
      ( [&] (const TaskInfo & ti) 
      {
        LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
        LocalHeap xheap = chelperheap.Split(ti.thread_nr, ti.nthreads);
        for (int i : sl)
        {
          HeapReset hr(lh);
          HeapReset hrx(xheap);
          Array<int> elnums_x(2, lh), elnums_per_x(2,lh), fnums1_x(6, lh), fnums2_x(6, lh), vnums1(8, lh), vnums2(8, lh);
          int facet_x = colfacets[i];
          int facet2_x = colfacets[i];
          // Horizontal edge - get facet elements w.r.t. first direction
          meshx->GetFacetElements (facet_x, elnums_x);
          int el1_x = elnums_x[0];
          auto & felx1 = spaces[0]->GetFE(ElementId(el1_x),lh);
          int ndofx1 = felx1.GetNDof();
          double ndofxinv = 1.0/ndofx1;          
          // The element facets:           
          fnums1_x = meshx->GetElFacets(ElementId(VOL,el1_x));
          int facnr_x1 = fnums1_x.Pos(facet_x);
          vnums1 = meshx->GetElVertices (ElementId(VOL,el1_x));
          const ElementTransformation & eltransx1 = meshx->GetTrafo(el1_x,lh);
          auto eltype1 = eltransx1.GetElementType();
          auto etfacet = ElementTopology::GetFacetType (eltype1, facnr_x1);
          Facet2ElementTrafo transform1(eltype1, vnums1);
          if(elnums_x.Size() < 2)
          {
            facet2_x = meshx->GetPeriodicFacet(facet_x);
            if(facet2_x > facet_x)
            {
              meshx->GetFacetElements (facet2_x, elnums_per_x);
              elnums_x.Append(elnums_per_x[0]);
            }
            else if(facet2_x < facet_x)
              continue;
          }
          if(elnums_x.Size() < 2)
          {
            meshx->GetFacetSurfaceElements(facet_x, elnums_x);
            int sel = elnums_x[0];
            ElementId sei(BND,sel);
            vnums2 = meshx->GetElVertices (sei);
            for(int j=0;j<nely;j++)
            {
              HeapReset hr(lh);
              ElementId ei1(j+el1_x*nely);               
              ElementTransformation & eltrans = tpfes->GetTrafo (ei1, lh);
              ElementTransformation & seltrans = meshx->GetTrafo (sei, lh);
              const FiniteElement & fel = tpfes->GetFE (ei1, lh);
              Array<int> dnums(fel.GetNDof(), lh);
              tpfes->GetDofNrs (ei1, dnums);
              FlatVector<double> elx(dnums.Size()*dimspace, lh), ely(dnums.Size()*dimspace, lh);
              x.GetIndirect(dnums, elx);
              if( facetboundaryintegrals != -1)
              {
                static_cast<TensorProductFacetBilinearFormIntegrator &>(*parts[facetboundaryintegrals]).ApplyFacetMatrix(fel,facnr_x1,eltrans,vnums1, seltrans, vnums2, elx, ely, lh);
                y.AddIndirect(dnums, ely);
              }
            }
            continue;
          }
          else
          {
            // TP Element number of the second element sharing the facet
            int el2_x = elnums_x[1];
            auto & felx2 = spaces[0]->GetFE(ElementId(el2_x),lh);
            int ndofx2 = felx2.GetNDof();
            const ElementTransformation & eltransx2 = meshx->GetTrafo(el2_x,lh);
            fnums2_x = meshx->GetElFacets(ElementId(VOL,el2_x));
            // Local position of second facet
            int facnr_x2 = fnums2_x.Pos(facet2_x);
            // vnums stores the elements vertex numbers (needed for facet2element trafo)
            vnums2 = meshx->GetElVertices (ElementId(VOL,el2_x));             
            // Prepare Integration Rules:
            int maxorderx = max2 (felx1.Order(), felx2.Order());
            auto eltype2 = eltransx2.GetElementType();
            const IntegrationRule & ir_facet = SelectIntegrationRule(etfacet, 2*maxorderx);
            const IntegrationRule & ir_volx1 = transform1(facnr_x1, ir_facet, lh);
            Facet2ElementTrafo transform2(eltype2, vnums2);
            const IntegrationRule & ir_volx2 = transform2(facnr_x2, ir_facet, lh);
            BaseMappedIntegrationRule & mirx1 = eltransx1(ir_volx1, lh);
            BaseMappedIntegrationRule & mirx2 = eltransx2(ir_volx2, lh);            
            mirx1.ComputeNormalsAndMeasure (eltype1, facnr_x1);
            FlatMatrix<> elvec_yslicemat(ndofx1+ndofx2,ndofyspace*dimspace,lh);
            Array<int> dnums_yslice((ndofx1+ndofx2)*ndofyspace,lh),dnums_yslice1(ndofx2*ndofyspace,lh);
            tpfes->GetSliceDofNrs(ElementId(el1_x),1,dnums_yslice,xheap);
            tpfes->GetSliceDofNrs(ElementId(el2_x),1,dnums_yslice1,xheap);
            dnums_yslice.Append(dnums_yslice1);
            x.GetIndirect (dnums_yslice, elvec_yslicemat.AsVector());
            static_cast<TensorProductFacetBilinearFormIntegrator &>(*parts[facetvolumeintegrals]).ApplyXFacetMatrix(felx1, eltransx1, felx2,eltransx2, elvec_yslicemat, &xheap, &mirx1,&mirx2,lh);
            int firstydof = 0;
            for(int j=0;j<nely;j++)
            {
              HeapReset hr(lh);
              ElementId elid(j+el1_x*nely);
              auto & tpfel = tpfes->GetFE(elid,lh);
              // int ndofy = spaces[1]->GetFE(ElementId(j),lh).GetNDof();
              int ndofy = ndofxinv*tpfel.GetNDof();
              IntRange dnumsy(firstydof,firstydof+dimspace*ndofy);
              firstydof+=dimspace*ndofy;
              const ElementTransformation & tptrafo = tpfes->GetTrafo(elid,lh);
              static_cast<TensorProductFacetBilinearFormIntegrator &>(*parts[facetvolumeintegrals]).ApplyYElementMatrix(tpfel,tptrafo,dnumsy,eltransx1.userdata,&mirx1,lh);
            }
            FlatMatrix<> elmat(ndofx1+ndofx2,ndofyspace*dimspace,lh);
            elmat = 0.0;
            static_cast<TensorProductFacetBilinearFormIntegrator &>(*parts[facetvolumeintegrals]).ApplyXFacetMatrixTrans(felx1,eltransx1,felx2,eltransx2,elmat,eltransx1.userdata,&mirx1,&mirx2,lh);
            //elvec_mat *= val;
            y.AddIndirect (dnums_yslice, elmat.AsVector());
          }
        }
      });
    }
    timerfac1.Stop();
    timerfac2.Start();
    for (FlatArray<int> colfacets : spaces[1]->FacetColoring())
    {
      SharedLoop2 sl(colfacets.Range());
      // task_manager -> CreateJob
      ParallelJob
      ( [&] (const TaskInfo & ti) 
      {
        LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
        LocalHeap yheap = chelperheap.Split(ti.thread_nr, ti.nthreads);
        for (int i : sl)
        {         
          HeapReset hr(lh);
          HeapReset hry(yheap);
          Array<int> elnums_y(2, lh), elnums_per_y(2,lh), fnums1_y(6, lh), fnums2_y(6, lh), vnums1(8, lh), vnums2(8, lh);
          int facet_y = colfacets[i];
          int facet2_y = colfacets[i];
          // Horizontal edge - get facet elements w.r.t. second direction
          meshy->GetFacetElements (facet_y, elnums_y);
          int el1_y = elnums_y[0];
          auto & fely1 = spaces[1]->GetFE(ElementId(el1_y),lh);
          int ndofy1 = fely1.GetNDof();
          double ndofyinv = 1.0/ndofy1;
          // The element facets:           
          fnums1_y = meshy->GetElFacets(ElementId(VOL,el1_y));
          int facnr_y1 = fnums1_y.Pos(facet_y);
          vnums1 = meshy->GetElVertices (ElementId(VOL,el1_y));
          const ElementTransformation & eltransy1 = meshy->GetTrafo(el1_y,lh);
          auto eltype1 = eltransy1.GetElementType();
          auto etfacet = ElementTopology::GetFacetType (eltype1, facnr_y1);
          Facet2ElementTrafo transform1(eltype1, vnums1);
          if(elnums_y.Size() < 2)
          {
            facet2_y = meshy->GetPeriodicFacet(facet_y);
            if(facet2_y > facet_y)
            {
              meshy->GetFacetElements (facet2_y, elnums_per_y);
              elnums_y.Append(elnums_per_y[0]);
            }
            else if(facet2_y < facet_y)
              continue;
          }
          if(elnums_y.Size() < 2)
          {
            meshy->GetFacetSurfaceElements(facet_y, elnums_y);
            int sel = elnums_y[0];
            ElementId sei(BND,sel);
            vnums2 = meshy->GetElVertices (sei);
            for(int j=0;j<nelx;j++)
            {
              HeapReset hr(lh);
              ElementId ei1(j*nely+el1_y);
              ElementTransformation & eltrans = tpfes->GetTrafo (ei1, lh);
              ElementTransformation & seltrans = meshy->GetTrafo (sei, lh);
              const FiniteElement & fel = tpfes->GetFE (ei1, lh);
              Array<int> dnums(fel.GetNDof(), lh);
               tpfes->GetDofNrs (ei1, dnums);
              FlatVector<double> elx(dnums.Size()*dimspace, lh), ely(dnums.Size()*dimspace, lh);
              x.GetIndirect(dnums, elx);
              if ( facetboundaryintegrals != -1)
              {
                static_cast<TensorProductFacetBilinearFormIntegrator &>(*parts[facetboundaryintegrals]).ApplyFacetMatrix(fel,facnr_y1+10,eltrans,vnums1, seltrans, vnums2, elx, ely, lh);
                y.AddIndirect(dnums, ely);
              }
            }
             continue;
          }
          else
          {
            // TP Element number of the second element sharing the facet
            int el2_y = elnums_y[1];
            auto & fely2 = spaces[1]->GetFE(ElementId(el2_y),lh);
            int ndofy2 = fely2.GetNDof();
            const ElementTransformation & eltransy2 = meshy->GetTrafo(ElementId(el2_y),lh);
            fnums2_y = meshy->GetElFacets(ElementId(VOL,el2_y));
            // Local position of second facet
            int facnr_y2 = fnums2_y.Pos(facet2_y);
            // vnums stores the elements vertex numbers (needed for facet2element trafo)
            vnums2 = meshy->GetElVertices (ElementId(VOL,el2_y));
            // Prepare Integration Rules:
            int maxordery = max2 (fely1.Order(), fely2.Order());
            auto eltype2 = eltransy2.GetElementType();
            const IntegrationRule & ir_facet = SelectIntegrationRule(etfacet, 2*maxordery);
            const IntegrationRule & ir_voly1 = transform1(facnr_y1, ir_facet, lh);
            Facet2ElementTrafo transform2(eltype2, vnums2);
            const IntegrationRule & ir_voly2 = transform2(facnr_y2, ir_facet, lh);
            BaseMappedIntegrationRule & miry1 = eltransy1(ir_voly1, lh);
            BaseMappedIntegrationRule & miry2 = eltransy2(ir_voly2, lh);            
            miry1.ComputeNormalsAndMeasure (eltype1, facnr_y1);
            FlatMatrix<> elvec_xslicemat(ndofy1+ndofy2,ndofxspace*dimspace,lh);
            Array<int> dnums_xslice((ndofy1+ndofy2)*ndofxspace,lh),dnums_xslice1(ndofy2*ndofxspace,lh);
            tpfes->GetSliceDofNrs(ElementId(el1_y),0,dnums_xslice,yheap);
            tpfes->GetSliceDofNrs(ElementId(el2_y),0,dnums_xslice1,yheap);
            dnums_xslice.Append(dnums_xslice1);
            x.GetIndirect (dnums_xslice, elvec_xslicemat.AsVector());
            static_cast<TensorProductFacetBilinearFormIntegrator &>(*parts[facetvolumeintegrals]).ApplyYFacetMatrix(fely1, eltransy1, fely2,eltransy2, elvec_xslicemat, &yheap, &miry1,&miry2,lh);
            int firstxdof = 0;
            for(int j=0;j<nelx;j++)
            {
              HeapReset hr(lh);
              ElementId elid(j*nely+el1_y);
              auto & tpfel = tpfes->GetFE(elid,lh);
              // int ndofx = spaces[0]->GetFE(ElementId(j),lh).GetNDof();
              int ndofx = ndofyinv*tpfel.GetNDof();
              IntRange dnumsx(firstxdof,firstxdof+dimspace*ndofx);
              firstxdof+=dimspace*ndofx;
              const ElementTransformation & tptrafo = tpfes->GetTrafo(elid,lh);
              static_cast<TensorProductFacetBilinearFormIntegrator &>(*parts[facetvolumeintegrals]).ApplyXElementMatrix(tpfel,tptrafo,dnumsx,eltransy1.userdata,&miry1,lh);
            }
            FlatMatrix<> elmat(ndofy1+ndofy2,dimspace*ndofxspace,lh);
            elmat = 0.0;
            static_cast<TensorProductFacetBilinearFormIntegrator &>(*parts[facetvolumeintegrals]).ApplyYFacetMatrixTrans(fely1,eltransy1,fely2,eltransy2,elmat,eltransy1.userdata,&miry1,&miry2,lh);
            //elvec_mat *= val;
            y.AddIndirect (dnums_xslice, elmat.AsVector());
          } // end: if elnums.Size != 1
        } // end: inner element loop
      });
    } // end: for( auto colfacets : spaces[0]->FacetColoring )   
    timerfac2.Stop();
  } 

  template <class SCAL>
  void S_BilinearForm<SCAL> :: AddMatrixTrans (double val,
                                               const BaseVector & x,
                                               BaseVector & y, LocalHeap & clh) const
  {
    if (geom_free_parts.Size())
      AddMatrixGF(val, x, y, true, clh);
    if (geom_free_parts.Size() == parts.Size()) return;

    if (!MixedSpaces())
      {
        for (auto vb : { VOL, BND, BBND, BBBND } )
          if (VB_parts[vb].Size())
            {
              IterateElements 
                (*fespace, vb, clh, 
                 [&] (FESpace::Element el, LocalHeap & lh)
                 {
                   // RegionTimer reg (timer_loop);                   
                   auto & fel = el.GetFE();
                   auto & trafo = el.GetTrafo();
                   auto dnums = el.GetDofs();
                   
                   FlatVector<SCAL> elvecx (dnums.Size() * fespace->GetDimension(), lh);
                   FlatVector<SCAL> elvecy (dnums.Size() * fespace->GetDimension(), lh);
                   
                   x.GetIndirect (dnums, elvecx);
                   this->fespace->TransformVec (el, elvecx, TRANSFORM_SOL);
                   
                   for (auto & bfi : VB_parts[vb])
                     {
                       if (!bfi->DefinedOn (el.GetIndex())) continue;
                       if (!bfi->DefinedOnElement (el.Nr())) continue;

                       auto & mapped_trafo = trafo.AddDeformation(bfi->GetDeformation().get(), lh);
                       
                       bfi->ApplyElementMatrixTrans (fel, mapped_trafo, elvecx, elvecy, 0, lh);

                       this->fespace->TransformVec (el, elvecy, TRANSFORM_RHS);
                       
                       elvecy *= val;
                       y.AddIndirect (dnums, elvecy, fespace->HasAtomicDofs());  // coloring
                     }
                 });
            } 
      }
    else // MixedSpaces
      {
        static Timer timer ("Apply Matrix Trans - mixed"); RegionTimer reg(timer);

        for (auto vb : { VOL, BND, BBND } )
          if (VB_parts[vb].Size())
            {        
              IterateElements 
                (*fespace, vb, clh, 
                 [&] (ElementId ei, LocalHeap & lh)
           
                 {
                   if (!fespace->DefinedOn (ei)) return;
                   if (!fespace2->DefinedOn (ei)) return;
                   const FiniteElement & fel1 = fespace->GetFE (ei, lh);
                   const FiniteElement & fel2 = fespace2->GetFE (ei, lh);
                   ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
             
                   Array<int> dnums1 (fel1.GetNDof(), lh);
                   fespace->GetDofNrs (ei, dnums1);
                   Array<int> dnums2 (fel2.GetNDof(), lh);
                   fespace2->GetDofNrs (ei, dnums2);
                   
                   FlatVector<SCAL> elvecx (dnums2.Size() * fespace2->GetDimension(), lh);
                   FlatVector<SCAL> elvecy (dnums1.Size() * fespace ->GetDimension(), lh);

                   x.GetIndirect (dnums2, elvecx);
                   this->fespace2->TransformVec (ei, elvecx, TRANSFORM_SOL);

                   for (auto & bfi : VB_parts[vb])
                     {
                       if (!bfi->DefinedOn (this->ma->GetElIndex (ei))) continue;
                       if (!bfi->DefinedOnElement (ei.Nr())) continue;                        

                       MixedFiniteElement fel(fel1, fel2);
                       bfi->ApplyElementMatrixTrans (fel, eltrans, elvecx, elvecy, 0, lh);
                       
                       this->fespace->TransformVec (ei, elvecy, TRANSFORM_RHS);
        
                       elvecy *= val;
                       y.AddIndirect (dnums1, elvecy);  // coloring	      
                     }
                 });
            }
      }
    
    // throw Exception ("AddMatrixTrans only supported for geom-free");
  }


  template <class SCAL>
  void S_BilinearForm<SCAL> :: AddMatrix1 (SCAL val,
                                           const BaseVector & x,
                                           BaseVector & y, LocalHeap & clh) const
  {
    if (geom_free_parts.Size())
      AddMatrixGF(val, x, y, false, clh);
    
    static Timer timer ("Apply Matrix");
    static Timer<> timervb[4] = { string("Apply Matrix - volume"),
                                string("Apply Matrix - boundary"),
                                string("Apply Matrix - cd2"), 
                                string("Apply Matrix - cd3") };

    // static Timer timer_loop ("Apply Matrix - all loop elmat");
    static Timer timer_applyelmat ("Apply Matrix - elmat");    

    static Timer timerDG ("Apply Matrix - DG");
    using TTimer = Timer<TNoTracing, TNoTiming>;
    static TTimer timerDGpar ("Apply Matrix - DG par");
    static TTimer timerDGapply ("Apply Matrix - DG par apply");
    static TTimer timerDG1 ("Apply Matrix - DG 1");
    static TTimer timerDG2 ("Apply Matrix - DG 2");
    static TTimer timerDG2a ("Apply Matrix - DG 2a");
    static TTimer timerDG2b ("Apply Matrix - DG 2b");
    static TTimer timerDG2c ("Apply Matrix - DG 2c");
    static TTimer timerDG3 ("Apply Matrix - DG 3");
    static TTimer timerDG4 ("Apply Matrix - DG 4");
    static TTimer timerDGfacet ("Apply Matrix - DG boundary");
    static TTimer timerDGfacet1 ("Apply Matrix - DG boundary 1");
    static TTimer timerDGfacet2 ("Apply Matrix - DG boundary 2");
    static Timer timerDGparallelfacets ("Apply Matrix - DG parallel facets");
    static Timer timerspecial("Apply Matrix - Special Elements");
    RegionTimer reg (timer);

    //     static int lh_size = 5000000;
    shared_ptr<TPHighOrderFESpace> tpfes = dynamic_pointer_cast<TPHighOrderFESpace>(fespace);
    if(tpfes)
    {
      AddMatrixTP(val,x,y,clh);
      return;
    }
    
    if (!MixedSpaces())
      {
        for (auto vb : { VOL, BND, BBND, BBBND } )
          if (VB_parts[vb].Size())
            {
              RegionTimer reg (timervb[vb]);
              
              IterateElements 
                (*fespace, vb, clh, 
                 [&] (FESpace::Element el, LocalHeap & lh)
                 {
                   // RegionTimer reg (timer_loop);                   
                   auto & fel = el.GetFE();
                   auto & trafo = el.GetTrafo();
                   auto dnums = el.GetDofs();
                                      
                   FlatVector<SCAL> elvecx (dnums.Size() * fespace->GetDimension(), lh);
                   FlatVector<SCAL> elvecy (dnums.Size() * fespace->GetDimension(), lh);

                   x.GetIndirect (dnums, elvecx);
                   this->fespace->TransformVec (el, elvecx, TRANSFORM_SOL);

                   for (auto & bfi : VB_parts[vb])
                     {
                       if (!bfi->DefinedOn (el.GetIndex())) continue;
                       if (!bfi->DefinedOnElement (el.Nr())) continue;

                       auto & mapped_trafo = trafo.AddDeformation(bfi->GetDeformation().get(), lh);

                       {
                         // RegionTimer reg (timer_applyelmat);
                         bfi->ApplyElementMatrix (fel, mapped_trafo, elvecx, elvecy, 0, lh);
                       }
                       
                       this->fespace->TransformVec (el, elvecy, TRANSFORM_RHS);
                       
                       elvecy *= val;
                       y.AddIndirect (dnums, elvecy, fespace->HasAtomicDofs());  // coloring
                     }
                 });
            } 
        

        {
          RegionTimer reg(timerDG);                    
          
          if ( (facetwise_skeleton_parts[VOL].Size() > 0) ||
               (facetwise_skeleton_parts[BND].Size() > 0) )
            
            for (auto colfacets : fespace->FacetColoring())
              {
                SharedLoop2 sl(colfacets.Size());

                ParallelJob
                  ( [&] (const TaskInfo & ti) 
                    {
                      LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
                      RegionTimer reg(timerDGpar);

                      Array<int> elnums(2, lh), elnums_per(2, lh), fnums1(6, lh), fnums2(6, lh),
                        vnums1(8, lh), vnums2(8, lh);

                  for (int i : sl)                
                     {
                       // timerDG1.Start();
                       HeapReset hr(lh);
                       int facet = colfacets[i];
                       int facet2 = colfacets[i];
                       ma->GetFacetElements (facet, elnums);
                       if (elnums.Size() == 0) continue; // coarse facets
                       int el1 = elnums[0];
                       ElementId ei1(VOL, el1);
                       fnums1 = ma->GetElFacets(ei1);
                       int facnr1 = fnums1.Pos(facet);
                       
                       // timerDG1.Stop();
                       if(elnums.Size() < 2)
                         {
                           if (ma->GetCommunicator().Size() > 1)
                             if (ma->GetDistantProcs (NodeId(NT_FACET, facet)).Size() > 0)
                               continue;

                           facet2 = ma->GetPeriodicFacet(facet);
                           if(facet2 > facet)
                             {
                               ma->GetFacetElements (facet2, elnums_per);
                               if (elnums_per.Size() > 1)
                                 throw Exception("DG-Apply failed due to invalid periodicity.");
                               elnums.Append(elnums_per[0]);
                             }
                           else if(facet2 < facet)
                             continue;
                         }
                       
                       if (elnums.Size()<2)
                         {
                           // RegionTimer reg(timerDGfacet);
                           
                           ma->GetFacetSurfaceElements (facet, elnums);
                           int sel = elnums[0];
                           ElementId sei(BND, sel);
                           
                           const FiniteElement & fel = fespace->GetFE (ei1, lh);
                           Array<int> dnums(fel.GetNDof(), lh);
                           vnums1 = ma->GetElVertices (ei1);
                           vnums2 = ma->GetElVertices (sei);
                           
                           ElementTransformation & eltrans = ma->GetTrafo (ei1, lh);
                           ElementTransformation & seltrans = ma->GetTrafo (sei, lh);
                           
                           fespace->GetDofNrs (ei1, dnums);
                           
                           for (auto & bfi : facetwise_skeleton_parts[BND])
                             {
                               if (!bfi->DefinedOn (seltrans.GetElementIndex())) continue;
                               if (!bfi->DefinedOnElement (facet)) continue;
                                         
                               FlatVector<SCAL> elx(dnums.Size()*this->fespace->GetDimension(), lh),
                                 ely(dnums.Size()*this->fespace->GetDimension(), lh);
                               x.GetIndirect(dnums, elx);
                               this->fespace->TransformVec (ei1, elx, TRANSFORM_SOL);
                               
                               auto & mapped_trafo = eltrans.AddDeformation(bfi->GetDeformation().get(), lh);
                               auto & mapped_strafo = seltrans.AddDeformation(bfi->GetDeformation().get(), lh);
                               bfi->ApplyFacetMatrix (fel,facnr1,mapped_trafo,vnums1, mapped_strafo, vnums2, elx, ely, lh);
                               this->fespace->TransformVec (ei1, ely, TRANSFORM_RHS);
                               y.AddIndirect(dnums, ely, fespace->HasAtomicDofs());
                             } //end for (numintegrators)
                           
                           continue;
                         } // end if boundary facet
                       
                       if (facetwise_skeleton_parts[VOL].Size() == 0)
                         continue;
                       
                       // timerDG2.Start();
                       // timerDG2a.Start();
                       // int el2 = elnums[1];
                       // ElementId ei2(VOL, el2);
                       ElementId ei2(VOL, elnums[1]);
                       
                       // fnums2 = ma->GetElFacets(ei2);
                       // int facnr2 = fnums2.Pos(facet2);
                       int facnr2 = ma->GetElFacets(ei2).Pos(facet2);

                       ElementTransformation & eltrans1 = ma->GetTrafo (ei1, lh);
                       ElementTransformation & eltrans2 = ma->GetTrafo (ei2, lh);
                       
                       const FiniteElement & fel1 = fespace->GetFE (ei1, lh);
                       const FiniteElement & fel2 = fespace->GetFE (ei2, lh);
                       // timerDG2a.Stop();
                       // timerDG2b.Start();                                 
                       Array<int> dnums1(fel1.GetNDof(), lh);
                       Array<int> dnums2(fel2.GetNDof(), lh);
                       fespace->GetDofNrs (ei1, dnums1);
                       fespace->GetDofNrs (ei2, dnums2);
                       vnums1 = ma->GetElVertices (ei1);
                       vnums2 = ma->GetElVertices (ei2);
                       
                       Array<int> dnums(fel1.GetNDof()+fel2.GetNDof(), lh);
                       dnums.Range(0, dnums1.Size()) = dnums1;
                       dnums.Range(dnums1.Size(), dnums.Size()) = dnums2;
                       FlatVector<SCAL> elx(dnums.Size()*fespace->GetDimension(), lh),
                         ely(dnums.Size()*fespace->GetDimension(), lh);
                       
                       x.GetIndirect(dnums, elx);
                       this->fespace->TransformVec (ei1, elx.Range(0, dnums1.Size()*fespace->GetDimension()), TRANSFORM_SOL);
                       this->fespace->TransformVec (ei2, elx.Range(dnums1.Size()*fespace->GetDimension(), dnums.Size()*fespace->GetDimension()), TRANSFORM_SOL);

                       RegionTimer reg2(timerDGapply);                     
                       for (auto & bfi : facetwise_skeleton_parts[VOL])                                   
                         {
                           if (!bfi->DefinedOn (ma->GetElIndex (ei1))) continue; 
                           if (!bfi->DefinedOn (ma->GetElIndex (ei2))) continue; 
                           if (!bfi->DefinedOnElement (facet) ) continue;
                           
                           auto & mapped_trafo1 = eltrans1.AddDeformation(bfi->GetDeformation().get(), lh);
                           auto & mapped_trafo2 = eltrans2.AddDeformation(bfi->GetDeformation().get(), lh);
                           bfi->ApplyFacetMatrix (fel1, facnr1, mapped_trafo1, vnums1,
                                                  fel2, facnr2, mapped_trafo2, vnums2, elx, ely, lh);
                           this->fespace->TransformVec (ei1, ely.Range(0, dnums1.Size()*fespace->GetDimension()), TRANSFORM_RHS);
                           this->fespace->TransformVec (ei2, ely.Range(dnums1.Size()*fespace->GetDimension(), dnums.Size()*fespace->GetDimension()), TRANSFORM_RHS);

                           y.AddIndirect(dnums, ely);
                         }
                     }
                 });
              }
          
                    

          if (elementwise_skeleton_parts.Size())
            IterateElements 
              (*fespace, VOL, clh, 
               [&] (ElementId ei1, LocalHeap & lh)
               {
                 {
                   int el1 = ei1.Nr();
                   Array<int> elnums(2, lh), elnums_per(2, lh), fnums1(6, lh), fnums2(6, lh),
                     vnums1(8, lh), vnums2(8, lh);
                   // RegionTimer reg1(timerDG1);
                   
                   fnums1 = ma->GetElFacets(ei1);
                   
                   for (int facnr1 : Range(fnums1))
                     {
                       HeapReset hr(lh);
                       int facet = fnums1[facnr1];
                       int facet2 = fnums1[facnr1];

                       ma->GetFacetElements(facet,elnums);
                       if (elnums.Size()<2) {
                         // #ifdef PARALLEL
                         auto comm = ma->GetCommunicator();
			 if( (comm.Size()>1) && (ma->GetDistantProcs (NodeId(StdNodeType(NT_FACET, ma->GetDimension()), fnums1[facnr1])).Size() > 0) )
			   continue;
                         // #endif
                         facet2 = ma->GetPeriodicFacet(fnums1[facnr1]);
                         if(facet2!=facet)
                           {
                             ma->GetFacetElements (facet2, elnums_per);
                             if (elnums_per.Size() > 1)
                               throw Exception("DG-Apply failed due to invalid periodicity.");
                             elnums.Append(elnums_per[0]);
                           }
                       }

                       if (elnums.Size()<2)
                         {
                           ma->GetFacetSurfaceElements (fnums1[facnr1], elnums);
                           int sel = elnums[0];
                           ElementId sei(BND, sel);
                           const FiniteElement & fel = fespace->GetFE (ei1, lh);
                           Array<int> dnums(fel.GetNDof(), lh);
                           vnums1 = ma->GetElVertices (ei1);
                           vnums2 = ma->GetElVertices (sei);     
                           
                           ElementTransformation & eltrans = ma->GetTrafo (ei1, lh);
                           ElementTransformation & seltrans = ma->GetTrafo (sei, lh);
                           
                           fespace->GetDofNrs (ei1, dnums);
                           if(fel.GetNDof() != dnums.Size())
                             {

                               (*testout) << "fel::GetNDof() = " << fel.GetNDof() << endl;
                               (*testout) << "dnums.Size() = " << dnums.Size() << endl;
                               (*testout) << "dnums = " << dnums << endl;
                               throw Exception ( string("Inconsistent number of degrees of freedom Surface fel::GetNDof() = ") + ToString(fel.GetNDof()) + string(" != dnums.Size() = ") + ToString(dnums.Size()) + string("!") );
                             }
                           
                           
                           for (auto & bfi : elementwise_skeleton_parts)
                             {
                               if (!bfi->DefinedOnElement (el1) ) continue;
                               FlatVector<SCAL> elx(dnums.Size()*fespace->GetDimension(), lh),
                                 ely(dnums.Size()*fespace->GetDimension(), lh);
                               x.GetIndirect(dnums, elx);
                               
                               // dynamic_cast<const FacetBilinearFormIntegrator&>(*bfi).
                               auto & mapped_trafo = eltrans.AddDeformation(bfi->GetDeformation().get(), lh);
                               auto & mapped_strafo = seltrans.AddDeformation(bfi->GetDeformation().get(), lh);
                               bfi->ApplyFacetMatrix (fel,facnr1,mapped_trafo,vnums1, mapped_strafo, vnums2, elx, ely, lh);
                               
                               y.AddIndirect(dnums, ely);
                               
                             } //end for (numintegrators)
                           
                           continue;
                         } // end if boundary facet

                       
                       // RegionTimer reg2(timerDG2);
                       int el2 = elnums[0] + elnums[1] - el1;
                       // T_ElementId<VOL,2> ei2(el2);
                       ElementId ei2(VOL, el2);
                       
                       fnums2 = ma->GetElFacets(ei2);
                       int facnr2 = fnums2.Pos(facet2);
                       
                       ElementTransformation & eltrans1 = ma->GetTrafo (ei1, lh);
                       ElementTransformation & eltrans2 = ma->GetTrafo (ei2, lh);
                       
                       const FiniteElement & fel1 = fespace->GetFE (ei1, lh);
                       const FiniteElement & fel2 = fespace->GetFE (ei2, lh);
                                 
                       Array<int> dnums1(fel1.GetNDof(), lh);
                       Array<int> dnums2(fel2.GetNDof(), lh);
                       fespace->GetDofNrs (ei1, dnums1);
                       fespace->GetDofNrs (ei2, dnums2);
                       
                       vnums1 = ma->GetElVertices (ei1);
                       vnums2 = ma->GetElVertices (ei2);
                       
                       if(fel1.GetNDof() != dnums1.Size() || ((elnums.Size()>1) && (fel2.GetNDof() != dnums2.Size() )))
                         {
                           cout << "facet, neighbouring fel(1): GetNDof() = " << fel1.GetNDof() << endl;
                           cout << "facet, neighbouring fel(2): GetNDof() = " << fel2.GetNDof() << endl;
                           cout << "facet, neighbouring fel(1): dnums.Size() = " << dnums1.Size() << endl;
                           cout << "facet, neighbouring fel(2): dnums.Size() = " << dnums2.Size() << endl;
                           throw Exception ( "Inconsistent number of degrees of freedom " );
                         }
                       
                       Array<int> dnums(fel1.GetNDof()+fel2.GetNDof(), lh);
                       /*
                         dnums.SetSize0();
                         dnums.Append(dnums1);
                         dnums.Append(dnums2);   
                       */
                       dnums.Range(0, dnums1.Size()) = dnums1;
                       dnums.Range(dnums1.Size(), dnums.Size()) = dnums2;
                       
                       FlatVector<SCAL> elx(dnums.Size()*fespace->GetDimension(), lh),
                         ely(dnums.Size()*fespace->GetDimension(), lh);
                       x.GetIndirect(dnums, elx);
                       
                       
                       for (auto & bfi : elementwise_skeleton_parts)
                         {
                           if (!bfi->DefinedOn (ma->GetElIndex (ei1))) continue; //TODO: treat as surface element
                           if (!bfi->DefinedOn (ma->GetElIndex (ei2))) continue; //TODO    
                           if (!bfi->DefinedOnElement (el1) ) continue;
                           // FacetBilinearFormIntegrator * fbfi = 
                           // dynamic_cast<FacetBilinearFormIntegrator*>(bfi.get());
                           
                           
                           auto & mapped_trafo1 = eltrans1.AddDeformation(bfi->GetDeformation().get(), lh);
                           auto & mapped_trafo2 = eltrans2.AddDeformation(bfi->GetDeformation().get(), lh);
                           bfi->ApplyFacetMatrix (fel1, facnr1, mapped_trafo1, vnums1,
                                                  fel2, facnr2, mapped_trafo2, vnums2, elx, ely, lh);
                           
                           /*
                             if (neighbor_testfunction)
                             {
                             lock_guard<mutex> guard(addelemfacbnd_mutex);
                             y.AddIndirect(dnums, ely);
                             }
                             else
                             {
                             y.AddIndirect(dnums1, ely.Range(0,dnums1.Size()));
                                       }
                           */
                           y.AddIndirect(dnums1, ely.Range(0,dnums1.Size()));
                           
                           if (bfi->GetDGFormulation().neighbor_testfunction)
                             {
                               int dim = this->fespace->GetDimension();
                               FlatVector<SCAL> swap_elx(elx.Size(), lh);
                               swap_elx.Range(0, dim*dnums2.Size()) = elx.Range(dim*dnums1.Size(), dim*dnums.Size());
                               swap_elx.Range(dim*dnums2.Size(), dim*dnums.Size()) = elx.Range(0, dim*dnums1.Size());
                               auto & mapped_trafo1 = eltrans1.AddDeformation(bfi->GetDeformation().get(), lh);
                               auto & mapped_trafo2 = eltrans2.AddDeformation(bfi->GetDeformation().get(), lh);
                               bfi->ApplyFacetMatrix (fel2, facnr2, mapped_trafo2, vnums2,
                                                      fel1, facnr1, mapped_trafo1, vnums1, swap_elx, ely, lh);
                               y.AddIndirect(dnums1, ely.Range(dim*dnums2.Size(), dim*dnums.Size()));
                             }
                         }
                     }
                 }                             
               });
        }

	
        // #ifdef PARALLEL
        auto comm = ma->GetCommunicator();
	if (comm.Size() > 1 && mpi_facet_parts.Size())
	  {
	    RegionTimer rt(timerDGparallelfacets);
        HeapReset hr(clh);
	    
	    //cout << "apply parallel DG facets, " << elementwise_skeleton_parts.Size() << " el-bound and " << facetwise_skeleton_parts[VOL].Size() << " facet parts" << ", " << mpi_facet_parts.Size() << " total parts " << endl;

            // int mrank = comm.Rank();
            int mnp = comm.Size();
	    Array<int> cnt_in(mnp), cnt_per(mnp);
	    if(!have_mpi_facet_data) {
	      os_per = Array<int>(mnp);
	      os_per = 0;
	    }
	    Array<MPI_Request> reqs;
	    Array<MPI_Request> reqr;
	    LocalHeap &lh(clh);
	    Array<int> elnums(2, lh), elnums2(2, lh), fnums(6, lh), vnums(8, lh);

	    size_t ne = ma->GetNE(VOL);
	    BitArray fine_facet(ma->GetNFacets());
	    fine_facet.Clear();
	    Array<int> elfacets;
	    for (int i = 0; i < ne; ++i) {
	      auto elfacets = ma->GetElFacets(ElementId(VOL,i));
	      for (auto f : elfacets) fine_facet.SetBit(f);
	    }
	    /** 
		We can have surf-els without elements. If we just skip these, 
		the order of facets can be jumbled between ranks.
		So we include facets from surf-els and jump to the identified facet
		(+ its local vol el) if on mpi-bnd.
	    **/
	    size_t nse = ma->GetNE(BND);
	    for (int i = 0; i < nse; ++i) {
	      auto selfacets = ma->GetElFacets(ElementId(BND, i));
	      for (auto f : selfacets) fine_facet.SetBit(f);
	    }

	    auto mpi_loop_range = (have_mpi_facet_data)?Range(1,3):Range(0,3);
	    
	    for(auto loop:mpi_loop_range) {
	      cnt_in = 0;
	      cnt_per = 0;
	      for(auto facet:Range(ma->GetNFacets())) {
		if(!fine_facet.Test(facet)) continue;
		NodeId facet_id(StdNodeType(NT_FACET, ma->GetDimension()), facet);
		auto fdps = ma->GetDistantProcs(facet_id);
		//skip non-mpi facets
		if (fdps.Size() == 0) continue;
		auto d = fdps[0];
		HeapReset hr(lh);

		ma->GetFacetSurfaceElements (facet, elnums2);
		ma->GetFacetElements(facet, elnums);
		bool periodic_facet = elnums2.Size()!=0;
		if (periodic_facet) { // dont double up on facets!
		  auto facet2 = ma->GetPeriodicFacet(facet);
		  if (facet>facet2) continue;
		}
		if (periodic_facet && !elnums.Size()) // use the identified facet!
		  {
		    facet = ma->GetPeriodicFacet(facet);
		    ma->GetFacetElements(facet, elnums);
		    ma->GetFacetSurfaceElements (facet, elnums2);
		  }

		ElementId eiv(VOL, elnums[0]);

		fnums = ma->GetElFacets(eiv);
		int facetnr = fnums.Pos(facet);

		const FiniteElement & fel = fespace->GetFE (eiv, lh);
		
		
		Array<int> dnums(fel.GetNDof(), lh);
		vnums = ma->GetElVertices(eiv);

		ElementTransformation & eltrans = ma->GetTrafo (eiv, lh);
		fespace->GetDofNrs (eiv, dnums);

		for(auto igt:mpi_facet_parts) {

		  FlatVector<SCAL> elx(dnums.Size()*this->fespace->GetDimension(), lh);
		  x.GetIndirect(dnums, elx);
		  FlatVector<SCAL> trace_values;
		  dynamic_cast<const FacetBilinearFormIntegrator*>(igt.get())->  
		    CalcTraceValues(fel,facetnr,eltrans,vnums, trace_values, elx, lh);
		  if (loop == 0) {
		    os_per[d] += trace_values.Size();
		    if(periodic_facet) cnt_per[d] += trace_values.Size();
		    else cnt_in[d] += trace_values.Size();
		  }
		  else if (loop == 1) {
		    auto offset = periodic_facet ? (os_per[d] + cnt_per[d]) : cnt_in[d];
		    FlatVector<SCAL> tmp(trace_values.Size(), &( send_table[d][offset] ));
		    tmp = trace_values;
		    if(periodic_facet) cnt_per[d] += trace_values.Size();
		    else cnt_in[d] += trace_values.Size();
		  }
		  else {
		    auto offset = periodic_facet ? (os_per[d] + cnt_per[d]) : cnt_in[d];
		    FlatVector<SCAL> trace_other(trace_values.Size(), &( recv_table[d][offset] ));
		    if(periodic_facet) cnt_per[d] += trace_values.Size();
		    else cnt_in[d] += trace_values.Size();

		    FlatVector<SCAL> ely(dnums.Size()*this->fespace->GetDimension(), lh);
		    dynamic_cast<const FacetBilinearFormIntegrator*>(igt.get())->  
		      ApplyFromTraceValues(fel,facetnr,eltrans,vnums, trace_other,  elx, ely, lh);
			
		    y.AddIndirect(dnums, ely);
		  }
		}		
	      }
	      
	      if(loop==0) {
		send_table = Table<SCAL> (os_per);
		recv_table = Table<SCAL> (os_per);
		os_per = cnt_in;
		for(auto r:send_table)
		  r = -1;
		for(auto r:recv_table)
		  r = -2;
		have_mpi_facet_data = true;
	      }
	      else if(loop==1) {
		for(auto dp:Range(mnp))
		  if(send_table[dp].Size()) {
		    reqs.Append(comm.ISend(send_table[dp], dp, MPI_TAG_SOLVE));
		    reqr.Append(comm.IRecv(recv_table[dp], dp, MPI_TAG_SOLVE));
		  }
		MyMPI_WaitAll(reqr);
	      }
	    }
	    if(reqs.Size()) MyMPI_WaitAll(reqs);
	  }
        // #endif

        static mutex specelmutex;
        if (specialelements.Size())
          {
            RegionTimer regt(timerspecial);
            ParallelForRange(IntRange(specialelements.Size()), [&](IntRange r)
            {
              Array<int> dnums;
              LocalHeap lh = clh.Split();
              for(auto i : r)
                {
                  HeapReset hr(lh);
                  const SpecialElement & el = *specialelements[i];
                  el.GetDofNrs (dnums);
                  FlatVector<SCAL> elvecx (dnums.Size() * fespace->GetDimension(), lh);
                  FlatVector<SCAL> elvecy (dnums.Size() * fespace->GetDimension(), lh);
                  x.GetIndirect (dnums, elvecx);
                  el.Apply (elvecx, elvecy, lh);
                  elvecy *= val;
                  {
                    lock_guard<mutex> lock(specelmutex);
                    y.AddIndirect (dnums, elvecy);
                  }
                }
            });
          }
        /*
              }
            catch (LocalHeapOverflow lhex)
              {
                lh_size *= 2;
                atempt++;
                cerr << "Trying automatic heapsize increase to " << lh_size << endl;
              }
          }
                */
      }
    else // MixedSpaces
      {
        static Timer timer ("Apply Matrix - mixed"); RegionTimer reg(timer);

        for (auto vb : { VOL, BND, BBND } )
          if (VB_parts[vb].Size())
            {        
              IterateElements 
                (*fespace2, vb, clh, 
                 [&] (ElementId ei, LocalHeap & lh)
           
                 {
                   if (!fespace->DefinedOn (ei)) return;
                   if (!fespace2->DefinedOn (ei)) return;
                   const FiniteElement & fel1 = fespace->GetFE (ei, lh);
                   const FiniteElement & fel2 = fespace2->GetFE (ei, lh);
                   ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
             
                   Array<int> dnums1 (fel1.GetNDof(), lh);
                   fespace->GetDofNrs (ei, dnums1);
                   Array<int> dnums2 (fel2.GetNDof(), lh);
                   fespace2->GetDofNrs (ei, dnums2);
                   
                   FlatVector<SCAL> elvecx (dnums1.Size() * fespace->GetDimension(), lh);
                   FlatVector<SCAL> elvecy (dnums2.Size() * fespace2->GetDimension(), lh);

                   x.GetIndirect (dnums1, elvecx);
                   this->fespace->TransformVec (ei, elvecx, TRANSFORM_SOL);

                   for (auto & bfi : VB_parts[vb])
                     {
                       if (!bfi->DefinedOn (this->ma->GetElIndex (ei))) continue;
                       if (!bfi->DefinedOnElement (ei.Nr())) continue;                        

                       MixedFiniteElement fel(fel1, fel2);
                       bfi->ApplyElementMatrix (fel, eltrans, elvecx, elvecy, 0, lh);
                       
                       this->fespace2->TransformVec (ei, elvecy, TRANSFORM_RHS);
        
                       elvecy *= val;
                       y.AddIndirect (dnums2, elvecy);  // coloring	      
                     }
                 });
            }

        {
          RegionTimer reg(timerDG);

          if ( (facetwise_skeleton_parts[VOL].Size() > 0) ||
               (facetwise_skeleton_parts[BND].Size() > 0) )

            for (auto colfacets : fespace->FacetColoring())
              {
                SharedLoop2 sl(colfacets.Size());

                ParallelJob
                  ( [&] (const TaskInfo & ti)
                    {
                      LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
                      RegionTimer reg(timerDGpar);

                      Array<int> elnums(2, lh), elnums_per(2, lh), fnums1(6, lh), fnums2(6, lh),
                        vnums1(8, lh), vnums2(8, lh);

                  for (int i : sl)
                     {
                       // timerDG1.Start();
                       HeapReset hr(lh);
                       int facet = colfacets[i];
                       int facet2 = colfacets[i];
                       ma->GetFacetElements (facet, elnums);
                       if (elnums.Size() == 0) continue; // coarse facets
                       int el1 = elnums[0];
                       ElementId ei1(VOL, el1);
                       fnums1 = ma->GetElFacets(ei1);
                       int facnr1 = fnums1.Pos(facet);

                       // timerDG1.Stop();
                       if(elnums.Size() < 2)
                         {
                           if (ma->GetCommunicator().Size() > 1)
                             if (ma->GetDistantProcs (NodeId(NT_FACET, facet)).Size() > 0)
                               continue;

                           facet2 = ma->GetPeriodicFacet(facet);
                           if(facet2 > facet)
                             {
                               ma->GetFacetElements (facet2, elnums_per);
                               if (elnums_per.Size() > 1)
                                 throw Exception("DG-Apply failed due to invalid periodicity.");
                               elnums.Append(elnums_per[0]);
                             }
                           else if(facet2 < facet)
                             continue;
                         }

                       if (elnums.Size()<2)
                         {
                           // RegionTimer reg(timerDGfacet);

                           ma->GetFacetSurfaceElements (facet, elnums);
                           int sel = elnums[0];
                           ElementId sei(BND, sel);

                           const FiniteElement & fel1 = fespace->GetFE (ei1, lh);
                           const FiniteElement & fel2 = fespace2->GetFE (ei1, lh);
                           MixedFiniteElement fel(fel1, fel2);

                           vnums1 = ma->GetElVertices (ei1);
                           vnums2 = ma->GetElVertices (sei);

                           ElementTransformation & eltrans = ma->GetTrafo (ei1, lh);
                           ElementTransformation & seltrans = ma->GetTrafo (sei, lh);

                           Array<int> dnums1 (fel1.GetNDof(), lh);
                           fespace->GetDofNrs (ei1, dnums1);
                           Array<int> dnums2 (fel2.GetNDof(), lh);
                           fespace2->GetDofNrs (ei1, dnums2);

                           for (auto & bfi : facetwise_skeleton_parts[BND])
                             {
                               if (!bfi->DefinedOn (seltrans.GetElementIndex())) continue;
                               if (!bfi->DefinedOnElement (facet)) continue;

                               FlatVector<SCAL> elx(dnums1.Size()*this->fespace->GetDimension(), lh),
                                 ely(dnums2.Size()*this->fespace->GetDimension(), lh);
                               x.GetIndirect(dnums1, elx);
                               this->fespace->TransformVec (ei1, elx, TRANSFORM_SOL);

                               auto & mapped_trafo = eltrans.AddDeformation(bfi->GetDeformation().get(), lh);
                               auto & mapped_strafo = seltrans.AddDeformation(bfi->GetDeformation().get(), lh);

                               bfi->ApplyFacetMatrix (fel,facnr1,mapped_trafo,vnums1, mapped_strafo, vnums2, elx, ely, lh);
                               this->fespace2->TransformVec (ei1, ely, TRANSFORM_RHS);
                               y.AddIndirect(dnums2, ely, fespace2->HasAtomicDofs());
                             } //end for (numintegrators)

                           continue;
                         } // end if boundary facet

                       if (facetwise_skeleton_parts[VOL].Size() == 0)
                         continue;

                       // timerDG2.Start();
                       // timerDG2a.Start();
                       // int el2 = elnums[1];
                       // ElementId ei2(VOL, el2);
                       ElementId ei2(VOL, elnums[1]);

                       // fnums2 = ma->GetElFacets(ei2);
                       // int facnr2 = fnums2.Pos(facet2);
                       int facnr2 = ma->GetElFacets(ei2).Pos(facet2);

                       ElementTransformation & eltrans1 = ma->GetTrafo (ei1, lh);
                       ElementTransformation & eltrans2 = ma->GetTrafo (ei2, lh);

                       const FiniteElement & fel1_trial = fespace->GetFE (ei1, lh);
                       const FiniteElement & fel1_test = fespace2->GetFE (ei1, lh);
                       MixedFiniteElement fel1(fel1_trial, fel1_test);
                       const FiniteElement & fel2_trial = fespace->GetFE (ei2, lh);
                       const FiniteElement & fel2_test = fespace2->GetFE (ei2, lh);
                       MixedFiniteElement fel2(fel2_trial, fel2_test);

                       // timerDG2a.Stop();
                       // timerDG2b.Start();
                       Array<int> dnums1_trial(fel1_trial.GetNDof(), lh), dnums2_trial(fel2_trial.GetNDof(), lh);
                       Array<int> dnums1_test(fel1_test.GetNDof(), lh), dnums2_test(fel2_trial.GetNDof(), lh);
                       fespace->GetDofNrs (ei1, dnums1_trial);
                       fespace->GetDofNrs (ei2, dnums2_trial);
                       fespace2->GetDofNrs (ei1, dnums1_test);
                       fespace2->GetDofNrs (ei2, dnums2_test);
                       vnums1 = ma->GetElVertices (ei1);
                       vnums2 = ma->GetElVertices (ei2);

                       Array<int> dnums_trial(fel1_trial.GetNDof()+fel2_trial.GetNDof(), lh);
                       dnums_trial.Range(0, dnums1_trial.Size()) = dnums1_trial;
                       dnums_trial.Range(dnums1_trial.Size(), dnums_trial.Size()) = dnums2_trial;
                       Array<int> dnums_test(fel1_test.GetNDof()+fel2_test.GetNDof(), lh);
                       dnums_test.Range(0, dnums1_test.Size()) = dnums1_test;
                       dnums_test.Range(dnums1_test.Size(), dnums_test.Size()) = dnums2_test;

                       FlatVector<SCAL> elx(dnums_trial.Size()*fespace->GetDimension(), lh),
                         ely(dnums_test.Size()*fespace2->GetDimension(), lh);

                       x.GetIndirect(dnums_trial, elx);
                       this->fespace->TransformVec (ei1, elx.Range(0, dnums1_trial.Size()*fespace->GetDimension()), TRANSFORM_SOL);
                       this->fespace->TransformVec (ei2, elx.Range(dnums1_trial.Size()*fespace->GetDimension(), dnums_trial.Size()*fespace->GetDimension()), TRANSFORM_SOL);

                       RegionTimer reg2(timerDGapply);
                       for (auto & bfi : facetwise_skeleton_parts[VOL])
                         {
                           if (!bfi->DefinedOn (ma->GetElIndex (ei1))) continue;
                           if (!bfi->DefinedOn (ma->GetElIndex (ei2))) continue;
                           if (!bfi->DefinedOnElement (facet) ) continue;

                           auto & mapped_trafo1 = eltrans1.AddDeformation(bfi->GetDeformation().get(), lh);
                           auto & mapped_trafo2 = eltrans2.AddDeformation(bfi->GetDeformation().get(), lh);
                           bfi->ApplyFacetMatrix (fel1, facnr1, mapped_trafo1, vnums1,
                                                  fel2, facnr2, mapped_trafo2, vnums2, elx, ely, lh);
                           this->fespace2->TransformVec (ei1, ely.Range(0, dnums1_test.Size()*fespace2->GetDimension()), TRANSFORM_RHS);
                           this->fespace2->TransformVec (ei2, ely.Range(dnums1_test.Size()*fespace2->GetDimension(), dnums_test.Size()*fespace2->GetDimension()), TRANSFORM_RHS);

                           y.AddIndirect(dnums_test, ely);
                         }
                     }
                 });
              }

          if (elementwise_skeleton_parts.Size())
              throw Exception("Elementwise skeleton not yet implemented for Mixed Spaces + Apply!");
          auto comm = ma->GetCommunicator();
          if (comm.Size() > 1 && mpi_facet_parts.Size())
              throw Exception("MPI facets not yet implemented for Mixed Spaces + Apply!");
        }
      }
  }


  template <class SCAL>
  void S_BilinearForm<SCAL> :: AddMatrixGF (SCAL val,
                                            const BaseVector & x,
                                            BaseVector & y,
                                            bool transpose,
                                            LocalHeap & lh) const
  {
    static Timer t("BilinearForm::Apply - geomfree");
    static Timer tgetx("BilinearForm::Apply - get x");    
    static Timer tx("BilinearForm::Apply - transform x");
    static Timer ty("BilinearForm::Apply - transform y");
    static Timer taddy("BilinearForm::Apply - add y");        
    static Timer tgf("BilinearForm::Apply - geomfree gridfunction");
    static Timer tgfmult("BilinearForm::Apply - geomfree gridfunction - mult");
    static Timer tm("BilinearForm::Apply - geomfree mult");
    static Timer teval("BilinearForm::Apply - evaluate");
    RegionTimer reg(t);

    // classify:
    auto fesx = GetTrialSpace();
    auto fesy = GetTestSpace();
    if (transpose) Swap (fesx, fesy);
    auto ma = GetMeshAccess();

    /*
    Array<short> classnr(ma->GetNE());
    ma->IterateElements
      (VOL, lh, [&] (auto el, LocalHeap & llh)
       {
         classnr[el.Nr()] = 
           SwitchET<ET_TRIG,ET_TET>
           (el.GetType(),
            [el] (auto et) { return ET_trait<et.ElementType()>::GetClassNr(el.Vertices()); });
       });
    
    TableCreator<size_t> creator;
    for ( ; !creator.Done(); creator++)
      for (auto i : Range(classnr))
        creator.Add (classnr[i], i);
    Table<size_t> table = creator.MoveTable();
    */
    
    const Table<size_t> & table = ma->GetElementsOfClass();
    // auto dof_tablex = fesx->CreateDofTable(VOL);



#ifdef OLDPARALLEL
    
    for (auto elclass_inds : table)
      {
        if (elclass_inds.Size() == 0) continue;
        HeapReset hr(lh);
        
        ElementId ei(VOL,elclass_inds[0]);
        auto & felx = GetTrialSpace()->GetFE (ei, lh);
        auto & fely = GetTestSpace()->GetFE (ei, lh);
        auto & trafo = GetTrialSpace()->GetMeshAccess()->GetTrafo(ei, lh);

        Matrix<> melx(elclass_inds.Size(), felx.GetNDof()*fesx->GetDimension());
        Matrix<> mely(elclass_inds.Size(), fely.GetNDof()*fesy->GetDimension());
        mely = 0.0;

        {
          RegionTimer reg(tgetx);
          ParallelForRange
            (elclass_inds.Size(), [&] (IntRange myrange)
             {
               Array<DofId> dofnr;
               for (auto i : myrange)
                 {
                   fesx->GetDofNrs( { VOL, elclass_inds[i] }, dofnr);
                   x.GetIndirect(dofnr, melx.Row(i));
                   // x.GetIndirect(dof_tablex[elclass_inds[i]], melx.Row(i));
                 }
             });
        }
        
        for (auto bfi1 : geom_free_parts)
          {
            auto bfi = dynamic_pointer_cast<SymbolicBilinearFormIntegrator> (bfi1);
            
            auto & trial_proxies = bfi->TrialProxies();
            auto & test_proxies = bfi->TestProxies();
            auto & gridfunction_cfs = bfi->GridFunctionCoefficients();
            auto & cf = bfi->GetCoefficientFunction();
            
            VorB element_vb = bfi->ElementVB();
            Facet2ElementTrafo f2el(felx.ElementType(), bfi->ElementVB());


            
            Array<Matrix<>> mgfs;
            {
              RegionTimer reg(tgf);            
              for (CoefficientFunction * cf : gridfunction_cfs)
                {
                  auto gfcf = dynamic_cast<GridFunctionCoefficientFunction*> (cf);
                  if (!gfcf)
                    {
                      cout << "no gf" << endl;
                      continue;
                    }
                  
                  auto fes = gfcf->GetGridFunction().GetFESpace();
                  auto & felgf = fes->GetFE(ei, lh);
                  Matrix<> mgf(elclass_inds.Size(), felgf.GetNDof()*fes->GetDimension());
                  auto & vec = gfcf->GetGridFunction().GetVector();
                  
                  /*
                  auto diffop = gfcf->GetDifferentialOperator(trafo.VB());
                  
                  FlatMatrix<SIMD<double>> bmat(felgf.GetNDof()*fes->GetDimension()*
                                                diffop->Dim(), // ??? right Dim
                                                simd_ir.Size(), lh);
                  FlatMatrix<double> hbmat(felgf.GetNDof()*fes->GetDimension(),
                                           diffop->Dim()*    // right Dim ???
                                           simd_ir.Size()*SIMD<double>::Size(),
                                           &bmat(0)[0]);
                  bmat = SIMD<double> (0.0);
                  diffop->CalcMatrix(felgf, simd_mir, bmat);
                  
                  Matrix<SIMD<double>> hmgfxi(elclass_inds.Size(), diffop->Dim()*simd_ir.Size());
                  FlatMatrix<> hhmgfxi(hmgfxi.Height(), hmgfxi.Width()*SIMD<double>::Size(), &hmgfxi(0)[0]);                      
                  */  
                  ParallelForRange
                    (elclass_inds.Size(), [&] (IntRange myrange) {
                      Array<DofId> dofnr;
                      for (auto i : myrange)
                        {
                          fes->GetDofNrs( { VOL, elclass_inds[i] }, dofnr);
                          vec.GetIndirect(dofnr, mgf.Row(i));                    
                        }
                    });
                      
                  mgfs.Append (move(mgf));
                }
            }
            

            
            for (auto facet : Range(f2el.GetNFacets()))
              {
                const SIMD_IntegrationRule & simd_ir = (element_vb == VOL) ? 
                  bfi->Get_SIMD_IntegrationRule (felx, lh)
                  :
                  f2el(facet,
                       bfi->GetSIMDIntegrationRule(f2el.FacetType (facet),
                                                   felx.Order()+fely.Order()+bfi->GetBonusIntegrationOrder()),
                       lh);
                
                auto & simd_mir = trafo(simd_ir, lh);

                Array<Matrix<SIMD<double>>> melxi;
                {
                  RegionTimer regx(tx);
                  for (auto proxynr : Range(trial_proxies))
                    {
                      auto proxy = trial_proxies[proxynr];
                      FlatMatrix<SIMD<double>> bmatx(melx.Width()*proxy->Dimension(),
                                                     simd_ir.Size(), lh);
                      FlatMatrix<double> hbmatx(melx.Width(),
                                                proxy->Dimension()*simd_ir.Size()*SIMD<double>::Size(),
                                                &bmatx(0)[0]);

                      bmatx = SIMD<double> (0.0);
                      proxy->Evaluator()->CalcMatrix(felx, simd_mir, bmatx);
                      Matrix<SIMD<double>> hmelxi(elclass_inds.Size(), proxy->Dimension()*simd_ir.Size());
                      FlatMatrix<> hhmelxi(hmelxi.Height(), hmelxi.Width()*SIMD<double>::Size(), &hmelxi(0)[0]);
                      ParallelForRange (elclass_inds.Size(), [&] (IntRange myrange) {
                          hhmelxi.Rows(myrange) = melx.Rows(myrange) * hbmatx;
                        });
                      melxi.Append (std::move(hmelxi));
                      tx.AddFlops (elclass_inds.Size()*hbmatx.Height()*hbmatx.Width());
                    }
                }

                Array<Matrix<SIMD<double>>> mgfxi;
                {
                  RegionTimer reg(tgf);
                  int cntgf = 0;
                  for (CoefficientFunction * cf : gridfunction_cfs)
                    {
                      auto gfcf = dynamic_cast<GridFunctionCoefficientFunction*> (cf);
                      if (!gfcf)
                        {
                          cout << "no gf" << endl;
                          continue;
                        }
                      
                      auto fes = gfcf->GetGridFunction().GetFESpace();
                      auto & felgf = fes->GetFE(ei, lh);
                      auto diffop = gfcf->GetDifferentialOperator(trafo.VB());
                      
                      FlatMatrix<SIMD<double>> bmat(felgf.GetNDof()*fes->GetDimension()*
                                                    diffop->Dim(), // ??? right Dim
                                                    simd_ir.Size(), lh);
                      FlatMatrix<double> hbmat(felgf.GetNDof()*fes->GetDimension(),
                                               diffop->Dim()*    // right Dim ???
                                               simd_ir.Size()*SIMD<double>::Size(),
                                               &bmat(0)[0]);
                      bmat = SIMD<double> (0.0);
                      diffop->CalcMatrix(felgf, simd_mir, bmat);
                      
                      // Matrix<> mgf(elclass_inds.Size(), felgf.GetNDof()*fes->GetDimension());
                      Matrix<SIMD<double>> hmgfxi(elclass_inds.Size(), diffop->Dim()*simd_ir.Size());
                      FlatMatrix<> hhmgfxi(hmgfxi.Height(), hmgfxi.Width()*SIMD<double>::Size(), &hmgfxi(0)[0]);                      
                      // auto & vec = gfcf->GetGridFunction().GetVector();
                      
                      ParallelForRange
                        (elclass_inds.Size(), [&] (IntRange myrange) {
                          /*
                          Array<DofId> dofnr;
                          for (auto i : myrange)
                            {
                              fes->GetDofNrs( { VOL, elclass_inds[i] }, dofnr);
                              vec.GetIndirect(dofnr, mgf.Row(i));                    
                            }
                          */
                          // RegionTracer rt(TaskManager::GetThreadId(), tgfmult);                          
                          // hhmgfxi.Rows(myrange) = mgf.Rows(myrange) * hbmat;
                          hhmgfxi.Rows(myrange) = mgfs[cntgf].Rows(myrange) * hbmat;
                        });
                      
                      mgfxi.Append (move(hmgfxi));
                      cntgf++;
                    }
                }
                
                Array<Matrix<SIMD<double>>> melyi;
                for (auto proxy : test_proxies)
                  melyi.Append (Matrix<SIMD<double>>(elclass_inds.Size(),
                                                     proxy->Dimension()*simd_ir.Size()));

                {
                RegionTimer reg(teval);
                ParallelForRange
                  (elclass_inds.Size(), [&] (IntRange myrange)
                   {
                     LocalHeap llh = lh.Split();
                     ProxyUserData ud(trial_proxies.Size(), gridfunction_cfs.Size(), llh);
                     
                     auto & trafo = GetTrialSpace()->GetMeshAccess()->GetTrafo(ei, llh);
                     auto & simd_mir = trafo(simd_ir, llh);
                     /*
                     if (element_vb == BND)
                       {
                         simd_mir.ComputeNormalsAndMeasure (felx.ElementType(), facet);
                         cout << "need it ? " << endl;
                       }
                     */
                     const_cast<ElementTransformation&>(trafo).userdata = &ud;
                     
                     ud.fel = &felx;
                     // for (CoefficientFunction * cf : gridfunction_cfs)
                     // ud.AssignMemory (cf, simd_ir.GetNIP(), cf->Dimension(), llh);
                     
                     for (auto i : myrange)
                       {
                         for (int proxynr : Range(trial_proxies))
                           {
                             auto proxy = trial_proxies[proxynr];
                             ud.AssignMemory (proxy, FlatMatrix<SIMD<double>> (proxy->Dimension(), simd_ir.Size(), &melxi[proxynr](i,0)));
                           }

                         for (int cfnr : Range(gridfunction_cfs))
                           {
                             CoefficientFunction * cf = gridfunction_cfs[cfnr];
                             ud.AssignMemory (cf, FlatMatrix<SIMD<double>> (cf->Dimension(), simd_ir.Size(), &mgfxi[cfnr](i,0)));
                           }
                             
                         /*
                         int cfnr = 0;
                         for (CoefficientFunction * cf : gridfunction_cfs)
                           {
                             auto mat = ud.GetAMemory(cf);
                             ud.SetComputed(cf, true);                    
                             auto & bigmat = mgfxi[cfnr];
                             FlatVector<> sp(bigmat.Width(), &mat(0));
                             sp = bigmat.Row(i);
                             cfnr++;
                           }
                         */
                         
                         for (auto proxynr : Range(test_proxies))
                           {
                             auto proxy = test_proxies[proxynr];
                             FlatMatrix<SIMD<double>> simd_proxyvalues(proxy->Dimension(), simd_ir.Size(), &melyi[proxynr](i,0));
                             {
                               // RegionTracer rt(TaskManager::GetThreadId(), teval);
                               for (int k = 0; k < proxy->Dimension(); k++)
                                 {
                                   ud.testfunction = proxy;
                                   ud.test_comp = k;
                                   cf -> Evaluate (simd_mir, simd_proxyvalues.Rows(k,k+1));
                                 }
                             }
                           }
                       }
                   });
                }
                {
                  RegionTimer regy(ty);            
                  for (auto proxynr : Range(test_proxies))
                    {
                      auto proxy = test_proxies[proxynr];
                      FlatMatrix<SIMD<double>> bmaty(mely.Width()*proxy->Dimension(),
                                                     simd_ir.Size(), lh);
                      FlatMatrix<double> hbmaty(mely.Width(),
                                                proxy->Dimension()*simd_ir.Size()*SIMD<double>::Size(),
                                                &bmaty(0)[0]);

                      bmaty = SIMD<double> (0.0);                      
                      proxy->Evaluator()->CalcMatrix(fely, simd_mir, bmaty);
                      for (size_t i = 0; i < bmaty.Height(); i++)
                        {
                          auto row = bmaty.Row(i);
                          for (size_t j = 0; j < row.Size(); j++)
                            row(j) *= simd_mir[j].GetWeight(); 
                        }
                      
                      FlatMatrix<> hmely(melyi[proxynr].Height(),
                                         melyi[proxynr].Width()*SIMD<double>::Size(),
                                         &melyi[proxynr](0,0)[0]);
                      
                      ParallelForRange (elclass_inds.Size(), [&] (IntRange myrange) {
                          mely.Rows(myrange) += hmely.Rows(myrange) * Trans(hbmaty);
                        });
                    }
                }
              }
            /*
            FlatVector<SCAL> elys(fely.GetNDof()*fesy->GetDimension(), lh);
            for (auto i : Range(elclass_inds))
              {
                fesy->GetDofNrs({VOL, elclass_inds[i]}, dofnr);
                elys = val * mely.Row(i);
                y.AddIndirect(dofnr, elys, true);     // atomic add
              }
            */
            {
              RegionTimer reg(taddy);
              bool needs_atomic = fesy->ElementColoring().Size() > 1;
              ParallelForRange (elclass_inds.Size(), [&] (IntRange myrange) {
                  Array<DofId> dofnr;
                  Vector<SCAL> elys(fely.GetNDof()*fesy->GetDimension());
                  for (auto i : myrange)
                    {
                      fesy->GetDofNrs({VOL, elclass_inds[i]}, dofnr);
                      elys = val * mely.Row(i);
                      y.AddIndirect(dofnr, elys, needs_atomic);     // atomic add
                    }
                });
            }
          }
#endif



    for (auto elclass_inds : table)
      ParallelForRange (elclass_inds.Size(), [&] (IntRange myrange) {

          LocalHeap & clh = lh;
          LocalHeap lh = clh.Split();
          if (myrange.Size() == 0) return;

          auto myinds = elclass_inds.Range(myrange);
        
          ElementId ei(VOL, myinds[0]);
          auto & felx = GetTrialSpace()->GetFE (ei, lh);
          auto & fely = GetTestSpace()->GetFE (ei, lh);
          auto & trafo = GetTrialSpace()->GetMeshAccess()->GetTrafo(ei, lh);
          
          Matrix<> melx(myrange.Size(), felx.GetNDof()*fesx->GetDimension());
          Matrix<> mely(myrange.Size(), fely.GetNDof()*fesy->GetDimension());
          mely = 0.0;
          auto tid = TaskManager::GetThreadId();
          
          {
            RegionTracer rt(tid, tgetx);
            Array<DofId> dofnr;
            for (auto i : Range(myinds))
              {
                fesx->GetDofNrs( { VOL, myinds[i] }, dofnr);
                x.GetIndirect(dofnr, melx.Row(i));
              }
          }

          
          for (auto bfi1 : geom_free_parts)
            {
              auto bfi = dynamic_pointer_cast<SymbolicBilinearFormIntegrator> (bfi1);
              
              auto & trial_proxies = bfi->TrialProxies();
              auto & test_proxies = bfi->TestProxies();
              auto & gridfunction_cfs = bfi->GridFunctionCoefficients();
              auto & cf = bfi->GetCoefficientFunction();
              
              VorB element_vb = bfi->ElementVB();
              Facet2ElementTrafo f2el(felx.ElementType(), bfi->ElementVB());
              
              Array<Matrix<>> mgfs;
              {
                RegionTracer reg(tid, tgf);            
                for (CoefficientFunction * cf : gridfunction_cfs)
                  {
                    auto gfcf = dynamic_cast<GridFunctionCoefficientFunction*> (cf);
                    auto fes = gfcf->GetGridFunction().GetFESpace();
                    auto & felgf = fes->GetFE(ei, lh);
                    Matrix<> mgf(myinds.Size(), felgf.GetNDof()*fes->GetDimension());
                    auto & vec = gfcf->GetGridFunction().GetVector();
                    
                    Array<DofId> dofnr;
                    for (auto i : Range(myinds))
                      {
                        fes->GetDofNrs( { VOL, myinds[i] }, dofnr);
                        vec.GetIndirect(dofnr, mgf.Row(i));                    
                      }
                    mgfs.Append (move(mgf));
                  }
              }

              
              for (auto facet : Range(f2el.GetNFacets()))
                {
                  const SIMD_IntegrationRule & simd_ir = (element_vb == VOL) ? 
                    bfi->Get_SIMD_IntegrationRule (felx, lh)
                    :
                    f2el(facet,
                         bfi->GetSIMDIntegrationRule(f2el.FacetType (facet),
                                                     felx.Order()+fely.Order()+bfi->GetBonusIntegrationOrder()),
                         lh);
                  
                  auto & simd_mir = trafo(simd_ir, lh);
                  
                  Array<Matrix<SIMD<double>>> melxi;
                  {
                    RegionTracer regx(tid, tx);
                    for (auto proxynr : Range(trial_proxies))
                      {
                        auto proxy = trial_proxies[proxynr];
                        FlatMatrix<SIMD<double>> bmatx(melx.Width()*proxy->Dimension(),
                                                       simd_ir.Size(), lh);
                        FlatMatrix<double> hbmatx(melx.Width(),
                                                  proxy->Dimension()*simd_ir.Size()*SIMD<double>::Size(),
                                                  (double*)&bmatx(0));
                        
                        bmatx = SIMD<double> (0.0);
                        proxy->Evaluator()->CalcMatrix(felx, simd_mir, bmatx);
                        Matrix<SIMD<double>> hmelxi(myinds.Size(), proxy->Dimension()*simd_ir.Size());
                        FlatMatrix<> hhmelxi(hmelxi.Height(), hmelxi.Width()*SIMD<double>::Size(), (double*)&hmelxi(0));

                        hhmelxi = melx * hbmatx;

                        melxi.Append (std::move(hmelxi));
                      }
                  }

                  Array<Matrix<SIMD<double>>> mgfxi;
                  {
                    RegionTracer reg(tid, tgf);
                    int cntgf = 0;
                    for (CoefficientFunction * cf : gridfunction_cfs)
                      {
                        auto gfcf = dynamic_cast<GridFunctionCoefficientFunction*> (cf);
                      
                        auto fes = gfcf->GetGridFunction().GetFESpace();
                        auto & felgf = fes->GetFE(ei, lh);
                        auto diffop = gfcf->GetDifferentialOperator(trafo.VB());
                        
                        FlatMatrix<SIMD<double>> bmat(felgf.GetNDof()*fes->GetDimension()*
                                                      diffop->Dim(), // ??? right Dim
                                                      simd_ir.Size(), lh);
                        FlatMatrix<double> hbmat(felgf.GetNDof()*fes->GetDimension(),
                                                 diffop->Dim()*    // right Dim ???
                                                 simd_ir.Size()*SIMD<double>::Size(),
                                                 (double*)&bmat(0));
                        bmat = SIMD<double> (0.0);
                        diffop->CalcMatrix(felgf, simd_mir, bmat);
                        
                        Matrix<SIMD<double>> hmgfxi(myinds.Size(), diffop->Dim()*simd_ir.Size());
                        FlatMatrix<> hhmgfxi(hmgfxi.Height(), hmgfxi.Width()*SIMD<double>::Size(), (double*)&hmgfxi(0));      
                        
                        hhmgfxi = mgfs[cntgf] * hbmat;
                        
                        mgfxi.Append (move(hmgfxi));
                        cntgf++;
                      }
                  }
                
                  Array<Matrix<SIMD<double>>> melyi;
                  for (auto proxy : test_proxies)
                    melyi.Append (Matrix<SIMD<double>>(myinds.Size(),
                                                       proxy->Dimension()*simd_ir.Size()));
                  
                  {
                    RegionTracer reg(tid, teval);

                    ProxyUserData ud(trial_proxies.Size(), gridfunction_cfs.Size(), lh);
                     
                     auto & trafo = GetTrialSpace()->GetMeshAccess()->GetTrafo(ei, lh);
                     auto & simd_mir = trafo(simd_ir, lh);
                     const_cast<ElementTransformation&>(trafo).userdata = &ud;
                     
                     ud.fel = &felx;
                     
                     for (auto i : Range(myinds))
                       {
                         for (int proxynr : Range(trial_proxies))
                           {
                             auto proxy = trial_proxies[proxynr];
                             ud.AssignMemory (proxy, FlatMatrix<SIMD<double>> (proxy->Dimension(), simd_ir.Size(), &melxi[proxynr](i,0)));
                           }
                         
                         for (int cfnr : Range(gridfunction_cfs))
                           {
                             CoefficientFunction * cf = gridfunction_cfs[cfnr];
                             ud.AssignMemory (cf, FlatMatrix<SIMD<double>> (cf->Dimension(), simd_ir.Size(), &mgfxi[cfnr](i,0)));
                           }
                             
                         for (auto proxynr : Range(test_proxies))
                           {
                             auto proxy = test_proxies[proxynr];
                             FlatMatrix<SIMD<double>> simd_proxyvalues(proxy->Dimension(), simd_ir.Size(), &melyi[proxynr](i,0));
                             {
                               // RegionTracer rt(TaskManager::GetThreadId(), teval);
                               for (int k = 0; k < proxy->Dimension(); k++)
                                 {
                                   ud.testfunction = proxy;
                                   ud.test_comp = k;
                                   cf -> Evaluate (simd_mir, simd_proxyvalues.Rows(k,k+1));
                                 }
                             }
                           }
                       }
                  }
                  {
                    RegionTracer regy(tid, ty);            
                    for (auto proxynr : Range(test_proxies))
                      {
                        auto proxy = test_proxies[proxynr];
                        FlatMatrix<SIMD<double>> bmaty(mely.Width()*proxy->Dimension(),
                                                       simd_ir.Size(), lh);
                        FlatMatrix<double> hbmaty(mely.Width(),
                                                  proxy->Dimension()*simd_ir.Size()*SIMD<double>::Size(),
                                                  (double*)&bmaty(0));

                        bmaty = SIMD<double> (0.0);                      
                        proxy->Evaluator()->CalcMatrix(fely, simd_mir, bmaty);
                        for (size_t i = 0; i < bmaty.Height(); i++)
                          {
                            auto row = bmaty.Row(i);
                            for (size_t j = 0; j < row.Size(); j++)
                              row(j) *= simd_mir[j].GetWeight(); 
                          }
                        
                        FlatMatrix<> hmely(melyi[proxynr].Height(),
                                           melyi[proxynr].Width()*SIMD<double>::Size(),
                                           (double*)&melyi[proxynr](0,0));
                        
                        mely += hmely * Trans(hbmaty);
                      }
                  }
                } // for facet
            } // for integrator   
                  
          {
            RegionTracer reg(tid, taddy);
            bool needs_atomic = fesy->ElementColoring().Size() > 1;
            
            Array<DofId> dofnr;
            Vector<SCAL> elys(fely.GetNDof()*fesy->GetDimension());
            
            for (auto i : Range(myinds))
              {
                fesy->GetDofNrs({VOL, myinds[i]}, dofnr);
                elys = val * mely.Row(i);
                y.AddIndirect(dofnr, elys, needs_atomic);     // atomic add
              }
          }
        });

        



        
        /*
        Matrix<> temp_x(elclass_inds.Size(), !transpose ? elmat.Width() : elmat.Height());
        Matrix<> temp_y(elclass_inds.Size(), !transpose ? elmat.Height() : elmat.Width());

        ParallelForRange
          (Range(elclass_inds),
           [&] (IntRange myrange)
           {
             Array<DofId> dofs;
             
             int tid = TaskManager::GetThreadId();
             {
               RegionTimer r(tx;
               auto fvx = x.FVDouble();
               for (auto i : myrange)
                 {
                   fesx->GetDofNrs(ElementId(VOL,elclass_inds[i]), dofs);
                   x.GetIndirect(dofs, temp_x.Row(i));
                 }
             }
             
             {
               RegionTimer r(tm;
               RegionTracer rt(tid, tm);
               NgProfiler::AddThreadFlops(tm, tid, elmat.Height()*elmat.Width()*myrange.Size());
               if (!transpose)
                 temp_y.Rows(myrange) = temp_x.Rows(myrange) * Trans(elmat);
               else
                 temp_y.Rows(myrange) = temp_x.Rows(myrange) * elmat;
             }
             {
               RegionTimer r(ty;
               for (auto i : myrange)
                 {
                   fesy->GetDofNrs(ElementId(VOL,elclass_inds[i]), dofs);
                   y.AddIndirect(dofs, temp_y.Row(i), true);
                   // todo: check if atomic is actually needed for space
                 }
             }
           });
        */
      }
  
  
  


  template <class SCAL>
  void S_BilinearForm<SCAL> :: ApplyLinearizedMatrixAdd1 (SCAL val,
                                                          const BaseVector & lin,
                                                          const BaseVector & x,
                                                          BaseVector & y, LocalHeap & lh) const
  {
    if (!MixedSpaces())
      
      {
        Array<int> dnums;
      
        int ne = ma->GetNE();
        int dim = GetFESpace()->GetDimension(); 
        //LocalHeap lh (2000000, "biform-ApplyLinearized");

        bool hasbound = false;
        bool hasinner = false;

        for (int j = 0; j < NumIntegrators(); j++)
          {
            const BilinearFormIntegrator & bfi = *GetIntegrator(j);
            if (bfi.BoundaryForm())
              hasbound = true;
            else
              hasinner = true;
          }


        if (hasinner)
          for (int i = 0; i < ne; i++)
            {
              HeapReset hr(lh);
              ElementId ei(VOL, i);
              const FiniteElement & fel = fespace->GetFE (ei, lh);
              ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
              fespace->GetDofNrs (ei, dnums);
          
              FlatVector<SCAL> elveclin (dnums.Size() * dim, lh);
              FlatVector<SCAL> elvecx (dnums.Size() * dim, lh);
              FlatVector<SCAL> elvecy (dnums.Size() * dim, lh);

              lin.GetIndirect (dnums, elveclin);
              fespace->TransformVec (ei, elveclin, TRANSFORM_SOL);
              
              x.GetIndirect (dnums, elvecx);
              fespace->TransformVec (ei, elvecx, TRANSFORM_SOL);

              for (int j = 0; j < NumIntegrators(); j++)
                {
                  const BilinearFormIntegrator & bfi = *parts[j];

                  if (bfi.BoundaryForm()) continue;
                  if (!bfi.DefinedOn (ma->GetElIndex (ei))) continue;
                  if (!bfi.DefinedOnElement(ei.Nr())) continue;

                  bfi.ApplyLinearizedElementMatrix (fel, eltrans, elveclin, elvecx, elvecy, lh);

                  fespace->TransformVec (ei, elvecy, TRANSFORM_RHS);

                  elvecy *= val;

                  y.AddIndirect (dnums, elvecy);
                }
            }

        int nse = ma->GetNSE();
        if (hasbound)
          for (int i = 0; i < nse; i++)
            {
              HeapReset hr(lh);
              ElementId sei(BND, i);
              const FiniteElement & fel = fespace->GetFE (sei, lh);
              ElementTransformation & eltrans = ma->GetTrafo (sei, lh);
              fespace->GetDofNrs (sei, dnums);
            
              FlatVector<SCAL> elveclin (dnums.Size() * dim, lh);
              FlatVector<SCAL> elvecx (dnums.Size() * dim, lh);
              FlatVector<SCAL> elvecy (dnums.Size() * dim, lh);
            
              lin.GetIndirect (dnums, elveclin);
              fespace->TransformVec (sei, elveclin, TRANSFORM_SOL);
              x.GetIndirect (dnums, elvecx);
              fespace->TransformVec (sei, elvecx, TRANSFORM_SOL);
          
              for (int j = 0; j < NumIntegrators(); j++)
                {
                  const BilinearFormIntegrator & bfi = *parts[j];
                
                  if (!bfi.BoundaryForm()) continue;
                  if (!bfi.DefinedOn (eltrans.GetElementIndex())) continue;
                  if (!bfi.DefinedOnElement(sei.Nr())) continue;
              
                  bfi.ApplyLinearizedElementMatrix (fel, eltrans, elveclin, elvecx, elvecy, lh);
                  fespace->TransformVec (sei, elvecy, TRANSFORM_RHS);
                  elvecy *= val;
                  y.AddIndirect (dnums, elvecy);
                }
            }

        for (int i = 0; i < specialelements.Size(); i++)
          {
            HeapReset hr(lh);
            const SpecialElement & el = *specialelements[i];
            el.GetDofNrs (dnums);
          
            FlatVector<SCAL> elvecx (dnums.Size() * dim, lh);
            FlatVector<SCAL> elvecy (dnums.Size() * dim, lh);

            x.GetIndirect (dnums, elvecx);
            el.Apply (elvecx, elvecy, lh);

            elvecy *= val;
            y.AddIndirect (dnums, elvecy);
            lh.CleanUp();
          }
      }
    else
      { 
        cout << "apply not implemented for mixed" << endl;
      }
  }







  template <class SCAL>
  double S_BilinearForm<SCAL> :: Energy (const BaseVector & x, LocalHeap & lh) const
  {
    static Timer t("BilinearForm::Energy"); RegionTimer reg(t);
    
    double energy = 0.0;

    if (!MixedSpaces())
      {
        //LocalHeap lh (2000000, "biform-energy", true);

        for (auto vb : { VOL, BND, BBND, BBBND })
          if (VB_parts[vb].Size())
            IterateElements 
              (*fespace, vb, lh, 
               [&] (FESpace::Element ei, LocalHeap & lh)
               {
                 const FiniteElement & fel = fespace->GetFE (ei, lh);
                 ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
                 
                 FlatArray<int> dnums = ei.GetDofs();
                 FlatVector<SCAL> elvecx (dnums.Size()*GetFESpace()->GetDimension(), lh);
                 
                 x.GetIndirect (dnums, elvecx);
                 fespace->TransformVec (ei, elvecx, TRANSFORM_SOL);
                 
                 double energy_T = 0;
                 
                 for (auto & bfi : VB_parts[vb])
                   {
                     if (!bfi->DefinedOn (ei.GetIndex())) continue;
                     if (!bfi->DefinedOnElement(ei.Nr())) continue;
                     energy_T += bfi->Energy (fel, eltrans, elvecx, lh);
                   }
                 
                 AtomicAdd(energy, energy_T);
               });
        
        
        /*
        bool hasbound = false;
        bool hasinner = false;

        for (int j = 0; j < NumIntegrators(); j++)
          {
            const BilinearFormIntegrator & bfi = *GetIntegrator(j);
            if (bfi.VB()==BND)
              hasbound = true;
            else if(bfi.VB()==VOL)
              hasinner = true;
	    else
	      throw Exception("Energy not implemented for BBND objects yet!");
          }

        if (hasinner)
          IterateElements 
            (*fespace, VOL, lh, 
             [&] (FESpace::Element ei, LocalHeap & lh)
             {
               const FiniteElement & fel = fespace->GetFE (ei, lh);
               ElementTransformation & eltrans = ma->GetTrafo (ei, lh);

               FlatArray<int> dnums = ei.GetDofs();
               FlatVector<SCAL> elvecx (dnums.Size()*GetFESpace()->GetDimension(), lh);
               
               x.GetIndirect (dnums, elvecx);
               fespace->TransformVec (ei, elvecx, TRANSFORM_SOL);
               
               double energy_T = 0;

               for (auto bfi : parts)
                 {
                   if (bfi->BoundaryForm()) continue;
		   if (!bfi->DefinedOn (ei.GetIndex())) continue;
                   energy_T += bfi->Energy (fel, eltrans, elvecx, lh);
                 }

               energy += energy_T;
             });

        int nse = ma->GetNSE();
        Array<int> dnums;
        if (hasbound)
          for (int i = 0; i < nse; i++)
            {
              HeapReset hr(lh);
              ElementId sei(BND, i);
              const FiniteElement & fel = fespace->GetFE (sei, lh);
              ElementTransformation & eltrans = ma->GetTrafo (sei, lh);
              fespace->GetDofNrs (sei, dnums);
            
              FlatVector<SCAL> elvecx (dnums.Size() * GetFESpace()->GetDimension(), lh);
              x.GetIndirect (dnums, elvecx);
              fespace->TransformVec (sei, elvecx, TRANSFORM_SOL);
          
              for (int j = 0; j < NumIntegrators(); j++)
                {
                  const BilinearFormIntegrator & bfi = *parts[j];
		  if (!bfi.DefinedOn (ma->GetElIndex (sei))) continue;
                
                  if (!bfi.BoundaryForm()) continue;
                  energy += bfi.Energy (fel, eltrans, elvecx, lh);
                }
            }
        */
        Array<int> dnums;
        for (int i = 0; i < specialelements.Size(); i++)
          {
            HeapReset hr(lh);

            const SpecialElement & el = *specialelements[i];
            el.GetDofNrs (dnums);

            FlatVector<SCAL> elvecx (dnums.Size() * GetFESpace()->GetDimension(), lh);
            x.GetIndirect (dnums, elvecx);
          
            AtomicAdd(energy, el.Energy (elvecx, lh));
          }
      }
    return energy;
  }


  template <class SCAL>
  void S_BilinearForm<SCAL> :: 
  AddDiagElementMatrix (FlatArray<int> dnums1,
                        FlatVector<SCAL> diag,
                        bool inner_element, int elnr,
                        LocalHeap & lh)
  {
    throw Exception ("Baseclass::AddDiagElementMatrix");
  }
 

  template <class TSCAL>
  void S_BilinearForm<TSCAL>::LapackEigenSystem(FlatMatrix<TSCAL> & elmat, LocalHeap & lh) const 
  {
    if (!this->symmetric || this->fespace->IsComplex())
      {
        Vector<Complex> lami(elmat.Height());
        Matrix<TSCAL> evecs(elmat.Height());
        FlatMatrix<TSCAL> elmat_save(elmat.Height(), elmat.Width(), lh);
        elmat_save = elmat;
#ifdef LAPACK
        LapackEigenValues (elmat_save, lami, evecs);
        (*testout) << "lami = " 
                   << endl << lami << endl << "evecs: " << endl << evecs << endl;
#endif
      }
    else
      {
        Vector<TSCAL> lami(elmat.Height());
        Matrix<TSCAL> evecs(elmat.Height());
#ifdef LAPACK
        LapackEigenValuesSymmetric (elmat, lami, evecs);
#else
        CalcEigenSystem (elmat, lami, evecs);
#endif
        (*testout) << "lami = " 
                   << endl << lami << endl << "evecs: " << endl << evecs << endl;
      }
  }




  template class S_BilinearForm<double>;
  template class S_BilinearForm<Complex>;

  /*
  template <class TM, class TV>
  T_BilinearForm<TM,TV>::
  T_BilinearForm (shared_ptr<FESpace> afespace, const string & aname, const Flags & flags)
    : S_BilinearForm<TSCAL> (afespace, aname, flags)
  { }

  template <class TM, class TV>
  T_BilinearForm<TM,TV>::
  T_BilinearForm (shared_ptr<FESpace> afespace, 
                  shared_ptr<FESpace> afespace2,
                  const string & aname,
                  const Flags & flags)
    : S_BilinearForm<TSCAL> (afespace, afespace2, aname, flags) 
  { } 

  template <class TM, class TV>
  T_BilinearForm<TM,TV>:: ~T_BilinearForm () { }
  */

  
  template <class TM, class TV>
    shared_ptr<BilinearForm>  T_BilinearForm<TM,TV>::
    GetLowOrderBilinearForm()
  {
    if (this->low_order_bilinear_form)
      return this->low_order_bilinear_form;
    
    auto lospace = this->fespace->LowOrderFESpacePtr();
    if (!lospace) return nullptr;

    cout << IM(3) << "creating low order biform on demand" << endl;
    this->low_order_bilinear_form = 
      make_shared<T_BilinearForm<TM,TV>> (lospace, this->name+string(" low-order"), this->flags);

    for (auto igt : this->parts)
      this->low_order_bilinear_form->AddIntegrator(igt);      

    if (this->mats.Size())
      {
        LocalHeap lh(10*1000*1000);
        this->low_order_bilinear_form -> Assemble(lh);
      }
    return this->low_order_bilinear_form;
  }

  template <class TM, class TV>
  void T_BilinearForm<TM,TV>::
  AllocateMatrix ()
  {
    if (this->mats.Size() == this->ma->GetNLevels())
      return;

    MatrixGraph graph = this->GetGraph (this->ma->GetNLevels()-1, false);

    auto spmat = make_shared<SparseMatrix<TM,TV,TV>> (graph, 1);
    this->GetMemoryTracer().Track(*spmat, "mymatrix");
    mymatrix = spmat; // .get();
    
    if (this->spd) spmat->SetSPD();
    shared_ptr<BaseMatrix> mat = spmat;

    if (this->GetFESpace()->IsParallel())
      mat = make_shared<ParallelMatrix> (mat, this->GetTrialSpace()->GetParallelDofs(),
					 this->GetTestSpace()->GetParallelDofs());

    this->mats.SetSize(this->ma->GetNLevels());
    this->mats.Last() = mat;

    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size()-1; i++)
        this->mats[i].reset();

    this->AllocateInternalMatrices();
  }



  template <class TM, class TV>
  void T_BilinearForm<TM,TV>::
  CleanUpLevel ()
  {
    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size(); i++)
        this->mats[i].reset();
  }



  template <class TM, class TV>
  unique_ptr<BaseVector> T_BilinearForm<TM, TV>::
  CreateRowVector() const
  {
    auto afespace = this->GetTrialSpace();
    
    if (afespace->IsParallel())
      return make_unique<ParallelVVector<TV>> (afespace->GetParallelDofs());
    else
      return make_unique<VVector<TV>> (afespace->GetNDof());
    // return make_unique<S_BaseVectorPtr<TSCAL>> (afespace->GetNDof(), sizeof(TV)/sizeof(TSCAL));
  }
  
  template <class TM, class TV>
    unique_ptr<BaseVector> T_BilinearForm<TM, TV>::
  CreateColVector() const
  {
    auto afespace = this->GetTestSpace();
    if (afespace->IsParallel())
      return make_unique<ParallelVVector<TV>> (afespace->GetNDof(), afespace->GetParallelDofs());
    else
      return make_unique<VVector<TV>> (afespace->GetNDof());
  }



  /*
  template <class TM, class TS>
  inline void AddPartOfElementMatrix(TM & dest, const FlatMatrix<TS> & source,
                                     const int start1, const int start2)
  {
    int hi = Height(dest);
    int wi = Width(dest);
          
    for (int k = 0; k < hi; k++)
      for (int l = 0; l < wi; l++)
        dest(k,l) += source(start1*hi+k, start2*wi+l);
  }
  
  template <>
  inline void AddPartOfElementMatrix(double & dest, 
                                     const FlatMatrix<double> & source,
                                     const int start1, const int start2)
  {
    dest += source(start1, start2);
  }

  template <>
  inline void AddPartOfElementMatrix(Complex & dest, 
                                     const FlatMatrix<Complex> & source,
                                     const int start1, const int start2)
  {
    dest += source(start1, start2);
  }
  */
  
  
  template <class TM, class TV>
  void T_BilinearForm<TM,TV>::
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<TSCAL> elmat,
                    ElementId id, bool addatomic, 
                    LocalHeap & lh) 
  {
    mymatrix -> TMATRIX::AddElementMatrix (dnums1, dnums2, elmat, addatomic || this->fespace->HasAtomicDofs());
  }


  template <class TM, class TV>
  T_BilinearFormSymmetric<TM,TV> :: 
  T_BilinearFormSymmetric (shared_ptr<FESpace> afespace, const string & aname,
                           const Flags & flags)
    : S_BilinearForm<TSCAL> (afespace, aname, flags) 
  {
    /*
    if (this->fespace->LowOrderFESpacePtr())
      {
        this->low_order_bilinear_form = 
          make_shared<T_BilinearFormSymmetric<TM,TV>> 
          (this->fespace->LowOrderFESpacePtr(), aname+string(" low-order"), flags);
      }
    */
  }

  template <class TM, class TV>
  T_BilinearFormSymmetric<TM,TV> :: 
  ~T_BilinearFormSymmetric ()
  { ; }

  template <class TM, class TV>
    shared_ptr<BilinearForm>  T_BilinearFormSymmetric<TM,TV>::
    GetLowOrderBilinearForm()
  {
    if (this->low_order_bilinear_form)
      return this->low_order_bilinear_form;
    
    auto lospace = this->fespace->LowOrderFESpacePtr();
    if (!lospace) return nullptr;

    cout << IM(3) << "creating low order biform on demand" << endl;
    this->low_order_bilinear_form = 
      make_shared<T_BilinearFormSymmetric<TM,TV>> (lospace, this->name+string(" low-order"), this->flags);

    for (auto igt : this->parts)
      this->low_order_bilinear_form->AddIntegrator(igt);      

    if (this->mats.Size())
      {
        LocalHeap lh(10*1000*1000);
        this->low_order_bilinear_form -> Assemble(lh);
      }
    return this->low_order_bilinear_form;
  }


  
  template <class TM, class TV>
  void T_BilinearFormSymmetric<TM,TV> :: 
  AllocateMatrix ()
  {
    if (this->mats.Size() == this->ma->GetNLevels())
      return;

    MatrixGraph graph = this->GetGraph (this->ma->GetNLevels()-1, true);

    auto spmat = make_shared<SparseMatrixSymmetric<TM,TV>> (graph, 1);
    mymatrix = spmat; // .get();
    this->GetMemoryTracer().Track(*spmat, "mymatrix");
    
    if (this->spd) spmat->SetSPD();
    shared_ptr<BaseMatrix> mat = spmat;

    if (this->GetFESpace()->IsParallel())
      mat = make_shared<ParallelMatrix> (mat, this->GetTrialSpace()->GetParallelDofs(),
					 this->GetTestSpace()->GetParallelDofs());
    
    this->mats.Append (mat);

    // delete graph;

    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size()-1; i++)
        this->mats[i].reset();

    this->AllocateInternalMatrices();    
  }


  template <class TM, class TV>
  void T_BilinearFormSymmetric<TM,TV>::
  CleanUpLevel ()
  {
    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size(); i++)
        this->mats[i].reset();
  }


  template <class TM, class TV>
  unique_ptr<BaseVector> T_BilinearFormSymmetric<TM, TV>::
  CreateRowVector() const
  {
    auto fes = this->GetTrialSpace();

    if ( fes->IsParallel() )
      return make_unique<ParallelVVector<TV>> (fes->GetParallelDofs());
    else
      return make_unique<VVector<TV>> (fes->GetNDof());
  }

  template <class TM, class TV>
  unique_ptr<BaseVector> T_BilinearFormSymmetric<TM, TV>::
  CreateColVector() const
  {
    auto fes = this->GetTestSpace();

    if ( fes->IsParallel() )
      return make_unique<ParallelVVector<TV>> (fes->GetParallelDofs());
    else
      return make_unique<VVector<TV>> (fes->GetNDof());
  }
  


  template <class TM, class TV>
  void T_BilinearFormSymmetric<TM,TV> :: 
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<TSCAL> elmat,
                    ElementId id, bool addatomic,
                    LocalHeap & lh) 
  {
    mymatrix -> TMATRIX::AddElementMatrixSymmetric (dnums1, elmat, addatomic || this->fespace->HasAtomicDofs());
  }






  template <class TM, class TV>
  void T_BilinearFormSymmetric<TM,TV>::LapackEigenSystem(FlatMatrix<TSCAL> & elmat, LocalHeap & lh) const
  {
    if ( !this->fespace->IsComplex() )
      {
        Vector<TSCAL> lami(elmat.Height());
        Matrix<TSCAL> evecs(elmat.Height());
#ifdef LAPACK
        LapackEigenValuesSymmetric (elmat, lami, evecs);
#else
        CalcEigenSystem (elmat, lami, evecs);
#endif
        (*testout) << "lami = " << endl << lami << endl << "evecs: " << endl << evecs << endl;
      }
    else
      {
        Vector<Complex> lami(elmat.Height());
        Matrix<TSCAL> evecs(elmat.Height());
        FlatMatrix<TSCAL> elmat_save(elmat.Height(), elmat.Width(), lh);
        elmat_save = elmat;
#ifdef LAPACK
        LapackEigenValues (elmat_save, lami, evecs);
        (*testout) << "LAPACK NS for complex symmetric problem \nlami = " 
                   << endl << lami << endl << "evecs: " << endl << evecs << endl;
#endif
      }
  }



  
  /* *********************  T_BilinearFormDynBlocks ************** */


  template <typename TSCAL>
    shared_ptr<BilinearForm>  T_BilinearFormDynBlocks<TSCAL>::
    GetLowOrderBilinearForm()
  {
    if (this->low_order_bilinear_form)
      return this->low_order_bilinear_form;
    
    auto lospace = this->fespace->LowOrderFESpacePtr();
    if (!lospace) return nullptr;

    cout << IM(3) << "creating low order biform on demand" << endl;
    this->low_order_bilinear_form = 
      make_shared<T_BilinearFormDynBlocks<TSCAL>> (lospace, lospace, this->name+string(" low-order"), this->flags);

    for (auto igt : this->parts)
      this->low_order_bilinear_form->AddIntegrator(igt);      

    if (this->mats.Size())
      {
        LocalHeap lh(10*1000*1000);
        this->low_order_bilinear_form -> Assemble(lh);
      }
    return this->low_order_bilinear_form;
  }

  template <typename TSCAL>
  void T_BilinearFormDynBlocks<TSCAL>::
  AllocateMatrix ()
  {
    if (this->mats.Size() == this->ma->GetNLevels())
      return;

    MatrixGraph graph = this->GetGraph (this->ma->GetNLevels()-1, false);

    auto spmat = make_shared<SparseBlockMatrix<TSCAL>>
      (graph, blockheight, blockwidth, true);
    
    this->GetMemoryTracer().Track(*spmat, "mymatrix");
    mymatrix = spmat; // .get();
    
    if (this->spd) spmat->SetSPD();
    shared_ptr<BaseMatrix> mat = spmat;

    if (this->GetFESpace()->IsParallel())
      mat = make_shared<ParallelMatrix> (mat, this->GetTrialSpace()->GetParallelDofs(),
					 this->GetTestSpace()->GetParallelDofs());

    this->mats.SetSize(this->ma->GetNLevels());
    this->mats.Last() = mat;

    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size()-1; i++)
        this->mats[i].reset();

    this->AllocateInternalMatrices();
  }


  template <typename TSCAL>
  void T_BilinearFormDynBlocks<TSCAL>::
  CleanUpLevel ()
  {
    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size(); i++)
        this->mats[i].reset();
  }



  template <typename TSCAL>
  unique_ptr<BaseVector> T_BilinearFormDynBlocks<TSCAL>::
  CreateRowVector() const
  {
    auto afespace = this->GetTrialSpace();
    
    if (afespace->IsParallel())
      return make_unique<S_ParallelBaseVectorPtr<TSCAL>> (afespace->GetNDof(),
                                                          blockwidth,
                                                          afespace->GetParallelDofs(),
                                                          CUMULATED);
    else
      return make_unique<S_BaseVectorPtr<TSCAL>> (afespace->GetNDof(), blockwidth);
  }
  
  template <typename TSCAL>
    unique_ptr<BaseVector> T_BilinearFormDynBlocks<TSCAL>::
  CreateColVector() const
  {
    auto afespace = this->GetTestSpace();
    if (afespace->IsParallel())
      return make_unique<S_ParallelBaseVectorPtr<TSCAL>> (afespace->GetNDof(),
                                                         blockheight,
                                                         afespace->GetParallelDofs(),
                                                         DISTRIBUTED);
    else
      return make_unique<S_BaseVectorPtr<TSCAL>> (afespace->GetNDof(),
                                                  blockheight);
  }


  
  template <typename TSCAL>
  void T_BilinearFormDynBlocks<TSCAL>::
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<TSCAL> elmat,
                    ElementId id, bool addatomic, 
                    LocalHeap & lh) 
  {
    mymatrix -> TMATRIX::AddElementMatrix (dnums1, dnums2, elmat, addatomic || this->fespace->HasAtomicDofs());
  }







  


  template <class TM>
  T_BilinearFormDiagonal<TM> :: 
  T_BilinearFormDiagonal (shared_ptr<FESpace> afespace, const string & aname,
                          const Flags & flags)
    : S_BilinearForm<TSCAL> (afespace, aname, flags) 
  { 
    this->diagonal = 1;

    if (this->fespace->LowOrderFESpacePtr())
      {
        this->low_order_bilinear_form = 
          make_shared<T_BilinearFormSymmetric<TM>> 
          (this->fespace->LowOrderFESpacePtr(), aname+string(" low-order"), flags);
        this->low_order_bilinear_form -> SetDiagonal (0);
      }
  }

  template <class TM>
  T_BilinearFormDiagonal<TM> :: 
  ~T_BilinearFormDiagonal ()
  { ; }


  template <class TM>
  void T_BilinearFormDiagonal<TM> :: 
  AllocateMatrix ()
  {
    if (this->mats.Size() == this->ma->GetNLevels())
      return;

    size_t ndof = this->fespace->GetNDof();

    mymatrix = make_shared<DiagonalMatrix<TM>> (ndof);
    
    shared_ptr<BaseMatrix> mat = mymatrix;
    if (this->GetFESpace()->IsParallel())
      mat = make_shared<ParallelMatrix> (mat, this->GetTrialSpace()->GetParallelDofs(),
                                         this->GetTestSpace()->GetParallelDofs());
    this->mats.Append (mat);
    
    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size()-1; i++)
        this->mats[i].reset();
  }


  template <class TM>
  unique_ptr<BaseVector> T_BilinearFormDiagonal<TM> :: 
  CreateRowVector() const
  {
    auto fes = this->GetTrialSpace();
    if ( fes->IsParallel() )
      return make_unique<ParallelVVector<TV_COL>> (fes->GetParallelDofs());
    else
      return make_unique<VVector<TV_COL>> (fes->GetNDof());
  }
  
  template <class TM>
  unique_ptr<BaseVector> T_BilinearFormDiagonal<TM> :: 
  CreateColVector() const
  {
    auto fes = this->GetTestSpace();
    if ( fes->IsParallel() )
      return make_unique<ParallelVVector<TV_COL>> (fes->GetParallelDofs());
    else
      return make_unique<VVector<TV_COL>> (fes->GetNDof());
  }


  ///
  template <class TM>
  void T_BilinearFormDiagonal<TM> :: 
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<TSCAL> elmat,
                    ElementId id, bool addatomic, 
                    LocalHeap & lh) 
  {
    if (addatomic) throw Exception ("atomic add for DiagonalMatrix not implemented");
    
    for (int i = 0; i < dnums1.Size(); i++)
      if (IsRegularDof(dnums1[i]))
        {
          TM & mij = (*mymatrix)(dnums1[i]); // , dnums1[i]);
          int hi = Height (mij);
          int wi = Width (mij);

          for (int k = 0; k < hi; k++)
            for (int l = 0; l < wi; l++)
              mij(k,l) += elmat(i*hi+k, i*wi+l);
        }
  }



  ///
  template <> void T_BilinearFormDiagonal<double>::
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<double> elmat,
                    ElementId id, bool addatomic,
                    LocalHeap & lh) 
  {
    if (addatomic) throw Exception ("atomic add for DiagonalMatrix not implemented");
    
    for (int i = 0; i < dnums1.Size(); i++)
      if (IsRegularDof(dnums1[i]))
        (*mymatrix)(dnums1[i]) += elmat(i, i);
  }




  ///
  template <> void T_BilinearFormDiagonal<Complex>::
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<Complex> elmat,
                    ElementId id, bool addatomic,
                    LocalHeap & lh) 
  {
    if (addatomic) throw Exception ("atomic add for DiagonalMatrix not implemented");
    
    for (int i = 0; i < dnums1.Size(); i++)
      if (IsRegularDof(dnums1[i]))
        (*mymatrix)(dnums1[i]) += elmat(i, i);
  }


  template <class TM>
  void T_BilinearFormDiagonal<TM> :: 
  AddDiagElementMatrix (FlatArray<int> dnums1,
                        FlatVector<TSCAL> diag,
                        bool inner_element, int elnr,
                        LocalHeap & lh) 
  {
    throw Exception ("generic AddDiagElementMatrix not implemented");
    /*
      TMATRIX & mat = dynamic_cast<TMATRIX&> (*this->mats.Last());
    
      for (int i = 0; i < dnums1.Size(); i++)
      if (dnums1[i] != -1)
      {
      TM mij = diag(dnums1[i]);
      int hi = Height (mij);
      int wi = Width (mij);
          
      for (int k = 0; k < hi; k++)
      for (int l = 0; l < wi; l++)
      mij(k,l) += elmat(i*hi+k, i*wi+l);
      }
    */
  }



  ///
  template <> void T_BilinearFormDiagonal<double>::
    AddDiagElementMatrix (FlatArray<int> dnums1,
                          FlatVector<double> diag,
                          bool inner_element, int elnr,
                          LocalHeap & lh) 
  {
    for (int i = 0; i < dnums1.Size(); i++)
      if (IsRegularDof(dnums1[i]))
        (*mymatrix)(dnums1[i]) += diag(i);
  }

  ///
  template <> void T_BilinearFormDiagonal<Complex>::
    AddDiagElementMatrix (FlatArray<int> dnums1,
                          FlatVector<Complex> diag,
                          bool inner_element, int elnr,
                          LocalHeap & lh) 
  {
    for (int i = 0; i < dnums1.Size(); i++)
      if (IsRegularDof(dnums1[i]))
        (*mymatrix)(dnums1[i]) += diag(i);
  }














  ComponentBilinearForm :: ComponentBilinearForm (shared_ptr<BilinearForm> abase_blf, int acomp, int ancomp)
    : BilinearForm( (*dynamic_pointer_cast<CompoundFESpace> (abase_blf->GetFESpace()))[acomp], "comp-lf", Flags()), 
      base_blf(abase_blf), comp(acomp) // , ncomp(ancomp)
  { 
    ;
  }

  BilinearForm & ComponentBilinearForm :: AddIntegrator (shared_ptr<BilinearFormIntegrator> bfi)
  {
    auto block_bfi = make_shared<CompoundBilinearFormIntegrator> (bfi, comp);
    block_bfi->SetDefinedOn (bfi->GetDefinedOn());
    base_blf -> AddIntegrator (block_bfi);
    return *this;
  }








  template <int CBSIZE>
  shared_ptr<BilinearForm> CreateBilinearForm1 (int cb_size,
                                                shared_ptr<FESpace> space, const string & name, const Flags & flags)
  {
    if (CBSIZE == cb_size)
      return make_shared<T_BilinearFormSymmetric<double, Vec<CBSIZE,Complex>>> (space, name, flags);
    else
      return shared_ptr<BilinearForm> (CreateBilinearForm1<CBSIZE-1> (cb_size, space, name, flags));
  }

  template <> 
  shared_ptr<BilinearForm> CreateBilinearForm1<1> (int cb_size,
                                         shared_ptr<FESpace> space, const string & name, const Flags & flags)
  {
    return make_shared<T_BilinearFormSymmetric<double, Complex>> (space, name, flags);
  }

  template <> 
  shared_ptr<BilinearForm> CreateBilinearForm1<0> (int cb_size,
                                                   shared_ptr<FESpace> space, const string & name, const Flags & flags)
  {
    throw Exception ("Illegal cacheblocksize" + ToString (cb_size));
  }
  

  shared_ptr<BilinearForm> CreateBilinearForm (shared_ptr<FESpace> space,
                                               const string & name,
                                               const Flags & flags)
  {
    BilinearForm * bf = NULL;

    if (flags.GetDefineFlag ("ebe")){
      if ( space->IsComplex() )
        return make_shared<ElementByElement_BilinearForm<Complex>> (space, name, flags);
      else 
        return make_shared<ElementByElement_BilinearForm<double>> (space, name, flags);
    }

    if (flags.GetDefineFlag ("nonassemble"))
      {
        if ( space->IsComplex() )
          return make_shared<S_BilinearFormNonAssemble<Complex>> (space, name, flags);
        else 
          return make_shared<S_BilinearFormNonAssemble<double>> (space, name, flags);
      }

    
    // bool symmetric_storage = flags.GetDefineFlag ("symmetric") || flags.GetDefineFlag ("spd");
    // if (flags.GetDefineFlag("nonsym") || flags.GetDefineFlag("nonsym_storage")) symmetric_storage = false;
    bool symmetric_storage = false;
    if (flags.GetDefineFlagX("nonsym_storage").IsFalse() || flags.GetDefineFlagX("symmetric_storage").IsTrue())
      symmetric_storage = true;
    
    
    if (symmetric_storage)
      {

        if (space->IsComplex() && flags.GetDefineFlag ("real"))
          {
            if(flags.NumFlagDefined("cacheblocksize"))
              {
                int cacheblocksize = int(flags.GetNumFlag("cacheblocksize", 1)); 
                /*
#if MAX_CACHEBLOCKS >= 2 
#ifdef GOLD
                if (cacheblocksize == 2)
                  return new T_BilinearFormSymmetric<double,MD<2>> (*space, name, flags);
#endif
#endif
                */
                return CreateBilinearForm1<MAX_CACHEBLOCKS> (cacheblocksize, space, name, flags);
                /*
                  switch(int(flags.GetNumFlag("cacheblocksize",1)))
                  {
                  #if MAX_CACHEBLOCKS >= 2
                  case 2:
                  return new T_BilinearFormSymmetric<double,Vec<2,Complex> > (*space, name, flags);
                  #endif
                  #if MAX_CACHEBLOCKS >= 3
                  case 3:
                  return new T_BilinearFormSymmetric<double,Vec<3,Complex> > (*space, name, flags);
                  case 4:
                  return new T_BilinearFormSymmetric<double,Vec<4,Complex> > (*space, name, flags);
                  #endif
                  #if MAX_CACHEBLOCKS >= 5
                  case 5:
                  return new T_BilinearFormSymmetric<double,Vec<5,Complex> > (*space, name, flags);
                  case 6:
                  return new T_BilinearFormSymmetric<double,Vec<6,Complex> > (*space, name, flags);
                  case 7:
                  return new T_BilinearFormSymmetric<double,Vec<7,Complex> > (*space, name, flags);
                  case 8:
                  return new T_BilinearFormSymmetric<double,Vec<8,Complex> > (*space, name, flags);
                  case 9:
                  return new T_BilinearFormSymmetric<double,Vec<9,Complex> > (*space, name, flags);
                  case 10:
                  return new T_BilinearFormSymmetric<double,Vec<10,Complex> > (*space, name, flags);
                  case 11:
                  return new T_BilinearFormSymmetric<double,Vec<11,Complex> > (*space, name, flags);
                  case 12:
                  return new T_BilinearFormSymmetric<double,Vec<12,Complex> > (*space, name, flags);
                  case 13:
                  return new T_BilinearFormSymmetric<double,Vec<13,Complex> > (*space, name, flags);
                  case 14:
                  return new T_BilinearFormSymmetric<double,Vec<14,Complex> > (*space, name, flags);
                  case 15:
                  return new T_BilinearFormSymmetric<double,Vec<15,Complex> > (*space, name, flags);
                  #endif
                  }
                */
              }
            else
              return make_shared<T_BilinearFormSymmetric<double,Complex>> (space, name, flags);
          }

        if(flags.NumFlagDefined("cacheblocksize"))
          {
            CreateSymMatObject4 (bf, T_BilinearFormSymmetric, 
                                 space->GetDimension(),
                                 int(flags.GetNumFlag("cacheblocksize",1)),
                                 space->IsComplex(),   
                                 space, name, flags);
          }
        else
          {
            //bf = CreateSymMatObject<T_BilinearFormSymmetric, BilinearForm> //, const FESpace, const string, const Flags>
            //  (space->GetDimension(), space->IsComplex(), *space, name, flags);

            CreateSymMatObject3(bf, T_BilinearFormSymmetric,
                                space->GetDimension(), space->IsComplex(),
                                space, name, flags);
          }
      }
    else if (flags.GetDefineFlag ("diagonal"))
      {
	/*
        bf = CreateSymMatObject<T_BilinearFormDiagonal, BilinearForm> //, const FESpace, const string, const Flags>
          (space->GetDimension(), space->IsComplex(), *space, name, flags);
	*/
        
        CreateSymMatObject3 (bf, T_BilinearFormDiagonal, 
                             space->GetDimension(), space->IsComplex(),   
                             space, name, flags);
      }
    else
      {
        if (space->IsComplex() && flags.GetDefineFlag ("real"))
          {
            
            if(flags.NumFlagDefined("cacheblocksize"))
              {
                switch(int(flags.GetNumFlag("cacheblocksize",1)))
                  {
#if MAX_CACHEBLOCKS >= 2
                  case 2:
                    return new T_BilinearForm<double,Vec<2,Complex> > (*space, name, flags);
#endif
#if MAX_CACHEBLOCKS >= 3
                  case 3:
                    return new T_BilinearForm<double,Vec<3,Complex> > (*space, name, flags);
                  case 4:
                    return new T_BilinearForm<double,Vec<4,Complex> > (*space, name, flags);
#endif
#if MAX_CACHEBLOCKS >= 5
                  case 5:
                    return new T_BilinearForm<double,Vec<5,Complex> > (*space, name, flags);
                  case 6:
                    return new T_BilinearForm<double,Vec<6,Complex> > (*space, name, flags);
                  case 7:
                    return new T_BilinearForm<double,Vec<7,Complex> > (*space, name, flags);
                  case 8:
                    return new T_BilinearForm<double,Vec<8,Complex> > (*space, name, flags);
                  case 9:
                    return new T_BilinearForm<double,Vec<9,Complex> > (*space, name, flags);
                  case 10:
                    return new T_BilinearForm<double,Vec<10,Complex> > (*space, name, flags);
                  case 11:
                    return new T_BilinearForm<double,Vec<11,Complex> > (*space, name, flags);
                  case 12:
                    return new T_BilinearForm<double,Vec<12,Complex> > (*space, name, flags);
                  case 13:
                    return new T_BilinearForm<double,Vec<13,Complex> > (*space, name, flags);
                  case 14:
                    return new T_BilinearForm<double,Vec<14,Complex> > (*space, name, flags);
                  case 15:
                    return new T_BilinearForm<double,Vec<15,Complex> > (*space, name, flags);
#endif
                  }
              }
            else
              return make_shared<T_BilinearForm<double,Complex>> (space, name, flags);
          }
        
        if(flags.NumFlagDefined("cacheblocksize"))
          {
            CreateSymMatObject4 (bf, T_BilinearForm, 
                                 space->GetDimension(),
                                 int(flags.GetNumFlag("cacheblocksize",1)),
                                 space->IsComplex(),   
                                 space, name, flags);
          }
        else
          {
            //bf = CreateSymMatObject<T_BilinearForm, BilinearForm> 
            //  (space->GetDimension(), space->IsComplex(), *space, name, flags);

            if (space->GetDimension() > MAX_SYS_DIM)
              {
                return make_shared<T_BilinearFormDynBlocks<double>> (space, name, flags);
              }
            else
              CreateSymMatObject3 (bf, T_BilinearForm,
                                   space->GetDimension(), space->IsComplex(),
                                   space, name, flags);
          }
      }

    if (!bf)
      throw Exception (string ("Could not create BilinearForm, space-dimension is ")+
                       ToString (space->GetDimension()) +
                       "\neither define MAX_SYS_DIM with higher value and recompile,"
                       "\nor set flag 'nonassemble'");
    return shared_ptr<BilinearForm> (bf);
  }



 


  void BilinearForm :: GalerkinProjection ()
  {
    auto prol = fespace->GetProlongation();

    if ( !low_order_bilinear_form )
      for (int finelevel = ma->GetNLevels()-1; finelevel>0; finelevel--)
        {
          auto prolMat = prol->CreateProlongationMatrix (finelevel);
          
          if (prolMat)					  
            mats[finelevel-1] = dynamic_cast< const BaseSparseMatrix& >(GetMatrix(finelevel)).
              Restrict(*prolMat,dynamic_pointer_cast<BaseSparseMatrix>(GetMatrixPtr(finelevel-1)));
        }
  }

  
  BilinearFormApplication :: 
  BilinearFormApplication (shared_ptr<BilinearForm> abf,
                           LocalHeap & alh)
    : bf (abf), lh(alh)
  {
    ;
  }

  void BilinearFormApplication :: 
  Mult (const BaseVector & v, BaseVector & prod) const
  {
    static Timer t("BilinearFormApplication"); RegionTimer r(t);
    v.Cumulate();

    prod = 0;
    bf -> AddMatrix (1, v, prod, lh);
    
    prod.SetParallelStatus (DISTRIBUTED);
  }

  void BilinearFormApplication :: 
  MultAdd (double val, const BaseVector & v, BaseVector & prod) const
  {
    v.Cumulate();
    prod.Distribute();

    bf -> AddMatrix (val, v, prod, lh);
  }

  void BilinearFormApplication :: 
  MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const
  {
    v.Cumulate();
    prod.Distribute();

    bf -> AddMatrix (val, v, prod, lh);
  }

  void BilinearFormApplication :: 
  MultTransAdd (double val, const BaseVector & v, BaseVector & prod) const
  {
    v.Cumulate();
    prod.Distribute();

    bf -> AddMatrixTrans (val, v, prod, lh);
  }

  shared_ptr<BilinearForm> CreateBilinearForm (shared_ptr<FESpace> space,
                                               shared_ptr<FESpace> space2,
                                               const string & name,
                                               const Flags & flags)
  {
    if (flags.GetDefineFlag ("nonassemble"))
      {
        if ( space->IsComplex() )
          return make_shared<S_BilinearFormNonAssemble<Complex>> (space, space2, name, flags);
        else 
          return make_shared<S_BilinearFormNonAssemble<double>> (space, space2, name, flags);
      }
    else
      {
        if ( space->IsComplex() )
          return make_shared<T_BilinearForm<Complex>> (space, space2, name, flags);
        else 
          return make_shared<T_BilinearForm<double>> (space, space2, name, flags);
      }
    // throw Exception ("cannot craeate mixes-space without nonassemble - flag");
  }
  
  /*
    void BilinearFormApplication :: 
    MultTransAdd (double val, const BaseVector & v, BaseVector & prod) const 
    {
    bf -> AddMatrix (val, v, prod);
    }
  */

  AutoVector BilinearFormApplication :: 
  CreateVector () const
  {
    return bf -> CreateRowVector();
  }

  AutoVector BilinearFormApplication :: 
  CreateRowVector () const
  {
    return bf -> CreateRowVector();
  }
  
  AutoVector BilinearFormApplication :: 
  CreateColVector () const
  {
    return bf -> CreateColVector();
  }
  
  
  LinearizedBilinearFormApplication ::
  LinearizedBilinearFormApplication (shared_ptr<BilinearForm> abf,
                                     const BaseVector * aveclin,
                                     LocalHeap & alh)
    : BilinearFormApplication (abf, alh), veclin(aveclin)
  {
    ;
  }

  void  LinearizedBilinearFormApplication :: 
  Mult (const BaseVector & v, BaseVector & prod) const
  {
    prod = 0;
    bf->ApplyLinearizedMatrixAdd (1, *veclin, v, prod, lh);
  }

  void LinearizedBilinearFormApplication :: 
  MultAdd (double val, const BaseVector & v, BaseVector & prod) const
  {
    bf->ApplyLinearizedMatrixAdd (val, *veclin, v, prod, lh);
  }

  void LinearizedBilinearFormApplication :: 
  MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const
  {
    bf->ApplyLinearizedMatrixAdd (val, *veclin, v, prod, lh);
  }

  template <class SCAL>
  S_BilinearFormNonAssemble<SCAL> :: 
  S_BilinearFormNonAssemble (shared_ptr<FESpace> afespace, const string & aname,
                                 const Flags & flags)
    : S_BilinearForm<SCAL> (afespace, aname, flags)
  { ; }

  template <class SCAL>
  S_BilinearFormNonAssemble<SCAL> :: 
  S_BilinearFormNonAssemble (shared_ptr<FESpace> afespace, shared_ptr<FESpace> afespace2,
                             const string & aname, const Flags & flags)
    : S_BilinearForm<SCAL> (afespace, afespace2, aname, flags)
  { ; }

  template <class SCAL>
  unique_ptr<BaseVector> S_BilinearFormNonAssemble<SCAL> :: 
  CreateRowVector() const
  {
    auto afespace = this->GetTrialSpace();
    if (afespace->IsParallel())
      return make_unique<S_ParallelBaseVectorPtr<SCAL>> (afespace->GetNDof(), afespace->GetDimension(), afespace->GetParallelDofs(), DISTRIBUTED);
    return make_unique<S_BaseVectorPtr<SCAL>> (afespace->GetNDof(), afespace->GetDimension());
  }
  
  template <class SCAL>
  unique_ptr<BaseVector> S_BilinearFormNonAssemble<SCAL> :: 
  CreateColVector() const
  {
    auto afespace = this->GetTestSpace();
    if (afespace->IsParallel())
      return make_unique<S_ParallelBaseVectorPtr<SCAL>> (afespace->GetNDof(), afespace->GetDimension(), afespace->GetParallelDofs(), DISTRIBUTED);
    return make_unique<S_BaseVectorPtr<SCAL>> (afespace->GetNDof(), afespace->GetDimension());
  }
  
  
  template <class SCAL>
  ElementByElement_BilinearForm<SCAL> :: 
  ElementByElement_BilinearForm (shared_ptr<FESpace> afespace, const string & aname,
                                 const Flags & flags)
    : S_BilinearForm<SCAL> (afespace, aname, flags)
  { ; }

  template <class SCAL>
  ElementByElement_BilinearForm<SCAL> :: ~ElementByElement_BilinearForm ()
  { ; }



  
  template <class SCAL>
  void ElementByElement_BilinearForm<SCAL> :: AllocateMatrix ()
  {
    auto fespace = this->fespace;
    this->mats.Append (make_shared<ElementByElementMatrix<SCAL>> (fespace->GetNDof(), this->ma->GetNE()+this->ma->GetNSE() ));
  }


  template<class SCAL>
    unique_ptr<BaseVector> ElementByElement_BilinearForm<SCAL> :: CreateRowVector() const
  {
    return make_unique<VVector<SCAL>> (this->GetTrialSpace()->GetNDof());
  }
  template<class SCAL>
    unique_ptr<BaseVector> ElementByElement_BilinearForm<SCAL> :: CreateColVector() const
  {
    return make_unique<VVector<SCAL>> (this->GetTestSpace()->GetNDof());
  }

  template<class SCAL>
  void ElementByElement_BilinearForm<SCAL> :: 
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<SCAL> elmat,
                    ElementId id, bool addatomic,
                    LocalHeap & lh)
  {
    int nr = id.Nr();
    if (id.IsBoundary()) nr += this->ma->GetNE();

    if (addatomic) throw Exception ("atomic add for EBE Matrix not implemented");    
    dynamic_cast<ElementByElementMatrix<SCAL>&> (this->GetMatrix()).AddElementMatrix (nr, dnums1, dnums2, elmat);
  }
  

  template class ElementByElement_BilinearForm<double>;
  template class ElementByElement_BilinearForm<Complex>;

  template class T_BilinearForm<double,double>;
  template class T_BilinearFormSymmetric<double,double>;

  template class T_BilinearForm<Complex,Complex>;
}
