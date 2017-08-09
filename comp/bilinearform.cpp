#include <comp.hpp>
#include <multigrid.hpp>

#include <parallelngs.hpp>

namespace ngcomp
{
  // dummy function header 
  void CalcEigenSystem (FlatMatrix<Complex> & elmat, 
                        FlatVector<Complex> & lami, 
                        FlatMatrix<Complex> & evecs)
  { ; }


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
    SetEliminateInternal (flags.GetDefineFlag ("eliminate_internal"));
    SetKeepInternal (eliminate_internal && 
                     !flags.GetDefineFlag ("nokeep_internal"));
    SetStoreInner (flags.GetDefineFlag ("store_inner"));
    precompute = flags.GetDefineFlag ("precompute");
    checksum = flags.GetDefineFlag ("checksum");
    spd = flags.GetDefineFlag ("spd");
    if (spd) symmetric = true;
  }


  BilinearForm :: 
  BilinearForm (shared_ptr<FESpace> afespace,
                shared_ptr<FESpace> afespace2, const string & aname,
                const Flags & flags)
    : NGS_Object(afespace->GetMeshAccess(), flags, aname), fespace(afespace), fespace2(afespace2)
  {
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
    if (flags.GetDefineFlag ("eliminate_internal")) SetEliminateInternal (1);
    SetKeepInternal (eliminate_internal && 
                     !flags.GetDefineFlag ("nokeep_internal"));
    if (flags.GetDefineFlag ("store_inner")) SetStoreInner (1);

    precompute = flags.GetDefineFlag ("precompute");
    checksum = flags.GetDefineFlag ("checksum");
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
    bfi = FixDimension(bfi, fespace->GetSpacialDimension());
    
    if (symmetric && !bfi->IsSymmetric())
      throw Exception (string ("Adding non-symmetric integrator to symmetric bilinear-form\n")+
                       string ("bfi is ")+bfi->Name());

    parts.Append (bfi);

    if(bfi->SkeletonForm())
      {
        auto dgform = bfi -> GetDGFormulation();
        if (dgform.element_boundary) {
          elementwise_skeleton_parts.Append(bfi);
#ifdef PARALLEL
	  auto fbfi = dynamic_pointer_cast<FacetBilinearFormIntegrator> (bfi);
	  if (!fbfi)  throw Exception ("not a FacetBFI");
	  mpi_facet_parts.Append(fbfi);
#endif
	}
	else
          {
            if (bfi->VB() > 1) throw Exception ("skeletonform makes sense only for VOL or BND");
            // facetwise_skeleton_parts[bfi->VB()].Append(bfi);
            auto fbfi = dynamic_pointer_cast<FacetBilinearFormIntegrator> (bfi);
            if (!fbfi)  throw Exception ("not a FacetBFI");
            facetwise_skeleton_parts[bfi->VB()] += fbfi;
#ifdef PARALLEL
	    if(bfi->VB()==VOL) { //BND-integrators are per definition not mpi!!
	      mpi_facet_parts.Append(fbfi);
	    }
#endif
          }
      }
    else
      VB_parts[bfi->VB()].Append(bfi);
    
    if (low_order_bilinear_form)
      low_order_bilinear_form -> AddIntegrator (parts.Last());

    return *this;
  }


  /*
    void BilinearForm :: AddIndependentIntegrator (BilinearFormIntegrator * bfi,
    const int master_surface,
    const int slave,
    const bool deletable)
    {
    independent_parts.Append (bfi);
    independent_parts_deletable.Append(deletable);
    Vec<2,int> indices(master_surface,slave); 
    independent_meshindex.Append(indices);
    if (low_order_bilinear_form)
    low_order_bilinear_form -> AddIndependentIntegrator (independent_parts.Last(),
    master_surface,slave,
    deletable);
    }
  */

  BilinearForm :: ~BilinearForm ()
  {
    ; // delete low_order_bilinear_form;
  }

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

  void BilinearForm :: SetPreconditioner (Preconditioner * pre)
  {
    // cout << "SetPreconditioner, type fes = " << typeid(*fespace).name() << ", type pre = " << typeid(*pre).name() << endl;
    if (preconditioners.Contains(pre))
      throw Exception (string("preconditioner ")+typeid(*pre).name()+ " already registered in bfa");
    if (pre->GetFlags().GetDefineFlag("not_register_for_auto_update"))
      throw Exception (string("'not_register_for_auto_update' set, but preconditioner") + typeid(*pre).name() +" registers anyway");
    preconditioners.Append (pre);
  }


  MatrixGraph * BilinearForm :: GetGraph (int level, bool symmetric)
  {
    static Timer timer ("BilinearForm::GetGraph");
    RegionTimer reg (timer);

    int ndof = fespace->GetNDof();
    int nf = ma->GetNFacets();
    int neV = ma->GetNE(VOL);
    int neB = ma->GetNE(BND);
    int neBB = ma->GetNE(BBND);
    const Array<SpecialElement*> & specialelements = fespace->GetSpecialElements();
    int nspe = specialelements.Size();

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
	    int nre = ma->GetNE(vb);
	    ParallelForRange (Range(nre), [&](IntRange r)
			      {
				Array<DofId> dnums;
				for (auto i : r)
				  {
				    auto eid = ElementId(vb,i);
				    if (!fespace->DefinedOn (vb,ma->GetElIndex(eid))) continue;
				    
				    if (vb == VOL && eliminate_internal)
				      fespace->GetDofNrs (i, dnums, EXTERNAL_DOF);
				    else
				      fespace->GetDofNrs (eid, dnums);
				    int shift = (vb==VOL) ? 0 : ((vb==BND) ? neV : neV+neB);
				    for (int d : dnums)
				      if (d != -1) creator.Add (shift+i, d);
                              }
			      });
	  }

        
        for (int i = 0; i < specialelements.Size(); i++)
          {
            specialelements[i]->GetDofNrs (dnums);
            for (int d : dnums)
              if (d != -1) creator.Add (neV+neB+neBB+i, d);
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
                if (dnums_dg[j] != -1 && (j==0 || (dnums_dg[j] != dnums_dg[j-1]) ))
                  creator.Add (neV+neB+neBB+nspe+i, dnums_dg[j]);
            }
        }

      }
    
    MatrixGraph * graph;

    if (!fespace2)
      {
        graph = new MatrixGraph (ndof, *creator.GetTable(), *creator.GetTable(), symmetric);
        delete creator.GetTable();
      }
    else
      {
        TableCreator<int> creator2(maxind);
        for ( ; !creator2.Done(); creator2++)
          {
	    for(VorB vb  : {VOL, BND, BBND})
	      {
		int nre = ma->GetNE(vb);
		for (int i = 0; i < nre; i++)
		  {
		    auto eid = ElementId(vb,i);
		    if (!fespace2->DefinedOn (vb,ma->GetElIndex(eid))) continue;
		    
		    if (vb == VOL && eliminate_internal)
		      fespace2->GetDofNrs (i, dnums, EXTERNAL_DOF);
		    else
		      fespace2->GetDofNrs (eid, dnums);
		    
		    int shift = (vb==VOL) ? 0 : ((vb==BND) ? neV : neV + neB);
		    for (int d : dnums)
		      if (d != -1) creator2.Add (shift+i, d);
		  }
		
		/*
		// just not tested ...
		for (int i = 0; i < specialelements.Size(); i++)
		{
                specialelements[i]->GetDofNrs (dnums);
		
                for (int j = 0; j < dnums.Size(); j++)
		if (dnums[j] != -1)
		creator2.Add (ne+nse+i, dnums[j]);
		}
		
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
		if (!fespace2->DefinedOn (ma->GetElIndex(elnr))) continue;
		fespace2->GetDofNrs (elnr, dnums);
		for (int j = 0; j < dnums.Size(); j++)
		if (dnums[j] != -1)
		creator2.Add (ne+nse+nspe+i, dnums[j]);
		}
                }
		*/
	      }
	  }

        graph = new MatrixGraph (fespace2->GetNDof(), *creator2.GetTable(), *creator.GetTable(), symmetric);        
        delete creator.GetTable();
        delete creator2.GetTable();
      }

    graph -> FindSameNZE();
    return graph;
  }






  void BilinearForm :: Assemble (LocalHeap & lh)
  {
    if (mats.Size() == ma->GetNLevels())
      return;


    if (nonassemble)
      {
        mats.Append (make_shared<BilinearFormApplication> (shared_ptr<BilinearForm>(this, NOOP_Deleter))); 
      
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

            auto vecf = mats.Last()->CreateVector();
            auto vecu = mats.Last()->CreateVector();
          
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
      
        auto vecf = mats.Last()->CreateVector();
        auto vecu = mats.Last()->CreateVector();
      
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

        int nze = dynamic_cast<const BaseSparseMatrix &> (*mats.Last()) . NZE();
        cout << "NZE = " << nze << ", MFLOP = " << double (nze * steps) / time * 1e-6 << endl;
        cout << "type = " << typeid(*mats.Last()).name() << endl;
      }

    if (galerkin)
      GalerkinProjection();
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



  void BilinearForm :: PrintReport (ostream & ost) const
  {
    ost << "on space " << GetFESpace()->GetName() << endl
        << "symmetric   = " << symmetric << endl
        << "multilevel  = " << multilevel << endl
        << "nonassemble = " << nonassemble << endl
        << "printelmat = " << printelmat << endl
        << "elmatev    = " << elmat_ev << endl
        << "eliminate_internal = " << eliminate_internal << endl
        << "keep_internal = " << keep_internal << endl
        << "store_inner = " << store_inner << endl
        << "integrators: " << endl;
  
    for (int i = 0; i < parts.Size(); i++)
      ost << "  " << parts[i]->Name() << endl;
  }


  void BilinearForm :: MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    if (low_order_bilinear_form)
      low_order_bilinear_form -> MemoryUsage (mu);
  
    int olds = mu.Size();

    for (int i = 0; i < mats.Size(); i++)
      if (mats[i]) mats[i]->MemoryUsage (mu);

    for (int i = olds; i < mu.Size(); i++)
      mu[i]->AddName (string(" bf ")+GetName());
  }


  template <class SCAL>
  S_BilinearForm<SCAL> :: ~S_BilinearForm() { ; }




  template <class SCAL>
  void S_BilinearForm<SCAL> :: DoAssemble (LocalHeap & clh)
  {
    static Timer mattimer("Matrix assembling");

    static Timer mattimer_checkintegrators("Matrix assembling check integrators");
    static Timer mattimer1a("Matrix assembling initializing");
    static Timer mattimer_finalize("Matrix assembling finalize matrix");
    static Timer mattimer_VB[] = { Timer("Matrix assembling vol"),
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
		 if(bfi->DefinedOn(el.GetIndex()))
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
            
	    if (fespace->specialelements.Size())
	      loopsteps += fespace->specialelements.Size(); 
	    
            size_t gcnt = 0; //global count (for all cases)
            
            for (auto pre : preconditioners)
              pre -> InitLevel(fespace->GetFreeDofs());

	    mattimer1a.Stop();

	    for (VorB vb : { VOL, BND, BBND })
	      {
		if (!VB_parts[vb].Size()) continue;
                
                RegionTimer reg(mattimer_VB[vb]);
		size_t ne = ma->GetNE(vb);
                
                if (vb==VOL && diagonal)
                  {
                    double prevtime = WallTime();
                    Array<DofId> dnums;

                    for (int i = 0; i < ne; i++)
                      {
                        HeapReset hr(clh);
                        ElementId id(vb, i);
                        gcnt++;
			
                        if (clock()-prevtime > 0.1 * CLOCKS_PER_SEC)
                          {
                            cout << IM(3) << "\rassemble element " << i << "/" << ne << flush;
                            ma->SetThreadPercentage ( 100.0*gcnt / (loopsteps) );
                            prevtime = clock();
                          }
                        
                        // clh.CleanUp(heapp);
			
                        if (!fespace->DefinedOn (vb,ma->GetElIndex (id))) continue;
			
                        const FiniteElement & fel = fespace->GetFE (id, clh);
                        ElementTransformation & eltrans = ma->GetTrafo (id, clh);
                        auto ei = ElementId(vb,i);
                        fespace->GetDofNrs (ei, dnums);
                        
                        FlatVector<SCAL> sum_diag(dnums.Size()*fespace->GetDimension(), clh);
                        sum_diag = 0;
                        
                        for (auto & bfip : VB_parts[vb])
                          {
                            const BilinearFormIntegrator & bfi = *bfip;
                            if (!bfi.DefinedOn (ma->GetElIndex (id))) continue;
                            if (!bfi.DefinedOnElement (i)) continue;
                            
                            FlatVector<double> diag;
                            try
                              {
                                bfi.CalcElementMatrixDiag (fel, eltrans, diag, clh);
				
                                if (printelmat) 
                                  {
                                    testout->precision(8);
                                    (*testout) << "elnum= " << ElementId(vb,i) << endl;
                                    (*testout) << "integrator " << bfi.Name() << endl;
                                    (*testout) << "dnums = " << endl << dnums << endl;
                                    (*testout) << "diag-elmat = " << endl << diag << endl;
                                  }
                              }
                            catch (Exception & e)
                              {
                                e.Append (string("in Assemble Element Matrix, bfi = ") + 
                                          bfi.Name() + string("\n"));
                                throw;
                              }
                            catch (exception & e)
                              {
                                throw (Exception (string(e.what()) +
                                                  string("in Assemble Element Matrix, bfi = ") + 
                                                  bfi.Name() + string("\n")));
                              }
                            
                            sum_diag += diag;
                          }
			
                        AddDiagElementMatrix (dnums, sum_diag, 1, i, clh);

                        if (check_unused)
                          for (int j = 0; j < dnums.Size(); j++)
                            if (dnums[j] != -1)
                              useddof[dnums[j]] = true;
                      }
                    cout << IM(3) << "\rassemble element " << ne << "/" << ne << endl;
                  }
                else // not diagonal
                  {
                    ProgressOutput progress(ma,string("assemble ") + ToString(vb) + string(" element"), ma->GetNE(vb));
                    if (vb == VOL && eliminate_internal && keep_internal)
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
                             cout << "fel:GetNDof() = " << fel.GetNDof() << endl;
                             cout << "dnums.Size() = " << dnums.Size() << endl;
                             
                             *testout << "Info from finite element: " << endl;
                             fel.Print (*testout);
                             (*testout) << "fel:GetNDof() = " << fel.GetNDof() << endl;
                             (*testout) << "dnums.Size() = " << dnums.Size() << endl;
                             (*testout) << "dnums = " << dnums << endl;
                             throw Exception ( "Inconsistent number of degrees of freedom " );
                           }
                         
                         int elmat_size = dnums.Size()*fespace->GetDimension();
                         FlatMatrix<SCAL> sum_elmat(elmat_size, lh);
                         sum_elmat = 0;
			 bool elem_has_integrator = false;
                         for (auto & bfip : VB_parts[vb])
                           {
                             const BilinearFormIntegrator & bfi = *bfip;
                             if (!bfi.DefinedOn (el.GetIndex())) continue;                        
                             if (!bfi.DefinedOnElement (el.Nr())) continue;                        

                             elem_has_integrator = true;

                             if (printelmat || elmat_ev)
                               {
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
                                           if (d == -1) *testout << "0 ";
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
                             else
                               {
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
                           }
                         
                         
                         if (!elem_has_integrator) return;
                         
                         fespace->TransformMat (el, sum_elmat, TRANSFORM_MAT_LEFT_RIGHT);
			 
                         if (elmat_ev)
                           {
                             (*testout) << "sum matrix:" << endl;
                             LapackEigenSystem(sum_elmat, lh);
                           }
                         
                         if (vb == VOL && eliminate_internal)
                           {
                             // static Timer statcondtimer("static condensation", 1);
                             // RegionTimer regstat (statcondtimer);
                             
                             Array<int> idofs1(dnums.Size(), lh);
                             
                             fespace->GetDofNrs (el.Nr(), idofs1, LOCAL_DOF);
                             for (int j = 0; j < idofs1.Size(); j++)
                               idofs1[j] = dnums.Pos(idofs1[j]);
                             
                             if (printelmat) 
                               {
                                 lock_guard<mutex> guard(printelmat_mutex);
                                 *testout << "eliminate internal" << endl;
                                 *testout << "idofs1 = " << idofs1 << endl;
                               }
                             
                             if (idofs1.Size())
                               {
                                 HeapReset hr (lh);
				 
                                 int size = sum_elmat.Height();
                                 int dim = size / dnums.Size();
				 
                                 int sizei = dim * idofs1.Size();
                                 int sizeo = size - sizei;
				 
                                 FlatArray<int> idofs (sizei, lh);
                                 FlatArray<int> odofs (sizeo, lh);
				     
                                 for (int j = 0, k = 0; j < idofs1.Size(); j++)
                                   for (int jj = 0; jj < dim; jj++)
                                     idofs[k++] = dim*idofs1[j]+jj;
                                 
                                 for (int j = 0, k = 0; j < size; j++)
                                   if (!idofs.Contains(j))
                                     odofs[k++] = j;
                                 
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
                                 if (!keep_internal) 
                                   {
                                     LapackAInvBt (d, b);    // b <--- b d^-1
                                     LapackMultAddABt (b, c, -1, a);                                 
                                   }
                                 else
                                   {
                                     Array<int> idnums1(dnums.Size(), lh), 
                                       ednums1(dnums.Size(), lh);
                                     fespace->GetDofNrs(el.Nr(),idnums1,LOCAL_DOF);
                                     fespace->GetDofNrs(el.Nr(),ednums1,EXTERNAL_DOF);
                                     
                                     Array<int> idnums(dim*idnums1.Size(), lh);
                                     Array<int> ednums(dim*ednums1.Size(), lh);
                                     idnums.SetSize0(); 
                                     ednums.SetSize0();
                                     for (size_t j = 0; j < idnums1.Size(); j++)
                                       idnums += dim*IntRange(idnums1[j], idnums1[j]+1);
                                     for (size_t j = 0; j < ednums1.Size(); j++)
                                       ednums += dim * IntRange(ednums1[j], ednums1[j]+1);
                                     
                                     if (store_inner)
                                       innermatrix ->AddElementMatrix(el.Nr(),idnums,idnums,d);
                                     
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
                                     
                                     
                                     LapackInverse (d);
                                     FlatMatrix<SCAL> he (sizei, sizeo, lh);
                                     he = 0.0;
                                     he -= d * Trans(c) | Lapack;
                                     harmonicext ->AddElementMatrix(el.Nr(),idnums,ednums,he);
                                     if (!symmetric)
                                       {
                                         FlatMatrix<SCAL> het (sizeo, sizei, lh);
                                         het = 0.0;
                                         LapackMultAddAB (b, d, -1, het);
                                         //  het = -1.0 * b * d;
                                         static_cast<ElementByElementMatrix<SCAL>*>(harmonicexttrans.get())
                                           ->AddElementMatrix(el.Nr(),ednums,idnums,het);
                                       }
                                     
                                     innersolve ->AddElementMatrix(el.Nr(),idnums,idnums,d);
                                     a += b * he | Lapack;
                                     
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
                                   dnums[idofs1[k]] = -1;
                               }
                           }
                         if (printelmat)
                           {
                             lock_guard<mutex> guard(printelmat_mutex);
                             *testout<< "elem " << el << ", elmat = " << endl << sum_elmat << endl;
                           }
                         
                         AddElementMatrix (dnums, dnums, sum_elmat, el, lh);
			 
                         for (auto pre : preconditioners)
                           pre -> AddElementMatrix (dnums, sum_elmat, el, lh);

                         if (check_unused)
                           {
                             if (printelmat)
                               *testout << "set these as useddof: " << dnums << endl;
                             for (auto d : dnums)
                               if (d != -1) useddof[d] = true;
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
              { // inner facets
                // auto vb = VOL;
                size_t ne = ma->GetNE(VOL);
                BitArray fine_facet(nf);
                fine_facet.Clear();
                Array<int> elfacets;
                
                for (int i = 0; i < ne; ++i)
                  {
                    auto elfacets = ma->GetElFacets(ElementId(VOL,i));
                    for (auto f : elfacets) fine_facet.Set(f);
                  }
                
                int cnt = 0;
                // better: use facet coloring
                ParallelForRange
                  ( IntRange(nf), [&] ( IntRange r )
                    {
                      LocalHeap lh = clh.Split();
                      
                      Array<int> elnums_per(2, lh);

                      Array<int> dnums, dnums1, dnums2, elnums, fnums, vnums1, vnums2;
                      for (int i : r)
                        {
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
                          
                          {
                            lock_guard<mutex> guard(printmatasstatus2_mutex);
                            cnt++;
                            gcnt++;
                            if (cnt % 10 == 0)
                              cout << "\rassemble inner facet element " << cnt << "/" << nf << flush;
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
                          for (auto & bfi : facetwise_skeleton_parts[VOL])
                            {
                              // if (!bfi->SkeletonForm()) continue;
                              // if (bfi->VB() == BND) continue;
                              if (!bfi->DefinedOn (ma->GetElIndex (ei1))) continue;
                              if (!bfi->DefinedOn (ma->GetElIndex (ei2))) continue;
                              if (!bfi->DefinedOnElement(i)) continue;                        
                              
                              /*
                                for (int k = 0; k < dnums.Size(); k++)
                                if (dnums[k] != -1)
                                useddof[dnums[k]] = true;
                              */
                              if (check_unused)
                                for (auto d : dnums)
                                  if (d != -1) useddof[d] = true;
                              
                              int elmat_size = (dnums1.Size()+dnums2.Size())*fespace->GetDimension();
                              FlatMatrix<SCAL> elmat(elmat_size, lh);
                              
                              shared_ptr<FacetBilinearFormIntegrator> fbfi = 
                                dynamic_pointer_cast<FacetBilinearFormIntegrator>(bfi);
                              
                              if (fbfi)
                                {
                                  fbfi->CalcFacetMatrix (fel1,facnr1,eltrans1,vnums1,
                                                         fel2,facnr2,eltrans2,vnums2, elmat, lh);
                                }
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
                                  fbfi->CalcFacetMatrix (cfel1[comp], facnr1, eltrans1, vnums1,
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
                              
                              // *testout << "elmat : \n" << elmat << endl;
                              
                              fespace->TransformMat (ei1, elmat.Rows(0,dnums1.Size()), TRANSFORM_MAT_LEFT);
                              fespace->TransformMat (ei2, elmat.Rows(dnums1.Size(),dnums2.Size()), TRANSFORM_MAT_LEFT);
                              fespace->TransformMat (ei1, elmat.Cols(0,dnums1.Size()), TRANSFORM_MAT_RIGHT);
                              fespace->TransformMat (ei2, elmat.Cols(dnums1.Size(),dnums2.Size()), TRANSFORM_MAT_RIGHT);
                              
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
                              
                              {
                                lock_guard<mutex> guard(addelemfacin_mutex);
                                AddElementMatrix (compressed_dnums, compressed_dnums, compressed_elmat, ElementId(BND,i), lh);
                              }
                            }
                        }
                    });
                cout << "\rassemble inner facet element " << nf << "/" << nf << endl;
                
              }

            
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
                                     cout << "Surface fel:GetNDof() = " << fel.GetNDof() << endl;
                                     cout << "dnums.Size() = " << dnums.Size() << endl;
                                     
                                     (*testout) << "fel:GetNDof() = " << fel.GetNDof() << endl;
                                     (*testout) << "dnums.Size() = " << dnums.Size() << endl;
                                     (*testout) << "dnums = " << dnums << endl;
                                     throw Exception ( "Inconsistent number of degrees of freedom " );
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

                                     if (check_unused)
                                       for (auto d : dnums)
                                         if (d != -1) useddof[d] = true;
                                     
                                     int elmat_size = dnums.Size()*fespace->GetDimension();
                                     FlatMatrix<SCAL> elmat(elmat_size, lh);
                                     
                                     dynamic_cast<const FacetBilinearFormIntegrator&>(*bfi).  
                                       CalcFacetMatrix (fel,facnr1,eltrans,vnums1, seltrans, vnums2, elmat, lh);
                                     
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
                                       lock_guard<mutex> guard(addelemfacbnd_mutex);
                                       AddElementMatrix (dnums, dnums, elmat, ElementId(VOL,el1), lh);
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
                                 cout << "\rassemble inner facet element " << cnt << "/" << nf << flush;
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
                                     if (d != -1) useddof[d] = true;
                                 
                                 int elmat_size = (dnums1.Size()+dnums2.Size())*fespace->GetDimension();
                                 FlatMatrix<SCAL> elmat(elmat_size, lh);
                                 
                                 shared_ptr<FacetBilinearFormIntegrator> fbfi = 
                                   dynamic_pointer_cast<FacetBilinearFormIntegrator>(bfi);
                                 
                                 if (fbfi)
                                   {
                                     fbfi->CalcFacetMatrix (fel1,facnr1,eltrans1,vnums1,
                                                            fel2,facnr2,eltrans2,vnums2, elmat, lh);
                                   }
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
                                 
                                 // *testout << "elmat : \n" << elmat << endl;
                                 
                                 fespace->TransformMat (ei1, elmat.Rows(0,dnums1.Size()), TRANSFORM_MAT_LEFT);
                                 fespace->TransformMat (ei2, elmat.Rows(dnums1.Size(),dnums2.Size()), TRANSFORM_MAT_LEFT);
                                 fespace->TransformMat (ei1, elmat.Cols(0,dnums1.Size()), TRANSFORM_MAT_RIGHT);
                                 fespace->TransformMat (ei2, elmat.Cols(dnums1.Size(),dnums2.Size()), TRANSFORM_MAT_RIGHT);
                                 
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
                                   lock_guard<mutex> guard(addelemfacin_mutex);
                                   AddElementMatrix (compressed_dnums, compressed_dnums, compressed_elmat,
                                                     ElementId(BND,el1), lh);
                                 }
                               }
                           }
                       }                             
                   });
                cout << "\rassemble inner facet element " << nf << "/" << nf << endl;
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
                              
                          if (!fespace->DefinedOn (BND,ma->GetElIndex (sei))) continue;
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
                              cout << "Surface fel:GetNDof() = " << fel.GetNDof() << endl;
                              cout << "dnums.Size() = " << dnums.Size() << endl;
				  
                              (*testout) << "fel:GetNDof() = " << fel.GetNDof() << endl;
                              (*testout) << "dnums.Size() = " << dnums.Size() << endl;
                              (*testout) << "dnums = " << dnums << endl;
                              throw Exception ( "Inconsistent number of degrees of freedom " );
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
                                  if (dnums[k] != -1)
                                    useddof[dnums[k]] = true;
                                  
                              int elmat_size = dnums.Size()*fespace->GetDimension();
                              FlatMatrix<SCAL> elmat(elmat_size, lh);
				  
                              // original version did not compile on MacOS V
                              const FacetBilinearFormIntegrator & fbfi = 
                                dynamic_cast<const FacetBilinearFormIntegrator&>(*bfi);  
                              fbfi.CalcFacetMatrix (fel,facnr,eltrans,vnums, seltrans, svnums, elmat, lh);
				  
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
                                lock_guard<mutex> guard(addelemfacbnd_mutex);
                                AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), lh);
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
            
            
            int nspecel = 0;
            ParallelForRange( IntRange(fespace->specialelements.Size()), [&] ( IntRange r )
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
                      cout << "\rassemble special element " << nspecel << "/" << fespace->specialelements.Size() << flush;
                    ma->SetThreadPercentage ( 100.0*(gcnt) / (loopsteps) );
                  }
                  
                  const SpecialElement & el = *fespace->specialelements[i];
                
                  el.GetDofNrs (dnums);
                
                  FlatMatrix<SCAL> elmat(dnums.Size(), lh);
                  el.Assemble (elmat, lh);
                  
                  {
                    lock_guard<mutex> guard(printmatspecel2_mutex);
                    if (check_unused)
                      for (int j = 0; j < dnums.Size(); j++)
                        if (dnums[j] != -1)
                          useddof[dnums[j]] = true;
                    
                    AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), lh);
                  }
                  
                  assembledspecialelements = true;
                  lh.CleanUp();
                }
            });
            if(assembledspecialelements) cout << "\rassemble special element " 
                                              << fespace->specialelements.Size() << "/" << fespace->specialelements.Size() << endl;

            
            
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
                    AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), clh);
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
                      AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), clh);
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
            
            int MASK = eliminate_internal ? EXTERNAL_DOF : ANY_DOF;
            bool first_time = true;

            if (MyMPI_GetNTasks() == 1 && check_unused)
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

        else

          {
            // mixed spaces

            cout << "assemble mixed bilinearform" << endl;
      
            BaseMatrix & mat = GetMatrix();
            mat = 0.0;


            bool hasbound = false;
            
            if (VB_parts[VOL].Size())
              IterateElements 
                (*fespace, VOL, clh,          // coloring for 1 space is enough
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
                   for (auto & bfi : VB_parts[VOL])
                     {
                       MixedFiniteElement fel(fel1, fel2);
                       bfi->CalcElementMatrix (fel, eltrans, elmat, lh);
                       /*
                        fespace->Transform (i, true, elmat, TRANSFORM_MAT_RIGHT);
                        fespace2->Transform (i, true, elmat, TRANFORM_MAT_LEFT);
                       */
                       AddElementMatrix (dnums2, dnums1, elmat, ei, lh);
                     }
                 });


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

                       /*
                        fespace->Transform (i, true, elmat, TRANSFORM_MAT_RIGHT);
                        fespace2->Transform (i, true, elmat, TRANFORM_MAT_LEFT);
                       */

                       AddElementMatrix (dnums2, dnums1, elmat, ei, lh);
                     }
                 });

            if (print) *testout << "mat = " << mat << endl;
            
            cout << endl;
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
                    fespace->GetDofNrs (i, dnums, LOCAL_DOF);            
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
                     fespace->GetDofNrs (ei.Nr(), idofs, LOCAL_DOF);
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
    // static Timer timerVB[] = { timervol, timerbound, timerbbound };

    static mutex addelmatboundary1_mutex;

    RegionTimer reg (timer);

    if (this->mats.Size() < this->ma->GetNLevels())
      AllocateMatrix();

    // timestamp = ++global_timestamp;
    timestamp = GetNextTimeStamp();

    try
      {
        int ndof = fespace->GetNDof();
        Array<bool> useddof(ndof);
        useddof = false;
      
        BaseMatrix & mat = GetMatrix();
        mat = 0.0;
      
        cout << IM(3) << "Assemble linearization" << endl;
      
        Array<int> dnums;

        for (auto pre : preconditioners)
          pre -> InitLevel(fespace->GetFreeDofs());

        if (VB_parts[VOL].Size())
          {
            RegionTimer reg(timervol);
            ProgressOutput progress (ma, "assemble element", ma->GetNE());

            if (eliminate_internal && keep_internal)
              {
                size_t ndof = fespace->GetNDof();
                size_t ne = ma->GetNE();
                harmonicext = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
                if (!symmetric)
                  harmonicexttrans = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
                else
                  harmonicexttrans = make_shared<Transpose>(*harmonicext);
                innersolve = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
                if (store_inner)
                  innermatrix = make_shared<ElementByElementMatrix<SCAL>>(ndof, ne);
              }
            
            IterateElements 
              (*fespace, VOL, clh,  [&] (FESpace::Element el, LocalHeap & lh)
               {
                 const FiniteElement & fel = fespace->GetFE (el, lh);
                 ElementTransformation & eltrans = ma->GetTrafo (el, lh);
                 
                 Array<int> dnums(fel.GetNDof(), lh);
                 fespace->GetDofNrs (el, dnums);
                 
                 if(fel.GetNDof() != dnums.Size())
                   {
                     cout << "fel::GetNDof() = " << fel.GetNDof() << endl;
                     cout << "dnums.Size() = " << dnums.Size() << endl;
                   }

                 for (auto d : dnums) 
                   if (d != -1) useddof[d] = true;
                 
                 FlatMatrix<SCAL> sum_elmat(dnums.Size()*fespace->GetDimension(), lh);
                 FlatMatrix<SCAL> elmat(dnums.Size()*fespace->GetDimension(), lh);
                 sum_elmat = 0;
                 
                 FlatVector<SCAL> elveclin (dnums.Size()*fespace->GetDimension(), lh);
                 lin.GetIndirect (dnums, elveclin);
                 fespace->TransformVec (el, elveclin, TRANSFORM_SOL);

                 for (auto & bfi : VB_parts[VOL])
                   {
                     if (!bfi->DefinedOn (el.GetIndex())) continue;
                     if (!bfi->DefinedOnElement (el.Nr())) continue;
                     
                     try
                       {
                         bfi->CalcLinearizedElementMatrix (fel, eltrans, elveclin, elmat, lh);
                         
                         if (printelmat) 
                           {
                             testout->precision(8);
                             (*testout) << "elnum= " << el.Nr() << endl;
                             (*testout) << "eltype " << fel.ElementType() << endl;
                             (*testout) << "integrator " << bfi->Name() << endl;
                             (*testout) << "dnums = " << endl << dnums << endl;
                             (*testout) << "elveclin = " << endl << elveclin << endl;
                             (*testout) << "elmat = " << endl << elmat << endl;
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

                 
                 fespace->TransformMat (el, sum_elmat, TRANSFORM_MAT_LEFT_RIGHT);
                 
                 if (printelmat)
                   *testout << "summat = " << sum_elmat << endl;
                 

                 if (eliminate_internal)
                   {
                     static Timer statcondtimer("static condensation", 2);
                     RegionTimer regstat (statcondtimer);
                     
                     ArrayMem<int,100> idofs, idofs1, odofs;
                     int i = el.Nr();

                     fespace->GetDofNrs (i, idofs1, LOCAL_DOF);
                     for (int j = 0; j < idofs1.Size(); j++)
                       idofs1[j] = dnums.Pos(idofs1[j]);
                          
                     if (printelmat) 
                       {
                         *testout << "eliminate internal" << endl;
                         *testout << "idofs1 = " << idofs1 << endl;
                       }
                     
                     
                     if (idofs1.Size())
                       {
                         HeapReset hr (lh);
                         
                         int size = sum_elmat.Height();
                         int dim = size / dnums.Size();
                         
                         idofs.SetSize (0);
                         for (int j = 0; j < idofs1.Size(); j++)
                           for (int jj = 0; jj < dim; jj++)
                             idofs.Append (dim*idofs1[j]+jj);
                         
                         int sizei = idofs.Size();
                         int sizeo = size - sizei;
                              
                         odofs.SetSize (0);
                         for (int j = 0; j < size; j++)
                           if (!idofs.Contains(j))
                             odofs.Append(j);
                         
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
                         if (!keep_internal)
                           {
                             LapackAInvBt (d, b);
                             LapackMultAddABt (b, c, -1, a);
                           }
                         else
                           {
                             ArrayMem<int,50> idnums1, idnums;
                             ArrayMem<int,50> ednums1, ednums;
                             
                             fespace->GetDofNrs(i,idnums1,LOCAL_DOF);
                             fespace->GetDofNrs(i,ednums1,EXTERNAL_DOF);
                             
                             for (int j = 0; j < idnums1.Size(); j++)
                               idnums += dim*IntRange(idnums1[j], idnums1[j]+1);
                             for (int j = 0; j < ednums1.Size(); j++)
                               ednums += dim * IntRange(ednums1[j], ednums1[j]+1);
                             
                             if (store_inner)
                               static_cast<ElementByElementMatrix<SCAL>*>(innermatrix.get())
                                      ->AddElementMatrix(i,idnums,idnums,d);
                             
                             LapackInverse (d);
                             FlatMatrix<SCAL> he (sizei, sizeo, lh);
                             he=0.0;
                             LapackMultAddABt (d, c, -1, he);
                             // he = -1.0 * d * Trans(c);
                             static_cast<ElementByElementMatrix<SCAL>*>(harmonicext.get())
                               ->AddElementMatrix(i,idnums,ednums,he);
                                  
                             if (!symmetric)
                               {
                                 FlatMatrix<SCAL> het (sizeo, sizei, lh);
                                 het = 0.0;
                                 LapackMultAddAB (b, d, -1, het);
                                 //  het = -1.0 * b * d;
                                 static_cast<ElementByElementMatrix<SCAL>*>(harmonicexttrans.get())
                                   ->AddElementMatrix(i,ednums,idnums,het);
                               }
                             static_cast<ElementByElementMatrix<SCAL>*>(innersolve.get())
                               ->AddElementMatrix(i,idnums,idnums,d);
                             
                             LapackMultAddAB (b, he, 1.0, a);
                           }                                 
                         
                         
                         if (printelmat)
                           *testout << "schur matrix = " << a << endl;
                         
                         sum_elmat.Rows(odofs).Cols(odofs) = a;
                         
                         for (int k = 0; k < idofs1.Size(); k++)
                           dnums[idofs1[k]] = -1;
                       }
                   }
                 

                 AddElementMatrix (dnums, dnums, sum_elmat, el, lh);

                 for (auto pre : preconditioners)
                   pre -> AddElementMatrix (dnums, sum_elmat, el, lh);
               });
            
            progress.Done();
          }

      
        if (VB_parts[BND].Size())
          {
            RegionTimer reg(timerbound);
            ProgressOutput progress (ma, "assemble surface element", ma->GetNE(BND));

            IterateElements 
              (*fespace, BND, clh,  [&] (FESpace::Element el, LocalHeap & lh)
               {
                 progress.Update();
                 
                 const FiniteElement & fel = fespace->GetFE (el, lh);
                 ElementTransformation & eltrans = ma->GetTrafo (el, lh);
                 Array<int> dnums(fel.GetNDof(), lh);
                 fespace->GetDofNrs (el, dnums);

                 for (auto d : dnums)
                   if (d != -1) useddof[d] = true;

                 FlatVector<SCAL> elveclin (dnums.Size()*fespace->GetDimension(), lh);
                 FlatMatrix<SCAL> elmat (dnums.Size()*fespace->GetDimension(), lh);
                 FlatMatrix<SCAL> sum_elmat (dnums.Size()*fespace->GetDimension(), lh);
                 sum_elmat = SCAL(0.0);
                  
                 lin.GetIndirect (dnums, elveclin);
                 fespace->TransformVec (el, elveclin, TRANSFORM_SOL);

                 for (auto & bfi : VB_parts[BND])
                   {
                     if (!bfi->DefinedOn (el.GetIndex())) continue;
                     if (!bfi->DefinedOnElement (el.Nr())) continue;
                     
                     bfi->CalcLinearizedElementMatrix (fel, eltrans, elveclin, elmat, lh);
                     
                     fespace->TransformMat (el, elmat, TRANSFORM_MAT_LEFT_RIGHT);
                     
                     sum_elmat += elmat;
                        
                     if (printelmat) 
                       {
                         lock_guard<mutex> guard(addelmatboundary1_mutex);
                         testout->precision(8);
                         (*testout) << "surface-elnum= " << el.Nr() << endl;
                         (*testout) << "eltype " << fel.ElementType() << endl;
                         (*testout) << "boundary = " << ma->GetMaterial (el) << endl;
                         (*testout) << "integrator " << bfi->Name() << endl;
                         (*testout) << "dnums = " << endl << dnums << endl;
                         (*testout) << "elveclin = " << endl << elveclin << endl;
                         (*testout) << "elmat = " << endl << elmat << endl;
                       }
                   }
                 
                 AddElementMatrix (dnums, dnums, sum_elmat, el, lh);
                 
                 for (auto pre : preconditioners)
                   pre -> AddElementMatrix (dnums, sum_elmat, el, lh);
               });
            progress.Done();
          }
        

        if (facetwise_skeleton_parts[VOL].Size())
          throw Exception ("CalcLinearization for skeleton-VOL not implemented");
        if (elementwise_skeleton_parts.Size())
          throw Exception ("CalcLinearization for skeleton-VOL not implemented");
        

        if (facetwise_skeleton_parts[BND].Size())
          {
            int ne = ma->GetNE(BND);
            
            ProgressOutput progress (ma, "assemble skelton element", ne);
            
            ParallelForRange
              ( IntRange(ne), [&] ( IntRange r )
                {
                  LocalHeap lh = clh.Split();
                  Array<int> fnums, elnums, vnums, svnums, dnums;
                  
                  for (int i : r)
                    {
                      progress.Update();                      
                      HeapReset hr(lh);
                      ElementId sei(BND, i);
                      
                      if (!fespace->DefinedOn (BND,ma->GetElIndex (sei))) continue;
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
                          cout << "Surface fel:GetNDof() = " << fel.GetNDof() << endl;
                          cout << "dnums.Size() = " << dnums.Size() << endl;
                          
                          (*testout) << "fel:GetNDof() = " << fel.GetNDof() << endl;
                          (*testout) << "dnums.Size() = " << dnums.Size() << endl;
                          (*testout) << "dnums = " << dnums << endl;
                          throw Exception ( "Inconsistent number of degrees of freedom " );
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
                          
                          for (int k = 0; k < dnums.Size(); k++)
                            if (dnums[k] != -1)
                              useddof[dnums[k]] = true;
                          
                          int elmat_size = dnums.Size()*fespace->GetDimension();
                          FlatMatrix<SCAL> elmat(elmat_size, lh);
                          FlatVector<SCAL> elveclin(elmat_size, lh);

                          lin.GetIndirect (dnums, elveclin);
                          fespace->TransformVec (ei, elveclin, TRANSFORM_SOL);                     
                          
                          // original version did not compile on MacOS V
                          const FacetBilinearFormIntegrator & fbfi = 
                            dynamic_cast<const FacetBilinearFormIntegrator&>(*bfi);  
                          fbfi.CalcLinearizedFacetMatrix (fel,facnr,eltrans,vnums, seltrans, svnums, elveclin, elmat, lh);
                          
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
                            lock_guard<mutex> guard(addelmatboundary1_mutex);                            
                            AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), lh);
                          }
                        }//end for (numintegrators)
                    }//end for nse                  
                });//end of parallel
            progress.Done();
          } // if facetwise_skeleton_parts[BND].size
        
        



        
        if (fespace->specialelements.Size())
          cout << "special elements: " << fespace->specialelements.Size() << endl;

        for (int i = 0; i < fespace->specialelements.Size(); i++)
          {
            HeapReset hr(clh);
            const SpecialElement & el = *fespace->specialelements[i];
            el.GetDofNrs (dnums);
          
            for (int j = 0; j < dnums.Size(); j++)
              if (dnums[j] != -1)
                useddof[dnums[j]] = true;
          
            FlatMatrix<SCAL> elmat;
            el.Assemble (elmat, clh);
          
            AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), clh);
          }
      
      
        // add eps to avoid empty lines
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
                AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), clh);
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
                  AddElementMatrix (dnums, dnums, elmat, ElementId(BND,i), clh);
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
    bool hasbound = false;
    bool hasinner = false;
    bool hasskeletonbound = false;
    bool hasskeletoninner = false;
    int volumeintegrals = -1;
    for(int j=0;j<parts.Size();j++)
    {
      const BilinearFormIntegrator & bfi = *GetIntegrator(j);
      if (bfi.BoundaryForm())
        if (bfi.SkeletonForm())
          hasskeletonbound = true;
        else
          hasbound = true;
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
    bool needs_facet_loop = false;
    bool needs_element_boundary_loop = false;
    bool neighbor_testfunction = false;
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
            needs_facet_loop = true;
            facetvolumeintegrals = j;
          }
          if (!dgform.element_boundary && parts[j]->BoundaryForm())
          {
            needs_facet_loop = true;
            facetboundaryintegrals = j;
          }
          if (dgform.element_boundary)
          {
            throw Exception("Element boundary formulation is not implemented for tensor product spaces, please reformulate as skeleton integrals");
            needs_element_boundary_loop = true;
          }
        }
      }
      // do we need locks for neighbor - testfunctions ?
      for (int j = 0; j < NumIntegrators(); j++)
        if (parts[j] -> SkeletonForm())
        {
          auto dgform = parts[j] -> GetDGFormulation();
          if (dgform.neighbor_testfunction)
            neighbor_testfunction = true;
        }
    }
    else
      return;
      
    if(facetvolumeintegrals == -1 && facetboundaryintegrals == -1)
      return;
    auto & nels = tpfes->GetNels();
    auto & nfacets = tpfes->GetNFacets();
    timerfac1.Start();
    for (FlatArray<int> colfacets : spaces[0]->FacetColoring())
    {
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
          // Horzontal edge - get facet elements w.r.t. first direction
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
          // Horzontal edge - get facet elements w.r.t. second direction
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
  void S_BilinearForm<SCAL> :: AddMatrix1 (SCAL val,
                                           const BaseVector & x,
                                           BaseVector & y, LocalHeap & clh) const
  {
    static Timer timer ("Apply Matrix");
    static Timer timervol ("Apply Matrix - volume");
    static Timer timerbound ("Apply Matrix - boundary");
    static Timer timerDG ("Apply Matrix - DG");
    constexpr int tlevel = 4;
    static Timer timerDGpar ("Apply Matrix - DG par", tlevel);
    static Timer timerDGapply ("Apply Matrix - DG par apply", tlevel);
    static Timer timerDG1 ("Apply Matrix - DG 1", tlevel);
    static Timer timerDG2 ("Apply Matrix - DG 2", tlevel);
    static Timer timerDG2a ("Apply Matrix - DG 2a", tlevel);
    static Timer timerDG2b ("Apply Matrix - DG 2b", tlevel);
    static Timer timerDG2c ("Apply Matrix - DG 2c", tlevel);
    static Timer timerDG3 ("Apply Matrix - DG 3", tlevel);
    static Timer timerDG4 ("Apply Matrix - DG 4", tlevel);
    static Timer timerDGfacet ("Apply Matrix - DG boundary", tlevel);
    static Timer timerDGfacet1 ("Apply Matrix - DG boundary 1", tlevel);
    static Timer timerDGfacet2 ("Apply Matrix - DG boundary 2", tlevel);
    static Timer timerDGparallelfacets ("Apply Matrix - DG parallel facets");
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
        // bool hasbound = VB_parts[BND].Size();
        // bool hasinner = VB_parts[VOL].Size();
        // bool hasskeletonbound = VB_skeleton_parts[BND].Size();
        // bool hasskeletoninner = VB_skeleton_parts[VOL].Size();

        /*
        bool hasbound = false;
        bool hasinner = false;
        bool hasskeletonbound = false;
        bool hasskeletoninner = false;

        for (int j = 0; j < NumIntegrators(); j++)
          {
            const BilinearFormIntegrator & bfi = *GetIntegrator(j);
            if (bfi.BoundaryForm())
              if (bfi.SkeletonForm())
                hasskeletonbound = true;
              else
                hasbound = true;
            else
              if (bfi.SkeletonForm())
                hasskeletoninner = true; 
              else
                hasinner = true;
          }
        */      

        for (auto vb : { VOL, BND, BBND } )
          if (VB_parts[vb].Size())
            {
              RegionTimer reg (timervol);
              
              IterateElements 
                (*fespace, vb, clh, 
                 [&] (FESpace::Element el, LocalHeap & lh)
                 {
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
                       bfi->ApplyElementMatrix (fel, trafo, elvecx, elvecy, 0, lh);
                       
                       this->fespace->TransformVec (el, elvecy, TRANSFORM_RHS);
                       
                       elvecy *= val;
                       y.AddIndirect (dnums, elvecy, fespace->HasAtomicDofs());  // coloring	                      
                     }
                 });
            } // hasinner
        
        
        /*
        if (hasbound)
          {
            RegionTimer reg (timerbound);

		    // LocalHeap clh (lh_size*TaskManager::GetMaxThreads(), "biform-AddMatrix - Heap");
                    
		    IterateElements 
		      (*fespace, BND, clh, 
		       [&] (ElementId ei, LocalHeap & lh)
                       
		       {
                         HeapReset hr(lh);
                         
                         const FiniteElement & fel = fespace->GetFE (ei, lh);
                         ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
			 Array<int> dnums (fel.GetNDof(), lh);
                         fespace->GetDofNrs (ei, dnums);
                         
                         ApplyElementMatrix(x,y,val,dnums,eltrans,ei.Nr(),1,cnt,lh,&fel);
                       });
                  }
        */

        

        // if (hasskeletonbound||hasskeletoninner)
        {
          RegionTimer reg(timerDG);                    
          
                    /*
                    bool needs_facet_loop = false;
                    bool needs_element_boundary_loop = false;
                    
                    for (int j = 0; j < NumIntegrators(); j++)
                      if (parts[j] -> SkeletonForm())
                        {
                          auto dgform = parts[j] -> GetDGFormulation();
                          if (!dgform.element_boundary)   //  && !parts[j]->BoundaryForm())  // too much !!!!
                            needs_facet_loop = true;
                          if (dgform.element_boundary)
                            needs_element_boundary_loop = true;
                            // throw Exception ("No BilinearFormApplication-Implementation for Facet-Integrators yet");
                        }
                    */

                    
                    /*
                    // do we need locks for neighbor - testfunctions ?
                    bool neighbor_testfunction = false;
                    for (int j = 0; j < NumIntegrators(); j++)
                      if (parts[j] -> SkeletonForm())
                        {
                          auto dgform = parts[j] -> GetDGFormulation();
                          if (dgform.neighbor_testfunction)
                            neighbor_testfunction = true;
                        }
                    */

                    // if (needs_facet_loop && !fespace->UsesDGCoupling())
                    // throw Exception ("skeleton-form needs \"dgjumps\" : True flag for FESpace");

                    // facet-loop
          if ( (facetwise_skeleton_parts[VOL].Size() > 0) ||
               (facetwise_skeleton_parts[BND].Size() > 0) )
            
            for (auto colfacets : fespace->FacetColoring())
              {
              /*
              ParallelForRange
                (colfacets.Size(), [&] (IntRange r)
                 {
                   RegionTimer reg(timerDGpar);
                   LocalHeap lh = clh.Split();
                   Array<int> elnums(2, lh), elnums_per(2, lh), fnums1(6, lh), fnums2(6, lh), vnums1(8, lh), vnums2(8, lh);
                   
                   for (int i : r)
              */
                SharedLoop2 sl(colfacets.Size());

                // task_manager -> CreateJob
                ParallelJob
                  ( [&] (const TaskInfo & ti) 
                    {
                      LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);
                      // ArrayMem<int,100> temp_dnums;
                      RegionTimer reg(timerDGpar);
                      // LocalHeap lh = clh.Split();
                      Array<int> elnums(2, lh), elnums_per(2, lh), fnums1(6, lh), fnums2(6, lh), vnums1(8, lh), vnums2(8, lh);

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
#ifdef PARALLEL
			   if( (ma->GetDistantProcs (NodeId(StdNodeType(NT_FACET, ma->GetDimension()), facet)).Size() > 0) && (MyMPI_GetNTasks()>1) )
			     continue;
#endif
                           facet2 = ma->GetPeriodicFacet(facet);
                           if(facet2 > facet)
                             {
                               ma->GetFacetElements (facet2, elnums_per);
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
                                         
                               FlatVector<SCAL> elx(dnums.Size()*this->fespace->GetDimension(), lh),
                                 ely(dnums.Size()*this->fespace->GetDimension(), lh);
                               x.GetIndirect(dnums, elx);
                               
                               bfi->ApplyFacetMatrix (fel,facnr1,eltrans,vnums1, seltrans, vnums2, elx, ely, lh);
                               y.AddIndirect(dnums, ely, fespace->HasAtomicDofs());
                             } //end for (numintegrators)
                           
                           continue;
                         } // end if boundary facet
                       
                       // timerDG2.Start();
                       // timerDG2a.Start();
                       int el2 = elnums[1];
                       ElementId ei2(VOL, el2);
                       
                       fnums2 = ma->GetElFacets(ei2);
                       int facnr2 = fnums2.Pos(facet2);
                       
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

                       RegionTimer reg2(timerDGapply);                     
                       for (auto & bfi : facetwise_skeleton_parts[VOL])                                   
                         {
                           if (!bfi->DefinedOn (ma->GetElIndex (ei1))) continue; 
                           if (!bfi->DefinedOn (ma->GetElIndex (ei2))) continue; 
                           
                           bfi->ApplyFacetMatrix (fel1, facnr1, eltrans1, vnums1,
                                                  fel2, facnr2, eltrans2, vnums2, elx, ely, lh);

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
                       
                       ma->GetFacetElements(fnums1[facnr1],elnums);
                       if (elnums.Size()<2) {
#ifdef PARALLEL
			 if( (ma->GetDistantProcs (NodeId(StdNodeType(NT_FACET, ma->GetDimension()), fnums1[facnr1])).Size() > 0) && (MyMPI_GetNTasks()>1) )
			   continue;
#endif
                         if(ma->GetPeriodicFacet(fnums1[facnr1])!=fnums1[facnr1])
                           {
                             ma->GetFacetElements (ma->GetPeriodicFacet(fnums1[facnr1]), elnums_per);
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
                               cout << "Surface fel:GetNDof() = " << fel.GetNDof() << endl;
                               cout << "dnums.Size() = " << dnums.Size() << endl;
                               
                               (*testout) << "fel:GetNDof() = " << fel.GetNDof() << endl;
                               (*testout) << "dnums.Size() = " << dnums.Size() << endl;
                               (*testout) << "dnums = " << dnums << endl;
                               throw Exception ( "Inconsistent number of degrees of freedom " );
                             }
                           
                           
                           // for (int j = 0; j < NumIntegrators(); j++)
                           for (auto & bfi : elementwise_skeleton_parts)
                             {
                               // const BilinearFormIntegrator & bfi = *parts[j];
                               
                               // if (bfi.BoundaryForm()) continue;
                               // if (!bfi.SkeletonForm()) continue;
                               // if (!bfi.GetDGFormulation().element_boundary) continue;
                               /*
                                 int elmat_size = dnums.Size()*fespace->GetDimension();
                                 FlatMatrix<SCAL> elmat(elmat_size, lh);
                                 
                                 dynamic_cast<const FacetBilinearFormIntegrator&>(bfi).  
                                 CalcFacetMatrix (fel,facnr1,eltrans,vnums1, seltrans, elmat, lh);
                                 
                                 fespace->TransformMat (el1, false, elmat, TRANSFORM_MAT_LEFT_RIGHT);
                               */
                               
                               FlatVector<SCAL> elx(dnums.Size()*fespace->GetDimension(), lh),
                                 ely(dnums.Size()*fespace->GetDimension(), lh);
                               x.GetIndirect(dnums, elx);
                               
                               dynamic_cast<const FacetBilinearFormIntegrator&>(*bfi).
                                 ApplyFacetMatrix (fel,facnr1,eltrans,vnums1, seltrans, vnums2, elx, ely, lh);
                               
                               y.AddIndirect(dnums, ely);
                               
                             } //end for (numintegrators)
                           
                           continue;
                         } // end if boundary facet

                       
                       // RegionTimer reg2(timerDG2);
                       int el2 = elnums[0] + elnums[1] - el1;
                       // T_ElementId<VOL,2> ei2(el2);
                       ElementId ei2(VOL, el2);
                       
                       fnums2 = ma->GetElFacets(ei2);
                       int facnr2 = fnums2.Pos(ma->GetPeriodicFacet(fnums1[facnr1]));

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
                       
                       
                       //  for (int j = 0; j < NumIntegrators(); j++)
                       for (auto & bfi : elementwise_skeleton_parts)
                         {
                           // shared_ptr<BilinearFormIntegrator> bfi = parts[j];
                           // BilinearFormIntegrator * bfi = parts[j].get();
                           
                           // if (!bfi->SkeletonForm()) continue;
                           // if (bfi->BoundaryForm()) continue;
                           // if (!bfi->GetDGFormulation().element_boundary) continue;                                     
                           if (!bfi->DefinedOn (ma->GetElIndex (ei1))) continue; //TODO: treat as surface element
                           if (!bfi->DefinedOn (ma->GetElIndex (ei2))) continue; //TODO    
                           
                           FacetBilinearFormIntegrator * fbfi = 
                             dynamic_cast<FacetBilinearFormIntegrator*>(bfi.get());
                           
                           fbfi->ApplyFacetMatrix (fel1, facnr1, eltrans1, vnums1,
                                                   fel2, facnr2, eltrans2, vnums2, elx, ely, lh);
                           
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
                           
                           if (fbfi->GetDGFormulation().neighbor_testfunction)
                             {
                               int dim = this->fespace->GetDimension();
                               FlatVector<SCAL> swap_elx(elx.Size(), lh);
                               swap_elx.Range(0, dim*dnums2.Size()) = elx.Range(dim*dnums1.Size(), dim*dnums.Size());
                               swap_elx.Range(dim*dnums2.Size(), dim*dnums.Size()) = elx.Range(0, dim*dnums1.Size());
                               fbfi->ApplyFacetMatrix (fel2, facnr2, eltrans2, vnums2,
                                                       fel1, facnr1, eltrans1, vnums1, swap_elx, ely, lh);
                               y.AddIndirect(dnums1, ely.Range(dim*dnums1.Size(), dim*dnums.Size()));
                             }
                         }
                     }
                 }                             
               });
        }

	
#ifdef PARALLEL
	if( (MyMPI_GetNTasks()>1) &&
	    (mpi_facet_parts.Size()) )
	  {
	    RegionTimer rt(timerDGparallelfacets);
	    
	    //cout << "apply parallel DG facets, " << elementwise_skeleton_parts.Size() << " el-bound and " << facetwise_skeleton_parts[VOL].Size() << " facet parts" << ", " << mpi_facet_parts.Size() << " total parts " << endl;

	    /**
	       TODO:
		 - convert ranks to ngs_comm-ranks and send via that comm
		 - make integrators that are not defined everywhere working
		   (this requires checking DefinedOn on el-indices on both sides of the facet!
	     **/
	    
	    MPI_Comm mcomm = ma->GetCommunicator();
	    int mrank, mnp;
	    MPI_Comm_rank(mcomm, &mrank);
	    MPI_Comm_size(mcomm, &mnp);
	    Array<int> cnt(mnp);
	    Array<MPI_Request> reqs;
	    Array<MPI_Request> reqr;
	    LocalHeap &lh(clh);
	    Array<int> elnums(2, lh), fnums(6, lh), vnums(8, lh);

	    size_t ne = ma->GetNE(VOL);
	    BitArray fine_facet(ma->GetNFacets());
	    fine_facet.Clear();
	    Array<int> elfacets;
	    for (int i = 0; i < ne; ++i) {
	      auto elfacets = ma->GetElFacets(ElementId(VOL,i));
	      for (auto f : elfacets) fine_facet.Set(f);
	    }

	    auto mpi_loop_range = (have_mpi_facet_data)?Range(1,3):Range(0,3);
	    
	    for(auto loop:mpi_loop_range) {
	      cnt = 0;
	      for(auto facet:Range(ma->GetNFacets())) {
		NodeId facet_id(StdNodeType(NT_FACET, ma->GetDimension()), facet);
		if(!fine_facet.Test(facet)) continue;
		auto fdps = ma->GetDistantProcs(facet_id);
		//skip non-mpi facets
		if (fdps.Size() == 0)
		  continue;
		auto d = fdps[0];
		HeapReset hr(lh);
		
		ma->GetFacetElements(facet, elnums);

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
		    cnt[d] += trace_values.Size();
		  }
		  else if (loop == 1) {
		    FlatVector<SCAL> tmp(trace_values.Size(), &( send_table[d][cnt[d]] ));
		    tmp = trace_values;

		    cnt[d] += trace_values.Size();
		  }
		  else {
		    FlatVector<SCAL> trace_other(trace_values.Size(), &( recv_table[d][cnt[d]] ));
		    cnt[d]+= trace_values.Size();

		    FlatVector<SCAL> ely(dnums.Size()*this->fespace->GetDimension(), lh);
		    dynamic_cast<const FacetBilinearFormIntegrator*>(igt.get())->  
		      ApplyFromTraceValues(fel,facetnr,eltrans,vnums, trace_other,  elx, ely, lh);
			
		    y.AddIndirect(dnums, ely);
		  }
		}		
	      }
	      
	      if(loop==0) {
		send_table = Table<SCAL> (cnt);
		recv_table = Table<SCAL> (cnt);
		for(auto r:send_table)
		  r = -1;
		for(auto r:recv_table)
		  r = -2;
		have_mpi_facet_data = true;
	      }
	      else if(loop==1) {
		for(auto dp:Range(mnp))
		  if(send_table[dp].Size()) {
		    reqs.Append(MyMPI_ISend(send_table[dp], dp, MPI_TAG_SOLVE, mcomm));
		    reqr.Append(MyMPI_IRecv(recv_table[dp], dp, MPI_TAG_SOLVE, mcomm));
		  }
		MyMPI_WaitAll(reqr);
	      }
	    }
	  }	    
#endif
        


        
        if (fespace->specialelements.Size())
          {
            // LocalHeap lh(lh_size, "biform-AddMatrix (c)");
            Array<int> dnums;
            // ElementTransformation * dummy_eltrans = NULL;
            for (int i = 0; i < fespace->specialelements.Size(); i++)
              {
                HeapReset hr(clh);
                const SpecialElement & el = *fespace->specialelements[i];
                el.GetDofNrs (dnums);

                FlatVector<SCAL> elvecx (dnums.Size() * fespace->GetDimension(), clh);
                FlatVector<SCAL> elvecy (dnums.Size() * fespace->GetDimension(), clh);
                
                x.GetIndirect (dnums, elvecx);

                el.Apply (elvecx, elvecy, clh);
                elvecy *= val;
                y.AddIndirect (dnums, elvecy);

                // ApplyElementMatrix(x,y,val,dnums,*dummy_eltrans,i,2,cnt,lh,NULL,&el);
              }
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

                   //for (int j = 0; j < this->NumIntegrators(); j++)
                   for (auto & bfi : VB_parts[vb])
                     {
                       // BilinearFormIntegrator & bfi = *this->parts[j];
                       // if (bfi.SkeletonForm()) continue;
                       // if (bfi.BoundaryForm()) continue;
                       if (!bfi->DefinedOn (this->ma->GetElIndex (ei))) continue;

                       MixedFiniteElement fel(fel1, fel2);
                       bfi->ApplyElementMatrix (fel, eltrans, elvecx, elvecy, 0, lh);
                       
                       this->fespace->TransformVec (ei, elvecy, TRANSFORM_RHS);
        
                       elvecy *= val;
                       y.AddIndirect (dnums2, elvecy);  // coloring	      
                     }
                 });
            }
        // cout << "apply not implemented for mixed" << endl;
      }
  }






  template <class SCAL>
  void S_BilinearForm<SCAL> :: ApplyLinearizedMatrixAdd1 (SCAL val,
                                                          const BaseVector & lin,
                                                          const BaseVector & x,
                                                          BaseVector & y) const
  {
    if (!MixedSpaces())
      
      {
        Array<int> dnums;
      
        int ne = ma->GetNE();
        int dim = GetFESpace()->GetDimension(); 
        LocalHeap lh (2000000, "biform-ApplyLinearized");

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
              
                  bfi.ApplyLinearizedElementMatrix (fel, eltrans, elveclin, elvecx, elvecy, lh);
                  fespace->TransformVec (sei, elvecy, TRANSFORM_RHS);
                  elvecy *= val;
                  y.AddIndirect (dnums, elvecy);
                }
            }

        for (int i = 0; i < fespace->specialelements.Size(); i++)
          {
            HeapReset hr(lh);
            const SpecialElement & el = *fespace->specialelements[i];
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
  double S_BilinearForm<SCAL> :: Energy (const BaseVector & x) const
  {
    static Timer t("BilinearForm::Energy"); RegionTimer reg(t);
    
    atomic<double> energy(0.0);

    if (!MixedSpaces())
      {
        LocalHeap lh (2000000, "biform-energy", true);

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

        for (int i = 0; i < fespace->specialelements.Size(); i++)
          {
            HeapReset hr(lh);

            const SpecialElement & el = *fespace->specialelements[i];
            el.GetDofNrs (dnums);

            FlatVector<SCAL> elvecx (dnums.Size() * GetFESpace()->GetDimension(), lh);
            x.GetIndirect (dnums, elvecx);
          
            energy += el.Energy (elvecx, lh);
          }
      }
    return energy;
  }


  template <class SCAL>
  void S_BilinearForm<SCAL> :: 
  AddDiagElementMatrix (const Array<int> & dnums1,
                        const FlatVector<SCAL> & diag,
                        bool inner_element, int elnr,
                        LocalHeap & lh)
  {
    throw Exception ("Baseclass::AddDiagElementMatrix");
  }
 



  template class S_BilinearForm<double>;
  template class S_BilinearForm<Complex>;

  template <class TM, class TV>
  T_BilinearForm<TM,TV>::
  T_BilinearForm (shared_ptr<FESpace> afespace, const string & aname, const Flags & flags)
    : S_BilinearForm<TSCAL> (afespace, aname, flags)
  { 
    if (this->fespace->LowOrderFESpacePtr())
      this->low_order_bilinear_form = 
        make_shared<T_BilinearForm<TM,TV>> 
        (this->fespace->LowOrderFESpacePtr(), aname+string(" low-order"), flags);
  }

  template <class TM, class TV>
  T_BilinearForm<TM,TV>::
  T_BilinearForm (shared_ptr<FESpace> afespace, 
                  shared_ptr<FESpace> afespace2,
                  const string & aname,
                  const Flags & flags)
    : S_BilinearForm<TSCAL> (afespace, afespace2, aname, flags) 
  {
    ;
  }

  template <class TM, class TV>
  T_BilinearForm<TM,TV>::
  ~T_BilinearForm ()
  {
    /*
    for (int i = 0; i < this->mats.Size(); i++)
      {
      delete this->mats[i];
        this->mats[i].reset();
      }
    */
  }



  template <class TM, class TV>
  void T_BilinearForm<TM,TV>::
  AllocateMatrix ()
  {
    if (this->mats.Size() == this->ma->GetNLevels())
      return;

    MatrixGraph * graph = this->GetGraph (this->ma->GetNLevels()-1, false);

    auto spmat = make_shared<SparseMatrix<TM,TV,TV>> (*graph, 1);
    mymatrix = spmat.get();
    
    if (this->spd) spmat->SetSPD();
    shared_ptr<BaseMatrix> mat = spmat;


#ifdef PARALLEL
    if ( this->GetFESpace()->IsParallel() )
      mat = make_shared<ParallelMatrix> (mat, &this->GetFESpace()->GetParallelDofs());
#endif
    this->mats.Append (mat);

    delete graph;

    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size()-1; i++)
        this->mats[i].reset();
  }


  template <class TM, class TV>
  void T_BilinearForm<TM,TV>::
  CleanUpLevel ()
  {
    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size(); i++)
        {
          /*
          delete this->mats[i];
          this->mats[i] = 0;
          */
          this->mats[i].reset();
        }
  }



  template <class TM, class TV>
  shared_ptr<BaseVector> T_BilinearForm<TM, TV>::
  CreateVector() const
  {
    auto afespace = this->fespace;
#ifdef PARALLEL
    if ( afespace->IsParallel() )
      return make_shared<ParallelVVector<TV>> (afespace->GetNDof(), &afespace->GetParallelDofs());
    else
#endif
      return make_shared<VVector<TV>> (afespace->GetNDof());
  }





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
    
  
  template <class TM, class TV>
  void T_BilinearForm<TM,TV>::
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<TSCAL> elmat,
                    ElementId id,
                    LocalHeap & lh) 
  {
    /*
    BaseMatrix * hmat = this->mats.Last().get();
    
#ifdef PARALLEL
    ParallelMatrix * parmat = dynamic_cast<ParallelMatrix*> (hmat);
    if (parmat) hmat = &parmat->GetMatrix();
#endif   

    TMATRIX & mat = dynamic_cast<TMATRIX&> (*hmat);
    mat.AddElementMatrix (dnums1, dnums2, elmat, this->fespace->HasAtomicDofs());
    */
    mymatrix -> TMATRIX::AddElementMatrix (dnums1, dnums2, elmat, this->fespace->HasAtomicDofs());
  }


  template <class TM, class TV>
  void T_BilinearForm<TM,TV>::LapackEigenSystem(FlatMatrix<TSCAL> & elmat, LocalHeap & lh) const 
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


  /*
  template <class SCAL>
  void S_BilinearForm<SCAL> :: ApplyElementMatrix(const BaseVector & x,
                                                  BaseVector & y,
                                                  const SCAL & val,
                                                  const Array<int> & dnums,
                                                  const ElementTransformation & eltrans,
                                                  const int elnum,
                                                  const int type,
                                                  int & cnt,
                                                  LocalHeap & lh,
                                                  const FiniteElement * fel,
                                                  const SpecialElement * sel) const
  {
    FlatVector<SCAL> elvecx (dnums.Size() * this->fespace->GetDimension(), lh);
    FlatVector<SCAL> elvecy (dnums.Size() * this->fespace->GetDimension(), lh);

    x.GetIndirect (dnums, elvecx);

    if(type == 0 || type == 1)
      {
        this->fespace->TransformVec (ElementId(VorB(type), elnum), elvecx, TRANSFORM_SOL);

        for (int j = 0; j < this->NumIntegrators(); j++)
          {
            BilinearFormIntegrator & bfi = *this->parts[j];
            if (bfi.SkeletonForm()) continue;
            if (type == 0 && bfi.BoundaryForm()) continue;
            if (type == 0 && !bfi.DefinedOn (this->ma->GetElIndex (elnum))) continue;
            if (type == 1 && !bfi.BoundaryForm()) continue;
            if (type == 1 && !bfi.DefinedOn (this->ma->GetSElIndex (elnum))) continue;
            
            if (this->precompute)
              bfi.ApplyElementMatrix (*fel, eltrans, elvecx, elvecy, 
                                      this->precomputed_data[elnum*this->NumIntegrators()+j], lh);
            else
              bfi.ApplyElementMatrix (*fel, eltrans, elvecx, elvecy, 0, lh);

            this->fespace->TransformVec (ElementId(VorB(type), elnum), elvecy, TRANSFORM_RHS);

            elvecy *= val;
            y.AddIndirect (dnums, elvecy, fespace->HasAtomicDofs());  // coloring	    
          }
      }
    else if (type == 2)
      {
        sel->Apply (elvecx, elvecy, lh);
        elvecy *= val;
        y.AddIndirect (dnums, elvecy);
      }
                      
  }
  */

  template <class TM, class TV>
  T_BilinearFormSymmetric<TM,TV> :: 
  T_BilinearFormSymmetric (shared_ptr<FESpace> afespace, const string & aname,
                           const Flags & flags)
    : S_BilinearForm<TSCAL> (afespace, aname, flags) 
  {
    if (this->fespace->LowOrderFESpacePtr())
      {
        this->low_order_bilinear_form = 
          make_shared<T_BilinearFormSymmetric<TM,TV>> 
          (this->fespace->LowOrderFESpacePtr(), aname+string(" low-order"), flags);
      }
  }

  template <class TM, class TV>
  T_BilinearFormSymmetric<TM,TV> :: 
  ~T_BilinearFormSymmetric ()
  {
    /*
    for (int i = 0; i < this->mats.Size(); i++)
      delete this->mats[i];
    */
  }


  template <class TM, class TV>
  void T_BilinearFormSymmetric<TM,TV> :: 
  AllocateMatrix ()
  {
    if (this->mats.Size() == this->ma->GetNLevels())
      return;

    MatrixGraph * graph = this->GetGraph (this->ma->GetNLevels()-1, true);

    auto spmat = make_shared<SparseMatrixSymmetric<TM,TV>> (*graph, 1);
    mymatrix = spmat.get();
    
    if (this->spd) spmat->SetSPD();
    shared_ptr<BaseMatrix> mat = spmat;

#ifdef PARALLEL
    if ( this->GetFESpace()->IsParallel() )
      mat = make_shared<ParallelMatrix> (mat, &this->GetFESpace()->GetParallelDofs());
#endif
    this->mats.Append (mat);

    delete graph;

    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size()-1; i++)
        this->mats[i].reset();
  }


  template <class TM, class TV>
  void T_BilinearFormSymmetric<TM,TV>::
  CleanUpLevel ()
  {
    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size(); i++)
        {
          /*
          delete this->mats[i];
          this->mats[i] = 0;
          */
          this->mats[i].reset();
        }
  }


  template <class TM, class TV>
  shared_ptr<BaseVector> T_BilinearFormSymmetric<TM, TV>::
  CreateVector() const
  {
    auto afespace = this->fespace;
#ifdef PARALLEL
    if ( afespace->IsParallel() )
      return make_shared<ParallelVVector<TV>> (afespace->GetNDof(), &afespace->GetParallelDofs());
    else
#endif
      // return new VVector<TV> (afespace->GetNDof());
      return make_shared<VVector<TV>> (afespace->GetNDof());
  }



  template <class TM, class TV>
  void T_BilinearFormSymmetric<TM,TV> :: 
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<TSCAL> elmat,
                    ElementId id, 
                    LocalHeap & lh) 
  {
    /*
    BaseMatrix * hmat = this->mats.Last().get();

#ifdef PARALLEL
    ParallelMatrix * parmat = dynamic_cast<ParallelMatrix*> (hmat);
    if (parmat) hmat = &parmat->GetMatrix();
#endif   

    TMATRIX & mat = dynamic_cast<TMATRIX&> (*hmat);

    mat.AddElementMatrix (dnums1, elmat, this->fespace->HasAtomicDofs());
    */
    mymatrix -> TMATRIX::AddElementMatrix (dnums1, elmat, this->fespace->HasAtomicDofs());
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

#ifdef NONE
  template <class TM, class TV>
  void T_BilinearFormSymmetric<TM,TV>::ApplyElementMatrix(const BaseVector & x,
                                                          BaseVector & y,
                                                          const TSCAL & val,
                                                          const Array<int> & dnums,
                                                          const ElementTransformation & eltrans,
                                                          const int elnum,
                                                          const int type,
                                                          int & cnt,
                                                          LocalHeap & lh,
                                                          const FiniteElement * fel,
                                                          const SpecialElement * sel) const
  {
    FlatVector<typename mat_traits<TV>::TSCAL> elvecx (dnums.Size() * this->GetFESpace()->GetDimension(), lh);
    FlatVector<typename mat_traits<TV>::TSCAL> elvecy (dnums.Size() * this->GetFESpace()->GetDimension(), lh);
                      
    x.GetIndirect (dnums, elvecx);

    if(type == 0 || type == 1)
      {
        this->fespace->TransformVec (elnum, (type == 1), elvecx, TRANSFORM_SOL);

        for (int j = 0; j < this->NumIntegrators(); j++)
          {
            BilinearFormIntegrator & bfi = *this->parts[j];
            if (bfi.SkeletonForm()) continue;
            if (type == 0 && bfi.BoundaryForm()) continue;
            if (type == 0 && !bfi.DefinedOn (this->ma->GetElIndex (elnum))) continue;
            if (type == 1 && !bfi.BoundaryForm()) continue;
            
            
            static Timer elementtimer ("Element matrix application", 2);
            elementtimer.Start();
            
            if (this->precompute)
              // bfi.ApplyElementMatrix (*fel, eltrans, elvecx, elvecy, this->precomputed_data[cnt++], lh);
              bfi.ApplyElementMatrix (*fel, eltrans, elvecx, elvecy, 
                                      this->precomputed_data[elnum*this->NumIntegrators()+j], lh);
            else
              bfi.ApplyElementMatrix (*fel, eltrans, elvecx, elvecy, 0, lh);
            
            elementtimer.Stop();
            
            /*
              testout->precision (12);
              (*testout) << "el " << i << ", dom = " << ma->GetElIndex(i) << ",integrator = " << typeid(bfi).name() << endl
              << "elx = " << elvecx 
              << "ely = " << elvecy << endl;
            */
            BilinearForm::GetFESpace()->TransformVec (elnum, (type == 1), elvecy, TRANSFORM_RHS);
            elvecy *= val;
            y.AddIndirect (dnums, elvecy);
          }
      }
    else if (type == 2)
      {
        sel->Apply (elvecx, elvecy, lh);
        elvecy *= val;
        y.AddIndirect (dnums, elvecy);
      }
                      
  }
#endif






  template <class TM>
  T_BilinearFormDiagonal<TM> :: 
  T_BilinearFormDiagonal (shared_ptr<FESpace> afespace, const string & aname,
                          const Flags & flags)
    : S_BilinearForm<TSCAL> (afespace, aname, flags) 
  { 
    this->diagonal = 1;
    cout << " !!!!!!!!!!!!!!!!!! allocated diagonal matrix !!!!!!!!!!!!!" << endl;

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
  {
    /*
    for (int i = 0; i < this->mats.Size(); i++)
      delete this->mats[i];
    */
  }

  ///
  template <class TM>
  void T_BilinearFormDiagonal<TM> :: 
  AllocateMatrix ()
  {
    if (this->mats.Size() == this->ma->GetNLevels())
      return;

    int ndof = this->fespace->GetNDof();
    MatrixGraph * graph = new MatrixGraph (ndof, 1);
    for (int i = 0; i < ndof; i++)
      graph->CreatePosition (i, i);

    // graphs.Append (graph);
    this->mats.Append (make_shared<SparseMatrixSymmetric<TM>> (*graph, 1));
    delete graph;

    if (!this->multilevel || this->low_order_bilinear_form)
      for (int i = 0; i < this->mats.Size()-1; i++)
        {
          /*
          delete this->mats[i];
          this->mats[i] = 0;
          */
          this->mats[i].reset();
        }
  }


  template <class TM>
  shared_ptr<BaseVector> T_BilinearFormDiagonal<TM> :: 
  CreateVector() const
  {
    auto afespace = this->fespace;
#ifdef PARALLEL
    if ( afespace->IsParallel() )
      return make_shared<ParallelVVector<TV_COL>> (afespace->GetNDof(), &afespace->GetParallelDofs());
    else
#endif
      // return new VVector<TV_COL> (afespace->GetNDof());
      return make_shared<VVector<TV_COL>> (afespace->GetNDof());
  }

  ///
  template <class TM>
  void T_BilinearFormDiagonal<TM> :: 
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<TSCAL> elmat,
                    ElementId id, 
                    LocalHeap & lh) 
  {
    TMATRIX & mat = dynamic_cast<TMATRIX&> (*this->mats.Last());

    for (int i = 0; i < dnums1.Size(); i++)
      if (dnums1[i] != -1)
        {
          TM & mij = mat(dnums1[i], dnums1[i]);
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
                    ElementId id, 
                    LocalHeap & lh) 
  {
    TMATRIX & mat = dynamic_cast<TMATRIX&> (GetMatrix());

    for (int i = 0; i < dnums1.Size(); i++)
      if (dnums1[i] != -1)
        mat(dnums1[i], dnums1[i]) += elmat(i, i);
  }




  ///
  template <> void T_BilinearFormDiagonal<Complex>::
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<Complex> elmat,
                    ElementId id, 
                    LocalHeap & lh) 
  {
    TMATRIX & mat = dynamic_cast<TMATRIX&> (GetMatrix()); 

    for (int i = 0; i < dnums1.Size(); i++)
      if (dnums1[i] != -1)
        mat(dnums1[i], dnums1[i]) += elmat(i, i);
  }


#ifdef VERYOLD
  template <class TM >
  void T_BilinearFormDiagonal<TM>::ApplyElementMatrix(const BaseVector & x,
                                                      BaseVector & y,
                                                      const TSCAL & val,
                                                      const Array<int> & dnums,
                                                      const ElementTransformation & eltrans,
                                                      const int elnum,
                                                      const int type,
                                                      int & cnt,
                                                      LocalHeap & lh,
                                                      const FiniteElement * fel,
                                                      const SpecialElement * sel) const
  {
    typedef typename mat_traits<TM>::TV_ROW TV_ROW;
    typedef typename mat_traits<TM>::TV_COL TV_COL;

    FlatVector<typename mat_traits<TV_ROW>::TSCAL> elvecx (dnums.Size() * this->GetFESpace()->GetDimension(), lh);
    FlatVector<typename mat_traits<TV_COL>::TSCAL> elvecy (dnums.Size() * this->GetFESpace()->GetDimension(), lh);
                      
    x.GetIndirect (dnums, elvecx);

    if(type == 0 || type == 1)
      {
        ElementId ei(VorB(type), elnum);
        this->fespace->TransformVec (ei, elvecx, TRANSFORM_SOL);

        for (int j = 0; j < this->NumIntegrators(); j++)
          {
            BilinearFormIntegrator & bfi = *this->parts[j];
            if (bfi.SkeletonForm()) continue;
            if (type == 0 && bfi.BoundaryForm()) continue;
            if (type == 0 && !bfi.DefinedOn (this->ma->GetElIndex (elnum))) continue;
            if (type == 1 && !bfi.BoundaryForm()) continue;
            
            
            static Timer elementtimer ("Element matrix application");
            elementtimer.Start();
            
            if (this->precompute)
              bfi.ApplyElementMatrix (*fel, eltrans, elvecx, elvecy, this->precomputed_data[cnt++], lh);
            else
              bfi.ApplyElementMatrix (*fel, eltrans, elvecx, elvecy, 0, lh);

            elementtimer.Stop();
            
            /*
              testout->precision (12);
              (*testout) << "el " << i << ", dom = " << ma->GetElIndex(i) << ",integrator = " << typeid(bfi).name() << endl
              << "elx = " << elvecx 
              << "ely = " << elvecy << endl;
            */
            BilinearForm::GetFESpace()->TransformVec (ei, elvecy, TRANSFORM_RHS);
        
            elvecy *= val;
            y.AddIndirect (dnums, elvecy);
          }
      }
    else if (type == 2)
      {
        sel->Apply (elvecx, elvecy, lh);
        elvecy *= val;
        y.AddIndirect (dnums, elvecy);
      }
                      
  }
#endif


  template <class TM>
  void T_BilinearFormDiagonal<TM> :: 
  AddDiagElementMatrix (const Array<int> & dnums1,
                        const FlatVector<TSCAL> & diag,
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
  AddDiagElementMatrix (const Array<int> & dnums1,
                        const FlatVector<double> & diag,
                        bool inner_element, int elnr,
                        LocalHeap & lh) 
  {
    TMATRIX & mat = dynamic_cast<TMATRIX&> (GetMatrix());

    for (int i = 0; i < dnums1.Size(); i++)
      if (dnums1[i] != -1)
        mat(dnums1[i], dnums1[i]) += diag(i);
  }

  ///
  template <> void T_BilinearFormDiagonal<Complex>::
  AddDiagElementMatrix (const Array<int> & dnums1,
                        const FlatVector<Complex> & diag,
                        bool inner_element, int elnr,
                        LocalHeap & lh) 
  {
    TMATRIX & mat = dynamic_cast<TMATRIX&> (GetMatrix()); 

    for (int i = 0; i < dnums1.Size(); i++)
      if (dnums1[i] != -1)
        mat(dnums1[i], dnums1[i]) += diag(i);
  }














  ComponentBilinearForm :: ComponentBilinearForm (shared_ptr<BilinearForm> abase_blf, int acomp, int ancomp)
    : BilinearForm( (*dynamic_pointer_cast<CompoundFESpace> (abase_blf->GetFESpace()))[acomp], "comp-lf", Flags()), 
      base_blf(abase_blf), comp(acomp) // , ncomp(ancomp)
  { 
    ;
  }

  BilinearForm & ComponentBilinearForm :: AddIntegrator (shared_ptr<BilinearFormIntegrator> bfi)
  {
    cout << "adding a block-bfi integrator" << endl;
    auto block_bfi = make_shared<CompoundBilinearFormIntegrator> (bfi, comp);
    block_bfi->SetDefinedOn (bfi->GetDefinedOn());
    base_blf -> AddIntegrator (block_bfi);
    cout << "comp is defined on : " << block_bfi->GetDefinedOn() << endl;
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

    
    bool symmetric_storage = flags.GetDefineFlag ("symmetric") || flags.GetDefineFlag ("spd");
    if (flags.GetDefineFlag("nonsym") || flags.GetDefineFlag("nonsym_storage")) symmetric_storage = false;
    
    
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
    SparseMatrix<double>* prolMat = NULL;

    if ( !low_order_bilinear_form )
      for (int i = GetNLevels()-1; i>0; i--)
        {
          prolMat = prol->CreateProlongationMatrix( i );
          
          /*
            GetMatrix( i-1 ) = 
            *( dynamic_cast< const BaseSparseMatrix& >( GetMatrix( i ) ).
            Restrict( *prolMat, &( dynamic_cast< BaseSparseMatrix& >
            ( GetMatrix( i-1 ) ) ) ) );
          */
          mats[i-1] =
            shared_ptr<BaseMatrix>
            ( dynamic_cast< const BaseSparseMatrix& >(GetMatrix(i) ).
              Restrict( *prolMat, &( dynamic_cast< BaseSparseMatrix& >
                                     (GetMatrix( i-1 ) ) ) ) );
          
          delete prolMat;
        }
  }

  
  BilinearFormApplication :: 
  BilinearFormApplication (shared_ptr<BilinearForm> abf)
    : bf (abf)
  {
    ;
  }

  void BilinearFormApplication :: 
  Mult (const BaseVector & v, BaseVector & prod) const
  {
    v.Cumulate();

    prod = 0;

    bool done = false;
    static int lh_size = 10*1000*1000;
    
    while(!done && lh_size < 1000*1000*1000)
      {
        try
          {
            LocalHeap lh(lh_size*TaskManager::GetMaxThreads(), "biform-AddMatrix - Heap");
            bf -> AddMatrix (1, v, prod, lh);
            done = true;
          }            
        catch (LocalHeapOverflow lhex)
          {
            lh_size *= 5;
            cerr << "Trying automatic heapsize increase to " << lh_size << endl;
          }
      }    

    prod.SetParallelStatus (DISTRIBUTED);
  }

  void BilinearFormApplication :: 
  MultAdd (double val, const BaseVector & v, BaseVector & prod) const
  {
    v.Cumulate();
    prod.Distribute();

    // bf -> AddMatrix (val, v, prod);

    bool done = false;
    static int lh_size = 10*1000*1000;
    
    while(!done && lh_size < 1000*1000*1000)
      {
        try
          {
            LocalHeap lh(lh_size*TaskManager::GetMaxThreads(), "biform-AddMatrix - Heap");
            bf -> AddMatrix (val, v, prod, lh);
            done = true;
          }            
        catch (LocalHeapOverflow lhex)
          {
            lh_size *= 5;
            cerr << "Trying automatic heapsize increase to " << lh_size << endl;
          }
      }    
    
  }

  void BilinearFormApplication :: 
  MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const
  {
    v.Cumulate();
    prod.Distribute();

    // bf -> AddMatrix (val, v, prod);

    bool done = false;
    static int lh_size = 10*1000*1000;
    
    while(!done && lh_size < 1000*1000*1000)
      {
        try
          {
            LocalHeap lh(lh_size*TaskManager::GetMaxThreads(), "biform-AddMatrix - Heap");
            bf -> AddMatrix (val, v, prod, lh);
            done = true;
          }            
        catch (LocalHeapOverflow lhex)
          {
            lh_size *= 5;
            cerr << "Trying automatic heapsize increase to " << lh_size << endl;
          }
      }    
    
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
    return bf -> CreateVector();
  }

  LinearizedBilinearFormApplication ::
  LinearizedBilinearFormApplication (shared_ptr<BilinearForm> abf,
                                     const BaseVector * aveclin)
    : BilinearFormApplication (abf), veclin(aveclin)
  {
    ;
  }

  void  LinearizedBilinearFormApplication :: 
  Mult (const BaseVector & v, BaseVector & prod) const
  {
    prod = 0;
    bf->ApplyLinearizedMatrixAdd (1, *veclin, v, prod);
  }

  void LinearizedBilinearFormApplication :: 
  MultAdd (double val, const BaseVector & v, BaseVector & prod) const
  {
    bf->ApplyLinearizedMatrixAdd (val, *veclin, v, prod);
  }

  void LinearizedBilinearFormApplication :: 
  MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const
  {
    bf->ApplyLinearizedMatrixAdd (val, *veclin, v, prod);
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
  shared_ptr<BaseVector> ElementByElement_BilinearForm<SCAL> :: CreateVector() const
  {
    return make_shared<VVector<SCAL>> (this->fespace->GetNDof());
  }

  template<class SCAL>
  void ElementByElement_BilinearForm<SCAL> :: 
  AddElementMatrix (FlatArray<int> dnums1,
                    FlatArray<int> dnums2,
                    BareSliceMatrix<SCAL> elmat,
                    ElementId id,
                    LocalHeap & lh)
  {
    int nr = id.Nr();
    if (id.IsBoundary()) nr += this->ma->GetNE();
    dynamic_cast<ElementByElementMatrix<SCAL>&> (this->GetMatrix()).AddElementMatrix (nr, dnums1, dnums2, elmat);
  }
  

  template class ElementByElement_BilinearForm<double>;
  template class ElementByElement_BilinearForm<Complex>;

  template class T_BilinearForm<double,double>;
  template class T_BilinearFormSymmetric<double,double>;
}
