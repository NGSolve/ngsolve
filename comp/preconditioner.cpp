#include <comp.hpp>

#include <multigrid.hpp>
//using namespace ngmg;

#include <solve.hpp>

#include <parallelngs.hpp>

namespace ngcomp
{
  using namespace ngmg;
  using namespace ngcomp;
  using namespace ngparallel;

  /*
  Preconditioner :: Preconditioner ()
  {
    test = false;
    timing = false;
    print = false;
    laterupdate = false;
    testresult_ok = testresult_min = testresult_max = NULL;
  }
 

  Preconditioner :: Preconditioner (const Flags & aflags)
    :flags(aflags) 
  {
    test = flags.GetDefineFlag ("test");
    timing = flags.GetDefineFlag ("timing");
    print = flags.GetDefineFlag ("print");
    laterupdate = flags.GetDefineFlag ("laterupdate");
    testresult_ok = testresult_min = testresult_max = NULL;
  }
  */  

  Preconditioner :: Preconditioner (const PDE * const apde, const Flags & aflags, const string aname)
    : NGS_Object(apde->GetMeshAccess(), aname), flags(aflags)
  {
    test = flags.GetDefineFlag ("test");
    timing = flags.GetDefineFlag ("timing");
    print = flags.GetDefineFlag ("print");
    laterupdate = flags.GetDefineFlag ("laterupdate");
    testresult_ok = testresult_min = testresult_max = NULL;
    // using lapack for testing EV
    uselapack = flags.GetDefineFlag ("lapacktest");
    if ( uselapack ) test = 1;

    if(test)
      {
	string testresult_ok_name = flags.GetStringFlag ("testresultok","");
	string testresult_min_name = flags.GetStringFlag ("testresultmin","");
	string testresult_max_name = flags.GetStringFlag ("testresultmax","");
	
	if(testresult_ok_name != "") testresult_ok = &(const_cast<PDE*>(apde)->GetVariable(testresult_ok_name));
	if(testresult_min_name != "") testresult_min = &(const_cast<PDE*>(apde)->GetVariable(testresult_min_name));
	if(testresult_max_name != "") testresult_max = &(const_cast<PDE*>(apde)->GetVariable(testresult_max_name));
      }

    on_proc = int ( flags.GetNumFlag("only_on", -1));
  }
  
  Preconditioner :: ~Preconditioner ()
  {
    ;
  }

  void Preconditioner :: Test () const
  {
    cout << "Compute eigenvalues" << endl;
    const BaseMatrix & amat = GetAMatrix();
    const BaseMatrix & pre = GetMatrix();
    int eigenretval;

    if ( !this->uselapack )
      {
        EigenSystem eigen (amat, pre);
        eigen.SetPrecision(1e-30);
        eigen.SetMaxSteps(1000); 
        
        eigen.SetPrecision(1e-15);
        eigenretval = eigen.Calc();
        eigen.PrintEigenValues (*testout);
        (cout) << " Min Eigenvalue : "  << eigen.EigenValue(1) << endl; 
        (cout) << " Max Eigenvalue : " << eigen.MaxEigenValue() << endl; 
        (cout) << " Condition   " << eigen.MaxEigenValue()/eigen.EigenValue(1) << endl; 
        (*testout) << " Min Eigenvalue : "  << eigen.EigenValue(1) << endl; 
        (*testout) << " Max Eigenvalue : " << eigen.MaxEigenValue() << endl; 
        
        if(testresult_ok) *testresult_ok = eigenretval;
        if(testresult_min) *testresult_min = eigen.EigenValue(1);
        if(testresult_max) *testresult_max = eigen.MaxEigenValue();
        
        
        //    (*testout) << " Condition   " << eigen.MaxEigenValue()/eigen.EigenValue(1) << endl; 
        //    for (int i = 1; i < min2 (10, eigen.NumEigenValues()); i++)
        //      cout << "cond(i) = " << eigen.MaxEigenValue() / eigen.EigenValue(i) << endl;
        (*testout) << " Condition   " << eigen.MaxEigenValue()/eigen.EigenValue(1) << endl;
        
      }

    else
      {
#ifdef LAPACK   
        int n = amat.Height();
	int n_elim = 0;
	BitArray internaldofs (n);
	internaldofs.Clear();

	for ( int i = 0; i < n; i++ )
	  {
	    const FlatArray<const int> rowindices = 
	      dynamic_cast<const BaseSparseMatrix&>(amat).GetRowIndices(i);
	    if ( rowindices.Size() <= 1 )
	      internaldofs.Set(i);
	    else
	      n_elim++;
	  }

        Matrix<Complex> mat(n_elim), mat2(n_elim), ev(n_elim);
        BaseVector & v1 = *amat.CreateVector();
        BaseVector & v2 = *amat.CreateVector();
        FlatVector<Complex> fv1 = v1.FVComplex();
        FlatVector<Complex> fv2 = v2.FVComplex();
        
	int i_elim = 0, j_elim = 0;
        for (int i = 0; i < n; i++)
          {
	    if ( internaldofs.Test(i) ) continue;
	    j_elim = 0;
            fv1 = 0.0;
            fv1(i) = 1.0;
            v2 = amat * v1;
            v1 = pre * v2;

            for (int j = 0; j < n; j++)
	      {
		if ( internaldofs.Test(j) ) continue;
		mat(j_elim,i_elim) = fv1(j);
		j_elim++;
	      }
	    i_elim++;
          }

        mat2 = Complex(0.0);
        for (int i = 0; i < n_elim; i++)
          mat2(i,i) = 1.0;
        
        cout << "call lapack" << endl;
        Vector<Complex> lami(n_elim);
        // LapackEigenValues (mat, lami);
        LaEigNSSolve (n_elim, &mat(0,0), &mat2(0,0), &lami(0), 1, &ev(0,0), 0, 'B');

        ofstream out ("eigenvalues.out");
        for (int i = 0; i < n_elim; i++)
          out << lami(i).real() << " " << lami(i).imag() << "\n";
#endif
	;
      }
  }



  void Preconditioner :: Timing () const
  {
    cout << "Timing Preconditioner ... " << flush;
    const BaseMatrix & amat = GetAMatrix();
    const BaseMatrix & pre = GetMatrix();

    clock_t starttime;
    double time;
    starttime = clock();

    BaseVector & vecf = *pre.CreateVector();
    BaseVector & vecu = *pre.CreateVector();

    vecf = 1;
    int steps = 0;
    do
      {
	vecu = pre * vecf;
	steps++;
	time = double(clock() - starttime) / CLOCKS_PER_SEC;
      }
    while (time < 2.0);

    cout << " 1 step takes " 
	 << time / steps
	 << " seconds" << endl;



    starttime = clock();
    steps = 0;
    do
      {
	vecu = amat * vecf;
	steps++;
	time = double(clock() - starttime) / CLOCKS_PER_SEC;
      }
    while (time < 2.0);

    cout << ", 1 matrix takes " 
	 << time / steps
	 << " seconds" << endl;

  }





  MGPreconditioner :: MGPreconditioner (PDE * pde, Flags & aflags, const string aname)
    : Preconditioner (pde,aflags,aname)
  {
    // cout << "*** MGPreconditioner constructor called" << endl;
    mgtest = flags.GetDefineFlag ("mgtest");
    mgfile = flags.GetStringFlag ("mgfile","mgtest.out"); 
    mgnumber = int(flags.GetNumFlag("mgnumber",1)); 
   
 

    const MeshAccess & ma = pde->GetMeshAccess();
    bfa = pde->GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL));

    const LinearForm * lfconstraint = 
      pde->GetLinearForm (flags.GetStringFlag ("constraint", ""),1);
    const FESpace & fes = bfa->GetFESpace();

    

    const BilinearForm * lo_bfa = bfa;
    const FESpace * lo_fes = &fes;

    if (&bfa->GetLowOrderBilinearForm())
      {
	lo_bfa = &bfa->GetLowOrderBilinearForm();
	lo_fes = &fes.LowOrderFESpace();
      }
    else if (id == 0 && ntasks > 1 )
      {
	lo_bfa = bfa;
	lo_fes = &fes;
      }
    
    Smoother * sm = NULL;
    //    const char * smoother = flags.GetStringFlag ("smoother", "point");
    smoothertype = flags.GetStringFlag ("smoother", "point");

    if (smoothertype == "point")
      {
	sm = new GSSmoother (ma, *lo_bfa);
      }
    else if (smoothertype == "line")
      {
	sm = new AnisotropicSmoother (ma, *lo_bfa);
      }
    else if (smoothertype == "block") 
      {
	if (!lfconstraint)
	  sm = new BlockSmoother (ma, *lo_bfa, flags);
	else
	  sm = new BlockSmoother (ma, *lo_bfa, *lfconstraint, flags);
      }
    /*
    else if (smoothertype == "potential")
      {
	sm = new PotentialSmoother (ma, *lo_bfa);
      }
    */
    else
      cerr << "Unknown Smoother " << smoothertype << endl;

    if (!sm)
      throw Exception ("smoother could not be allocated"); 

    const Prolongation * prol = lo_fes->GetProlongation();

    mgp = new MultigridPreconditioner (ma, *lo_fes, *lo_bfa, sm, 
				       const_cast<Prolongation*> (prol));
    mgp->SetOwnProlongation (0);
    mgp->SetSmoothingSteps (int(flags.GetNumFlag ("smoothingsteps", 1)));
    mgp->SetCycle (int(flags.GetNumFlag ("cycle", 1)));
    mgp->SetIncreaseSmoothingSteps (int(flags.GetNumFlag ("increasesmoothingsteps", 1)));
    mgp->SetCoarseSmoothingSteps (int(flags.GetNumFlag ("coarsesmoothingsteps", 1)));
    mgp->SetUpdateAll( flags.GetDefineFlag( "updateall" ) );

    MultigridPreconditioner::COARSETYPE ct = MultigridPreconditioner::EXACT_COARSE;
    const char * coarse = flags.GetStringFlag ("coarsetype", "direct");
    if (strcmp (coarse, "smoothing") == 0)
      ct = MultigridPreconditioner::SMOOTHING_COARSE;
    else if (strcmp (coarse, "cg") == 0)
      ct = MultigridPreconditioner::CG_COARSE;
    mgp->SetCoarseType (ct);
    

    coarse_pre = 
      pde->GetPreconditioner (flags.GetStringFlag ("coarseprecond", ""), 1);

    if (coarse_pre)
      mgp->SetCoarseType (MultigridPreconditioner::USER_COARSE);

    finesmoothingsteps = int (flags.GetNumFlag ("finesmoothingsteps", 1));

    tlp = 0;
    inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
  }

 
  MGPreconditioner :: ~MGPreconditioner()
  {
    delete mgp;
    delete tlp;
  }

  void MGPreconditioner :: FreeSmootherMem ()
  {
    if(mgp) mgp->FreeMem();
    if(tlp) tlp->FreeMem();
  }
  

  void MGPreconditioner :: Update ()
  {
    const BilinearForm * lo_bfa = &bfa->GetLowOrderBilinearForm();

    INVERSETYPE invtype, loinvtype;
    invtype = dynamic_cast<const BaseSparseMatrix & > (bfa->GetMatrix()).SetInverseType (inversetype);
    if (lo_bfa)
      loinvtype = dynamic_cast<const BaseSparseMatrix & > (lo_bfa->GetMatrix()) .SetInverseType (inversetype);


    mgp->Update();

    if (coarse_pre)
      {
	mgp->SetCoarseGridPreconditioner (&coarse_pre->GetMatrix());
	mgp->SetOwnCoarseGridPreconditioner(false);
      }

    if (&bfa->GetLowOrderBilinearForm() || ntasks > 1)
      {
	delete tlp;

	Smoother * fine_smoother = NULL;

	fine_smoother = new BlockSmoother (bfa->GetMeshAccess(), *bfa, flags);

	tlp = new TwoLevelMatrix (&bfa->GetMatrix(),
				  mgp,
				  fine_smoother,
				  bfa->GetMeshAccess().GetNLevels()-1);
	
	tlp -> SetSmoothingSteps (finesmoothingsteps);
	tlp -> Update();
      }
    else
      tlp = 0;

    if (timing) Timing();
    if (test) Test();
    if (mgtest) MgTest();

    dynamic_cast<const BaseSparseMatrix & > (bfa->GetMatrix()).SetInverseType ( invtype );
    if (lo_bfa)
      dynamic_cast<const BaseSparseMatrix & > (lo_bfa->GetMatrix()) .SetInverseType ( loinvtype );
  }
  
  
  void MGPreconditioner :: CleanUpLevel () 
  { 
    if (&bfa->GetLowOrderBilinearForm())
      {
	delete tlp;
	tlp = 0;
      }    
  }


  const ngla::BaseMatrix & MGPreconditioner :: GetMatrix() const
  {
    if (tlp)
      return *tlp;
    else
      return *mgp; 
  }

  void MGPreconditioner::PrintReport (ostream & ost)
  {
    ost << "Multigrid preconditioner" << endl
	<< "bilinear-form = " << bfa->GetName() << endl
	<< "smoothertype = " << smoothertype << endl;
  }


  void MGPreconditioner::MemoryUsage (Array<MemoryUsageStruct*> & mu) const
  {
    int olds = mu.Size();

    if (&GetMatrix()) GetMatrix().MemoryUsage (mu);;

    for (int i = olds; i < mu.Size(); i++)
      mu[i]->AddName (string(" mgpre ")); // +GetName());
  }
    



  void MGPreconditioner :: MgTest () const
  {
    cout << "Compute eigenvalues" << endl;
    const BaseMatrix & amat = GetAMatrix();
    const BaseMatrix & pre = GetMatrix();
    
    int eigenretval;
    
    EigenSystem eigen (amat, pre);
    eigen.SetPrecision(1e-30);
    eigen.SetMaxSteps(1000); 
    eigenretval = eigen.Calc();
    eigen.PrintEigenValues (*testout);
    (cout) << " Min Eigenvalue : "  << eigen.EigenValue(mgnumber) << endl; 
    (cout) << " Max Eigenvalue : " << eigen.MaxEigenValue() << endl; 
    (cout) << " Condition   " << eigen.MaxEigenValue()/eigen.EigenValue(mgnumber) << endl; 
    (*testout) << " Min Eigenvalue : "  << eigen.EigenValue(mgnumber) << endl; 
    (*testout) << " Max Eigenvalue : " << eigen.MaxEigenValue() << endl;
    (*testout) << " Condition   " << eigen.MaxEigenValue()/eigen.EigenValue(mgnumber) << endl;
    static ofstream condout (mgfile.c_str());

    // double cond;

    condout << bfa->GetFESpace().GetNDof() << "\t" << bfa->GetFESpace().GetOrder() << "\t" << eigen.EigenValue(mgnumber) << "\t" << eigen.MaxEigenValue() << "\t" 
	    << eigen.MaxEigenValue()/eigen.EigenValue(mgnumber) <<  "\t" << endl;
    
    if(testresult_ok) *testresult_ok = eigenretval;
    if(testresult_min) *testresult_min = eigen.EigenValue(mgnumber);
    if(testresult_max) *testresult_max = eigen.MaxEigenValue();

  }


  DirectPreconditioner :: DirectPreconditioner (PDE * pde, Flags & aflags, const string aname)
    : Preconditioner(pde,aflags,aname)
  {
    bfa = pde->GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL));
    inverse = NULL;
    inversetype = flags.GetStringFlag("inverse", "sparsecholesky");
  }

  DirectPreconditioner :: ~DirectPreconditioner()
  {
    if ( inverse )
      delete inverse;
  }

  void DirectPreconditioner :: Update ()
  {
    if ( inverse )
      delete inverse;
    try
      {
// 	inverse = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix())
// 	  .InverseMatrix();
	const BaseSparseMatrix & amatrix = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix());

	amatrix.SetInverseType ( inversetype );
	if ( this->on_proc == -1  || this->on_proc == id )
	  {
	    const BitArray * freedofs = bfa->GetFESpace().GetFreeDofs();

	    if (freedofs)
	      inverse = amatrix.InverseMatrix(freedofs);  // change to const BitArray *
	    else
	      inverse = amatrix.InverseMatrix();

	    // delete freedofs;

	    if (print)
	      (*testout) << "inverse = " << endl << (*inverse) << endl;
	  }
      }
    catch (exception & e)
      {
	throw Exception ("DirectPreconditioner: needs a sparse matrix (or has memory problems)");
      }
    // (*testout) << "mat = " << bfa->GetMatrix() << endl;
    // (*testout) << "inv = " << (*inverse) << endl;

    //    if (test) Test();
  }

  void DirectPreconditioner :: CleanUpLevel ()
  {
    if ( inverse )
      delete inverse;
    inverse = 0;
  }

  const ngla::BaseMatrix & DirectPreconditioner :: GetMatrix() const
  {
    return *inverse;
  }







  DNDDPreconditioner :: DNDDPreconditioner (PDE * apde, Flags & aflags, const string aname)
    : Preconditioner(apde,aflags,aname)  {
    pde = apde;
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL));
    inverse = NULL;
  }

  DNDDPreconditioner :: ~DNDDPreconditioner()
  {
    delete inverse;
  }

  void DNDDPreconditioner :: Update ()
  {
    cout << "Dirichlet-Neumann DD Preconditioner" << endl;

    int i, j;

    const MeshAccess & ma = pde->GetMeshAccess();

    int np = bfa->GetFESpace().GetNDof();
    int ne = ma.GetNE();
    Array<int> domain(np), dnums;
    const FESpace & space = bfa->GetFESpace();

    for (i = 0; i < np; i++)
      domain[i] = 1;

    for (i = 0; i < ne; i++)
      {
	space.GetDofNrs (i, dnums);
	int elind = ma.GetElIndex(i);
	if (elind == 2 || elind == 4)
	  for (j = 0; j < dnums.Size(); j++)
	    domain[dnums[j]] = 2;
      }

    delete inverse;
    //  inverse = bfa->GetMatrix().InverseMatrix();
    // inverse = new SparseSystemLDL<SysMatrixC1d,SysVectorC1d> (dynamic_cast<const SparseSystemMatrix<SysMatrixC1d,SysVectorC1d>&> (bfa->GetMatrix()), NULL, &domain);
  }

  const ngla::BaseMatrix & DNDDPreconditioner :: GetMatrix() const
  {
    return *inverse;
  }







  LocalPreconditioner :: LocalPreconditioner (PDE * pde, Flags & aflags, const string aname)
    : Preconditioner (pde,aflags,aname),
      coarse_pre(0)
  {
    bfa = pde->GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL));
    block = flags.GetDefineFlag ("block");
    jacobi = NULL;
    locprectest = flags.GetDefineFlag ("mgtest");
    locprecfile = flags.GetStringFlag ("mgfile","locprectest.out"); 

    string smoother = flags.GetStringFlag("smoother","");
    if ( smoother == "block" )
      block = true;

    // coarse-grid preconditioner only used in parallel!!
    ct = "NO_COARSE";
    const char * coarse = flags.GetStringFlag ("coarsetype", "nocoarse");
    if (strcmp (coarse, "smoothing") == 0)
      ct = "SMOOTHING_COARSE";
    else if (strcmp (coarse, "direct") == 0)
      ct = "DIRECT_COARSE";
    

    coarse_pre = 
      pde->GetPreconditioner (flags.GetStringFlag ("coarseprecond", ""), 1);
    if (coarse_pre)
      ct = "USER_COARSE";

  }
   

  LocalPreconditioner :: ~LocalPreconditioner()
  {
    delete jacobi;
  }

  void LocalPreconditioner :: Update ()
  {
    cout << "Update Local Preconditioner" << flush;
    delete jacobi;
    
    // const BaseSparseMatrix& amatrix = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix());
    // 	if ( inversetype != "none" )
    // 	amatrix.SetInverseType ( inversetype );
	
    int blocktype = int (flags.GetNumFlag ( "blocktype", -1));

    bool parallel = (this->on_proc == -1);
    if ( !parallel && id != this->on_proc )
      {
	jacobi = 0; 
	return;
      }

//     BaseBlockJacobiPrecond::COARSE_TYPE bbct;
//     switch ( ct  )
//       {
//       case NO_COARSE:
// 	bbct = NO_COARSE; break;
//       case DIRECT_COARSE:
// 	bbct = DIRECT_COARSE; break;
//       }

    if ( block && blocktype == -1 ) blocktype = 0;
    if ( blocktype >= 0 )
      {
	// new: blocktypes, specified in fespace
	if (bfa->UsesEliminateInternal())
	  flags.SetFlag("eliminate_internal");
	Table<int> * blocks = bfa->GetFESpace().CreateSmoothingBlocks(flags);
	jacobi = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix())
	  .CreateBlockJacobiPrecond(*blocks, 0, coarse_pre, parallel);
	dynamic_cast<BaseBlockJacobiPrecond&> (*jacobi) . InitCoarseType(ct);
      }
    else if (block)
      {
	cout << "\nFlag block deprecated: use -blocktype=<typeno> instead" << endl;
      }
    else
      {
	jacobi = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix())
	  .CreateJacobiPrecond(bfa->GetFESpace().GetFreeDofs());
      }


    
    if (test) Test();
    if(locprectest) LocPrecTest(); 
    //    cout << " done" << endl;
  }

  const ngla::BaseMatrix & LocalPreconditioner :: GetMatrix() const
  {
    return *jacobi;
  }

 void LocalPreconditioner :: LocPrecTest () const
  {
    cout << "Compute eigenvalues" << endl;
    const BaseMatrix & amat = GetAMatrix();
    const BaseMatrix & pre = GetMatrix();
    
    int eigenretval;
    
    EigenSystem eigen (amat, pre);
    eigen.SetPrecision(1e-30);
    eigen.SetMaxSteps(1000); 
    eigenretval = eigen.Calc();
    eigen.PrintEigenValues (*testout);
    (cout) << " Min Eigenvalue : "  << eigen.EigenValue(1) << endl; 
    (cout) << " Max Eigenvalue : " << eigen.MaxEigenValue() << endl; 
    (cout) << " Condition   " << eigen.MaxEigenValue()/eigen.EigenValue(1) << endl; 
    (*testout) << " Min Eigenvalue : "  << eigen.EigenValue(1) << endl; 
    (*testout) << " Max Eigenvalue : " << eigen.MaxEigenValue() << endl;
    (*testout) << " Condition   " << eigen.MaxEigenValue()/eigen.EigenValue(1) << endl;
    static ofstream condout (locprecfile.c_str());

    // double cond;

    condout << bfa->GetFESpace().GetNDof() << "\t" << bfa->GetFESpace().GetOrder() << "\t" << eigen.EigenValue(1) << "\t" << eigen.MaxEigenValue() << "\t" 
	    << eigen.MaxEigenValue()/eigen.EigenValue(1) <<  "\t" << endl;
    
    if(testresult_ok) *testresult_ok = eigenretval;
    if(testresult_min) *testresult_min = eigen.EigenValue(1);
    if(testresult_max) *testresult_max = eigen.MaxEigenValue();

  }




 
 
#ifdef OLD
  VEFC_Preconditioner :: VEFC_Preconditioner (PDE * pde, Flags & aflags, const string aname)
    : Preconditioner (pde,aflags,aname)
  {
    bfa = pde->GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL));
    jacobi = NULL;
  }

  VEFC_Preconditioner :: ~VEFC_Preconditioner()
  {
    delete jacobi;
  }

  void VEFC_Preconditioner :: Update ()
  {
    cout << "Update VEFC_ Preconditioner" << flush;
    //  delete jacobi;

    
    jacobi = new VEFC_Matrix (bfa->GetMatrix(),
			      bfa->GetFESpace());
    /*
    Table<int> * blocks = bfa->GetFESpace().CreateSmoothingBlocks();

    jacobi = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix())
      .CreateBlockJacobiPrecond(*blocks);
    */

    if (test) Test();
  }
#endif

  const ngla::BaseMatrix & VEFC_Preconditioner :: GetMatrix() const
  {
    return *jacobi;
  }


  TwoLevelPreconditioner :: TwoLevelPreconditioner (PDE * apde, Flags & aflags, const string aname)
    : Preconditioner(apde,aflags,aname)
  {
    pde = apde;
    bfa = pde->GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL));
    cpre = pde->GetPreconditioner (flags.GetStringFlag ("coarsepreconditioner", NULL));
    smoothingsteps = int (flags.GetNumFlag ("smoothingsteps", 1));
    premat = NULL;
  }

  TwoLevelPreconditioner :: ~TwoLevelPreconditioner()
  {
    delete premat;
  }

  void TwoLevelPreconditioner :: Update ()
  {
    /*
      delete premat;

      Smoother * smoother = new GSSmoother (pde->GetMeshAccess(), *bfa);
      //  Smoother * smoother = new EBESmoother (pde->GetMeshAccess(), *bfa);

      premat = new TwoLevelMatrix (&bfa->GetMatrix(),
      &cpre->GetMatrix(),
      smoother,
      pde->GetMeshAccess().GetNLevels());

      //			       bfa->GetMatrix().CreateJacobiPrecond());
      premat -> SetSmoothingSteps (smoothingsteps);
      premat -> Update();
    */

    /*
      cout << "2-Level Preconditioner" << endl;
      EigenSystem eigen (bfa->GetMatrix(), *premat);
      eigen.Calc();
      eigen.PrintEigenValues (cout);
    */
  }



  ComplexPreconditioner :: ComplexPreconditioner (PDE * apde, Flags & aflags, const string aname)
    : Preconditioner(apde,aflags,aname)
  {
    dim = int (flags.GetNumFlag ("dim", 1));
    cm = 0;
    creal = apde->GetPreconditioner (flags.GetStringFlag ("realpreconditioner",""));
  }

  ComplexPreconditioner :: ~ComplexPreconditioner()
  { 
    ;
  }

  void ComplexPreconditioner :: Update ()
  {
    delete cm;
    switch (dim)
      {
      case 1:
	cm = new Real2ComplexMatrix<double,Complex> (&creal->GetMatrix());
	break;
      case 2:
	cm = new Real2ComplexMatrix<Vec<2,double>,Vec<2,Complex> > (&creal->GetMatrix());
	break;
      case 3:
	cm = new Real2ComplexMatrix<Vec<3,double>,Vec<3,Complex> > (&creal->GetMatrix());
	break;
      case 4:
	cm = new Real2ComplexMatrix<Vec<4,double>,Vec<4,Complex> > (&creal->GetMatrix());
	break;
      default:
	cout << "Error: dimension " << dim << " for complex preconditioner not supported!" << endl;
      }
  }



  ChebychevPreconditioner :: ChebychevPreconditioner (PDE * apde, Flags & aflags, const string aname)
    : Preconditioner(apde,aflags,aname)
  {
    steps = int (flags.GetNumFlag ("steps",10.));
    cm = 0;
    csimple = apde->GetPreconditioner (flags.GetStringFlag ("csimple",""));
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform",""));
    test = flags.GetDefineFlag ("test");
  }

  ChebychevPreconditioner :: ~ChebychevPreconditioner()
  { 
    ;
  }

  void ChebychevPreconditioner :: Update ()
  {
    delete cm;

    cout << "Compute eigenvalues csimple" << endl;
    const BaseMatrix & amat = bfa->GetMatrix();
    const BaseMatrix & pre = csimple->GetMatrix();

    EigenSystem eigen (amat, pre);
    eigen.SetPrecision(1e-30);
    eigen.SetMaxSteps(1000); 
    eigen.Calc();

    double lmin = eigen.EigenValue(1); 
    double lmax = eigen.MaxEigenValue(); 
    
    (*testout) << " Min Eigenvalue csimple: "  << eigen.EigenValue(1) << endl; 
    (*testout) << " Max Eigenvalue csimple : " << eigen.MaxEigenValue() << endl; 
    (cout) << " Min Eigenvalue csimple: "  << eigen.EigenValue(1) << endl; 
    (cout)<< " Max Eigenvalue csimple: " << eigen.MaxEigenValue() << endl; 
    (*testout) << " Condition csimple  " << eigen.MaxEigenValue()/eigen.EigenValue(1) << endl; 
    (cout) << " Condition csimple" << eigen.MaxEigenValue()/eigen.EigenValue(1) << endl; 
    eigen.PrintEigenValues(cout); 
  
    cm = new ChebyshevIteration(amat, pre, steps); 
    cm->SetBounds(1-lmax,1-lmin); 
    if(test) Test(); 

  }





  ////////////////////////////////////////////////////////////////////////////////
  // added 08/19/2003, FB

  NonsymmetricPreconditioner :: NonsymmetricPreconditioner (PDE * apde, Flags & aflags, const string aname)
    : Preconditioner(apde,aflags,aname)
  {
    dim = int (flags.GetNumFlag ("dim", 1));
    cm = 0;
    cbase = apde->GetPreconditioner (flags.GetStringFlag ("basepreconditioner",""));
  }

  NonsymmetricPreconditioner :: ~NonsymmetricPreconditioner()
  { 
    ;
  }

  void NonsymmetricPreconditioner :: Update ()
  {
    delete cm;
    switch (dim)
      {
      case 2:
	cm = new Small2BigNonSymMatrix<double, Vec<2,double> > (&cbase->GetMatrix());
	break;
      case 4:
	cm = new Small2BigNonSymMatrix<Vec<2,double>, Vec<4,double> > (&cbase->GetMatrix());
	break;
      case 6:
	cm = new Small2BigNonSymMatrix<Vec<3,double>, Vec<6,double> > (&cbase->GetMatrix());
	break;
       case 8:
 	cm = new Small2BigNonSymMatrix<Vec<4,double>, Vec<8,double> > (&cbase->GetMatrix());
 	break;
//       case 2:
// 	cm = new Sym2NonSymMatrix<Vec<2,double> > (&cbase->GetMatrix());
// 	break;
//       case 4:
// 	cm = new Sym2NonSymMatrix<Vec<4,double> > (&cbase->GetMatrix());
// 	break;
//       case 6:
// 	cm = new Sym2NonSymMatrix<Vec<6,double> > (&cbase->GetMatrix());
// 	break;
//        case 8:
//  	cm = new Sym2NonSymMatrix<Vec<8,double> > (&cbase->GetMatrix());
//  	break;
      default:
	cout << "Error: dimension " << dim;
	cout << " for nonsymmetric preconditioner not supported!" << endl;
      }
  }





  ////////////////////////////////////////////////////////////////////////////////


  CommutingAMGPreconditioner :: CommutingAMGPreconditioner (PDE * apde, Flags & aflags, const string aname)
    : Preconditioner (apde,aflags,aname), pde(apde)
  {
    bfa = pde->GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    while (&bfa->GetLowOrderBilinearForm())
      bfa = &bfa->GetLowOrderBilinearForm();

    coefse = pde->GetCoefficientFunction (flags.GetStringFlag ("coefse", ""),1);    
    coefe = pde->GetCoefficientFunction (flags.GetStringFlag ("coefe", ""),1);    
    coeff = pde->GetCoefficientFunction (flags.GetStringFlag ("coeff", ""),1);    

    hcurl = dynamic_cast<const NedelecFESpace*> (&bfa->GetFESpace()) != 0;
    levels = int (flags.GetNumFlag ("levels", 10));
    coarsegrid = flags.GetDefineFlag ("coarsegrid");

    amg = NULL;
  }

  CommutingAMGPreconditioner :: ~CommutingAMGPreconditioner ()
  {
    delete amg;
  }

  void CommutingAMGPreconditioner :: Update ()
  {
    cout << "Update amg" << endl;

#ifdef PARALLEL
    if ( this->on_proc != id && this->on_proc != -1)
      {
	amg = NULL;
	return;
      }
#endif

    const MeshAccess & ma = pde->GetMeshAccess();

    int nedge = ma.GetNEdges(); 
    int nface = ma.GetNFaces(); 
    int nel = ma.GetNE();

    if (coarsegrid && ma.GetNLevels() > 1)
      return;
      

    delete amg;

    cout << "get edges" << endl;

    Array<INT<2> > e2v (nedge);
    for (int i = 0; i < nedge; i++)
      {
	// if (ma.GetClusterRepEdge (i) >= 0)
	if (1)
	  ma.GetEdgePNums (i, e2v[i][0], e2v[i][1]);
	else
	  e2v[i][0] = e2v[i][1] = -1;
      }

    cout << "get faces" << endl;

    Array<INT<4> > f2v (nface);
    Array<int> fpnums;
    for (int i = 0; i < nface; i++)
      {
	if (ma.GetClusterRepFace (i) >= 0)
	  {
	    ma.GetFacePNums (i, fpnums);
	    
	    f2v[i][3] = -1;
	    for (int j = 0; j < fpnums.Size(); j++)
	      f2v[i][j] = fpnums[j];
	  }
	else
	  {
	    for (int j = 0; j < 4; j++)
	      f2v[i][j] = -1;
	  }
      }
    //(*testout) << "f2v " << f2v << endl;

    cout << "compute hi" << endl;

    // compute edge and face weights:

    Array<double> weighte(nedge), weightf(nface);
    Array<double> hi;    // edge length
    Array<double> ai;    // area of face


    hi.SetSize(nedge);
    for (int i = 0; i < nedge; i++)
      {
	Vec<3> p1, p2, v;
	if (e2v[i][0] != -1 && e2v[i][1] != -1)
	  {
	    ma.GetPoint (e2v[i][0], p1);
	    ma.GetPoint (e2v[i][1], p2);
	    v = p2 - p1;
	    hi[i] = L2Norm (v);
	  }
	else
	  hi[i] = 1;
      }


    if (hcurl)
      {
	//cout << "comptue ai" << endl;
	ai.SetSize(nface);
	for (int i = 0; i < nface; i++)
	  {
	    Vec<3> p1, p2, p3, v1, v2, vn;
	    ma.GetPoint (f2v[i][0], p1);
	    ma.GetPoint (f2v[i][1], p2);
	    ma.GetPoint (f2v[i][2], p3);
	    v1 = p2 - p1;
	    v2 = p3 - p1;
	    vn = Cross (v1, v2);
	    ai[i] = L2Norm (vn);
	  }
      }

    Array<int> ednums(12), edorient(12);
    Array<int> fanums(12), faorient(12);
    LocalHeap lh (10000, "CommutingAMG");
    ElementTransformation eltrans;
    IntegrationPoint ip(0, 0, 0, 0);


    cout << "compute weights" << endl;

    weighte = 0;
    weightf = 0;
    for (int i = 0; i < nel; i++)
      {
	lh.CleanUp();
	ma.GetElEdges (i, ednums);

	ma.GetElementTransformation (i, eltrans, lh);
	SpecificIntegrationPoint<3,3> sip(ip, eltrans, lh);

	double vol = ma.ElementVolume (i);


	double vale = Evaluate (*coefe, sip);

	for (int j = 0; j < ednums.Size(); j++)
	  weighte[ednums[j]] += vale * vol / sqr (hi[ednums[j]]);

	if (hcurl)
	  {
	    ma.GetElFaces (i, fanums);
	    double valf = Evaluate (*coeff, sip);
	    for (int j = 0; j < fanums.Size(); j++)
	      weightf[fanums[j]] += valf * vol / sqr (ai[fanums[j]]);
	  }
      }

    int nsel = ma.GetNSE();
    if (coefse)
      for (int i = 0; i < nsel; i++)
	{
	  lh.CleanUp();
	  ma.GetSElEdges (i, ednums);

	  ma.GetSurfaceElementTransformation (i, eltrans, lh);
	  SpecificIntegrationPoint<2,3> sip(ip, eltrans, lh);

	  double vol = ma.SurfaceElementVolume (i);
	  double vale = Evaluate (*coefse, sip);

	  for (int j = 0; j < ednums.Size(); j++)
	    weighte[ednums[j]] += vale * vol / sqr (hi[ednums[j]]);
	}



    clock_t starttime, endtime;
    starttime = clock();

    Array< Vec<3> > vertices;
    /*
    int nv = ma.GetNV();
    vertices.SetSize(nv);
    for (int i = 0; i < nv; i++)
      ma.GetPoint(i, vertices[i]);
    */
    CommutingAMG * amgmat;
    if (hcurl)
      amgmat = new AMG_HCurl (bfa->GetMatrix(), vertices, e2v, f2v, weighte, weightf, levels);
    else
      amgmat = new AMG_H1 (bfa->GetMatrix(), e2v, weighte, levels);

    endtime = clock();
    cout << "AMG coarsening time = " << double(endtime - starttime)/CLOCKS_PER_SEC << " sec" << endl;

    amgmat->ComputeMatrices (dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix()));

    endtime = clock();
    cout << "AMG projection time = " << double(endtime - starttime)/CLOCKS_PER_SEC << " sec" << endl;

    cout << "Total NZE = " << amgmat->NZE() << endl;

    amg = amgmat;

    /*
    ChebyshevIteration * cheby =
      new ChebyshevIteration (bfa->GetMatrix(), *amgmat, 10);
    cheby -> SetBounds (0, 0.98);
    amg = cheby;
    */

    cout << "matrices done" << endl;

    if (timing) Timing();
    if (test) Test();
  }


  void CommutingAMGPreconditioner :: CleanUpLevel () 
  { 
    if (coarsegrid)
      return;

    delete amg;
    amg = 0;
  }







  PreconditionerClasses::PreconditionerInfo::
  PreconditionerInfo (const string & aname,
	       Preconditioner* (*acreator)(const PDE & pde, const Flags & flags))
    : name(aname), creator(acreator) 
  {
    ;
  }
  
  PreconditionerClasses :: PreconditionerClasses ()
  {
    ;
  }

  PreconditionerClasses :: ~PreconditionerClasses()
  {
    for(int i=0; i<prea.Size(); i++)
      delete prea[i];
  }
  
  void PreconditionerClasses :: 
  AddPreconditioner (const string & aname,
		     Preconditioner* (*acreator)(const PDE & pde, 
						 const Flags & flags))
  {
    prea.Append (new PreconditionerInfo(aname, acreator));
  }

  const PreconditionerClasses::PreconditionerInfo * 
  PreconditionerClasses::GetPreconditioner(const string & name)
  {
    for (int i = 0; i < prea.Size(); i++)
      {
	if (name == prea[i]->name)
	  return prea[i];
      }
    return 0;
  }

  void PreconditionerClasses :: Print (ostream & ost) const
  {
    ost << endl << "Preconditioners:" << endl;
    ost <<         "---------" << endl;
    ost << setw(20) << "Name" << endl;
    for (int i = 0; i < prea.Size(); i++)
      ost << setw(20) << prea[i]->name << endl;
  }

 
  PreconditionerClasses & GetPreconditionerClasses ()
  {
    static PreconditionerClasses fecl;
    return fecl;
  }




}




