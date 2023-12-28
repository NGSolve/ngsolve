// #include <comp.hpp>
#include "bilinearform.hpp"
#include "preconditioner.hpp"
#include <multigrid.hpp>
// #include <solve.hpp>
#include <parallelngs.hpp>


namespace ngcomp
{
  using namespace ngmg;

  
  Preconditioner :: Preconditioner (shared_ptr<BilinearForm> bfa, const Flags & aflags,
                                    const string aname)
    : NGS_Object(bfa? bfa->GetMeshAccess(): nullptr, aflags, aname), bf(bfa)
  {
    test = flags.GetDefineFlag ("test");
    timing = flags.GetDefineFlag ("timing");
    print = flags.GetDefineFlag ("print");
    laterupdate = flags.GetDefineFlag ("laterupdate");
    testresult_ok = testresult_min = testresult_max = NULL;
    // using lapack for testing EV
    uselapack = flags.GetDefineFlag ("lapacktest");
    if ( uselapack ) test = 1;

    on_proc = int ( flags.GetNumFlag("only_on", -1));
    if (!flags.GetDefineFlag ("not_register_for_auto_update"))
      {
        bfa->SetPreconditioner(this);
        is_registered = true;
      }

  }
  

  
  Preconditioner :: ~Preconditioner ()
  {
    if (auto bfp = bf.lock(); is_registered && bfp)
      bfp->UnsetPreconditioner(this);
  }

  void Preconditioner :: Test () const
  {
    cout << IM(1) << "Compute eigenvalues" << endl;
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
        cout << IM(1) << " Min Eigenvalue : "  << eigen.EigenValue(1) << endl; 
        cout << IM(1) << " Max Eigenvalue : " << eigen.MaxEigenValue() << endl; 
        cout << IM(1) << " Condition   " << eigen.MaxEigenValue()/eigen.EigenValue(1) << endl; 
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
	    FlatArray rowindices = 
	      dynamic_cast<const BaseSparseMatrix&>(amat).GetRowIndices(i);
	    if ( rowindices.Size() <= 1 )
	      internaldofs.SetBit(i);
	    else
	      n_elim++;
	  }

        Matrix<Complex> mat(n_elim), mat2(n_elim), ev(n_elim);
        BaseVector & v1 = *amat.CreateColVector();
        BaseVector & v2 = *amat.CreateColVector();
        FlatVector<Complex> fv1 = v1.FVComplex();
        // FlatVector<Complex> fv2 = v2.FVComplex();
        
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
    cout << IM(1) << "Timing Preconditioner ... " << flush;
    const BaseMatrix & amat = GetAMatrix();
    const BaseMatrix & pre = GetMatrix();

    clock_t starttime;
    double time;
    starttime = clock();

    BaseVector & vecf = *pre.CreateColVector();
    BaseVector & vecu = *pre.CreateColVector();

    vecf = 1;
    int steps = 0;
    do
      {
	vecu = pre * vecf;
	steps++;
	time = double(clock() - starttime) / CLOCKS_PER_SEC;
      }
    while (time < 2.0);

    cout << IM(1) << " 1 step takes " << time / steps << " seconds" << endl;


    starttime = clock();
    steps = 0;
    do
      {
	vecu = amat * vecf;
	steps++;
	time = double(clock() - starttime) / CLOCKS_PER_SEC;
      }
    while (time < 2.0);

    cout << IM(1) << ", 1 matrix takes "
	 << time / steps << " seconds" << endl;
  }


  void Preconditioner :: ThrowPreconditionerNotReady() const
  {
    throw Exception("Preconditioner not ready: either update manually, or use update from assembling");
  }







  MGPreconditioner :: MGPreconditioner (shared_ptr<BilinearForm> abfa, const Flags & aflags, const string aname)
    : Preconditioner (abfa,aflags,aname)
  {
    // cout << "*** MGPreconditioner constructor called" << endl;
    mgtest = flags.GetDefineFlag ("mgtest");
    mgfile = flags.GetStringFlag ("mgfile","mgtest.out"); 
    mgnumber = int(flags.GetNumFlag("mgnumber",1)); 
    

    shared_ptr<MeshAccess> ma = abfa->GetMeshAccess();
    bfa = abfa;
    // bfa -> SetPreconditioner (this);


    // shared_ptr<LinearForm> lfconstraint = pde.GetLinearForm (flags.GetStringFlag ("constraint", ""),1);
    shared_ptr<FESpace> fes = bfa->GetFESpace();
    

    shared_ptr<BilinearForm> lo_bfa = bfa;
    shared_ptr<FESpace> lo_fes = fes;

    if (bfa->GetLowOrderBilinearForm())
      {
	lo_bfa = bfa->GetLowOrderBilinearForm();
	lo_fes = fes->LowOrderFESpacePtr();
      }
    /*
    else if (id == 0 && ntasks > 1 )  // not supported anymore
      {
	lo_bfa = bfa;
	lo_fes = &fes;
      }
    */

    shared_ptr<Smoother> sm = nullptr;
    //    const char * smoother = flags.GetStringFlag ("smoother", "point");
    smoothertype = flags.GetStringFlag ("smoother", "point");

    if (smoothertype == "point")
      {
	sm = make_shared<GSSmoother> (*ma, *lo_bfa);
      }
    else if (smoothertype == "line")
      {
	sm = make_shared<AnisotropicSmoother> (*ma, *lo_bfa);
      }
    else if (smoothertype == "block") 
      {
	// if (!lfconstraint)
        sm = make_shared<BlockSmoother> (*ma, *lo_bfa, flags);
          /*
            else
            sm = new BlockSmoother (*ma, *lo_bfa, *lfconstraint, flags);
          */
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

    auto prol = lo_fes->GetProlongation();

    mgp = make_shared<MultigridPreconditioner> (lo_bfa, sm, prol);
    mgp->SetSmoothingSteps (int(flags.GetNumFlag ("smoothingsteps", 1)));
    mgp->SetCycle (int(flags.GetNumFlag ("cycle", 1)));
    mgp->SetIncreaseSmoothingSteps (int(flags.GetNumFlag ("increasesmoothingsteps", 1)));
    mgp->SetCoarseSmoothingSteps (int(flags.GetNumFlag ("coarsesmoothingsteps", 1)));
    mgp->SetUpdateAll( flags.GetDefineFlag( "updateall" ) );
    mgp->SetHarmonicExtensionProlongation (flags.GetDefineFlag("he_prolongation"));    
    mgp->SetUpdateAlways(flags.GetDefineFlag("updatealways"));

    MultigridPreconditioner::COARSETYPE ct = MultigridPreconditioner::EXACT_COARSE;
    const string & coarse = flags.GetStringFlag ("coarsetype", "direct");
    if (coarse == "smoothing")
      ct = MultigridPreconditioner::SMOOTHING_COARSE;
    else if (coarse == "cg")
      ct = MultigridPreconditioner::CG_COARSE;
    mgp->SetCoarseType (ct);
    

    // coarse_pre = pde.GetPreconditioner (flags.GetStringFlag ("coarseprecond", ""), 1);
    // if (coarse_pre) mgp->SetCoarseType (MultigridPreconditioner::USER_COARSE);

    finesmoothingsteps = int (flags.GetNumFlag ("finesmoothingsteps", 1));

    tlp = 0;
    inversetype = flags.GetStringFlag("inverse", GetInverseName (default_inversetype));
    GetMemoryTracer().Track(*mgp, "MultiGridPreconditioner");
  }


  void MGPreconditioner :: Update ()
  {
    static Timer t("MGPreconditioner::Update"); RegionTimer reg(t);
    
    shared_ptr<BilinearForm> lo_bfa = bfa->GetLowOrderBilinearForm();

    INVERSETYPE invtype, loinvtype = default_inversetype;
    invtype = dynamic_cast<const BaseSparseMatrix & > (bfa->GetMatrix()).SetInverseType (inversetype);
    if (lo_bfa)
      loinvtype = dynamic_cast<const BaseSparseMatrix & > (lo_bfa->GetMatrix()) .SetInverseType (inversetype);


    mgp->Update();

    if (coarse_pre)
      {
	mgp->SetCoarseGridPreconditioner (shared_ptr<BaseMatrix> (const_cast<BaseMatrix*>(&coarse_pre->GetMatrix()), NOOP_Deleter));
      }

    if (bfa->GetLowOrderBilinearForm()) //  || ntasks > 1) not supported anymore
      {
        static Timer t("MGPreconditioner::Update - fine precond"); RegionTimer reg(t);
        auto fine_smoother = make_shared<BlockSmoother> (*bfa->GetMeshAccess(), *bfa, flags);
        GetMemoryTracer().Track(*fine_smoother, "FineSmoother");
        tlp = make_shared<TwoLevelMatrix> (&bfa->GetMatrix(),
                                           &*mgp,
                                           fine_smoother,
                                           bfa->GetMeshAccess()->GetNLevels()-1);
        tlp -> SetSmoothingSteps (finesmoothingsteps);
        if (bfa->GetFESpace()->LowOrderEmbedding())
          tlp->SetEmbedding (bfa->GetFESpace()->LowOrderEmbedding());
	tlp -> Update();
      }
    else
      tlp = nullptr;

    if (timing) Timing();
    if (test) Test();
    if (mgtest) MgTest();

    dynamic_cast<const BaseSparseMatrix & > (bfa->GetMatrix()).SetInverseType ( invtype );
    if (lo_bfa)
      dynamic_cast<const BaseSparseMatrix & > (lo_bfa->GetMatrix()) .SetInverseType ( loinvtype );
  }
  
  
  void MGPreconditioner :: CleanUpLevel () 
  { 
    if (bfa->GetLowOrderBilinearForm())
      {
	tlp = nullptr;
      }    
  }


  const ngla::BaseMatrix & MGPreconditioner :: GetMatrix() const
  {
    if (tlp)
      return *tlp;
    else
      return *mgp; 
  }

  void MGPreconditioner::PrintReport (ostream & ost) const
  {
    ost << "Multigrid preconditioner" << endl
	<< "bilinear-form = " << bfa->GetName() << endl
	<< "smoothertype = " << smoothertype << endl;
  }


  Array<MemoryUsage> MGPreconditioner::GetMemoryUsage () const
  {
    auto mu = GetMatrix().GetMemoryUsage ();;

    for (int i = 0; i < mu.Size(); i++)
      mu[i].AddName (string(" mgpre ")); // +GetName());
    return mu;
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

    condout << bfa->GetFESpace()->GetNDof() << "\t" << bfa->GetFESpace()->GetOrder() << "\t" << eigen.EigenValue(mgnumber) << "\t" << eigen.MaxEigenValue() << "\t" 
	    << eigen.MaxEigenValue()/eigen.EigenValue(mgnumber) <<  "\t" << endl;
    
    if(testresult_ok) *testresult_ok = eigenretval;
    if(testresult_min) *testresult_min = eigen.EigenValue(mgnumber);
    if(testresult_max) *testresult_max = eigen.MaxEigenValue();

  }

  void MGPreconditioner :: SetDirectSolverCluster(shared_ptr<Array<int>> cluster)
  {
    if(tlp)
      {
        if(auto blocksmoother = dynamic_cast<BlockSmoother*>(&tlp->GetSmoother()))
          blocksmoother->SetDirectSolverCluster(cluster);
      }
    else
      if(auto blocksmoother = dynamic_cast<BlockSmoother*>(&mgp->GetSmoother()))
        blocksmoother->SetDirectSolverCluster(cluster);
  }

  void MGPreconditioner::SetCoarsePreconditioner(shared_ptr<Preconditioner> prec)
  {
      coarse_pre = prec;
      mgp->SetCoarseType (MultigridPreconditioner::USER_COARSE);
  }

  
  // ****************************** DirectPreconditioner **************************


  class NGS_DLL_HEADER DirectPreconditioner : public Preconditioner
  {
    shared_ptr<BilinearForm> bfa;
    shared_ptr<BaseMatrix> inverse;
    string inversetype;

  public:
    DirectPreconditioner (shared_ptr<BilinearForm> abfa, const Flags & aflags,
			  const string aname = "directprecond")
      : Preconditioner(abfa,aflags,aname), bfa(abfa)
    {
      // bfa -> SetPreconditioner (this);
      inversetype = flags.GetStringFlag("inverse", GetInverseName (default_inversetype));
    }

    ///
    virtual ~DirectPreconditioner()
    {
      ; //delete inverse;
    }
    
    virtual void FinalizeLevel (const BaseMatrix * mat) 
    {
      Update();
    }
    
    ///
    virtual void Update ()
    {
      // delete inverse;
      if (GetTimeStamp() == bfa->GetTimeStamp()) return;
      timestamp = bfa->GetTimeStamp();
      
      cout << IM(3) << "Update Direct Solver Preconditioner" << flush;
      
      try
	{                                          
          auto have_sparse_fact = dynamic_pointer_cast<SparseFactorization> (inverse);
          if (have_sparse_fact && have_sparse_fact -> SupportsUpdate())
            {
              if (have_sparse_fact->GetAMatrix() == bfa->GetMatrixPtr())
                {
                  // cout << "have the same matrix, can update factorization" << endl;
                  have_sparse_fact->Update();
                  return;
                }
            }
          
	  bfa->GetMatrix().SetInverseType (inversetype);
	  shared_ptr<BitArray> freedofs = 
	    bfa->GetFESpace()->GetFreeDofs (bfa->UsesEliminateInternal());
	  inverse = bfa->GetMatrix().InverseMatrix(freedofs);
	}
      catch (exception & e)
	{
	  throw Exception (string("caught exception in DirectPreconditioner: \n") +
                           e.what() + 
                           "\nneeds a sparse matrix (or has memory problems)");
	}
      GetMemoryTracer().Track(*inverse, "Inverse");
    }

    virtual void CleanUpLevel ()
    {
      // delete inverse;
      inverse = nullptr;
    }

    virtual const BaseMatrix & GetMatrix() const
    {
      if (!inverse)
        ThrowPreconditionerNotReady();        
      return *inverse;
    }

    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa->GetMatrix(); 
    }

    virtual const char * ClassName() const
    {
      return "Direct Preconditioner"; 
    }
  };







  // ****************************** LocalPreconditioner *******************************


  /**
     Local (Block-Jacobi or Block-Gauss-Seidel) preconditioner
  */
  class LocalPreconditioner : public Preconditioner
  {
  protected:
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<BaseMatrix> jacobi;
    ///
    bool block;
    bool locprectest; 
    string locprecfile; 

    string ct;
    shared_ptr<Preconditioner> coarse_pre;
    function<shared_ptr<Table<DofId>>(FESpace&)> blockcreator;
  public:
    ///
    LocalPreconditioner (shared_ptr<BilinearForm> bfa, const Flags & aflags,
			 const string aname = "localprecond");
    ///
    virtual ~LocalPreconditioner() { ; }
    ///
    virtual bool IsComplex() const { return jacobi->IsComplex(); }
    
    ///
    virtual void FinalizeLevel (const BaseMatrix * mat) 
    {
      cout << IM(3) << "Update Local Preconditioner" << flush;
      timestamp = bfa->GetTimeStamp();
      int blocktype = int (flags.GetNumFlag ( "blocktype", -1));
      
      // if (MyMPI_GetNTasks() != 1) return;
      bool parallel = (this->on_proc == -1);
      /*
        if ( !parallel && id != this->on_proc )
        {
	jacobi = 0; 
	return;
        }
      */

//     BaseBlockJacobiPrecond::COARSE_TYPE bbct;
//     switch ( ct  )
//       {
//       case NO_COARSE:
// 	bbct = NO_COARSE; break;
//       case DIRECT_COARSE:
// 	bbct = DIRECT_COARSE; break;
//       }

      if (blockcreator)
        {
          shared_ptr<Table<int>> blocks = blockcreator(*bfa->GetFESpace());
          // cout << "created blocks: " << *blocks << endl;
          jacobi = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix())
            .CreateBlockJacobiPrecond(blocks, 0, parallel, bfa->GetFESpace()->GetFreeDofs());
          return;
        }
      
      if ( block && blocktype == -1 ) blocktype = 0;
      if ( blocktype >= 0 )
        {
          // new: blocktypes, specified in fespace
          if (bfa->UsesEliminateInternal())
            flags.SetFlag("eliminate_internal");
          shared_ptr<Table<int>> blocks = bfa->GetFESpace()->CreateSmoothingBlocks(flags);
          jacobi = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix())
            .CreateBlockJacobiPrecond(blocks, 0, parallel, bfa->GetFESpace()->GetFreeDofs());
        }
      else if (block)
        {
          cout << "\nFlag block deprecated: use -blocktype=<typeno> instead" << endl;
        }
      else
        {
          shared_ptr<BaseMatrix> mat = bfa->GetMatrixPtr();
#ifdef PARALLEL
          if (dynamic_pointer_cast<ParallelMatrix> (mat))
            mat = dynamic_pointer_cast<ParallelMatrix> (mat)->GetMatrix();
#endif
          jacobi = dynamic_pointer_cast<BaseSparseMatrix> (mat)
            -> CreateJacobiPrecond(bfa->GetFESpace()->GetFreeDofs(bfa->UsesEliminateInternal()));
        }
    }

    virtual void Update ()
    {
      if (GetTimeStamp() < bfa->GetTimeStamp())
        FinalizeLevel (&bfa->GetMatrix());
      if (test) Test();
      if(locprectest) LocPrecTest(); 
    }


    ///
    virtual const BaseMatrix & GetMatrix() const
    {
      if (!jacobi)
        ThrowPreconditionerNotReady();
      return *jacobi;
    }
    
    virtual shared_ptr<BaseMatrix> GetMatrixPtr()
    {
      if (!jacobi)
        ThrowPreconditionerNotReady();
      return jacobi;
    }

    ///
    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa->GetMatrix(); 
    }
    ///
    virtual const char * ClassName() const
    { return "Local Preconditioner"; }
    void LocPrecTest () const;
  };





  LocalPreconditioner :: 
  LocalPreconditioner  (shared_ptr<BilinearForm> abfa, const Flags & aflags,
                        const string aname)

    : Preconditioner (abfa,aflags,aname), bfa(abfa)
  {
    // bfa -> SetPreconditioner (this);

    block = flags.GetDefineFlag ("block");
    locprectest = flags.GetDefineFlag ("mgtest");
    locprecfile = flags.GetStringFlag ("mgfile","locprectest.out"); 

    string smoother = flags.GetStringFlag("smoother","");
    if ( smoother == "block" )
      block = true;

    // coarse-grid preconditioner only used in parallel!!
    ct = "NO_COARSE";

    if(flags.AnyFlagDefined("blockcreator"))
      {
        blockcreator = std::any_cast<function<shared_ptr<Table<DofId>>(const FESpace&)>>(flags.GetAnyFlag("blockcreator"));
        cout << IM(3) << "local pre, got blockcreator" << endl;
      }
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

    condout << bfa->GetFESpace()->GetNDof() << "\t" << bfa->GetFESpace()->GetOrder() << "\t" << eigen.EigenValue(1) << "\t" << eigen.MaxEigenValue() << "\t" 
	    << eigen.MaxEigenValue()/eigen.EigenValue(1) <<  "\t" << endl;
    
    if(testresult_ok) *testresult_ok = eigenretval;
    if(testresult_min) *testresult_min = eigen.EigenValue(1);
    if(testresult_max) *testresult_max = eigen.MaxEigenValue();

  }







  // ****************************** TwoLevelPreconditioner *******************************


  TwoLevelPreconditioner :: ~TwoLevelPreconditioner()
  {
    delete premat;
  }

  const BaseMatrix & TwoLevelPreconditioner :: GetMatrix() const
  {
    return *premat;
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
	cm = new Real2ComplexMatrix<double,Complex> (creal->GetMatrixPtr());
	break;
      case 2:
	cm = new Real2ComplexMatrix<Vec<2,double>,Vec<2,Complex> > (creal->GetMatrixPtr());
	break;
      case 3:
	cm = new Real2ComplexMatrix<Vec<3,double>,Vec<3,Complex> > (creal->GetMatrixPtr());
	break;
      case 4:
	cm = new Real2ComplexMatrix<Vec<4,double>,Vec<4,Complex> > (creal->GetMatrixPtr());
	break;
      default:
	cout << "Error: dimension " << dim << " for complex preconditioner not supported!" << endl;
      }
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


  CommutingAMGPreconditioner :: ~CommutingAMGPreconditioner ()
  {
    delete amg;
  }

  void CommutingAMGPreconditioner :: Update ()
  {
    cout << "Update amg" << endl;

#ifdef PARALLEL
    auto comm = ma->GetCommunicator();
    if ( this->on_proc != comm.Rank() && this->on_proc != -1)
      {
	amg = NULL;
	return;
      }
#endif

    shared_ptr<MeshAccess> ma = nullptr; //pde->GetMeshAccess();

    size_t nedge = ma->GetNEdges(); 
    size_t nface = ma->GetNFaces(); 
    size_t nel = ma->GetNE();

    if (coarsegrid && ma->GetNLevels() > 1)
      return;
      

    delete amg;

    // cout << "get edges" << endl;

    Array<INT<2> > e2v (nedge);
    e2v = INT<2> (-1, -1);
    for (size_t i = 0; i < nedge; i++)
      if (ma->GetClusterRepEdge (i) >= 0)
	e2v[i] = ma->GetEdgePNums (i);

    // cout << "get faces" << endl;

    Array<INT<4> > f2v (nface);
    for (int i = 0; i < nface; i++)
      {
	if (ma->GetClusterRepFace (i) >= 0)
	  {
	    auto fpnums = ma->GetFacePNums (i);
	    
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
      if (e2v[i][0] != -1)
	{
	  Vec<3> v = 
	    ma->GetPoint<3> (e2v[i][0]) -
	    ma->GetPoint<3> (e2v[i][1]);
	  hi[i] = L2Norm (v);
	}
    
    if (hcurl)
      {
	ai.SetSize(nface);
	for (int i = 0; i < nface; i++)
	  {
	    Vec<3> p1 = ma->GetPoint<3> (f2v[i][0]);
	    Vec<3> p2 = ma->GetPoint<3> (f2v[i][1]);
	    Vec<3> p3 = ma->GetPoint<3> (f2v[i][2]);
	    Vec<3> vn = Cross (Vec<3> (p2-p1), Vec<3> (p3-p1));
	    ai[i] = L2Norm (vn);
	  }
      }

    // Array<int> fanums(12);
    LocalHeap lh (10000, "CommutingAMG");
    IntegrationPoint ip(0, 0, 0, 0);

    cout << "compute weights" << endl;

    weighte = 0;
    weightf = 0;
    for (int i = 0; i < nel; i++)
      {
        ElementId ei(VOL,i);
	HeapReset hr(lh);
	auto ednums = ma->GetElEdges (ei);

	ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
	MappedIntegrationPoint<3,3> sip(ip, eltrans);

	double vol = ma->ElementVolume (i);
	double vale = ngfem::Evaluate (*coefe, sip);

	for (int j = 0; j < ednums.Size(); j++)
	  weighte[ednums[j]] += vale * vol / sqr (hi[ednums[j]]);

	if (hcurl)
	  {
	    auto fanums = ma->GetElFaces (ei);
	    double valf = ngfem::Evaluate (*coeff, sip);
	    for (int j = 0; j < fanums.Size(); j++)
	      weightf[fanums[j]] += valf * vol / sqr (ai[fanums[j]]);
	  }
      }

    int nsel = ma->GetNSE();
    if (coefse)
      for (int i = 0; i < nsel; i++)
	{
          ElementId sei(BND,i);
	  HeapReset hr(lh);
	  auto ednums = ma->GetElEdges (sei);
	  ElementTransformation & eltrans = ma->GetTrafo (sei, lh);

	  MappedIntegrationPoint<2,3> sip(ip, eltrans);

	  double vol = ma->SurfaceElementVolume (i);
	  double vale = ngfem::Evaluate (*coefse, sip);

	  for (int j = 0; j < ednums.Size(); j++)
	    weighte[ednums[j]] += vale * vol / sqr (hi[ednums[j]]);
	}


    Timer timer_coarse("AMG - Coarsening time");

    timer_coarse.Start();
    Array< Vec<3> > vertices;

    int nv = ma->GetNV();
    vertices.SetSize(nv);
    for (int i = 0; i < nv; i++)
      ma->GetPoint(i, vertices[i]);

    CommutingAMG * amgmat;
    if (hcurl)
      amgmat = new AMG_HCurl (bfa->GetMatrix(), vertices, e2v, f2v, weighte, weightf, levels);
    else
      amgmat = new AMG_H1 (bfa->GetMatrix(), e2v, weighte, levels);
    timer_coarse.Stop();

    cout << "AMG coarsening time = " << timer_coarse.GetTime() << " sec" << endl;

    Timer timer_proj("AMG - projection time");
    timer_proj.Start();
    amgmat->ComputeMatrices (dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix()));
    timer_proj.Stop();

    cout << "AMG projection time = " << timer_proj.GetTime() << " sec" << endl;
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
                      function<shared_ptr<Preconditioner>(shared_ptr<BilinearForm>,const Flags &,const string)> acreatorbf,
                      DocInfo adocinfo)
    : name(aname), creatorbf(acreatorbf), docinfo(adocinfo)
  {
    ;
  }

  /*
  PreconditionerClasses :: PreconditionerClasses () { }
  PreconditionerClasses :: ~PreconditionerClasses() { }
  */
  
  void PreconditionerClasses :: 
  AddPreconditioner (const string & aname,
                     function<shared_ptr<Preconditioner>(shared_ptr<BilinearForm>,const Flags &,const string)> acreatorbf,
                     DocInfo docinfo)
  {
    prea.Append (make_unique<PreconditionerInfo>(aname, acreatorbf, docinfo));
  }

  const PreconditionerClasses::PreconditionerInfo * 
  PreconditionerClasses::GetPreconditioner(const string & name)
  {
    for (int i = 0; i < prea.Size(); i++)
      if (name == prea[i]->name)
        return prea[i].get();
    return nullptr;
  }

  void PreconditionerClasses :: Print (ostream & ost) const
  {
    ost << endl << "Preconditioners:" << endl;
    ost <<         "---------" << endl;
    ost << setw(20) << "Name" << endl;
    for (int i = 0; i < prea.Size(); i++)
      ost << setw(20) << prea[i]->name << endl;
  }

  void PreconditionerClasses :: Cleanup ()
  {
      prea.DeleteAll();
  }

 
  PreconditionerClasses & GetPreconditionerClasses ()
  {
    static PreconditionerClasses fecl;
    return fecl;
  }

  RegisterPreconditioner<MGPreconditioner> registerMG("multigrid");
  RegisterPreconditioner<DirectPreconditioner> registerDirect("direct");
  RegisterPreconditioner<LocalPreconditioner> registerlocal("local");

}




