#ifdef LAPACK

#include <solve.hpp>


/*

  Solves for the generalized EVP 

  A u = lam M u

  Uses a simultaneous preconditioned inverse iteration

 */


namespace ngsolve
{
  class NumProcEVP_AM : public NumProc
  {
  protected:
    shared_ptr<BilinearForm> bfa;
    shared_ptr<BilinearForm> bfm;
    shared_ptr<GridFunction> gfu;
    shared_ptr<Preconditioner> pre;

    int maxsteps, nr, maxnewton;
    double prec;
    bool print;
    bool coarse;
    string variable;

  public:
    NumProcEVP_AM (shared_ptr<PDE> apde, const Flags & flags);
    virtual ~NumProcEVP_AM();

    virtual void Do(LocalHeap & lh);

    virtual string GetClassName () const
    {
      return "Eigenvalue Problem";
    }

    virtual void PrintReport (ostream & ost) const
    {
      ost << GetClassName() << endl
	  << "Bilinear-form A = " << bfa->GetName() << endl
	  << "Bilinear-form M = " << bfm->GetName() << endl
	  << "Gridfunction  = " << gfu->GetName() << endl
	  << "Preconditioner = " << ((pre) ? pre->ClassName() : "None") << endl;
    }
    
    static void PrintDoc (ostream & ost)
    {
      ost << 
	"\n\nNumproc evpAM:\n"			\
	"--------------\n"						\
	"Solves a generalized symmetric eigenvalue problem\n\n" \
	"Required flags:\n" 
	"-bilinearforma=<bfname>\n" 
	"    bilinear-form providing the stiffness matrix\n"	\
	"-bilinearformm=<bfname>\n" 
	"    bilinear-form providing the mass matrix\n"	\
	"-gridfunction=<gfname>\n" 
	"    a gridfunction. multidim=xx defines num of calculated ev\n"	\
	"-preconditioner=<prename>\n" \
	"    a preconditioner for the stiffness matrix\n" \
	"-maxsteps=n\n" 
	  << endl;
    }
  };


  
  NumProcEVP_AM :: NumProcEVP_AM (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearforma", ""));
    bfm = apde->GetBilinearForm (flags.GetStringFlag ("bilinearformm", ""));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
    pre = apde->GetPreconditioner (flags.GetStringFlag ("preconditioner", ""));
    maxsteps = int(flags.GetNumFlag ("maxsteps", 200));
    variable = flags.GetStringFlag ("variable", "eigenvalue");

    maxnewton = int(flags.GetNumFlag ("maxnewton", 0));
    nr = int(flags.GetNumFlag ("nr", 0));
  }

  NumProcEVP_AM :: ~NumProcEVP_AM()
  {
    ;
  }


  void NumProcEVP_AM :: Do(LocalHeap & lh)
  {
    cout << "solve  evp" << endl;

    // simultaneous iteration without coarse grid

    const BaseMatrix & mata = bfa->GetMatrix();
    const BaseMatrix & matm = bfm->GetMatrix();
    const BaseMatrix & matpre = pre->GetMatrix();

    int num = gfu -> GetMultiDim();
    cout << "num = " << num << endl;
	

    BaseVector & vecu = gfu->GetVector(0);
	
    auto d = vecu.CreateVector();
    auto w = vecu.CreateVector();
    auto au = vecu.CreateVector();
    auto mu = vecu.CreateVector();
    
    Vector<double> lami(2*num);

    
    cout.precision(16);

    Array<shared_ptr<BaseVector>> hvec(2*num);
    for (int i = 0; i < 2*num; i++)
      hvec[i] = vecu.CreateVector();

    // random initial guess, in the range of the precond
    for (int i = 0; i < num; i++)
      {
	hvec[0]->SetRandom();
	gfu->GetVector(i) = matpre * *hvec[0];
      }
	

    Matrix<double> ha(2*num, 2*num);
    Matrix<double> hm(2*num, 2*num);
    Matrix<double> evecs(2*num, 2*num);


    for (int it = 0; it < maxsteps; it++)
      {
	cout << "it = " << it << " " << endl;
	
	for (int ei = 0; ei < num; ei++)
	  {
	    BaseVector & vecu = gfu->GetVector (ei);

	    *au = mata * vecu;
	    *mu = matm * vecu;
	    
	    double lam = 
	      InnerProduct (vecu, *au) /
	      InnerProduct (vecu, *mu);
	    
	    *d = *au - lam * *mu;
	    *w = matpre * *d;
	    double norm = InnerProduct (*w, *mu = matm* *w);

	    *hvec[ei] = vecu;
	    *hvec[num+ei] = (1.0 / sqrt(norm)) * *w;
	  }
	    
	    
	for (int ei = 0; ei < 2*num; ei++)
	  {
	    *au = mata * *hvec[ei];
	    *mu = matm * *hvec[ei];
		
	    for (int ej  = 0; ej < 2*num; ej++)
	      {
		ha(ei, ej) = InnerProduct (*au, *hvec[ej]);
		hm(ei, ej) = InnerProduct (*mu, *hvec[ej]);
	      } 
	  }

        cout << "ha = " << endl << ha << endl;
        cout << "hm = " << endl << hm << endl;

#ifdef LAPACK
	LapackEigenValuesSymmetric (ha, hm, lami, evecs);
#else
        cerr << "This numproc needs LAPACK" << endl;
#endif
	
	for (int ei = 0; ei < num; ei++)
	  cout << "lam(" << ei << ") = " << lami(ei) << endl;
	    
	for (int ei = 0; ei < num; ei++)
	  {
	    *au = 0;
	    for (int i = 0; i < 2*num; i++)
	      *au += evecs(ei, i) * *hvec[i];
	    gfu->GetVector (ei) = *au;
	  }
      }

    
    for (int i = 0; i < num; i++)
      {
        stringstream vn;
        vn << variable << i;
        shared_ptr<PDE> (pde)->AddVariable (vn.str(), lami(i));
      }










#ifdef A_MULTIGRID_EVP_SOLVER


    if (!bfa->GetFESpace().IsComplex())
      {
	int level = ma.GetNLevels()-1;
	cout << "level = " << level;
	
	if (level <= 0) return;
	
	const BilinearForm * bfac = 
	  bfa->HasLowOrderBilinearForm() ?  
	  &bfa->GetLowOrderBilinearForm() : bfa;

	const BilinearForm * bfmc = 
	  bfm->HasLowOrderBilinearForm() ? 
	  &bfm->GetLowOrderBilinearForm() : bfm;

	const SparseMatrix<double> & mata = 
	  dynamic_cast<const SparseMatrix<double>&> (bfa->GetMatrix());
	const SparseMatrix<double> & matac = 
	  dynamic_cast<const SparseMatrix<double>&> (bfac->GetMatrix(0));
	const SparseMatrix<double> & matm = 
	  dynamic_cast<const SparseMatrix<double>&> (bfm->GetMatrix());
	const SparseMatrix<double> & matmc = 
	  dynamic_cast<const SparseMatrix<double>&> (bfmc->GetMatrix(0));
	
	BaseVector & vecu = gfu->GetVector();
	
	const BaseMatrix & matpre = pre->GetMatrix();
	
	// int n = matac.Height();
	int nc = matac.Height();
	// nc = 0;
	
	int i, j, k, l;
	
	BaseVector & d = *vecu.CreateVector();
	BaseVector & w = *vecu.CreateVector();
	BaseVector & au = *vecu.CreateVector();
	BaseVector & mu = *vecu.CreateVector();
	
	BaseVector & auc = *matac.CreateVector();
	BaseVector & muc = *matac.CreateVector();
	
	Vector<double> lami(nc+2*num);
	
	
	vecu.SetRandom();
	
	const Prolongation & prol = *gfu->GetFESpace().GetProlongation();
	
	cout.precision(16);
	
	Array<BaseVector*> vecui(2*num);
	Array<BaseVector*> vecunew(num);
	
	
	for (i = 0; i < num; i++)
	  {
	    stringstream strname;
	    strname << "gfu" << i;
	    string name = strname.str();
	    
	    Flags flags;
	    GridFunction * gf = CreateGridFunction (&gfu->GetFESpace(), name, flags);
	    gf->Update();
	    gf->Visualize (name);
	    
	    vecui[i] = &gf->GetVector(); // vecu.CreateVector();
	    vecui[i] -> SetRandom();
	  }
	
	for (i = num; i < 2*num; i++)
	  {
	    vecui[i] = vecu.CreateVector();
	    vecui[i] -> SetRandom();
	  }
	
	for (i = 0; i < num; i++)
	  {
	    vecunew[i] = vecu.CreateVector();
	  }
	
	Matrix<double> ha(nc+2*num, nc+2*num);
	Matrix<double> hm(nc+2*num, nc+2*num);


	/*
	FlatVector<> fu = vecu.FVDouble();
	if (nc > 0)
	  {
	    Matrix<> matac2(nc);
	    Matrix<> matmc2(nc);
	    for (int i = 0; i < nc; i++)
	      {
		fu = 0;
		fu(i) = 1;
		for (l = 1; l <= level; l++)
		  prol.ProlongateInline (l, vecu);
		au = mata * vecu;
		mu = matm * vecu;
		
		for (int j = 0; j < nc; j++)
		  {
		    fu = 0;
		    fu(j) = 1;
		    for (l = 1; l <= level; l++)
		      prol.ProlongateInline (l, vecu);
		    
		    matac2(i,j) = InnerProduct (au, vecu);
		    matmc2(i,j) = InnerProduct (mu, vecu);
		  }
	      }
	    (*testout) << "matac = " << endl << matac << endl;
	    (*testout) << "matac2 = " << endl << matac2 << endl;
	    (*testout) << "matmc = " << endl << matmc << endl;
	    (*testout) << "matmc2 = " << endl << matmc2 << endl;

	    for (int i = 0; i < nc; i++)
	      for (int j = 0; j < nc; j++)
		{
                  matac2(i,j) -= matac(i,j);
		  matmc2(i,j) -= matmc(i,j);
		}

	    (*testout) << "diffa = " << endl << matac2 << endl;
	    (*testout) << "diffm = " << endl << matmc2 << endl;

	  }
	*/


	
	// large space 
	for (int it = 0; it < maxsteps; it++)
	  {
	    cout << "it = " << it << " " << endl;
	    
	    // coarse grid matrix
	    for (i = 0; i < nc; i++)
	      for (j = 0; j < nc; j++)
		{
		  ha(i,j) = matac(i,j);
		  hm(i,j) = matmc(i,j);
		}
	    
	    
	    for (int ei = 0; ei < num; ei++)
	      {
		au = mata * *vecui[ei];
		mu = matm * *vecui[ei];
		
		double lam = 
		  InnerProduct (*vecui[ei], au) /
		  InnerProduct (*vecui[ei], mu);
		
		d = au - lam * mu;
		w = matpre * d;
		double norm = InnerProduct (w, mu=matm*w);
		*vecui[num+ei] = (1.0 / sqrt(norm)) * w;
	      }
	    
	    
	    for (int ei = 0; ei < 2*num; ei++)
	      {
		au = mata * *vecui[ei];
		mu = matm * *vecui[ei];
		
		for (int ej  = 0; ej < 2*num; ej++)
		  {
		    ha(nc+ei, nc+ej) = InnerProduct (au, *vecui[ej]);
		    hm(nc+ei, nc+ej) = InnerProduct (mu, *vecui[ej]);
		  }
		
		for (l = level; l >= 1; l--)
		  {
		    prol.RestrictInline (l, au);
		    prol.RestrictInline (l, mu);
		  }
		
		FlatVector<double> fau = au.FVDouble();
		FlatVector<double> fmu = mu.FVDouble();
		
		for (i = 0; i < nc; i++)
		  {
		    ha(i, nc+ei) = ha(nc+ei, i) = fau(i);
		    hm(i, nc+ei) = hm(nc+ei, i) = fmu(i);
		  }
	      }
	    
	    (*testout) << "ha = " << endl << ha << endl;
	    (*testout) << "hm = " << endl << hm << endl;
	   
#ifdef LAPACK 
	    // cout << "lapack ..." << flush;
	    LapackGHEPEPairs (nc+2*num, &ha(0,0), &hm(0,0), &lami(0));
	    // cout << " done " << flush;
#else
            cerr << "This numproc needs LAPACK" << endl;	    
#endif


	    for (int ei = 0; ei < num; ei++)
	      cout << "lam(" << ei << ") = " << lami(ei) << endl;
	    
	    for (int ei = 0; ei < num; ei++)
	      {
		FlatVector<double> fau = au.FVDouble();
		for (i = 0; i < nc; i++)
		  fau(i) = ha(ei,i);
		
		if (nc > 0)
		  for (l = 1; l <= level; l++)
		    prol.ProlongateInline (l, au);
		else
		  au = 0;
		
		for (i = 0; i < 2*num; i++)
		  au += ha(ei,nc+i) * *vecui[i];
		
		*vecunew[ei] = au;
	      }
	    
	    for (int ei = 0; ei < num; ei++)
	      *vecui[ei] = *vecunew[ei];



	    /*
	    // coarse space, u+w

	    Matrix<double> ha1(nc+num, nc+num);
	    Matrix<double> hm1(nc+num, nc+num);
	    Matrix<double> ha2(2*num, 2*num);
	    Matrix<double> hm2(2*num, 2*num);

	    for (int it = 0; it < maxsteps; it++)
	    {
	    cout << "it = " << it << " ";

	    // coarse grid matrix
	    for (i = 0; i < nc; i++)
	    for (j = 0; j < nc; j++)
	    {
	    ha1(i,j) = matac(i,j);
	    hm1(i,j) = matmc(i,j);
	    }

	    for (int ei = 0; ei < num; ei++)
	    {
	    au = mata * *vecui[ei];
	    mu = matm * *vecui[ei];
	    
	    for (int ej  = 0; ej < num; ej++)
	    {
	    ha1(nc+ei, nc+ej) = InnerProduct (au, *vecui[ej]);
	    hm1(nc+ei, nc+ej) = InnerProduct (mu, *vecui[ej]);
	    }
	    
	    for (l = level; l >= 1; l--)
	    {
	    prol.RestrictInline (l, au);
	    prol.RestrictInline (l, mu);
	    }

	    FlatVector<double> fau = au.FVDouble();
	    FlatVector<double> fmu = mu.FVDouble();

	    for (i = 0; i < nc; i++)
	    {
	    ha1(i, nc+ei) = ha1(nc+ei, i) = fau(i);
	    hm1(i, nc+ei) = hm1(nc+ei, i) = fmu(i);
	    }
	    }

#ifdef LAPACK
	    cout << "lapack ..." << flush;
	    LapackGHEPEPairs (nc+num, &ha1(0,0), &hm1(0,0), &lami(0));
	    cout << " done " << flush;
#else
            cout << "This numproc needs LAPACK!" << endl;
#endif

	
	    for (int ei = 0; ei < num; ei++)
	    cout << "lam(" << ei << ") = " << lami(ei) << endl;
	
	    for (int ei = 0; ei < num; ei++)
	    {
	    FlatVector<double> fau = au.FVDouble();
	    for (i = 0; i < nc; i++)
	    fau(i) = ha1(ei,i);

	    for (l = 1; l <= level; l++)
	    prol.ProlongateInline (l, au);

	    for (i = 0; i < num; i++)
	    au += ha1(ei,nc+i) * *vecui[i];

	    *vecunew[ei] = au;
	    }

	    for (int ei = 0; ei < num; ei++)
	    *vecui[ei] = *vecunew[ei];



	    // u+w
	    for (int ei = 0; ei < num; ei++)
	    {
	    au = mata * *vecui[ei];
	    mu = matm * *vecui[ei];

	    double lam = 
	    InnerProduct (*vecui[ei], au) /
	    InnerProduct (*vecui[ei], mu);
	      
	    d = au - lam * mu;
	    w = matpre * d;
	    *vecui[num+ei] = 1.0 / sqrt (InnerProduct (d,w)) * w;
	    }


	    for (int ei = 0; ei < 2*num; ei++)
	    {
	    au = mata * *vecui[ei];
	    mu = matm * *vecui[ei];
	    
	    for (int ej  = 0; ej < 2*num; ej++)
	    {
	    ha2(ei, ej) = InnerProduct (au, *vecui[ej]);
	    hm2(ei, ej) = InnerProduct (mu, *vecui[ej]);
	    }
	    }

	    (*testout) << "ha2 = " << endl << ha2 
	    << "hm2 = " << endl << hm2 << endl;

#ifdef LAPACK
	    cout << "lapack ..." << flush;
	    LapackGHEPEPairs (2*num, &ha2(0,0), &hm2(0,0), &lami(0));
	    cout << " done " << flush;
#else
            cerr << "This numproc needs LAPACK" << endl;	
#endif

	    for (int ei = 0; ei < num; ei++)
	    cout << "lam(" << ei << ") = " << lami(ei) << endl;
	
	    for (int ei = 0; ei < num; ei++)
	    {
	    FlatVector<double> fau = au.FVDouble();
	    for (i = 0; i < nc; i++)
	    fau(i) = ha2(ei,i);

	    au = 0;
	    for (i = 0; i < 2*num; i++)
	    au += ha2(ei,i) * *vecui[i];

	    *vecunew[ei] = au;
	    }

	    for (int ei = 0; ei < num; ei++)
	    *vecui[ei] = *vecunew[ei];
	    }
	    */




	    /*
	      Matrix<double> ha2(2,2);
	      Matrix<double> hm2(2,2);

	      for (int it = 0; it < maxsteps; it++)
	      {
	      cout << "it = " << it << " ";
	      au = mata * vecu;
	      mu = matm * vecu;

	      double uau = InnerProduct (au, vecu);
	      double umu = InnerProduct (mu, vecu);
	
	      double lam = uau / umu;

	      cout << "lam = " << lam << endl;
	      (*testout) << "lam = " << lam << endl;
	      (*testout) << ", uau = " << uau << ", umu = " << umu;

	      d = au - lam * mu;
	      w = matpre * d;

	      ha2(0,0) = uau;
	      hm2(0,0) = umu;

	      au = mata * w;
	      mu = matm * w;

	      ha2(1,0) = ha2(0,1) = InnerProduct (au, vecu);
	      hm2(1,0) = hm2(0,1) = InnerProduct (mu, vecu);
	      ha2(1,1) = InnerProduct (au, w);
	      hm2(1,1) = InnerProduct (mu, w);

	      (*testout) << ", waw = " << ha2(1,1) << ", wmw = " << hm2(1,1) << endl;

	      if (hm2(1,1) < 1e-18) break;

#ifdef LAPACK
	      cout << "lapack ..." << flush;
	      LapackGHEPEPairs (2, &ha2(0,0), &hm2(0,0), &lami(0));
	      cout << " done " << flush;
#else
              cerr << "This numproc needs LAPACK! " << endl;
#endif
	
	      vecu *= ha2(0,0);
	      vecu += ha2(0,1) * w;
	      }
	    */
	  }
      }
#endif

  }


  static RegisterNumProc<NumProcEVP_AM> npevpam("evpAM");
}
#endif
