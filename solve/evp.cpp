#ifdef LAPACK

#include "../include/solve.hpp"


namespace ngsolve
{


  /* ***************************** Numproc EVP *************************** */

  ///
  class NumProcEVP : public NumProc
  {
  protected:
    shared_ptr<BilinearForm> bfa;
    shared_ptr<BilinearForm> bfm;

    shared_ptr<GridFunction> gfu;
    shared_ptr<Preconditioner> pre;
    int num;

    double prec, shift, shifti;
    bool print;

    string filename;

    enum SOLVER { DENSE, ARNOLDI };
    SOLVER solver;

  public:
    NumProcEVP (shared_ptr<PDE> apde, const Flags & flags);
    virtual ~NumProcEVP() { ; }

    virtual void Do(LocalHeap & lh);

    virtual string GetClassName () const
    {
      return " Eigenvalue Solver";
    }

    virtual void PrintReport (ostream & ost) const
    {
      ost << GetClassName() << endl
	  << "Bilinear-form A = " << bfa->GetName() << endl
	  << "Bilinear-form M = " << bfm->GetName() << endl
	  << "Gridfunction  = " << gfu->GetName() << endl;
    }
  };


  NumProcEVP :: NumProcEVP (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    shared_ptr<PDE> spde (pde);
    bfa = spde->GetBilinearForm (flags.GetStringFlag ("bilinearforma", ""));
    bfm = spde->GetBilinearForm (flags.GetStringFlag ("bilinearformm", ""));
    gfu = spde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
    pre = spde->GetPreconditioner (flags.GetStringFlag ("preconditioner", ""),1);
    num = int(flags.GetNumFlag ("num", 500));
    shift = flags.GetNumFlag ("shift",1); 
    shifti = flags.GetNumFlag ("shifti",0); 

    filename = flags.GetStringFlag ("filename","eigen.out"); 

    solver = ARNOLDI;
    if (flags.GetDefineFlag("dense")) solver = DENSE;
  }



  void NumProcEVP :: Do(LocalHeap & lh)
  {
    if (solver == ARNOLDI)
      {
        cout << "new version using linalg - Arnoldi" << endl;

        if (bfa->GetFESpace()->IsComplex())
          {
            cout << "complex evp" << endl;
            
            Arnoldi<Complex> arnoldi (bfa->GetMatrixPtr(), bfm->GetMatrixPtr(), 
                                      bfa->GetFESpace()->GetFreeDofs() );
            arnoldi.SetShift (Complex(shift,shifti));
            
            int nev = gfu->GetMultiDim();
            Array<shared_ptr<BaseVector>> evecs(nev);

            Array<Complex> lam(nev);
            if (pre)
              arnoldi.Calc (num, lam, nev, evecs, pre->GetMatrixPtr());
            else
              arnoldi.Calc (num, lam, nev, evecs, 0);
            
            ofstream eigenout(filename.c_str());
            eigenout.precision(16);
            for (int i = 0; i < lam.Size(); ++i)
              eigenout << lam[i].real() << "\t" << lam[i].imag() << "\t"
                       << sqrt(lam[i]).real() << "\t" << sqrt(lam[i]).imag() 
                       << endl;
            
            for (int i = 0; i < nev; i++)
              gfu->GetVector(i) = *evecs[i];

            cout << "lam = " << endl << lam << endl;
          }
        else
          {
            cout << "real evp" << endl;
            Arnoldi<double> arnoldi (bfa->GetMatrixPtr(), bfm->GetMatrixPtr(), 
                                     bfa->GetFESpace()->GetFreeDofs() );
            arnoldi.SetShift (shift);
            
            int nev = gfu->GetMultiDim();
            Array<shared_ptr<BaseVector>> evecs(nev);
            // for (int i = 0; i  < nev; i++)
            // evecs[i] = &gfu->GetVector(i);
            Array<Complex> lam(nev);
            if (pre)
              arnoldi.Calc (num, lam, nev, evecs, pre->GetMatrixPtr());
            else
              arnoldi.Calc (num, lam, nev, evecs, 0);
            
            ofstream eigenout(filename.c_str());
            eigenout.precision(16);
            for (int i = 0; i < lam.Size(); ++i)
              eigenout << lam[i].real() << "\t" << lam[i].imag() << endl;
            
            for (int i = 0; i < nev; i++)
              gfu->GetVector(i) = *evecs[i];

            cout << "lam = " << endl << lam << endl;
          }
        
        return;
      }

    int dim = bfa->GetFESpace()->GetDimension();
    int size = bfa->GetMatrix().Height();
    bool iscomplex = bfa->GetFESpace()->IsComplex();

    
    if (solver == DENSE)
      {
        cout << "solve evp with dense Lapack" << endl;

	int n = dim*size;

	if (!iscomplex)
	  {
	    VVector<> hv(n), hva(n), hvm(n);
	    Matrix<> mata(n), matm(n);
	    
	    for (int i = 0; i < n; i++)
	      {
		hv = 0.0;
		hv(i) = 1.0;
		hva = bfa->GetMatrix() * hv;
		hvm = bfm->GetMatrix() * hv;
		
		for (int j = 0; j < n; j++)
		  {
		    mata(j,i) = hva(j);
		    matm(j,i) = hvm(j);
		  }
	      }
	    cout << "has matrices" << endl;
	    
	    Vector<> lami(n);
	    Matrix<> evecs(n);
	    
	    LapackGHEPEPairs (n, &mata(0,0), &matm(0,0), &lami(0));
	    cout.precision(12);
	    cout << "lami = " << lami << endl;
	    
	    // LaEigNSSolve (n, &mata(0,0), &matm(0,0), &lami(0), 1, &evecs(0,0), 0, 'N');
	    
	    cout << "eigensystem done" << endl;
	  }
	else
	  {
	    cout << "dense complex solver" << endl;
	    VVector<Complex> hv(n), hva(n), hvm(n);
	    Matrix<Complex> mata(n), matm(n);
	    
	    for (int i = 0; i < n; i++)
	      {
		hv = 0.0;
		hv(i) = 1.0;
		hva = bfa->GetMatrix() * hv;
		hvm = bfm->GetMatrix() * hv;
		
		for (int j = 0; j < n; j++)
		  {
		    mata(j,i) = hva(j);
		    matm(j,i) = hvm(j);
		  }
	      }
	    cout << "has matrices" << endl;
	    
	    *testout << "mata = " << mata << endl;
	    *testout << "matm = " << matm << endl;

	    Vector<Complex> lami(n);
	    Matrix<Complex> evecs(n);
	    
            cout << "calling lapack ... " << flush;
	    LaEigNSSolve (n, &mata(0,0), &matm(0,0), &lami(0), 1, &evecs(0,0), 0, 'N');
            cout << "done" << endl;

	    for (int i = 0; i < min2 (gfu->GetMultiDim(), n); i++)
	      {
		FlatVector<Complex> vu = gfu->GetVector(i).FVComplex();
		for (int j = 0; j < n; j++)
		  vu(j) = evecs(i,j);
	      }


	    cout.precision(12);
	    cout << "lami = " << lami << endl;

	    ofstream ev ("eigenvalues.out");
	    for (int i = 0; i < n; i++)
	      {
		Complex lam = lami(i);
		ev << i << " " << lam.real() << " " << lam.imag();
                lam = sqrt(lami(i));
		ev  << " " << lam.real() << " " << lam.imag() << endl;
	      }
	    cout << "eigensystem done" << endl;
	  }
      }
    
    else if (solver == ARNOLDI)

      {
        //  old version 
	int n = dim*size;
	// VVector<Complex> hv(n), hva(n), hvm(n), hwa(n), hwm(n), w(n);

	// Arnoldi:
	int m = min2 (num, n);

	Matrix<Complex> matV (m,n), matVt(n,m);
	Matrix<Complex> matH(m), matHt(m), matI(m);
	
	const BaseSparseMatrix&  mata1 =
	  dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix());
	
	BaseSparseMatrix&  mata =
	  dynamic_cast<BaseSparseMatrix&> (*mata1.CreateMatrix());
	
	const BaseSparseMatrix& matm = 
	  dynamic_cast<const BaseSparseMatrix&> (bfm->GetMatrix()); 

        
        BaseVector & hv = *mata.CreateColVector();
        BaseVector & hva = *mata.CreateColVector();
        BaseVector & hvm = *mata.CreateColVector();
        BaseVector & hwa = *mata.CreateColVector();
        BaseVector & hwm = *mata.CreateColVector();
        BaseVector & w = *mata.CreateColVector();

        FlatVector<Complex> hvc = hv.FVComplex();
        FlatVector<Complex> hvmc = hvm.FVComplex();


	mata.AsVector() = mata1.AsVector() - shift*matm.AsVector();  
	auto inva = mata.InverseMatrix(bfa->GetFESpace()->GetFreeDofs());	    
	
	// BaseMatrix * inva = dynamic_cast<const BaseSparseMatrix&> (bfa->GetMatrix()).InverseMatrix();

	hv.SetRandom();
	matV = Complex(0.0);
	matH = Complex(0.0);
	matI = Complex(0.0);
	for (int i = 0; i < m; i++)
	  matI(i,i) = 1.0;

	Complex len = sqrt (S_InnerProduct<Complex> (hv, hv));
	hv /= len;


	for (int i = 0; i < m; i++)
	  {
	    cout << "i = " << i << endl;
	    for (int j = 0; j < n; j++)
	      matV(i,j) = hvc(j);
	    
	    hva = bfm->GetMatrix() * hv;
	    hvm = *inva * hva;

	    for (int j = 0; j <= i; j++)
	      {
		Complex sum = 0.0;
		for (int k = 0; k < n; k++)
		  sum += hvmc(k) * matV(j,k);
		matH(j,i) = sum;
		for (int k = 0; k < n; k++)
		  hvmc(k) -= sum * matV(j,k);
	      }
	    
	    hv = hvm;
	    Complex len = sqrt (S_InnerProduct<Complex> (hv, hv));
	    if (i<m-1)
	      matH(i+1,i) = len; 

	    //  cout << ", len = " << len << endl;
	    hv /= len;
	  }

	
	//	(*testout) << "matH = " << endl << matH << endl;
	//	(*testout) << "V^T V = " << endl << matV * Trans (matV) << endl;
	cout << "has Hessenberg" << endl;

	Vector<Complex> lami(m);
	Matrix<Complex> evecs(m);
	Matrix<Complex> levecs(n,m);
	matHt = Trans (matH);

	evecs = Complex (0.0);
	lami = Complex (0.0);

	// LaEigNSSolve (mata.Height(), &mata(0,0), &matm(0,0), &lami(0), 1, &evecs(0,0), 0, 'N');
	// LaEigNSSolve (matH.Height(), &matHt(0,0), &matI(0,0), &lami(0), 1, &evecs(0,0), 0, 'N');

	LapackHessenbergEP (matH.Height(), &matHt(0,0), &lami(0), &evecs(0,0));

	cout << "hessenberg is back" << endl;

	// levecs = Trans(matV) * Trans(evecs);
	matVt = Trans (matV);
	levecs = Trans(evecs * Trans (matVt));
	
	for (int i = 0; i < m; i++)
	  lami(i) =  1.0 / lami(i) + shift;

	cout << "eigensystem done" << endl;


	// sort by angle
	BitArray used(m);
	Array<int> reorder(m);

	used.Set();
	for (int i = 0; i < m; i++)
	  reorder[i] = i;
	/*
	used.Clear();
	for (int i = 0; i < m; i++)
	  {
	    int mini = -1;
	    double val = 0;
	    for (int j = 0; j < m; j++)
	      if (!used[j] && lami(j).real() < 10000 && lami(j) != Complex(100,100))
		{
		  if (mini < 0 ||
		      lami(j).imag() / lami(j).real() > val)
		    {
		      mini = j;
		      val = lami(j).imag() / lami(j).real();
		    }
		}
	    if (mini != -1)
	      {
		reorder.Append (mini);
		used.Set(mini);
	      }
	  }
	for (int i = 0; i < m; i++)
	  {
	    if (!used[i] && lami(i) != Complex(100,100))
	      reorder.Append(i);
	  }
	*/

        

	ofstream out ("eigenvalues.dat");
	//	for (int i = 0; i < reorder.Size(); i++)
	//	  out << i << " " << lami(reorder[i]).real() << " " << lami(reorder[i]).imag() << endl;
	
	/*
	  (*testout) << "mata = " << bfa->GetMatrix() << endl
	  << " = " << endl << mata << endl;
	*/

	for (int i = 0; i < min2 (size_t(gfu->GetMultiDim()), reorder.Size()); i++)
	  {
	    FlatVector<Complex> vu = gfu->GetVector(i).FVComplex();
	    for (int j = 0; j < n; j++)
	      vu(j) = levecs(j, reorder[i]);
	  }



	// compute dlam/dalpha
	ofstream out_deriv("eigenvalues_deriv.dat");
	// BaseVector & hv = *gfu->GetVector(0).CreateVector();
          
#ifdef DERIVxxx
	Complex halpha = alpha;
	for (int i = 0; i < min2 (gfu->GetMultiDim(), reorder.Size()); i++)
	  {
	    double eps = 0.001;
	    BaseVector & vu = gfu->GetVector(i);
	    Complex lam = lami(reorder[i]);

	    hv = 0.0;
	    alpha = halpha + eps;
	    bfa -> AddMatrix (Complex (1/(2*eps)), vu, hv);
	    Complex aua =  S_InnerProduct<Complex> (vu, hv);
	    alpha = halpha - eps;
	    bfa -> AddMatrix (-Complex (1/(2*eps)), vu, hv);
	    Complex dlama = S_InnerProduct<Complex> (vu, hv);

	    hv = 0.0;
	    alpha = halpha + eps;
	    bfm -> AddMatrix (-lam * Complex (1/(2*eps)), vu, hv);
	    alpha = halpha - eps;
	    bfm -> AddMatrix (lam * Complex (1/(2*eps)), vu, hv);
	    Complex dlamm = S_InnerProduct<Complex> (vu, hv);

	    hv = bfm->GetMatrix() * vu;
	    Complex umu =  S_InnerProduct<Complex> (vu, hv);
	    
	    out_deriv << i << "  " 
		      << lam.real() << "  " << lam.imag() << "     " 
		      << " aua = " << aua 
		      << " da = " << dlama
		      << " dm = " << dlamm
		      << " diff = " << dlama+dlamm
		      << " dl/da = " << (dlama+dlamm)/umu << endl;
	  }
	alpha = halpha;
#endif

	
	int postit = 0;
	// postiterate:
	Mat<2,2,Complex> sma, smm;
	Vec<2,Complex> smlam;

	if (postit < 1)
	  // for (int i = 0; i < min2 (gfu->GetMultiDim(), reorder.Size()); i++)
          for (int i = 0; i < reorder.Size(); i++)
	    {
	      Complex omegai = sqrt (lami(reorder[i]));
	      out << i << "  " << omegai.real() << "  " << omegai.imag() << endl;
	    }

	for (int i = 0; i < min2 (size_t(gfu->GetMultiDim()), reorder.Size()); i++)
	  {
	    for (int k = 0; k < postit; k++)
	      {
		//	cout << "lam = " << lami(reorder[i]) << endl;
	       	BaseVector & vu = gfu->GetVector(i);
		//SZ hva = bfa->GetMatrix() * vu;
		hva = mata * vu;
	       	hvm = bfm->GetMatrix() * vu;
		hv = hva - lami(reorder[i]) * hvm;
		// cout << "def = " << hv.L2Norm() << endl;
		w = (*inva) * hv;
		
		//SZ	hwa = bfa->GetMatrix() * w;
		hwa = mata * w;
		hwm = bfm->GetMatrix() * w;

		sma(0,0) = S_InnerProduct<Complex> (hva, vu);
		sma(0,1) = sma(1,0) = S_InnerProduct<Complex> (hva, w);
		sma(0,0) = S_InnerProduct<Complex> (hwa, w);

		smm(0,0) = S_InnerProduct<Complex> (hvm, vu);
		smm(0,1) = sma(1,0) = S_InnerProduct<Complex> (hvm, w);
		smm(0,0) = S_InnerProduct<Complex> (hwm, w);

		LaEigNSSolve(2, &sma(0,0), &smm(0,0), &smlam(0), 0, 0, 0, 'N');
		//	cout << "smlami = " << smlam << endl;
		
		out << i << " " 
		    << lami(reorder[i]).real() << " " << lami(reorder[i]).imag() 
		    << "  ev1 = " << smlam(0) << " ev2 = " << smlam(1) << endl;
	      }
	  }
	
	
      }
    
  }


  
  static RegisterNumProc<NumProcEVP> npinitevp("evp");

}
#endif
