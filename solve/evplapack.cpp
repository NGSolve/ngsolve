#ifdef LAPACK
#include "../include/solve.hpp"

namespace ngfem
{
  extern Complex alpha;
}

namespace ngsolve
{
  using namespace ngsolve;


  /* ***************************** Numproc EVP-Lapack *************************** */

  ///
  class NumProcEVPLapack : public NumProc
  {
  protected:
    BilinearForm * bfa;
    BilinearForm * bfm;

    GridFunction * gfu;
    Preconditioner * pre;
    int num;

    double prec, shift;
    bool print;

    string eigenoutstr;
    enum SOLVER { DENSE, ARNOLDI };
    ///
    SOLVER solver;

  public:
    NumProcEVPLapack (PDE & apde, const Flags & flags);
    virtual ~NumProcEVPLapack();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcEVPLapack (pde, flags);
    }

    virtual void Do(LocalHeap & lh);

    virtual string GetClassName () const
    {
      return " Eigenvalue Problem (Lapack)";
    }

    virtual void PrintReport (ostream & ost)
    {
      ost << GetClassName() << endl
	  << "Bilinear-form A = " << bfa->GetName() << endl
	  << "Bilinear-form M = " << bfm->GetName() << endl
	  << "Gridfunction  = " << gfu->GetName() << endl;
    }
  };


  NumProcEVPLapack :: NumProcEVPLapack (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearforma", ""));
    bfm = pde.GetBilinearForm (flags.GetStringFlag ("bilinearformm", ""));
    gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
    pre = pde.GetPreconditioner (flags.GetStringFlag ("preconditioner", ""),1);
    num = int(flags.GetNumFlag ("num", 500));
    shift = flags.GetNumFlag ("shift",1.0); 

    eigenoutstr = flags.GetStringFlag ("eigenout","eigen.out"); 

    solver = ARNOLDI;
    if (flags.GetDefineFlag("dense")) solver = DENSE;
  }

  NumProcEVPLapack :: ~NumProcEVPLapack()
  {
    ;
  }


  void NumProcEVPLapack :: Do(LocalHeap & lh)
  {
    cout << "solve evp with Lapack" << endl;

    if (solver == ARNOLDI)
      {

        cout << "new version using linalg - Arnoldi" << endl;
        
        
        if (bfa->GetFESpace().IsComplex())
          {
            cout << "complex evp" << endl;

            Arnoldi<Complex> arnoldi (bfa->GetMatrix(), bfm->GetMatrix(), bfa->GetFESpace().GetFreeDofs() );
            arnoldi.SetShift (Complex(shift,0));
            
            int nev = gfu->GetMultiDim();
            Array<BaseVector*> evecs(nev);
            Array<Complex> lam(nev);
            if (pre)
              arnoldi.Calc (num, lam, nev, evecs, &pre->GetMatrix());
            else
              arnoldi.Calc (num, lam, nev, evecs, 0);
            
            ofstream eigenout(eigenoutstr.c_str());
            for (int i = 0; i < lam.Size(); ++i)
              eigenout << lam[i].real() << "\t" << lam[i].imag() << "\t"
                       << sqrt(lam[i]).real() << "\t" << sqrt(lam[i]).imag() 
                       << endl;
            
            for (int i = 0; i < nev; i++)
              gfu->GetVector(i) = *evecs[i];
          }
        else
          {
            cout << "real evp" << endl;
            Arnoldi<double> arnoldi (bfa->GetMatrix(), bfm->GetMatrix(), bfa->GetFESpace().GetFreeDofs() );
            arnoldi.SetShift (shift);
            
            int nev = gfu->GetMultiDim();
            Array<BaseVector*> evecs(nev);
            Array<Complex> lam(nev);
            if (pre)
              arnoldi.Calc (num, lam, nev, evecs, &pre->GetMatrix());
            else
              arnoldi.Calc (num, lam, nev, evecs, 0);
            
            ofstream eigenout(eigenoutstr.c_str());
            for (int i = 0; i < lam.Size(); ++i)
              eigenout << lam[i].real() << "\t" << lam[i].imag() << endl;
            
            for (int i = 0; i < nev; i++)
              gfu->GetVector(i) = *evecs[i];
            // cout << "lam = " << endl << lam << endl;
          }
        
        return;
      }

    int dim = bfa->GetFESpace().GetDimension();
    int size = bfa->GetMatrix().Height();
    bool iscomplex = bfa->GetFESpace().IsComplex();

    
    if (solver == DENSE)
      {
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
	    
	    LaEigNSSolve (n, &mata(0,0), &matm(0,0), &lami(0), 1, &evecs(0,0), 0, 'N');

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
		Complex lam =  (lami(i));
		ev << i << " " << lam.real() << " " << lam.imag() << endl;
	      }
	    cout << "eigensystem done" << endl;
	  }
      }
  }
  




  


  namespace
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    { 
      GetNumProcs().AddNumProc ("evplapack", NumProcEVPLapack::Create, NumProcEVPLapack::PrintDoc);
    }
    
    
    Init init;
    
  }
  


}
#endif
