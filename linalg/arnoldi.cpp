/**************************************************************************/
/* File:   arnoldi.cpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   August 2011                                                    */
/**************************************************************************/

/* 

Arnoldi Eigenvalue Solver
  
*/ 

// #include <la.hpp>
#include "vvector.hpp"
#include "arnoldi.hpp"
#include "cg.hpp"

namespace ngla
{
  
  template <typename SCAL>
  void Arnoldi<SCAL>::Calc (int inumval, Array<Complex> & lam, int inumev, 
                            Array<shared_ptr<BaseVector>> & hevecs, 
                            shared_ptr<BaseMatrix> pre) const
  {
    size_t numval = inumval;
    size_t numev = inumev;
    
    static Timer t("arnoldi");    
    static Timer t2("arnoldi - orthogonalize");    
    static Timer t3("arnoldi - compute large vectors");

    RegionTimer reg(t);

    auto hv  = a->CreateColVector();
    auto hv2 = a->CreateColVector();
    auto hva = a->CreateColVector();
    auto hvm = a->CreateColVector();
   
    size_t n = hv.template FV<SCAL>().Size();
    size_t m = min2 (numval, n);


    Matrix<SCAL> matH(m);
    Array<shared_ptr<BaseVector>> abv(m);
    for (int i = 0; i < m; i++)
      abv[i] = a->CreateColVector();

    auto mat_shift = a->CreateMatrix();
    mat_shift->AsVector() = a->AsVector() - shift*b->AsVector();  
    shared_ptr<BaseMatrix> inv;
    if (!pre)
      {
        if (inversetype)
          mat_shift -> SetInverseType(*inversetype);
        inv = mat_shift->InverseMatrix (freedofs);
      }
    else
      {
        auto itso = make_shared<GMRESSolver<double>> (mat_shift, pre);
        itso->SetPrintRates(1);
        itso->SetMaxSteps(2000);
        inv = itso;
      }

    hv.SetRandom();
    hv.SetParallelStatus (CUMULATED);
    FlatVector<SCAL> fv = hv.template FV<SCAL>();
    if (freedofs)
      for (int i = 0; i < hv.Size(); i++)
	if (! (*freedofs)[i] ) fv(i) = 0;

    t2.Start();
    // matV = SCAL(0.0);   why ?
    matH = SCAL(0.0);

    *hv2 = *hv;
    SCAL len = sqrt (S_InnerProduct<SCAL> (*hv, *hv2)); // parallel
    *hv /= len;
    
    for (int i = 0; i < m; i++)
      {
	cout << IM(1) << "\ri = " << i << "/" << m << flush;
	/*
	for (int j = 0; j < n; j++)
	  matV(i,j) = hv.FV<SCAL>()(j);
	*/
	*abv[i] = *hv;

	*hva = *b * *hv;
	*hvm = *inv * *hva;

	for (int j = 0; j <= i; j++)
	  {
            /*
            SCAL sum = 0.0;
	    for (int k = 0; k < n; k++)
	      sum += hvm.FV<SCAL>()(k) * matV(j,k);
	    matH(j,i) = sum;
	    for (int k = 0; k < n; k++)
	      hvm.FV<SCAL>()(k) -= sum * matV(j,k);
            */
            /*
            SCAL sum = 0.0;
            FlatVector<SCAL> abvj = abv[j] -> FV<SCAL>();
            FlatVector<SCAL> fv_hvm = hvm.FV<SCAL>();
	    for (int k = 0; k < n; k++)
	      sum += fv_hvm(k) * abvj(k);
	    matH(j,i) = sum;
	    for (int k = 0; k < n; k++)
	      fv_hvm(k) -= sum * abvj(k);
            */

	    matH(j,i) = S_InnerProduct<SCAL> (*hvm, *abv[j]);
	    *hvm -= matH(j,i) * *abv[j];
	  }
		
	*hv = *hvm;
	*hv2 = *hv;
	SCAL len = sqrt (S_InnerProduct<SCAL> (*hv, *hv2));
	if (i<m-1) matH(i+1,i) = len; 
	
	*hv /= len;
      }
      
    t2.Stop();
    t2.AddFlops (double(n)*m*m);
    cout << IM(3) << "n = " << n << ", m = " << m << " n*m*m = " << n*m*m << endl;
    cout << IM(3) << "\ri = " << m << "/" << m << endl;

	    
    Vector<Complex> lami(m);
    Matrix<Complex> evecs(m);    
    Matrix<Complex> matHt(m);

    matHt = Trans (matH);
    
    evecs = Complex (0.0);
    lami = Complex (0.0);

    cout << IM(5) << "Solve Hessenberg evp with Lapack ... " << flush;
    LapackHessenbergEP (matH.Height(), &matHt(0,0), &lami(0), &evecs(0,0));
    cout << IM(5) << "done" << endl;
	    
    for (int i = 0; i < m; i++)
      lami(i) =  1.0 / lami(i) + shift;

    lam.SetSize (m);
    for (int i = 0; i < m; i++)
      lam[i] = lami(i);

    t3.Start();
    if (numev>0)
      {
	int nout = min2 (numev, m); 
	hevecs.SetSize(nout);
	for (int i = 0; i< nout; i++)
	  {
            if (a->IsComplex())
              hevecs[i] = a->CreateColVector();
            else // real biform and system-vecors not yet supported
              hevecs[i] =  make_shared<VVector<Complex>> (a->Height());
            
	    *hevecs[i] = 0;
	    for (int j = 0; j < m; j++)
	      *hevecs[i] += evecs(i,j) * *abv[j];
	    // hevecs[i]->FVComplex() = Trans(matV)*evecs.Row(i);
	  }
      }
    t3.Stop();
  } 
	

  template class Arnoldi<double>;
  template class Arnoldi<Complex>;


}
