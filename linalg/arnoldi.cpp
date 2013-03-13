/**************************************************************************/
/* File:   arnoldi.cpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   August 2011                                                    */
/**************************************************************************/

/* 

Arnoldi Eigenvalue Solver
  
*/ 

#include <la.hpp>

namespace ngla
{
  
  template <typename SCAL>
  void Arnoldi<SCAL>::Calc (int numval, Array<Complex> & lam, int numev, Array<BaseVector*> & hevecs) const
  { 
    BaseVector & hv  = *a.CreateVector();
    BaseVector & hv2 = *a.CreateVector();
    BaseVector & hva = *a.CreateVector();
    BaseVector & hvm = *a.CreateVector();
   
    int n = hv.Size()*hv.EntrySize();    // matrix size

    /*
    SCAL x = 0.;
    Complex y = 0.;
    if ( typeid(x).name() ==  typeid(y).name()) n=n/2;
    */
    if (typeid(SCAL) == typeid(Complex) ) n=n/2;


    // int m = min2 (numval, n);
    int m = numval;


    // Matrix<SCAL> matV (m,n);
    Matrix<SCAL> matH(m);
    Array<BaseVector*> abv(m);
    for (int i = 0; i < m; i++)
      abv[i] = a.CreateVector();

    BaseMatrix & mat_shift = *a.CreateMatrix();
    mat_shift.AsVector() = a.AsVector() - shift*b.AsVector();  
    BaseMatrix * inv = mat_shift.InverseMatrix (freedofs);	  

    hv.SetRandom();
    hv.SetParallelStatus (CUMULATED);
    // hv = 1;
    FlatVector<double> fv = hv.FV<double>();
    if (freedofs)
      for (int i = 0; i < hv.Size(); i++)
	if (! (*freedofs)[i] ) fv(i) = 0;

    // matV = SCAL(0.0);   why ?
    matH = SCAL(0.0);

    hv2 = hv;
    SCAL len = sqrt (S_InnerProduct<SCAL> (hv, hv2)); // parallel
    hv /= len;

    for (int i = 0; i < m; i++)
      {
	cout << IM(1) << "i = " << i << endl;
	/*
	for (int j = 0; j < n; j++)
	  matV(i,j) = hv.FV<SCAL>()(j);
	*/
	*abv[i] = hv;

	hva = b * hv;
	hvm = *inv * hva;

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
	    matH(j,i) = S_InnerProduct<SCAL> (hvm, *abv[j]);
	    hvm -= matH(j,i) * *abv[j];
	  }
		
	hv = hvm;

	hv2 = hv;
	SCAL len = sqrt (S_InnerProduct<SCAL> (hv, hv2));
	if (i<m-1) matH(i+1,i) = len; 
	
	hv /= len;
      }
      
    delete inv;
	    
    cout << "has Hessenberg" << endl;
    *testout << "hessenberg = " << endl << matH << endl;
	    
    Vector<Complex> lami(m);
    Matrix<Complex> evecs(m);    
    Matrix<Complex> matHt(m);

    matHt = Trans (matH);
	    
    evecs = Complex (0.0);
    lami = Complex (0.0);
	
    LapackHessenbergEP (matH.Height(), &matHt(0,0), &lami(0), &evecs(0,0));
    
	    
    for (int i = 0; i < m; i++)
      lami(i) =  1.0 / lami(i) + shift;

    lam.SetSize (m);
    for (int i = 0; i < m; i++)
      lam[i] = lami(i);
    
    if (numev>0)
      {
	int nout = min2 (numev, m); 
	hevecs.SetSize(nout);
	for (int i = 0; i< nout; i++)
	  {
	    hevecs[i] = a.CreateVector();
	    *hevecs[i] = 0;
	    for (int j = 0; j < m; j++)
	      *hevecs[i] += evecs(i,j) * *abv[j];
	    // hevecs[i]->FVComplex() = Trans(matV)*evecs.Row(i);
	  }
      }
  } 
	

  template class Arnoldi<double>;
  template class Arnoldi<Complex>;


}
