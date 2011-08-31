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
    BaseVector & hv = *a.CreateVector();
    BaseVector & hva = *a.CreateVector();
    BaseVector & hvm = *a.CreateVector();
    // BaseVector & hwa = *a.CreateVector();
    // BaseVector & hwm = *a.CreateVector();
    // BaseVector & w = *a.CreateVector();
    
    int n = hv.Size()*hv.EntrySize();    // matrix size
    int m = min2 (numval, n);

    Matrix<SCAL> matV (m,n), matVt(n,m);
    Matrix<SCAL> matH(m), matI(m);



    const BaseSparseMatrix&  mata1 =
      dynamic_cast<const BaseSparseMatrix&> (a);
    
    BaseSparseMatrix&  mata =
      dynamic_cast<BaseSparseMatrix&> (*mata1.CreateMatrix());
    
    const BaseSparseMatrix& matm = 
      dynamic_cast<const BaseSparseMatrix&> (b);
    
    mata.AsVector() = mata1.AsVector() - shift*matm.AsVector();  
    
    BaseMatrix * inva = mata.InverseMatrix();	  
    

    hv.SetRandom();
    matV = SCAL(0.0);
    matH = SCAL(0.0);
    matI = SCAL(0.0);
    for (int i = 0; i < m; i++)
      matI(i,i) = 1.0;

    SCAL len = sqrt (S_InnerProduct<SCAL> (hv, hv));
    hv /= len;



    for (int i = 0; i < m; i++)
      {
	cout << "i = " << i << endl;
	for (int j = 0; j < n; j++)
	  matV(i,j) = hv.FV<SCAL>()(j);
		
	hva = matm * hv;
	hvm = *inva * hva;

	for (int j = 0; j <= i; j++)
	  {
	    SCAL sum = 0.0;
	    for (int k = 0; k < n; k++)
	      sum += hvm.FV<SCAL>()(k) * matV(j,k);
	    matH(j,i) = sum;
	    for (int k = 0; k < n; k++)
	      hvm.FV<SCAL>()(k) -= sum * matV(j,k);
	  }
		
	hv = hvm;
	
	SCAL len = sqrt (S_InnerProduct<SCAL> (hv, hv));
	if (i<m-1) matH(i+1,i) = len; 
	
	hv /= len;
      }
	    
    cout << "has Hessenberg" << endl;
	    
    Vector<Complex> lami(m);
    Matrix<Complex> evecs(m);
    Matrix<Complex> levecs(n,m);
    Matrix<Complex> matHt(m);

    matHt = Trans (matH);
	    
    evecs = Complex (0.0);
    lami = Complex (0.0);
	    
    LapackHessenbergEP (matH.Height(), &matHt(0,0), &lami(0), &evecs(0,0));
	    
    // levecs = Trans(matV) * Trans(evecs);
    matVt = Trans (matV);
    levecs = Trans(evecs * Trans (matVt));
	    
    for (int i = 0; i < m; i++)
      lami(i) =  1.0 / lami(i) + shift;

    lam.SetSize (m);
    for (int i = 0; i < m; i++)
      lam[i] = lami(i);



    ofstream out ("eigenvalues.dat");
    out.precision (12);
    for (int i = 0; i < m; i++) 
      out << i << "  " << lami[i].real() << "  " << lami[i].imag() << endl;
    cout << "wrote eigenvalues to file 'eigenvalues.dat'" << endl;

    cout << "eigensystem done" << endl;

    /*
    for (int i = 0; i < min2 (gfu->GetMultiDim(), m); i++)
      {
	FlatVector<Complex> vu = gfu->GetVector(i).FVComplex();
	for (int j = 0; j < n; j++)
	  vu(j) = levecs(j, i);
      }
    */
  } 
	

  template class Arnoldi<double>;
  template class Arnoldi<Complex>;




}
