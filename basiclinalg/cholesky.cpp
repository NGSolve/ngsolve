#include <bla.hpp>

namespace ngbla
{
  using namespace ngbla;

/* 
   See Stoer, Einf. in die Num. Math, S 146
*/

  // Compute A = L D L^T decomposition
  // A_{ij} = \sum_k L_ik D_kk L_jk
  // L .. lower left factor, row-wise storage

  template <class T>
  void FlatCholeskyFactors<T> :: 
  Factor (const FlatMatrix<T> & a)
  {
    n = a.Height();
    bool dots = (n > 1000);

    //    diag = new T[n*(n+1)/2];
    lfact = diag+n;

    T x;
    
    for (int i = 0; i < n; i++)
      {
	if (dots && i % 10 == 0)
	  (cout) << "." << flush;
	
	for (int j = i; j < n; j++)
	  {
	    x = a(j,i);

	    T * pik = PRow (i);
	    T * pjk = PRow (j);
	    
	    for (int k = 0; k < i; k++)
	      x -= pjk[k] * diag[k] * Trans (pik[k]);

	    if (i == j)
	      {
		diag[i] = x;
	      }
	    else
	      {
		T invd;
		CalcInverse (diag[i], invd);
		// pjk[i] = invd * x;
		pjk[i] = x * invd;
	      }
	  }
      }

    for (int i = 0; i < n; i++)
      {
	T invd;
	CalcInverse (diag[i], invd);
	diag[i] = invd;
      }
    
    if (dots) (cout) << endl;

    /*
    (*testout) << "cholesky factors = " << endl;
    Print (*testout);
    cout << "experimental Cholesky!!!" << endl;
    diag[n-1] = 0;
    */
  }

  /*  
  template <class T>
  FlatCholeskyFactors<T> :: ~CholeskyFactors ()
  {
    delete [] diag;
  }
  */

  
  /*  gone to the header
  template <class T>  
  void FlatCholeskyFactors<T> :: 
  Mult (FlatVector<TV> x, FlatVector<TV> y) const
  {
    TV sum, val;
    const T *pj;

    for (int i = 0; i < n; i++)
      y(i) = x(i);
    
    for (int i = 0; i < n; i++)
      {
	sum = y(i);

	pj = PRow(i);
	for (int j = 0; j < i; ++j)
	  sum -= pj[j] * y(j);

	y(i) = sum;
      }

    for (int i = 0; i < n; i++)
      {
	sum = diag[i] * y(i);
	y(i) = sum;
      }

    for (int i = n-1; i >= 0; i--)
      {
	pj = PRow(i);
	val = y(i);
	for (int j = 0; j < i; ++j)
	  y(j) -= pj[j] * val;
      }
  }
  */
  
  template <class T>
  ostream & FlatCholeskyFactors<T> :: Print (ostream & ost) const
  {
    const T * pj;
    
    ost << "Diag: " << endl;
    for (int i = 0; i < n; i++)
      ost << i << ": " << diag[i] << endl;
    
    for (int i = 0; i < n; i ++)
      {
	ost << i << ": ";
	pj = PRow(i);
	for (int j = 0; j < i; ++j, ++pj)
	  ost << *pj << "  ";
	ost << endl;
      }
    return ost;
  }

  template class FlatCholeskyFactors<double>;
  template class FlatCholeskyFactors<Complex>;
#if MAX_SYS_DIM >= 1
  template class FlatCholeskyFactors<Mat<1,1,double> >;
  template class FlatCholeskyFactors<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class FlatCholeskyFactors<Mat<2,2,double> >;
  template class FlatCholeskyFactors<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class FlatCholeskyFactors<Mat<3,3,double> >;
  template class FlatCholeskyFactors<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class FlatCholeskyFactors<Mat<4,4,double> >;
  template class FlatCholeskyFactors<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class FlatCholeskyFactors<Mat<5,5,double> >;
  template class FlatCholeskyFactors<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class FlatCholeskyFactors<Mat<6,6,double> >;
  template class FlatCholeskyFactors<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class FlatCholeskyFactors<Mat<7,7,double> >;
  template class FlatCholeskyFactors<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class FlatCholeskyFactors<Mat<8,8,double> >;
  template class FlatCholeskyFactors<Mat<8,8,Complex> >;
#endif
}
