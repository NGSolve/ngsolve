#include <bla.hpp>

namespace ngbla
{
  using namespace ngstd;
  using namespace ngbla;

  template <class TM>
  void CheckPos (const TM & m)
  {;}

  void CheckPos (const double & m)
  {
    if (m <= 0)
      {
	cout << "diag is " << m << endl;
	throw Exception ("diag is <= 0");
      }
  }



  // Compute A = L D L^T decomposition
  // A_{ij} = \sum_k L_ik D_kk L_jk
  // L .. lower left factor, columne-wise storage

  template <class T>
  void FlatBandCholeskyFactors<T> :: 
  Factor (const FlatSymBandMatrix<T> & a )
  {
    static int timer = NgProfiler::CreateTimer ("Band Cholesky");
    NgProfiler::RegionTimer reg (timer);

    T x;
    ArrayMem<T, 100> help(n);

    for (int i = 0; i < n; i++)
      {
        int mink = max2(0, i-bw+1);
        int ki = Index(i, mink);

        for (int k = mink; k < i; k++, ki++)
          help[k] = mem[k] * Trans (mem[ki]);

	int maxj = min2(n, i+bw);
	for (int j = i; j < maxj; j++)
	  {
	    x = a(j,i);

	    int mink = max2(0, j-bw+1);
	    int ki = Index(i, mink);
	    int kj = Index(j, mink);

            NgProfiler::AddFlops (timer, i-mink);

	    for (int k = mink; k < i; k++, ki++, kj++)
              x -= mem[kj] * help[k];
            // x -= mem[kj] * mem[k] * Trans (mem[ki]);

	    if (i == j)
	      {
		mem[i] = x;
	      }
	    else
	      {
		T invd;
		CalcInverse (mem[i], invd);
		(*this)(j,i) = x * invd;
	      }
	  }
      }

    for (int i = 0; i < n; i++)
      {
	T invd;
	CalcInverse (mem[i], invd);
	mem[i] = invd;
      }
  }
  
  
  /*
  template <class T>  
  void FlatBandCholeskyFactors<T> :: 
  Mult (const FlatVector<TV> & x, FlatVector<TV> & y) const
  {
    const TV * hx = &x(0);
    TV * hy = &y(0);
    const T * hm = &mem[0];

    for (int i = 0; i < n; i++)
      hy[i] = hx[i];

    int i, jj = n;
    for (i = 0; i < bw-1; i++)
      {
	TV sum = TSCAL(0.0);  

	for (int j = 0; j < i; j++, jj++)
	  sum += hm[jj] * hy[j];

	hy[i] -= sum;
      }

    for (  ; i < n; i++)
      {
	TV sum = TSCAL(0.0);

	for (int j = i-bw+1; j < i; j++, jj++)
	  sum += hm[jj] * hy[j];

	hy[i] -= sum;
      }

    for (int i = 0; i < n; i++)
      {
	TV sum = mem[i] * hy[i];
	hy[i] = sum;
      }


    // jj = n + (n-1) * (bw-1) - bw*(bw-1)/2;   
    for (i = n-1; i >= bw-1; i--)
      {
	jj -= bw-1;
	TV val = hy[i];

	int firstj = i-bw+1;
	for (int j = 0; j < bw-1; j++)
	  hy[firstj+j] -= Trans (mem[jj+j]) * val;
      }
    
    for (  ; i >= 0; i--) 
      {
	jj -= i;
	TV val = hy[i];

	for (int j = 0; j < i; j++)
	  hy[j] -= Trans (mem[jj+j]) * val;
      }
  }
  */

  
  template <class T>
  ostream & FlatBandCholeskyFactors<T> :: Print (ostream & ost) const
  {
    ost << "Diag: " << endl;
    for (int i = 0; i < n; i++)
      ost << i << ": " << mem[i] << endl;
    
    for (int i = 0; i < n; i ++)
      {
	ost << i << ": ";
	for (int j = max2(0, i-bw+1); j < i; j++)
	  ost << (*this)(i,j) << "  ";
	ost << endl;
      }
    return ost;
  }

  template class FlatBandCholeskyFactors<double>;
  template class FlatBandCholeskyFactors<Complex>;
#if MAX_SYS_DIM >= 1
  template class FlatBandCholeskyFactors<Mat<1,1,double> >;
  template class FlatBandCholeskyFactors<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class FlatBandCholeskyFactors<Mat<2,2,double> >;
  template class FlatBandCholeskyFactors<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class FlatBandCholeskyFactors<Mat<3,3,double> >;
  template class FlatBandCholeskyFactors<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class FlatBandCholeskyFactors<Mat<4,4,double> >;
  template class FlatBandCholeskyFactors<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class FlatBandCholeskyFactors<Mat<5,5,double> >;
  template class FlatBandCholeskyFactors<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class FlatBandCholeskyFactors<Mat<6,6,double> >;
  template class FlatBandCholeskyFactors<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class FlatBandCholeskyFactors<Mat<7,7,double> >;
  template class FlatBandCholeskyFactors<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class FlatBandCholeskyFactors<Mat<8,8,double> >;
  template class FlatBandCholeskyFactors<Mat<8,8,Complex> >;
#endif
}
