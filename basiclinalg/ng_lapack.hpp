#ifndef FILE_NG_LAPACK
#define FILE_NG_LAPACK

// #include <mkl_cblas.h>

/****************************************************************************/
/* File:   ng_lapack.hpp                                                    */
/* Author: Joachim Schoeberl                                                */
/* Date:   21. Nov. 2004                                                    */
/****************************************************************************/



namespace ngbla 
{

  class T_Lapack { };
  static constexpr T_Lapack Lapack;

  template <typename TA>
  class LapackExpr : public Expr<LapackExpr<TA> >
  {
    const TA & a;
  public:
    LapackExpr (const TA & aa) : a(aa) { ; }
    const TA & A() const { return a; }
    size_t Height() const { return a.Height(); }
    size_t Width() const { return a.Width(); }
  };
  
  template <typename TA>
  INLINE LapackExpr<TA> operator| (const Expr<TA> & a, T_Lapack /* tl */)
  {
    return LapackExpr<TA> (a.Spec());
  }

  
  template <typename TOP, typename T, typename TB>
  class assign_trait<TOP, T, LapackExpr<TB>, int>
  {
  public:
    static INLINE T & Assign (MatExpr<T> & self, const Expr<LapackExpr<TB>> & v)
    {
    #ifdef LAPACK
      if constexpr (std::is_same_v<TOP,typename MatExpr<T>::As>)
                     LapackMultAdd (v.Spec().A().A(), v.Spec().A().B(), 1.0, self.Spec(), 0.0);
      if constexpr (std::is_same_v<TOP,typename MatExpr<T>::AsAdd>)
                     LapackMultAdd (v.Spec().A().A(), v.Spec().A().B(), 1.0, self.Spec(), 1.0);
      if constexpr (std::is_same_v<TOP,typename MatExpr<T>::AsSub>)
                     LapackMultAdd (v.Spec().A().A(), v.Spec().A().B(), -1.0, self.Spec(), 1.0);
      return self.Spec();
    #else // LAPACK
      throw Exception("No Lapack");
    #endif // LAPACK
    }
  };  


  



  
#ifdef LAPACK

  extern "C" {
#ifdef MKL_ILP64
    typedef long int integer;
#else
    typedef int integer;
#endif
    // typedef char logical;
    typedef integer logical;
    typedef float real;
    typedef double doublereal;
    typedef Complex doublecomplex;
    typedef complex<float> singlecomplex;


// Windows SDK defines VOID in the file WinNT.h
#ifndef VOID
    typedef void VOID;
#endif

    typedef int ftnlen;
    typedef int L_fp;  // ?


#include "clapack.h"
  }


  // Interface to lapack functions
  NGS_DLL_HEADER int sgemm(char *transa, char *transb, integer *m, integer *
		  n, integer *k, real *alpha, real *a, integer *lda, 
		  real *b, integer *ldb, real *beta, real *c__, 
		  integer *ldc);

  NGS_DLL_HEADER int dgemm(char *transa, char *transb, integer *m, integer *
		  n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
		  doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
		  integer *ldc);

  NGS_DLL_HEADER int sgemm(char *transa, char *transb, integer *m, integer *
		  n, integer *k, real *alpha, real *a, integer *lda, 
		  real *b, integer *ldb, real *beta, real *c__, 
		  integer *ldc);

  NGS_DLL_HEADER int zgemm(char *transa, char *transb, integer *m, integer *
		    n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, 
		    doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *
		    c__, integer *ldc);

  NGS_DLL_HEADER int dger(integer *m, integer *n, doublereal *alpha,
                   doublereal *x, integer *incx, doublereal *y, integer *incy,
                   doublereal *a, integer *lda);
  
  NGS_DLL_HEADER int dgetrf(integer* n, integer* m, double* a, integer* lda, integer* ipiv, integer* info);
  NGS_DLL_HEADER int dgetri(integer* n, double* a, integer* lda, integer* ipiv,
                            double* hwork, integer* lwork, integer* info);
  NGS_DLL_HEADER int dgetrs(char *trans, integer *n, integer *nrhs, 
                            doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
                            ldb, integer *info);
  

  inline int gemm(char *transa, char *transb, integer *m, integer *
      n, integer *k, real *alpha, real *a, integer *lda,
      real *b, integer *ldb, real *beta, real *c__,
      integer *ldc)
  {
    return sgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
  }

  inline int gemm(char *transa, char *transb, integer *m, integer *
      n, integer *k, doublereal *alpha, doublereal *a, integer *lda,
      doublereal *b, integer *ldb, doublereal *beta, doublereal *c__,
      integer *ldc)
  {
    return dgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
  }

  inline int gemm(char *transa, char *transb, integer *m, integer *
      n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda,
      doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *
      c__, integer *ldc)
  {
    return zgemm (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
  }

  
  // BLAS 1

  inline double LapackDot (FlatVector<double> x, FlatVector<double> y)
  {
    integer n = x.Size();
    integer incx = 1;
    integer incy = 1;
    return ddot_ (&n, &x(0), &incx, &y(0), &incy);
  }



  // BLAS 2

  /*
  extern "C" 
  void dgemv_ (char & trans, int & m, int & n, double & alpha, double & a, int & lda, double & x, int & incx,
               double & beta, double & y, int & incy);
  */

  inline void LapackMultAx (ngbla::FlatMatrix<double> a,
                            ngbla::FlatVector<double> x,
                            ngbla::FlatVector<double> y)
  {
    char trans = 'T';
    integer m = a.Width();
    integer n = a.Height();
    double alpha = 1;
    integer lda = max(size_t(1), a.Width());
    integer incx = 1;
    double beta = 0;
    integer incy = 1;
    dgemv_ (&trans, &m, &n, &alpha, &a(0,0), &lda, &x(0), &incx, &beta, &y(0), &incy);
  }

  inline void LapackMultAx (ngbla::FlatMatrix<Complex> a,
                            ngbla::FlatVector<Complex> x,
                            ngbla::FlatVector<Complex> y)
  {
    char trans = 'T';
    integer m = a.Width();
    integer n = a.Height();
    Complex alpha(1,0);
    integer lda = max(size_t(1), a.Width());
    integer incx = 1;
    Complex beta(0, 0);
    integer incy = 1;
    zgemv_ (&trans, &m, &n, &alpha, &a(0,0), &lda, &x(0), &incx, &beta, &y(0), &incy);
  }



  inline void LapackMultAtx (ngbla::FlatMatrix<double> a,
                             ngbla::FlatVector<double> x,
                             ngbla::FlatVector<double> y)
  {
    char trans = 'N';
    integer m = a.Width();
    integer n = a.Height();
    double alpha = 1;
    integer lda = max(size_t(1), a.Width());
    integer incx = 1;
    double beta = 0;
    integer incy = 1;
    dgemv_ (&trans, &m, &n, &alpha, &a(0,0), &lda, &x(0), &incx, &beta, &y(0), &incy);
  }


  inline void LapackAddxyt (ngbla::FlatMatrix<double> a,
                            double fac,
                            ngbla::FlatVector<double> x,
                            ngbla::FlatVector<double> y)
  {
    integer m = a.Width();
    integer n = a.Height();
    double alpha = fac;
    integer lda = max(size_t(1), a.Width());
    integer incx = 1;
    integer incy = 1;

    dger (&m, &n, &alpha, &y(0), &incx, &x(0), &incy, &a(0,0), &lda);
  }





  template <typename SCAL>
  inline void BASE_LapackMult (SliceMatrix<SCAL> a, bool transa, 
			       SliceMatrix<SCAL> b, bool transb, 
			       SliceMatrix<SCAL> c)
  {
    char transa_ = transa ? 'T' : 'N';
    char transb_ = transb ? 'T' : 'N'; 

    integer n = c.Width();
    integer m = c.Height();
    if (n == 0 || m == 0) return;
    integer k = transa ? a.Height() : a.Width();
    SCAL alpha = 1.0;
    SCAL beta = 0;
    integer lda = max(size_t(1), a.Dist());
    integer ldb = max(size_t(1), b.Dist());
    integer ldc = max(size_t(1), c.Dist());

    gemm (&transb_, &transa_, &n, &m, &k, &alpha, 
	  &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }

  template <typename SCAL>
  inline void BASE_LapackMultAdd (SliceMatrix<SCAL> a, bool transa, 
				  SliceMatrix<SCAL> b, bool transb, 
				  SCAL aalpha,
				  SliceMatrix<SCAL> c,
				  SCAL abeta)
  {
    char transa_ = transa ? 'T' : 'N';
    char transb_ = transb ? 'T' : 'N'; 

    integer n = c.Width();
    integer m = c.Height();
    if (n == 0 || m == 0) return;
    integer k = transa ? a.Height() : a.Width();
    SCAL alpha = aalpha;
    SCAL beta = abeta;
    integer lda = max(size_t(1), a.Dist());
    integer ldb = max(size_t(1), b.Dist());
    integer ldc = max(size_t(1), c.Dist());

    gemm (&transb_, &transa_, &n, &m, &k, &alpha, 
	  &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }



  template <typename TA, typename TB, typename TC>
  inline void LapackMult (const TA & a, 
			  const TB & b, 
			  const TC & c)
  { BASE_LapackMult<typename TC::TSCAL> (a, false, b, false, c); }

  template <typename TA, typename TB, typename TC>
  inline void LapackMult (const TA & a, 
			  TransExpr<TB> b, 
			  const TC & c)
  { BASE_LapackMult<typename TC::TSCAL> (a, false, b.A(), true, c); }

  template <typename TA, typename TB, typename TC>
  inline void LapackMult (TransExpr<TA> a, 
			  const TB & b, 
			  const TC & c)
  { BASE_LapackMult<typename TC::TSCAL> (a.A(), true, b, false, c); }

  template <typename TA, typename TB, typename TC>
  inline void LapackMult (TransExpr<TA> a, 
			  TransExpr<TB> b, 
			  const TC & c)
  { BASE_LapackMult<typename TC::TSCAL> (a.A(), true, b.A(), true, c); }

  

  /*
  inline void LapackMult (SliceMatrix<double> a, 
			  SliceMatrix<double> b, 
			  SliceMatrix<double> c)
  { BASE_LapackMult<double> (a, false, b, false, c); }

  template <typename TA>
  inline void LapackMult (SliceMatrix<double> a, 
			  TransExpr<TA> b, 
			  SliceMatrix<double> c)
  { BASE_LapackMult<double> (a, false, b.A(), true, c); }

  template <typename TA>
  inline void LapackMult (TransExpr<TA> a, 
			  SliceMatrix<double> b, 
			  SliceMatrix<double> c)
  { BASE_LapackMult<double> (a.A(), true, b, false, c); }

  template <typename TA, typename TB>
  inline void LapackMult (TransExpr<TA> a, 
			  TransExpr<TB> b,
			  SliceMatrix<double> c)
  { BASE_LapackMult<double> (a.A(), true, b.A(), true, c); }
  */


  /*
  inline void LapackMult (SliceMatrix<Complex> a, 
			  SliceMatrix<Complex> b, 
			  SliceMatrix<Complex> c)
  { BASE_LapackMult<Complex> (a, false, b, false, c); }

  template <typename TA>
  inline void LapackMult (SliceMatrix<Complex> a, 
			  TransExpr<TA> b, 
			  SliceMatrix<Complex> c)
  { BASE_LapackMult<Complex> (a, false, b.A(), true, c); }

  template <typename TA>
  inline void LapackMult (TransExpr<TA> a, 
			  SliceMatrix<Complex> b, 
			  SliceMatrix<Complex> c)
  { BASE_LapackMult<Complex> (a.A(), true, b, false, c); }
  
  template <typename TA, typename TB>
  inline void LapackMult (TransExpr<TA> a, 
			  TransExpr<TB> b,
			  SliceMatrix<Complex> c)
  { BASE_LapackMult<Complex> (a.A(), true, b.A(), true, c); }
  */



  inline void LapackMultAdd (SliceMatrix<float> a, 
			     SliceMatrix<float> b, 
			     float alpha,
			     SliceMatrix<float> c,
			     float beta)
  { BASE_LapackMultAdd<float> (a, false, b, false, alpha, c, beta); }

  inline void LapackMultAdd (SliceMatrix<double> a, 
			     SliceMatrix<double> b, 
			     double alpha,
			     SliceMatrix<double> c,
			     double beta)
  { BASE_LapackMultAdd<double> (a, false, b, false, alpha, c, beta); }

  inline void LapackMultAdd (SliceMatrix<double,ColMajor> a, 
			     SliceMatrix<double> b, 
			     double alpha,
			     SliceMatrix<double> c,
			     double beta)
  { BASE_LapackMultAdd<double> (Trans(a), true, b, false, alpha, c, beta); }

  inline void LapackMultAdd (SliceMatrix<double> a, 
			     SliceMatrix<double,ColMajor> b, 
			     double alpha,
			     SliceMatrix<double> c,
			     double beta)
  { BASE_LapackMultAdd<double> (a, false, Trans(b), true, alpha, c, beta); }

  inline void LapackMultAdd (SliceMatrix<double,ColMajor> a, 
			     SliceMatrix<double,ColMajor> b, 
			     double alpha,
			     SliceMatrix<double> c,
			     double beta)
  { BASE_LapackMultAdd<double> (Trans(a), true, Trans(b), true, alpha, c, beta); }



  /*
  template <typename TA>
  inline void LapackMultAdd (SliceMatrix<double> a, 
			     TransExpr<TA> b, 
			     double alpha,			     
			     SliceMatrix<double> c,
			     double beta = 1.0)
  { BASE_LapackMultAdd<double> (a, false, b.A(), true, alpha, c, beta); }

  template <typename TA>
  inline void LapackMultAdd (TransExpr<TA> a, 
			     SliceMatrix<double> b, 
			     double alpha,
			     SliceMatrix<double> c,
			     double beta = 1.0)
  { BASE_LapackMultAdd<double> (a.A(), true, b, false, alpha, c, beta); }
  
  template <typename TA, typename TB>
  inline void LapackMultAdd (TransExpr<TA> a, 
			     TransExpr<TB> b,
			     double alpha,
			     SliceMatrix<double> c,
			     double beta = 1.0)
  { BASE_LapackMultAdd<double> (a.A(), true, b.A(), true, alpha, c, beta); }
  */


  inline void LapackMultAdd (SliceMatrix<Complex> a, 
			     SliceMatrix<Complex> b, 
			     Complex alpha,
			     SliceMatrix<Complex> c,
			     Complex beta = 1.0)
  { BASE_LapackMultAdd<Complex> (a, false, b, false, alpha, c, beta); }

  inline void LapackMultAdd (SliceMatrix<Complex,ColMajor> a, 
			     SliceMatrix<Complex> b, 
			     Complex alpha,
			     SliceMatrix<Complex> c,
			     Complex beta)
  { BASE_LapackMultAdd<Complex> (Trans(a), true, b, false, alpha, c, beta); }

  inline void LapackMultAdd (SliceMatrix<Complex> a, 
			     SliceMatrix<Complex,ColMajor> b, 
			     Complex alpha,
			     SliceMatrix<Complex> c,
			     Complex beta)
  { BASE_LapackMultAdd<Complex> (a, false, Trans(b), true, alpha, c, beta); }

  inline void LapackMultAdd (SliceMatrix<Complex,ColMajor> a, 
			     SliceMatrix<Complex,ColMajor> b, 
			     Complex alpha,
			     SliceMatrix<Complex> c,
			     Complex beta)
  { BASE_LapackMultAdd<Complex> (Trans(a), true, Trans(b), true, alpha, c, beta); }


  template <typename TA>
  inline void LapackMultAdd (SliceMatrix<Complex> a, 
			     TransExpr<TA> b, 
			     Complex alpha,
			     SliceMatrix<Complex> c,
			     Complex beta = 1.0)
  { BASE_LapackMultAdd<Complex> (a, false, b.A(), true, alpha, c, beta); }

  template <typename TA>
  inline void LapackMultAdd (TransExpr<TA> a, 
			     SliceMatrix<Complex> b, 
			     Complex alpha,
			     SliceMatrix<Complex> c,
			     Complex beta)
  { BASE_LapackMultAdd<Complex> (a.A(), true, b, false, alpha, c, beta); }
  
  template <typename TA, typename TB>
  inline void LapackMultAdd (TransExpr<TA> a, 
			     TransExpr<TB> b,
			     Complex alpha,
			     SliceMatrix<Complex> c,
			     Complex beta)
  { BASE_LapackMultAdd<Complex> (a.A(), true, b.A(), true, alpha, c, beta); }



  template <typename TA, typename TB, typename Talpha, typename Tc, typename Tbeta>
  inline void LapackMultAdd (MinusExpr<TA> a, 
			     const TB & b, 
			     Talpha alpha,			     
			     const Tc & c,
			     Tbeta beta)
  { LapackMultAdd (a.A(), b, -alpha, c, beta); }



  


  // we don't have a lapack function for that 
  void LapackMultAdd (SliceMatrix<double> a, 
                      SliceMatrix<Complex,ColMajor> b, 
                      Complex alpha,
                      SliceMatrix<Complex> c,
                      Complex beta);
  // ---> moved implementation to ng_blas
  /*
  {
    if (beta == 0.0)
      c = alpha * (a * b);
    else
      {
        c *= beta;
        c += alpha * (a * b);
      }
    // BASE_LapackMultAdd<double> (Trans(a), true, Trans(b), true, alpha, c, beta);
  }
  */

  template <typename TA, typename TB, typename ALPHA, typename BETA, typename TC>
  inline void LapackMultAdd (TA a, 
                             TB b, 
                             ALPHA alpha,
                             SliceMatrix<TC,ColMajor> c,
                             BETA beta)
  {
    LapackMultAdd (Trans(b), Trans(a), alpha, Trans(c), beta);
  }
  

  /*

    // old-style function names, new driver
  inline void LapackMultABt (SliceMatrix<double> a, 
			     SliceMatrix<double> b, 
			     SliceMatrix<double> c)
  { BASE_LapackMult (a, false, b, true, c); }

  inline void LapackMultABt (SliceMatrix<Complex> a, 
			     SliceMatrix<Complex> b, 
			     SliceMatrix<Complex> c)
  { BASE_LapackMult (a, false, b, true, c); }

  inline void LapackMultAB (SliceMatrix<double> a, 
			    SliceMatrix<double> b, 
			    SliceMatrix<double> c)
  { BASE_LapackMult (a, false, b, false, c); }

  inline void LapackMultAB (SliceMatrix<Complex> a, 
			    SliceMatrix<Complex> b, 
			    SliceMatrix<Complex> c)
  { BASE_LapackMult (a, false, b, false, c); }


  inline void LapackMultAtB (SliceMatrix<double> a, 
			     SliceMatrix<double> b, 
			     SliceMatrix<double> c)
  { BASE_LapackMult (a, true, b, false, c); }

  inline void LapackMultAtB (SliceMatrix<Complex> a, 
			     SliceMatrix<Complex> b, 
			     SliceMatrix<Complex> c)
  { BASE_LapackMult (a, true, b, false, c); }


  inline void LapackMultAtBt (SliceMatrix<double> a, 
			      SliceMatrix<double> b, 
			      SliceMatrix<double> c)
  { BASE_LapackMult (a, true, b, true, c); }
  inline void LapackMultAtBt (SliceMatrix<Complex> a, 
			      SliceMatrix<Complex> b, 
			      SliceMatrix<Complex> c)
  { BASE_LapackMult (a, true, b, true, c); }

  */






    // old functions for compatibility 
  inline void LapackMultABt (SliceMatrix<double> a, 
			     SliceMatrix<double> b, 
			     SliceMatrix<double> c)
  {
    char transa = 'T';
    char transb = 'N';
    integer m = c.Height();
    integer n = c.Width();
    integer k = a.Width();
    double alpha = 1.0;
    double beta = 0;
    integer lda = max(size_t(1), a.Dist());
    integer ldb = max(size_t(1), b.Dist());
    integer ldc = max(size_t(1), c.Dist());

    dgemm (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }

 
  inline void LapackMultAB (SliceMatrix<double> a, 
			    SliceMatrix<double> b, 
			    SliceMatrix<double> c)
  {
    char transa = 'N';
    char transb = 'N';
    integer m = c.Height();
    integer n = c.Width();
    integer k = a.Width();
    double alpha = 1.0;
    double beta = 0;
    integer lda = max(size_t(1), a.Dist());
    integer ldb = max(size_t(1), b.Dist());
    integer ldc = max(size_t(1), c.Dist());
    dgemm (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }

  inline void LapackMultAtB (ngbla::FlatMatrix<double> a, 
                             ngbla::FlatMatrix<double> b,
                             ngbla::SliceMatrix<double> c)
  {
    char transa = 'N';
    char transb = 'T';
    integer m = c.Height();  // changed n,m
    integer n = c.Width(); 
    integer k = a.Height();
    double alpha = 1.0;
    double beta = 0;
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer ldc = max(size_t(1), c.Dist()); // c.Width();

    dgemm (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }

  inline void LapackMultAtBt (ngbla::FlatMatrix<double> a, 
                              ngbla::FlatMatrix<double> b,
                              ngbla::FlatMatrix<double> c)
  {
    char transa = 'T';
    char transb = 'T';
    integer m = c.Height();
    integer n = c.Width();
    integer k = a.Height();
    double alpha = 1.0;
    double beta = 0;
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer ldc = max(size_t(1), c.Width());

    dgemm (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }

  inline void LapackMultABt (ngbla::FlatMatrix<ngbla::Complex> a, 
                             ngbla::FlatMatrix<ngbla::Complex> b,
                             ngbla::FlatMatrix<ngbla::Complex> c)
  {
    //   c = a * Trans (b);

    char transa = 'T';
    char transb = 'N';
    integer m = c.Height();
    integer n = c.Width();
    integer k = a.Width();
    Complex alpha(1,0); // double alpha[2] =  { 1.0, 0.0 };
    Complex beta(0,0);  // double beta[2] = { 0.0, 0.0 };
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer ldc = max(size_t(1), c.Width());

    zgemm (&transa, &transb, &n, &m, &k, &alpha,  
            &b(0,0), &ldb, 
            &a(0,0), &lda, &beta, 
            &c(0,0), &ldc);

  }



  inline void LapackMultAddABt (ngbla::FlatMatrix<double> a, 
                                ngbla::FlatMatrix<double> b,
                                double fac,
                                ngbla::FlatMatrix<double> c)
  {
    char transa = 'T';
    char transb = 'N';
    integer m = c.Height();
    integer n = c.Width();
    integer k = a.Width();
    double alpha = fac;
    double beta = 1.0;
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer ldc = max(size_t(1), c.Width());

    dgemm (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }


  inline void LapackMultAddAB (ngbla::FlatMatrix<double> a, 
                               ngbla::FlatMatrix<double> b,
                               double fac,
                               ngbla::FlatMatrix<double> c)
  {
    char transa = 'N';
    char transb = 'N';
    integer m = c.Height();
    integer n = c.Width();
    integer k = a.Width();
    double alpha = fac;
    double beta = 1.0;
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer ldc = max(size_t(1), c.Width());

    dgemm (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }

  inline void LapackMultAddAtB (ngbla::FlatMatrix<double> a, 
				ngbla::FlatMatrix<double> b,
				double fac,
				ngbla::FlatMatrix<double> c)
  {
    char transa = 'N';
    char transb = 'T';
    integer m = c.Height();  // changed n,m
    integer n = c.Width(); 
    integer k = a.Height();
    double alpha = fac;
    double beta = 1.0;
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer ldc = max(size_t(1), c.Width());

    dgemm (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }


  inline void LapackMultAddABt (ngbla::FlatMatrix<ngbla::Complex> a, 
                                ngbla::FlatMatrix<ngbla::Complex> b,
                                double fac,
                                ngbla::FlatMatrix<ngbla::Complex> c)
  { 
    // c += fac * a * Trans (b);
    char transa = 'T';
    char transb = 'N';
    integer m = c.Height();
    integer n = c.Width();
    integer k = a.Width();
    Complex alpha(fac, 0);
    Complex beta(1,0);
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer ldc = max(size_t(1), c.Width());

    zgemm (&transa, &transb, &n, &m, &k, &alpha, 
            &b(0,0), &ldb, 
            &a(0,0), &lda, &beta, 
            &c(0,0), &ldc);
  }


 
  inline void LapackMultAddAtB (ngbla::FlatMatrix<ngbla::Complex> a, 
                                ngbla::FlatMatrix<ngbla::Complex> b,
                                double fac,
                                ngbla::FlatMatrix<ngbla::Complex> c)
  { 
    // c += fac * a * Trans (b);
    char transa = 'N';
    char transb = 'T';
    integer m = c.Height();
    integer n = c.Width();
    integer k = a.Height();
    Complex alpha(fac, 0); // double alpha[2] = { fac, 0 };
    Complex beta(1,0);     // double beta[2] = { 1.0, 0 };
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer ldc = max(size_t(1), c.Width());

    zgemm (&transa, &transb, &n, &m, &k, &alpha, 
            &b(0,0), &ldb, 
            &a(0,0), &lda, &beta, 
            &c(0,0), &ldc);
  }




  inline void LapackMultAddAB (ngbla::FlatMatrix<ngbla::Complex> a, 
                               ngbla::FlatMatrix<ngbla::Complex> b,
                               double fac,
                               ngbla::FlatMatrix<ngbla::Complex> c)
  {
    char transa = 'N';
    char transb = 'N';
    integer m = c.Height();
    integer n = c.Width();
    integer k = a.Width();
    Complex alpha(fac, 0); // double alpha[2] = { fac, 0 };
    Complex beta(1,0);    // double beta[2] = { 1.0, 0 };
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer ldc = max(size_t(1), c.Width());

    zgemm (&transa, &transb, &n, &m, &k, &alpha, 
            &b(0,0), &ldb, 
            &a(0,0), &lda, &beta, 
            &c(0,0), &ldc);
  }




  /*
  extern "C"
  void dgetri_ (int & n, double & a, int & lda,
                int & ipiv, double & work, int & lwork, int & info);

  extern "C"
  void dgetrf_ (int & n, int & m, double & a, int & lda,
                int & ipiv, int & info);

  extern "C"
  void dgetrs_ (char & trans, int & n, int & nrhs, 
                double & a, int & lda, int & ipiv, 
                double & b, int & ldb, int & info);



  extern "C"
  void dpotrf_ (char & uplo, int & n, double & a, int & lda,
                int & info);

  extern "C"
  void dpotrs_ (char & uplo, int & n, int & nrhs, 
                double & a, int & lda,
                double & b, int & ldb, int & info);
  */


  template <ORDERING ORD>
  class LapackLU
  {
    Matrix <double, ORD> a;
    ArrayMem<integer,100> ipiv;
    
  public:
    LapackLU (Matrix<double,ORD> _a)
      : a(std::move(_a)), ipiv(a.Height())
    {
      integer m = a.Height();
      if (m == 0) return;
      integer n = a.Width();
      integer lda = a.Dist();

      integer info;
      dgetrf(&n, &m, &a(0,0), &lda, &ipiv[0], &info);
    }
    
    template <typename Db>
    void Solve (VectorView<double,Db> b) const
    {
      /*
      int dgetrs_(char *trans, integer *n, integer *nrhs, 
                  doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
                  ldb, integer *info);
      */
      char transa =  (ORD == ColMajor) ? 'N' : 'T';
      integer n = a.Height();
      integer nrhs = 1;
      integer lda = a.Dist();
      integer ldb = b.Size();
      integer info;
      dgetrs(&transa, &n, &nrhs, a.Data(), &lda, ipiv.Data(), b.Data(), &ldb, &info);
    }
    
    Matrix <double,ORD> Inverse() &&
    {
      double hwork;
      integer lwork = -1;
      integer n = a.Height();      
      integer lda = a.Dist();
      integer info;
      dgetri(&n, &a(0,0), &lda, &ipiv[0], &hwork, &lwork, &info);
      lwork = integer(hwork);
      ArrayMem<double,1000> work(lwork);
      dgetri(&n, &a(0,0), &lda, &ipiv[0], &work[0], &lwork, &info);
      return std::move(a);
    }
  };


  

  inline void LapackInverse (ngbla::SliceMatrix<double> a)
  {
    integer m = a.Height();
    if (m == 0) return;
    integer n = a.Width();
    integer lda = max(size_t(1), a.Dist());

    ArrayMem<integer,100> ipiv(n);
    integer info;
    
    dgetrf(&n, &m, &a(0,0), &lda, &ipiv[0], &info);

    double hwork;
    integer lwork = -1;
    dgetri(&n, &a(0,0), &lda, &ipiv[0], &hwork, &lwork, &info);
    lwork = integer(hwork);

    ArrayMem<double,1000> work(lwork);
    dgetri(&n, &a(0,0), &lda, &ipiv[0], &work[0], &lwork, &info);
  }



  inline void LapackInverseSPD (ngbla::SliceMatrix<double> a)
  {
    integer n = a.Width();
    if (n == 0) return;
    integer lda = max(size_t(1), a.Dist());

    integer info;
    char uplo = 'U';

    dpotrf_ (&uplo, &n, &a(0,0), &lda, &info);
    dpotri_ (&uplo, &n, &a(0,0), &lda, &info);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < i; j++)
	a(j,i) = a(i,j);
  }









  /*
    Compoutes B <--- B A^{-1}    (trans = 'N')
    Compoutes B <--- B A^{-T}    (trans = 'T')
    Compoutes B <--- B A^{-H}    (trans = 'H')

    // trans = 'N': solve A x = b^T
    // trans = 'T': solve A^T x = b^T
    // trans = 'C': solve A^H x = b^T
  */
  inline void LapackAInvBt (ngbla::FlatMatrix<double> a, ngbla::FlatMatrix<double> b, char trans = 'N')
  {
    integer m = a.Height();
    integer n = a.Width();
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer nrhs = b.Height();

    ArrayMem<integer,100> ipiv(n);
    integer info;
    // char uplo = 'L';

    dgetrf_ (&n, &m, &a(0,0), &lda, &ipiv[0], &info);
    dgetrs_ (&trans, &n, &nrhs, &a(0,0), &lda, &ipiv[0], &b(0,0), &ldb, &info);

    /*
    // symmetric, non-spd
    int lwork = -1;
    double hwork[1] = { 1.0 };
    dsytrf_ (&uplo, &n, &a(0,0), &lda, &ipiv[0], &hwork[0], &lwork, &info);
    lwork = int(hwork[0]);
    ArrayMem<double, 1000> work(lwork);
    dsytrf_ (&uplo, &n, &a(0,0), &lda, &ipiv[0], &work[0], &lwork, &info);
    dsytrs_ (&uplo, &n, &nrhs, &a(0,0), &lda, &ipiv[0], &b(0,0), &ldb, &info);
    */

    /*
    // spd
    dpotrf_ (&uplo, &n, &a(0,0), &lda, &info);
    dpotrs_ (&uplo, &n, &nrhs, &a(0,0), &lda, &b(0,0), &ldb, &info);
    */
  }




  /*
  extern "C"
  void zgetri_ (int & n, double & a, int & lda,
                int & ipiv, double & work, int & lwork, int & info);

  extern "C"
  void zgetrf_ (int & n, int & m, double & a, int & lda,
                int & ipiv, int & info);

  extern "C"
  void zgetrs_ (char & trans, int & m, int & nrhs, double & a, int & lda,
                int & ipiv, 
                double & b, int & ldb, 
                int & info);
  */


  inline void LapackInverse (ngbla::SliceMatrix<ngbla::Complex> a)
  {
    integer m = a.Height();
    if (m == 0) return;

    integer n = a.Width();
    integer lda = a.Dist();
    integer * ipiv = new integer[n];
    integer lwork = 100*n;
    Complex * work = new Complex[lwork];
    integer info;
  
    // std::cout << "a = " << std::endl << a << std::endl;
    zgetrf_ (&n, &m, &a(0,0), &lda, ipiv, &info);
    // std::cout << "factors = " << std::endl << a << std::endl;

    if (info != 0)
      {
	std::cout << "ZGETRF::info = " << info << std::endl;
	// *testout << "ZGETRF::info = " << info << std::endl;
	// *testout << "a = " << endl << a << endl;
      }
    zgetri_ (&n,  &a(0,0), &lda, ipiv, work, &lwork, &info);
    if (info != 0)
      std::cout << "ZGETRI::info = " << info << std::endl;
  
    delete [] work;
    delete [] ipiv;

  }


  /*
    trans = 'N': solve A x = b^T
    trans = 'T': solve A^T x = b^T
    trans = 'C': solve A^H x = b^T
  */

  inline void LapackAInvBt (ngbla::FlatMatrix<ngbla::Complex> a, ngbla::FlatMatrix<ngbla::Complex> b, char trans = 'N')
  {
    integer m = a.Height();
    integer n = a.Width();
    integer lda = max(size_t(1), a.Width());
    integer ldb = max(size_t(1), b.Width());
    integer nrhs = b.Height();
    integer * ipiv = new integer[n];
    integer lwork = 100*n;
    double * work = new double[2*lwork];
    integer info;
  

    zgetrf_ (&n, &m,&a(0,0), &lda, ipiv, &info);
    zgetrs_ (&trans, &n, &nrhs, &a(0,0), &lda, ipiv, &b(0,0), &ldb, &info);

    delete [] work;
    delete [] ipiv;
    //   std::cerr << "complex LapackAInvBt not implemented" << std::endl;
  }


  /*
  extern "C"
  void dsyev_(char & jobz, char & uplo, int & n , double & A , int & lda, double & w,
              double & work, int & lwork, int & info); 

  extern "C"
  void zgeev_( char *jobvl, char *jobvr, int *n, std::complex<double> *A, int * lda, std::complex<double>* lami, 
               std::complex<double> * vl, int * nvl, std::complex<double> * vr, int * nvr, 
               std::complex<double> * work, int * lwork, double * rwork, int * info);
  */

  //  extern "C"
  // void dgeev_( char *jobvl, char *jobvr, int *n, double *A, int * lda, double* lami_re, double *lami_im, 
  // double * vl, int * nvl, double * vr, int * nvr, 
  // double * work, int * lwork,  /* double * rwork, */ int * info);



  NGS_DLL_HEADER
  void LapackEigenValuesSymmetric (ngbla::FlatMatrix<double> a,
                                   ngbla::FlatVector<double> lami,
                                   ngbla::FlatMatrix<double> evecs = ngbla::FlatMatrix<double>(0,0));
  /*
  {
    char jobz, uplo = 'U'; 
    integer n = a.Height();
    integer lwork=(n+2)*n+1;
 
    double* work = new double[lwork];
    integer info; 
 
    double * matA;

    if ( evecs.Height() )
      {
        // eigenvectors are calculated
        evecs = a;
        jobz = 'V';
        matA = &evecs(0,0);
      }
    else
      {
        // only eigenvalues are calculated, matrix a is destroyed!!
        jobz = 'N';
        matA = &a(0,0);
      }
    dsyev_(&jobz, &uplo , &n , matA, &n, &lami(0), work, &lwork, &info); 

    if (info)
      std::cerr << "LapackEigenValuesSymmetric, info = " << info << std::endl;

    delete [] work; 
  }
  */



  inline void LapackEigenValues (ngbla::FlatMatrix<double> a,
                                 ngbla::FlatVector<ngbla::Complex> lami, 
                                 ngbla::FlatMatrix<double> eveci )
  {
    char jobvr = 'V' , jobvl= 'N';

    integer n = a.Height();
    integer nvl = 1; 
    integer nvr = eveci.Width() ; 
  
    double * vl = 0; 
    double * vr;//  = new std::complex<double> [nvr*n];
    double * lami_re = new double[n], * lami_im = new double[n];

    integer lwork = 8*n; 
    double * work = new double[lwork]; 
    double *rwork = new double[8*n];  
    integer info = 0;
  
    if ( eveci.Width() )
      {
        vr = &eveci(0,0);
      }
    else
      {
        nvr = n;
        vr =  new double [nvr*n];
      }

    dgeev_(&jobvl, &jobvr, &n, &a(0,0), &n, lami_re, lami_im, vl, &nvl, vr, &nvr, work, &lwork, /* rwork, */ &info);
  
    if(info != 0) 	
      {
        std::cout << "**** Error in zggev_, info = " << info << " *****" << std::endl; 
        return;
      }
    
    for ( size_t i = 0; i < lami.Size(); i++ )
      lami(i) = ngbla::Complex (lami_re[i], lami_im[i]);

    delete[] work; 
    delete[] rwork;
    if ( !eveci.Width() )
      delete[] vr;
    delete [] lami_re;
    delete [] lami_im;
  }

  inline void LapackEigenValues (ngbla::FlatMatrix<ngbla::Complex> a,
                                 ngbla::FlatVector<ngbla::Complex> lami, 
                                 ngbla::FlatMatrix<ngbla::Complex> eveci )
  {
    char jobvr = 'V' , jobvl= 'N';

    integer n = a.Height();
    integer nvl = 1; 
    integer nvr = eveci.Width() ; 
  
    Complex * vl = 0;
    Complex * vr;//  = new std::complex<double> [nvr*n];
  
    integer lwork = 8*n; 
    Complex * work = new Complex [lwork];
    double *rwork = new double[8*n];  
    integer info = 0;
  
    if ( eveci.Width() )
      {
        vr = &eveci(0,0);
      }
    else
      {
        nvr = n;
        vr =  new Complex [nvr*n];
      }

    zgeev_(&jobvl, &jobvr, &n, &a(0,0), &n, &lami(0), vl, &nvl, vr, &nvr, work, &lwork, rwork, &info);
    //  alpha, beta, &vl, &nvl, vr, &nvr,  
    // 	     work , &lwork, rwork,  &info);
  
    if(info != 0) 	
      {
        std::cout << "**** Error in zggev_, info = " << info << " *****" << std::endl; 
        return;
      }
    
    delete[] work; 
    delete[] rwork;
    if ( !eveci.Width() )
      delete[] vr;
  }





  inline void LapackEigenValuesSymmetric (ngbla::FlatMatrix<ngbla::Complex> a,
                                          ngbla::FlatVector<ngbla::Complex> lami, 
                                          ngbla::FlatMatrix<ngbla::Complex> eveci = ngbla::FlatMatrix<ngbla::Complex> (0,0) )
  {
    //   std::cerr << "complex evp not implemented" << std::endl;
    LapackEigenValues ( a, lami, eveci );
  }



  // Solve complex generalized eigenvalue problem (QZ) 
  /*
  extern "C"
  void zggev_(char *jobvl,char* jobvr,int* N,std::complex<double>* A, int* lda,std::complex<double>* B, int* ldb, std::complex<double>* alpha, std::complex<double>* beta, std::complex<double>* vl, int* ldvl, std::complex<double>* vr, int* ldvr, std::complex<double>* work, int* lwork, double* rwork, int* info); 
  */


  // Attention A,B are overwritten !!! 
  inline void LapackEigenValues (ngbla::FlatMatrix<ngbla::Complex> a,
                                 ngbla::FlatMatrix<ngbla::Complex> b,
                                 ngbla::FlatVector<ngbla::Complex> lami)

  {
    integer n = a.Height();
    // integer evecs_bool = 0;
    // std::complex<double> * evecs, * dummy;

    char jobvr = 'N', jobvl= 'N';
    // bool balancing = 0; 
  
    Complex * alpha= new Complex[n];
    Complex * beta = new Complex[n];
    Complex vl=0.;
  
    integer nvl = 1; 
    Complex * vr = NULL;
  
    Complex * work = new Complex[8*n];
    integer lwork = 8*n; 
    double *rwork = new double[8*n];  
  
    integer nvr = n ; 
  
    //std::complex<double>  * A1,*B1; 
    integer i; 
  
    // char job=balance_type; // Permute and Scale in Balancing 
    // integer ihi,ilo; 
    double * lscale, *rscale; 
    lscale = new double[n]; 
    rscale = new double[n]; 
    double * work2; 
    work2 = new double[6*n];
  
    // char side = 'R'; 
  
    integer info = 0;

    // integer ii; 
    // ii=0; 
   
    // if(balancing) zggbal_(&job,&n, A, &n , B, &n, &ilo, &ihi,  lscale, rscale, work2, &info) ; 
  
    // if(info == 0 ) 
    if (1)
      {  
        zggev_(&jobvl, &jobvr, &n, &a(0,0), &n, &b(0,0), &n, alpha, beta, &vl, &nvl, vr, &nvr,  
               work , &lwork, rwork,  &info);
      
        if(info==0) 	
          {
            /*
              if(jobvr == 'V' && balancing) 
              {
	      zggbak_(&job, &side, &n, &ilo, &ihi, lscale, rscale, &n, vr, &n,&info)  ;
	      
	      if(info!=0)
              { 
              std::cout << "*****  Error in zggbak_ **** " << endl; 
              return; 
              }
              }
            */
          } 
        else 
          {
            std::cout << "**** Error in zggev_, info = " << info << " *****" << std::endl; 
            // return;
          }
      }	
    else 
      {
        std::cout << "**** Error in zggbal_ **** " << std::endl; 
        return;
      }
    
    delete [] work; 
    delete [] rwork;
  
    delete [] lscale; 
    delete [] rscale; 
    delete [] work2;   
 
    for(i=0;i<n;i++)
      {
        if(abs(beta[i]) >= 1.e-30) 
          lami[i]=Complex(alpha[i]/beta[i]);
        else 
          {
            lami[i] = Complex(100.,100.);
          }
      } 
  
    /*  
        std:: complex<double> resid[n]; 
        double error;  
        for(int k=0; k<n;k++)
        {
	for(i=0; i<n;i++) 
        { 
        resid[i] = 0.; 
        for(int j=0;j<n;j++)
        {
        resid[i] += A[j*n + i]*evecs[k*n+j];
        resid[i] -=B[j*n+i]*evecs[k*n+j]*lami[k];
        }
    
        error = abs(resid[i]) * abs(resid[i]) ;
        } 
	error = sqrt(error); 
	cout << " lami (i) " << lami[k] <<  "\t" <<  alpha[k] << "\t"  << beta[k]  << endl; 
	cout << " error lapack " << k << "  \t " << error << endl; 
        }
    */
    delete [] alpha; 
    delete [] beta; 
 
  }   










  /*
  extern "C"
  void dsygv_(int & itype, char & jobzm, char & uplo , int & n1, double & A , int & n, 
              double & B, int  & nb, double & lami, double & work, int & lwork, int & info); 
  */

  inline void LapackEigenValuesSymmetric (ngbla::FlatMatrix<double> a,
                                          ngbla::FlatMatrix<double> b,
                                          ngbla::FlatVector<double> lami,
                                          ngbla::FlatMatrix<double> evecs = ngbla::FlatMatrix<double>(0,0))
  {
    char jobz = 'N' , uplo = 'U'; 
    integer n = a.Height();

    integer lwork=(n+2)*n+1;
    double* work = new double[lwork];
  
    integer info; 
    integer itype =1; 

    if ( evecs.Height() )
      jobz = 'V';
    else
      jobz = 'N';

    dsygv_(&itype, &jobz, &uplo , &n , &a(0,0), &n, &b(0,0), &n, 
           &lami(0), work, &lwork, &info); 

    if ( evecs.Height() )
      evecs = a;

    if (info)
      std::cerr << "LapackEigenValuesSymmetric, info = " << info << std::endl;
  
    delete [] work; 
  }

  // A = U * diag(S) * V
  NGS_DLL_HEADER void LapackSVD (SliceMatrix<double, ColMajor> A,
                         SliceMatrix<double, ColMajor> U,
                         SliceMatrix<double, ColMajor> V,
                         FlatVector<double> S,                         
                         bool all);


  // A = U * diag(S) * V
  inline void LapackSVD (SliceMatrix<double> A,
                         SliceMatrix<double> U,
                         SliceMatrix<double> V,
                         FlatVector<double> S,                         
                         bool all)
  {
    LapackSVD (Trans(A), Trans(V), Trans(U), S, all);
  }
  

  // A = U * diag(S) * V
  NGS_DLL_HEADER void LapackSVD (SliceMatrix<Complex, ColMajor> A,
                         SliceMatrix<Complex, ColMajor> U,
                         SliceMatrix<Complex, ColMajor> V,
                         FlatVector<double> S,                         
                         bool all);


  // A = U * diag(S) * V
  inline void LapackSVD (SliceMatrix<Complex> A,
                         SliceMatrix<Complex> U,
                         SliceMatrix<Complex> V,
                         FlatVector<double> S,                         
                         bool all)
  {
    LapackSVD (Trans(A), Trans(V), Trans(U), S, all);
  }
  
  

  

#else

  typedef int integer;

  inline void LapackMultAtx (ngbla::FlatMatrix<double> a,
                             ngbla::FlatVector<double> x,
                             ngbla::FlatVector<double> y)
  { y = Trans (a) * x; }


  inline void LapackAddxyt (ngbla::FlatMatrix<double> a,
                            ngbla::FlatVector<double> x,
                            ngbla::FlatVector<double> y)
  { y += Trans (a) * x; }


  template <typename TA, typename TB>
  inline void LapackMult (const TA & a, const TB & b, 
			  ngbla::SliceMatrix<double> c)
  { c = a * b; }

  template <typename TA, typename TB>
  inline void LapackMult (const TA & a, const TB & b, 
			  ngbla::SliceMatrix<Complex> c)
  { c = a * b; }

  template <typename TA, typename TB>
  inline void LapackMultAdd (const TA & a,
			     const TB & b, 
			     double alpha,
			     SliceMatrix<double> c,
			     double beta)
  { c *= beta; c += alpha * a * b; }

  template <typename TA, typename TB>
  inline void LapackMultAdd (const TA & a,
			     const TB & b, 
			     Complex alpha,
			     SliceMatrix<Complex> c,
			     Complex beta)
  { c *= beta; c += alpha * a * b; }



  inline void LapackMultABt (ngbla::FlatMatrix<double> a, 
                             ngbla::FlatMatrix<double> b,
                             ngbla::FlatMatrix<double> c)
  { c = a * Trans (b); }

  inline void LapackMultAtB (ngbla::FlatMatrix<double> a, 
                             ngbla::FlatMatrix<double> b,
                             ngbla::FlatMatrix<double> c)
  { c = Trans(a) * b; }


  inline void LapackMultAB (ngbla::FlatMatrix<double> a, 
                            ngbla::FlatMatrix<double> b,
                            ngbla::FlatMatrix<double> c)
  { c = a * b; }


  inline void LapackMultABt (ngbla::FlatMatrix<ngbla::Complex> a, 
                             ngbla::FlatMatrix<ngbla::Complex> b,
                             ngbla::FlatMatrix<ngbla::Complex> c)
  { c = a * Trans (b); }

  inline void LapackMultAtB (ngbla::FlatMatrix<Complex> a, 
                             ngbla::FlatMatrix<Complex> b,
                             ngbla::FlatMatrix<Complex> c)
  { c = Trans(a) * b; }



  inline void LapackMultAddAB (ngbla::FlatMatrix<double> a, 
                                ngbla::FlatMatrix<double> b,
                                double fac,
                                ngbla::FlatMatrix<double> c)
  { c += fac * a * b; }

  inline void LapackMultAddABt (ngbla::FlatMatrix<double> a, 
                                ngbla::FlatMatrix<double> b,
                                double fac,
                                ngbla::FlatMatrix<double> c)
  { c += fac * a * Trans (b); }

  inline void LapackMultAddAtB (ngbla::FlatMatrix<double> a, 
                                ngbla::FlatMatrix<double> b,
                                double fac,
                                ngbla::FlatMatrix<double> c)
  { c += fac * Trans(a) * b; }



  inline void LapackMultAddAB (ngbla::FlatMatrix<ngbla::Complex> a, 
                                ngbla::FlatMatrix<ngbla::Complex> b,
                                double fac,
                                ngbla::FlatMatrix<ngbla::Complex> c)

  { c += fac * a * b; }

  inline void LapackMultAddABt (ngbla::FlatMatrix<ngbla::Complex> a, 
                                ngbla::FlatMatrix<ngbla::Complex> b,
                                double fac,
                                ngbla::FlatMatrix<ngbla::Complex> c)

  { c += fac * a * Trans (b); }


  
  inline void LapackInverse (ngbla::FlatMatrix<double> a)
  { 
    CalcInverse (a);
    /*
    ngbla::Matrix<> hm(a.Height());
    CalcInverse (a, hm);
    a = hm;
    */
  }

  inline void LapackInverse (ngbla::FlatMatrix<ngbla::Complex> a)
  { 
    CalcInverse (a);
    /*
    // std::cerr << "sorry, Inverse not available without LAPACK" << std::endl;
    ngbla::Matrix<Complex> hm(a.Height());
    CalcInverse (a, hm);
    a = hm;
    */
  }

  inline void LapackAInvBt (ngbla::FlatMatrix<double> a, ngbla::FlatMatrix<double> b, char trans = 'N')
  {
    LapackInverse (a);
    ngbla::Matrix<> hb (b.Height(), b.Width());
    if (trans == 'T')
      hb = b * Trans(a);
    else
      hb = b * a;
    b = hb;
  }

  inline void LapackAInvBt (ngbla::FlatMatrix<Complex> a, ngbla::FlatMatrix<Complex> b, char trans = 'N')
  {
    LapackInverse (a);
    ngbla::Matrix<Complex> hb (b.Height(), b.Width());
    if (trans == 'T')
      hb = b * Trans(a);
    else
      hb = b * a;
    b = hb;
  }


  inline void LapackEigenValuesSymmetric (ngbla::FlatMatrix<double> a,
                                          ngbla::FlatVector<double> lami)
  { 
    std::cerr << "sorry, EVP not available without LAPACK" << std::endl;
  }

  inline void LapackEigenValuesSymmetric (ngbla::FlatMatrix<double> a,
                                          ngbla::FlatMatrix<double> b,
                                          ngbla::FlatVector<double> lami)
  { 
    std::cerr << "sorry, EVP not available without LAPACK" << std::endl;
  }


  inline void LapackEigenValuesSymmetric (ngbla::FlatMatrix<ngbla::Complex> a,
                                          ngbla::FlatVector<ngbla::Complex> lami)
  {
    std::cerr << "sorry, EVP not available without LAPACK" << std::endl;
  }


#endif




  // several LAPACK eigenvalue solvers  

#ifdef LAPACK

  void LaEigNSSolve(int n, double * A, double * B, std::complex<double> * lami, int evecs_bool, double *evecs_re, double *evecs_im, char balance_type);
  void LaEigNSSolve(int n, std::complex<double> * A, std::complex<double> * B, std::complex<double> * lami, int evecs_bool, std::complex<double> *evecs, std::complex<double> *dummy, char balance_type);

  void LapackSSEP(int n, double* A, double* lami, double* evecs); 

  void LapackHessenbergEP (int n, std::complex<double> * H, std::complex<double> * lami, std::complex<double> * evecs); 
  void LapackGHEP(int n, double* A, double* B,  double* lami) ; 
  int LapackGHEPEPairs(int n, double* A, double* B, double* lami);
  int LapackGHEPEPairs(int n, std::complex<double>* A, std::complex<double>* B, double* lami);
  // A,B overwritten in A eigenvectors z^H B z = 1  
  
  //void LaEigNSSolve(const LaGenMatDouble &A, LaVectorDouble &eigvals);
  //void LaEigNSSolveIP(LaGenMatDouble &A, LaVectorDouble &eigvals);

  void LaEigNSSolveTest();
  void LaLinearSolveComplex(int n, std::complex<double> * A, std::complex<double> * F); 
  void LaLinearSolve(int n, double * A, double * F); 
  void LaLinearSolveRHS(int n, double * A, double * F); 


  void LaEigNSSolveX(int n, std::complex<double> * A, std::complex<double> * B, std::complex<double> * lami, int evecs_bool, std::complex<double> * evecs, std::complex<double> * dummy, char balance_type); 
  void LaEigNSSolveX(int n, double * A, double * B, std::complex<double> * lami, int evecs_bool, double * evecs, double * dummy, char balance_type);

#else
  
  inline void LapackHessenbergEP (int n, std::complex<double> * H, std::complex<double> * lami, std::complex<double> * evecs)
  {
    cerr << "Sorry, HessebergEP not available without Lapack" << endl;
  }
  

#endif

}


#endif
