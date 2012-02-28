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
  // Interface to lapack functions


#ifdef LAPACK

  extern "C" {
    typedef char logical;
    typedef int integer;
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


  
  // BLAS 1

  inline double LapackDot (FlatVector<double> x, FlatVector<double> y)
  {
    int n = x.Size();
    int incx = 1;
    int incy = 1;
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
    int m = a.Width();
    int n = a.Height();
    double alpha = 1;
    int lda = a.Width();
    int incx = 1;
    double beta = 0;
    int incy = 1;
    dgemv_ (&trans, &m, &n, &alpha, &a(0,0), &lda, &x(0), &incx, &beta, &y(0), &incy);
  }

  inline void LapackMultAx (ngbla::FlatMatrix<Complex> a,
                            ngbla::FlatVector<Complex> x,
                            ngbla::FlatVector<Complex> y)
  {
    char trans = 'T';
    int m = a.Width();
    int n = a.Height();
    Complex alpha(1,0);
    int lda = a.Width();
    int incx = 1;
    Complex beta(0, 0);
    int incy = 1;
    zgemv_ (&trans, &m, &n, &alpha, &a(0,0), &lda, &x(0), &incx, &beta, &y(0), &incy);
  }



  inline void LapackMultAtx (ngbla::FlatMatrix<double> a,
                             ngbla::FlatVector<double> x,
                             ngbla::FlatVector<double> y)
  {
    char trans = 'N';
    int m = a.Width();
    int n = a.Height();
    double alpha = 1;
    int lda = a.Width();
    int incx = 1;
    double beta = 0;
    int incy = 1;
    dgemv_ (&trans, &m, &n, &alpha, &a(0,0), &lda, &x(0), &incx, &beta, &y(0), &incy);
  }


  /*
  extern "C"
  void dger_ (int & m, int & n, double & alpha, double & x, int & incx, double & y, int & incy, double & a, int & lda);
  */

  inline void LapackAddxyt (ngbla::FlatMatrix<double> a,
                            double fac,
                            ngbla::FlatVector<double> x,
                            ngbla::FlatVector<double> y)
  {
    int m = a.Width();
    int n = a.Height();
    double alpha = fac;
    int lda = a.Width();
    int incx = 1;
    int incy = 1;

    dger_ (&m, &n, &alpha, &y(0), &incx, &x(0), &incy, &a(0,0), &lda);
  }






  inline void LapackMultABt (FlatMatrix<double> a, FlatMatrix<double> b, FlatMatrix<double> c)
  {
    char transa = 'T';
    char transb = 'N';
    int m = c.Height();
    int n = c.Width();
    int k = a.Width();
    double alpha = 1.0;
    double beta = 0;
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();

    dgemm_ (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }
  

  inline void LapackMultAB (FlatMatrix<double> a, FlatMatrix<double> b, FlatMatrix<double> c)
  {
    char transa = 'N';
    char transb = 'N';
    int m = c.Height();
    int n = c.Width();
    int k = a.Width();
    double alpha = 1.0;
    double beta = 0;
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();
    dgemm_ (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }


  inline void LapackMultAtB (ngbla::FlatMatrix<double> a, 
                             ngbla::FlatMatrix<double> b,
                             ngbla::FlatMatrix<double> c)
  {
    char transa = 'N';
    char transb = 'T';
    int m = c.Height();  // changed n,m
    int n = c.Width(); 
    int k = a.Height();
    double alpha = 1.0;
    double beta = 0;
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();

    dgemm_ (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }


  inline void LapackMultAtBt (ngbla::FlatMatrix<double> a, 
                              ngbla::FlatMatrix<double> b,
                              ngbla::FlatMatrix<double> c)
  {
    char transa = 'T';
    char transb = 'T';
    int m = c.Height();
    int n = c.Width();
    int k = a.Height();
    double alpha = 1.0;
    double beta = 0;
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();

    dgemm_ (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }








  inline void LapackMultABt (ngbla::FlatMatrix<ngbla::Complex> a, 
                             ngbla::FlatMatrix<ngbla::Complex> b,
                             ngbla::FlatMatrix<ngbla::Complex> c)
  {
    //   c = a * Trans (b);

    char transa = 'T';
    char transb = 'N';
    int m = c.Height();
    int n = c.Width();
    int k = a.Width();
    Complex alpha(1,0); // double alpha[2] =  { 1.0, 0.0 };
    Complex beta(0,0);  // double beta[2] = { 0.0, 0.0 };
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();

    zgemm_ (&transa, &transb, &n, &m, &k, &alpha,  
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
    int m = c.Height();
    int n = c.Width();
    int k = a.Width();
    double alpha = fac;
    double beta = 1.0;
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();

    dgemm_ (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }


  inline void LapackMultAddAB (ngbla::FlatMatrix<double> a, 
                               ngbla::FlatMatrix<double> b,
                               double fac,
                               ngbla::FlatMatrix<double> c)
  {
    char transa = 'N';
    char transb = 'N';
    int m = c.Height();
    int n = c.Width();
    int k = a.Width();
    double alpha = fac;
    double beta = 1.0;
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();

    dgemm_ (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }

  inline void LapackMultAddAtB (ngbla::FlatMatrix<double> a, 
				ngbla::FlatMatrix<double> b,
				double fac,
				ngbla::FlatMatrix<double> c)
  {
    char transa = 'N';
    char transb = 'T';
    int m = c.Height();  // changed n,m
    int n = c.Width(); 
    int k = a.Height();
    double alpha = fac;
    double beta = 1.0;
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();

    dgemm_ (&transa, &transb, &n, &m, &k, &alpha, &b(0,0), &ldb, &a(0,0), &lda, &beta, &c(0,0), &ldc);
  }


  inline void LapackMultAddABt (ngbla::FlatMatrix<ngbla::Complex> a, 
                                ngbla::FlatMatrix<ngbla::Complex> b,
                                double fac,
                                ngbla::FlatMatrix<ngbla::Complex> c)
  { 
    // c += fac * a * Trans (b);
    char transa = 'T';
    char transb = 'N';
    int m = c.Height();
    int n = c.Width();
    int k = a.Width();
    Complex alpha(fac, 0);
    Complex beta(1,0);
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();

    zgemm_ (&transa, &transb, &n, &m, &k, &alpha, 
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
    int m = c.Height();
    int n = c.Width();
    int k = a.Height();
    Complex alpha(fac, 0); // double alpha[2] = { fac, 0 };
    Complex beta(1,0);     // double beta[2] = { 1.0, 0 };
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();

    zgemm_ (&transa, &transb, &n, &m, &k, &alpha, 
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
    int m = c.Height();
    int n = c.Width();
    int k = a.Width();
    Complex alpha(fac, 0); // double alpha[2] = { fac, 0 };
    Complex beta(1,0);    // double beta[2] = { 1.0, 0 };
    int lda = a.Width();
    int ldb = b.Width();
    int ldc = c.Width();

    zgemm_ (&transa, &transb, &n, &m, &k, &alpha, 
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



  inline void LapackInverse (ngbla::FlatMatrix<double> a)
  {
    int m = a.Height();
    if (m == 0) return;
    int n = a.Width();
    int lda = a.Width();
    int * ipiv = new int[n];
    int lwork = 100*n;
    double * work = new double[lwork];
    int info;

    dgetrf_ (&n, &m, &a(0,0), &lda, ipiv, &info);
    dgetri_ (&n, &a(0,0), &lda, ipiv, work, &lwork, &info);

    delete [] work;
    delete [] ipiv;
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
    int m = a.Height();
    int n = a.Width();
    int lda = a.Width();
    int ldb = b.Width();
    int nrhs = b.Height();
    //int * ipiv = new int[n];
    ArrayMem<int,100> ipiv(n);
    int info;
    // char uplo = 'L';

    dgetrf_ (&n, &m, &a(0,0), &lda, &ipiv[0], &info);
    dgetrs_ (&trans, &n, &nrhs, &a(0,0), &lda, &ipiv[0], &b(0,0), &ldb, &info);

    /*
      dpotrf_ (uplo, n, a(0,0), lda, info);
      dpotrs_ (uplo, n, nrhs, a(0,0), lda, b(0,0), ldb, info);
    */

    // delete [] ipiv;
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


  inline void LapackInverse (ngbla::FlatMatrix<ngbla::Complex> a)
  {
    int m = a.Height();
    if (m == 0) return;

    int n = a.Width();
    int lda = a.Width();
    int * ipiv = new int[n];
    int lwork = 100*n;
    Complex * work = new Complex[lwork];
    int info;
  
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
    int m = a.Height();
    int n = a.Width();
    int lda = a.Width();
    int ldb = b.Width();
    int nrhs = b.Height();
    int * ipiv = new int[n];
    int lwork = 100*n;
    double * work = new double[2*lwork];
    int info;
  

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



  inline void LapackEigenValuesSymmetric (ngbla::FlatMatrix<double> a,
                                          ngbla::FlatVector<double> lami,
                                          ngbla::FlatMatrix<double> evecs = ngbla::FlatMatrix<double>(0,0)){
    char jobz, uplo = 'U'; 
    int n = a.Height();
    int lwork=(n+2)*n+1;
 
    double* work = new double[lwork];
    int info; 
 
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



  inline void LapackEigenValues (ngbla::FlatMatrix<double> a,
                                 ngbla::FlatVector<ngbla::Complex> lami, 
                                 ngbla::FlatMatrix<double> eveci )
  {
    char jobvr = 'V' , jobvl= 'N';

    int n = a.Height();
    int nvl = 1; 
    int nvr = eveci.Width() ; 
  
    double * vl = 0; 
    double * vr;//  = new std::complex<double> [nvr*n];
    double * lami_re = new double[n], * lami_im = new double[n];

    int lwork = 8*n; 
    double * work = new double[lwork]; 
    double *rwork = new double[8*n];  
    int info = 0;
  
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
    
    for ( int i = 0; i < lami.Size(); i++ )
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

    int n = a.Height();
    int nvl = 1; 
    int nvr = eveci.Width() ; 
  
    std::complex<double> * vl = 0; 
    std::complex<double> * vr;//  = new std::complex<double> [nvr*n];
  
    int lwork = 8*n; 
    std::complex<double> * work = new std::complex<double>[lwork]; 
    double *rwork = new double[8*n];  
    int info = 0;
  
    if ( eveci.Width() )
      {
        vr = &eveci(0,0);
      }
    else
      {
        nvr = n;
        vr =  new std::complex<double> [nvr*n];
      }

    zgeev_(&jobvl, &jobvr, &n, (std::complex<double>*)(void*)&a(0,0), &n, (std::complex<double>*)(void*)&lami(0), vl, &nvl, vr, &nvr, work, &lwork, rwork, &info);
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
    int n = a.Height();
    // int evecs_bool = 0;
    // std::complex<double> * evecs, * dummy;

    char jobvr = 'N', jobvl= 'N';
    // bool balancing = 0; 
  
    std::complex<double> * alpha= new std::complex<double>[n];
    std::complex<double> * beta = new std::complex<double>[n]; 
    std::complex<double> vl=0.; 
  
    int nvl = 1; 
    std::complex<double> * vr ;
  
    std::complex<double> * work = new std::complex<double>[8*n]; 
    int lwork = 8*n; 
    double *rwork = new double[8*n];  
  
    int nvr = n ; 
  
    //std::complex<double>  * A1,*B1; 
    int i; 
  
    // char job=balance_type; // Permute and Scale in Balancing 
    // int ihi,ilo; 
    double * lscale, *rscale; 
    lscale = new double[n]; 
    rscale = new double[n]; 
    double * work2; 
    work2 = new double[6*n];
  
    // char side = 'R'; 
  
    int info = 0;

    // int ii; 
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
          lami[i]=std::complex<double>(alpha[i]/beta[i]);     
        else 
          {
            lami[i] = std::complex<double>(100.,100.);
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
    int n = a.Height();

    int lwork=(n+2)*n+1;
    double* work = new double[lwork];
  
    int info; 
    int itype =1; 

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


#else


  inline void LapackMultAtx (ngbla::FlatMatrix<double> a,
                             ngbla::FlatVector<double> x,
                             ngbla::FlatVector<double> y)
  { y = Trans (a) * x; }


  inline void LapackAddxyt (ngbla::FlatMatrix<double> a,
                            ngbla::FlatVector<double> x,
                            ngbla::FlatVector<double> y)
  { y += Trans (a) * x; }




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
    ngbla::Matrix<> hm(a.Height());
    CalcInverse (a, hm);
    a = hm;
  }

  inline void LapackInverse (ngbla::FlatMatrix<ngbla::Complex> a)
  { 
    // std::cerr << "sorry, Inverse not available without LAPACK" << std::endl;
    ngbla::Matrix<Complex> hm(a.Height());
    CalcInverse (a, hm);
    a = hm;
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
