#ifdef LAPACK
// #include <iostream>
// #include <complex>
// using namespace std;

#include "bla.hpp"
// #include "LapackInterface.hpp" 

/*
  #include <mkl_lapack.h>
  #define dggev_ dggev
  #define dsyev_ dsyev
  #define dggbak_ dggbak
  #define dggbal_ dggbal
  #define zggev_ zggev
  #define zggbal_ zggbal
  #define zggev_ zggev
*/

// #include "LapackGEP.hpp" 




namespace ngbla
{


  int sgemm(char *transa, char *transb, integer *m, integer *
		  n, integer *k, real *alpha, real *a, integer *lda, 
		  real *b, integer *ldb, real *beta, real *c__, 
		  integer *ldc)
  {
    return sgemm_ (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
  }


  int dgemm(char *transa, char *transb, integer *m, integer *
		  n, integer *k, doublereal *alpha, doublereal *a, integer *lda, 
		  doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
		  integer *ldc)
  {
    return dgemm_ (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
  }

  int zgemm(char *transa, char *transb, integer *m, integer *
		    n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, 
		    doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *
		    c__, integer *ldc)
  {
    return zgemm_ (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c__, ldc);
  }

  int dger(integer *m, integer *n, doublereal *alpha,
                   doublereal *x, integer *incx, doublereal *y, integer *incy,
                   doublereal *a, integer *lda)
  {
    return dger_ (m, n, alpha, y, incx, x, incy, a, lda);
  }

  int dgetri(integer* n, double* a, integer* lda, integer* ipiv,
             double* hwork, integer* lwork, integer* info)
  {
    return dgetri_(n,a,lda,ipiv,hwork,lwork,info);
  }

  int dgetrf(integer* n, integer* m, double* a, integer* lda, integer* ipiv, integer* info)
  {
    return dgetrf_(n,m,a,lda,ipiv,info);
  }

  int dgetrs(char *trans, integer *n, integer *nrhs, 
             doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
             ldb, integer *info)
  {
    return dgetrs_(trans, n, nrhs, a, lda, ipiv, b, ldb, info);
  }
  
  
}



namespace ngbla
{

  void LaEigNSSolveTest()
  {   

    double *  A= new double[16];
    double *  B= new double[16]; 
 
    
    int i;   
    for(i=0;i<16;i++) {A[i]=0.;} 
  
    A[0]=1.; 
    A[4]=2.; 
    A[1]=2.;
    A[5]=5.; 
    A[15]=1.; 
    A[10]=0.2; 
    
    for(i=0;i<16;i++) B[i] = A[i]; 
  
    char jobvr = 'V';
    // char jobvl= 'N';
    integer n = 4;
   
    integer info;
   
    double * wr = new double[n]; 

    double* work = new double[16*n]; 
    integer lwork = 16*n; 

    char uplo = 'U' ; 

    dsyev_(&jobvr,&uplo , &n , A , &n, wr, work, &lwork, &info); 

    //  dggev_(&jobvl, &jobvr, &n, A, &n, B, &n, wr, wi, beta, &vl, &nvl, &vr, &nvl, 
    // work , &lwork, &info); 




    double *ev = new double[4];
    double *v2 = new double[4]; 
  
    /*cout << " Matrix B (original) " << endl; 
      for(i=0;i<n;i++) 
      {
      for(j=0;j<n;j++) 
      cout << B[i*n+j] << "\t" ; 
	  
      cout << endl; 
      }

      cout << " Matrix A (original) " << endl; 
      for(i=0;i<n;i++) 
      {
      for(j=0;j<n;j++) 
      cout << A[i*n+j] << "\t" ; 
	  
      cout << endl; 
      }

      for(i=0;i<n;i++)
      {
      for(j=0;j<n;j++) { ev[j] = A[i*n+j]; v2[j] = 0.; }
      cout << "lami" <<  (wr[i]) << endl;
      cout << " Residuum " << endl;  
      for(j=0;j<n;j++)
      {
      for(k=0;k<n;k++)
      {
      v2[j] += B[j*n+k]*ev[k];
      }
      cout <<  v2[j] - wr[i]*ev[j] << "\t" ; 
      }
      cout << endl; 
      } 
    */    
    delete[] A; 
    delete[] B; 
    delete[] wr; 
    delete[] work; 
    delete[] ev; 
    delete[] v2; 
  }
       
 


  //Attention: A,B are overwritten !! 
  void LaEigNSSolve(int hn, double * A, double * B, std::complex<double> * lami, int evecs_bool, double * evecs_re, double * evecs_im, char balance_type )
  {   
    integer n = hn;
    char jobvr , jobvl= 'N';
    bool balancing = 0; 

    if ( balance_type == 'B' || balance_type == 'S' || balance_type == 'P' )
      balancing =1; 

  
    integer info;
   
    double * wi= new double[n], * wr = new double[n]; 
    double * beta = new double[n]; 
    double vl=0; 
   
    integer nvl = 1; 
    integer nvr ; 

    if (evecs_bool)
      {
        jobvr = 'V'; 
        nvr = n; 
      }
    else 
      { 
        nvr=1; 
        jobvr = 'N';
      }

    double *vr = new double[nvr*nvr]; 
  
    double* work = new double[16*n]; 
    integer lwork = 16*n; 
    int i,j;
   
    char job=balance_type; // Permute and Scale in Balancing 
    integer ihi,ilo; 
    double * lscale, *rscale; 
    lscale = new double[n]; 
    rscale = new double[n]; 

    char side = 'R'; 
  

    if(balancing) 
      dggbal_(&job,&n, A,&n , B, &n, &ilo, &ihi,  lscale, rscale, work, &info) ; 
    else info =0; 

    if(info ==0 ) 
      { 
        dggev_(&jobvl, &jobvr, &n, A, &n, B, &n, wr, wi, beta, &vl, &nvl, vr, &nvr, 
               work , &lwork, &info);
      
        if(info==0) 	
          { 
            if(jobvr == 'V' && balancing) 
              {
                dggbak_(&job, &side, &n, &ilo, &ihi, lscale, rscale, &n, vr, &n,&info)  ; 
                if(info!=0)
                  {
                    cout << " Error in dggbak_ :: info  " << info << endl; 
                    return;
                  }
              }
          }
        else 
          {
            cout << " Error in dggev_ :: info  " << info << endl; 
            return;
          }  
      }
    else 
      {
        cout << " Error in dggbal_ :: info " << info << endl; 
        return; 
      }

    delete[] lscale; 
    delete[] rscale;  
  
    for(i=0;i<n;i++)
      {
        if (fabs(beta[i])>1e-30)  // not infinite eigenvalue 
          lami[i]=std::complex<double>(wr[i]/beta[i],wi[i]/beta[i]);
        else 
          {
            lami[i]=std::complex<double>(100.,100.); 
          }
      }
  
    if(evecs_bool)
      {
      
        for(i=0;i<n;i++)
          {
	 
            if( imag(lami[i])== 0. || beta[i] == 0.) //real eigenvalue -> real eigenvector in i-th line
              { 
                for(j=0;j<n;j++) 
                  {
                    // evecs[i*n+j]= std::complex<double>(vr[i*n + j],0.);
                    evecs_re[i*n+j]= vr[i*n+j];
                    evecs_im[i*n+j]= 0.; 
                  }
              } 
            else // conjugate complex eigenvectors 
              {
                for(j=0;j<n;j++)
                  {
                    // evecs[i*n+j]= std::complex<double>(vr[i*n+j],vr[(i+1)*n+j]);
                    // evecs[(i+1)*n+j]=std::complex<double>(vr[i*n+j],-vr[(i+1)*n+j]);
                    evecs_re[i*n+j]= vr[i*n+j];
                    evecs_re[(i+1)*n+j]=vr[i*n+j];
                    evecs_im[i*n+j]= vr[(i+1)*n+j];
                    evecs_im[(i+1)*n+j]=-vr[(i+1)*n+j];
                  }
                i++; 
              } 
          }
      }
   
 
    delete[] wi;
    delete[] wr; 
    delete[] beta; 
    delete[] work; 
    delete[] vr;


 
  }

  // Attention A,B are overwritten !!! 
  void LaEigNSSolve(int hn, std::complex<double> * A, std::complex<double> * B, std::complex<double> * lami, int evecs_bool, std::complex<double> * evecs, std::complex<double> * dummy, char balance_type)
  {
    integer n = hn;
    balance_type = 'N'; 
  
  
    std::complex<double> * at = new std::complex<double> [n*n];
    std::complex<double> * bt = new std::complex<double> [n*n];
  
    for(int i=0;i<n;i++)
      for(int   j=0;j<n;j++)
        at[j*n+i] = A[i*n+j];
      
    for(int i=0;i<n;i++)
      for(int   j=0;j<n;j++)
        bt[j*n+i] = B[i*n+j];
  
  
    char jobvr , jobvl= 'N';
    bool balancing = 0; 
  
    balance_type = 'P'; 
    if ( balance_type == 'B' || balance_type == 'S' || balance_type == 'P' )
      balancing =1; 
  
    balancing = 0; 
  
    std::complex<double> * alpha= new std::complex<double>[n];
    std::complex<double> * beta = new std::complex<double>[n]; 
    std::complex<double> vl=0.; 
  
    integer nvl = 1; 
    std::complex<double> * vr = 0;
  
    std::complex<double> * work = new std::complex<double>[8*n]; 
    integer lwork = 8*n; 
    double *rwork = new double[8*n];  
  
    integer nvr = n ; 
    if (evecs_bool) 
      {
        jobvr = 'V'; 
        vr = evecs; 
      }
    else jobvr = 'N'; 
  
    //std::complex<double>  * A1,*B1; 
    int i; 
  
    char job=balance_type; // Permute and Scale in Balancing 
    integer ihi,ilo; 
   
    double * lscale = new double[n]; 
    double * rscale = new double[n]; 
   
    double * work2 = new double[6*n];
  
    char side = 'R'; 
  
    integer info = 0;

    if(balancing) zggbal_(&job,&n, at, &n , bt, &n, &ilo, &ihi,  lscale, rscale, work2, &info) ; 
  
    if(info == 0 ) 
      {  
        zggev_(&jobvl, &jobvr, &n, at, &n, bt, &n, alpha, beta, &vl, &nvl, vr, &nvr,  
               work , &lwork, rwork,  &info);
      
        if(info==0) 	
          // if(info>=0) 	
          {
            if(jobvr == 'V' && balancing) 
              {
                zggbak_(&job, &side, &n, &ilo, &ihi, lscale, rscale, &n, vr, &n,&info)  ;
            
                if(info!=0)
                  { 
                    cout << "*****  Error in zggbak_ **** " << endl; 
                    return; 
                  }
              }
          } 
        else 
          {
            cout << "**** Error in zggev_, info = " << info << " *****" << endl; 
            return;
          }
      }	
    else 
      {
        cout << "**** Error in zggbal_ **** " << endl; 
        return;
      }
    
    delete[] work; 
    delete[] rwork;
  
    delete[] lscale; 
    delete[] rscale; 
    delete[] work2;   
 
  
 
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
    delete[] alpha; 
    delete[] beta; 
    delete[] at;
    delete[] bt; 
 
  }   


  // Attention A,B are overwritten !!! 
  void LaEigNSSolveX(int hn, std::complex<double> * A, std::complex<double> * B, std::complex<double> * lami, int evecs_bool, std::complex<double> * evecs, std::complex<double> * dummy, char balance_type)
  {
  
    integer n = hn;
    std::complex<double> * at = new std::complex<double> [n*n];
    std::complex<double> * bt = new std::complex<double> [n*n];
  
    for(int i=0;i<n;i++)
      for(int   j=0;j<n;j++)
        at[j*n+i] = A[i*n+j];
      
    for(int i=0;i<n;i++)
      for(int   j=0;j<n;j++)
        bt[j*n+i] = B[i*n+j];
  
  
  
  
  
  
    char jobvr , jobvl= 'N';
 
    if(balance_type != 'B' && balance_type != 'S' && balance_type != 'P') 
      balance_type = 'N'; 
  
    std::complex<double> * alpha= new std::complex<double>[n];
    std::complex<double> * beta = new std::complex<double>[n]; 
    std::complex<double> vl=0.; 
  
    integer nvl = 1; 
    std::complex<double> * vr = 0;
  
    std::complex<double> * work = new std::complex<double>[20*n]; 
    integer lwork = 20*n; 
    double *rwork = new double[20*n];  
  
    integer nvr = n ; 
    if (evecs_bool) 
      {
        jobvr = 'V'; 
        vr = evecs; 
      }
    else jobvr = 'N'; 
  
    //std::complex<double>  * A1,*B1; 
    int i; 
  
    // char job=balance_type; // Permute and Scale in Balancing 
    integer ihi,ilo; 
   
    double * lscale = new double[4*n]; 
    double * rscale = new double[4*n]; 
   
  
    integer * iwork = new integer[n+2];
  
  
    // char side = 'R'; 
    char sense = 'N'; // no recoprocal condition number computed !! 
    double rconde, rcondv; // not referenced if sense == N 
    logical bwork; // not referenced 
  
    integer info = 0;


    double abnrm, bbnrm; // 1-norm of A and B 
    
    zggevx_(&balance_type,&jobvl, &jobvr, &sense, &n, at, &n, bt, &n, alpha, beta, &vl, &nvl, vr, &nvr,  
            &ilo, &ihi, lscale, rscale, &abnrm, &bbnrm, &rconde, &rcondv, work , &lwork, rwork,  iwork, &bwork, &info);
    
    if(info!=0)       
      // if(info>=0)  
      {
        cout << " ****** INFO " << info << endl;  
        cout << "*****  Error in zggevx_ **** " << endl; 
        //throw (); 
        return; 
      }
    
    for(int i=0;i<n;i++) 
      cout << " i " << i << " alpha " << alpha[i] << " beta " << beta[i] << endl; 
        
        
        
    delete[] iwork;      
    delete[] work; 
    delete[] rwork;
  
    delete[] lscale; 
    delete[] rscale;    
 
  
 
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
    delete[] alpha; 
    delete[] beta; 
    delete[] at;
    delete[] bt; 
 
  }   

  void LaEigNSSolveX(int n, double * A, double * B, std::complex<double> * lami, int evecs_bool, double * evecs, double * dummy, char balance_type)
  {
    cout << "LaEigNSSolveX not implemented for double" << endl;
    // throw();  
  }

  void LapackSSEP(int hn, double* A, double* lami, double* evecs)  
  {
    integer n = hn;

    for(int i=0;i<n*n;i++) evecs[i] = A[i]; 
    
    char jobzm = 'V' , uplo = 'U'; 
  
    integer lwork=2*n*n; 
 
    double* work = new double[lwork];
  
    integer info; 
 
    dsyev_(&jobzm,&uplo , &n , evecs , &n, lami, work, &lwork, &info); 

    delete[] work; 
  }


  /*
  extern "C"
  void zhseqr_(char & job, char & compz, int & n, 
               int & ilo, int & ihi, complex<double> & h, int & ldh, 
               complex<double> & w, complex<double> & z, int & ldz, 
               complex<double> & work, int & lwork, int & info);

  extern "C"
  void zhsein_ (char & side, char & eigsrc, char & initv,
                int * select, 
		int & n, complex<double> & h, int & ldh, 
		complex<double> & w,
                complex<double> & vl, int & ldvl,
                complex<double> & vr, int & ldvr,
		int & mm, int & m, 
                complex<double> & work, double & rwork,
                int & ifaill, int & ifailr, int & info);
  */

  void LapackHessenbergEP (int hn, std::complex<double> * A, std::complex<double> * lami, std::complex<double> * evecs)
  {
    integer n = hn;
    integer lwork = 2 * n * n;  // or 6 n ?
    complex<double> * work = new complex<double>[lwork];

    complex<double> * hA = new complex<double>[n*n];
    for (int i = 0; i < n*n; i++)  { hA[i] = A[i]; }

    logical * select = new logical[n];
    for (int i = 0; i < n; i++) select[i] = 1; // 'V';


    complex<double> vl;
    //  complex<double> * vl = new complex<double>[n*n];
    //  complex<double> * vr = new complex<double>[n*n];

    integer info;

    // cout << "calls zhseqr" << endl;

    char job = 'E', compz = 'N';
    integer ilo = 1, ihi = n, ldh = n, ldz = n;
    zhseqr_(&job, &compz, &n, &ilo, &ihi, hA, &ldh, lami, evecs, &ldz, work, &lwork, &info);
    //  zhseqr_('S', 'I', n, 1, n, *A, n, *lami, *evecs, n, *work, lwork, info);


    if (info)
      cout << "error in eigensolver, info = " << info << endl;

    // for (int i = 0; i < n; i++)
    // cout << "ev(" << i << ") = " << lami[i] << endl;

    for (int i = 0; i < n*n; i++)  { hA[i] = A[i]; }
    double * rwork = new double[n];

    integer m = 0;
    char side = 'R', eigsrc = 'Q', initv = 'N';
    n = hn;
    ldh = n; 
    integer ldvl = n, ldvr = n, mm = n;
    integer * ifaill = new integer[n];
    integer * ifailr = new integer[n];

    for (int i = 0; i < n*n; i++)
      evecs[i] = -1.0;
    // cout << "call zhsein" << endl;

    zhsein_ (&side, &eigsrc, &initv, select, &n, A, &ldh, lami, &vl, &ldvl, evecs, &ldvr,
             &mm, &m, work, rwork, ifaill, ifailr, &info);

    if (info)
      cout << "error in eigensolver, info = " << info << endl;
    
    /*
    cout << "m  = " << m << endl;
    cout << "info = " << info << endl;
    cout << "m = " << m << endl;
    cout << "rwork[0] = " << rwork[0] << endl;
    cout << "evecs(0,0) = " << evecs[0] << endl;
    for (int i = 0; i < n; i++)
      cout << "ifail[" << i << "] = " << ifailr[i] << endl;
    // cout << "ifaill = " << ifaill << endl;
    //   cout << "ifailr = " << ifailr << endl;
    cout << "info = " << info << endl;
    */

    delete[] select;
    delete[] hA;
    delete[] rwork;
    delete[] work;
    // cout << "hessenberg complete" << endl;
  }



  void LapackGHEP(int hn, double* A, double* B,  double* lami)  
  {
    integer n = hn;
    double *B1 = new double[n*n]; 
    double *A1 = new double[n*n]; 
  
    for(int i=0;i<n*n;i++)
      { 
        A1[i] = A[i];
        B1[i] = B[i]; 
      } 
 
    char jobzm = 'V' , uplo = 'U'; 
  
    integer lwork=16*n; 
 
    double* work = new double[lwork];
  
    integer info; 
    integer itype =1; 


    dsygv_(&itype,&jobzm,&uplo , &n , A1 , &n, B1, &n, lami, work, &lwork, &info); 
  

    delete[] A1; 
    delete[] B1; 
    delete[] work; 
  }
       
  int LapackGHEPEPairs(int hn, double* A, double* B,  double* lami)  
  {
    char jobzm = 'V' , uplo = 'U'; 
    
    integer n = hn;
    integer lwork=4*n; 
 
    double* work = new double[lwork];
  
    integer info; 
    integer itype =1; 
    integer lda = n;
    integer ldb = n;

    dsygv_(&itype,&jobzm,&uplo , &n , A , &lda, B, &ldb, lami, work, &lwork, &info); 

    if(info != 0) 
      {
        cout << "LapackGHEPEPairs Info " << info << endl;  
        cout << "n = " << n << endl; 
      }
  
  
    delete [] work;  
    return(info); 
  }

  int LapackGHEPEPairs(int hn, complex<double>* A, complex<double>* B,  double* lami)  
  {
    integer n = hn;
    char jobzm = 'V' , uplo = 'U'; 
  
    integer lwork=8*n; 

    std::complex<double>* work = new std::complex<double>[lwork];
    double* rwork = new double[lwork]; 
  
    integer info; 
    integer itype =1; 
    integer lda = n;
    integer ldb = n;

    cout << " zhegv " << endl; 

    cout << " A s " << endl; 
    for(int j=0;j<n;j++)
      {
	for(int k=0;k<n;k++)	
	  cout <<  A[j*n+k] << " \t " ; 
	cout << endl; 
      }
    cout << " M " << endl; 
    for(int j=0;j<n;j++)
      {
	for(int k=0;k<n;k++)	
	  cout <<  B[j*n+k] << " \t " ; 
	cout << endl; 
      }


    zhegv_(&itype,&jobzm,&uplo , &n , A , &lda, B, &ldb, lami, work, &lwork, rwork, &info); 

    cout << " ... is back " << endl; 
    if(info != 0) 
      {
        cout << "LapackGHEPEPairs Info " << info << endl;  
        cout << "n = " << n << endl; 
      
      
      }
  
    delete [] work; 
    delete [] rwork;
  
    return(info); 
  }
       



  void LaLinearSolveComplex(int hn, std::complex<double> * A, std::complex<double> * F)
  {
    // Solve Ax=F
    // A on exit LU factorization 
    // F is overwritten by solution x 

    integer n = hn;
    integer nrhs =1; 
    integer *ipiv; 
    ipiv = new integer[n]; 
    integer info; 

  
    zgesv_(&n, &nrhs, A, &n, ipiv, F, &n, &info ); 

    if(info!=0) 
      cout << " ***** Error in LapackGEP.cpp LaLinearSolveComplex : info =  " <<  info << endl; 
    delete[] ipiv; 

    return; 
  } 


  void LaLinearSolve(int hn, double * A, double * F)
  {
    // Invert
    // A on exit LU factorization 
    // F is overwritten by solution x 

    integer n = hn;
    integer nrhs = n; 
    integer *ipiv; 
    ipiv = new integer[n*n]; 
    integer info; 

  
    dgesv_(&n, &nrhs, A, &n, ipiv, F, &n, &info ); 

    if(info!=0) 
      cout << " ***** Error in LapackGEP.cpp LaLinearSolveComplex : info =  " <<  info << endl; 
    delete[] ipiv; 

    return; 
  } 


  void LaLinearSolveRHS(int hn, double * A, double * F)
  {
    // Solve linear system A x = F for 1 rhs 
    // A on exit LU factorization 
    // F is overwritten by solution x 


    integer n = hn;
    integer nrhs = 1; 
    integer *ipiv; 
    ipiv = new integer[n]; 
    integer info; 

 

    dgesv_(&n, &nrhs, A, &n, ipiv, F, &n, &info ); 




    if(info!=0) 
      cout << " ***** Error in LapackGEP.cpp LaLinearSolveComplex : info =  " <<  info << endl; 
    delete[] ipiv; 

    return; 
  } 
}



#endif
