notwendig ?

// old style interface

#ifndef LAEVINTERFACE_HH
#define LAEVINTERFACE_HH

#include "arch.hpp" 


#ifdef OWN_LAPACK
extern "C"
{
// ********************* Eigen Solve Routines ***************************

  // Solve generalized double eigenvalue problem (QZ) 
  void dggev_(char *jobvl,char *jobvr,int *N,double *A,int *lda,double *B,
		       int *ldb,double *alphaR,double *alphaI,double* beta,
		       double* vl, int *ldvl,double *vr,int *ldvr,
		       double *work, int *lwork,int *info); 

  
  // Solve complex generalized eigenvalue problem (QZ) 
  void zggev_(char *jobvl,char* jobvr,int* N,complex<double>* A, int* lda,complex<double>* B, int* ldb, complex<double>* alpha, complex<double>* beta, complex<double>* vl, int* ldvl, complex<double>* vr, int* ldvr, complex<double>* work, int* lwork, double* rwork, int* info); 
  
  // Expert routine for solving generalized eigenvalue problem (QZ) 
  //  zggevx 
  //  dggevx 

  //Solve double standard symmetric eigenvlaue problem 
  void dsyev_(char *jobz, char *uplo, int* n , double *A , int* lda, double* w,
	      double* work, int* lwork, int* info); 


  //Balancing and permuting matrix pencil routines 
  void dggbal_(char *job, int* n, double* A, int* lda, double* B, int* ldb, 
	       int* ilo, int* ihi, double* lscale, double* rscale, 
	       double* work, int* info); 
  void dggbak_(char* job, char* side, int* n, int* ilo, int* ihi, double* lscale,
	       double* rscale, int* m, double* v, int* ldv, int* info); 
  void zggbal_(char *job, int* n, std::complex<double> *A, int* lda, 
	       std::complex<double> *B, int* ldb, int* ilo, int* ihi, 
	       double* lscale, double* rscale, double* work, int* info); 
  void zggbak_(char* job, char* side, int* n, int* ilo, int* ihi, double* lscale,
	       double* rscale, int* m, std::complex<double> *v, int* ldv, 
	       int* info);
  void dsygv_(int* itype, char *jobzm, char* uplo , int* n1, double* A , int * n, double* B, int *nb, double * lami, double * work, int * lwork, int * info); 
  void zhegv_(int* itype, char *jobzm, char* uplo , int* n1, std::complex<double>* A , int * n, std::complex<double>* B, int *nb, double * lami, std::complex<double> * work, int * lwork,double * rwork, int * info); 
  
  // void zsygv_(int* itype, char *jobzm, char* uplo , int* n1, std::complex<double>* A , int * n, std::complex<double>* B, int *nb, std::complex<double> * lami, std::complex<double> * work, int * lwork, int * info); 

  // General Linear Complex System Solver
  void zgesv_(int * n, int * nrhs, std::complex<double> * A, int * lda, int * ipiv, std::complex<double> * B, int * ldb, int * info );
 void dgesv_(int * n, int * nrhs, double * A, int * lda, int * ipiv, double * B, int * ldb, int * info );

 void zggevx_(char * balance, char * jobvl,  char * jobvr, char * sense, int * n, complex<double> * at,  int * n2, complex<double> * bt, int * n3, complex<double> * alpha, complex<double> * beta, complex<double> * vl, int * nvl, complex<double> * vr, int * nvr,  
          int * ilo, int * ihi, double * lscale, double *  rscale, double * abnrm, double * bbnrm, double * rconde, double * rcondv, complex<double> * work , int * lwork, double * rwork, int* iwork, bool * bwork, int * info);
 
 
}

#endif
#endif 
