#ifndef FILE_NGBLA
#define FILE_NGBLA


#ifdef USE_BLAS

extern "C"
{
  // #include <cblas.h>
#include <mkl_cblas.h>
}
#endif


#include <ngstd.hpp>


/// namespace for basic linear algebra
namespace ngbla
{
  using namespace std;
  using namespace ngstd;

  using ngstd::CArray;

  typedef std::complex<double> Complex;
  inline double fabs (Complex v) { return std::abs (v); }

  inline bool IsComplex(double v) { return false; }
  inline bool IsComplex(Complex v) { return true; }

#include "expr.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "cholesky.hpp"
#include "symmetricmatrix.hpp"
#include "bandmatrix.hpp"


  /// Computes eigenvalues and vectors of the symmetric matrix mat.
  extern void CalcEigenSystem (const FlatMatrix<double> & mat,
			       FlatVector<double> & lami,
			       FlatMatrix<double> & eigenvecs);



  /// Computes the Schur Complement.
  extern void CalcSchurComplement (const FlatMatrix<double> a, 
				   FlatMatrix<double> s,
				   const BitArray & used,
				   LocalHeap & lh);
}

#include "ng_lapack.hpp"

#endif
