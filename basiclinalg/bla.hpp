#ifndef FILE_NGBLA
#define FILE_NGBLA

// #ifdef USE_BLAS
// extern "C"
// {
// #include <mkl_cblas.h>
// }
// #endif


#include <ngstd.hpp>


/// namespace for basic linear algebra
namespace ngbla
{
  using namespace std;
  using namespace ngstd;

  using ngstd::CArray;

#ifdef USE_MYCOMPLEX
  typedef MyComplex<double> Complex;
  using std::fabs;
  inline double fabs (Complex v) { return ngstd::abs (v); }
#else
  // typedef std::complex<double> Complex;
  using std::fabs;
  inline double fabs (Complex v) { return std::abs (v); }
#endif

  using ngcore::AtomicAdd;
  inline void AtomicAdd (Complex & x, Complex y)
  {
    auto real = y.real();
    ngstd::AtomicAdd (reinterpret_cast<double(&)[2]>(x)[0], real);
    auto imag = y.imag();
    ngstd::AtomicAdd (reinterpret_cast<double(&)[2]>(x)[1], imag);
  }

  inline bool IsComplex(double v) { return false; }
  inline bool IsComplex(Complex v) { return true; }
}


namespace ngstd
{
  inline Archive & operator& (Archive & ar, Complex & c)
  {
    double hr, hi;
    if (ar.Output())
      { hr = c.real(); hi = c.imag(); }
    ar & hr & hi;
    if (ar.Input())
      c = Complex(hr, hi);
    return ar;
  }
}


#ifdef PARALLEL
namespace ngcore
{
  template <> struct MPI_typetrait<ngbla::Complex> {
    static MPI_Datatype MPIType ()  { return MPI_CXX_DOUBLE_COMPLEX; }
      // return MPI_C_DOUBLE_COMPLEX;   // no MPI_SUM defined ??
      // return MPI_DOUBLE_COMPLEX;
  };
}
#endif




#include "expr.hpp"
#include "vector.hpp"
#include "matrix.hpp"
#include "avector.hpp"
#include "ngblas.hpp"
#include "cholesky.hpp"
#include "symmetricmatrix.hpp"
#include "bandmatrix.hpp"
#include "triangular.hpp"
#include "householder.hpp"
#include "tensor.hpp"

namespace ngbla
{

  /// Computes eigenvalues and vectors of the symmetric matrix mat.
  extern NGS_DLL_HEADER void CalcEigenSystem (const FlatMatrix<double> & mat,
			       FlatVector<double> & lami,
			       FlatMatrix<double> & eigenvecs);



  /// Computes the Schur Complement.
  extern NGS_DLL_HEADER void CalcSchurComplement (const FlatMatrix<double> a, 
				   FlatMatrix<double> s,
				   const BitArray & used,
				   LocalHeap & lh);
}

#include "ng_lapack.hpp"

#endif
