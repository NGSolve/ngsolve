#ifndef COMPLEX_WRAPPER
#define COMPLEX_WRAPPER

#include <complex>

#ifdef USE_MYCOMPLEX
#include <mycomplex.hpp>
#endif




// #ifdef __clang__
#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
namespace std
{
  // avoid expensive call to complex mult by using the grammar school implementation
  INLINE std::complex<double> operator* (std::complex<double> a, std::complex<double> b)
  {
    return std::complex<double> (a.real()*b.real()-a.imag()*b.imag(),
                                 a.real()*b.imag()+a.imag()*b.real());
  }
}
#endif



namespace ngcore
{
#ifdef USE_MYCOMPLEX
  typedef ngstd::MyComplex<double> Complex;
  using std::fabs;
  inline double fabs (Complex v) { return ngstd::abs (v); }
#else
  typedef std::complex<double> Complex;
  using std::fabs;
  inline double fabs (Complex v) { return std::abs (v); }
#endif
}



/// namespace for basic linear algebra
namespace ngbla
{
  using ngcore::Complex;

  using ngcore::AtomicAdd;
  inline void AtomicAdd (Complex & x, Complex y)
  {
    auto real = y.real();
    ngcore::AtomicAdd (reinterpret_cast<double(&)[2]>(x)[0], real);
    auto imag = y.imag();
    ngcore::AtomicAdd (reinterpret_cast<double(&)[2]>(x)[1], imag);
  }

  inline bool IsComplex(double v) { return false; }
  inline bool IsComplex(Complex v) { return true; }
}


namespace ngstd
{
  using ngcore::Complex;
  INLINE Complex IfPos (Complex a, Complex b, Complex c)
  {
    return a.real() > 0 ? b : c;
  }
}



#endif
