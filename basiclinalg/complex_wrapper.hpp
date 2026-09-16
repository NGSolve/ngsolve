#ifndef COMPLEX_WRAPPER
#define COMPLEX_WRAPPER



#include <complex>





namespace ngcore
{
  using Complex32 = std::complex<float>;
  
  typedef std::complex<double> Complex;
  using std::fabs;
  inline double fabs (Complex v) { return std::abs (v); }
  
}


/// namespace for basic linear algebra
namespace ngbla
{
  using ngcore::Complex;
  using ngcore::Complex32;

  using ngcore::AtomicAdd;
  inline void AtomicAdd (Complex32 & x, Complex32 y)
  {
    auto real = y.real();
    auto imag = y.imag();
    ngcore::AtomicAdd (reinterpret_cast<float(&)[2]>(x)[0], real);
    ngcore::AtomicAdd (reinterpret_cast<float(&)[2]>(x)[1], imag);
  }

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
