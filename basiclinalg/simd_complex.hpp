#ifndef FILE_SIMD_COMPLEX
#define FILE_SIMD_COMPLEX

/**************************************************************************/
/* File:   simd.hpp                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   06. Nov. 16                                                    */
/**************************************************************************/

#include <core/simd.hpp>

namespace ngcore
{

  INLINE SIMD<mask64> Mask128 (int64_t nr)
  {
#ifdef __CUDA_ARCH__
    return 2*nr;

#else // __CUDA_ARCH__
  #ifdef __AVX512F__

    return _mm512_cmpgt_epi64_mask(_mm512_set1_epi64(nr),
                                   _mm512_set_epi64(3, 3, 2, 2, 1, 1, 0, 0));

  #elif defined(__AVX__)

    return my_mm256_cmpgt_epi64(_mm256_set1_epi64x(nr),
                                _mm256_set_epi64x(1, 1, 0, 0));

  #elif defined NETGEN_ARCH_AMD64
    return _mm_cmpgt_epi32(_mm_set1_epi32(nr),
                           _mm_set_epi32(0, 0, 0, 0));
  #else
    return 2*nr;
  #endif
#endif // __CUDA_ARCH__
  }


  
  template <int N>
  class SIMD<Complex, N>
  {
    SIMD<double,N> re, im;
  public:
    SIMD () = default;
    SIMD (SIMD<double,N> _r, SIMD<double,N> _i = 0.0) : re(_r), im(_i) { ; }
    SIMD (Complex c) : re(c.real()), im(c.imag()) { ; }
    SIMD (double d) : re(d), im(0.0) { ; }
    static constexpr int Size() { return SIMD<double>::Size(); }

    auto real() const { return re; }
    auto imag() const { return im; }
    auto & real() { return re; }
    auto & imag() { return im; }


    // Numbers in SIMD structure are not necessarily in same order as in memory
    // for instance:
    // [x0,y0,x1,y1,x2,y2,x3,y3] -> [x0,x2,x1,x3,y0,y2,y1,y3]
    void LoadFast (Complex * p)
    {
      SIMD<double,N> c1((double*)p);
      SIMD<double,N> c2((double*)(p+SIMD<double>::Size()/2));
      std::tie(re,im) = Unpack(c1,c2);
    }

    void StoreFast (Complex * p)
    {
      auto [h1,h2] = Unpack(re,im);
      h1.Store((double*)p);
      h2.Store((double*)(p+SIMD<double,N>::Size()/2));
    }

    void LoadFast (Complex * p, int nr)
    {
      SIMD<double> c1((double*)p, Mask128(nr));
      SIMD<double> c2((double*)(p+SIMD<double>::Size()/2), Mask128(nr-SIMD<double>::Size()/2));
      std::tie(re,im) = Unpack(c1,c2);
    }

    void StoreFast (Complex * p, int nr)
    {
      auto [h1,h2] = Unpack(re,im);
      h1.Store((double*)p, Mask128(nr));
      h2.Store((double*)(p+SIMD<double>::Size()/2), Mask128(nr-SIMD<double>::Size()/2));
    }
  };

  // templatize all with N ???
  template <int N>
  inline SIMD<double, N> Real(SIMD<Complex, N> a) { return a.real(); }
  template <int N>
  inline SIMD<double, N> Imag(SIMD<Complex, N> a) { return a.imag(); }

  template <int N> INLINE auto operator+ (SIMD<Complex, N> a, SIMD<Complex, N> b)
  { return SIMD<Complex, N> (a.real()+b.real(), a.imag()+b.imag()); }
  template <int N> INLINE auto operator+ (SIMD<Complex, N> a, SIMD<double, N> b)
  { return SIMD<Complex, N> (a.real()+b, a.imag()); }
  template <int N> INLINE auto operator+ (SIMD<double, N> b, SIMD<Complex, N> a)
  { return SIMD<Complex, N> (a.real()+b, a.imag()); }
  template <int N> INLINE auto operator+ (SIMD<Complex, N> a, double b)
  { return SIMD<Complex, N> (a.real()+b, a.imag()); }
  template <int N> INLINE auto operator+ (double b, SIMD<Complex, N> a)
  { return SIMD<Complex, N> (a.real()+b, a.imag()); }
  template <int N> INLINE auto operator+ (SIMD<Complex, N> a, Complex b)
  { return a + SIMD<Complex, N> (b); }
  template <int N> INLINE auto operator+ (Complex b, SIMD<Complex, N> a)
  { return a + SIMD<Complex, N> (b); }    

  
  template <int N> INLINE auto operator- (SIMD<Complex, N> a, SIMD<Complex, N> b)
  { return SIMD<Complex, N> (a.real()-b.real(), a.imag()-b.imag()); }
  

  template <int N>
  INLINE auto operator- (SIMD<Complex, N> a)
  { return SIMD<Complex, N>(-a.real(), -a.imag()); }

  
  template <int N> INLINE auto operator* (SIMD<Complex, N> a, SIMD<Complex, N> b)
  { return SIMD<Complex, N> (a.real()*b.real()-a.imag()*b.imag(),
                             a.real()*b.imag()+a.imag()*b.real()); }
  template <int N> INLINE auto operator* (SIMD<double, N> a, SIMD<Complex, N> b)
  { return SIMD<Complex, N> (a*b.real(), a*b.imag()); }
  template <int N> INLINE auto operator* (SIMD<Complex, N> b, SIMD<double, N> a)
  { return SIMD<Complex, N> (a*b.real(), a*b.imag()); }
  
  template <int N> INLINE auto operator* (double a, SIMD<Complex, N> b)
  { return SIMD<Complex, N> (a*b.real(), a*b.imag()); }
  template <int N> INLINE auto operator* (SIMD<Complex, N> b, double a)
  { return SIMD<Complex, N> (a*b.real(), a*b.imag()); }
  template <int N> INLINE auto operator* (Complex a, SIMD<Complex, N> b)
  { return SIMD<Complex, N> (a)*b; }
  template <int N> INLINE auto operator* (SIMD<Complex, N> b, Complex a)
  { return SIMD<Complex, N> (a)*b; }
  template <int N> INLINE auto operator* (SIMD<double, N> b, Complex a)
  { return SIMD<Complex, N> (a.real()*b, a.imag()*b); }
  template <int N> INLINE auto operator* (Complex a, SIMD<double, N> b)
  { return SIMD<Complex, N> (a.real()*b, a.imag()*b); }

  
  template <int N>
  INLINE SIMD<Complex, N> & operator+= (SIMD<Complex, N> & a, SIMD<Complex, N> b)
  { a.real()+=b.real(); a.imag()+=b.imag(); return a; }
  template <int N>
  INLINE SIMD<Complex, N> & operator-= (SIMD<Complex, N> & a, SIMD<Complex, N> b)
  { a.real()-=b.real(); a.imag()-=b.imag(); return a; }

  template <int N>
  INLINE SIMD<Complex, N> & operator*= (SIMD<Complex, N> & a, double b)
  { a.real()*= b; a.imag() *= b; return a; }  
  template <int N>
  INLINE SIMD<Complex, N> & operator*= (SIMD<Complex, N> & a, Complex b)
  { a = a*SIMD<Complex, N>(b); return a; }
  template <int N>
  INLINE SIMD<Complex, N> & operator*= (SIMD<Complex, N> & a, SIMD<double> b)
  { a.real()*= b; a.imag() *= b; return a; }
  template <int N>
  INLINE SIMD<Complex, N> & operator*= (SIMD<Complex, N> & a, SIMD<Complex, N> b)
  { a = a*b; return a; }

  template <int N>
  INLINE SIMD<Complex, N> Inv (SIMD<Complex, N> a)
  {
    SIMD<double> n2 = a.real()*a.real()+a.imag()*a.imag();
    return SIMD<Complex, N> (a.real()/n2, -a.imag()/n2);
  }
  
  template <int N> INLINE auto operator/ (double a, SIMD<Complex, N> b)
  { return a * Inv(b); }

  template <int N> INLINE auto operator/ (Complex a, SIMD<Complex, N> b)
  { return a * Inv(b); }

  template <int N> INLINE auto operator/ (SIMD<Complex, N> a, SIMD<Complex, N> b)
  { return a * Inv(b); }

  template <int N> INLINE auto operator/ (SIMD<Complex, N> a, SIMD<double, N> b)
  { return a * (1.0/b); }
  
  template <int N>
  INLINE Complex HSum (SIMD<Complex, N> sc)
  {
    auto [re,im] = HSum(sc.real(), sc.imag());
    return Complex(re,im);
  }

  template <int N>
  INLINE auto HSum (SIMD<Complex, N> sc1, SIMD<Complex, N> sc2)
  {
    auto [re1,im1,re2,im2] = HSum(sc1.real(), sc1.imag(), sc2.real(), sc2.imag());
    return std::tuple(Complex(re1,im1), Complex(re2,im2));
  }
  
  template <int N>
  ostream & operator<< (ostream & ost, SIMD<Complex, N> c)
  {
    ost << c.real() << ", " << c.imag();
    return ost;
  }


  template <typename FUNC, int N>
  SIMD<Complex, N> SIMDComplexWrapper (SIMD<Complex, N> x, FUNC f)
  {
    Complex hx[N];
    x.StoreFast(hx);
    for (auto & hxi : hx) hxi = f(hxi);
    x.LoadFast(hx);
    return x;
  }

  template <int N>
  inline SIMD<Complex, N> cos (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return cos(c); }); }
  
  template <int N>
  inline SIMD<Complex, N> sin (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return sin(c); }); }

  template <int N>
  inline SIMD<Complex, N> tan (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return tan(c); }); }

  template <int N>
  inline SIMD<Complex, N> atan (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return atan(c); }); }

  template <int N>
  inline SIMD<Complex, N> acos (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return acos(c); }); }

  template <int N>
  inline SIMD<Complex, N> asin (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return asin(c); }); }

  template <int N>
  inline SIMD<Complex, N> cosh (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return cosh(c); }); }
  
  template <int N>
  inline SIMD<Complex, N> sinh (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return sinh(c); }); }
  
  template <int N>
  inline SIMD<Complex, N> exp (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return exp(c); }); }

  template <int N>
  inline SIMD<Complex, N> log (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return log(c); }); }

  template <int N>
  inline SIMD<Complex, N> sqrt (SIMD<Complex, N> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return sqrt(c); }); }

  template <int N>
  inline SIMD<Complex, N> Conj (SIMD<Complex, N> x)
  { return SIMD<Complex, N> (x.real(), -x.imag()); } 

  template <int N>
  INLINE SIMD<Complex, N> IfPos (SIMD<Complex, N> a, SIMD<Complex, N> b, SIMD<Complex, N> c)
  {
    return SIMD<Complex, N> (IfPos (a.real(), b.real(), c.real()),
                             IfPos (a.real(), b.imag(), c.imag()));
  }
}


namespace ngbla
{
  template <typename T> struct is_scalar_type;
  template <int N>
  struct is_scalar_type<ngcore::SIMD<double,N>> { static constexpr bool value = true; };

  template <typename T> struct is_scalar_type;
  template <int N>
  struct is_scalar_type<ngcore::SIMD<ngcore::Complex,N>> { static constexpr bool value = true; };
}


#endif
