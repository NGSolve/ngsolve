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
  using std::tie;

  typedef std::complex<double> Complex;

  INLINE SIMD<mask64> Mask128 (int64_t nr)
  {
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
  }


  
  template <>
  class SIMD<Complex>
  {
    SIMD<double> re, im;
  public:
    SIMD () = default;
    SIMD (SIMD<double> _r, SIMD<double> _i = 0.0) : re(_r), im(_i) { ; }
    SIMD (Complex c) : re(c.real()), im(c.imag()) { ; }
    SIMD (double d) : re(d), im(0.0) { ; }
    static constexpr int Size() { return SIMD<double>::Size(); }

    SIMD<double> real() const { return re; }
    SIMD<double> imag() const { return im; }
    SIMD<double> & real() { return re; }
    SIMD<double> & imag() { return im; }


    // Numbers in SIMD structure are not necessarily in same order as in memory
    // for instance:
    // [x0,y0,x1,y1,x2,y2,x3,y3] -> [x0,x2,x1,x3,y0,y2,y1,y3]
    void LoadFast (Complex * p)
    {
      SIMD<double> c1((double*)p);
      SIMD<double> c2((double*)(p+SIMD<double>::Size()/2));
      tie(re,im) = Unpack(c1,c2);
    }

    void StoreFast (Complex * p)
    {
      SIMD<double> h1, h2;
      tie(h1,h2) = Unpack(re,im);
      h1.Store((double*)p);
      h2.Store((double*)(p+SIMD<double>::Size()/2));
    }

    void LoadFast (Complex * p, int nr)
    {
      SIMD<double> c1((double*)p, Mask128(nr));
      SIMD<double> c2((double*)(p+SIMD<double>::Size()/2), Mask128(nr-SIMD<double>::Size()/2));
      tie(re,im) = Unpack(c1,c2);
    }

    void StoreFast (Complex * p, int nr)
    {
      SIMD<double> h1, h2;
      tie(h1,h2) = Unpack(re,im);
      h1.Store((double*)p, Mask128(nr));
      h2.Store((double*)(p+SIMD<double>::Size()/2), Mask128(nr-SIMD<double>::Size()/2));
    }
  };

  inline SIMD<double> Real(SIMD<Complex> a) { return a.real(); }
  inline SIMD<double> Imag(SIMD<Complex> a) { return a.imag(); }

  INLINE SIMD<Complex> operator+ (SIMD<Complex> a, SIMD<Complex> b)
  { return SIMD<Complex> (a.real()+b.real(), a.imag()+b.imag()); }
  INLINE SIMD<Complex> & operator+= (SIMD<Complex> & a, SIMD<Complex> b)
  { a.real()+=b.real(); a.imag()+=b.imag(); return a; }
  INLINE SIMD<Complex> operator- (SIMD<Complex> a, SIMD<Complex> b)
  { return SIMD<Complex> (a.real()-b.real(), a.imag()-b.imag()); }
  INLINE SIMD<Complex> operator- (SIMD<Complex> a)
  { return SIMD<Complex> (-a.real(), -a.imag()); }
  INLINE SIMD<Complex> operator* (SIMD<Complex> a, SIMD<Complex> b)
    { return SIMD<Complex> (a.real()*b.real()-a.imag()*b.imag(),
                            a.real()*b.imag()+a.imag()*b.real()); }
  INLINE SIMD<Complex> operator* (SIMD<double> a, SIMD<Complex> b)
  { return SIMD<Complex> (a*b.real(), a*b.imag()); }
  INLINE SIMD<Complex> operator* (SIMD<Complex> b, SIMD<double> a)
  { return SIMD<Complex> (a*b.real(), a*b.imag()); }
  INLINE SIMD<Complex> operator* (double a, SIMD<Complex> b)
  { return SIMD<Complex> (a*b.real(), a*b.imag()); }
  INLINE SIMD<Complex> & operator*= (SIMD<Complex> & a, double b)
  { a.real()*= b; a.imag() *= b; return a; }  
  INLINE SIMD<Complex> & operator*= (SIMD<Complex> & a, Complex b)
  { a = a*SIMD<Complex>(b); return a; }
  INLINE SIMD<Complex> & operator*= (SIMD<Complex> & a, SIMD<double> b)
  { a.real()*= b; a.imag() *= b; return a; }
  INLINE SIMD<Complex> & operator*= (SIMD<Complex> & a, SIMD<Complex> b)
  { a = a*b; return a; }

  INLINE SIMD<Complex> Inv (SIMD<Complex> a)
  {
    SIMD<double> n2 = a.real()*a.real()+a.imag()*a.imag();
    return SIMD<Complex> (a.real()/n2, -a.imag()/n2);
  }
  
  INLINE SIMD<Complex> operator/ (SIMD<Complex> a, SIMD<Complex> b)
  {
    return a * Inv(b);
  }

  
  INLINE Complex HSum (SIMD<Complex> sc)
  {
    double re, im;
    std::tie(re,im) = HSum(sc.real(), sc.imag());
    return Complex(re,im);
  }

  INLINE auto HSum (SIMD<Complex> sc1, SIMD<Complex> sc2)
  {
    // double re1, im1, re2, im2;
    // std::tie(re1,im1,re2,im2) = HSum(sc1.real(), sc1.imag(), sc2.real(), sc2.imag());
    auto [re1,im1,re2,im2] = HSum(sc1.real(), sc1.imag(), sc2.real(), sc2.imag());
    return make_tuple(Complex(re1,im1), Complex(re2,im2));
  }
  
  INLINE ostream & operator<< (ostream & ost, SIMD<Complex> c)
  {
    ost << c.real() << ", " << c.imag();
    return ost;
  }


  template <typename FUNC>
  SIMD<Complex> SIMDComplexWrapper (SIMD<Complex> x, FUNC f)
  {
    Complex hx[SIMD<double>::Size()];
    x.StoreFast(hx);
    for (auto & hxi : hx) hxi = f(hxi);
    x.LoadFast(hx);
    return x;
  }

  inline SIMD<Complex> cos (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return cos(c); }); }
  
  inline SIMD<Complex> sin (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return sin(c); }); }

  inline SIMD<Complex> tan (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return tan(c); }); }

  inline SIMD<Complex> atan (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return atan(c); }); }

  inline SIMD<Complex> acos (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return acos(c); }); }

  inline SIMD<Complex> asin (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return asin(c); }); }

  inline SIMD<Complex> cosh (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return cosh(c); }); }
  
  inline SIMD<Complex> sinh (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return sinh(c); }); }
  
  inline SIMD<Complex> exp (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return exp(c); }); }

  inline SIMD<Complex> log (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return log(c); }); }

  inline SIMD<Complex> sqrt (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return sqrt(c); }); }

  inline SIMD<Complex> Conj (SIMD<Complex> x)
  { return SIMD<Complex> (x.real(), -x.imag()); } 

  INLINE Complex IfPos (Complex a, Complex b, Complex c)
  {
    return Complex (IfPos (a.real(), b.real(), c.real()),
                    IfPos (a.real(), b.imag(), c.imag()));
  }

  INLINE SIMD<Complex> IfPos (SIMD<Complex> a, SIMD<Complex> b, SIMD<Complex> c)
  {
    return SIMD<Complex> (IfPos (a.real(), b.real(), c.real()),
                          IfPos (a.real(), b.imag(), c.imag()));
  }
}



#endif
