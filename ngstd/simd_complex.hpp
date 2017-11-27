#ifndef FILE_SIMD_COMPLEX
#define FILE_SIMD_COMPLEX

/**************************************************************************/
/* File:   simd.hpp                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   06. Nov. 16                                                    */
/**************************************************************************/

namespace ngstd
{
  typedef std::complex<double> Complex;
    
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

#ifdef __AVX512F__
    void Load (Complex * p)
    {
      __m512d c1 = _mm512_loadu_pd((double*)p);
      __m512d c2 = _mm512_loadu_pd((double*)(p+4));
      re = _mm512_unpacklo_pd(c1,c2);
      im = _mm512_unpackhi_pd(c1,c2);
    }
    void Store (Complex * p) const
    {
      _mm512_storeu_pd ((double*)p, _mm512_unpacklo_pd(re.Data(),im.Data()));
      _mm512_storeu_pd ((double*)(p+4), _mm512_unpackhi_pd(re.Data(),im.Data()));
    }
    void Load (Complex * p, size_t mask)
    {
      __mmask8 mask1 = _mm512_cmpgt_epi64_mask(_mm512_set1_epi64(mask&3),
                                                 _mm512_set_epi64(3, 3, 2, 2, 1, 1, 0, 0));
      __mmask8 mask2 = _mm512_cmpgt_epi64_mask(_mm512_set1_epi64(mask&3),
                                                 _mm512_set_epi64(7, 7, 6, 6, 5, 5, 4, 4));

      __m512d c1 = _mm512_mask_load_pd(_mm512_setzero_pd(), mask1, (double*)p);
      __m512d c2 = _mm512_mask_load_pd(_mm512_setzero_pd(), mask2, (double*)(p+4));
      re = _mm512_unpacklo_pd(c1,c2);
      im = _mm512_unpackhi_pd(c1,c2);
    }
    void Store (Complex * p, size_t mask) const
    {
      __mmask8 mask1 = _mm512_cmpgt_epi64_mask(_mm512_set1_epi64(mask&3),
                                                 _mm512_set_epi64(3, 3, 2, 2, 1, 1, 0, 0));
      __mmask8 mask2 = _mm512_cmpgt_epi64_mask(_mm512_set1_epi64(mask&3),
                                                 _mm512_set_epi64(7, 7, 6, 6, 5, 5, 4, 4));
      _mm512_mask_store_pd ((double*)p, mask1,  _mm512_unpacklo_pd(re.Data(), im.Data()));
      _mm512_mask_store_pd ((double*)(p+4), mask2, _mm512_unpackhi_pd(re.Data(), im.Data()));
    }

#elif  defined (__AVX__)
    void Load (Complex * p)
    {
      __m256d c1 = _mm256_loadu_pd((double*)p);
      __m256d c2 = _mm256_loadu_pd((double*)(p+2));
      re = _mm256_unpacklo_pd(c1,c2);
      im = _mm256_unpackhi_pd(c1,c2);
    }
    void Store (Complex * p) const
    {
      _mm256_storeu_pd ((double*)p, _mm256_unpacklo_pd(re.Data(),im.Data()));
      _mm256_storeu_pd ((double*)(p+2), _mm256_unpackhi_pd(re.Data(),im.Data()));
    }

    void Load (Complex * p, size_t mask)
    {
      __m256i mask1 = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(mask&3),
                                           _mm256_set_epi64x(1, 1, 0, 0));
      __m256i mask2 = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(mask&3),

                                           _mm256_set_epi64x(3, 3, 2, 2));
      
      __m256d c1 = _mm256_maskload_pd((double*)p, mask1);
      __m256d c2 = _mm256_maskload_pd((double*)(p+2), mask2);
      re = _mm256_unpacklo_pd(c1,c2);
      im = _mm256_unpackhi_pd(c1,c2);
    }
    void Store (Complex * p, size_t mask) const
    {
      __m256i mask1 = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(mask&3),
                                           _mm256_set_epi64x(1, 1, 0, 0));
      __m256i mask2 = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(mask&3),
                                           _mm256_set_epi64x(3, 3, 2, 2));

      _mm256_maskstore_pd ((double*)p, mask1, _mm256_unpacklo_pd(re.Data(),im.Data()));
      _mm256_maskstore_pd ((double*)(p+2), mask2, _mm256_unpackhi_pd(re.Data(),im.Data()));
    }
#else
    void Load (Complex * p)
    {
      Complex c = *p;
      re = c.real();
      im = c.imag();
    }
    void Store (Complex * p) const
    {
      *p = Complex(re.Data(), im.Data());
    }
#endif
  };
    
  INLINE SIMD<Complex> operator+ (SIMD<Complex> a, SIMD<Complex> b)
  { return SIMD<Complex> (a.real()+b.real(), a.imag()+b.imag()); }
  INLINE SIMD<Complex> & operator+= (SIMD<Complex> & a, SIMD<Complex> b)
  { a.real()+=b.real(); a.imag()+=b.imag(); return a; }
  INLINE SIMD<Complex> operator- (SIMD<Complex> a, SIMD<Complex> b)
  { return SIMD<Complex> (a.real()-b.real(), a.imag()-b.imag()); }
  INLINE SIMD<Complex> operator* (SIMD<Complex> a, SIMD<Complex> b)
    { return SIMD<Complex> (a.real()*b.real()-a.imag()*b.imag(),
                            a.real()*b.imag()+a.imag()*b.real()); }
  INLINE SIMD<Complex> operator* (SIMD<double> a, SIMD<Complex> b)
  { return SIMD<Complex> (a*b.real(), a*b.imag()); }
  INLINE SIMD<Complex> operator* (SIMD<Complex> b, SIMD<double> a)
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
    double re1, im1, re2, im2;
    std::tie(re1,im1,re2,im2) = HSum(sc1.real(), sc1.imag(), sc2.real(), sc2.imag());
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
    x.Store(hx);
    for (auto & hxi : hx) hxi = f(hxi);
    x.Load(hx);
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
  
  inline SIMD<Complex> exp (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return exp(c); }); }

  inline SIMD<Complex> log (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return log(c); }); }

  inline SIMD<Complex> sqrt (SIMD<Complex> x)
  { return SIMDComplexWrapper (x, [](Complex c) { return sqrt(c); }); }

  inline SIMD<Complex> Conj (SIMD<Complex> x)
  { return SIMD<Complex> (x.real(), -x.imag()); } 

}



#endif
