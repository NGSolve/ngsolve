#ifndef NG_SIMD_AVX_HPP
#define NG_SIMD_AVX_HPP

/**************************************************************************/
/* File:   simd_avx.hpp                                                       */
/* Author: Joachim Schoeberl, Matthias Hochsteger                         */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

#include <immintrin.h>

namespace ngstd
{
#if defined(__AVX2__)
  INLINE __m256i my_mm256_cmpgt_epi64 (__m256i a, __m256i b)
  {
    return _mm256_cmpgt_epi64 (a,b);
  }
#else
  INLINE __m256i my_mm256_cmpgt_epi64 (__m256i a, __m256i b)
  {
    __m128i rlo = _mm_cmpgt_epi64(_mm256_extractf128_si256(a, 0),
                                  _mm256_extractf128_si256(b, 0));
    __m128i rhi = _mm_cmpgt_epi64(_mm256_extractf128_si256(a, 1),
                                  _mm256_extractf128_si256(b, 1));
    return _mm256_insertf128_si256 (_mm256_castsi128_si256(rlo), rhi, 1);
  }
#endif


  template <> 
  class SIMD<mask64,4>
  {
    __m256i mask;
  public:
    SIMD (size_t i)
      : mask(my_mm256_cmpgt_epi64(_mm256_set1_epi64x(i),
                                  _mm256_set_epi64x(3, 2, 1, 0)))
    { ; }
    SIMD (__m256i _mask) : mask(_mask) { ; }    
    SIMD (__m256d _mask) : mask(_mm256_castpd_si256(_mask)) { ; }    
    __m256i Data() const { return mask; }
    static constexpr int Size() { return 4; }    
    static SIMD<mask64, 4> GetMaskFromBits (unsigned int i);
  };

  static SIMD<mask64, 4> masks_from_4bits[16] = {
    _mm256_set_epi64x (0,0,0,0), _mm256_set_epi64x (0,0,0,-1),
    _mm256_set_epi64x (0,0,-1,0), _mm256_set_epi64x (0,0,-1,-1),
    _mm256_set_epi64x (0,-1,0,0), _mm256_set_epi64x (0,-1,0,-1),
    _mm256_set_epi64x (0,-1,-1,0), _mm256_set_epi64x (0,-1,-1,-1),
    _mm256_set_epi64x (-1,0,0,0), _mm256_set_epi64x (-1,0,0,-1),
    _mm256_set_epi64x (-1,0,-1,0), _mm256_set_epi64x (-1,0,-1,-1),
    _mm256_set_epi64x (-1,-1,0,0), _mm256_set_epi64x (-1,-1,0,-1),
    _mm256_set_epi64x (-1,-1,-1,0), _mm256_set_epi64x (-1,-1,-1,-1)
  };

  INLINE SIMD<mask64, 4> SIMD<mask64, 4> :: GetMaskFromBits (unsigned int i)
  {
    return masks_from_4bits[i & 15];
  }

  template<>
  class SIMD<int64_t,4> 
  {
    __m256i data;
    
  public:
    static constexpr int Size() { return 4; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (int64_t val) { data = _mm256_set1_epi64x(val); }
    SIMD (int64_t v0, int64_t v1, int64_t v2, int64_t v3) { data = _mm256_set_epi64x(v3,v2,v1,v0); }
    // SIMD (SIMD<double,2> v0, SIMD<double,2> v1) : SIMD(v0[0], v0[1], v1[0], v1[1]) { ; }
    SIMD (__m256i _data) { data = _data; }

    INLINE auto operator[] (int i) const { return ((int64_t*)(&data))[i]; }
    INLINE __m256i Data() const { return data; }
    INLINE __m256i & Data() { return data; }

    SIMD<int64_t,2> Lo() const { return _mm256_extractf128_si256(data, 0); }
    SIMD<int64_t,2> Hi() const { return _mm256_extractf128_si256(data, 1); }
    static SIMD FirstInt(int n0=0) { return { n0+0, n0+1, n0+2, n0+3 }; }
  };


  INLINE SIMD<int64_t,4> operator-(SIMD<int64_t,4> a) { return _mm256_sub_epi64(_mm256_setzero_si256(), a.Data()); }

#ifdef __AVX2__
  INLINE SIMD<int64_t,4> operator+ (SIMD<int64_t,4> a, SIMD<int64_t,4> b) { return _mm256_add_epi64(a.Data(),b.Data()); }
  INLINE SIMD<int64_t,4> operator- (SIMD<int64_t,4> a, SIMD<int64_t,4> b) { return _mm256_sub_epi64(a.Data(),b.Data()); }
#else
  INLINE SIMD<int64_t,4> operator+ (SIMD<int64_t,4> a, SIMD<int64_t,4> b) {
      auto lo_sum = _mm256_extractf128_si256(a.Data(), 0) + _mm256_extractf128_si256(b.Data(), 0);
      auto hi_sum = _mm256_extractf128_si256(a.Data(), 1) + _mm256_extractf128_si256(b.Data(), 1);
      return _mm256_set_m128i(hi_sum,lo_sum);
  }
  INLINE SIMD<int64_t,4> operator- (SIMD<int64_t,4> a, SIMD<int64_t,4> b) {
      auto lo_sub = _mm256_extractf128_si256(a.Data(), 0) - _mm256_extractf128_si256(b.Data(), 0);
      auto hi_sub = _mm256_extractf128_si256(a.Data(), 1) - _mm256_extractf128_si256(b.Data(), 1);
      return _mm256_set_m128i(hi_sub,lo_sub);
  }
#endif

  template<>
  class SIMD<double,4>
  {
    __m256d data;
    
  public:
    static constexpr int Size() { return 4; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (double val) { data = _mm256_set1_pd(val); }
    SIMD (int val)    { data = _mm256_set1_pd(val); }
    SIMD (size_t val) { data = _mm256_set1_pd(val); }
    SIMD (double v0, double v1, double v2, double v3) { data = _mm256_set_pd(v3,v2,v1,v0); }
    SIMD (SIMD<double,2> v0, SIMD<double,2> v1) : SIMD(v0[0], v0[1], v1[0], v1[1]) { ; }
    SIMD (double const * p) { data = _mm256_loadu_pd(p); }
    SIMD (double const * p, SIMD<mask64,4> mask) { data = _mm256_maskload_pd(p, mask.Data()); }
    SIMD (__m256d _data) { data = _data; }

    void Store (double * p) { _mm256_storeu_pd(p, data); }
    void Store (double * p, SIMD<mask64,4> mask) { _mm256_maskstore_pd(p, mask.Data(), data); }    
    
    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>                                                                    SIMD (const T & func)
    {   
      data = _mm256_set_pd(func(3), func(2), func(1), func(0));              
    }   
    
    INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    INLINE double & operator[] (int i) { return ((double*)(&data))[i]; }
    // [[deprecated("don't write to individual elments of SIMD")]]                
    // INLINE double & operator[] (int i) { return ((double*)(&data))[i]; }
    template <int I>
    double Get() const { return ((double*)(&data))[I]; }
    INLINE __m256d Data() const { return data; }
    INLINE __m256d & Data() { return data; }

    SIMD<double,2> Lo() const { return _mm256_extractf128_pd(data, 0); }
    SIMD<double,2> Hi() const { return _mm256_extractf128_pd(data, 1); }

    operator tuple<double&,double&,double&,double&> ()
    { return tuple<double&,double&,double&,double&>((*this)[0], (*this)[1], (*this)[2], (*this)[3]); }
  };

  INLINE auto Unpack (SIMD<double,4> a, SIMD<double,4> b)
  {
    return make_tuple(SIMD<double,4>(_mm256_unpacklo_pd(a.Data(),b.Data())),
                      SIMD<double,4>(_mm256_unpackhi_pd(a.Data(),b.Data())));
  }
  
  INLINE SIMD<double,4> operator- (SIMD<double,4> a) { return _mm256_xor_pd(a.Data(), _mm256_set1_pd(-0.0)); }
  INLINE SIMD<double,4> operator+ (SIMD<double,4> a, SIMD<double,4> b) { return _mm256_add_pd(a.Data(),b.Data()); }
  INLINE SIMD<double,4> operator- (SIMD<double,4> a, SIMD<double,4> b) { return _mm256_sub_pd(a.Data(),b.Data()); }
  INLINE SIMD<double,4> operator* (SIMD<double,4> a, SIMD<double,4> b) { return _mm256_mul_pd(a.Data(),b.Data()); }
  INLINE SIMD<double,4> operator/ (SIMD<double,4> a, SIMD<double,4> b) { return _mm256_div_pd(a.Data(),b.Data()); }
  INLINE SIMD<double,4> operator* (double a, SIMD<double,4> b) { return _mm256_set1_pd(a)*b.Data(); }
  INLINE SIMD<double,4> operator* (SIMD<double,4> b, double a) { return _mm256_set1_pd(a)*b.Data(); }

  INLINE SIMD<double,4> sqrt (SIMD<double,4> a) { return _mm256_sqrt_pd(a.Data()); }
  INLINE SIMD<double,4> floor (SIMD<double,4> a) { return _mm256_floor_pd(a.Data()); }
  INLINE SIMD<double,4> ceil (SIMD<double,4> a) { return _mm256_ceil_pd(a.Data()); }
  INLINE SIMD<double,4> fabs (SIMD<double,4> a) { return _mm256_max_pd(a.Data(), -a.Data()); }

  INLINE SIMD<mask64,4> operator<= (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_LE_OQ); }
  INLINE SIMD<mask64,4> operator< (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_LT_OQ); }
  INLINE SIMD<mask64,4> operator>= (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_GE_OQ); }
  INLINE SIMD<mask64,4> operator> (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_GT_OQ); }
  INLINE SIMD<mask64,4> operator== (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_EQ_OQ); }
  INLINE SIMD<mask64,4> operator!= (SIMD<double,4> a , SIMD<double,4> b)
  { return _mm256_cmp_pd (a.Data(), b.Data(), _CMP_NEQ_OQ); }

  INLINE SIMD<mask64,4> operator<= (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  _mm256_xor_si256(_mm256_cmpgt_epi64(a.Data(),b.Data()),_mm256_set1_epi32(-1)); }
  INLINE SIMD<mask64,4> operator< (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  my_mm256_cmpgt_epi64(b.Data(),a.Data()); }
  INLINE SIMD<mask64,4> operator>= (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  _mm256_xor_si256(_mm256_cmpgt_epi64(b.Data(),a.Data()),_mm256_set1_epi32(-1)); }
  INLINE SIMD<mask64,4> operator> (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  my_mm256_cmpgt_epi64(a.Data(),b.Data()); }
  INLINE SIMD<mask64,4> operator== (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  _mm256_cmpeq_epi64(a.Data(),b.Data()); }
  INLINE SIMD<mask64,4> operator!= (SIMD<int64_t,4> a , SIMD<int64_t,4> b)
  { return  _mm256_xor_si256(_mm256_cmpeq_epi64(a.Data(),b.Data()),_mm256_set1_epi32(-1)); }

#ifdef __AVX2__
  INLINE SIMD<mask64,4> operator&& (SIMD<mask64,4> a, SIMD<mask64,4> b)
  { return _mm256_and_si256 (a.Data(), b.Data()); }
  INLINE SIMD<mask64,4> operator|| (SIMD<mask64,4> a, SIMD<mask64,4> b)
  { return _mm256_or_si256 (a.Data(), b.Data()); }
  INLINE SIMD<mask64,4> operator! (SIMD<mask64,4> a)
  { return _mm256_xor_si256 (a.Data(), _mm256_cmpeq_epi64(a.Data(),a.Data())); }
#else //AVX2 is a superset of AVX. Without it, it is necessary to reinterpret the types
  INLINE SIMD<mask64,4> operator&& (SIMD<mask64,4> a, SIMD<mask64,4> b)
  { return _mm256_castpd_si256(_mm256_and_pd (_mm256_castsi256_pd(a.Data()),_mm256_castsi256_pd( b.Data()))); }
  INLINE SIMD<mask64,4> operator|| (SIMD<mask64,4> a, SIMD<mask64,4> b)
  { return _mm256_castpd_si256(_mm256_or_pd (_mm256_castsi256_pd(a.Data()), _mm256_castsi256_pd(b.Data()))); }
  INLINE SIMD<mask64,4> operator! (SIMD<mask64,4> a)
  { return _mm256_castpd_si256(_mm256_xor_pd (_mm256_castsi256_pd(a.Data()),_mm256_castsi256_pd( _mm256_cmpeq_epi64(a.Data(),a.Data())))); }
#endif
  INLINE SIMD<double,4> If (SIMD<mask64,4> a, SIMD<double,4> b, SIMD<double,4> c)
  { return _mm256_blendv_pd(c.Data(), b.Data(), _mm256_castsi256_pd(a.Data())); }
  
  INLINE SIMD<double,4> IfPos (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    auto cp = _mm256_cmp_pd (a.Data(), _mm256_setzero_pd(), _CMP_GT_OS);
    return _mm256_blendv_pd(c.Data(), b.Data(), cp);
  }

  INLINE SIMD<double,4> IfZero (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    auto cp = _mm256_cmp_pd (a.Data(), _mm256_setzero_pd(), _CMP_EQ_OS);
    return _mm256_blendv_pd(c.Data(), b.Data(), cp);
  }

  INLINE double HSum (SIMD<double,4> sd)
  {
    // __m128d hv = _mm_add_pd (_mm256_extractf128_pd(sd.Data(),0), _mm256_extractf128_pd(sd.Data(),1));
    __m128d hv = (sd.Lo()+sd.Hi()).Data();
    return _mm_cvtsd_f64 (_mm_hadd_pd (hv, hv));
  }

  INLINE auto HSum (SIMD<double,4> sd1, SIMD<double,4> sd2)
  {
    __m256d hv = _mm256_hadd_pd(sd1.Data(), sd2.Data());
    __m128d hv2 = _mm_add_pd (_mm256_extractf128_pd(hv,0), _mm256_extractf128_pd(hv,1));
    return SIMD<double,2>(_mm_cvtsd_f64 (hv2),  _mm_cvtsd_f64(_mm_shuffle_pd (hv2, hv2, 3)));
  }

  INLINE auto HSum (SIMD<double,4> v1, SIMD<double,4> v2, SIMD<double,4> v3, SIMD<double,4> v4)
  {
    __m256d hsum1 = _mm256_hadd_pd (v1.Data(), v2.Data());
    __m256d hsum2 = _mm256_hadd_pd (v3.Data(), v4.Data());
    SIMD<double,4> hsum = _mm256_add_pd (_mm256_permute2f128_pd (hsum1, hsum2, 1+2*16),
                                         _mm256_blend_pd (hsum1, hsum2, 12));
    return hsum;
    // return make_tuple(hsum[0], hsum[1], hsum[2], hsum[3]);
  }

  
  INLINE SIMD<int64_t,4> If (SIMD<mask64,4> a, SIMD<int64_t,4> b, SIMD<int64_t,4> c)
  { return _mm256_castpd_si256(_mm256_blendv_pd(_mm256_castsi256_pd(c.Data()), _mm256_castsi256_pd(b.Data()),
                                                _mm256_castsi256_pd(a.Data()))); }

}

#endif // NG_SIMD_AVX_HPP

