#ifndef NG_SIMD_AVX_HPP
#define NG_SIMD_AVX_HPP

/**************************************************************************/
/* File:   simd_avx.hpp                                                       */
/* Author: Joachim Schoeberl, Matthias Hochsteger                         */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

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
    int64_t operator[] (int i) const { return ((int64_t*)(&mask))[i]; }    
  };

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
    INLINE auto & operator[] (int i) { return ((int64_t*)(&data))[i]; }
    INLINE __m256i Data() const { return data; }
    INLINE __m256i & Data() { return data; }

    SIMD<int64_t,2> Lo() const { return _mm256_extractf128_si256(data, 0); }
    SIMD<int64_t,2> Hi() const { return _mm256_extractf128_si256(data, 1); }
    static SIMD FirstInt() { return { 0, 1, 2, 3 }; }
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
}
#endif // NG_SIMD_AVX_HPP

