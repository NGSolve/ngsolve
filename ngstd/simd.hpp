#ifndef FILE_SIMD
#define FILE_SIMD

/**************************************************************************/
/* File:   simd.hpp                                                       */
/* Author: Joachim Schoeberl, Matthias Hochsteger                         */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

// #include <immintrin.h>

#ifdef WIN32
#ifndef AVX_OPERATORS_DEFINED
#define AVX_OPERATORS_DEFINED
INLINE __m128d operator- (__m128d a) { return _mm_xor_pd(a, _mm_set1_pd(-0.0)); }
INLINE __m128d operator+ (__m128d a, __m128d b) { return _mm_add_pd(a,b); }
INLINE __m128d operator- (__m128d a, __m128d b) { return _mm_sub_pd(a,b); }
INLINE __m128d operator* (__m128d a, __m128d b) { return _mm_mul_pd(a,b); }
INLINE __m128d operator/ (__m128d a, __m128d b) { return _mm_div_pd(a,b); }
INLINE __m128d operator* (double a, __m128d b) { return _mm_set1_pd(a)*b; }
INLINE __m128d operator* (__m128d b, double a) { return _mm_set1_pd(a)*b; }

INLINE __m128d operator+= (__m128d &a, __m128d b) { return a = a+b; }
INLINE __m128d operator-= (__m128d &a, __m128d b) { return a = a-b; }
INLINE __m128d operator*= (__m128d &a, __m128d b) { return a = a*b; }
INLINE __m128d operator/= (__m128d &a, __m128d b) { return a = a/b; }

INLINE __m256d operator- (__m256d a) { return _mm256_xor_pd(a, _mm256_set1_pd(-0.0)); }
INLINE __m256d operator+ (__m256d a, __m256d b) { return _mm256_add_pd(a,b); }
INLINE __m256d operator- (__m256d a, __m256d b) { return _mm256_sub_pd(a,b); }
INLINE __m256d operator* (__m256d a, __m256d b) { return _mm256_mul_pd(a,b); }
INLINE __m256d operator/ (__m256d a, __m256d b) { return _mm256_div_pd(a,b); }
INLINE __m256d operator* (double a, __m256d b) { return _mm256_set1_pd(a)*b; }
INLINE __m256d operator* (__m256d b, double a) { return _mm256_set1_pd(a)*b; }

INLINE __m256d operator+= (__m256d &a, __m256d b) { return a = a+b; }
INLINE __m256d operator-= (__m256d &a, __m256d b) { return a = a-b; }
INLINE __m256d operator*= (__m256d &a, __m256d b) { return a = a*b; }
INLINE __m256d operator/= (__m256d &a, __m256d b) { return a = a/b; }

/*Some integer arithmetic operations are not defined. Check
https://github.com/vectorclass/version2/blob/master/vectori128.h
for possible implementations*/

INLINE __m128i operator-(__m128i a) { return _mm_sub_epi64(_mm_setzero_si128(), a); }
INLINE __m128i operator+ (__m128i a, __m128i b) { return _mm_add_epi64(a,b); }
INLINE __m128i operator- (__m128i a, __m128i b) { return _mm_sub_epi64(a,b); }

INLINE __m128i operator+= (__m128i &a, __m128i b) { return a = a+b; }
INLINE __m128i operator-= (__m128i &a, __m128i b) { return a = a-b; }

INLINE __m256i operator-(__m256i a) { return _mm256_sub_epi64(_mm256_setzero_si256(), a); }

#ifdef __AVX2__
INLINE __m256i operator+ (__m256i a, __m256i b) { return _mm256_add_epi64(a,b); }
INLINE __m256i operator- (__m256i a, __m256i b) { return _mm256_sub_epi64(a,b); }
#else
INLINE __m256i operator+ (__m256i a, __m256i b) {
  auto lo_sum = _mm256_extractf128_si256(a, 0) + _mm256_extractf128_si256(b, 0);
  auto hi_sum = _mm256_extractf128_si256(a, 1) + _mm256_extractf128_si256(b, 1);
  return _mm256_set_m128i(hi_sum,lo_sum);
}
INLINE __m256i operator- (__m256i a, __m256i b) {
  auto lo_sub = _mm256_extractf128_si256(a, 0) - _mm256_extractf128_si256(b, 0);
  auto hi_sub = _mm256_extractf128_si256(a, 1) - _mm256_extractf128_si256(b, 1);
  return _mm256_set_m128i(hi_sub,lo_sub);
}
#endif
INLINE __m256i operator+= (__m256i &a, __m256i b) { return a = a+b; }
INLINE __m256i operator-= (__m256i &a, __m256i b) { return a = a-b; }

#endif
#endif



namespace ngstd
{

  // MSVC does not define SSE. It's always present on 64bit cpus
#if (defined(_M_AMD64) || defined(_M_X64) || defined(__AVX__))
#ifndef __SSE__
#define __SSE__
#endif
#ifndef __SSE2__
#define __SSE2__
#endif
#endif


  
  constexpr int GetDefaultSIMDSize() {
#if defined __AVX512F__
    return 8;
#elif defined __AVX__
    return 4;
#elif defined __SSE__
    return 2;
#else
    return 1;
#endif
  }
  
#if defined __AVX512F__
    typedef __m512 tAVX;
    typedef __m512d tAVXd;
#elif defined __AVX__
    typedef __m256 tAVX;
    typedef __m256d tAVXd; 
#elif defined __SSE__
    typedef __m128 tAVX;
    typedef __m128d tAVXd; 
#endif

  template <typename T, int N=GetDefaultSIMDSize()> class SIMD;

  

  /*
  template <typename T>
  class AlignedAlloc
  {
  protected:
    static void * aligned_malloc(size_t s)
    {
      // Assume 16 byte alignment of standard library
#ifdef __SSE__        
      if(alignof(T)>16)
        return  _mm_malloc(s, alignof(T));
      else
#endif          
        // return malloc(s);
        return new char[s];
    }
    
      static void aligned_free(void *p)
      {
#ifdef __SSE__        
        if(alignof(T)>16)
          _mm_free(p);
        else
#endif
          // free(p);
          delete (char*)p;
      }

  public:
    void * operator new (size_t s, void *p) { return p; }
    void * operator new (size_t s) { return aligned_malloc(s); }
    void * operator new[] (size_t s) { return aligned_malloc(s); }
    void operator delete (void * p) { aligned_free(p); }
    void operator delete[] (void * p) { aligned_free(p); }
  };


#else
  */  
  // it's only a dummy without AVX
  template <typename T>
  class AlignedAlloc { ; };

  // #endif

#ifndef __AVX__
  INLINE __m128i my_mm_cmpgt_epi64(__m128i a, __m128i b) {
    auto  res_lo = _mm_cvtsi128_si64(a)  > _mm_cvtsi128_si64(b) ? -1:0;
    auto  res_hi = _mm_cvtsi128_si64(_mm_srli_si128(a,8)) > _mm_cvtsi128_si64(_mm_srli_si128(b,8)) ? -1 : 0;
    return _mm_set_epi64x(res_hi,res_lo);
  }
#else
  INLINE __m128i my_mm_cmpgt_epi64(__m128i a, __m128i b) {
    return _mm_cmpgt_epi64(a,b);
  }
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
#endif

  // typedef int64_t mask64;
  class mask64;
  
  template <> 
  class SIMD<mask64,1>
  {
    int64_t mask;
  public:
    SIMD (size_t i)
      : mask(i > 0 ? -1 : 0) { ; }
    bool Data() const { return mask; }
    static constexpr int Size() { return 1; }    
    auto operator[] (int /* i */) const { return mask; }
  };


#ifdef __SSE__
  template <> 
  class SIMD<mask64,2>
  {
    __m128i mask;
  public:
    SIMD (int i)
      : mask(_mm_cmpgt_epi32(_mm_set1_epi32(i),
                             _mm_set_epi32(1, 1, 0, 0)))
    { ; }
    SIMD (__m128i _mask) : mask(_mask) { ; }
    __m128i Data() const { return mask; }
    static constexpr int Size() { return 2; }    
    int64_t operator[] (int i) const { return ((int64_t*)(&mask))[i]; }    
  };
#endif
  
  
#ifdef __AVX__
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
#else
  template <> 
  class SIMD<mask64,4>
  {
    SIMD<mask64,2> lo, hi;
  public:
    SIMD (int i) : lo(i), hi(i-2) { ; } 
    SIMD<mask64,2> Lo() const { return lo; }
    SIMD<mask64,2> Hi() const { return hi; }
    static constexpr int Size() { return 4; }    
    // mask64 operator[] (int i) const { return ((mask64*)(&mask))[i]; }    
  };
  
#endif


#ifdef __AVX512F__
  template <> 
  class SIMD<mask64,8>
  {
    __mmask8 mask;
  public:
    SIMD (size_t i)
      : mask(_mm512_cmpgt_epi64_mask(_mm512_set1_epi64(i),
                                     _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)))
    { ; }
    SIMD (int i)
      : mask(_mm512_cmpgt_epi64_mask(_mm512_set1_epi64(i),
                                     _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)))
    { ; }
    SIMD (int64_t i)
      : mask(_mm512_cmpgt_epi64_mask(_mm512_set1_epi64(i),
                                     _mm512_set_epi64(7, 6, 5, 4, 3, 2, 1, 0)))
    { ; }
    SIMD (__mmask8 _mask) : mask(_mask) { ; }        
    __mmask8 Data() const { return mask; }
    static constexpr int Size() { return 8; }    
    // mask64 operator[] (int i) const { return ((mask64*)(&mask))[i]; }    
  };
#endif


  
  template<>
  class SIMD<int64_t,1>
  {
    int64_t data;
    
  public:
    static constexpr int Size() { return 1; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;
    SIMD (int64_t val) { data = val; }

    int64_t operator[] (int i) const { return ((int64_t*)(&data))[i]; }
    auto Data() const { return data; }
    auto & Data() { return data; }
  };

  
  template<>
  class SIMD<double,1>
  {
    double data;
    
  public:
    static constexpr int Size() { return 1; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;
    SIMD (double val) { data = val; }
    SIMD (int val)    { data = val; }
    SIMD (size_t val) { data = val; }
    SIMD (double const * p) { data = *p; }
    SIMD (double const * p, SIMD<mask64,1> mask) { data = mask.Data() ? *p : 0.0; }
    
    template <typename T, typename std::enable_if<std::is_convertible<T,std::function<double(int)>>::value,int>::type = 0>
    SIMD (const T & func)
    {
      data = func(0);
    }
    
    template <typename T, typename std::enable_if<std::is_convertible<T,std::function<double(int)>>::value,int>::type = 0>
    SIMD & operator= (const T & func)
    {
      data = func(0);
      return *this;
    }

    void Store (double * p) { *p = data; }
    void Store (double * p, SIMD<mask64,1> mask) { if (mask.Data()) *p = data; }
    
    double operator[] (int i) const { return ((double*)(&data))[i]; }
    double Data() const { return data; }
    double & Data() { return data; }
  };
  

#ifdef __SSE__

  template<>
  class SIMD<int64_t,2> 
  {
    __m128i data;
    
  public:
    static constexpr int Size() { return 2; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD (int64_t v0, int64_t v1) { data = _mm_set_epi64x(v1,v0); }
    
    SIMD & operator= (const SIMD &) = default;

    SIMD (int64_t val) { data = _mm_set1_epi64x(val); }
    SIMD (__m128i _data) { data = _data; }

    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<int64_t(int)>>::value, int>::type = 0>
    SIMD (const T & func)
    {   
      data = _mm_set_epi64(func(1), func(0));             
    }
    
    INLINE auto operator[] (int i) const { return ((int64_t*)(&data))[i]; }
    INLINE auto & operator[] (int i) { return ((int64_t*)(&data))[i]; }
    INLINE __m128i Data() const { return data; }
    INLINE __m128i & Data() { return data; }
    static SIMD FirstInt() { return { 0, 1 }; }    
  };


  
  template<>
  class alignas(16) SIMD<double,2> : public AlignedAlloc<SIMD<double,2>>
  {
    __m128d data;
    
  public:
    static constexpr int Size() { return 2; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD (double v0, double v1) { data = _mm_set_pd(v1,v0); }
    
    SIMD & operator= (const SIMD &) = default;

    SIMD (double val) { data = _mm_set1_pd(val); }
    SIMD (int val)    { data = _mm_set1_pd(val); }
    SIMD (size_t val) { data = _mm_set1_pd(val); }

    SIMD (double const * p) { data = _mm_loadu_pd(p); }
    SIMD (double const * p, SIMD<mask64,2> mask)
      {
#ifdef __AVX__
        data = _mm_maskload_pd(p, mask.Data());
#else
        // this versions segfaults if p points to the last allowed element
        // happened on Mac with the new SparseCholesky-factorization
        // data = _mm_and_pd(_mm_castsi128_pd(mask.Data()), _mm_loadu_pd(p));
        data = _mm_set_pd (mask[1] ? p[1] : 0.0, mask[0] ? p[0] : 0.0);        
#endif
      }
    SIMD (__m128d _data) { data = _data; }

    void Store (double * p) { _mm_storeu_pd(p, data); }
    void Store (double * p, SIMD<mask64,2> mask)
    {
#ifdef __AVX__
      _mm_maskstore_pd(p, mask.Data(), data);
#else
      /*
      _mm_storeu_pd (p, _mm_or_pd (_mm_and_pd(_mm_castsi128_pd(mask.Data()), data),
                                   _mm_andnot_pd(_mm_castsi128_pd(mask.Data()), _mm_loadu_pd(p))));
      */
      if (mask[0]) p[0] = (*this)[0];
      if (mask[1]) p[1] = (*this)[1];
#endif
    }    
    
    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>                                                                    SIMD (const T & func)
    {   
      data = _mm_set_pd(func(1), func(0));              
    }   
    
    INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    INLINE double & operator[] (int i) { return ((double*)(&data))[i]; }
    INLINE __m128d Data() const { return data; }
    INLINE __m128d & Data() { return data; }

    operator tuple<double&,double&> ()
    { return tuple<double&,double&>((*this)[0], (*this)[1]); }
  };

  INLINE auto Unpack (SIMD<double,2> a, SIMD<double,2> b)
  {
    return make_tuple(SIMD<double,2>(_mm_unpacklo_pd(a.Data(),b.Data())),
                      SIMD<double,2>(_mm_unpackhi_pd(a.Data(),b.Data())));
  }
  
#endif

  
  

#ifdef __AVX__

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


  
  template<>
  class /* alignas(32) */ SIMD<double,4> : public AlignedAlloc<SIMD<double,4>>
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
    // [[deprecated("don't write to individual elments of SIMD")]]                
    // INLINE double & operator[] (int i) { return ((double*)(&data))[i]; }
    template <int I>
    double Get() const { return ((double*)(&data))[I]; }
    INLINE __m256d Data() const { return data; }
    INLINE __m256d & Data() { return data; }

    SIMD<double,2> Lo() const { return _mm256_extractf128_pd(data, 0); }
    SIMD<double,2> Hi() const { return _mm256_extractf128_pd(data, 1); }

    // operator tuple<double&,double&,double&,double&> ()
    // { return tuple<double&,double&,double&,double&>((*this)[0], (*this)[1], (*this)[2], (*this)[3]); }
  };

  INLINE auto Unpack (SIMD<double,4> a, SIMD<double,4> b)
  {
    return make_tuple(SIMD<double,4>(_mm256_unpacklo_pd(a.Data(),b.Data())),
                      SIMD<double,4>(_mm256_unpackhi_pd(a.Data(),b.Data())));
  }
  

#else // AVX

  template<>
  class /* alignas(32) */ SIMD<int64_t,4> : public AlignedAlloc<SIMD<int64_t,4>>
  {
    SIMD<int64_t,2> data[2];
    
  public:
    static constexpr int Size() { return 4; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;
    
    SIMD (int64_t val) : data{val,val} { ; }
    SIMD (int64_t v0, int64_t v1, int64_t v2, int64_t v3)
    {
      data[0] = SIMD<int64_t,2>(v0,v1);
      data[1] = SIMD<int64_t,2>(v2,v3);
    }    
    SIMD (SIMD<int64_t,2> lo, SIMD<int64_t,2> hi) : data{lo,hi} { ; }

    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<int64_t(int)>>::value, int>::type = 0>
      SIMD (const T & func)
    {   
      data = _mm256_set_epi64x(func(3), func(2), func(1), func(0));             
    }
    
    auto Lo() const { return data[0]; }
    auto & Lo() { return data[0]; }
    auto Hi() const { return data[1]; }
    auto & Hi() { return data[1]; }
    
    INLINE int64_t operator[] (int i) const { return ((int64_t*)(&data[0]))[i]; }
    INLINE int64_t & operator[] (int i) { return ((int64_t*)(&data[0]))[i]; }

    operator tuple<int64_t&,int64_t&,int64_t&,int64_t&> ()
    { return tuple<int64_t&,int64_t&,int64_t&,int64_t&>((*this)[0], (*this)[1], (*this)[2], (*this)[3]); }

    static SIMD FirstInt() { return { 0, 1, 2, 3 }; }
  };
  
  template<>
  class /* alignas(32) */ SIMD<double,4> : public AlignedAlloc<SIMD<double,4>>
  {
    SIMD<double,2> data[2];
    
  public:
    static constexpr int Size() { return 4; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD (SIMD<double,2> lo, SIMD<double,2> hi) : data{lo,hi} { ; }
    SIMD (double v0, double v1, double v2, double v3)
      {
        data[0] = SIMD<double,2>(v0,v1);
        data[1] = SIMD<double,2>(v2,v3);
      }
    
    SIMD & operator= (const SIMD &) = default;

    SIMD (double val) : data{val,val} { ; }
    SIMD (int val)    : data{val,val} { ; } 
    SIMD (size_t val) : data{val,val} { ; } 

    SIMD (double const * p) : data{p,p+2} { ; }
    SIMD (double const * p, SIMD<mask64,4> mask)
      {
        data[0] = SIMD<double,2> (p, mask.Lo());
        data[1] = SIMD<double,2> (p+2, mask.Hi());
      }

    void Store (double * p) { data[0].Store(p); data[1].Store(p+2); }
    void Store (double * p, SIMD<mask64,4> mask)
    {
      data[0].Store(p, mask.Lo());
      data[1].Store(p+2, mask.Hi());
    }    

    /*
    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>                                                                    SIMD (const T & func)
    {   
      data[0] =  = _mm256_set_pd(func(3), func(2), func(1), func(0));              
    }   
    */

    auto Lo() const { return data[0]; }
    auto & Lo() { return data[0]; }
    auto Hi() const { return data[1]; }
    auto & Hi() { return data[1]; }
    
    INLINE double operator[] (int i) const { return ((double*)(&data[0]))[i]; }
    INLINE double & operator[] (int i) { return ((double*)(&data[0]))[i]; }
    // INLINE __m256d Data() const { return data; }
    // INLINE __m256d & Data() { return data; }

    operator tuple<double&,double&,double&,double&> ()
    { return tuple<double&,double&,double&,double&>((*this)[0], (*this)[1], (*this)[2], (*this)[3]); }
  };


  INLINE auto Unpack (SIMD<double,4> a, SIMD<double,4> b)
  {
    SIMD<double,2> h1,h2,h3,h4;
    tie(h1,h2) = Unpack(a.Lo(), b.Lo());
    tie(h3,h4) = Unpack(a.Hi(), b.Hi());
    return make_tuple(SIMD<double,4> (h1,h3), SIMD<double,4>(h2,h4));
  }  
  
#endif



  
#ifdef __AVX512F__
  template<>
  class /* alignas(64) */ SIMD<double,8> : public AlignedAlloc<SIMD<double,8>>
  {
    __m512d data;
  public:
    static constexpr int Size() { return 8; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (double val) { data = _mm512_set1_pd(val); }
    SIMD (int val)    { data = _mm512_set1_pd(val); }
    SIMD (size_t val) { data = _mm512_set1_pd(val); }
    SIMD (double const * p) { data = _mm512_loadu_pd(p); }
    SIMD (double const * p, SIMD<mask64,8> mask)
      { data = _mm512_mask_loadu_pd(_mm512_setzero_pd(), mask.Data(), p); }
    SIMD (__m512d _data) { data = _data; }
    
    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>
      SIMD (const T & func)
    {   
      data = _mm512_set_pd(func(7), func(6), func(5), func(4), func(3), func(2), func(1), func(0));              
    }

    void Store (double * p) { _mm512_storeu_pd(p, data); }
    void Store (double * p, SIMD<mask64,8> mask) { _mm512_mask_storeu_pd(p, mask.Data(), data); }    
    
    /*
    template <typename T>
    SIMD (const T & val)
    {
//       SIMD_function(val, std::is_convertible<T, std::function<double(int)>>());
      SIMD_function(val, has_call_operator<T>::value);
    }
    
    template <typename T>
    SIMD & operator= (const T & val)
    {
//       SIMD_function(val, std::is_convertible<T, std::function<double(int)>>());
      SIMD_function(val, has_call_operator<T>::value);
      return *this;
    }
    */
    
    template <typename Function>
    void SIMD_function (const Function & func, std::true_type)
    {
      /*
      data = _mm512_set_pd(func(7), func(6), func(5), func(4),
                           func(3), func(2), func(1), func(0));
      */
      data = (__m512){ func(7), func(6), func(5), func(4),
                       func(3), func(2), func(1), func(0) };
                       
      
    }
    
    // not a function
    void SIMD_function (double const * p, std::false_type)
    {
      data = _mm512_loadu_pd(p);
    }
    
    void SIMD_function (double val, std::false_type)
    {
      data = _mm512_set1_pd(val);
    }
    
    void SIMD_function (__m512d _data, std::false_type)
    {
      data = _data;
    }
    
    INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    INLINE double & operator[] (int i) { return ((double*)(&data))[i]; }
    INLINE __m512d Data() const { return data; }
    INLINE __m512d & Data() { return data; }
  };

  template<>
  class SIMD<int64_t,8> 
  {
    __m512i data;
    
  public:
    static constexpr int Size() { return 8; }
    SIMD () {}
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (int64_t val) { data = _mm512_set1_epi64(val); }
    SIMD (int64_t v0, int64_t v1, int64_t v2, int64_t v3, int64_t v4, int64_t v5, int64_t v6, int64_t v7) { data = _mm512_set_epi64(v7,v6,v5,v4,v3,v2,v1,v0); }
    SIMD (__m512i _data) { data = _data; }

    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<int64_t(int)>>::value, int>::type = 0>
      SIMD (const T & func)
    {   
      data = _mm512_set_epi64(func(7), func(6), func(5), func(4), func(3), func(2), func(1), func(0));              
    }
    

    INLINE auto operator[] (int i) const { return ((int64_t*)(&data))[i]; }
    INLINE auto & operator[] (int i) { return ((int64_t*)(&data))[i]; }
    INLINE __m512i Data() const { return data; }
    INLINE __m512i & Data() { return data; }
    static SIMD FirstInt() { return { 0, 1, 2, 3, 4, 5, 6, 7 }; }
    // SIMD<double,4> Lo() const { return _mm512_extracti64x4_epi64(data, 0); }
    // SIMD<double,4> Hi() const { return _mm512_extracti64x4_epi64(data, 1); }
  };
#endif
  

  //integer arithmetic operators
    template <int N>
  INLINE SIMD<int64_t,N> operator+ (SIMD<int64_t,N> a, SIMD<int64_t,N> b) { return a.Data()+b.Data(); }
#ifndef __AVX__
  INLINE SIMD<int64_t,4> operator+ (SIMD<int64_t,4> a, SIMD<int64_t,4> b) { return { a.Lo()+b.Lo(), a.Hi()+b.Hi() }; }
#endif
  
  template <int N>
  INLINE SIMD<int64_t,N> operator+ (SIMD<int64_t,N> a, int64_t b) { return a+SIMD<int64_t,N>(b); }
  template <int N>
  INLINE SIMD<int64_t,N> operator+ (int64_t a, SIMD<int64_t,N> b) { return SIMD<int64_t,N>(a)+b; }
  template <int N>  
  INLINE SIMD<int64_t,N> operator- (SIMD<int64_t,N> a, SIMD<int64_t,N> b) { return a.Data()-b.Data(); }
#ifndef __AVX__
  INLINE SIMD<int64_t,4> operator- (SIMD<int64_t,4> a, SIMD<int64_t,4> b) { return { a.Lo()-b.Lo(), a.Hi()-b.Hi() }; }
#endif

  template <int N>  
  INLINE SIMD<int64_t,N> operator- (int64_t a, SIMD<int64_t,N> b) { return SIMD<int64_t,N>(a)-b; }
  template <int N>  
  INLINE SIMD<int64_t,N> operator- (SIMD<int64_t,N> a, int64_t b) { return a-SIMD<int64_t,N>(b); }
  template <int N>  
  INLINE SIMD<int64_t,N> operator- (SIMD<int64_t,N> a) { return -a.Data(); }
#ifndef __AVX__
  INLINE SIMD<int64_t,4> operator- (SIMD<int64_t,4> a) { return { -a.Lo(), -a.Hi() }; }
#endif
  
  
  template <int N>  
  INLINE SIMD<int64_t,N> & operator+= (SIMD<int64_t,N> & a, SIMD<int64_t,N> b) { a=a+b; return a; }
  template <int N>  
  INLINE SIMD<int64_t,N> & operator+= (SIMD<int64_t,N> & a, int64_t b) { a+=SIMD<int64_t,N>(b); return a; }
  template <int N>  
  INLINE SIMD<int64_t,N> & operator-= (SIMD<int64_t,N> & a, SIMD<int64_t,N> b) { a.Data()-=b.Data(); return a; }
  template <int N>  
  INLINE SIMD<int64_t,N> & operator-= (SIMD<int64_t,N> & a, int64_t b) { a-=SIMD<int64_t,N>(b); return a; }

  
  template <int N>
  INLINE SIMD<double,N> operator+ (SIMD<double,N> a, SIMD<double,N> b) { return a.Data()+b.Data(); }
#ifndef __AVX__
  INLINE SIMD<double,4> operator+ (SIMD<double,4> a, SIMD<double,4> b) { return { a.Lo()+b.Lo(), a.Hi()+b.Hi() }; }
#endif
  
  template <int N>
  INLINE SIMD<double,N> operator+ (SIMD<double,N> a, double b) { return a+SIMD<double,N>(b); }
  template <int N>
  INLINE SIMD<double,N> operator+ (double a, SIMD<double,N> b) { return SIMD<double,N>(a)+b; }
  template <int N>  
  INLINE SIMD<double,N> operator- (SIMD<double,N> a, SIMD<double,N> b) { return a.Data()-b.Data(); }
#ifndef __AVX__
  INLINE SIMD<double,4> operator- (SIMD<double,4> a, SIMD<double,4> b) { return { a.Lo()-b.Lo(), a.Hi()-b.Hi() }; }
#endif

  template <int N>  
  INLINE SIMD<double,N> operator- (double a, SIMD<double,N> b) { return SIMD<double,N>(a)-b; }
  template <int N>  
  INLINE SIMD<double,N> operator- (SIMD<double,N> a, double b) { return a-SIMD<double,N>(b); }
  template <int N>  
  INLINE SIMD<double,N> operator- (SIMD<double,N> a) { return -a.Data(); }
#ifndef __AVX__
  INLINE SIMD<double,4> operator- (SIMD<double,4> a) { return { -a.Lo(), -a.Hi() }; }
#endif
  
  template <int N>  
  INLINE SIMD<double,N> operator* (SIMD<double,N> a, SIMD<double,N> b) { return a.Data()*b.Data(); }
#ifndef __AVX__
  INLINE SIMD<double,4> operator* (SIMD<double,4> a, SIMD<double,4> b) { return { a.Lo()*b.Lo(), a.Hi()*b.Hi() }; }
#endif
  template <int N>  
  INLINE SIMD<double,N> operator* (double a, SIMD<double,N> b) { return SIMD<double,N>(a)*b; }
  template <int N>  
  INLINE SIMD<double,N> operator* (SIMD<double,N> b, double a) { return SIMD<double,N>(a)*b; }
  template <int N>  
  INLINE SIMD<double,N> operator/ (SIMD<double,N> a, SIMD<double,N> b) { return a.Data()/b.Data(); }
  template <int N>  
  INLINE SIMD<double,N> operator/ (SIMD<double,N> a, double b) { return a/SIMD<double,N>(b); }
  template <int N>  
  INLINE SIMD<double,N> operator/ (double a, SIMD<double,N> b) { return SIMD<double,N>(a)/b; }
  template <int N>  
  INLINE SIMD<double,N> & operator+= (SIMD<double,N> & a, SIMD<double,N> b) { a=a+b; return a; }
  template <int N>  
  INLINE SIMD<double,N> & operator+= (SIMD<double,N> & a, double b) { a+=SIMD<double,N>(b); return a; }
  template <int N>  
  INLINE SIMD<double,N> & operator-= (SIMD<double,N> & a, SIMD<double,N> b) { a.Data()-=b.Data(); return a; }
  template <int N>  
  INLINE SIMD<double,N> & operator-= (SIMD<double,N> & a, double b) { a-=SIMD<double,N>(b); return a; }
  template <int N>  
  INLINE SIMD<double,N> & operator*= (SIMD<double,N> & a, SIMD<double,N> b) { a=a*b; return a; }
  template <int N>  
  INLINE SIMD<double,N> & operator*= (SIMD<double,N> & a, double b) { a*=SIMD<double,N>(b); return a; }
  template <int N>  
  INLINE SIMD<double,N> & operator/= (SIMD<double,N> & a, SIMD<double,N> b) { a.Data()/=b.Data(); return a; }

  template <int N>    
  INLINE SIMD<double,N> L2Norm2 (SIMD<double,N> a) { return a.Data()*a.Data(); }
  template <int N>
  INLINE SIMD<double,N> Trans (SIMD<double,N> a) { return a; }


  INLINE double IfPos (double a, double b, double c) { return a>0 ? b : c; }
  INLINE double IfZero (double a, double b, double c) { return a==0. ? b : c; }
  
#ifdef __SSE__

  
  INLINE __m128d my_mm_hadd_pd(__m128d a, __m128d b) {
#if defined(__SSE3__) || defined(__AVX__)
    return _mm_hadd_pd(a,b); 
#else
    return _mm_add_pd( _mm_unpacklo_pd(a,b), _mm_unpackhi_pd(a,b) );
#endif
  }

  INLINE SIMD<double,2> sqrt (SIMD<double,2> a) { return _mm_sqrt_pd(a.Data()); }
  INLINE SIMD<double,2> fabs (SIMD<double,2> a) { return _mm_max_pd(a.Data(), -a.Data()); }
  using std::floor;
  INLINE SIMD<double,2> floor (SIMD<double,2> a)
  { return ngstd::SIMD<double,2>([&](int i)->double { return floor(a[i]); } ); }
  using std::ceil;  
  INLINE SIMD<double,2> ceil (SIMD<double,2> a) 
  { return ngstd::SIMD<double,2>([&](int i)->double { return ceil(a[i]); } ); }

  INLINE SIMD<mask64,2> operator<= (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmple_pd(a.Data(),b.Data())); }
  INLINE SIMD<mask64,2> operator< (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmplt_pd(a.Data(),b.Data())); }
  INLINE SIMD<mask64,2> operator>= (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmpge_pd(a.Data(),b.Data())); }
  INLINE SIMD<mask64,2> operator> (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmpgt_pd(a.Data(),b.Data())); }
  INLINE SIMD<mask64,2> operator== (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmpeq_pd(a.Data(),b.Data())); }
  INLINE SIMD<mask64,2> operator!= (SIMD<double,2> a , SIMD<double,2> b)
  { return _mm_castpd_si128( _mm_cmpneq_pd(a.Data(),b.Data())); }

  INLINE SIMD<mask64,2> operator<= (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  _mm_xor_si128(_mm_cmpgt_epi64(a.Data(),b.Data()),_mm_set1_epi32(-1)); }
  INLINE SIMD<mask64,2> operator< (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  my_mm_cmpgt_epi64(b.Data(),a.Data()); }
  INLINE SIMD<mask64,2> operator>= (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  _mm_xor_si128(_mm_cmpgt_epi64(b.Data(),a.Data()),_mm_set1_epi32(-1)); }
  INLINE SIMD<mask64,2> operator> (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  my_mm_cmpgt_epi64(a.Data(),b.Data()); }
  INLINE SIMD<mask64,2> operator== (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  _mm_cmpeq_epi64(a.Data(),b.Data()); }
  INLINE SIMD<mask64,2> operator!= (SIMD<int64_t,2> a , SIMD<int64_t,2> b)
  { return  _mm_xor_si128(_mm_cmpeq_epi64(a.Data(),b.Data()),_mm_set1_epi32(-1)); }

  
  
 INLINE SIMD<mask64,2> operator&& (SIMD<mask64,2> a, SIMD<mask64,2> b)
  { return _mm_castpd_si128(_mm_and_pd (_mm_castsi128_pd(a.Data()),_mm_castsi128_pd( b.Data()))); }
  INLINE SIMD<mask64,2> operator|| (SIMD<mask64,2> a, SIMD<mask64,2> b)
  { return _mm_castpd_si128(_mm_or_pd (_mm_castsi128_pd(a.Data()), _mm_castsi128_pd(b.Data()))); }
  INLINE SIMD<mask64,2> operator! (SIMD<mask64,2> a)
  { return _mm_castpd_si128(_mm_xor_pd (_mm_castsi128_pd(a.Data()),_mm_castsi128_pd( _mm_cmpeq_epi64(a.Data(),a.Data())))); }
#ifdef __SSE4_1__
  INLINE SIMD<double,2> If (SIMD<mask64,2> a, SIMD<double,2> b, SIMD<double,2> c)
  { return _mm_blendv_pd(c.Data(), b.Data(), _mm_castsi128_pd(a.Data())); }
#else
  INLINE SIMD<double,2> If (SIMD<mask64,2> a, SIMD<double,2> b, SIMD<double,2> c)
  {
    return _mm_or_pd(
                      _mm_andnot_pd(_mm_castsi128_pd(a.Data()),c.Data()),
                      _mm_and_pd(b.Data(),_mm_castsi128_pd(a.Data()))
                      );}
#endif

  
  INLINE SIMD<double,2> IfPos (SIMD<double,2> a, SIMD<double,2> b, SIMD<double,2> c)
  { return ngstd::SIMD<double,2>([&](int i)->double { return a[i]>0 ? b[i] : c[i]; }); }
  INLINE SIMD<double,2> IfZero (SIMD<double,2> a, SIMD<double,2> b, SIMD<double,2> c)
  { return ngstd::SIMD<double,2>([&](int i)->double { return a[i]==0. ? b[i] : c[i]; }); }

  
  INLINE double HSum (SIMD<double,2> sd)
  {
    return _mm_cvtsd_f64 (my_mm_hadd_pd (sd.Data(), sd.Data()));
  }

  INLINE auto HSum (SIMD<double,2> sd1, SIMD<double,2> sd2)
  {
    __m128d hv2 = my_mm_hadd_pd(sd1.Data(), sd2.Data());
    return SIMD<double,2> (hv2);
    // return SIMD<double,2>(_mm_cvtsd_f64 (hv2),  _mm_cvtsd_f64(_mm_shuffle_pd (hv2, hv2, 3)));
  }

  INLINE auto HSum (SIMD<double,2> v1, SIMD<double,2> v2, SIMD<double,2> v3, SIMD<double,2> v4)
  {
    SIMD<double,2> hsum1 = my_mm_hadd_pd (v1.Data(), v2.Data());
    SIMD<double,2> hsum2 = my_mm_hadd_pd (v3.Data(), v4.Data());
    return SIMD<double,4> (hsum1, hsum2);
  }

  INLINE SIMD<int64_t, 2> If(SIMD<mask64, 2> a, SIMD<int64_t, 2> b,
                             SIMD<int64_t, 2> c) {
    return _mm_or_si128(
                        _mm_andnot_si128(a.Data(),c.Data()),
                        _mm_and_si128(b.Data(),a.Data())
                        );
  }
#endif

  
#ifdef __AVX__
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

  INLINE SIMD<double,4> LoadIndirect (double * p, SIMD<int64_t,4> index)
  { return _mm256_i64gather_pd (p, index.Data(), 8); }


#else  // AVX

  INLINE SIMD<double,4> If (SIMD<mask64,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    return SIMD<double,4> (If(a.Lo(), b.Lo(), c.Lo()), If(a.Hi(), b.Hi(), c.Hi()));
  }

  
#endif  
  

#ifdef __AVX512F__
  INLINE SIMD<double,8> sqrt (SIMD<double,8> a) { return _mm512_sqrt_pd(a.Data()); }
  INLINE SIMD<double,8> floor (SIMD<double,8> a) { return _mm512_floor_pd(a.Data()); }
  INLINE SIMD<double,8> ceil (SIMD<double,8> a) { return _mm512_ceil_pd(a.Data()); }  
  INLINE SIMD<double,8> fabs (SIMD<double,8> a) { return _mm512_max_pd(a.Data(), -a.Data()); }

  INLINE SIMD<mask64,8> operator<= (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_LE_OQ); }
  INLINE SIMD<mask64,8> operator< (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_LT_OQ); }
  INLINE SIMD<mask64,8> operator>= (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_GE_OQ); }
  INLINE SIMD<mask64,8> operator> (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_GT_OQ); }
  INLINE SIMD<mask64,8> operator== (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_EQ_OQ); }
  INLINE SIMD<mask64,8> operator!= (SIMD<double,8> a , SIMD<double,8> b)
  { return _mm512_cmp_pd_mask (a.Data(), b.Data(), _CMP_NEQ_OQ); }

  INLINE SIMD<mask64,8> operator<= (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(), _MM_CMPINT_LE); }
  INLINE SIMD<mask64,8> operator< (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(), _MM_CMPINT_LT); }
  INLINE SIMD<mask64,8> operator>= (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(),  _MM_CMPINT_NLT); }
  INLINE SIMD<mask64,8> operator> (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(), _MM_CMPINT_NLE); }
  INLINE SIMD<mask64,8> operator== (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(), _MM_CMPINT_EQ); }
  INLINE SIMD<mask64,8> operator!= (SIMD<int64_t,8> a , SIMD<int64_t,8> b)
  { return _mm512_cmp_epi64_mask (a.Data(), b.Data(), _MM_CMPINT_NE); }
  
  INLINE SIMD<mask64,8> operator&& (SIMD<mask64,8> a, SIMD<mask64,8> b)
  { return (__mmask8)(a.Data() & b.Data()); }
  INLINE SIMD<mask64,8> operator|| (SIMD<mask64,8> a, SIMD<mask64,8> b)
  { return (__mmask8)(a.Data() | b.Data()); }
  INLINE SIMD<mask64,8> operator! (SIMD<mask64,8> a)
  { return (__mmask8)(~a.Data()); }
  
  INLINE SIMD<double,8> If (SIMD<mask64,8> a, SIMD<double,8> b, SIMD<double,8> c)
  { return _mm512_mask_blend_pd(a.Data(), c.Data(), b.Data()); }
  
  INLINE SIMD<double,8> IfPos (SIMD<double,8> a, SIMD<double> b, SIMD<double> c)
  {
    auto k = _mm512_cmp_pd_mask(a.Data(),_mm512_setzero_pd(), _CMP_GT_OS);
    return _mm512_mask_blend_pd(k,c.Data(),b.Data());
  }
  INLINE SIMD<double,8> IfZero (SIMD<double,8> a, SIMD<double,8> b, SIMD<double,8> c)
  {
    auto k = _mm512_cmp_pd_mask(a.Data(),_mm512_setzero_pd(), _CMP_EQ_OS);
    return _mm512_mask_blend_pd(k,c.Data(),b.Data());
  }


   INLINE SIMD<int64_t,8> If (SIMD<mask64,8> a, SIMD<int64_t,8> b, SIMD<int64_t,8> c)
  { return _mm512_mask_blend_epi64(a.Data(), c.Data(), b.Data()); }

  INLINE auto Unpack (SIMD<double,8> a, SIMD<double,8> b)
  {
    return make_tuple(SIMD<double,8>(_mm512_unpacklo_pd(a.Data(),b.Data())),
                      SIMD<double,8>(_mm512_unpackhi_pd(a.Data(),b.Data())));
  }
  
   
  INLINE double HSum (SIMD<double,8> sd)
  {
    SIMD<double,4> low = _mm512_extractf64x4_pd(sd.Data(),0);
    SIMD<double,4> high = _mm512_extractf64x4_pd(sd.Data(),1);
    return HSum(low)+HSum(high);
  }

  INLINE auto HSum (SIMD<double,8> sd1, SIMD<double,8> sd2)
  {
    return SIMD<double,2>(HSum(sd1), HSum(sd2));
  }

  INLINE SIMD<double,4> HSum (SIMD<double,8> v1, SIMD<double,8> v2, SIMD<double,8> v3, SIMD<double,8> v4)
  {
    /*
    SIMD<double,4> high1 = _mm512_extractf64x4_pd(v1.Data(),1);
    SIMD<double,4> high2 = _mm512_extractf64x4_pd(v2.Data(),1);
    SIMD<double,4> high3 = _mm512_extractf64x4_pd(v3.Data(),1);
    SIMD<double,4> high4 = _mm512_extractf64x4_pd(v4.Data(),1);
    SIMD<double,4> low1 = _mm512_extractf64x4_pd(v1.Data(),0);
    SIMD<double,4> low2 = _mm512_extractf64x4_pd(v2.Data(),0);
    SIMD<double,4> low3 = _mm512_extractf64x4_pd(v3.Data(),0);
    SIMD<double,4> low4 = _mm512_extractf64x4_pd(v4.Data(),0);
    return HSum(low1,low2,low3,low4) + HSum(high1,high2,high3,high4);
    */
    
    SIMD<double> lo,hi;
    tie(lo,hi) = Unpack(v1, v2);
    SIMD<double> sum01 = lo+hi;
    tie(lo,hi) = Unpack(v3, v4);
    SIMD<double> sum23 = lo+hi;
    // sum01  b a b a b a b a
    // sum23  d c d c d c d c
    // __m512 perm = _mm512_permutex2var_pd (sum01.Data(), _mm512_set_epi64(1,2,3,4,5,6,7,8), sum23.Data());
    __m256d ab =  _mm512_extractf64x4_pd(sum01.Data(),0) + _mm512_extractf64x4_pd(sum01.Data(),1);
    __m256d cd =  _mm512_extractf64x4_pd(sum23.Data(),0) + _mm512_extractf64x4_pd(sum23.Data(),1);
    return _mm256_add_pd (_mm256_permute2f128_pd (ab, cd, 1+2*16), _mm256_blend_pd (ab, cd, 12));
  }
  
#endif



  


  /*  
  INLINE SIMD<double> operator+ (SIMD<double> a, SIMD<double> b) { return a.Data()+b.Data(); }
  INLINE SIMD<double> operator- (SIMD<double> a, SIMD<double> b) { return a.Data()-b.Data(); }
  INLINE SIMD<double> operator- (SIMD<double> a) { return -a.Data(); }
  INLINE SIMD<double> operator* (SIMD<double> a, SIMD<double> b) { return a.Data()*b.Data(); }
  INLINE SIMD<double> operator/ (SIMD<double> a, SIMD<double> b) { return a.Data()/b.Data(); }
  INLINE SIMD<double> operator* (double a, SIMD<double> b) { return SIMD<double>(a)*b; }
  INLINE SIMD<double> operator* (SIMD<double> b, double a) { return SIMD<double>(a)*b; }
  INLINE SIMD<double> operator+= (SIMD<double> & a, SIMD<double> b) { return a.Data()+=b.Data(); }
  INLINE SIMD<double> operator-= (SIMD<double> & a, SIMD<double> b) { return a.Data()-=b.Data(); }
  INLINE SIMD<double> operator*= (SIMD<double> & a, SIMD<double> b) { return a.Data()*=b.Data(); }
  INLINE SIMD<double> operator/= (SIMD<double> & a, SIMD<double> b) { return a.Data()/=b.Data(); }
  */
  
  INLINE SIMD<double,1> sqrt (SIMD<double,1> a) { return std::sqrt(a.Data()); }
  INLINE SIMD<double,1> floor (SIMD<double,1> a) { return std::floor(a.Data()); }
  INLINE SIMD<double,1> ceil (SIMD<double,1> a) { return std::ceil(a.Data()); }
  INLINE SIMD<double,1> fabs (SIMD<double,1> a) { return std::fabs(a.Data()); }
  INLINE SIMD<double,1> L2Norm2 (SIMD<double,1> a) { return a.Data()*a.Data(); }
  INLINE SIMD<double,1> Trans (SIMD<double,1> a) { return a; }
  INLINE SIMD<double,1> IfPos (SIMD<double,1> a, SIMD<double,1> b, SIMD<double,1> c)
  {
    return (a.Data() > 0) ? b : c;
  }
  INLINE SIMD<double,1> IfZero (SIMD<double,1> a, SIMD<double,1> b, SIMD<double,1> c)
  {
    return (a.Data() == 0.) ? b : c;
  }

  INLINE double HSum (SIMD<double,1> sd)
  { return sd.Data(); }
  INLINE auto HSum (SIMD<double,1> sd1, SIMD<double,1> sd2)
  { return SIMD<double,2>(sd1.Data(), sd2.Data()); }
  INLINE auto HSum (SIMD<double,1> sd1, SIMD<double,1> sd2, SIMD<double,1> sd3, SIMD<double,1> sd4)
  { return SIMD<double,4>(sd1.Data(), sd2.Data(), sd3.Data(), sd4.Data()); }


  

  INLINE void SIMDTranspose (SIMD<double,4> a1, SIMD<double,4> a2, SIMD <double,4> a3, SIMD<double,4> a4,
                             SIMD<double,4> & b1, SIMD<double,4> & b2, SIMD<double,4> & b3, SIMD<double,4> & b4)
  {
    SIMD<double,4> h1,h2,h3,h4;
    tie(h1,h2) = Unpack(a1,a2);
    tie(h3,h4) = Unpack(a3,a4);
    b1 = SIMD<double,4> (h1.Lo(), h3.Lo());
    b2 = SIMD<double,4> (h2.Lo(), h4.Lo());
    b3 = SIMD<double,4> (h1.Hi(), h3.Hi());
    b4 = SIMD<double,4> (h2.Hi(), h4.Hi());
  }


  template <int NUM, typename FUNC>
  INLINE void Iterate2 (FUNC f)
  {
    if constexpr (NUM > 1) Iterate2<NUM-1> (f);
    if constexpr (NUM >= 1) f(IC<NUM-1>());
  }

  
  template <typename T, int N>
  ostream & operator<< (ostream & ost, SIMD<T,N> simd)
  {
    /*
    ost << simd[0];
    for (int i = 1; i < simd.Size(); i++)
      ost << " " << simd[i];
    */
    Iterate2<simd.Size()> ([&] (auto I) {
        if (I.value != 0) ost << " ";
        ost << get<I.value>(simd);
      });
    return ost;
  }

  using std::exp;
  template <int N>
  INLINE ngstd::SIMD<double,N> exp (ngstd::SIMD<double,N> a) {
    return ngstd::SIMD<double>([a](int i)->double { return exp(a[i]); } );
  }

  using std::log;
  template <int N>  
  INLINE ngstd::SIMD<double,N> log (ngstd::SIMD<double,N> a) {
    return ngstd::SIMD<double,N>([a](int i)->double { return log(a[i]); } );
  }

  using std::pow;
  template <int N>    
  INLINE ngstd::SIMD<double,N> pow (ngstd::SIMD<double,N> a, double x) {
    return ngstd::SIMD<double,N>([a,x](int i)->double { return pow(a[i],x); } );
  }

  template <int N>
  INLINE ngstd::SIMD<double,N> pow (ngstd::SIMD<double,N> a, ngstd::SIMD<double,N> b) {
    return ngstd::SIMD<double,N>([a,b](int i)->double { return pow(a[i],b[i]); } );
  }

  using std::sin;
  template <int N>      
  INLINE ngstd::SIMD<double,N> sin (ngstd::SIMD<double,N> a) {
    return ngstd::SIMD<double,N>([a](int i)->double { return sin(a[i]); } );
  }
  
  using std::cos;
  template <int N>        
  INLINE ngstd::SIMD<double,N> cos (ngstd::SIMD<double,N> a) {
    return ngstd::SIMD<double,N>([a](int i)->double { return cos(a[i]); } );
  }

  using std::tan;
  template <int N>        
  INLINE ngstd::SIMD<double,N> tan (ngstd::SIMD<double,N> a) {
    return ngstd::SIMD<double,N>([a](int i)->double { return tan(a[i]); } );
  }

  using std::atan;
  template <int N>          
  INLINE ngstd::SIMD<double,N> atan (ngstd::SIMD<double,N> a) {
    return ngstd::SIMD<double,N>([a](int i)->double { return atan(a[i]); } );
  }

  using std::atan2;
  template <int N>          
  INLINE ngstd::SIMD<double,N> atan2 (ngstd::SIMD<double,N> y, ngstd::SIMD<double,N> x) {
    return ngstd::SIMD<double,N>([y,x](int i)->double { return atan2(y[i], x[i]); } );
  }

  using std::acos;
  template <int N>          
  INLINE ngstd::SIMD<double,N> acos (ngstd::SIMD<double,N> a) {
    return ngstd::SIMD<double,N>([a](int i)->double { return acos(a[i]); } );
  }

  using std::asin;
  template <int N>          
  INLINE ngstd::SIMD<double,N> asin (ngstd::SIMD<double,N> a) {
    return ngstd::SIMD<double,N>([a](int i)->double { return asin(a[i]); } );
  }

  using std::sinh;
  template <int N>      
  INLINE ngstd::SIMD<double,N> sinh (ngstd::SIMD<double,N> a) {
    return ngstd::SIMD<double,N>([a](int i)->double { return sinh(a[i]); } );
  }
  
  using std::cosh;
  template <int N>        
  INLINE ngstd::SIMD<double,N> cosh (ngstd::SIMD<double,N> a) {
    return ngstd::SIMD<double,N>([a](int i)->double { return cosh(a[i]); } );
  }
  

  template <int D, typename T>
  class MultiSIMD : public AlignedAlloc<MultiSIMD<D,T>>
  {
    SIMD<T> head;
    MultiSIMD<D-1,T> tail;
  public:
    MultiSIMD () = default;
    MultiSIMD (const MultiSIMD & ) = default;
    MultiSIMD (T v) : head(v), tail(v) { ; } 
    MultiSIMD (SIMD<T> _head, MultiSIMD<D-1,T> _tail)
      : head(_head), tail(_tail) { ; }
    template <typename ... ARGS>
    MultiSIMD (SIMD<T> _v0, SIMD<T> _v1, ARGS ... args)
      : head(_v0), tail(_v1, args...) { ; }
    SIMD<T> Head() const { return head; }
    MultiSIMD<D-1,T> Tail() const { return tail; }
    SIMD<T> & Head() { return head; }
    MultiSIMD<D-1,T> & Tail() { return tail; }

    template <int NR>
    SIMD<T> Get() const { return NR==0 ? head : tail.template Get<NR-1>(); }
    template <int NR>
    SIMD<T> & Get() { return NR==0 ? head : tail.template Get<NR-1>(); }
    auto MakeTuple () { return tuple_cat(tuple<SIMD<T>&> (head), tail.MakeTuple()); }
    // not yet possible for MSVC
    // operator auto () { return MakeTuple(); }
  };

  template <typename T>
  class MultiSIMD<2,T> : public AlignedAlloc<MultiSIMD<2,T>>
  {
    SIMD<T> v0, v1;
  public:
    MultiSIMD () = default;
    MultiSIMD (const MultiSIMD & ) = default;
    MultiSIMD (T v) : v0(v), v1(v) { ; } 
    MultiSIMD (SIMD<T> _v0, SIMD<T> _v1) : v0(_v0), v1(_v1) { ; }
    
    SIMD<T> Head() const { return v0; }
    SIMD<T> Tail() const { return v1; }
    SIMD<T> & Head() { return v0; }
    SIMD<T> & Tail() { return v1; } 

    template <int NR>
    SIMD<T> Get() const { return NR==0 ? v0 : v1; }
    template <int NR>
    SIMD<T> & Get() { return NR==0 ? v0 : v1; }

    auto MakeTuple () { return tuple<SIMD<T>&,SIMD<T>&>(v0, v1); }
    operator tuple<SIMD<T>&,SIMD<T>&> () { return MakeTuple(); }
  };

  template <int D> INLINE MultiSIMD<D,double> operator+ (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()+b.Head(), a.Tail()+b.Tail()); }
  template <int D> INLINE MultiSIMD<D,double> operator+ (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a+b.Head(), a+b.Tail()); }
  template <int D> INLINE MultiSIMD<D,double> operator+ (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> (a+b.Head(), a+b.Tail()); }
  
  template <int D> INLINE MultiSIMD<D,double> operator- (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()-b.Head(), a.Tail()-b.Tail()); }
  template <int D> INLINE MultiSIMD<D,double> operator- (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a-b.Head(), a-b.Tail()); }
  template <int D> INLINE MultiSIMD<D,double> operator- (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> (b.Head()-a, b.Tail()-a); }
  template <int D> INLINE MultiSIMD<D,double> operator- (MultiSIMD<D,double> a)
  { return MultiSIMD<D,double> (-a.Head(), -a.Tail()); }
  template <int D> INLINE MultiSIMD<D,double> operator* (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()*b.Head(), a.Tail()*b.Tail()); }
  template <int D> INLINE MultiSIMD<D,double> operator/ (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()/b.Head(), a.Tail()/b.Tail()); }
  template <int D> INLINE MultiSIMD<D,double> operator* (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> ( a*b.Head(), a*b.Tail()); }
  template <int D> INLINE MultiSIMD<D,double> operator* (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> ( a*b.Head(), a*b.Tail()); }  

  template <int D> INLINE MultiSIMD<D,double> & operator+= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b) 
  { a.Head()+=b.Head(); a.Tail()+=b.Tail(); return a; }
  template <int D> INLINE MultiSIMD<D,double> operator-= (MultiSIMD<D,double> & a, double b)
  { a.Head()-=b; a.Tail()-=b; return a; }
  template <int D> INLINE MultiSIMD<D,double> operator-= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b)
  { a.Head()-=b.Head(); a.Tail()-=b.Tail(); return a; }
  template <int D> INLINE MultiSIMD<D,double> & operator*= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b)
  { a.Head()*=b.Head(); a.Tail()*=b.Tail(); return a; }
  template <int D> INLINE MultiSIMD<D,double> & operator*= (MultiSIMD<D,double> & a, double b)
  { a.Head()*=b; a.Tail()*=b; return a; }
  // INLINE MultiSIMD<double> operator/= (MultiSIMD<double> & a, MultiSIMD<double> b) { return a.Data()/=b.Data(); }


  template <int D, typename T>
  ostream & operator<< (ostream & ost, MultiSIMD<D,T> multi)
  {
    ost << multi.Head() << " " << multi.Tail();
    return ost;
  }

  INLINE SIMD<double> HVSum (SIMD<double> a) { return a; }
  template <int D>
  INLINE SIMD<double> HVSum (MultiSIMD<D,double> a) { return a.Head() + HVSum(a.Tail()); }

  template <int D> INLINE double HSum (MultiSIMD<D,double> a) { return HSum(HVSum(a)); }
  template <int D> INLINE auto HSum (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return HSum(HVSum(a), HVSum(b)); }





  template <typename T1, typename T2, typename T3>
  // a*b+c
  INLINE auto FMA(T1 a, T2 b, T3 c)
  {
    return a*b+c;
  }

#ifdef __AVX512F__
  INLINE SIMD<double,8> FMA (SIMD<double,8> a, SIMD<double,8> b, SIMD<double,8> c)
  {
    return _mm512_fmadd_pd (a.Data(), b.Data(), c.Data());
  }
  INLINE SIMD<double,8> FMA (const double & a, SIMD<double,8> b, SIMD<double,8> c)
  {
    return _mm512_fmadd_pd (_mm512_set1_pd(a), b.Data(), c.Data());    
  }
#endif
#ifdef __FMA__
  INLINE SIMD<double,4> FMA (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> c)
  {
    return _mm256_fmadd_pd (a.Data(), b.Data(), c.Data());
  }
  INLINE SIMD<double,4> FMA (const double & a, SIMD<double,4> b, SIMD<double,4> c)
  {
    return _mm256_fmadd_pd (_mm256_set1_pd(a), b.Data(), c.Data());
  }
#endif

  // update form of fma
  template <int N>
  void FMAasm (SIMD<double,N> a, SIMD<double,N> b, SIMD<double,N> & sum)
  {
    sum = FMA(a,b,sum);
  }

#if defined(__FMA__) && !defined(__AVX512F__)
  // make sure to use the update-version of fma
  // important in matrix kernels using 12 sum-registers, 3 a-values and updated b-value
  // avx512 has enough registers, and gcc seems to use only the first 16 z-regs
  INLINE void FMAasm (SIMD<double,4> a, SIMD<double,4> b, SIMD<double,4> & sum)
  {
    asm ("vfmadd231pd %[a], %[b], %[sum]"
         : [sum] "+x" (sum.Data())
         : [a] "x" (a.Data()), [b] "x" (b.Data())
         );
  }
#endif


  
  
  template <int D>
  INLINE MultiSIMD<D,double> FMA(MultiSIMD<D,double> a, MultiSIMD<D,double> b, MultiSIMD<D,double> c)
  {
    return MultiSIMD<D,double> (FMA (a.Head(), b.Head(), c.Head()), FMA (a.Tail(), b.Tail(), c.Tail()));
  }
  
  template <int D>
  INLINE MultiSIMD<D,double> FMA(const double & a, MultiSIMD<D,double> b, MultiSIMD<D,double> c)
  {
    return MultiSIMD<D,double> (FMA (a, b.Head(), c.Head()), FMA (a, b.Tail(), c.Tail()));
  }


#if defined __AVX512F__

  INLINE SIMD<mask64, 8> GetMaskFromBits (unsigned int i)
  {
    return SIMD<mask64, 8>(__mmask8(i));
  }

#elif defined __AVX__

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

  INLINE SIMD<mask64, 4> GetMaskFromBits (unsigned int i)
  {
    return masks_from_4bits[i & 15];
  }

#elif defined __SSE__

  static SIMD<mask64, 2> masks_from_2bits[4] = {
    _mm_set_epi32 (0,0,0,0), _mm_set_epi32 (0,0,-1,0),
    _mm_set_epi32 (-1,0,0,0), _mm_set_epi32 (-1,0,-1,0),
  };

  INLINE SIMD<mask64, 2> GetMaskFromBits (unsigned int i)
  {
    return masks_from_2bits[i & 3];
  }

#endif

  template <int i, typename T, int N>
  T get(SIMD<T,N> a) { return a[i]; }
  
}

namespace std
{
  // structured binding support
  template <typename T, int N >
  struct tuple_size<ngstd::SIMD<T,N>> : std::integral_constant<std::size_t, N> {};
  template<size_t N, typename T, int M> struct tuple_element<N,ngstd::SIMD<T,M>> { using type = T; };
}



#endif
