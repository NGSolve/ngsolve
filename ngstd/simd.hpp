#ifndef FILE_SIMD
#define FILE_SIMD

/**************************************************************************/
/* File:   simd.hpp                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

#include <immintrin.h>

#ifdef WIN32
inline __m128d operator- (__m128d a) { return _mm_xor_pd(a, _mm_set1_pd(-0.0)); }
inline __m128d operator+ (__m128d a, __m128d b) { return _mm_add_pd(a,b); }
inline __m128d operator- (__m128d a, __m128d b) { return _mm_sub_pd(a,b); }
inline __m128d operator* (__m128d a, __m128d b) { return _mm_mul_pd(a,b); }
inline __m128d operator/ (__m128d a, __m128d b) { return _mm_div_pd(a,b); }
inline __m128d operator* (double a, __m128d b) { return _mm_set1_pd(a)*b; }
inline __m128d operator* (__m128d b, double a) { return _mm_set1_pd(a)*b; }

inline __m128d operator+= (__m128d &a, __m128d b) { return a = a+b; }
inline __m128d operator-= (__m128d &a, __m128d b) { return a = a-b; }
inline __m128d operator*= (__m128d &a, __m128d b) { return a = a*b; }
inline __m128d operator/= (__m128d &a, __m128d b) { return a = a/b; }

inline __m256d operator- (__m256d a) { return _mm256_xor_pd(a, _mm256_set1_pd(-0.0)); }
inline __m256d operator+ (__m256d a, __m256d b) { return _mm256_add_pd(a,b); }
inline __m256d operator- (__m256d a, __m256d b) { return _mm256_sub_pd(a,b); }
inline __m256d operator* (__m256d a, __m256d b) { return _mm256_mul_pd(a,b); }
inline __m256d operator/ (__m256d a, __m256d b) { return _mm256_div_pd(a,b); }
inline __m256d operator* (double a, __m256d b) { return _mm256_set1_pd(a)*b; }
inline __m256d operator* (__m256d b, double a) { return _mm256_set1_pd(a)*b; }

inline __m256d operator+= (__m256d &a, __m256d b) { return a = a+b; }
inline __m256d operator-= (__m256d &a, __m256d b) { return a = a-b; }
inline __m256d operator*= (__m256d &a, __m256d b) { return a = a*b; }
inline __m256d operator/= (__m256d &a, __m256d b) { return a = a/b; }
#endif

namespace ngstd
{
  template <typename T> class SIMD;

  template <typename T>
  struct has_call_operator
  {
      template <typename C> static std::true_type check( decltype( sizeof(&C::operator() )) ) { return std::true_type(); }
      template <typename> static std::false_type check(...) { return std::false_type(); }
      typedef decltype( check<T>(sizeof(char)) ) type;
      static constexpr type value = type();
  };

#ifdef __AVX2__
  
  template<>
  class alignas(32) SIMD<double>
  {
    __m256d data;
    
  public:
    static constexpr int Size() { return 4; }
    SIMD () = default;
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (double val)
    {
      data = _mm256_set1_pd(val);
    }
    
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
    
    template <typename Function>
    void SIMD_function (const Function & func, std::true_type)
    {
      data = _mm256_set_pd(func(3), func(2), func(1), func(0));
    }
    
    // not a function
    void SIMD_function (double const * p, std::false_type)
    {
      data = _mm256_loadu_pd(p);
    }
    
    void SIMD_function (double val, std::false_type)
    {
      data = _mm256_set1_pd(val);
    }
    
    void SIMD_function (__m256d _data, std::false_type)
    {
      data = _data;
    }
    
    double operator[] (int i) const { return ((double*)(&data))[i]; }
    __m256d Data() const { return data; }
    __m256d & Data() { return data; }
  };
  
  
  inline SIMD<double> operator+ (SIMD<double> a, SIMD<double> b) { return a.Data()+b.Data(); }
  inline SIMD<double> operator- (SIMD<double> a, SIMD<double> b) { return a.Data()-b.Data(); }
  inline SIMD<double> operator- (SIMD<double> a) { return -a.Data(); }
  inline SIMD<double> operator* (SIMD<double> a, SIMD<double> b) { return a.Data()*b.Data(); }
  inline SIMD<double> operator/ (SIMD<double> a, SIMD<double> b) { return a.Data()/b.Data(); }
  inline SIMD<double> operator* (double a, SIMD<double> b) { return SIMD<double>(a)*b; }
  inline SIMD<double> operator* (SIMD<double> b, double a) { return SIMD<double>(a)*b; }
  inline SIMD<double> operator+= (SIMD<double> & a, SIMD<double> b) { return a.Data()+=b.Data(); }
  inline SIMD<double> operator-= (SIMD<double> & a, SIMD<double> b) { return a.Data()-=b.Data(); }
  inline SIMD<double> operator*= (SIMD<double> & a, SIMD<double> b) { return a.Data()*=b.Data(); }
  inline SIMD<double> operator/= (SIMD<double> & a, SIMD<double> b) { return a.Data()/=b.Data(); }

  inline SIMD<double> sqrt (SIMD<double> a) { return _mm256_sqrt_pd(a.Data()); }
  inline SIMD<double> fabs (SIMD<double> a) { return _mm256_max_pd(a.Data(), -a.Data()); }
  inline SIMD<double> L2Norm2 (SIMD<double> a) { return a.Data()*a.Data(); }
  inline SIMD<double> Trans (SIMD<double> a) { return a; }
  inline SIMD<double> IfPos (SIMD<double> a, SIMD<double> b, SIMD<double> c)
  {
    auto cp = _mm256_cmp_pd (a.Data(), _mm256_setzero_pd(), _CMP_GT_OS);
    return _mm256_blendv_pd(c.Data(), b.Data(), cp);
  }

  inline double HSum (SIMD<double> sd)
  {
    __m128d hv = _mm_add_pd (_mm256_extractf128_pd(sd.Data(),0), _mm256_extractf128_pd(sd.Data(),1));
    return _mm_cvtsd_f64 (_mm_hadd_pd (hv, hv));
  }
  
  


#else

  template<>
  class SIMD<double>
  {
    double data;
    
  public:
    static constexpr int Size() { return 1; }
    SIMD () = default;
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;
    
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
    
    template <typename Function>
    void SIMD_function (const Function & func, std::true_type)
    {
      data = func(0);
    }
    
    // not a function
    void SIMD_function (double const * p, std::false_type)
    {
      data = *p;
    }
    
    void SIMD_function (double val, std::false_type)
    {
      data = val;
    }
    
    double operator[] (int i) const { return ((double*)(&data))[i]; }
    double Data() const { return data; }
    double & Data() { return data; }
  };
  
  
  inline SIMD<double> operator+ (SIMD<double> a, SIMD<double> b) { return a.Data()+b.Data(); }
  inline SIMD<double> operator- (SIMD<double> a, SIMD<double> b) { return a.Data()-b.Data(); }
  inline SIMD<double> operator- (SIMD<double> a) { return -a.Data(); }
  inline SIMD<double> operator* (SIMD<double> a, SIMD<double> b) { return a.Data()*b.Data(); }
  inline SIMD<double> operator/ (SIMD<double> a, SIMD<double> b) { return a.Data()/b.Data(); }
  inline SIMD<double> operator* (double a, SIMD<double> b) { return SIMD<double>(a)*b; }
  inline SIMD<double> operator* (SIMD<double> b, double a) { return SIMD<double>(a)*b; }
  inline SIMD<double> operator+= (SIMD<double> & a, SIMD<double> b) { return a.Data()+=b.Data(); }
  inline SIMD<double> operator-= (SIMD<double> & a, SIMD<double> b) { return a.Data()-=b.Data(); }
  inline SIMD<double> operator*= (SIMD<double> & a, SIMD<double> b) { return a.Data()*=b.Data(); }
  inline SIMD<double> operator/= (SIMD<double> & a, SIMD<double> b) { return a.Data()/=b.Data(); }

  inline SIMD<double> sqrt (SIMD<double> a) { return std::sqrt(a.Data()); }
  inline SIMD<double> fabs (SIMD<double> a) { return std::fabs(a.Data()); }
  inline SIMD<double> L2Norm2 (SIMD<double> a) { return a.Data()*a.Data(); }
  inline SIMD<double> Trans (SIMD<double> a) { return a; }
  inline SIMD<double> IfPos (SIMD<double> a, SIMD<double> b, SIMD<double> c)
  {
    return (a.Data() > 0) ? b : c;
  }

  inline double HSum (SIMD<double> sd)
  { return sd.Data(); }

#endif





  
  
  template <typename T>
  ostream & operator<< (ostream & ost, SIMD<T> simd)
  {
    ost << simd[0];
    for (int i = 1; i < simd.Size(); i++)
      ost << " " << simd[i];
    return ost;
  }
}


#endif
