#ifndef FILE_MYSIMD
#define FILE_MYSIMD

/**************************************************************************/
/* File:   mysimd.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

#include <immintrin.h>

#ifdef WIN32 
#ifndef AVX_OPERATORS_DEFINED
#define AVX_OPERATORS_DEFINED
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
#endif



namespace netgen
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

#ifdef __AVX__

  template <typename T>
  class AlignedAlloc
  {
  public:
    void * operator new (size_t s, void *p) { return p; }
    void * operator new (size_t s) { return  _mm_malloc(s, alignof(T)); }
    void * operator new[] (size_t s) { return  _mm_malloc(s, alignof(T)); }
    void operator delete (void * p) { _mm_free(p); }
    void operator delete[] (void * p) { _mm_free(p); }
  };
  
  template<>
  class alignas(32) SIMD<double> : public AlignedAlloc<SIMD<double>>
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

    SIMD (__m256d adata)
      : data(adata)
      { ; } 

    /*
    template <typename T>
    SIMD (const T & val)
    {
//       SIMD_function(val, std::is_convertible<T, std::function<double(int)>>());
      SIMD_function(val, has_call_operator<T>::value);
    }
    */

    /*
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
    
    inline double operator[] (int i) const { return ((double*)(&data))[i]; }
    inline __m256d Data() const { return data; }
    inline __m256d & Data() { return data; }
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

  using std::sqrt;
  using std::fabs;
  
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

  // it's only a dummy without AVX
  template <typename T>
  class AlignedAlloc { ; };
  
  template<>
  class SIMD<double>
  {
    double data;
    
  public:
    static constexpr int Size() { return 1; }
    SIMD () = default;
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (double val)
      : data(val) { ; } 

    /*
    template <typename T>
    SIMD (const T & val)
    {
//       SIMD_function(val, std::is_convertible<T, std::function<double(int)>>());
      SIMD_function(val, has_call_operator<T>::value);
    }
    */

    
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

  using std::sqrt;
  using std::fabs;
  
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

  /*
  using std::exp;
  inline netgen::SIMD<double> exp (netgen::SIMD<double> a) {
    return netgen::SIMD<double>([&](int i)->double { return exp(a[i]); } );
}

  using std::log;
inline netgen::SIMD<double> log (netgen::SIMD<double> a) {
  return netgen::SIMD<double>([&](int i)->double { return log(a[i]); } );
}

  using std::pow;
inline netgen::SIMD<double> pow (netgen::SIMD<double> a, double x) {
  return netgen::SIMD<double>([&](int i)->double { return pow(a[i],x); } );
}
  */


  template <int D, typename T>
  class MultiSIMD
  {
    SIMD<T> head;
    MultiSIMD<D-1,T> tail;
  public:
    MultiSIMD () = default;
    MultiSIMD (const MultiSIMD & ) = default;
    MultiSIMD (T v) : head(v), tail(v) { ; } 
    MultiSIMD (SIMD<T> _head, MultiSIMD<D-1,T> _tail)
      : head(_head), tail(_tail) { ; }
    SIMD<T> Head() const { return head; }
    MultiSIMD<D-1,T> Tail() const { return tail; }
    SIMD<T> & Head() { return head; }
    MultiSIMD<D-1,T> & Tail() { return tail; }

    template <int NR>
    SIMD<T> Get() const { return NR==0 ? head : tail.template Get<NR-1>(); }
    template <int NR>
    SIMD<T> & Get() { return NR==0 ? head : tail.template Get<NR-1>(); }
  };

  template <typename T>
  class MultiSIMD<2,T>
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
  };

  template <int D> inline MultiSIMD<D,double> operator+ (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()+b.Head(), a.Tail()+b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator+ (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a+b.Head(), a+b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator+ (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> (a+b.Head(), a+b.Tail()); }
  
  template <int D> inline MultiSIMD<D,double> operator- (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()-b.Head(), a.Tail()-b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator- (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a-b.Head(), a-b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator- (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> (b.Head()-a, b.Tail()-a); }
  template <int D> inline MultiSIMD<D,double> operator- (MultiSIMD<D,double> a)
  { return MultiSIMD<D,double> (-a.Head(), -a.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator* (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()*b.Head(), a.Tail()*b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator/ (MultiSIMD<D,double> a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> (a.Head()/b.Head(), a.Tail()/b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator* (double a, MultiSIMD<D,double> b)
  { return MultiSIMD<D,double> ( a*b.Head(), a*b.Tail()); }
  template <int D> inline MultiSIMD<D,double> operator* (MultiSIMD<D,double> b, double a)
  { return MultiSIMD<D,double> ( a*b.Head(), a*b.Tail()); }  

  template <int D> inline MultiSIMD<D,double> & operator+= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b) 
  { a.Head()+=b.Head(); a.Tail()+=b.Tail(); return a; }
  template <int D> inline MultiSIMD<D,double> operator-= (MultiSIMD<D,double> & a, double b)
  { a.Head()-=b; a.Tail()-=b; return a; }
  template <int D> inline MultiSIMD<D,double> operator-= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b)
  { a.Head()-=b.Head(); a.Tail()-=b.Tail(); return a; }
  template <int D> inline MultiSIMD<D,double> & operator*= (MultiSIMD<D,double> & a, MultiSIMD<D,double> b)
  { a.Head()*=b.Head(); a.Tail()*=b.Tail(); return a; }
  template <int D> inline MultiSIMD<D,double> & operator*= (MultiSIMD<D,double> & a, double b)
  { a.Head()*=b; a.Tail()*=b; return a; }
  // inline MultiSIMD<double> operator/= (MultiSIMD<double> & a, MultiSIMD<double> b) { return a.Data()/=b.Data(); }


  template <int D, typename T>
  ostream & operator<< (ostream & ost, MultiSIMD<D,T> multi)
  {
    ost << multi.Head() << " " << multi.Tail();
    return ost;
  }

}

#endif
