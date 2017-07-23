#ifndef FILE_SIMD
#define FILE_SIMD

/**************************************************************************/
/* File:   simd.hpp                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 16                                                    */
/**************************************************************************/

#include <immintrin.h>

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
#endif
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

#ifdef __AVX__

#ifdef __AVX512F__
  template<>
  class alignas(64) SIMD<double> 
  {
    __m512d data;
    
  public:
    static constexpr int Size() { return 8; }
    SIMD () = default;
    SIMD (const SIMD &) = default;
    SIMD & operator= (const SIMD &) = default;

    SIMD (double val)
    {
      data = _mm512_set1_pd(val);
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
      /*
      data = _mm512_set_pd(func(7), func(6), func(5), func(4),
                           func(3), func(2), func(1), func(0));
      */
      data = (__m512){ func(7), func(6), func(5), func(4),
                       func(3), func(2), func(1), func(0));
                       
      
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
    INLINE __m512d Data() const { return data; }
    INLINE __m512d & Data() { return data; }
  };
 
#else

  /*
  template <size_t ALIGN = 64>
  class AlignedAlloc
  {
  public:
    void * operator new (size_t s) { return  _mm_malloc(s, ALIGN); }
    void * operator new[] (size_t s) { return  _mm_malloc(s, ALIGN); }
    void operator delete (void * p) { _mm_free(p); }
    void operator delete[] (void * p) { _mm_free(p); }
  };
  */

  template <typename T>
  class AlignedAlloc
  {
    protected:
      static void * aligned_malloc(size_t s)
      {
        // Assume 16 byte alignment of standard library
        if(alignof(T)<=16)
            return malloc(s);
        else
            return  _mm_malloc(s, alignof(T));
      }

      static void aligned_free(void *p)
      {
        if(alignof(T)<=16)
            free(p);
        else
            _mm_free(p);
      }

  public:
    void * operator new (size_t s, void *p) { return p; }
    void * operator new (size_t s) { return aligned_malloc(s); }
    void * operator new[] (size_t s) { return aligned_malloc(s); }
    void operator delete (void * p) { aligned_free(p); }
    void operator delete[] (void * p) { aligned_free(p); }
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

    SIMD (double val) { data = _mm256_set1_pd(val); }
    SIMD (int val)    { data = _mm256_set1_pd(val); }
    SIMD (size_t val) { data = _mm256_set1_pd(val); }

    SIMD (double const * p) { data = _mm256_loadu_pd(p); }
    SIMD (__m256d _data) { data = _data; }
    
    template<typename T, typename std::enable_if<std::is_convertible<T, std::function<double(int)>>::value, int>::type = 0>                                                                    SIMD (const T & func)
    {   
      data = _mm256_set_pd(func(3), func(2), func(1), func(0));              
    }   
    
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
    
    /*
    void * operator new (size_t s) { return  _mm_malloc(s, 64); }
    void * operator new[] (size_t s) { return  _mm_malloc(s, 64); }
    void operator delete (void * p) { _mm_free(p); }
    void operator delete[] (void * p) { _mm_free(p); }
    */

    /*
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
    */
    
    INLINE double operator[] (int i) const { return ((double*)(&data))[i]; }
    INLINE __m256d Data() const { return data; }
    INLINE __m256d & Data() { return data; }
  };
#endif

  
  
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

  INLINE SIMD<double> L2Norm2 (SIMD<double> a) { return a.Data()*a.Data(); }
  INLINE SIMD<double> Trans (SIMD<double> a) { return a; }

#ifdef __AVX512F__
  INLINE SIMD<double> sqrt (SIMD<double> a) { return _mm512_sqrt_pd(a.Data()); }
  INLINE SIMD<double> fabs (SIMD<double> a) { return _mm512_max_pd(a.Data(), -a.Data()); }
  INLINE SIMD<double> IfPos (SIMD<double> a, SIMD<double> b, SIMD<double> c)
  {
    /*
    auto cp = _mm512_cmp_pd (a.Data(), _mm512_setzero_pd(), _CMP_GT_OS);
    return _mm512_blendv_pd(c.Data(), b.Data(), cp);
    */
    throw Exception ("IfPos missing for AVX512");
  }

  INLINE double HSum (SIMD<double> sd)
  {
    throw Exception ("HSum missing for AVX512");    
    // __m128d hv = _mm_add_pd (_mm256_extractf128_pd(sd.Data(),0), _mm256_extractf128_pd(sd.Data(),1));
    // return _mm_cvtsd_f64 (_mm_hadd_pd (hv, hv));
  }

#else
  INLINE SIMD<double> sqrt (SIMD<double> a) { return _mm256_sqrt_pd(a.Data()); }
  INLINE SIMD<double> fabs (SIMD<double> a) { return _mm256_max_pd(a.Data(), -a.Data()); }
  INLINE SIMD<double> IfPos (SIMD<double> a, SIMD<double> b, SIMD<double> c)
  {
    auto cp = _mm256_cmp_pd (a.Data(), _mm256_setzero_pd(), _CMP_GT_OS);
    return _mm256_blendv_pd(c.Data(), b.Data(), cp);
  }

  INLINE double HSum (SIMD<double> sd)
  {
    __m128d hv = _mm_add_pd (_mm256_extractf128_pd(sd.Data(),0), _mm256_extractf128_pd(sd.Data(),1));
    return _mm_cvtsd_f64 (_mm_hadd_pd (hv, hv));
  }

  INLINE auto HSum (SIMD<double> sd1, SIMD<double> sd2)
  {
    __m256d hv = _mm256_hadd_pd(sd1.Data(), sd2.Data());
    __m128d hv2 = _mm_add_pd (_mm256_extractf128_pd(hv,0), _mm256_extractf128_pd(hv,1));
    return make_tuple(_mm_cvtsd_f64 (hv2),  _mm_cvtsd_f64(_mm_shuffle_pd (hv2, hv2, 3)));
  }

  INLINE auto HSum (SIMD<double> v1, SIMD<double> v2, SIMD<double> v3, SIMD<double> v4)
  {
    __m256d hsum1 = _mm256_hadd_pd (v1.Data(), v2.Data());
    __m256d hsum2 = _mm256_hadd_pd (v3.Data(), v4.Data());
    __m256d hsum = _mm256_add_pd (_mm256_permute2f128_pd (hsum1, hsum2, 1+2*16),
                                  _mm256_blend_pd (hsum1, hsum2, 12));
    // return hsum;
    return make_tuple(hsum[0], hsum[1], hsum[2], hsum[3]);
  }
  
#endif  
  


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
    SIMD (double val) { data = val; }
    SIMD (int val)    { data = val; }
    SIMD (size_t val) { data = val; }
    SIMD (double const * p) { data = *p; }
    
    
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
    
    double operator[] (int i) const { return ((double*)(&data))[i]; }
    double Data() const { return data; }
    double & Data() { return data; }
  };
  
  
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

  INLINE SIMD<double> sqrt (SIMD<double> a) { return std::sqrt(a.Data()); }
  INLINE SIMD<double> fabs (SIMD<double> a) { return std::fabs(a.Data()); }
  INLINE SIMD<double> L2Norm2 (SIMD<double> a) { return a.Data()*a.Data(); }
  INLINE SIMD<double> Trans (SIMD<double> a) { return a; }
  INLINE SIMD<double> IfPos (SIMD<double> a, SIMD<double> b, SIMD<double> c)
  {
    return (a.Data() > 0) ? b : c;
  }

  INLINE double HSum (SIMD<double> sd)
  { return sd.Data(); }
  INLINE auto HSum (SIMD<double> sd1, SIMD<double> sd2)
  { return make_tuple(sd1.Data(), sd2.Data()); }
  INLINE auto HSum (SIMD<double> sd1, SIMD<double> sd2, SIMD<double> sd3, SIMD<double> sd4)
  { return make_tuple(sd1.Data(), sd2.Data(), sd3.Data(), sd4.Data()); }
#endif





  
  
  template <typename T>
  ostream & operator<< (ostream & ost, SIMD<T> simd)
  {
    ost << simd[0];
    for (int i = 1; i < simd.Size(); i++)
      ost << " " << simd[i];
    return ost;
  }

  using std::exp;
INLINE ngstd::SIMD<double> exp (ngstd::SIMD<double> a) {
  return ngstd::SIMD<double>([&](int i)->double { return exp(a[i]); } );
}

  using std::log;
INLINE ngstd::SIMD<double> log (ngstd::SIMD<double> a) {
  return ngstd::SIMD<double>([&](int i)->double { return log(a[i]); } );
}

  using std::pow;
INLINE ngstd::SIMD<double> pow (ngstd::SIMD<double> a, double x) {
  return ngstd::SIMD<double>([&](int i)->double { return pow(a[i],x); } );
}

  using std::sin;
INLINE ngstd::SIMD<double> sin (ngstd::SIMD<double> a) {
  return ngstd::SIMD<double>([&](int i)->double { return sin(a[i]); } );
}
  
  using std::cos;
INLINE ngstd::SIMD<double> cos (ngstd::SIMD<double> a) {
  return ngstd::SIMD<double>([&](int i)->double { return cos(a[i]); } );
}

  using std::tan;
INLINE ngstd::SIMD<double> tan (ngstd::SIMD<double> a) {
  return ngstd::SIMD<double>([&](int i)->double { return tan(a[i]); } );
}

  using std::atan;
INLINE ngstd::SIMD<double> atan (ngstd::SIMD<double> a) {
  return ngstd::SIMD<double>([&](int i)->double { return atan(a[i]); } );
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
    operator auto () { return MakeTuple(); }
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
    operator auto () { return MakeTuple(); }
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
  INLINE SIMD<double> FMA (SIMD<double> a, SIMD<double> b, SIMD<double> c)
  {
    return _mm512_fmadd_pd (a.Data(), b.Data(), c.Data());
  }
  INLINE SIMD<double> FMA (const double & a, SIMD<double> b, SIMD<double> c)
  {
    return _mm512_fmadd_pd (_mm256_set1_pd(a), b.Data(), c.Data());    
  }
#else
#ifdef __AVX2__
  INLINE SIMD<double> FMA (SIMD<double> a, SIMD<double> b, SIMD<double> c)
  {
    return _mm256_fmadd_pd (a.Data(), b.Data(), c.Data());
  }
  INLINE SIMD<double> FMA (const double & a, SIMD<double> b, SIMD<double> c)
  {
    return _mm256_fmadd_pd (_mm256_set1_pd(a), b.Data(), c.Data());
  }
#endif
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


#ifdef __AVX__
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
  
  class ExceptionNOSIMD : public Exception
  {
  public:
    using Exception :: Exception;
  };
  
}

#endif
