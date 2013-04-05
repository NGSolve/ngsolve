#ifndef FILE_NGBLA_MD
#define FILE_NGBLA_MD

/**************************************************************************/
/* File:   md.hpp                                                         */
/* Author: Joachim Schoeberl                                              */
/* Date:   03. Apr. 13                                                    */
/**************************************************************************/


/*
  Multiple-Data 
  implemeted by SSE, AVX, ....

*/



#if defined(__SSE3__) || defined (__AVX__)
      #ifndef WIN32
            #include <immintrin.h>
      #else
            #include <intrin.h>
      #endif
#endif


namespace ngbla
{

#ifdef __AVX__
  template <int D = 4> class MD;
#else
#ifdef __SSE3__
  template <int D = 2> class MD;
#else
  template <int D = 1> class MD;
#endif
#endif



#ifdef __AVX__

  template <> class MD<4>
  {
    __m256d data;

  public:
    enum { SIZE = 2 };
    MD () { ; }
    MD (__m256d at2) : data(at2) { ; }
    MD (double d1) { data = _mm256_set1_pd(d1); }

    MD & operator*= (const MD & v2)
    { data *= v2.data; return *this; }

    MD & operator+= (const MD & v2)
    { data += v2.data; return *this; }

    MD & operator-= (const MD & v2)
    { data -= v2.data; return *this; }
    
    __m256d Data() const { return data; }
    double operator[] (int i) const { return ((const double*)&data)[i]; }
    double & operator[] (int i) { return ((double*)&data)[i]; }
  };
#endif


#ifdef __SSE3__
  

#if defined(__INTEL_COMPILER) || defined(WIN32)
  inline __m128d operator+ (__m128d a, __m128d b) { return _mm_add_pd (a, b); }
  inline __m128d operator- (__m128d a, __m128d b) { return _mm_sub_pd (a, b); }
  inline __m128d operator* (__m128d a, __m128d b) { return _mm_mul_pd (a, b); }
  inline __m128d operator/ (__m128d a, __m128d b) { return _mm_div_pd (a, b); }
  inline __m128d& operator*= (__m128d& a, __m128d b) { return a=_mm_mul_pd (a, b); }
  inline __m128d& operator+= (__m128d& a, __m128d b) { return a=_mm_add_pd (a, b); }
  inline __m128d& operator-= (__m128d& a, __m128d b) { return a=_mm_sub_pd (a, b); }
#endif




  template <> class MD<2>
  {
    __m128d data;

  public:
    enum { SIZE = 2 };
    MD () { ; }
    MD (__m128d at2) : data(at2) { ; }
    // MD (int d1) { data = _mm_set1_pd(d1); }
    MD (double d1) { data = _mm_set1_pd(d1); }
    // MD (double d1, double d2) { data = _mm_set_pd(d2, d1); } 

    MD & operator*= (const MD & v2)
    { /*data = _mm_mul_pd (data, v2.data);*/data*=v2.data; return *this; }

    MD & operator+= (const MD & v2)
    { /*data = _mm_add_pd (data, v2.data);*/data+=v2.data; return *this; }

    MD & operator-= (const MD & v2)
    { /*data = _mm_sub_pd (data, v2.data);*/ data-=v2.data; return *this; }
    
    __m128d Data() const { return data; }
    double operator[] (int i) const { return ((const double*)&data)[i]; }
    double & operator[] (int i) { return ((double*)&data)[i]; }
  };
#endif

  template <> class MD<1>
  {
    double data;

  public:
    enum { SIZE = 1 };
    MD () { ; }
    MD (int d1) { data = d1; }
    MD (double d1) { data = d1; }

    MD & operator*= (const MD & v2)
    { data *= v2.data; return *this; }

    MD & operator+= (const MD & v2)
    { data += v2.data; return *this; }

    MD & operator-= (const MD & v2)
    { data -= v2.data; return *this; }

    double Data() const { return data; }
    double operator[] (int i) const { return ((const double*)&data)[i]; }
    double & operator[] (int i) { return ((double*)&data)[i]; }
  };







  template <int D>
  inline MD<D> operator+ (const MD<D> & a, const MD<D> & b)
  {
    return a.Data()+b.Data();
  }

  template <int D>
  inline MD<D> operator+ (double a, const MD<D> & b)
  {
    return MD<D>(a)+b;
  }

  template <int D>
  inline MD<D> operator+ (const MD<D> & a, double b)
  {
    return a+MD<D>(b);
  }




  template <int D>
  inline MD<D> operator- (const MD<D> & a, const MD<D> & b)
  {
    return a.Data()-b.Data();
  }

  template <int D>
  inline MD<D> operator- (const double & a, const MD<D> & b)
  {
    return MD<D>(a)-b;
  }

  template <int D>
  inline MD<D> operator- (const MD<D> & b, const double & a)
  {
    return b-MD<D>(a);
  }




  template <int D>
  inline MD<D> operator* (const MD<D> & a, const MD<D> & b)
  {
    return a.Data()*b.Data();
  }

  template <int D>
  inline MD<D> operator/ (const MD<D> & a, const MD<D> & b)
  {
    return a.Data()/b.Data();
  }

  template <int D>
  inline MD<D> operator* (const MD<D> & a, double b)
  {
    return b * a.Data();
  }

  template <int D>
  inline MD<D> operator* (double b, const MD<D> & a)
  {
    return MD<D>(b) * a;
  }
  
  template <int D>
  ALWAYS_INLINE inline double Sum (MD<D> a)
  {
    double sum = a[0];
    for (int j = 1; j < D; j++)
      sum += a[j];
    return sum;
  }

  ALWAYS_INLINE inline double Sum (MD<1> a)
  {
    return a[0];
  }



  template <int D>
  ALWAYS_INLINE inline double InnerProduct (const MD<D> & a, const MD<D> & b)
  {
    return Sum(a*b);
  }


  template <int D>
  inline ostream & operator<< (ostream & ost, const MD<D> & data)
  {
    ost << data[0];
    for (int j = 1; j < D; j++)
      ost << ", " << data[j];
    // ost << data[0] << ", " << data[1];
    return ost;
  }

  
  
  

  template <int D>
  class mat_traits<MD<D> >
  {
  public:
    typedef MD<D> TSCAL;
    typedef double TV_COL;
    typedef double TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = 0 };
  };


  template <int D>
  class mat_scale_type<MD<D>,double> 
  { 
  public: 
    typedef MD<D> TMAT; 
    typedef double TSCAL;
  };

  template <int D>
  class mat_prod_type<MD<D>,double> {public: typedef MD<D> TMAT;};
  template <int D>
  class mat_prod_type<double, MD<D> > {public: typedef MD<D> TMAT;};

  
}
#endif
