#include <bla.hpp>



namespace ngstd
{

#if defined(__AVX__)

#if defined __INTEL_COMPILER
  typedef const char * prefetch_ptr_t;
#elif defined WIN32
  // MSVC needs a char*
  typedef char * prefetch_ptr_t;
#else
  typedef void * prefetch_ptr_t;
#endif

  template<int HINT> INLINE void Prefetch (void * p) { _mm_prefetch (reinterpret_cast<prefetch_ptr_t>(p), _MM_HINT_T2); }
  template <> INLINE void Prefetch<0> (void * p) { _mm_prefetch (reinterpret_cast<prefetch_ptr_t>(p), _MM_HINT_T0); }
  template <> INLINE void Prefetch<1> (void * p) { _mm_prefetch (reinterpret_cast<prefetch_ptr_t>(p), _MM_HINT_T1); }
  
#endif

}
  

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
#endif // AVX_OPERATORS_DEFINED
#endif // WIN32



namespace ngbla
{


#if defined (__AVX__) && !defined(__AVX512F__)


  inline SIMD<double,4> operator+= (SIMD<double,4> & a, __m256d b) { return a += SIMD<double,4>(b); }
  inline SIMD<double,4> operator-= (SIMD<double,4> & a, __m256d b) { return a -= SIMD<double,4>(b); }  

  
  INLINE __m256d HAdd (__m256d v1, __m256d v2, __m256d v3, __m256d v4)
  {
    __m256d hsum1 = _mm256_hadd_pd (v1, v2);
    __m256d hsum2 = _mm256_hadd_pd (v3, v4);
    __m256d hsum = _mm256_add_pd (_mm256_permute2f128_pd (hsum1, hsum2, 1+2*16),
                                  _mm256_blend_pd (hsum1, hsum2, 12));
    return hsum;
  }


  
#include "matkernel.hpp"

  INLINE auto Transpose (__m256d a, __m256d b, __m256d c, __m256d d)
  {
    __m256d loab = _mm256_unpacklo_pd (a, b);
    __m256d hiab = _mm256_unpackhi_pd (a, b);
    __m256d locd = _mm256_unpacklo_pd (c, d);
    __m256d hicd = _mm256_unpackhi_pd (c, d);

    __m256d r0 = _mm256_permute2f128_pd(loab, locd, 32);
    __m256d r1 = _mm256_permute2f128_pd(hiab, hicd, 32);
    __m256d r2 = _mm256_permute2f128_pd(loab, locd, 49);
    __m256d r3 = _mm256_permute2f128_pd(hiab, hicd, 49);
    return make_tuple(r0, r1, r2, r3);
  }

  // b = Trans(a)
  void TransposeMatrix2 (SliceMatrix<> a, SliceMatrix<> b)
  {
    size_t h = a.Height();
    size_t w = a.Width();
    
    size_t dista = a.Dist();
    size_t distb = b.Dist();

    /*
      __m256i mask_width = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(w&3),
      _mm256_set_epi64x(3, 2, 1, 0));
    */
    __m256i mask_height = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(h&3),
                                               _mm256_set_epi64x(3, 2, 1, 0));
  
    size_t i = 0;
    double * pa1 = &a(0,0);
    double * pb1 = &b(0,0);
    for ( ; i+4 <= w; i+=4)
      {
        size_t j = 0;
        double * pa = pa1;
        double * pb = pb1;
        for ( ; j+4 <= h; j+=4)
          {
            __m256d a1 = _mm256_loadu_pd(pa);
            __m256d a2 = _mm256_loadu_pd(pa+dista);
            __m256d a3 = _mm256_loadu_pd(pa+2*dista);
            __m256d a4 = _mm256_loadu_pd(pa+3*dista);

            auto trans = Transpose(a1, a2, a3, a4);
          
            _mm256_storeu_pd(pb        , get<0> (trans));
            _mm256_storeu_pd(pb+distb  , get<1> (trans));
            _mm256_storeu_pd(pb+2*distb, get<2> (trans));
            _mm256_storeu_pd(pb+3*distb, get<3> (trans));
            pb += 4;
            pa += 4*dista;
          }

        __m256d a1 = _mm256_loadu_pd(pa);
        __m256d a2 = _mm256_loadu_pd(pa+((j+1<h) ? dista : 0));
        __m256d a3 = _mm256_loadu_pd(pa+((j+2<h) ? 2*dista : 0));
        __m256d a4 = _mm256_loadu_pd(pa+((j+3<h) ? 3*dista : 0));
      
        auto trans = Transpose(a1, a2, a3, a4);
      
        _mm256_maskstore_pd(pb        , mask_height, get<0> (trans));
        _mm256_maskstore_pd(pb+distb  , mask_height, get<1> (trans));
        _mm256_maskstore_pd(pb+2*distb, mask_height, get<2> (trans));
        _mm256_maskstore_pd(pb+3*distb, mask_height, get<3> (trans));
      
        pa1 += 4;
        pb1 += 4*distb;
      }

    for ( ; i < w; i++)
      {
        size_t j = 0;
        double * pa = pa1;
        double * pb = pb1;
        for ( ; j+4 <= h; j+=4)
          {
            __m256d valb = _mm256_set_pd(pa[3*dista], pa[2*dista], pa[dista], pa[0]);
            _mm256_storeu_pd(pb, valb);
            pb += 4;
            pa += 4*dista;
          }
        __m256d valb = _mm256_set_pd(pa[(j+3<h) ? 3*dista : 0],
                                     pa[(j+2<h) ? 2*dista : 0],
                                     pa[(j+1<h) ? dista : 0],
                                     pa[0]);
        _mm256_maskstore_pd(pb, mask_height, valb);
        pa1 += 1;
        pb1 += distb;
      }
  }



  void TransposeMatrix (SliceMatrix<> a, SliceMatrix<> b)
  {
    size_t h = a.Height();
    size_t i = 0;
    constexpr size_t bs = 64;
    
    for ( ; i+bs<=h; i+=bs)
      TransposeMatrix2(a.Rows(i,i+bs), b.Cols(i,i+bs));
    TransposeMatrix2(a.Rows(i,h), b.Cols(i,h));    
  }
  


  
  // namespace ngbla
  // {


  // n ... number of SIMDs
  /*
    INLINE void MyScal1x4 (int n, 
    __m256d * a1,
    __m256d * b1, __m256d * b2, __m256d * b3, __m256d * b4,
    __m256d & s1)
    {
    __m256d sum11 = _mm256_setzero_pd();
    __m256d sum12 = _mm256_setzero_pd();
    __m256d sum13 = _mm256_setzero_pd();
    __m256d sum14 = _mm256_setzero_pd();

    for (int i = 0; i < n; i++)
    {
    sum11 += a1[i] * b1[i];
    sum12 += a1[i] * b2[i];
    sum13 += a1[i] * b3[i];
    sum14 += a1[i] * b4[i];
    }
  
    s1 = HAdd (sum11, sum12, sum13, sum14);
    }
  */
  INLINE void MyScal1x4 (int n, 
                         SIMD<double> * a1,
                         SIMD<double> * b1, SIMD<double> * b2, SIMD<double> * b3, SIMD<double> * b4,
                         SIMD<double> & s1)
  {
    SIMD<double> sum11(0.0), sum12(0.0), sum13(0.0), sum14(0.0);

    for (int i = 0; i < n; i++)
      {
        sum11 += a1[i] * b1[i];
        sum12 += a1[i] * b2[i];
        sum13 += a1[i] * b3[i];
        sum14 += a1[i] * b4[i];
      }
  
    s1 = HAdd (sum11.Data(), sum12.Data(), sum13.Data(), sum14.Data());
  }

  template <int R = 0>
  INLINE auto MyScal1x4 (size_t n, double * a1, double * b1, double * b2, double * b3, double * b4)
  {
    SIMD<double> sum11(0.0), sum12(0.0), sum13(0.0), sum14(0.0);

    size_t i = 0;
#pragma nounroll      
    for ( ; i < 4*n; i+=4)
      {
        sum11 += SIMD<double> (_mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b1+i));
        sum12 += SIMD<double> (_mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b2+i));
        sum13 += SIMD<double> (_mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b3+i));
        sum14 += SIMD<double> (_mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b4+i));
      }

    if (R > 0)
      {
        __m256i mask;
        switch (R)
          {
          case 1: mask = _mm256_set_epi64x(0,0,0,-1); break;
          case 2: mask = _mm256_set_epi64x(0,0,-1,-1); break;
          case 3: mask = _mm256_set_epi64x(0,-1,-1,-1); break;
          default: ;
          }
        sum11 += SIMD<double> (_mm256_maskload_pd (a1+i, mask) * _mm256_maskload_pd (b1+i, mask));
        sum12 += SIMD<double> (_mm256_maskload_pd (a1+i, mask) * _mm256_maskload_pd (b2+i, mask));      
        sum13 += SIMD<double> (_mm256_maskload_pd (a1+i, mask) * _mm256_maskload_pd (b3+i, mask));
        sum14 += SIMD<double> (_mm256_maskload_pd (a1+i, mask) * _mm256_maskload_pd (b4+i, mask));
      }
    return SIMD<double> (HAdd (sum11.Data(), sum12.Data(), sum13.Data(), sum14.Data()));
  }


  // n ... number of SIMDs
  /*
    INLINE void MyScal2x4 (int n, 
    __m256d * a1, __m256d * a2,
    __m256d * b1, __m256d * b2, __m256d * b3, __m256d * b4,
    __m256d & s1, __m256d & s2)
    {
    __m256d sum11 = _mm256_setzero_pd();
    __m256d sum21 = _mm256_setzero_pd();

    __m256d sum12 = _mm256_setzero_pd();
    __m256d sum22 = _mm256_setzero_pd();

    __m256d sum13 = _mm256_setzero_pd();
    __m256d sum23 = _mm256_setzero_pd();

    __m256d sum14 = _mm256_setzero_pd();
    __m256d sum24 = _mm256_setzero_pd();

    for (int i = 0; i < n; i++)
    {
    sum11 += a1[i] * b1[i];
    sum21 += a2[i] * b1[i];
    sum12 += a1[i] * b2[i];
    sum22 += a2[i] * b2[i];
    sum13 += a1[i] * b3[i];
    sum23 += a2[i] * b3[i];
    sum14 += a1[i] * b4[i];
    sum24 += a2[i] * b4[i];
    }
  
    s1 = HAdd (sum11, sum12, sum13, sum14);
    s2 = HAdd (sum21, sum22, sum23, sum24);
    }
  */
  INLINE void MyScal2x4 (int n, 
                         SIMD<double> * a1, SIMD<double> * a2,
                         SIMD<double> * b1, SIMD<double> * b2, SIMD<double> * b3, SIMD<double> * b4,
                         SIMD<double> & s1, SIMD<double> & s2)
  {
    SIMD<double> sum11(0.0), sum12(0.0), sum13(0.0), sum14(0.0);
    SIMD<double> sum21(0.0), sum22(0.0), sum23(0.0), sum24(0.0);

#pragma nounroll    
    for (int i = 0; i < n; i++)
      {
        sum11 += a1[i] * b1[i];
        sum21 += a2[i] * b1[i];
        sum12 += a1[i] * b2[i];
        sum22 += a2[i] * b2[i];
        sum13 += a1[i] * b3[i];
        sum23 += a2[i] * b3[i];
        sum14 += a1[i] * b4[i];
        sum24 += a2[i] * b4[i];
      }
  
    s1 = HAdd (sum11.Data(), sum12.Data(), sum13.Data(), sum14.Data());
    s2 = HAdd (sum21.Data(), sum22.Data(), sum23.Data(), sum24.Data());
  }

  /*
    template <int R = 0>
    INLINE auto MyScal2x4 (size_t n, 
    SIMD<double> * a1, SIMD<double> * a2,
    SIMD<double> * b1, SIMD<double> * b2, SIMD<double> * b3, SIMD<double> * b4)
    {
    SIMD<double> sum11(0.0), sum12(0.0), sum13(0.0), sum14(0.0);
    SIMD<double> sum21(0.0), sum22(0.0), sum23(0.0), sum24(0.0);

    size_t i = 0;
    for ( ; i < n; i++)
    {
    sum11 += a1[i] * b1[i];
    sum21 += a2[i] * b1[i];
    sum12 += a1[i] * b2[i];
    sum22 += a2[i] * b2[i];
    sum13 += a1[i] * b3[i];
    sum23 += a2[i] * b3[i];
    sum14 += a1[i] * b4[i];
    sum24 += a2[i] * b4[i];
    }

    if (R > 0)
    {
    __m256i mask;
    switch (R)
    {
    case 1: mask = _mm256_set_epi64x(0,0,0,-1); break;
    case 2: mask = _mm256_set_epi64x(0,0,-1,-1); break;
    case 3: mask = _mm256_set_epi64x(0,-1,-1,-1); break;
    default: ;
    }
    sum11 += _mm256_maskload_pd ((double*)&a1[i].Data(), mask) * _mm256_maskload_pd ((double*)&b1[i].Data(), mask);
    sum12 += _mm256_maskload_pd ((double*)&a1[i].Data(), mask) * _mm256_maskload_pd ((double*)&b2[i].Data(), mask);      
    sum13 += _mm256_maskload_pd ((double*)&a1[i].Data(), mask) * _mm256_maskload_pd ((double*)&b3[i].Data(), mask);
    sum14 += _mm256_maskload_pd ((double*)&a1[i].Data(), mask) * _mm256_maskload_pd ((double*)&b4[i].Data(), mask);
    sum21 += _mm256_maskload_pd ((double*)&a2[i].Data(), mask) * _mm256_maskload_pd ((double*)&b1[i].Data(), mask);
    sum22 += _mm256_maskload_pd ((double*)&a2[i].Data(), mask) * _mm256_maskload_pd ((double*)&b2[i].Data(), mask);      
    sum23 += _mm256_maskload_pd ((double*)&a2[i].Data(), mask) * _mm256_maskload_pd ((double*)&b3[i].Data(), mask);
    sum24 += _mm256_maskload_pd ((double*)&a2[i].Data(), mask) * _mm256_maskload_pd ((double*)&b4[i].Data(), mask);
    }
    
    SIMD<double> s1 = HAdd (sum11.Data(), sum12.Data(), sum13.Data(), sum14.Data());
    SIMD<double> s2 = HAdd (sum21.Data(), sum22.Data(), sum23.Data(), sum24.Data());
    return make_tuple(s1,s2);
    }
  */



  // the best so far ...
  template <int R = 0>
  INLINE auto MyScal2x4 (size_t n, double * a1, double * a2,
                         double * b1, double * b2, double * b3, double * b4)
  {
    SIMD<double> sum11(0.0), sum12(0.0), sum13(0.0), sum14(0.0);
    SIMD<double> sum21(0.0), sum22(0.0), sum23(0.0), sum24(0.0);
    
    size_t n4 = 4*n;
#pragma nounroll    
    for (size_t i = 0; i < n4; i+=4)
      {
        sum11 += _mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b1+i);
        sum21 += _mm256_loadu_pd(a2+i) * _mm256_loadu_pd(b1+i);
        sum12 += _mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b2+i);
        sum22 += _mm256_loadu_pd(a2+i) * _mm256_loadu_pd(b2+i);
        sum13 += _mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b3+i);
        sum23 += _mm256_loadu_pd(a2+i) * _mm256_loadu_pd(b3+i);
        sum14 += _mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b4+i);
        sum24 += _mm256_loadu_pd(a2+i) * _mm256_loadu_pd(b4+i);
      }

    if (R > 0)
      {
        __m256i mask;
        switch (R)
          {
          case 1: mask = _mm256_set_epi64x(0,0,0,-1); break;
          case 2: mask = _mm256_set_epi64x(0,0,-1,-1); break;
          case 3: mask = _mm256_set_epi64x(0,-1,-1,-1); break;
          default: ;
          }
        sum11 += _mm256_maskload_pd (a1+n4, mask) * _mm256_maskload_pd (b1+n4, mask);
        sum12 += _mm256_maskload_pd (a1+n4, mask) * _mm256_maskload_pd (b2+n4, mask);      
        sum13 += _mm256_maskload_pd (a1+n4, mask) * _mm256_maskload_pd (b3+n4, mask);
        sum14 += _mm256_maskload_pd (a1+n4, mask) * _mm256_maskload_pd (b4+n4, mask);
        sum21 += _mm256_maskload_pd (a2+n4, mask) * _mm256_maskload_pd (b1+n4, mask);
        sum22 += _mm256_maskload_pd (a2+n4, mask) * _mm256_maskload_pd (b2+n4, mask);      
        sum23 += _mm256_maskload_pd (a2+n4, mask) * _mm256_maskload_pd (b3+n4, mask);
        sum24 += _mm256_maskload_pd (a2+n4, mask) * _mm256_maskload_pd (b4+n4, mask);
      }
    
    SIMD<double> s1 = HAdd (sum11.Data(), sum12.Data(), sum13.Data(), sum14.Data());
    SIMD<double> s2 = HAdd (sum21.Data(), sum22.Data(), sum23.Data(), sum24.Data());
    return make_tuple(s1,s2);
  }


  INLINE __m256i RestMask (int R)
  {
    switch (R)
      {
      case 0: return _mm256_set_epi64x(0,0,0,0); 
      case 1: return _mm256_set_epi64x(0,0,0,-1); 
      case 2: return _mm256_set_epi64x(0,0,-1,-1); 
      case 3: return _mm256_set_epi64x(0,-1,-1,-1); 
      default: ;
      }
    __assume(false);
    return _mm256_set_epi64x(0,0,0,0); 
  }


  template <int R = 0>
  INLINE auto MyScal3x4 (size_t n, double * pa1, double * pa2, double * pa3,
                         double * pb1, double * pb2, double * pb3, double * pb4)
  {
    SIMD<double> sum11(0.0), sum12(0.0), sum13(0.0), sum14(0.0);
    SIMD<double> sum21(0.0), sum22(0.0), sum23(0.0), sum24(0.0);
    SIMD<double> sum31(0.0), sum32(0.0), sum33(0.0), sum34(0.0);
    
    size_t n4 = 4*n;
#pragma nounroll    
    for (size_t i = 0; i < n4; i+=4)
      {
        __m256d a1 = _mm256_loadu_pd(pa1+i);
        __m256d a2 = _mm256_loadu_pd(pa2+i);
        __m256d a3 = _mm256_loadu_pd(pa3+i);
        
        __m256d b1 = _mm256_loadu_pd(pb1+i);
        sum11 += a1 * b1;
        sum21 += a2 * b1;
        sum31 += a3 * b1;
        __m256d b2 = _mm256_loadu_pd(pb2+i);
        sum12 += a1 * b2;
        sum22 += a2 * b2;
        sum32 += a3 * b2;
        __m256d b3 = _mm256_loadu_pd(pb3+i);
        sum13 += a1 * b3;
        sum23 += a2 * b3;
        sum33 += a3 * b3;
        __m256d b4 = _mm256_loadu_pd(pb4+i);
        sum14 += a1 * b4;
        sum24 += a2 * b4;
        sum34 += a3 * b4;
      }

    if (R > 0)
      {
        __m256d a1 = _mm256_maskload_pd (pa1+n4, RestMask(R));
        __m256d a2 = _mm256_maskload_pd (pa2+n4, RestMask(R));
        __m256d a3 = _mm256_maskload_pd (pa3+n4, RestMask(R));

        __m256d b1 = _mm256_maskload_pd (pb1+n4, RestMask(R));
        sum11 += a1 * b1;
        sum21 += a2 * b1;
        sum31 += a3 * b1;

        __m256d b2 = _mm256_maskload_pd (pb2+n4, RestMask(R));
        sum12 += a1 * b2;
        sum22 += a2 * b2;
        sum32 += a3 * b2;

        __m256d b3 = _mm256_maskload_pd (pb3+n4, RestMask(R));
        sum13 += a1 * b3;
        sum23 += a2 * b3;
        sum33 += a3 * b3;
      
        __m256d b4 = _mm256_maskload_pd (pb4+n4, RestMask(R));
        sum14 += a1 * b4;
        sum24 += a2 * b4;
        sum34 += a3 * b4;
      }
    
    SIMD<double> s1 = HAdd (sum11.Data(), sum12.Data(), sum13.Data(), sum14.Data());
    SIMD<double> s2 = HAdd (sum21.Data(), sum22.Data(), sum23.Data(), sum24.Data());
    SIMD<double> s3 = HAdd (sum31.Data(), sum32.Data(), sum33.Data(), sum34.Data());
    return make_tuple(s1,s2,s3);
  }

#ifndef __GNUC__
  template <int R = 0>
  INLINE auto MyScal3x4 (size_t n, double * pa, size_t da, double * pb, size_t db)
  {
    return MatKernelScalAB<3,4> (SIMD<double>::Size()*n+R, pa, da, pb, db);    
    return MyScal3x4<R> (n, pa, pa+da, pa+2*da, pb, pb+db, pb+2*db, pb+3*db);
  }
#else
  template <int R = 0>
  INLINE auto MyScal3x4 (size_t n, double * pa, size_t da, double * pb, size_t db)
  {
    return MatKernelScalAB<3,4> (SIMD<double>::Size()*n+R, pa, da, pb, db);

    
    SIMD<double> sum11(0.0), sum12(0.0), sum13(0.0), sum14(0.0);
    SIMD<double> sum21(0.0), sum22(0.0), sum23(0.0), sum24(0.0);
    SIMD<double> sum31(0.0), sum32(0.0), sum33(0.0), sum34(0.0);
    
    size_t n4 = 4*n;
#pragma nounroll    
    for (size_t i = 0; i < n4; i+=4)
      {
        /*
          __m256d a1 = _mm256_loadu_pd(pa);
          __m256d a2 = _mm256_loadu_pd(pa+da);
          __m256d a3 = _mm256_loadu_pd(pa+2*da);
        */
        __m256d a1, a2, a3;
        asm ("vmovupd (%[pa]), %[a1]\n\t"
             "vmovupd (%[pa],%[da]), %[a2]\n\t"
             "vmovupd (%[pa],%[da],2), %[a3]\n\t"
             "add $32,%[pa]"
             : [a1] "=&x" (a1), [a2] "=&x" (a2), [a3] "=&x" (a3),
               [pa] "+r" (pa)
             : [da] "r" (8*da) 
             );
        
        //__m256d b1 = _mm256_loadu_pd(pb);
        __m256d b1, b2, b3, b4;
        asm ("vmovupd (%[pb]), %[b1]"
             : [b1] "=&x" (b1) : [pb] "r" (pb), [db] "r" (8*db) );
        
        sum11 += a1 * b1;
        sum21 += a2 * b1;
        sum31 += a3 * b1;
        
        // __m256d b2 = _mm256_loadu_pd(pb+db);
        asm ("vmovupd (%[pb],%[db]), %[b2]"
             : [b2] "=&x" (b2) : [pb] "r" (pb), [db] "r" (8*db) );      
        
        sum12 += a1 * b2;
        sum22 += a2 * b2;
        sum32 += a3 * b2;
        // __m256d b3 = _mm256_loadu_pd(pb+2*db);
        asm ("vmovupd (%[pb],%[db],2), %[b3]"
             : [b3] "=&x" (b3) : [pb] "r" (pb), [db] "r" (8*db) );      
        
        sum13 += a1 * b3;
        sum23 += a2 * b3;
        sum33 += a3 * b3;
        // __m256d b4 = _mm256_loadu_pd(pb+3*db);
        asm ("vmovupd (%[pb],%[db3],8), %[b4]\n\t"
             "add $32,%[pb]"
             : [b4] "=&x" (b4), [pb] "+r" (pb) : [db3] "r" (3*db) );      
        
        sum14 += a1 * b4;
        sum24 += a2 * b4;
        sum34 += a3 * b4;
        // pa += 4;
        // pb += 4;
      }
    
    if (R > 0)
      {
        __m256d a1 = _mm256_maskload_pd (pa, RestMask(R));
        __m256d a2 = _mm256_maskload_pd (pa+da, RestMask(R));
        __m256d a3 = _mm256_maskload_pd (pa+2*da, RestMask(R));

        __m256d b1 = _mm256_maskload_pd (pb, RestMask(R));
        sum11 += a1 * b1;
        sum21 += a2 * b1;
        sum31 += a3 * b1;

        __m256d b2 = _mm256_maskload_pd (pb+db, RestMask(R));
        sum12 += a1 * b2;
        sum22 += a2 * b2;
        sum32 += a3 * b2;

        __m256d b3 = _mm256_maskload_pd (pb+2*db, RestMask(R));
        sum13 += a1 * b3;
        sum23 += a2 * b3;
        sum33 += a3 * b3;
      
        __m256d b4 = _mm256_maskload_pd (pb+3*db, RestMask(R));
        sum14 += a1 * b4;
        sum24 += a2 * b4;
        sum34 += a3 * b4;
      }
    
    SIMD<double> s1 = HAdd (sum11.Data(), sum12.Data(), sum13.Data(), sum14.Data());
    SIMD<double> s2 = HAdd (sum21.Data(), sum22.Data(), sum23.Data(), sum24.Data());
    SIMD<double> s3 = HAdd (sum31.Data(), sum32.Data(), sum33.Data(), sum34.Data());
    return make_tuple(s1,s2,s3);
  }


#endif
  

  /*
    template <int R = 0>
    INLINE auto MyScal2x4 (size_t n, double * a1, double * a2,
    double * b1, double * b2, double * b3, double * b4)
    {
    SIMD<double> sum11(0.0), sum12(0.0), sum13(0.0), sum14(0.0);
    SIMD<double> sum21(0.0), sum22(0.0), sum23(0.0), sum24(0.0);
    
    size_t i = 0;
    #pragma nounroll    
    for ( ; i < 4*n; i+=4)
    {
    sum11 += _mm256_loadu_pd(a1) * _mm256_loadu_pd(b1);
    sum21 += _mm256_loadu_pd(a2) * _mm256_loadu_pd(b1);
    sum12 += _mm256_loadu_pd(a1) * _mm256_loadu_pd(b2);
    sum22 += _mm256_loadu_pd(a2) * _mm256_loadu_pd(b2);
    sum13 += _mm256_loadu_pd(a1) * _mm256_loadu_pd(b3);
    sum23 += _mm256_loadu_pd(a2) * _mm256_loadu_pd(b3);
    sum14 += _mm256_loadu_pd(a1) * _mm256_loadu_pd(b4);
    sum24 += _mm256_loadu_pd(a2) * _mm256_loadu_pd(b4);
    a1 += 4; a2 += 4;
    b1 += 4; b2 += 4; b3 += 4; b4 += 4;
    }

    if (R > 0)
    {
    __m256i mask;
    switch (R)
    {
    case 1: mask = _mm256_set_epi64x(0,0,0,-1); break;
    case 2: mask = _mm256_set_epi64x(0,0,-1,-1); break;
    case 3: mask = _mm256_set_epi64x(0,-1,-1,-1); break;
    default: ;
    }
    sum11 += _mm256_maskload_pd (a1, mask) * _mm256_maskload_pd (b1, mask);
    sum12 += _mm256_maskload_pd (a1, mask) * _mm256_maskload_pd (b2, mask);      
    sum13 += _mm256_maskload_pd (a1, mask) * _mm256_maskload_pd (b3, mask);
    sum14 += _mm256_maskload_pd (a1, mask) * _mm256_maskload_pd (b4, mask);
    sum21 += _mm256_maskload_pd (a2, mask) * _mm256_maskload_pd (b1, mask);
    sum22 += _mm256_maskload_pd (a2, mask) * _mm256_maskload_pd (b2, mask);      
    sum23 += _mm256_maskload_pd (a2, mask) * _mm256_maskload_pd (b3, mask);
    sum24 += _mm256_maskload_pd (a2, mask) * _mm256_maskload_pd (b4, mask);
    }
    
    SIMD<double> s1 = HAdd (sum11.Data(), sum12.Data(), sum13.Data(), sum14.Data());
    SIMD<double> s2 = HAdd (sum21.Data(), sum22.Data(), sum23.Data(), sum24.Data());
    return make_tuple(s1,s2);
    }
  */


  /*
    template <int R = 0>
    INLINE auto MyScal2x4 (size_t n, 
    SIMD<double> * a1, SIMD<double> * a2,
    SIMD<double> * b1, SIMD<double> * b2, SIMD<double> * b3, SIMD<double> * b4)
    {
    SIMD<double> sum11(0.0), sum12(0.0), sum13(0.0), sum14(0.0);
    SIMD<double> sum21(0.0), sum22(0.0), sum23(0.0), sum24(0.0);

    for (size_t i = 0; i < n; i++)
    {
    sum11 += *a1 * *b1;
    sum21 += *a2 * *b1;
    sum12 += *a1 * *b2;
    sum22 += *a2 * *b2;
    sum13 += *a1 * *b3;
    sum23 += *a2 * *b3;
    sum14 += *a1 * *b4;
    sum24 += *a2 * *b4;
    a1++; a2++;
    b1++; b2++; b3++; b4++;
    }

    if (R > 0)
    {
    __m256i mask;
    switch (R)
    {
    case 1: mask = _mm256_set_epi64x(0,0,0,-1); break;
    case 2: mask = _mm256_set_epi64x(0,0,-1,-1); break;
    case 3: mask = _mm256_set_epi64x(0,-1,-1,-1); break;
    default: ;
    }
    sum11 += _mm256_maskload_pd ((double*)&a1->Data(), mask) * _mm256_maskload_pd ((double*)&b1->Data(), mask);
    sum12 += _mm256_maskload_pd ((double*)&a1->Data(), mask) * _mm256_maskload_pd ((double*)&b2->Data(), mask);      
    sum13 += _mm256_maskload_pd ((double*)&a1->Data(), mask) * _mm256_maskload_pd ((double*)&b3->Data(), mask);
    sum14 += _mm256_maskload_pd ((double*)&a1->Data(), mask) * _mm256_maskload_pd ((double*)&b4->Data(), mask);
    sum21 += _mm256_maskload_pd ((double*)&a2->Data(), mask) * _mm256_maskload_pd ((double*)&b1->Data(), mask);
    sum22 += _mm256_maskload_pd ((double*)&a2->Data(), mask) * _mm256_maskload_pd ((double*)&b2->Data(), mask);      
    sum23 += _mm256_maskload_pd ((double*)&a2->Data(), mask) * _mm256_maskload_pd ((double*)&b3->Data(), mask);
    sum24 += _mm256_maskload_pd ((double*)&a2->Data(), mask) * _mm256_maskload_pd ((double*)&b4->Data(), mask);
    }
    
    SIMD<double> s1 = HAdd (sum11.Data(), sum12.Data(), sum13.Data(), sum14.Data());
    SIMD<double> s2 = HAdd (sum21.Data(), sum22.Data(), sum23.Data(), sum24.Data());
    return make_tuple(s1,s2);
    }
  */

  // n ... number of doubles, must be even
  INLINE
  void MyScal4x4 (int n,
                  __m128d * pa1, __m128d * pa2, __m128d * pa3, __m128d * pa4,
                  __m128d * pb1, __m128d * pb2, __m128d * pb3, __m128d * pb4,
                  __m256d & s1, __m256d & s2, __m256d & s3, __m256d & s4)
  {
    __m256d sum11 = _mm256_setzero_pd();
    __m256d sum21 = _mm256_setzero_pd();
    __m256d sum31 = _mm256_setzero_pd();
    __m256d sum41 = _mm256_setzero_pd();
    __m256d sum12 = _mm256_setzero_pd();
    __m256d sum22 = _mm256_setzero_pd();
    __m256d sum32 = _mm256_setzero_pd();
    __m256d sum42 = _mm256_setzero_pd();

    int n2 = n/2;
    for (int i = 0; i < n2; i++)
      {
        __m256d a1 = _mm256_broadcast_pd(pa1+i);
        __m256d a2 = _mm256_broadcast_pd(pa2+i);
        __m256d a3 = _mm256_broadcast_pd(pa3+i);
        __m256d a4 = _mm256_broadcast_pd(pa4+i);

        __m256d b1 = _mm256_broadcast_pd(pb1+i);
        __m256d b2 = _mm256_broadcast_pd(pb2+i);
        __m256d b3 = _mm256_broadcast_pd(pb3+i);
        __m256d b4 = _mm256_broadcast_pd(pb4+i);
        __m256d mb1 = _mm256_blend_pd(b1, b3, 12);
        __m256d mb2 = _mm256_blend_pd(b2, b4, 12);
      
        sum11 += a1 * mb1;
        sum21 += a2 * mb1;
        sum31 += a3 * mb1;
        sum41 += a4 * mb1;
        sum12 += a1 * mb2;
        sum22 += a2 * mb2;
        sum32 += a3 * mb2;
        sum42 += a4 * mb2;
      }
    s1 = _mm256_hadd_pd(sum11, sum12);
    s2 = _mm256_hadd_pd(sum21, sum22);
    s3 = _mm256_hadd_pd(sum31, sum32);
    s4 = _mm256_hadd_pd(sum41, sum42);
  }
        

  // n ... number of full double-pairs  (2*n <= width)
  template <int R>
  INLINE auto MyScal4x4 (size_t n,
                         double * pa1, double * pa2, double * pa3, double * pa4,
                         double * pb1, double * pb2, double * pb3, double * pb4)
  {
    __m256d sum11 = _mm256_setzero_pd();
    __m256d sum21 = _mm256_setzero_pd();
    __m256d sum31 = _mm256_setzero_pd();
    __m256d sum41 = _mm256_setzero_pd();
    __m256d sum12 = _mm256_setzero_pd();
    __m256d sum22 = _mm256_setzero_pd();
    __m256d sum32 = _mm256_setzero_pd();
    __m256d sum42 = _mm256_setzero_pd();

    size_t n2 = 2*n;
#pragma nounroll        
    for (size_t i = 0; i < n2; i+=2)
      {
        __m256d a1 = _mm256_broadcast_pd((__m128d*)(pa1+i));
        __m256d a2 = _mm256_broadcast_pd((__m128d*)(pa2+i));
        __m256d a3 = _mm256_broadcast_pd((__m128d*)(pa3+i));
        __m256d a4 = _mm256_broadcast_pd((__m128d*)(pa4+i));

        __m256d b1 = _mm256_broadcast_pd((__m128d*)(pb1+i));
        __m256d b2 = _mm256_broadcast_pd((__m128d*)(pb2+i));
        __m256d b3 = _mm256_broadcast_pd((__m128d*)(pb3+i));
        __m256d b4 = _mm256_broadcast_pd((__m128d*)(pb4+i));
        __m256d mb1 = _mm256_blend_pd(b1, b3, 12);
        __m256d mb2 = _mm256_blend_pd(b2, b4, 12);
      
        sum11 += a1 * mb1;
        sum21 += a2 * mb1;
        sum31 += a3 * mb1;
        sum41 += a4 * mb1;
        sum12 += a1 * mb2;
        sum22 += a2 * mb2;
        sum32 += a3 * mb2;
        sum42 += a4 * mb2;
      }
    if ( (R==1) || (R==3) )
      {
        __m256d a1 = _mm256_blend_pd(_mm256_broadcast_sd(pa1+n2), _mm256_setzero_pd(), 10);
        __m256d a2 = _mm256_blend_pd(_mm256_broadcast_sd(pa2+n2), _mm256_setzero_pd(), 10);
        __m256d a3 = _mm256_blend_pd(_mm256_broadcast_sd(pa3+n2), _mm256_setzero_pd(), 10);
        __m256d a4 = _mm256_blend_pd(_mm256_broadcast_sd(pa4+n2), _mm256_setzero_pd(), 10);

        __m256d b1 = _mm256_blend_pd(_mm256_broadcast_sd(pb1+n2), _mm256_setzero_pd(), 10);
        __m256d b2 = _mm256_blend_pd(_mm256_broadcast_sd(pb2+n2), _mm256_setzero_pd(), 10);
        __m256d b3 = _mm256_blend_pd(_mm256_broadcast_sd(pb3+n2), _mm256_setzero_pd(), 10);
        __m256d b4 = _mm256_blend_pd(_mm256_broadcast_sd(pb4+n2), _mm256_setzero_pd(), 10);
        __m256d mb1 = _mm256_blend_pd(b1, b3, 12);
        __m256d mb2 = _mm256_blend_pd(b2, b4, 12);
      
        sum11 += a1 * mb1;
        sum21 += a2 * mb1;
        sum31 += a3 * mb1;
        sum41 += a4 * mb1;
        sum12 += a1 * mb2;
        sum22 += a2 * mb2;
        sum32 += a3 * mb2;
        sum42 += a4 * mb2;
      }
    
    __m256d s1 = _mm256_hadd_pd(sum11, sum12);
    __m256d s2 = _mm256_hadd_pd(sum21, sum22);
    __m256d s3 = _mm256_hadd_pd(sum31, sum32);
    __m256d s4 = _mm256_hadd_pd(sum41, sum42);
    return make_tuple (SIMD<double>(s1), SIMD<double>(s2), SIMD<double>(s3), SIMD<double>(s4));
  }
        

  // C += A * Trans(B)


  /*
    template <int R>
    INLINE void AddABt_Rest (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c)
    {
    size_t distc = c.Dist();
    size_t widthc = c.Width();
    size_t heightc = c.Height();

    size_t dista = a.Dist();
    size_t full_vwidtha = size_t(a.Width()) / 4;
    size_t full_vwidtha2 = size_t(a.Width()) / 2;
    size_t distb = b.Dist(); 

    __m256i mask_widthc = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(widthc&3),
    _mm256_set_epi64x(3, 2, 1, 0));
    
    double * pc = &c(0);
    double * pa = &a(0);
    size_t i = 0;

    for ( ; i < heightc-3; i += 4)
    {
    double * pc1 = pc;
    double * pc2 = pc + distc;
    double * pc3 = pc + 2*distc;
    double * pc4 = pc + 3*distc;
    pc = pc + 4*distc;
    double * pa1 = pa;
    double * pa2 = pa+dista;
    double * pa3 = pa+2*dista;
    double * pa4 = pa+3*dista;
    pa = pa + 4*dista;

    double * pb = &b(0);
    size_t j = 0;
    for ( ; j < widthc-3; j += 4)
    {
    double * pb1 = pb;
    double * pb2 = pb+distb;
    double * pb3 = pb+2*distb;
    double * pb4 = pb+3*distb;
    pb = pb + 4*distb;

    auto scal = MyScal4x4<R> (full_vwidtha2, pa1, pa2, pa3, pa4, pb1, pb2, pb3, pb4);
    auto s1 = _mm256_loadu_pd(pc1+j) + get<0>(scal).Data();
    auto s2 = _mm256_loadu_pd(pc2+j) + get<1>(scal).Data();
    auto s3 = _mm256_loadu_pd(pc3+j) + get<2>(scal).Data();
    auto s4 = _mm256_loadu_pd(pc4+j) + get<3>(scal).Data();
    _mm256_storeu_pd(pc1+j, s1);
    _mm256_storeu_pd(pc2+j, s2);
    _mm256_storeu_pd(pc3+j, s3);
    _mm256_storeu_pd(pc4+j, s4);
    }

    if (j < widthc)
    {
    double * pb1 = pb;
    double * pb2 = (j+1 < widthc) ? pb+distb : pb;
    double * pb3 = (j+2 < widthc) ? pb+2*distb : pb;
    double * pb4 = (j+3 < widthc) ? pb+3*distb : pb;
            
    auto scal = MyScal4x4<R> (full_vwidtha2, pa1, pa2, pa3, pa4, pb1, pb2, pb3, pb4);

    SIMD<double> s1 = get<0>(scal);
    SIMD<double> s2 = get<1>(scal);
    SIMD<double> s3 = get<2>(scal);
    SIMD<double> s4 = get<3>(scal);
    s1 += _mm256_maskload_pd (pc1+j, mask_widthc);
    s2 += _mm256_maskload_pd (pc2+j, mask_widthc);
    s3 += _mm256_maskload_pd (pc3+j, mask_widthc);
    s4 += _mm256_maskload_pd (pc4+j, mask_widthc);
    _mm256_maskstore_pd (pc1+j, mask_widthc, s1.Data());
    _mm256_maskstore_pd (pc2+j, mask_widthc, s2.Data());
    _mm256_maskstore_pd (pc3+j, mask_widthc, s3.Data());
    _mm256_maskstore_pd (pc4+j, mask_widthc, s4.Data());
    }
    }

    for ( ; i < heightc-1; i += 2)
    {
    double * pc1 = pc;
    double * pc2 = pc1 + distc;
    pc = pc2 + distc;
    double * pa1 = pa;
    double * pa2 = pa1 + dista;
    pa = pa2 + dista;

    double * pb = &b(0);
    size_t j = 0;
    for ( ; j < widthc-3; j += 4)
    {
    double * pb1 = pb;
    double * pb2 = pb + distb;
    double * pb3 = pb + 2*distb;
    double * pb4 = pb + 3*distb;
    pb = pb + 4*distb;
            
    auto scal = MyScal2x4<R> (full_vwidtha, pa1, pa2, pb1, pb2, pb3, pb4);
    auto s1 = _mm256_loadu_pd(pc1+j) + get<0>(scal).Data();
    auto s2 = _mm256_loadu_pd(pc2+j) + get<1>(scal).Data();
    _mm256_storeu_pd(pc1+j, s1);
    _mm256_storeu_pd(pc2+j, s2);
    }
    if (j < widthc)
    {
    double * pb1 = pb;
    double * pb2 = (j+1 < widthc) ? pb+distb : pb;
    double * pb3 = (j+2 < widthc) ? pb+2*distb : pb;
    double * pb4 = (j+3 < widthc) ? pb+3*distb : pb;
            
    auto scal = MyScal2x4<R> (full_vwidtha, pa1, pa2, pb1, pb2, pb3, pb4);

    SIMD<double> s1 = get<0>(scal);
    SIMD<double> s2 = get<1>(scal);
    s1 += _mm256_maskload_pd (pc1+j, mask_widthc);
    s2 += _mm256_maskload_pd (pc2+j, mask_widthc);
    _mm256_maskstore_pd (pc1+j, mask_widthc, s1.Data());
    _mm256_maskstore_pd (pc2+j, mask_widthc, s2.Data());
    }
    }

    if (i < heightc)
    {
    size_t j = 0;
    double * pb = &b(0);

    for ( ; j < widthc-3; j += 4)
    {
    double * pb1 = pb;
    double * pb2 = pb+distb;
    double * pb3 = pb+2*distb;
    double * pb4 = pb+3*distb;
    pb = pb + 4*distb;
            
    SIMD<double> s1 = MyScal1x4<R> (full_vwidtha, pa, pb1, pb2, pb3, pb4);
    s1 += _mm256_loadu_pd(pc+j);
    _mm256_storeu_pd(pc+j, s1.Data());
    }
    if (j < widthc)
    {
    double * pb1 = pb;
    double * pb2 = (j+1 < widthc) ? pb+distb : pb;
    double * pb3 = (j+2 < widthc) ? pb+2*distb : pb;
    double * pb4 = (j+3 < widthc) ? pb+3*distb : pb;

    SIMD<double> s1 = MyScal1x4<R> (full_vwidtha, pa, pb1, pb2, pb3, pb4);
    s1 += _mm256_maskload_pd (pc+j, mask_widthc);
    _mm256_maskstore_pd (pc+j, mask_widthc, s1.Data());
    }
    }
    }
  */


  
  template <int R>
  INLINE auto MyScal6x4 (size_t n,
                         double * pa1, double * pa2, double * pa3, double * pa4, double * pa5, double * pa6,
                         double * pb1, double * pb2, double * pb3, double * pb4)
  {
    __m256d sum11 = _mm256_setzero_pd();
    __m256d sum21 = _mm256_setzero_pd();
    __m256d sum31 = _mm256_setzero_pd();
    __m256d sum41 = _mm256_setzero_pd();
    __m256d sum51 = _mm256_setzero_pd();
    __m256d sum61 = _mm256_setzero_pd();
    __m256d sum12 = _mm256_setzero_pd();
    __m256d sum22 = _mm256_setzero_pd();
    __m256d sum32 = _mm256_setzero_pd();
    __m256d sum42 = _mm256_setzero_pd();
    __m256d sum52 = _mm256_setzero_pd();
    __m256d sum62 = _mm256_setzero_pd();

    // size_t n2 = 2*n;
#pragma nounroll        
    for (size_t i = 0; i < 2*n; i+=2)
      {
        __m256d b1 = _mm256_broadcast_pd((__m128d*)(pb1));
        __m256d b2 = _mm256_broadcast_pd((__m128d*)(pb2));
        __m256d b3 = _mm256_broadcast_pd((__m128d*)(pb3));
        __m256d b4 = _mm256_broadcast_pd((__m128d*)(pb4));
        __m256d mb1 = _mm256_blend_pd(b1, b3, 12);
        __m256d mb2 = _mm256_blend_pd(b2, b4, 12);

        __m256d a1 = _mm256_broadcast_pd((__m128d*)(pa1));        
        sum11 += a1 * mb1;
        sum12 += a1 * mb2;

        __m256d a2 = _mm256_broadcast_pd((__m128d*)(pa2));
        sum21 += a2 * mb1;
        sum22 += a2 * mb2;

        __m256d a3 = _mm256_broadcast_pd((__m128d*)(pa3));
        sum31 += a3 * mb1;
        sum32 += a3 * mb2;
        
        __m256d a4 = _mm256_broadcast_pd((__m128d*)(pa4));
        sum41 += a4 * mb1;
        sum42 += a4 * mb2;
        
        __m256d a5 = _mm256_broadcast_pd((__m128d*)(pa5));
        sum51 += a5 * mb1;
        sum52 += a5 * mb2;

        __m256d a6 = _mm256_broadcast_pd((__m128d*)(pa6));
        sum61 += a6 * mb1;
        sum62 += a6 * mb2;

        pa1 += 2;
        pa2 += 2;
        pa3 += 2;
        pa4 += 2;
        pa5 += 2;
        pa6 += 2;
        
        pb1 += 2;
        pb2 += 2;
        pb3 += 2;
        pb4 += 2;
      }
    if ( (R==1) || (R==3) )
      {
        __m256d b1 = _mm256_blend_pd(_mm256_broadcast_sd(pb1), _mm256_setzero_pd(), 10);
        __m256d b2 = _mm256_blend_pd(_mm256_broadcast_sd(pb2), _mm256_setzero_pd(), 10);
        __m256d b3 = _mm256_blend_pd(_mm256_broadcast_sd(pb3), _mm256_setzero_pd(), 10);
        __m256d b4 = _mm256_blend_pd(_mm256_broadcast_sd(pb4), _mm256_setzero_pd(), 10);
        __m256d mb1 = _mm256_blend_pd(b1, b3, 12);
        __m256d mb2 = _mm256_blend_pd(b2, b4, 12);
      
        __m256d a1 = _mm256_blend_pd(_mm256_broadcast_sd(pa1), _mm256_setzero_pd(), 10);
        sum11 += a1 * mb1;
        sum12 += a1 * mb2;
        
        __m256d a2 = _mm256_blend_pd(_mm256_broadcast_sd(pa2), _mm256_setzero_pd(), 10);
        sum21 += a2 * mb1;
        sum22 += a2 * mb2;

        __m256d a3 = _mm256_blend_pd(_mm256_broadcast_sd(pa3), _mm256_setzero_pd(), 10);
        sum31 += a3 * mb1;
        sum32 += a3 * mb2;

        __m256d a4 = _mm256_blend_pd(_mm256_broadcast_sd(pa4), _mm256_setzero_pd(), 10);
        sum41 += a4 * mb1;
        sum42 += a4 * mb2;

        __m256d a5 = _mm256_blend_pd(_mm256_broadcast_sd(pa5), _mm256_setzero_pd(), 10);
        sum51 += a5 * mb1;
        sum52 += a5 * mb2;

        __m256d a6 = _mm256_blend_pd(_mm256_broadcast_sd(pa6), _mm256_setzero_pd(), 10);
        sum61 += a6 * mb1;
        sum62 += a6 * mb2;
      }
    
    __m256d s1 = _mm256_hadd_pd(sum11, sum12);
    __m256d s2 = _mm256_hadd_pd(sum21, sum22);
    __m256d s3 = _mm256_hadd_pd(sum31, sum32);
    __m256d s4 = _mm256_hadd_pd(sum41, sum42);
    __m256d s5 = _mm256_hadd_pd(sum51, sum52);
    __m256d s6 = _mm256_hadd_pd(sum61, sum62);
    return make_tuple (SIMD<double>(s1), SIMD<double>(s2), SIMD<double>(s3),
                       SIMD<double>(s4), SIMD<double>(s5), SIMD<double>(s6));
  }

  template <int R>
  INLINE auto MyScal6x4 (size_t n,
                         double * pa, size_t da,
                         double * pb, size_t db)
  {
    return MyScal6x4<R> (n,
                         pa, pa+da, pa+2*da, pa+3*da, pa+4*da, pa+5*da,
                         pb, pb+db, pb+2*db, pb+3*db);
  }  

  template <int R, typename FUNC>
  INLINE void AddABt_Rest (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c,
                           FUNC func)
  {
    size_t distc = c.Dist();
    size_t widthc = c.Width();
    size_t heightc = c.Height();

    size_t dista = a.Dist();
    size_t full_vwidtha = size_t(a.Width()) / 4;
    // size_t full_vwidtha2 = size_t(a.Width()) / 2;
    size_t distb = b.Dist(); 

    __m256i mask_widthc = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(widthc&3),
                                               _mm256_set_epi64x(3, 2, 1, 0));
    
    double * pc = &c(0);
    double * pa = &a(0);
    size_t i = 0;

    /*
      for ( ; i+6 <= heightc; i += 6)
      {
      double * pc1 = pc;
      double * pc2 = pc + distc;
      double * pc3 = pc + 2*distc;
      double * pc4 = pc + 3*distc;
      double * pc5 = pc + 4*distc;
      double * pc6 = pc + 5*distc;
      pc = pc + 6*distc;
      double * pa1 = pa;
      double * pa2 = pa+dista;
      double * pa3 = pa+2*dista;
      double * pa4 = pa+3*dista;
      double * pa5 = pa+4*dista;
      double * pa6 = pa+5*dista;
      pa = pa + 6*dista;

      double * pb = &b(0);
      size_t j = 0;
      for ( ; j+4 <= widthc; j += 4, pb += 4*distb)
      {
      auto scal = MyScal6x4<R> (full_vwidtha2, pa1, dista, pb, distb);
      auto s1 = func (_mm256_loadu_pd(pc1+j), get<0>(scal).Data());
      auto s2 = func (_mm256_loadu_pd(pc2+j), get<1>(scal).Data());
      auto s3 = func (_mm256_loadu_pd(pc3+j), get<2>(scal).Data());
      auto s4 = func (_mm256_loadu_pd(pc4+j), get<3>(scal).Data());
      auto s5 = func (_mm256_loadu_pd(pc5+j), get<4>(scal).Data());
      auto s6 = func (_mm256_loadu_pd(pc6+j), get<5>(scal).Data());
      _mm256_storeu_pd(pc1+j, s1);
      _mm256_storeu_pd(pc2+j, s2);
      _mm256_storeu_pd(pc3+j, s3);
      _mm256_storeu_pd(pc4+j, s4);
      _mm256_storeu_pd(pc5+j, s5);
      _mm256_storeu_pd(pc6+j, s6);
      }

      if (j < widthc)
      {
      double * pb1 = pb;
      double * pb2 = (j+1 < widthc) ? pb+distb : pb;
      double * pb3 = (j+2 < widthc) ? pb+2*distb : pb;
      double * pb4 = (j+3 < widthc) ? pb+3*distb : pb;
            
      auto scal = MyScal6x4<R> (full_vwidtha2, pa1, pa2, pa3, pa4, pa5, pa6, pb1, pb2, pb3, pb4);

      auto s1 = func (_mm256_maskload_pd(pc1+j, mask_widthc), get<0>(scal).Data());
      auto s2 = func (_mm256_maskload_pd(pc2+j, mask_widthc), get<1>(scal).Data());
      auto s3 = func (_mm256_maskload_pd(pc3+j, mask_widthc), get<2>(scal).Data());
      auto s4 = func (_mm256_maskload_pd(pc4+j, mask_widthc), get<3>(scal).Data());
      auto s5 = func (_mm256_maskload_pd(pc5+j, mask_widthc), get<4>(scal).Data());
      auto s6 = func (_mm256_maskload_pd(pc6+j, mask_widthc), get<5>(scal).Data());            
      _mm256_maskstore_pd (pc1+j, mask_widthc, s1);
      _mm256_maskstore_pd (pc2+j, mask_widthc, s2);
      _mm256_maskstore_pd (pc3+j, mask_widthc, s3);
      _mm256_maskstore_pd (pc4+j, mask_widthc, s4);
      _mm256_maskstore_pd (pc5+j, mask_widthc, s5);
      _mm256_maskstore_pd (pc6+j, mask_widthc, s6);            
      }
      }
    */

    for ( ; i+3 <= heightc; i += 3)
      {
        double * pc1 = pc;
        double * pc2 = pc1 + distc;
        double * pc3 = pc2 + distc;
        pc = pc3 + distc;
        double * pa1 = pa;
        double * pa2 = pa1 + dista;
        double * pa3 = pa2 + dista;
        pa = pa3 + dista;

        double * pb = &b(0);
        size_t j = 0;
        for ( ; j+4 <= widthc; j += 4)
          {
            double * pb1 = pb;
            // double * pb2 = pb + distb;
            // double * pb3 = pb + 2*distb;
            // double * pb4 = pb + 3*distb;
            pb = pb + 4*distb;
            
            // auto scal = MyScal3x4<R> (full_vwidtha, pa1, pa2, pa3, pb1, pb2, pb3, pb4);
            auto scal = MyScal3x4<R> (full_vwidtha, pa1, dista, pb1, distb);
            auto s1 = func (_mm256_loadu_pd(pc1+j), get<0>(scal).Data());
            auto s2 = func (_mm256_loadu_pd(pc2+j), get<1>(scal).Data());
            auto s3 = func (_mm256_loadu_pd(pc3+j), get<2>(scal).Data());
            
            _mm256_storeu_pd(pc1+j, s1);
            _mm256_storeu_pd(pc2+j, s2);
            _mm256_storeu_pd(pc3+j, s3);
          }
        if (j < widthc)
          {
            double * pb1 = pb;
            double * pb2 = (j+1 < widthc) ? pb+distb : pb;
            double * pb3 = (j+2 < widthc) ? pb+2*distb : pb;
            double * pb4 = (j+3 < widthc) ? pb+3*distb : pb;
            
            auto scal = MyScal3x4<R> (full_vwidtha, pa1, pa2, pa3, pb1, pb2, pb3, pb4);
            auto s1 = func (_mm256_maskload_pd(pc1+j, mask_widthc), get<0>(scal).Data());
            auto s2 = func (_mm256_maskload_pd(pc2+j, mask_widthc), get<1>(scal).Data());
            auto s3 = func (_mm256_maskload_pd(pc3+j, mask_widthc), get<2>(scal).Data());
            _mm256_maskstore_pd (pc1+j, mask_widthc, s1);
            _mm256_maskstore_pd (pc2+j, mask_widthc, s2);
            _mm256_maskstore_pd (pc3+j, mask_widthc, s3);
          }
      }



    
    for ( ; i+2 <= heightc; i += 2)
      {
        double * pc1 = pc;
        double * pc2 = pc1 + distc;
        pc = pc2 + distc;
        double * pa1 = pa;
        double * pa2 = pa1 + dista;
        pa = pa2 + dista;

        double * pb = &b(0);
        size_t j = 0;
        for ( ; j+4 <= widthc; j += 4)
          {
            double * pb1 = pb;
            double * pb2 = pb + distb;
            double * pb3 = pb + 2*distb;
            double * pb4 = pb + 3*distb;
            pb = pb + 4*distb;
            
            auto scal = MyScal2x4<R> (full_vwidtha, pa1, pa2, pb1, pb2, pb3, pb4);
            auto s1 = func (_mm256_loadu_pd(pc1+j), get<0>(scal).Data());
            auto s2 = func (_mm256_loadu_pd(pc2+j), get<1>(scal).Data());
            
            _mm256_storeu_pd(pc1+j, s1);
            _mm256_storeu_pd(pc2+j, s2);
          }
        if (j < widthc)
          {
            double * pb1 = pb;
            double * pb2 = (j+1 < widthc) ? pb+distb : pb;
            double * pb3 = (j+2 < widthc) ? pb+2*distb : pb;
            double * pb4 = (j+3 < widthc) ? pb+3*distb : pb;
            
            auto scal = MyScal2x4<R> (full_vwidtha, pa1, pa2, pb1, pb2, pb3, pb4);
            auto s1 = func (_mm256_maskload_pd(pc1+j, mask_widthc), get<0>(scal).Data());
            auto s2 = func (_mm256_maskload_pd(pc2+j, mask_widthc), get<1>(scal).Data());

            _mm256_maskstore_pd (pc1+j, mask_widthc, s1);
            _mm256_maskstore_pd (pc2+j, mask_widthc, s2);
          }
      }

    if (i < heightc)
      {
        size_t j = 0;
        double * pb = &b(0);

        for ( ; j+4 <= widthc; j += 4)
          {
            double * pb1 = pb;
            double * pb2 = pb+distb;
            double * pb3 = pb+2*distb;
            double * pb4 = pb+3*distb;
            pb = pb + 4*distb;
            
            SIMD<double> scal = MyScal1x4<R> (full_vwidtha, pa, pb1, pb2, pb3, pb4);
            auto s1 = func (_mm256_loadu_pd(pc+j), scal.Data());
            _mm256_storeu_pd(pc+j, s1);
          }
        if (j < widthc)
          {
            double * pb1 = pb;
            double * pb2 = (j+1 < widthc) ? pb+distb : pb;
            double * pb3 = (j+2 < widthc) ? pb+2*distb : pb;
            double * pb4 = (j+3 < widthc) ? pb+3*distb : pb;

            SIMD<double> scal = MyScal1x4<R> (full_vwidtha, pa, pb1, pb2, pb3, pb4);
            auto s1 = func (_mm256_maskload_pd(pc+j, mask_widthc), scal.Data());
            _mm256_maskstore_pd (pc+j, mask_widthc, s1);
          }
      }
  }
  

  template <int R, typename FUNC>
  void AddABt_Rest1 (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c,
                     FUNC func)
  {
    constexpr size_t bs = 32;  // height b
    size_t hb = b.Height();
    size_t i = 0;
    for ( ; i+bs <= hb; i += bs)
      AddABt_Rest<R> (a, b.Rows(i,i+bs), c.Cols(i,i+bs), func);
    if (i < hb)
      AddABt_Rest<R> (a, b.Rows(i,hb), c.Cols(i,hb), func);      
  }

  template <int R, typename FUNC>
  void AddABt_Rest2 (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c,
                     FUNC func)
  {
    constexpr size_t bs = 96;  // height a
    size_t ha = a.Height();
    size_t i = 0;
    for ( ; i+bs <= ha; i += bs)
      AddABt_Rest1<R> (a.Rows(i,i+bs), b, c.Rows(i,i+bs), func);
    if (i < ha)
      AddABt_Rest1<R> (a.Rows(i,ha), b, c.Rows(i,ha), func);      
  }

  
  template <typename FUNC>
  void AddABt2 (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c,
                FUNC func)
  {
    constexpr size_t bs = 256; // inner-product loop
    size_t wa = a.Width();
    size_t i = 0;
    for ( ; i+bs <= wa; i += bs)
      AddABt_Rest2<0> (a.Cols(i, i+bs), b.Cols(i,i+bs), c, func);
    if (i == a.Width()) return;
    
    switch (wa & 3)
      {
      case 0: AddABt_Rest2<0>(a.Cols(i,wa), b.Cols(i,wa), c, func); break;
      case 1: AddABt_Rest2<1>(a.Cols(i,wa), b.Cols(i,wa), c, func); break;
      case 2: AddABt_Rest2<2>(a.Cols(i,wa), b.Cols(i,wa), c, func); break;
      case 3: AddABt_Rest2<3>(a.Cols(i,wa), b.Cols(i,wa), c, func); break;
      }
  }
  
  void XXAddABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    // c += a * Trans(b);
    // return;
    AddABt2 (a, b, c.AddSize(a.Height(),b.Height()), [] (auto c, auto ab) { return c+ab; });
  }

  void XXSubABt (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c)
  {
    // c -= a * Trans(b);
    // return;
    AddABt2 (a, b, c, [] (auto c, auto ab) { return c-ab; });
  }







#ifdef NONE
  void AddABtSymV1 (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c)
  {
    int i = 0;
    // clear overhead
    if (a.Width() != 4*a.VWidth())
      {
        int r = 4*a.VWidth()-a.Width();
        __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(r),
                                            _mm256_set_epi64x(0,1,2,3));

        __m256d zero = _mm256_setzero_pd();
        for (int i = 0; i < a.Height(); i++)
          _mm256_maskstore_pd((double*)&a.Get(i, a.VWidth()-1), mask, zero);
        for (int i = 0; i < b.Height(); i++)
          _mm256_maskstore_pd((double*)&b.Get(i, b.VWidth()-1), mask, zero);
      }
  
    if (a.VWidth() <= 0) return;
  
    for ( ; i+1 < c.Height(); i += 2)
      {
        int j = 0;
        for ( ; j < i; j += 4)
          {
            SIMD<double> s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            s1 += _mm256_loadu_pd(&c(i,j));
            s2 += _mm256_loadu_pd(&c(i+1,j));
            _mm256_storeu_pd(&c(i,j), s1.Data());
            _mm256_storeu_pd(&c(i+1,j), s2.Data());
          }
        if (j <= i)
          {
            SIMD<double> s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            __m128d s1l = _mm256_extractf128_pd(s1.Data(), 0);
            __m128d s2l = _mm256_extractf128_pd(s2.Data(), 0);
            s1l += _mm_loadu_pd(&c(i,j));
            s2l += _mm_loadu_pd(&c(i+1,j));
            _mm_storeu_pd(&c(i,j), s1l);
            _mm_storeu_pd(&c(i+1,j), s2l);
          }
      }
  
    if (i < c.Height())
      {
        int j = 0;
        for ( ; j+3 < c.Width(); j += 4)
          {
            SIMD<double> s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            s1 += _mm256_loadu_pd(&c(i,j));
            _mm256_storeu_pd(&c(i,j), s1.Data());
          }
        if (j < c.Width())
          {
            SIMD<double> s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
            for (int j2 = 0; j2 < c.Width()-j; j2++)
              c(i,j+j2) += ((double*)(&s1))[j2];
          }
      }
  }
#endif

  /*
    INLINE void AddABtSym (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c)
    {
    // clear overhead
    if (a.Width() != 4*a.VWidth())
    {
    int r = 4*a.VWidth()-a.Width();
    __m256i mask = _mm256_cmpgt_epi64(_mm256_set1_epi64x(r),
    _mm256_set_epi64x(0,1,2,3));

    __m256d zero = _mm256_setzero_pd();
    for (int i = 0; i < a.Height(); i++)
    _mm256_maskstore_pd((double*)&a.Get(i, a.VWidth()-1), mask, zero);
    for (int i = 0; i < b.Height(); i++)
    _mm256_maskstore_pd((double*)&b.Get(i, b.VWidth()-1), mask, zero);
    }
  
    if (a.VWidth() <= 0) return;
  
    int j = 0;
    for ( ; j < c.Width()-3; j += 4)
    {
    int i = j;
    double * pc = &c(i,j);
    for ( ; i < c.Height()-1; i += 2)
    {
    __m256d s1, s2;
    MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
    &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);

    // s1 += _mm256_loadu_pd(&c(i,j));
    // s2 += _mm256_loadu_pd(&c(i+1,j));
    // _mm256_storeu_pd(&c(i,j), s1);
    // _mm256_storeu_pd(&c(i+1,j), s2);


    s1 += _mm256_loadu_pd(pc);
    _mm256_storeu_pd(pc, s1);
    pc += c.Dist();
    s2 += _mm256_loadu_pd(pc);
    _mm256_storeu_pd(pc, s2);
    pc += c.Dist();          
    }
    if (i < c.Height())
    {
    __m256d s1;
    MyScal1x4 (a.VWidth(), &a.Get(i,0),
    &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1);

    s1 += _mm256_loadu_pd(pc);
    _mm256_storeu_pd(pc, s1);
    }
    }

    for ( ; j < c.Width(); j++)
    for (int i = j; i < c.Height(); i++)
    c(i,j) += InnerProduct(a.Row(i), b.Row(j));
    }
  */

  /*
  void AddABtSym (AFlatMatrix<double> a, AFlatMatrix<double> b, BareSliceMatrix<double> bc)
  {
    auto c = bc.AddSize(a.Height(), b.Height());
    // clear overhead
    if (a.Width() != 4*a.VWidth())
      {
        int r = 4*a.VWidth()-a.Width();
        __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(r),
                                            _mm256_set_epi64x(0,1,2,3));

        __m256d zero = _mm256_setzero_pd();
        for (int i = 0; i < a.Height(); i++)
          _mm256_maskstore_pd((double*)&a.Get(i, a.VWidth()-1), mask, zero);
        for (int i = 0; i < b.Height(); i++)
          _mm256_maskstore_pd((double*)&b.Get(i, b.VWidth()-1), mask, zero);
      }
  
    if (a.VWidth() <= 0) return;
  
    int j = 0;
    for ( ; j+3 < c.Width(); j += 4)
      {
        int i = j;
        double * pc = &c(i,j);
        for ( ; i+3 < c.Height(); i += 4)
          {
            __m256d s1, s2, s3, s4;
            MyScal4x4 (4*a.VWidth(),
                       (__m128d*)&a.Get(i,0), (__m128d*)&a.Get(i+1,0), (__m128d*)&a.Get(i+2,0), (__m128d*)&a.Get(i+3,0),
                       (__m128d*)&b.Get(j,0), (__m128d*)&b.Get(j+1,0), (__m128d*)&b.Get(j+2,0), (__m128d*)&b.Get(j+3,0), s1, s2, s3, s4);

            s1 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s1);
            pc += c.Dist();
            s2 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s2);
            pc += c.Dist();          
            s3 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s3);
            pc += c.Dist();
            s4 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s4);
            pc += c.Dist();          
          }

        if (i+1 < c.Height())
          {
            SIMD<double> s1, s2;
            MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);

            s1 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s1.Data());
            pc += c.Dist();
            s2 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s2.Data());
            pc += c.Dist();
            i += 2;
          }

        if (i < c.Height())
          {
            SIMD<double> s1;
            MyScal1x4 (a.VWidth(), &a.Get(i,0),
                       &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1);

            s1 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s1.Data());
          }
      }

    for ( ; j < c.Width(); j++)
      for (int i = j; i < c.Height(); i++)
        c(i,j) += InnerProduct(a.Row(i), b.Row(j));
  }
  */


  /*
    // now in ngblas, new version
  template <int R>
  void AddABtSymR (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> bc)
  {
    auto c = bc.AddSize(a.Height(), b.Height());
    if (a.Width() == 0) return;
  
    size_t j = 0;
    for ( ; j+3 < c.Width(); j += 4)
      {
        size_t i = j;
        double * pc = &c(i,j);
        for ( ; i+3 < c.Height(); i += 4)
          {
            SIMD<double> s1, s2, s3, s4;
            tie (s1,s2,s3,s4) = 
              MyScal4x4<R> (a.Width()/2,
                            &a(i,0), &a(i+1,0), &a(i+2,0), &a(i+3,0),
                            &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0));
            
            s1 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s1.Data());
            pc += c.Dist();
            s2 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s2.Data());
            pc += c.Dist();          
            s3 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s3.Data());
            pc += c.Dist();
            s4 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s4.Data());
            pc += c.Dist();          
          }

        if (i+1 < c.Height())
          {
            SIMD<double> s1, s2;
            tie (s1, s2) = 
              MyScal2x4<R> (a.Width()/4,
                         &a(i,0), &a(i+1,0),
                         &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0));
            s1 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s1.Data());
            pc += c.Dist();
            s2 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s2.Data());
            pc += c.Dist();
            i += 2;
          }

        if (i < c.Height())
          {
            SIMD<double> s1;
            s1 = MyScal1x4<R> (a.Width()/4, &a(i,0),
                               &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0));
            
            s1 += _mm256_loadu_pd(pc);
            _mm256_storeu_pd(pc, s1.Data());
          }
      }

    // tuning for rest 3
    if (j+3 == c.Width())
      {
        double sum00 = c(j  ,j  );
        double sum10 = c(j+1,j  );
        double sum11 = c(j+1,j+1);
        double sum20 = c(j+2,j  );
        double sum21 = c(j+2,j+1);
        double sum22 = c(j+2,j+2);
        for (size_t k = 0; k < a.Width(); k++)
          {
            sum00 += a(j  ,k) * b(j  ,k);
            sum10 += a(j+1,k) * b(j  ,k);
            sum11 += a(j+1,k) * b(j+1,k);
            sum20 += a(j+2,k) * b(j  ,k);
            sum21 += a(j+2,k) * b(j+1,k);
            sum22 += a(j+2,k) * b(j+2,k);
          }
        c(j  ,j  ) = sum00;
        c(j+1,j  ) = sum10;
        c(j+1,j+1) = sum11;
        c(j+2,j  ) = sum20;
        c(j+2,j+1) = sum21;
        c(j+2,j+2) = sum22;
        return;
      }
    for ( ; j < c.Width(); j++)
      for (size_t i = j; i < c.Height(); i++)
        c(i,j) += InnerProduct(a.Row(i), b.Row(j));
  }


  void AddABtSym (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> bc)
  {
    switch (a.Width() % 4)
      {
      case 0: AddABtSymR<0> (a, b, bc); break;
      case 1: AddABtSymR<1> (a, b, bc); break;
      case 2: AddABtSymR<2> (a, b, bc); break;
      case 3: AddABtSymR<3> (a, b, bc); break;
      default: __assume(false);
      }
  }
  */




  // n ... number of doubles
  INLINE
  void MyScal4x4 (int n,
                  double * pa1, double * pa2, double * pa3, double * pa4,
                  Complex * _pb1, Complex * _pb2, Complex * _pb3, Complex * _pb4,
                  __m256d & sum11, __m256d & sum21, __m256d & sum31, __m256d & sum41,
                  __m256d & sum12, __m256d & sum22, __m256d & sum32, __m256d & sum42)
  {
    __m128d * pb1 = reinterpret_cast<__m128d*> (_pb1);
    __m128d * pb2 = reinterpret_cast<__m128d*> (_pb2);
    __m128d * pb3 = reinterpret_cast<__m128d*> (_pb3);
    __m128d * pb4 = reinterpret_cast<__m128d*> (_pb4);


    sum11 = _mm256_setzero_pd();
    sum21 = _mm256_setzero_pd();
    sum31 = _mm256_setzero_pd();
    sum41 = _mm256_setzero_pd();
    sum12 = _mm256_setzero_pd();
    sum22 = _mm256_setzero_pd();
    sum32 = _mm256_setzero_pd();
    sum42 = _mm256_setzero_pd();

    for (int i = 0; i < n; i++)
      {
        __m256d a1 = _mm256_broadcast_sd(pa1+i);
        __m256d a2 = _mm256_broadcast_sd(pa2+i);
        __m256d a3 = _mm256_broadcast_sd(pa3+i);
        __m256d a4 = _mm256_broadcast_sd(pa4+i);

        __m256d b1 = _mm256_broadcast_pd(pb1+i);
        __m256d b2 = _mm256_broadcast_pd(pb2+i);
        __m256d b3 = _mm256_broadcast_pd(pb3+i);
        __m256d b4 = _mm256_broadcast_pd(pb4+i);
        __m256d mb1 = _mm256_blend_pd(b1, b2, 12);
        __m256d mb2 = _mm256_blend_pd(b3, b4, 12);
      
        sum11 += a1 * mb1;
        sum21 += a2 * mb1;
        sum31 += a3 * mb1;
        sum41 += a4 * mb1;
        sum12 += a1 * mb2;
        sum22 += a2 * mb2;
        sum32 += a3 * mb2;
        sum42 += a4 * mb2;
      }
    /*
      s1 = _mm256_hadd_pd(sum11, sum12);
      s2 = _mm256_hadd_pd(sum21, sum22);
      s3 = _mm256_hadd_pd(sum31, sum32);
      s4 = _mm256_hadd_pd(sum41, sum42);
    */
  }

  // n ... number of doubles
  template <typename TI>
  INLINE void MyScal2x4 (TI n,
                         double * pa1, double * pa2,
                         Complex * _pb1, Complex * _pb2, Complex * _pb3, Complex * _pb4,
                         __m256d & sum11, __m256d & sum21,
                         __m256d & sum12, __m256d & sum22)
  {
    __m128d * pb1 = reinterpret_cast<__m128d*> (_pb1);
    __m128d * pb2 = reinterpret_cast<__m128d*> (_pb2);
    __m128d * pb3 = reinterpret_cast<__m128d*> (_pb3);
    __m128d * pb4 = reinterpret_cast<__m128d*> (_pb4);
    
    sum11 = _mm256_setzero_pd();
    sum21 = _mm256_setzero_pd();
    sum12 = _mm256_setzero_pd();
    sum22 = _mm256_setzero_pd();
    
    for (TI i = 0; i < n; i++)
      {
        __m256d a1 = _mm256_broadcast_sd(pa1+i);
        __m256d a2 = _mm256_broadcast_sd(pa2+i);

        __m256d b1 = _mm256_broadcast_pd(pb1+i);
        __m256d b2 = _mm256_broadcast_pd(pb2+i);
        __m256d b3 = _mm256_broadcast_pd(pb3+i);
        __m256d b4 = _mm256_broadcast_pd(pb4+i);
        __m256d mb1 = _mm256_blend_pd(b1, b2, 12);
        __m256d mb2 = _mm256_blend_pd(b3, b4, 12);
      
        sum11 += a1 * mb1;
        sum21 += a2 * mb1;
        sum12 += a1 * mb2;
        sum22 += a2 * mb2;
      }
  }

  // n ... number of doubles
  INLINE
  void MyScal2x2 (int n,
                  double * pa1, double * pa2,
                  Complex * _pb1, Complex * _pb2,
                  __m256d & sum11, __m256d & sum21)
  {
    __m128d * pb1 = reinterpret_cast<__m128d*> (_pb1);
    __m128d * pb2 = reinterpret_cast<__m128d*> (_pb2);
    
    sum11 = _mm256_setzero_pd();
    sum21 = _mm256_setzero_pd();
    
    for (int i = 0; i < n; i++)
      {
        __m256d a1 = _mm256_broadcast_sd(pa1+i);
        __m256d a2 = _mm256_broadcast_sd(pa2+i);

        __m256d b1 = _mm256_broadcast_pd(pb1+i);
        __m256d b2 = _mm256_broadcast_pd(pb2+i);
        __m256d mb1 = _mm256_blend_pd(b1, b2, 12);

        sum11 += a1 * mb1;
        sum21 += a2 * mb1;
      }
  }



  // n ... number of doubles
  INLINE
  void MyScal1x4 (int n,
                  double * pa1,
                  Complex * _pb1, Complex * _pb2, Complex * _pb3, Complex * _pb4,
                  __m256d & sum11,
                  __m256d & sum12)
  {
    __m128d * pb1 = reinterpret_cast<__m128d*> (_pb1);
    __m128d * pb2 = reinterpret_cast<__m128d*> (_pb2);
    __m128d * pb3 = reinterpret_cast<__m128d*> (_pb3);
    __m128d * pb4 = reinterpret_cast<__m128d*> (_pb4);
    
    sum11 = _mm256_setzero_pd();
    sum12 = _mm256_setzero_pd();
    
    for (int i = 0; i < n; i++)
      {
        __m256d a1 = _mm256_broadcast_sd(pa1+i);

        __m256d b1 = _mm256_broadcast_pd(pb1+i);
        __m256d b2 = _mm256_broadcast_pd(pb2+i);
        __m256d b3 = _mm256_broadcast_pd(pb3+i);
        __m256d b4 = _mm256_broadcast_pd(pb4+i);
        __m256d mb1 = _mm256_blend_pd(b1, b2, 12);
        __m256d mb2 = _mm256_blend_pd(b3, b4, 12);
      
        sum11 += a1 * mb1;
        sum12 += a1 * mb2;
      }
  }


  void AddABtSym (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
  {  
    // static Timer t("AddABtSym, double-complex"); RegionTimer reg(t);
    // c += a * Trans(b);
    // return;
    /*
      cout << "a = " << endl << a << endl;
      cout << "b = " << endl << b << endl;
      cout << "c = " << endl << c << endl;
      cout << "c1 = " << endl << c+a*Trans(b) << endl;
    */
  
    if (a.Width() <= 0) return;
  
    int j = 0;
    for ( ; j+3 < c.Width(); j += 4)
      {
        int i = j;
        Complex * pc = &c(i,j);
        for ( ; i+3 < c.Height(); i += 4)
          {
            __m256d s11, s21, s31, s41;
            __m256d s12, s22, s32, s42;
            MyScal4x4 (a.Width(),
                       &a(i,0), &a(i+1,0), &a(i+2,0), &a(i+3,0),
                       &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0),
                       s11, s21, s31, s41,
                       s12, s22, s32, s42);

            double * pcd = reinterpret_cast<double*> (pc);
            s11 += _mm256_loadu_pd(pcd);
            _mm256_storeu_pd(pcd, s11);
            s12 += _mm256_loadu_pd(pcd+4);
            _mm256_storeu_pd(pcd+4, s12);
          
            pc += c.Dist();
            pcd = reinterpret_cast<double*> (pc);          
            s21 += _mm256_loadu_pd(pcd);
            _mm256_storeu_pd(pcd, s21);
            s22 += _mm256_loadu_pd(pcd+4);
            _mm256_storeu_pd(pcd+4, s22);

            pc += c.Dist();
            pcd = reinterpret_cast<double*> (pc);                    
            s31 += _mm256_loadu_pd(pcd);
            _mm256_storeu_pd(pcd, s31);
            s32 += _mm256_loadu_pd(pcd+4);
            _mm256_storeu_pd(pcd+4, s32);

            pc += c.Dist();
            pcd = reinterpret_cast<double*> (pc);                              
            s41 += _mm256_loadu_pd(pcd);
            _mm256_storeu_pd(pcd, s41);
            s42 += _mm256_loadu_pd(pcd+4);
            _mm256_storeu_pd(pcd+4, s42);
            pc += c.Dist();          
          }

        if (i+1 < c.Height())
          {
            __m256d s11, s21;
            __m256d s12, s22;
            MyScal2x4 (a.Width(), &a(i,0), &a(i+1,0),
                       &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0),
                       s11, s21, s12, s22);

            double * pcd = reinterpret_cast<double*> (pc);          
            s11 += _mm256_loadu_pd(pcd);
            _mm256_storeu_pd(pcd, s11);
            s12 += _mm256_loadu_pd(pcd+4);
            _mm256_storeu_pd(pcd+4, s12);
          
            pc += c.Dist();
            pcd = reinterpret_cast<double*> (pc);          
            s21 += _mm256_loadu_pd(pcd);
            _mm256_storeu_pd(pcd, s21);
            s22 += _mm256_loadu_pd(pcd+4);
            _mm256_storeu_pd(pcd+4, s22);
          
            pc += c.Dist();
            i += 2;
          }

        if (i < c.Height())
          {
            __m256d s11, s12;
            MyScal1x4 (a.Width(), &a(i,0),
                       &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0), s11, s12);
          
            double * pcd = reinterpret_cast<double*> (pc);                    
            s11 += _mm256_loadu_pd(pcd);
            _mm256_storeu_pd(pcd, s11);
            s12 += _mm256_loadu_pd(pcd+4);
            _mm256_storeu_pd(pcd+4, s12);
          }
      }

    for ( ; j < c.Width(); j++)
      for (int i = j; i < c.Height(); i++)
        c(i,j) += InnerProduct(a.Row(i), b.Row(j));

    /*
      cout << "c2 = " << endl << c << endl;
      if (c.Height() > 4)
      exit(1);
    */
  }

  void AddABtSym (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
  {
    c += a * Trans(b) | Lapack;
    // cout << "AddABTsym complex - complex not implemented" << endl;
  }




  void AddABt (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
  {
    /*
      int ah = a.Height();
      int bh = b.Height();
      int w = a.Width();

      if (bh > 128 && w > 128)
      {
      if (w > bh)
      {
      int w2 = w/2;
      w2 &= -4;
      AddABt (a.Cols(0, w2), b.Cols(0,w2), c);
      AddABt (a.Cols(w2, w), b.Cols(w2, w), c);
      return;
      }
      else
      {
      int h2 = bh/2;
      h2 &= -4;
      AddABt (a, b.Rows(0,h2), c.Cols(0,h2));
      AddABt (a, b.Rows(h2, bh), c.Cols(h2, bh));
      return;
      }
      }
    */
    /*
      if (a.Height() > 128)
      {
      int h2 = a.Height()/2;
      h2 &= -4;
      AddABt (a.Rows(0,h2), b, c.Rows(0,h2));
      AddABt (a.Rows(h2, a.Height()), b, c.Rows(h2, c.Height()));
      return;
      }
    */
  
    // c += a * Trans(b);
    // return;
  
    int i = 0;

    if (a.Width() <= 0) return;



    for ( ; i+3 < c.Height(); i += 4)
      {
        int j = 0;
        for ( ; j+3 < c.Width(); j += 4)
          {
            __m256d s11, s21, s31, s41, s12, s22, s32, s42;
            MyScal4x4 (a.Width(), &a(i,0), &a(i+1,0), &a(i+2,0), &a(i+3,0),
                       &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0),
                       s11, s21, s31, s41, s12, s22, s32, s42);
            s11 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j)));
            s21 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+1,j)));
            s31 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+2,j)));
            s41 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+3,j)));
            s12 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j+2)));
            s22 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+1,j+2)));
            s32 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+2,j+2)));
            s42 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+3,j+2)));
            
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j)), s11);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+1,j)), s21);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+2,j)), s31);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+3,j)), s41);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j+2)), s12);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+1,j+2)), s22);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+2,j+2)), s32);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+3,j+2)), s42);
          }
        
        for ( ; j+1 < c.Width(); j += 2)
          for (int i2 = i; i2 < i+4; i2 += 2)
            {
              __m256d s11, s21;
              MyScal2x2 (a.Width(), &a(i2,0), &a(i2+1,0),
                         &b(j,0), &b(j+1,0),
                         s11, s21);
              s11 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i2,j)));
              s21 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i2+1,j)));
              _mm256_storeu_pd(reinterpret_cast<double*>(&c(i2,j)), s11);
              _mm256_storeu_pd(reinterpret_cast<double*>(&c(i2+1,j)), s21);
            }

        for (  ; j < c.Width(); j++)
          for (int i2 = i; i2 < i+4; i2 ++)
            c(i2,j) += InnerProduct(a.Row(i2), b.Row(j));
      }


    
    for ( ; i+1 < c.Height(); i += 2)
      {
        int j = 0;
        for ( ; j+3 < c.Width(); j += 4)
          {
            __m256d s11, s21, s12, s22;
            MyScal2x4 (a.Width(), &a(i,0), &a(i+1,0),
                       &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0),
                       s11, s21, s12, s22);
            s11 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j)));
            s21 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+1,j)));
            s12 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j+2)));
            s22 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+1,j+2)));
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j)), s11);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+1,j)), s21);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j+2)), s12);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+1,j+2)), s22);
          }
        for ( ; j+1 < c.Width(); j += 2)
          {
            __m256d s11, s21;
            MyScal2x2 (a.Width(), &a(i,0), &a(i+1,0),
                       &b(j,0), &b(j+1,0),
                       s11, s21);
            s11 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j)));
            s21 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i+1,j)));
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j)), s11);
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i+1,j)), s21);
          }
        for (  ; j < c.Width(); j++)
          {
            c(i,j) += InnerProduct(a.Row(i), b.Row(j));
            c(i+1,j) += InnerProduct(a.Row(i+1), b.Row(j));
          }
        /*
          if (j < c.Width())
          {
          __m256d s1, s2;
          MyScal2x4 (a.VWidth(), &a.Get(i,0), &a.Get(i+1,0),
          &b.Get(j,0), &b.Get(j+1,0), &b.Get(j+2,0), &b.Get(j+3,0), s1, s2);
          for (int j2 = 0; j2 < c.Width()-j; j2++)
          {
          c(i,j+j2) += ((double*)(&s1))[j2];
          c(i+1,j+j2) += ((double*)(&s2))[j2];
          }
          }
        */
      }

  
    if (i < c.Height())
      {
        int j = 0;
        for ( ; j+3 < c.Width(); j += 4)
          {
            __m256d s11, s12;
            MyScal1x4 (a.Width(), &a(i,0),
                       &b(j,0), &b(j+1,0), &b(j+2,0), &b(j+3,0), s11, s12);
            s11 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j)));
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j)), s11);
            s12 += _mm256_loadu_pd(reinterpret_cast<double*>(&c(i,j+2)));
            _mm256_storeu_pd(reinterpret_cast<double*>(&c(i,j+2)), s12);
          }
        for ( ; j < c.Width(); j++)
          {
            c(i,j) += InnerProduct(a.Row(i), b.Row(j));            
          }
      }
  }

  void AddABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c)
  {
    c += a * Trans(b) | Lapack;  
    // cout << "addabt complex-complex not implemented" << endl;
  }






  // micro-kernels for A^t B
  // pa+i * pb+j   ->  C 3x16
  void MyScal3x16Trans (size_t ninner,
                        double * pa, size_t da,
                        double * pb, size_t db,
                        double * pc, size_t dc)
  {
    __assume (ninner > 0);
#ifndef __GNUC__    
    double * hpc = pc;
    __m256d sum11 = _mm256_loadu_pd(pc);
    __m256d sum12 = _mm256_loadu_pd(pc+4);
    __m256d sum13 = _mm256_loadu_pd(pc+8);
    __m256d sum14 = _mm256_loadu_pd(pc+12);
    pc += dc;
    __m256d sum21 = _mm256_loadu_pd(pc);
    __m256d sum22 = _mm256_loadu_pd(pc+4);
    __m256d sum23 = _mm256_loadu_pd(pc+8);
    __m256d sum24 = _mm256_loadu_pd(pc+12);
    pc += dc;
    __m256d sum31 = _mm256_loadu_pd(pc);
    __m256d sum32 = _mm256_loadu_pd(pc+4);
    __m256d sum33 = _mm256_loadu_pd(pc+8);
    __m256d sum34 = _mm256_loadu_pd(pc+12);
    pc += dc;
    pc = hpc;
#else
    __m256d sum11, sum12, sum13, sum14;
    __m256d sum21, sum22, sum23, sum24;
    __m256d sum31, sum32, sum33, sum34;

    asm (
         "vmovupd (%[pc]), %[sum11]\n\t"
         "vmovupd 32(%[pc]), %[sum12]\n\t"
         "vmovupd 64(%[pc]), %[sum13]\n\t"
         "vmovupd 96(%[pc]), %[sum14]\n\t"

         "vmovupd (%[pc],%[dc8]), %[sum21]\n\t"
         "vmovupd 32(%[pc],%[dc8]), %[sum22]\n\t"
         "vmovupd 64(%[pc],%[dc8]), %[sum23]\n\t"
         "vmovupd 96(%[pc],%[dc8]), %[sum24]\n\t"

         "vmovupd (%[pc],%[dc8],2), %[sum31]\n\t"
         "vmovupd 32(%[pc],%[dc8],2), %[sum32]\n\t"
         "vmovupd 64(%[pc],%[dc8],2), %[sum33]\n\t"
         "vmovupd 96(%[pc],%[dc8],2), %[sum34]\n\t"

         :
         [sum11] "+x" (sum11), [sum21] "+x" (sum21), [sum31] "+x" (sum31),
         [sum12] "+x" (sum12), [sum22] "+x" (sum22), [sum32] "+x" (sum32),
         [sum13] "+x" (sum13), [sum23] "+x" (sum23), [sum33] "+x" (sum33),
         [sum14] "+x" (sum14), [sum24] "+x" (sum24), [sum34] "+x" (sum34)
         : [pc] "r" (pc), [dc8]  "r" (8*dc)
         );
#endif
  
#pragma nounroll    
    for (size_t i = 0; i < ninner; i++, pa += da, pb += db)
      {
        __m256d a1 = _mm256_set1_pd(*pa);
        __m256d a2 = _mm256_set1_pd(pa[1]);
        __m256d a3 = _mm256_set1_pd(pa[2]);
        // Prefetch (pa+da,  _MM_HINT_T0);
        // Prefetch (pa+da,  _MM_HINT_T0);
        // Prefetch (pa+da+2,  _MM_HINT_T0);
        
        // Prefetch (pb+31,  _MM_HINT_T1);
        // Prefetch (pb+23,  _MM_HINT_T1);
        
        __m256d b1 = _mm256_loadu_pd(pb);
        sum11 -= a1 * b1;
        sum21 -= a2 * b1;
        sum31 -= a3 * b1;
        
        __m256d b2 = _mm256_loadu_pd(pb+4);
        sum12 -= a1 * b2;
        sum22 -= a2 * b2;
        sum32 -= a3 * b2;
        
        __m256d b3 = _mm256_loadu_pd(pb+8);
        sum13 -= a1 * b3;
        sum23 -= a2 * b3;
        sum33 -= a3 * b3;
        
        __m256d b4 = _mm256_loadu_pd(pb+12);
        sum14 -= a1 * b4;
        sum24 -= a2 * b4;
        sum34 -= a3 * b4;
      }

    _mm256_storeu_pd (pc, sum11);
    _mm256_storeu_pd (pc+4, sum12);
    _mm256_storeu_pd (pc+8, sum13);
    _mm256_storeu_pd (pc+12, sum14);
    pc += dc;
    _mm256_storeu_pd (pc, sum21);
    _mm256_storeu_pd (pc+4, sum22);
    _mm256_storeu_pd (pc+8, sum23);
    _mm256_storeu_pd (pc+12, sum24);
    pc += dc;
    _mm256_storeu_pd (pc, sum31);
    _mm256_storeu_pd (pc+4, sum32);
    _mm256_storeu_pd (pc+8, sum33);
    _mm256_storeu_pd (pc+12, sum34);
    pc += dc;
  }
  
  INLINE void MyScal3x4Trans (size_t ninner,
                              double * pa, size_t da,
                              double * pb, size_t db,
                              double * pc, size_t dc)
  {
    __assume (ninner > 0);        
    double * hpc = pc;
    SIMD<double> sum11 = _mm256_loadu_pd(pc);
    pc += dc;
    SIMD<double> sum21 = _mm256_loadu_pd(pc);
    pc += dc;
    SIMD<double> sum31 = _mm256_loadu_pd(pc);
    pc += dc;
    pc = hpc;
    
#pragma nounroll    
    for (size_t i = 0; i < ninner; i++, pa += da, pb += db)
      {
        __m256d a1 = _mm256_set1_pd(*pa);
        __m256d a2 = _mm256_set1_pd(pa[1]);
        __m256d a3 = _mm256_set1_pd(pa[2]);
        
        __m256d b1 = _mm256_loadu_pd(pb);
        sum11 -= a1 * b1;
        sum21 -= a2 * b1;
        sum31 -= a3 * b1;
      }

    _mm256_storeu_pd (pc, sum11.Data());
    pc += dc;
    _mm256_storeu_pd (pc, sum21.Data());
    pc += dc;
    _mm256_storeu_pd (pc, sum31.Data());
    pc += dc;
  }
  
  INLINE void MyScal3x4Trans (size_t ninner,
                              double * pa, size_t da,
                              double * pb, size_t db,
                              double * pc, size_t dc,
                              __m256i mask)
  {
    __assume (ninner > 0);    
    double * hpc = pc;
    SIMD<double> sum11 = _mm256_maskload_pd(pc, mask);
    pc += dc;
    SIMD<double> sum21 = _mm256_maskload_pd(pc, mask);
    pc += dc;
    SIMD<double> sum31 = _mm256_maskload_pd(pc, mask);
    pc += dc;
    pc = hpc;
    
#pragma nounroll    
    for (size_t i = 0; i < ninner; i++, pa += da, pb += db)
      {
        __m256d a1 = _mm256_set1_pd(*pa);
        __m256d a2 = _mm256_set1_pd(pa[1]);
        __m256d a3 = _mm256_set1_pd(pa[2]);
        
        __m256d b1 = _mm256_maskload_pd(pb, mask);
        sum11 -= a1 * b1;
        sum21 -= a2 * b1;
        sum31 -= a3 * b1;
      }

    _mm256_maskstore_pd (pc, mask, sum11.Data());
    pc += dc;
    _mm256_maskstore_pd (pc, mask, sum21.Data());
    pc += dc;
    _mm256_maskstore_pd (pc, mask, sum31.Data());
    pc += dc;
  }

  INLINE void MyScal1x16Trans (size_t ninner,
                               double * pa, size_t da,
                               double * pb, size_t db,
                               double * pc, size_t dc)
  {
    __assume (ninner > 0);    
    SIMD<double> sum11 = _mm256_loadu_pd(pc);
    SIMD<double> sum12 = _mm256_loadu_pd(pc+4);
    SIMD<double> sum13 = _mm256_loadu_pd(pc+8);
    SIMD<double> sum14 = _mm256_loadu_pd(pc+12);
    
#pragma nounroll    
    for (size_t i = 0; i < ninner; i++, pa += da, pb += db)
      {
        __m256d a1 = _mm256_set1_pd(*pa);
        __m256d b1 = _mm256_loadu_pd(pb);
        sum11 -= a1 * b1;
        
        __m256d b2 = _mm256_loadu_pd(pb+4);
        sum12 -= a1 * b2;
        
        __m256d b3 = _mm256_loadu_pd(pb+8);
        sum13 -= a1 * b3;
        
        __m256d b4 = _mm256_loadu_pd(pb+12);
        sum14 -= a1 * b4;
      }

    _mm256_storeu_pd (pc, sum11.Data());
    _mm256_storeu_pd (pc+4, sum12.Data());
    _mm256_storeu_pd (pc+8, sum13.Data());
    _mm256_storeu_pd (pc+12, sum14.Data());
  }
  
  INLINE void MyScal1x4Trans (size_t ninner,
                              double * pa, size_t da,
                              double * pb, size_t db,
                              double * pc, size_t dc)
  {
    __assume (ninner > 0);    
    SIMD<double> sum11 = _mm256_loadu_pd(pc);
#pragma nounroll    
    for (size_t i = 0; i < ninner; i++, pa += da, pb += db)
      {
        __m256d a1 = _mm256_set1_pd(*pa);
        __m256d b1 = _mm256_loadu_pd(pb);
        sum11 -= a1 * b1;
      }
    _mm256_storeu_pd (pc, sum11.Data());
  }
  
  INLINE void MyScal1x4Trans (size_t ninner,
                              double * pa, size_t da,
                              double * pb, size_t db,
                              double * pc, size_t dc,
                              __m256i mask)
  {
    SIMD<double> sum11 = _mm256_maskload_pd(pc, mask);
    __assume (ninner > 0);
#pragma nounroll    
    for (size_t i = 0; i < ninner; i++, pa += da, pb += db)
      {
        __m256d a1 = _mm256_set1_pd(*pa);
        __m256d b1 = _mm256_maskload_pd(pb, mask);
        sum11 -= a1 * b1;
      }
    _mm256_maskstore_pd (pc, mask, sum11.Data());
  }
  

  INLINE void MyScalx16Trans (size_t ninner, size_t wa,
                              double * pa, size_t da,
                              double * pb, size_t db,
                              double * pc, size_t dc)
  {
    size_t i = 0;
    constexpr auto HINT = 1;
    for ( ; i+3 <= wa; i+= 3, pa += 3, pc += 3*dc)
      {
        Prefetch<HINT> (pc+3*dc);
        Prefetch<HINT> (pc+3*dc+8);
        Prefetch<HINT> (pc+3*dc+16);
        Prefetch<HINT> (pc+4*dc);
        Prefetch<HINT> (pc+4*dc+8);
        Prefetch<HINT> (pc+4*dc+16);
        Prefetch<HINT> (pc+5*dc);
        Prefetch<HINT> (pc+5*dc+8);
        Prefetch<HINT> (pc+5*dc+16);
        MyScal3x16Trans (ninner, pa, da, pb, db, pc, dc);
      }
    for ( ; i < wa; i+= 1, pa += 1, pc += dc)
      MyScal1x16Trans (ninner, pa, da, pb, db, pc, dc);
  }
  INLINE void MyScalx4Trans (size_t ninner, size_t wa,
                             double * pa, size_t da,
                             double * pb, size_t db,
                             double * pc, size_t dc)
  {
    size_t i = 0;
    for ( ; i+3 <= wa; i+= 3, pa += 3, pc += 3*dc)
      MyScal3x4Trans (ninner, pa, da, pb, db, pc, dc);
    for ( ; i < wa; i+= 1, pa += 1, pc += dc)
      MyScal1x4Trans (ninner, pa, da, pb, db, pc, dc);
  }
  INLINE void MyScalx4Trans (size_t ninner, size_t wa,
                             double * pa, size_t da,
                             double * pb, size_t db,
                             double * pc, size_t dc,
                             __m256i mask)
  {
    size_t i = 0;
    for ( ; i+3 <= wa; i+= 3, pa += 3, pc += 3*dc)
      MyScal3x4Trans (ninner, pa, da, pb, db, pc, dc, mask);
    for ( ; i < wa; i+= 1, pa += 1, pc += dc)
      MyScal1x4Trans (ninner, pa, da, pb, db, pc, dc, mask);
  }

  void SubAtB1 (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c)
  {
    size_t wa = a.Width();
    size_t wb = b.Width();
    size_t da = a.Dist();
    size_t db = b.Dist();
    size_t dc = c.Dist();
    size_t ninner = a.Height();

    size_t j = 0;
    double * pa = &a(0,0);
    double * pb = &b(0,0);
    double * pc = &c(0,0);
    for ( ; j+16 <= wb; j+=16, pb += 16, pc += 16)
      {
        /*
          constexpr auto hint = _MM_HINT_T1;
          Prefetch (pc+0*dc, hint);
          Prefetch (pc+0*dc+8, hint);
          Prefetch (pc+0*dc+16, hint);
          Prefetch (pc+1*dc, hint);
          Prefetch (pc+1*dc+8, hint);
          Prefetch (pc+1*dc+16, hint);
          Prefetch (pc+2*dc, hint);
          Prefetch (pc+2*dc+8, hint);
          Prefetch (pc+2*dc+16, hint);
        */
        MyScalx16Trans (ninner, wa, pa, da, pb, db, pc, dc);
      }
    for ( ; j+4 <= wb; j+=4, pb += 4, pc += 4)
      MyScalx4Trans (ninner, wa, pa, da, pb, db, pc, dc);
    if (j < wb)
      {
        __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(wb-j),
                                            _mm256_set_epi64x(3,2,1,0));
        MyScalx4Trans (ninner, wa, pa, da, pb, db, pc, dc, mask);
      }
  }

  void SubAtB2 (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c)
  {
    constexpr size_t bs = 96;
    size_t wa = a.Width();
    size_t i = 0;
    for ( ; i+bs < wa; i += bs)
      SubAtB1 (a.Cols(i,i+bs), b, c.Rows(i,i+bs));
    if (i < wa)
      SubAtB1 (a.Cols(i,wa), b, c.Rows(i,wa));    
  }
  
  void SubAtBVersion3x16 (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c)
  {
    /*
      if (a.Height() < 16)
      {
      c -= Trans(a)*b;
      return;
      }
      if (c.Height() < 16)
      {
      c -= Trans(a)*b;
      return;
      }
      if (c.Width() < 16)
      {
      c -= Trans(a)*b;
      return;
      }
    */
    // static Timer tsub("avector::SubAtb - lapack");
    // RegionTimer reg(tsub);
    // tsub.AddFlops (size_t(a.Height()) * c.Height() * c.Width());
    
    // c -= Trans(a)*b | Lapack;
    // return;
    // cout << c.Height() << " x " << c.Width() << " x " << a.Height() << " = " << size_t(c.Height())*c.Width()*a.Height() << endl;
    constexpr size_t bs = 32;
    size_t ha = a.Height();
    size_t i = 0;
    for ( ; i+bs < ha; i += bs)
      SubAtB2 (a.Rows(i,i+bs), b.Rows(i,i+bs), c);
    if (i < ha)
      SubAtB2 (a.Rows(i,ha), b.Rows(i,ha), c);    
  }

  /*
    void SubAtB (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c)
    {
    Matrix<> hc = c;
    SubAtB1 (a, b, hc);
    c -= Trans(a) * b | Lapack;
    double err = L2Norm(hc-c);
    if (err > 1e-5)
    {
    cout << "c-hc: " << endl << c-hc << endl;
    exit(1);
    }
    // cout << "err = " << L2Norm(hc-c) << endl;
    }
  */


#else // AVX
  
  void TransposeMatrix(SliceMatrix<> a, SliceMatrix<> b)
  {
    b = Trans(a);    
  }
  
#endif // AVX


  // ////////////////////////// begin SubAtB  Version 4x12 ///////////////  

  
  constexpr size_t NA = 128;
  [[maybe_unused]]
  constexpr size_t NB = 96;
  constexpr size_t NK = 128;
  
#if defined (__AVX__) && !defined(__AVX512F__)

  // prefetch a row-major matrix
  void PreFetchMatrix (size_t h, size_t w, double * p, size_t dist)
  {
#pragma nounroll  
    for (size_t i = 0; i < h; i++, p += dist)
      for (size_t j = 0; j < w+8; j+=8) // cache line of 64 B i.e. 8 doubles
        Prefetch<1> (p+j);
  }

  // copy matrix (width is multiple of 4)
  void CopyMatrixIn (size_t h, size_t w,
                     double * ps, size_t dists,
                     double * pd, size_t distd)
  {
    __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(w&3),
                                        _mm256_set_epi64x(3,2,1,0));
    
    for (size_t i = 0; i < h; i++, pd += distd)
      {
        double * psnext = ps+dists;
        /*
          for (size_t j = 0; j+4 <= w; j+=4)
          {
          _mm256_storeu_pd(pd+j, _mm256_loadu_pd(ps+j));
          Prefetch (psnext+j,  _MM_HINT_T1);
          }
        */
        size_t j = 0;
        for ( ; j + 16 <= w; j+=16)
          {
            auto val1 = _mm256_loadu_pd(ps+j);
            auto val2 = _mm256_loadu_pd(ps+j+4);
            auto val3 = _mm256_loadu_pd(ps+j+8);
            auto val4 = _mm256_loadu_pd(ps+j+12);
            _mm256_storeu_pd (pd+j, val1);
            _mm256_storeu_pd (pd+j+4, val2);
            _mm256_storeu_pd (pd+j+8, val3);
            _mm256_storeu_pd (pd+j+12, val4);
            Prefetch<1> (psnext+j);
            Prefetch<1> (psnext+j+8);
          }
        for ( ; j +4 <= w; j+=4)
          _mm256_storeu_pd(pd+j, _mm256_loadu_pd(ps+j));
        _mm256_maskstore_pd (pd+j, mask, _mm256_maskload_pd(ps+j, mask));
        
        // for ( ; j < w; j++)
        // pd[j] = ps[j];
        // _mm256_storeu_pd(pd+j, _mm256_loadu_pd(ps+j));
      
        ps = psnext;
      }
  }

  // witdth is 12
  void CopyMatrixIn12 (size_t h, 
                       double * ps, size_t dists,
                       double * pd, size_t distd)
  {
    for (size_t i = 0; i < h; i++, ps += dists, pd += distd)
      {
        __m256d val0 = _mm256_loadu_pd (ps);
        __m256d val1 = _mm256_loadu_pd (ps+4);
        __m256d val2 = _mm256_loadu_pd (ps+8);
        _mm256_store_pd (pd, val0);
        _mm256_store_pd (pd+4, val1);
        _mm256_store_pd (pd+8, val2);
        Prefetch<1> (ps+2*dists);
        Prefetch<1> (ps+2*dists+8);
        // Prefetch (ps+4*dists, _MM_HINT_T2);
        // Prefetch (ps+4*dists+8, _MM_HINT_T2);
      }
  }

  
  // copy matrix (width is multiple of 4)
  void CopyMatrixOut (size_t h, size_t w,
                      double * ps, size_t dists,
                      double * pd, size_t distd)
  {
    __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(w&3),
                                        _mm256_set_epi64x(3,2,1,0));
    
    for (size_t i = 0; i < h; i++, pd += distd)
      {
        double * psnext = ps+dists;
        size_t j = 0;
        for ( ; j + 16 <= w; j+=16)
          {
            auto val1 = _mm256_loadu_pd(ps+j);
            auto val2 = _mm256_loadu_pd(ps+j+4);
            auto val3 = _mm256_loadu_pd(ps+j+8);
            auto val4 = _mm256_loadu_pd(ps+j+12);

            _mm256_storeu_pd (pd+j, val1);
            _mm256_storeu_pd (pd+j+4, val2);
            _mm256_storeu_pd (pd+j+8, val3);
            _mm256_storeu_pd (pd+j+12, val4);
            /*
              _mm256_stream_pd (pd+j, val1);
              _mm256_stream_pd (pd+j+4, val1);
              _mm256_stream_pd (pd+j+8, val1);
              _mm256_stream_pd (pd+j+12, val1);
            */
          }
        for ( ; j +4 <= w; j+=4)
          _mm256_storeu_pd(pd+j, _mm256_loadu_pd(ps+j));
        _mm256_maskstore_pd (pd+j, mask, _mm256_maskload_pd(ps+j, mask));        
        // for ( ; j < w; j++)
        // pd[j] = ps[j];
        ps = psnext;
      }
  }


  void CopyMatrixInScaleRows (size_t h, size_t w,
                              double * ps, size_t dists,
                              double * pd, size_t distd,
                              double * pscale, size_t distscale)
  {
    __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(w&3),
                                        _mm256_set_epi64x(3,2,1,0));
    
    for (size_t i = 0; i < h; i++, pd += distd, pscale += distscale)
      {
        double * psnext = ps+dists;
        __m256d scale = _mm256_set1_pd(*pscale);
        /*
          for (size_t j = 0; j+4 <= w; j+=4)
          {
          _mm256_storeu_pd(pd+j, _mm256_loadu_pd(ps+j));
          Prefetch (psnext+j,  _MM_HINT_T1);
          }
        */
        size_t j = 0;
        for ( ; j + 16 <= w; j+=16)
          {
            auto val1 = scale * _mm256_loadu_pd(ps+j);
            auto val2 = scale * _mm256_loadu_pd(ps+j+4);
            auto val3 = scale * _mm256_loadu_pd(ps+j+8);
            auto val4 = scale * _mm256_loadu_pd(ps+j+12);
            _mm256_storeu_pd (pd+j, val1);
            _mm256_storeu_pd (pd+j+4, val2);
            _mm256_storeu_pd (pd+j+8, val3);
            _mm256_storeu_pd (pd+j+12, val4);
            Prefetch<1> (psnext+j);
            Prefetch<1> (psnext+j+8);
          }
        for ( ; j +4 <= w; j+=4)
          _mm256_storeu_pd(pd+j, scale * _mm256_loadu_pd(ps+j));
        _mm256_maskstore_pd (pd+j, mask, scale * _mm256_maskload_pd(ps+j, mask));
        
        // for ( ; j < w; j++)
        // pd[j] = ps[j];
        // _mm256_storeu_pd(pd+j, _mm256_loadu_pd(ps+j));
      
        ps = psnext;
      }
  }

  // transpose matrix of avx-vectors (i.e. don't shuffle)
  void CopyMatrixInVTrans (size_t h, size_t w,
                           double * ps, size_t dists,
                           double * pd, size_t distd)
  {
    __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(w&3),
                                        _mm256_set_epi64x(3,2,1,0));
    size_t i = 0;

    /*
    for ( ; i+3 < h; i+=4, pd += 16)
      {
        double * psnext = ps+4*dists;
        size_t j = 0;
        double * pd2 = pd;
        
        for ( ; j +4 <= w; j+=4, pd2 += 4*distd)
          {
            // Prefetch (ps+4*dists+j, _MM_HINT_T1);
            // Prefetch (ps+5*dists+j, _MM_HINT_T1);
            // Prefetch (ps+6*dists+j, _MM_HINT_T1);
            // Prefetch (ps+7*dists+j, _MM_HINT_T1);
            _mm256_store_pd(pd2, _mm256_loadu_pd(ps+j));
            _mm256_store_pd(pd2+4, _mm256_loadu_pd(ps+dists+j));
            _mm256_store_pd(pd2+8, _mm256_loadu_pd(ps+2*dists+j));
            _mm256_store_pd(pd2+12, _mm256_loadu_pd(ps+3*dists+j));
          }

        _mm256_maskstore_pd (pd2, mask, _mm256_maskload_pd(ps+j, mask));
        _mm256_maskstore_pd (pd2+4, mask, _mm256_maskload_pd(ps+dists+j, mask));
        ps = psnext;
      }
    */
    
    for ( ; i+1 < h; i+=2, pd += 8)
      {
        double * psnext = ps+2*dists;
        size_t j = 0;
        double * pd2 = pd;
        
        for ( ; j+4 <= w; j+=4, pd2 += 4*distd)
          {
            // Prefetch (psnext+j, _MM_HINT_T0);
            // Prefetch (psnext+dists+j, _MM_HINT_T0);
            _mm256_store_pd(pd2, _mm256_loadu_pd(ps+j));
            _mm256_store_pd(pd2+4, _mm256_loadu_pd(ps+dists+j));
          }

        _mm256_maskstore_pd (pd2, mask, _mm256_maskload_pd(ps+j, mask));
        _mm256_maskstore_pd (pd2+4, mask, _mm256_maskload_pd(ps+dists+j, mask));
        ps = psnext;
      }

    for ( /* size_t i = 0 */; i < h; i++, pd += 4)
      {
        double * psnext = ps+dists;
        size_t j = 0;
        double * pd2 = pd;
        /*
        for ( ; j + 16 <= w; j+=16, pd2 += 16*distd)
          {
            // Prefetch (ps+j+32, _MM_HINT_T0);
            // Prefetch (ps+j+40, _MM_HINT_T0);
            Prefetch (psnext+j, _MM_HINT_T0);
            Prefetch (psnext+j+8, _MM_HINT_T0);
            auto val1 = _mm256_loadu_pd(ps+j);
            auto val2 = _mm256_loadu_pd(ps+j+4);
            auto val3 = _mm256_loadu_pd(ps+j+8);
            auto val4 = _mm256_loadu_pd(ps+j+12);

            _mm256_store_pd (pd2, val1);
            _mm256_store_pd (pd2+4*distd, val2);
            _mm256_store_pd (pd2+8*distd, val3);
            _mm256_store_pd (pd2+12*distd, val4);
          }
        */
        Prefetch<0> (psnext+j);
        Prefetch<0> (psnext+j+8);            
        for ( ; j +4 <= w; j+=4, pd2 += 4*distd)
          _mm256_store_pd(pd2, _mm256_loadu_pd(ps+j));
          // _mm256_stream_pd(pd2, _mm256_loadu_pd(ps+j));
        _mm256_maskstore_pd (pd2, mask, _mm256_maskload_pd(ps+j, mask));
        
        ps = psnext;
      }
  }

  
  void CopyMatrixInVTransScaleRows (size_t h, size_t w,
                                    double * ps, size_t dists,
                                    double * pd, size_t distd,
                                    double * pscale, size_t distscale)
  {
    __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(w&3),
                                        _mm256_set_epi64x(3,2,1,0));

    size_t i = 0;
    
    for ( ; i < h; i+=2, pd += 8, pscale += 2*distscale)
      {
        double * psnext = ps+2*dists;
        __m256d scale0 = _mm256_set1_pd(*pscale);
        __m256d scale1 = _mm256_set1_pd(*(pscale+distscale));

        size_t j = 0;
        double * pd2 = pd;

        for ( ; j+4 <= w; j+=4, pd2 += 4*distd)
          {
            _mm256_storeu_pd(pd2, scale0 * _mm256_loadu_pd(ps+j));
            _mm256_storeu_pd(pd2+4, scale1 * _mm256_loadu_pd(ps+dists+j));
          }
        _mm256_maskstore_pd (pd2, mask, scale0 * _mm256_maskload_pd(ps+j, mask));
        _mm256_maskstore_pd (pd2+4, mask, scale1 * _mm256_maskload_pd(ps+dists+j, mask));
        ps = psnext;
      }
      
    for ( ; i < h; i++, pd += 4, pscale += distscale)
      {
        double * psnext = ps+dists;
        __m256d scale = _mm256_set1_pd(*pscale);

        size_t j = 0;
        double * pd2 = pd;
        for ( ; j + 16 <= w; j+=16, pd2 += 16*distd)
          {
            Prefetch<0> (psnext+j);
            Prefetch<0> (psnext+j+8);            
            auto val1 = scale * _mm256_loadu_pd(ps+j);
            auto val2 = scale * _mm256_loadu_pd(ps+j+4);
            auto val3 = scale * _mm256_loadu_pd(ps+j+8);
            auto val4 = scale * _mm256_loadu_pd(ps+j+12);

            _mm256_store_pd (pd2, val1);
            _mm256_store_pd (pd2+4*distd, val2);
            _mm256_store_pd (pd2+8*distd, val3);
            _mm256_store_pd (pd2+12*distd, val4);
          }
        Prefetch<0> (psnext+j);
        Prefetch<0> (psnext+j+8);            
        for ( ; j +4 <= w; j+=4, pd2 += 4*distd)
          _mm256_storeu_pd(pd2, scale * _mm256_loadu_pd(ps+j));
        _mm256_maskstore_pd (pd2, mask, scale * _mm256_maskload_pd(ps+j, mask));
        
        ps = psnext;
      }
  }
  


  void CopyMatrix (SliceMatrix<> source,
                   SliceMatrix<> dest)
  {
    CopyMatrixIn (source.Height(), source.Width(),
                  &source(0,0), source.Dist(),
                  &dest(0,0), dest.Dist());
  }


  // micro-kernels for A^t B
  // pa+i * pb+j   ->  C 4x12
  INLINE
  void KernelScal4x12Trans (size_t ninner,
                            double * pa, size_t da,
                            double * pb, size_t db,
                            double * pc, size_t dc)
  {
    __assume (ninner > 0);
#ifndef __GNUC__    
    double * hpc = pc;
    __m256d sum11 = _mm256_loadu_pd(pc);
    __m256d sum12 = _mm256_loadu_pd(pc+4);
    __m256d sum13 = _mm256_loadu_pd(pc+8);
    pc += dc;
    __m256d sum21 = _mm256_loadu_pd(pc);
    __m256d sum22 = _mm256_loadu_pd(pc+4);
    __m256d sum23 = _mm256_loadu_pd(pc+8);
    pc += dc;
    __m256d sum31 = _mm256_loadu_pd(pc);
    __m256d sum32 = _mm256_loadu_pd(pc+4);
    __m256d sum33 = _mm256_loadu_pd(pc+8);
    pc += dc;
    __m256d sum41 = _mm256_loadu_pd(pc);
    __m256d sum42 = _mm256_loadu_pd(pc+4);
    __m256d sum43 = _mm256_loadu_pd(pc+8);
    pc = hpc;
#else
    __m256d sum11, sum12, sum13;
    __m256d sum21, sum22, sum23;
    __m256d sum31, sum32, sum33;
    __m256d sum41, sum42, sum43;

    asm (
         "vmovupd (%[pc]), %[sum11]\n\t"
         "vmovupd 32(%[pc]), %[sum12]\n\t"
         "vmovupd 64(%[pc]), %[sum13]\n\t"

         "vmovupd (%[pc],%[dc8]), %[sum21]\n\t"
         "vmovupd 32(%[pc],%[dc8]), %[sum22]\n\t"
         "vmovupd 64(%[pc],%[dc8]), %[sum23]\n\t"

         "vmovupd (%[pc],%[dc8],2), %[sum31]\n\t"
         "vmovupd 32(%[pc],%[dc8],2), %[sum32]\n\t"
         "vmovupd 64(%[pc],%[dc8],2), %[sum33]\n\t"

         "vmovupd (%[pc],%[dc3],8), %[sum41]\n\t"
         "vmovupd 32(%[pc],%[dc3],8), %[sum42]\n\t"
         "vmovupd 64(%[pc],%[dc3],8), %[sum43]\n\t"
         :
         [sum11] "+x" (sum11), [sum21] "+x" (sum21), [sum31] "+x" (sum31), [sum41] "+x" (sum41),
         [sum12] "+x" (sum12), [sum22] "+x" (sum22), [sum32] "+x" (sum32), [sum42] "+x" (sum42),
         [sum13] "+x" (sum13), [sum23] "+x" (sum23), [sum33] "+x" (sum33), [sum43] "+x" (sum43)
         : [pc] "r" (pc), [dc8]  "r" (8*dc), [dc3] "r" (3*dc)
         );
#endif
  
#pragma nounroll
    // #pragma unroll 2
    for (size_t i = 0; i < ninner; i++, pa += 4, pb += db)
      {
        __m256d b1 = _mm256_loadu_pd(pb);
        __m256d b2 = _mm256_loadu_pd(pb+4);
        __m256d b3 = _mm256_loadu_pd(pb+8);

        // Prefetch (pb+db,  _MM_HINT_T1);
        // Prefetch (pb+db+8,  _MM_HINT_T1);
        // Prefetch (pb+db+15,  _MM_HINT_T0);
        Prefetch<0> (pa+32);
        Prefetch<0> (pb+4*db);
        Prefetch<0> (pb+4*db+8);
        // Prefetch (pb+4*db+11,  _MM_HINT_T0);
        // Prefetch (pb+15,  _MM_HINT_T1);
        // Prefetch (pb+23,  _MM_HINT_T1);

        // Prefetch (pb,  _MM_HINT_T1);
        // Prefetch (pb+8,  _MM_HINT_T1);
        // Prefetch (pb+16,  _MM_HINT_T1);
        // Prefetch (pa+7,  _MM_HINT_T0);
        // Prefetch (pa+da,  _MM_HINT_T0);
        // Prefetch (pa+da+2,  _MM_HINT_T0);
        // Prefetch (pb+db+16,  _MM_HINT_T0);
        // Prefetch (pb+db+11,  _MM_HINT_T1);

        __m256d a1 = _mm256_set1_pd(*pa);
        sum11 -= a1 * b1;
        sum12 -= a1 * b2;
        sum13 -= a1 * b3;
        
        __m256d a2 = _mm256_set1_pd(pa[1]);
        sum21 -= a2 * b1;
        sum22 -= a2 * b2;
        sum23 -= a2 * b3;
        
        __m256d a3 = _mm256_set1_pd(pa[2]);
        sum31 -= a3 * b1;
        sum32 -= a3 * b2;
        sum33 -= a3 * b3;
        
        __m256d a4 = _mm256_set1_pd(pa[3]);
        sum41 -= a4 * b1;
        sum42 -= a4 * b2;
        sum43 -= a4 * b3;
      }

    _mm256_storeu_pd (pc, sum11);
    _mm256_storeu_pd (pc+4, sum12);
    _mm256_storeu_pd (pc+8, sum13);
    pc += dc;
    _mm256_storeu_pd (pc, sum21);
    _mm256_storeu_pd (pc+4, sum22);
    _mm256_storeu_pd (pc+8, sum23);
    pc += dc;
    _mm256_storeu_pd (pc, sum31);
    _mm256_storeu_pd (pc+4, sum32);
    _mm256_storeu_pd (pc+8, sum33);
    pc += dc;
    _mm256_storeu_pd (pc, sum41);
    _mm256_storeu_pd (pc+4, sum42);
    _mm256_storeu_pd (pc+8, sum43);
  }

  void KernelScal1x12Trans (size_t ninner,
                            double * pa, size_t da,
                            double * pb, size_t db,
                            double * pc, size_t dc)
  {
    __assume (ninner > 0);
#ifndef __GNUC__    
    __m256d sum11 = _mm256_loadu_pd(pc);
    __m256d sum12 = _mm256_loadu_pd(pc+4);
    __m256d sum13 = _mm256_loadu_pd(pc+8);
#else
    __m256d sum11, sum12, sum13;

    asm (
         "vmovupd (%[pc]), %[sum11]\n\t"
         "vmovupd 32(%[pc]), %[sum12]\n\t"
         "vmovupd 64(%[pc]), %[sum13]\n\t"
         :
         [sum11] "+x" (sum11),
         [sum12] "+x" (sum12),
         [sum13] "+x" (sum13)
         : [pc] "r" (pc)
         );
#endif
  
#pragma nounroll
    for (size_t i = 0; i < ninner; i++, pa += 4 /* da */, pb += db)
      {
        __m256d b1 = _mm256_loadu_pd(pb);
        __m256d b2 = _mm256_loadu_pd(pb+4);
        __m256d b3 = _mm256_loadu_pd(pb+8);

        __m256d a1 = _mm256_set1_pd(*pa);
        sum11 -= a1 * b1;
        sum12 -= a1 * b2;
        sum13 -= a1 * b3;
      }

    _mm256_storeu_pd (pc, sum11);
    _mm256_storeu_pd (pc+4, sum12);
    _mm256_storeu_pd (pc+8, sum13);
  }



  void KernelScal4x4Trans (size_t ninner,
                           double * pa, size_t da,
                           double * pb, size_t db,
                           double * pc, size_t dc)
  {
    __assume (ninner > 0);
#ifndef __GNUC__
    double * hpc = pc;    
    __m256d sum11 = _mm256_loadu_pd(pc);
    pc += dc;
    __m256d sum21 = _mm256_loadu_pd(pc);
    pc += dc;
    __m256d sum31 = _mm256_loadu_pd(pc);
    pc += dc;
    __m256d sum41 = _mm256_loadu_pd(pc);
    pc = hpc;
#else
    __m256d sum11;
    __m256d sum21;
    __m256d sum31;
    __m256d sum41;

    asm (
         "vmovupd (%[pc]), %[sum11]\n\t"
         "vmovupd (%[pc],%[dc8]), %[sum21]\n\t"
         "vmovupd (%[pc],%[dc8],2), %[sum31]\n\t"
         "vmovupd (%[pc],%[dc3],8), %[sum41]"
         : [sum11] "+x" (sum11), [sum21] "+x" (sum21), [sum31] "+x" (sum31), [sum41] "+x" (sum41)
         : [pc] "r" (pc), [dc8]  "r" (8*dc), [dc3] "r" (3*dc)
         );
#endif
  
#pragma nounroll
    for (size_t i = 0; i < ninner; i++, pa += 4 /* da */, pb += db)  // v-trans
      {
        __m256d b1 = _mm256_loadu_pd(pb);

        __m256d a1 = _mm256_set1_pd(*pa);
        sum11 -= a1 * b1;
        __m256d a2 = _mm256_set1_pd(pa[1]);
        sum21 -= a2 * b1;
        __m256d a3 = _mm256_set1_pd(pa[2]);
        sum31 -= a3 * b1;
        __m256d a4 = _mm256_set1_pd(pa[3]);
        sum41 -= a4 * b1;
      }

    _mm256_storeu_pd (pc, sum11);
    pc += dc;
    _mm256_storeu_pd (pc, sum21);
    pc += dc;
    _mm256_storeu_pd (pc, sum31);
    pc += dc;
    _mm256_storeu_pd (pc, sum41);
    pc += dc;
  }

  INLINE void KernelScal4x4Trans (size_t ninner,
                                  double * pa, size_t da,
                                  double * pb, size_t db,
                                  double * pc, size_t dc,
                                  __m256i mask)
  {
    __assume (ninner > 0);    
    double * hpc = pc;
    SIMD<double> sum11 = _mm256_maskload_pd(pc, mask);
    pc += dc;
    SIMD<double> sum21 = _mm256_maskload_pd(pc, mask);
    pc += dc;
    SIMD<double> sum31 = _mm256_maskload_pd(pc, mask);
    pc += dc;
    SIMD<double> sum41 = _mm256_maskload_pd(pc, mask);
    pc += dc;
    pc = hpc;
    
#pragma nounroll    
    for (size_t i = 0; i < ninner; i++, pa += 4 /* da */, pb += db)   // v-trans
      {
        __m256d a1 = _mm256_set1_pd(*pa);
        __m256d a2 = _mm256_set1_pd(pa[1]);
        __m256d a3 = _mm256_set1_pd(pa[2]);
        __m256d a4 = _mm256_set1_pd(pa[3]);
        
        __m256d b1 = _mm256_maskload_pd(pb, mask);
        sum11 -= a1 * b1;
        sum21 -= a2 * b1;
        sum31 -= a3 * b1;
        sum41 -= a4 * b1;        
      }

    _mm256_maskstore_pd (pc, mask, sum11.Data());
    pc += dc;
    _mm256_maskstore_pd (pc, mask, sum21.Data());
    pc += dc;
    _mm256_maskstore_pd (pc, mask, sum31.Data());
    pc += dc;
    _mm256_maskstore_pd (pc, mask, sum41.Data());
    pc += dc;
  }




  void KernelScal1x4Trans (size_t ninner,
                           double * pa, size_t da,
                           double * pb, size_t db,
                           double * pc, size_t dc)
  {
    __assume (ninner > 0);
    __m256d sum11 = _mm256_loadu_pd(pc);
  
#pragma nounroll
    for (size_t i = 0; i < ninner; i++, pa += 4 /* da */, pb += db)   // v-trans
      {
        __m256d b1 = _mm256_loadu_pd(pb);
        __m256d a1 = _mm256_set1_pd(*pa);
        sum11 -= a1 * b1;
      }

    _mm256_storeu_pd (pc, sum11);
  }

  INLINE void KernelScal1x4Trans (size_t ninner,
                                  double * pa, size_t da,
                                  double * pb, size_t db,
                                  double * pc, size_t dc,
                                  __m256i mask)
  {
    SIMD<double> sum11 = _mm256_maskload_pd(pc, mask);
    __assume (ninner > 0);
#pragma nounroll    
    for (size_t i = 0; i < ninner; i++, pa += 4 /* da */, pb += db)    // v-trans
      {
        __m256d a1 = _mm256_set1_pd(*pa);
        __m256d b1 = _mm256_maskload_pd(pb, mask);
        sum11 -= a1 * b1;
      }
    _mm256_maskstore_pd (pc, mask, sum11.Data());
  }

  // INLINE
  void KernelScalNx12Trans (
                            double * pa, size_t da,
                            double * pb, size_t db,
                            double * pc, size_t dc,
                            size_t ninner, size_t na
                            )
  {
    size_t i = 0;

    alignas (64) double memb[NK*12];
    CopyMatrixIn12 (ninner, pb, db, &memb[0], 12);
    // PreFetchMatrix (ninner, 12, pb+12, db);
    
    for ( ; i+4 <= na; i+=4, pc += 4*dc, pa += 4*da)    // v-trans
      {
        double * nextpc = pc+4*dc;
        for (size_t j = 0; j < 4; j++, nextpc+=dc)
          {
            Prefetch<0> (nextpc);
            Prefetch<0> (nextpc+8);
          }
        // KernelScal4x12Trans (ninner, pa, da, pb, db, pc, dc);
        KernelScal4x12Trans (ninner, pa, da, memb, 12, pc, dc);
      }
    for ( ; i+1 <= na; i+=1, pc += dc, pa += 1)         // v-trans
      // KernelScal1x12Trans (ninner, pa, da, pb, db, pc, dc);
      KernelScal1x12Trans (ninner, pa, da, memb, 12, pc, dc);
  }

  void KernelScalNx4Trans (
                            double * pa, size_t da,
                            double * pb, size_t db,
                            double * pc, size_t dc,
                            size_t ninner, size_t na
                            )
  {
    size_t i = 0;
    for ( ; i+4 <= na; i+=4, pc += 4*dc, pa += 4*da)    // v-trans
      KernelScal4x4Trans (ninner, pa, da, pb, db, pc, dc);
    for ( ; i+1 <= na; i+=1, pc += dc, pa += 1)         // v-trans
      KernelScal1x4Trans (ninner, pa, da, pb, db, pc, dc);
  }
  
  void KernelScalNx4Trans (
                            double * pa, size_t da,
                            double * pb, size_t db,
                            double * pc, size_t dc,
                            size_t ninner, size_t na,
                            __m256i mask)
  {
    size_t i = 0;
    for ( ; i+4 <= na; i+=4, pc += 4*dc, pa += 4*da)    // v-trans
      KernelScal4x4Trans (ninner, pa, da, pb, db, pc, dc, mask);
    for ( ; i+1 <= na; i+=1, pc += dc, pa += 1)         // v-trans
      KernelScal1x4Trans (ninner, pa, da, pb, db, pc, dc, mask);
  }

  void KernelScalNxMTrans (
                           double * pa, size_t da,
                           double * pb, size_t db,
                           double * pc, size_t dc,
                           size_t ninner, size_t na, size_t nb
                           )
  {
    size_t j = 0;
    for ( ; j+12 <= nb; j+=12 , pb += 12, pc += 12)
      KernelScalNx12Trans (pa, da, pb, db, pc, dc, ninner, na);
    for ( ; j+4 <= nb; j+=4, pb += 4, pc += 4)
      KernelScalNx4Trans (pa, da, pb, db, pc, dc, ninner, na);
    if (j < nb)
      {
        __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(nb-j),
                                            _mm256_set_epi64x(3,2,1,0));
        KernelScalNx4Trans (pa, da, pb, db, pc, dc, ninner, na, mask);
      }
  }

  
  INLINE
  void KernelScal4x12TransN (
                                    double * pa, size_t da,
                                    double * pb, size_t db,
                                    double * pc, size_t dc,
                                    size_t ninner, size_t nb
                                    )
  {
    size_t j = 0;
    for ( ; j+12 <= nb; j+=12, pb += 12, pc += 12)
      {
        /*
        // prefetching ????
        Prefetch (pc+15, _MM_HINT_T1);
        Prefetch (pc+23, _MM_HINT_T1);
        Prefetch (pc+dc+15, _MM_HINT_T1);
        Prefetch (pc+dc+23, _MM_HINT_T1);
        Prefetch (pc+2*dc+15, _MM_HINT_T1);
        Prefetch (pc+2*dc+23, _MM_HINT_T1);
        */
        KernelScal4x12Trans (ninner, pa, da, pb, db, pc, dc);
      }
    for ( ; j+4 <= nb; j+=4, pb += 4, pc += 4)
      KernelScal4x4Trans (ninner, pa, da, pb, db, pc, dc);
    if (j < nb)
      {
        __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(nb-j),
                                            _mm256_set_epi64x(3,2,1,0));
        KernelScal4x4Trans (ninner, pa, da, pb, db, pc, dc, mask);
      }
  }

  void KernelScal1x12TransN (
                             double * pa, size_t da,
                             double * pb, size_t db,
                             double * pc, size_t dc,
                             size_t ninner, size_t nb
                             )
  {
    size_t j = 0;
    for ( ; j+12 <= nb; j+=12, pb += 12, pc += 12)
      KernelScal1x12Trans (ninner, pa, da, pb, db, pc, dc);
    for ( ; j+4 <= nb; j+=4, pb += 4, pc += 4)
      KernelScal1x4Trans (ninner, pa, da, pb, db, pc, dc);
    if (j < nb)
      {
        __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(nb-j),
                                            _mm256_set_epi64x(3,2,1,0));
        KernelScal1x4Trans (ninner, pa, da, pb, db, pc, dc, mask);
      }
  }


  // ninner is now small (16) such that kernel runs at peek performance ...
  void KernelScal4x12TransNM (
                              double * pa, size_t da,
                              double * pb, size_t db,
                              double * pc, size_t dc,
                              size_t ninner, size_t na, size_t nb
                              )
  {
    // PreFetchMatrix (ninner, na, da, pa);
    // PreFetchMatrix (ninner, nb, db, pb);
    size_t i = 0;
    for ( ; i+4 <= na; i+=4, pc += 4*dc, pa += 4*da)    // v-trans
      KernelScal4x12TransN (pa, da, pb, db, pc, dc, ninner, nb);
    for ( ; i+1 <= na; i+=1, pc += dc, pa += 1)         // v-trans
      KernelScal1x12TransN (pa, da, pb, db, pc, dc, ninner, nb);

    /*
      alignas (64) double tempa[ninner*na+4];
      alignas (64) double tempb[ninner*nb+4];
      CopyMatrix (ninner, na, pa, da, &tempa[0], na);
      CopyMatrix (ninner, nb, pb, db, &tempb[0], nb);
      size_t i = 0;
      for ( ; i+4 <= na; i+=4, pc += 4*dc, pa += 4)
      KernelScal4x12TransN (&tempa[0], na, &tempb[0], nb, pc, dc, ninner, nb);
      for ( ; i+1 <= na; i+=1, pc += dc, pa += 1)
      KernelScal1x12TransN (&tempa[0], na, &tempb[0], nb, pc, dc, ninner, nb);
    */
  }



  // block-block mult (all matrices have fixed size)
  INLINE void BlockScal4x12Trans_BB_inline (
                                            double * pa, size_t da,
                                            double * pb, size_t db,
                                            double * pc, size_t dc,
                                            size_t na, size_t nb, size_t k
                                            )
  {
    /*
    // PreFetchMatrix (na, nb, dc, pc);  
    constexpr size_t bs = 16;
    size_t i = 0;
    for ( ; i+bs <= k; i+=bs, pa+=bs*da, pb+=bs*db)
    KernelScal4x12TransNM (pa, da, pb, db, pc, dc, bs, na, nb);
    if (i < k)
    KernelScal4x12TransNM (pa, da, pb, db, pc, dc, k-i, na, nb);
    */
    KernelScal4x12TransNM (pa, da, pb, db, pc, dc, k, na, nb);
    /*
    alignas(64) double tempc[NA*NB];
    CopyMatrixIn (na, nb, pc, dc, &tempc[0], NB);
    constexpr size_t bs = 128;
    size_t i = 0;
    for ( ; i+bs <= k; i+=bs, pa+=bs*4, pb+=bs*db)   // pa += bs*4*da
      KernelScal4x12TransNM (pa, da, pb, db, &tempc[0], NB, bs, na, nb);
    if (i < k)
      KernelScal4x12TransNM (pa, da, pb, db, &tempc[0], NB, k-i, na, nb);
    CopyMatrixOut (na, nb, &tempc[0], NB, pc, dc);
    */
  }

  void BlockScal4x12Trans_BB (
                              double * pa, size_t da,
                              double * pb, size_t db,
                              double * pc, size_t dc,
                              size_t na, size_t nb, size_t k
                              )
  {
    BlockScal4x12Trans_BB_inline (pa, da, pb, db, pc, dc, na, nb, k);
  }


  template <int NA, int NB, int K>
  void TBlockScal4x12Trans_BB (
                               double * pa, size_t da,
                               double * pb, size_t db,
                               double * pc, size_t dc
                               )
  {
    BlockScal4x12Trans_BB_inline (pa, da, pb, db, pc, dc, NA, NB, K);
  }

  template <> void TBlockScal4x12Trans_BB<64,48,64> (double * pa, size_t da,
                                                     double * pb, size_t db,
                                                     double * pc, size_t dc);

  void StdBlock  (
                  double * pa, size_t da,
                  double * pb, size_t db,
                  double * pc, size_t dc
                  )
  {
    // BlockScal4x12Trans_BB_inline (pa, da, pb, db, pc, dc, NA, NB, NK);
    alignas (64) double mema[NA*NK];
    CopyMatrixInVTrans (NK, NA, pa, da, &mema[0], NA);
    // BlockScal4x12Trans_BB_inline (&mema[0], NA, pb, db, pc, dc, NA, NB, NK);
  }






  // good: A is already fixed (128 x 128)
  void MySubAtB_BP (size_t na, size_t nb, size_t k,
                    double * pa, size_t da,
                    double * pb, size_t db,
                    double * pc, size_t dc)
  {
    SliceMatrix<> a(k, na, da, pa);
    SliceMatrix<> b(k, nb, db, pb);
    SliceMatrix<> c(na, nb, dc, pc);
    // c -= Trans(a)*b;
    // return;

    /*
      size_t i = 0;
      constexpr size_t bs = NB;
      for ( ; i+bs <= nb; i += bs, pb += bs, pc += bs)
      BlockScal4x12Trans_BB (pa, da, pb, db, pc, dc, na, bs, k);
      if (i < nb)
      BlockScal4x12Trans_BB (pa, da, pb, db, pc, dc, na, nb-i, k);    
    */
  
    alignas (64) double mema[NA*NK];
    CopyMatrixInVTrans (k, na, pa, da, &mema[0], NA);
    // SliceMatrix<> aloc(k,na, NA, &mema[0]);
    // aloc = a;

    KernelScalNxMTrans(mema, NA, pb, db, pc, dc,
                       k, na, nb);
    /*
    size_t i = 0;
    constexpr size_t bs = NB;
    for ( ; i+bs <= nb; i += bs, pb += bs, pc += bs)
      BlockScal4x12Trans_BB (mema, NA, pb, db, pc, dc, na, bs, k);
    if (i < nb)
      BlockScal4x12Trans_BB (mema, NA, pb, db, pc, dc, na, nb-i, k);  
    */
  }

  template <size_t na, size_t k>
  void MySubAtB_BP (size_t nb, 
                    double * pa, size_t da,
                    double * pb, size_t db,
                    double * pc, size_t dc)
  {
    /*
      SliceMatrix<> a(k, na, da, pa);
      SliceMatrix<> b(k, nb, db, pb);
      SliceMatrix<> c(na, nb, dc, pc);
      c -= Trans(a)*b;
      return;
    */
    /*
      size_t i = 0;
      constexpr size_t bs = NB;
      for ( ; i+bs <= nb; i += bs, pb += bs, pc += bs)
      BlockScal4x12Trans_BB_inline (pa, da, pb, db, pc, dc, na, bs,k);
      if (i < nb)
      BlockScal4x12Trans_BB (pa, da, pb, db, pc, dc, na, nb-i, k);    
    */
    alignas (64) double mema[na*k];
    CopyMatrixInVTrans (k, na, pa, da, &mema[0], na);

    size_t i = 0;
    constexpr size_t bs = NB;
    for ( ; i+bs <= nb; i += bs, pb += bs, pc += bs)
      BlockScal4x12Trans_BB_inline (mema, na, pb, db, pc, dc, na, bs,k);
    if (i < nb)
      BlockScal4x12Trans_BB_inline (mema, na, pb, db, pc, dc, na, nb-i, k);    
  }


  void MySubAtB_PM (size_t na, size_t nb, size_t k,
                    double * pa, size_t da,
                    double * pb, size_t db,
                    double * pc, size_t dc)
  {
    size_t i = 0;
    constexpr size_t bs = NK;
    for ( ; i+bs <= k; i += bs, pa += bs*da, pb += bs*db)
      MySubAtB_BP (na, nb, bs, pa, da, pb, db, pc, dc);
    if (i < k)
      MySubAtB_BP (na, nb, k-i, pa, da, pb, db, pc, dc);    
  }

  template <size_t na>
  void MySubAtB_PM (size_t nb, size_t k,
                    double * pa, size_t da,
                    double * pb, size_t db,
                    double * pc, size_t dc)
  {
    /*
      SliceMatrix<> a(k, na, da, pa);
      SliceMatrix<> b(k, nb, db, pb);
      SliceMatrix<> c(na, nb, dc, pc);
      c -= Trans(a)*b;
      return;
    */
    size_t i = 0;
    constexpr size_t bs = NK;
    for ( ; i+bs <= k; i += bs, pa += bs*da, pb += bs*db)
      MySubAtB_BP<na,bs> (nb, pa, da, pb, db, pc, dc);
    if (i < k)
      MySubAtB_BP (na, nb, k-i, pa, da, pb, db, pc, dc);    
  }


  // A ... k x na,  B ... k x nb,   C ... na x nb
  void SubAtB_MM (size_t na, size_t nb, size_t k,
                  double * pa, size_t da,
                  double * pb, size_t db,
                  double * pc, size_t dc)
  {
    size_t i = 0;
    constexpr size_t bs = NA;
    for ( ; i+bs < na; i += bs, pa += bs, pc += bs*dc)
      MySubAtB_PM<bs> (nb, k, pa, da, pb, db, pc, dc);
    // MySubAtB2 (bs, nb, k, pa, da, pb, db, pc, dc);
    if (i < na)
      MySubAtB_PM (na-i, nb, k, pa, da, pb, db, pc, dc);
  }
  
  void SubAtB (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c)
  {
    // c -= Trans(a) * b | Lapack;
    // return;
    SubAtB_MM (a.Width(), b.Width(), a.Height(),
               &a(0,0), a.Dist(), &b(0,0), b.Dist(), &c(0,0), c.Dist());
  }


  
  void MySubAtDB_BP (SliceMatrix<double> a,
                     SliceVector<double> diag,
                     SliceMatrix<double> b, SliceMatrix<double> c)
  {
    alignas (64) double mema[(NA+8)*NK];
    size_t na = a.Width();
    size_t nb = b.Width();
    size_t k = a.Height();
    
    CopyMatrixInVTransScaleRows (k, na,
                                 &a(0,0), a.Dist(), &mema[0], NA+8,
                                 &diag(0), diag.Dist());

    KernelScalNxMTrans(mema, NA+8, &b(0,0), b.Dist(), &c(0,0), c.Dist(),
                       k, na, nb);

    /*
    size_t i = 0;
    constexpr size_t bs = NB;
    for ( ; i+bs <= nb; i += bs)
      BlockScal4x12Trans_BB_inline (mema, NA+8, &b(size_t(0),i), b.Dist(), &c(size_t(0),i), c.Dist(), na, bs,k);
    if (i < nb)
      BlockScal4x12Trans_BB_inline (mema, NA+8, &b(size_t(0),i), b.Dist(), &c(size_t(0),i), c.Dist(), na, nb-i, k);    
    */
  }

#else

  void MySubAtDB_BP (SliceMatrix<double> a,
                     SliceVector<double> diag,
                     SliceMatrix<double> b, SliceMatrix<double> c)
  {
    alignas (64) double mema[NA*NK];
    FlatMatrix<> loca(a.Height(), a.Width(), mema);

    for (size_t i = 0; i < loca.Height(); i++)
      loca.Row(i) = diag(i)*a.Row(i);

    c -= Trans(loca) * b;
  }
  
#endif

  
  /*
  void MySubAtDB_PM (SliceMatrix<double> a,
                     SliceVector<double> diag,
                     SliceMatrix<double> b, SliceMatrix<double> c)
  {
    size_t k = a.Height();
    size_t i = 0;
    constexpr size_t bs = NK;
    for ( ; i+bs <= k; i += bs) 
      MySubAtDB_BP (a.Rows(i,i+bs), diag.Range(i,i+bs), b.Rows(i,i+bs), c);
    if (i < k)
      MySubAtDB_BP (a.Rows(i,k), diag.Range(i,k), b.Rows(i,k), c);      
  }

  void SubAtDB (SliceMatrix<double> a,
                SliceVector<double> diag,
                SliceMatrix<double> b, SliceMatrix<double> c)
  {
    size_t na = a.Width();
    size_t i = 0;
    constexpr size_t bs = NA;
    for ( ; i+bs <= na; i += bs)
      MySubAtDB_PM (a.Cols(i,i+bs), diag, b, c.Rows(i,i+bs));
    if (i < na)
      MySubAtDB_PM (a.Cols(i,na), diag, b, c.Rows(i,na));
  }
  */
  


  // ////////////////////////// end SubAtB  Version 4x12 ///////////////




  
#ifdef __AVX__
  
  
  // mat-mat product
  /*
  // b.Width <= 4
  INLINE
  void MultMatMat4(SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(b.Width()),
                                        _mm256_set_epi64x(3, 2, 1, 0));

    unsigned int da = a.Dist();
    int wa = a.Width();
    int r = 0;
    double * bpc = &c(0,0);
    unsigned int dc = c.Dist();
    double * ar = &a(0,0);
    for ( ; r+7 < a.Height(); r+=8)
      {
        __m256d sum1 = _mm256_setzero_pd();
        __m256d sum2 = _mm256_setzero_pd();
        __m256d sum3 = _mm256_setzero_pd();
        __m256d sum4 = _mm256_setzero_pd();
        __m256d sum5 = _mm256_setzero_pd();
        __m256d sum6 = _mm256_setzero_pd();
        __m256d sum7 = _mm256_setzero_pd();
        __m256d sum8 = _mm256_setzero_pd();

        for (int j = 0; j < wa; j++)
          {
            __m256d rb =  _mm256_blendv_pd(_mm256_setzero_pd(),
                                           _mm256_loadu_pd(&b(j,0)),
                                           _mm256_castsi256_pd(mask));
            double * arj = ar + j;
            double * arj4 = arj + 4*da;
            sum1 += _mm256_set1_pd(*arj) * rb;
            sum2 += _mm256_set1_pd(*(arj+da)) * rb;
            sum3 += _mm256_set1_pd(*(arj+2*da)) * rb;
            sum4 += _mm256_set1_pd(*(arj+3*da)) * rb;
            sum5 += _mm256_set1_pd(*(arj4)) * rb;
            sum6 += _mm256_set1_pd(*(arj4+da)) * rb;
            sum7 += _mm256_set1_pd(*(arj4+2*da)) * rb;
            sum8 += _mm256_set1_pd(*(arj4+3*da)) * rb;
          }

        _mm256_maskstore_pd(bpc, mask, sum1);
        _mm256_maskstore_pd(bpc+dc, mask, sum2);
        bpc += 2*dc;
        _mm256_maskstore_pd(bpc, mask, sum3);
        _mm256_maskstore_pd(bpc+dc, mask, sum4);
        bpc += 2*dc;
        _mm256_maskstore_pd(bpc, mask, sum5);
        _mm256_maskstore_pd(bpc+dc, mask, sum6);
        bpc += 2*dc;
        _mm256_maskstore_pd(bpc, mask, sum7);
        _mm256_maskstore_pd(bpc+dc, mask, sum8);
        bpc += 2*dc;
        ar += 8*da;
      }

    if (r+3 < a.Height())
      {
        __m256d sum1 = _mm256_setzero_pd();
        __m256d sum2 = _mm256_setzero_pd();
        __m256d sum3 = _mm256_setzero_pd();
        __m256d sum4 = _mm256_setzero_pd();
      
        for (int j = 0; j < wa; j++)
          {
            __m256d rb =  _mm256_blendv_pd(_mm256_setzero_pd(),
                                           _mm256_loadu_pd(&b(j,0)),
                                           _mm256_castsi256_pd(mask));

            double * arj = ar + j;
            sum1 += _mm256_set1_pd(*arj) * rb;
            sum2 += _mm256_set1_pd(*(arj+da)) * rb;
            sum3 += _mm256_set1_pd(*(arj+2*da)) * rb;
            sum4 += _mm256_set1_pd(*(arj+3*da)) * rb;
          }

        _mm256_maskstore_pd(bpc, mask, sum1);
        _mm256_maskstore_pd(bpc+dc, mask, sum2);
        bpc += 2*dc;
        _mm256_maskstore_pd(bpc, mask, sum3);
        _mm256_maskstore_pd(bpc+dc, mask, sum4);
        bpc += 2*dc;
        r += 4;
        ar += 4*da;
      }
    if (r < a.Height()-1)
      {
        __m256d sum1 = _mm256_setzero_pd();
        __m256d sum2 = _mm256_setzero_pd();
      
        for (int j = 0; j < wa; j++)
          {
            __m256d rb =  _mm256_blendv_pd(_mm256_setzero_pd(),
                                           _mm256_loadu_pd(&b(j,0)),
                                           _mm256_castsi256_pd(mask));
            double * arj = ar + j;
            sum1 += _mm256_set1_pd(*arj) * rb;
            sum2 += _mm256_set1_pd(*(arj+da)) * rb;
          }

        _mm256_maskstore_pd(bpc + 0*dc, mask, sum1);
        _mm256_maskstore_pd(bpc + 1*dc, mask, sum2);
        bpc += 2*dc;
        r += 2;
        ar += 2*da;
      }

    if (r < a.Height())
      {
        __m256d sum = _mm256_setzero_pd();
        for (int j = 0; j < wa; j++)
          {
            __m256d rb =  _mm256_loadu_pd(&b(j,0));
            double * arj = ar + j;
            sum += _mm256_set1_pd(*arj) * rb;

          }

        _mm256_maskstore_pd(bpc + 0*dc, mask, sum);
      }
  }
  */

    /*
  INLINE
  void MultMatMat8(SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    unsigned int da = a.Dist();
    int wa = a.Width();
    int r = 0;
    double * bpc = &c(0,0);
    unsigned int dc = c.Dist();
    double * ar = &a(0,0);
    for ( ; r+3 < a.Height(); r+=4)
      {
        __m256d sum11 = _mm256_setzero_pd();
        __m256d sum21 = _mm256_setzero_pd();
        __m256d sum31 = _mm256_setzero_pd();
        __m256d sum41 = _mm256_setzero_pd();
        __m256d sum12 = _mm256_setzero_pd();
        __m256d sum22 = _mm256_setzero_pd();
        __m256d sum32 = _mm256_setzero_pd();
        __m256d sum42 = _mm256_setzero_pd();

        for (int j = 0; j < wa; j++)
          {
            __m256d rb1 = _mm256_loadu_pd(&b(j,0));
            __m256d rb2 = _mm256_loadu_pd(&b(j,4));

            double * arj = ar + j;

            __m256d a1 = _mm256_set1_pd(*arj);
            __m256d a2 = _mm256_set1_pd(*(arj+da));
            __m256d a3 = _mm256_set1_pd(*(arj+2*da));
            __m256d a4 = _mm256_set1_pd(*(arj+3*da));

            sum11 += a1 * rb1;
            sum21 += a2 * rb1;
            sum31 += a3 * rb1;
            sum41 += a4 * rb1;
            sum12 += a1 * rb2;
            sum22 += a2 * rb2;
            sum32 += a3 * rb2;
            sum42 += a4 * rb2;
          }

        _mm256_storeu_pd(bpc, sum11);
        _mm256_storeu_pd(bpc+4, sum12);
        bpc += dc;
        _mm256_storeu_pd(bpc, sum21);
        _mm256_storeu_pd(bpc+4, sum22);
        bpc += dc;
        _mm256_storeu_pd(bpc, sum31);
        _mm256_storeu_pd(bpc+4, sum32);
        bpc += dc;
        _mm256_storeu_pd(bpc, sum41);
        _mm256_storeu_pd(bpc+4, sum42);
        bpc += dc;
        ar += 4*da;
      }

    for ( ; r < a.Height()-1; r+=2)
      {
        __m256d sum11 = _mm256_setzero_pd();
        __m256d sum21 = _mm256_setzero_pd();
        __m256d sum12 = _mm256_setzero_pd();
        __m256d sum22 = _mm256_setzero_pd();

        for (int j = 0; j < wa; j++)
          {
            __m256d rb1 = _mm256_loadu_pd(&b(j,0));
            __m256d rb2 = _mm256_loadu_pd(&b(j,4));

            double * arj = ar + j;

            __m256d a1 = _mm256_set1_pd(*arj);
            __m256d a2 = _mm256_set1_pd(*(arj+da));

            sum11 += a1 * rb1;
            sum21 += a2 * rb1;
            sum12 += a1 * rb2;
            sum22 += a2 * rb2;
          }

        _mm256_storeu_pd(bpc, sum11);
        _mm256_storeu_pd(bpc+4, sum12);
        bpc += dc;
        _mm256_storeu_pd(bpc, sum21);
        _mm256_storeu_pd(bpc+4, sum22);
        bpc += dc;
        ar += 2*da;
      }

    for ( ; r < a.Height(); r++)
      {
        __m256d sum11 = _mm256_setzero_pd();
        __m256d sum12 = _mm256_setzero_pd();

        for (int j = 0; j < wa; j++)
          {
            __m256d rb1 = _mm256_loadu_pd(&b(j,0));
            __m256d rb2 = _mm256_loadu_pd(&b(j,4));

            double * arj = ar + j;

            __m256d a1 = _mm256_set1_pd(*arj);
            sum11 += a1 * rb1;
            sum12 += a1 * rb2;
          }

        _mm256_storeu_pd(bpc, sum11);
        _mm256_storeu_pd(bpc+4, sum12);
        bpc += dc;
        ar += da;
      }
  }
    */



#ifdef NONE
  // c = a * Diag (d)
  void MultMatDiagMat(AFlatMatrixD a, AFlatVectorD diag, AFlatMatrixD c)
  {
    /*
      for (int i = 0; i < diag.Size(); i++)
      c.Col(i) = diag(i) * a.Col(i);
    */
    int rest = 4*diag.VSize() - diag.Size();
    int loops = diag.VSize();
    if (rest) loops--;
  
    for (int i = 0; i < c.Height(); i++)
      for (int j = 0; j < loops; j++)
        c.Get(i,j) = a.Get(i,j) * diag.Get(j);

    if (rest)
      {
        __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(4-rest),
                                            _mm256_set_epi64x(3, 2, 1, 0));

        __m256d md = _mm256_maskload_pd((double*)&diag.Get(loops), mask);

        for (int i = 0; i < c.Height(); i++)
          {
            __m256d ma = _mm256_maskload_pd((double*)&a.Get(i,loops), mask);
            __m256d prod = md * ma;
            _mm256_maskstore_pd((double*)&c.Get(i,loops), mask, prod);
          }
      }
  }
#endif
  


  // }

#endif






  
}
