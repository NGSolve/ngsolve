#include <bla.hpp>
using namespace ngbla;


#if defined(__AVX__)



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





INLINE __m256d HAdd (__m256d v1, __m256d v2, __m256d v3, __m256d v4)
{
  __m256d hsum1 = _mm256_hadd_pd (v1, v2);
  __m256d hsum2 = _mm256_hadd_pd (v3, v4);
  
  __m256d hsum = 
    _mm256_add_pd (_mm256_insertf128_pd (hsum1, 
                                         _mm256_extractf128_pd (hsum2, 0), 1),
                   _mm256_insertf128_pd (hsum2, 
                                         _mm256_extractf128_pd (hsum1, 1), 0));
  return hsum;
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
      sum11 += _mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b1+i);
      sum12 += _mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b2+i);
      sum13 += _mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b3+i);
      sum14 += _mm256_loadu_pd(a1+i) * _mm256_loadu_pd(b4+i);
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
      sum11 += _mm256_maskload_pd (a1+i, mask) * _mm256_maskload_pd (b1+i, mask);
      sum12 += _mm256_maskload_pd (a1+i, mask) * _mm256_maskload_pd (b2+i, mask);      
      sum13 += _mm256_maskload_pd (a1+i, mask) * _mm256_maskload_pd (b3+i, mask);
      sum14 += _mm256_maskload_pd (a1+i, mask) * _mm256_maskload_pd (b4+i, mask);
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


  void AddABt (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c)
  {
    if (a.Width() <= 0) return;
    switch (a.Width() & 3)
      {
      case 0: AddABt_Rest<0>(a,b,c); return;
      case 1: AddABt_Rest<1>(a,b,c); return;
      case 2: AddABt_Rest<2>(a,b,c); return;
      case 3: AddABt_Rest<3>(a,b,c); return;
      }
  }









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
  
    for ( ; i < c.Height()-1; i += 2)
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
        for ( ; j < c.Width()-3; j += 4)
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

  void AddABtSym (AFlatMatrix<double> a, AFlatMatrix<double> b, SliceMatrix<double> c)
  {
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
    for ( ; j < c.Width()-3; j += 4)
      {
        int i = j;
        double * pc = &c(i,j);
        for ( ; i < c.Height()-3; i += 4)
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

        if (i < c.Height()-1)
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
  for ( ; j < c.Width()-3; j += 4)
    {
      int i = j;
      Complex * pc = &c(i,j);
      for ( ; i < c.Height()-3; i += 4)
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

      if (i < c.Height()-1)
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



    for ( ; i < c.Height()-3; i += 4)
      {
        int j = 0;
        for ( ; j < c.Width()-3; j += 4)
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
        
        for ( ; j < c.Width()-1; j += 2)
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


    
    for ( ; i < c.Height()-1; i += 2)
      {
        int j = 0;
        for ( ; j < c.Width()-3; j += 4)
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
        for ( ; j < c.Width()-1; j += 2)
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
        for ( ; j < c.Width()-3; j += 4)
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





// mat-mat product

// b.Width <= 4
INLINE
void MultMatMat4(SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
{
  __m256i mask = my_mm256_cmpgt_epi64(_mm256_set1_epi64x(b.Width()),
				      _mm256_set_epi64x(3, 2, 1, 0));

  /*
  __m256i mask;
  switch (b.Width())
    {
    case 1:
      mask = _mm256_set_epi64x(0,0,0,-1); break;
    case 2:
      mask = _mm256_set_epi64x(0,0,-1,-1); break;
    case 3:
      mask = _mm256_set_epi64x(0,-1,-1,-1); break;
    case 4:
      mask = _mm256_set_epi64x(-1,-1,-1,-1); break;
    }
  */
  unsigned int da = a.Dist();
  int wa = a.Width();
  int r = 0;
  double * bpc = &c(0,0);
  unsigned int dc = c.Dist();
  double * ar = &a(0,0);
  for ( ; r < a.Height()-7; r+=8)
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

  if (r < a.Height()-3)
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

// b.Width() = 8
INLINE
void MultMatMat8(SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
{
  unsigned int da = a.Dist();
  int wa = a.Width();
  int r = 0;
  double * bpc = &c(0,0);
  unsigned int dc = c.Dist();
  double * ar = &a(0,0);
  for ( ; r < a.Height()-3; r+=4)
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





// c = a * b
void MultMatMat(SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
{
  int k = 0;
  for ( ; k < b.Width()-7; k += 8)
    MultMatMat8(a, b.Cols(k,k+8), c.Cols(k,k+8));
  for ( ; k < b.Width(); k += 4)
    {
      int end = min2(b.Width(), k+4);
      MultMatMat4(a, b.Cols(k,end), c.Cols(k,end));
    }
}



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



// }

#endif
