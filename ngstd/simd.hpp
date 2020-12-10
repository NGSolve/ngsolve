#ifndef NG_SIMD_HPP
#define NG_SIMD_HPP

#include "simd_base.hpp"

#if (defined(_M_AMD64) || defined(_M_X64) || defined(__SSE__))
#ifndef __SSE__
#define __SSE__
#endif
#include "simd_sse.hpp"
#endif

#ifdef __AVX__
#include "simd_avx.hpp"
#endif

#ifdef __AVX512F__
#include "simd_avx512.hpp"
#endif

namespace ngstd
{
  INLINE auto HSum (SIMD<double,2> v1, SIMD<double,2> v2, SIMD<double,2> v3, SIMD<double,2> v4)
  {
    SIMD<double,2> hsum1 = my_mm_hadd_pd (v1.Data(), v2.Data());
    SIMD<double,2> hsum2 = my_mm_hadd_pd (v3.Data(), v4.Data());
    return SIMD<double,4> (hsum1, hsum2);
  }


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

  template<int N>
  INLINE auto HSum (SIMD<double,N> s1, SIMD<double,N> s2)
  {
    return SIMD<double,2>(HSum(s1), HSum(s2));
  }

  template<int N>
  INLINE auto HSum (SIMD<double,N> s1, SIMD<double,N> s2, SIMD<double,N> s3, SIMD<double,N> s4 )
  {
    return SIMD<double,4>(HSum(s1), HSum(s2), HSum(s3), HSum(s4));
  }

  INLINE auto GetMaskFromBits( unsigned int i )
  {
    return SIMD<mask64>::GetMaskFromBits(i);
  }

// For Netgen interface
#if defined __AVX512F__
    typedef __m512 tAVX;
    typedef __m512d tAVXd;
#elif defined __AVX__
    typedef __m256 tAVX;
    typedef __m256d tAVXd; 
#elif (defined(_M_AMD64) || defined(_M_X64) || defined(__SSE__))
    typedef __m128 tAVX;
    typedef __m128d tAVXd; 
#endif

}

#endif // NG_SIMD_HPP
