#ifndef NG_SIMD_HPP
#define NG_SIMD_HPP

#include "simd_base.hpp"
#include "simd_sse.hpp"
#include "simd_avx.hpp"

namespace ngstd
{
  INLINE auto HSum (SIMD<double,2> v1, SIMD<double,2> v2, SIMD<double,2> v3, SIMD<double,2> v4)
  {
    SIMD<double,2> hsum1 = my_mm_hadd_pd (v1.Data(), v2.Data());
    SIMD<double,2> hsum2 = my_mm_hadd_pd (v3.Data(), v4.Data());
    return SIMD<double,4> (hsum1, hsum2);
  }


  template<int N>
  INLINE auto Unpack (SIMD<double,N> a, SIMD<double,N> b)
  {
    auto [a1,b1] = Unpack(a.Lo(), b.Lo());
    auto [a2,b2] = Unpack(a.Hi(), b.Hi());
    return make_tuple(SIMD<double,N>{ a1, b2 },
                      SIMD<double,N>{ b2, b2 });
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

}

#endif // NG_SIMD_HPP
