#include <bla.hpp>


/*
// do we have 32 vector-registers ?
#if defined(__AVX512F__) || defined(__arm64__)
constexpr bool reg32 = true;
#else
constexpr bool reg32 = false;
#endif
*/

namespace ngbla
{
  
  void CopyVector (BareVector<Complex> src, FlatVector<Complex> dest) NETGEN_NOEXCEPT
  {
    for (size_t i = 0; i < dest.Size(); i++)
      dest[i] = src[i];
  }
  
  void CopyVector (BareSliceVector<Complex> src, SliceVector<Complex> dest) NETGEN_NOEXCEPT
  {
    for (size_t i = 0; i < dest.Size(); i++)
      dest[i] = src[i];
  }

  void CopyVector (Complex alpha, BareVector<Complex> src, FlatVector<Complex> dest) NETGEN_NOEXCEPT
  {
    for (size_t i = 0; i < dest.Size(); i++)
      dest[i] = alpha * src[i];
  }
  
  void CopyVector (Complex alpha, BareSliceVector<Complex> src, SliceVector<Complex> dest) NETGEN_NOEXCEPT
  {
    for (size_t i = 0; i < dest.Size(); i++)
      dest[i] = alpha * src[i];
  }

  
  void AddVector (Complex alpha, BareVector<Complex> src, FlatVector<Complex> dest) NETGEN_NOEXCEPT
  {
    for (size_t i = 0; i < dest.Size(); i++)
      dest[i] += alpha * src[i];
  }
  
  void AddVector (Complex alpha, BareSliceVector<Complex> src, SliceVector<Complex> dest) NETGEN_NOEXCEPT
  {
    for (size_t i = 0; i < dest.Size(); i++)
      dest[i] += alpha * src[i];
  }

}

