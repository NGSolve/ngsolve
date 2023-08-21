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

  void SetVector (Complex val, FlatVector<Complex> dest) NETGEN_NOEXCEPT
  {
    for (size_t i = 0; i < dest.Size(); i++)
      dest[i] = val;
  }
  
  void SetVector (Complex val, SliceVector<Complex> dest) NETGEN_NOEXCEPT
  {
    for (size_t i = 0; i < dest.Size(); i++)
      dest[i] = val;
  }
  
  
  void CopyVector (BareVector<Complex> src, FlatVector<Complex> dest) NETGEN_NOEXCEPT
  {
    for (size_t i = 0; i < dest.Size(); i++)
      dest[i] = src[i];
  }
  
  void CopyVector (BareSliceVector<Complex> src, BareSliceVector<Complex> dest, size_t size) NETGEN_NOEXCEPT
  {
    for (size_t i = 0; i < size; i++)
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







  // **************** FlatVector ****************8

  template <typename TM, typename TV, typename FUNC>
  void NgGEMV (Complex s, BareSliceMatrix<TM,RowMajor> a, FlatVector<const TV> x, FlatVector<Complex> y, FUNC func) NETGEN_NOEXCEPT  
  {
    for (size_t i = 0; i < y.Size(); i++)
      {
        Complex sum = 0;
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), s*sum);
      }
  }
  template <typename TM, typename TV, typename FUNC>
  void NgGEMV (Complex s, BareSliceMatrix<TM,ColMajor> a, FlatVector<const TV> x, FlatVector<Complex> y, FUNC func) NETGEN_NOEXCEPT  
  {
    for (size_t i = 0; i < y.Size(); i++)
      {
        Complex sum = 0;
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), s*sum);
      }
  }

  template <>
  NGS_DLL_HEADER  
  void NgGEMV<false> (Complex s, BareSliceMatrix<double,RowMajor> a, FlatVector<const Complex> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<true> (Complex s, BareSliceMatrix<double,RowMajor> a, FlatVector<const Complex> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y+=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<false> (Complex s, BareSliceMatrix<double,ColMajor> a, FlatVector<const Complex> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<true> (Complex s, BareSliceMatrix<double,ColMajor> a, FlatVector<const Complex> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y+=sum; });
  }


  template <>
  NGS_DLL_HEADER    
  void NgGEMV<false> (Complex s, BareSliceMatrix<Complex,RowMajor> a, FlatVector<const double> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<true> (Complex s, BareSliceMatrix<Complex,RowMajor> a, FlatVector<const double> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y+=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<false> (Complex s, BareSliceMatrix<Complex,ColMajor> a, FlatVector<const double> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<true> (Complex s, BareSliceMatrix<Complex,ColMajor> a, FlatVector<const double> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y+=sum; });
  }


  
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<false> (Complex s, BareSliceMatrix<Complex,RowMajor> a, FlatVector<const Complex> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<true> (Complex s, BareSliceMatrix<Complex,RowMajor> a, FlatVector<const Complex> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y+=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<false> (Complex s, BareSliceMatrix<Complex,ColMajor> a, FlatVector<const Complex> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<true> (Complex s, BareSliceMatrix<Complex,ColMajor> a, FlatVector<const Complex> x, FlatVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y+=sum; });
  }


  // **************** SliceVector ****************8
  

  template <typename TM, typename FUNC>
  void NgGEMV (Complex s, BareSliceMatrix<TM,RowMajor> a, SliceVector<Complex> x, SliceVector<Complex> y, FUNC func) NETGEN_NOEXCEPT  
  {
    for (size_t i = 0; i < y.Size(); i++)
      {
        Complex sum = 0;
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), s*sum);
      }
  }
  template <typename TM, typename FUNC>
  void NgGEMV (Complex s, BareSliceMatrix<TM,ColMajor> a, SliceVector<Complex> x, SliceVector<Complex> y, FUNC func) NETGEN_NOEXCEPT  
  {
    for (size_t i = 0; i < y.Size(); i++)
      {
        Complex sum = 0;
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), s*sum);
      }
  }

  template <>
  NGS_DLL_HEADER    
  void NgGEMV<false> (Complex s, BareSliceMatrix<double,RowMajor> a, SliceVector<Complex> x, SliceVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<true> (Complex s, BareSliceMatrix<double,RowMajor> a, SliceVector<Complex> x, SliceVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y+=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<false> (Complex s, BareSliceMatrix<double,ColMajor> a, SliceVector<Complex> x, SliceVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<true> (Complex s, BareSliceMatrix<double,ColMajor> a, SliceVector<Complex> x, SliceVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y+=sum; });
  }

  template <>
  NGS_DLL_HEADER    
  void NgGEMV<false> (Complex s, BareSliceMatrix<Complex,RowMajor> a, SliceVector<Complex> x, SliceVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<true> (Complex s, BareSliceMatrix<Complex,RowMajor> a, SliceVector<Complex> x, SliceVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y+=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<false> (Complex s, BareSliceMatrix<Complex,ColMajor> a, SliceVector<Complex> x, SliceVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y=sum; });
  }
  template <>
  NGS_DLL_HEADER    
  void NgGEMV<true> (Complex s, BareSliceMatrix<Complex,ColMajor> a, SliceVector<Complex> x, SliceVector<Complex> y) NETGEN_NOEXCEPT
  {
    NgGEMV (s, a, x, y, [](Complex & y, Complex sum) { y+=sum; });
  }



  
}

