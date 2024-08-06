#include <bla.hpp>


// do we have 32 vector-registers ?
#if defined(__AVX512F__) || defined(__arm64__)
constexpr bool reg32 = true;
#else
constexpr bool reg32 = false;
#endif

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
    static Timer t("Complex MatVec"); RegionTimer reg(t);
    
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
    static Timer t("Complex MatVec, ColMajor"); RegionTimer reg(t);
    
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
    static Timer t("Complex MatVec, SliceVec"); RegionTimer reg(t);
    
    for (size_t i = 0; i < y.Size(); i++)
      {
        Complex sum = 0;
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), s*sum);
      }
  }


  
  template <size_t SX, typename FUNC>
  void NgGEMV_Short (Complex s, BareSliceMatrix<double,ColMajor> a, SliceVector<Complex> x, SliceVector<Complex> y, FUNC func) NETGEN_NOEXCEPT  
  {
    /*
    Vec<SX> hxr, hxi;
    for (size_t j = 0; j < SX; j++)
      {
        hxr(j) = x(j).real();
        hxi(j) = x(j).imag();
      }
    */
    size_t i = 0;
    constexpr size_t SW = SIMD<double>::Size();
    for ( ; i+SW <= y.Size(); i+= SW)
      {
        SIMD<double,SW> sumr = 0.0;
        SIMD<double,SW> sumi = 0.0;
        for (size_t j = 0; j < SX; j++)
          {
            SIMD<double,SW> aij(a.Addr(i,j));
            SIMD<double,SW> xr(x(j).real());
            SIMD<double,SW> xi(x(j).imag());
            sumr += aij*xr;
            sumi += aij*xi;
          }
        SIMD<double,SW> s_sumr = s.real()*sumr - s.imag()*sumi;
        SIMD<double,SW> s_sumi = s.real()*sumi + s.imag()*sumr;

        for (size_t k = 0; k < SW; k++)
          func(y(i+k), Complex(s_sumr[k], s_sumi[k]));
      }
           
    for ( ; i < y.Size(); i++)
      {
        Complex sum = 0;
        for (size_t j = 0; j < SX; j++)
          sum += a(i,j) * x(j);
        func(y(i), s*sum);
      }
  }

  void Test1 (Complex s, BareSliceMatrix<double,ColMajor> a, SliceVector<Complex> x, SliceVector<Complex> y)
  {
    NgGEMV_Short<2> (s, a, x, y,[](Complex & y, Complex sum) { y+=sum; });
  }
  
  template <typename FUNC>
  void NgGEMV (Complex s, BareSliceMatrix<double,ColMajor> a, SliceVector<Complex> x, SliceVector<Complex> y, FUNC func) NETGEN_NOEXCEPT  
  {
    // static Timer t("Complex MatVec, ColMajor, SliceVec"); RegionTimer reg(t);
    // t.AddFlops (2*x.Size()*y.Size());

    /*
    // dit not improve ...
    typedef void (*pmv)(Complex, BareSliceMatrix<double,ColMajor>, SliceVector<Complex>, SliceVector<Complex>, FUNC) NETGEN_NOEXCEPT;
    static pmv dispatch_matvec[12];
    static auto init_matvec = [] ()
    {
      Iterate<std::size(dispatch_matvec)> ([&] (auto i)
      { dispatch_matvec[i] = &NgGEMV_Short<i,FUNC>; });
      return 1;
    }();

    if (x.Size() < std::size(dispatch_matvec))
      {
        (*dispatch_matvec[x.Size()])(s, a, x, y, func);
        return;
      }
    */


    
    /*
    for (size_t i = 0; i < y.Size(); i++)
      {
        Complex sum = 0;
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), s*sum);
      }
    */

    /*
    size_t i = 0;
    for ( ; i+4 <= y.Size(); i+=4)
      {
        Complex sum0 = 0;
        Complex sum1 = 0;
        Complex sum2 = 0;
        Complex sum3 = 0;
        
        for (size_t j = 0; j < x.Size(); j++)
          {
            sum0 += a(i  ,j) * x(j);
            sum1 += a(i+1,j) * x(j);
            sum2 += a(i+2,j) * x(j);
            sum3 += a(i+3,j) * x(j);
          }
        
        func(y(i  ), s*sum0);
        func(y(i+1), s*sum1);
        func(y(i+2), s*sum2);
        func(y(i+3), s*sum3);
      }
    */
    /*
    size_t i = 0;
    for ( ; i+8 <= y.Size(); i+=8)
      {
        Complex sum0 = 0;
        Complex sum1 = 0;
        Complex sum2 = 0;
        Complex sum3 = 0;
        Complex sum4 = 0;
        Complex sum5 = 0;
        Complex sum6 = 0;
        Complex sum7 = 0;
        
        for (size_t j = 0; j < x.Size(); j++)
          {
            Complex xj = x(j);
            sum0 += a(i  ,j) * xj;
            sum1 += a(i+1,j) * xj;
            sum2 += a(i+2,j) * xj;
            sum3 += a(i+3,j) * xj;
            sum4 += a(i+4,j) * xj;
            sum5 += a(i+5,j) * xj;
            sum6 += a(i+6,j) * xj;
            sum7 += a(i+7,j) * xj;
          }
        
        func(y(i  ), s*sum0);
        func(y(i+1), s*sum1);
        func(y(i+2), s*sum2);
        func(y(i+3), s*sum3);
        func(y(i+4), s*sum4);
        func(y(i+5), s*sum5);
        func(y(i+6), s*sum6);
        func(y(i+7), s*sum7);
      }

    for ( ; i+2 <= y.Size(); i+=2)
      {
        Complex sum0 = 0;
        Complex sum1 = 0;

        for (size_t j = 0; j < x.Size(); j++)
          {
            Complex xj = x(j);
            sum0 += a(i  ,j) * xj;
            sum1 += a(i+1,j) * xj;
          }
        
        func(y(i  ), s*sum0);
        func(y(i+1), s*sum1);
      }

    for ( ; i < y.Size(); i++)
      {
        Complex sum = 0;
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), s*sum);
      }
    */

    constexpr int SW = SIMD<double>::Size();
    constexpr int MSW = (reg32 ? 8 : 4)*SIMD<double>::Size();    
    size_t i = 0;
    for ( ; i+MSW <= y.Size(); i+= MSW)
      {
        SIMD<double,MSW> sumr = 0.0;
        SIMD<double,MSW> sumi = 0.0;
        for (size_t j = 0; j < x.Size(); j++)
          {
            SIMD<double,MSW> aij(a.Addr(i,j));
            SIMD<double,MSW> xr(x(j).real());
            SIMD<double,MSW> xi(x(j).imag());
            sumr += aij*xr;
            sumi += aij*xi;
          }
        SIMD<double,MSW> s_sumr = s.real()*sumr - s.imag()*sumi;
        SIMD<double,MSW> s_sumi = s.real()*sumi + s.imag()*sumr;

        for (size_t k = 0; k < MSW; k++)
          func(y(i+k), Complex(s_sumr[k], s_sumi[k]));
      }
    
    for ( ; i+SW <= y.Size(); i+= SW)
      {
        SIMD<double,SW> sumr = 0.0;
        SIMD<double,SW> sumi = 0.0;
        for (size_t j = 0; j < x.Size(); j++)
          {
            SIMD<double,SW> aij(a.Addr(i,j));
            SIMD<double,SW> xr(x(j).real());
            SIMD<double,SW> xi(x(j).imag());
            sumr += aij*xr;
            sumi += aij*xi;
          }
        SIMD<double,SW> s_sumr = s.real()*sumr - s.imag()*sumi;
        SIMD<double,SW> s_sumi = s.real()*sumi + s.imag()*sumr;

        for (size_t k = 0; k < SW; k++)
          func(y(i+k), Complex(s_sumr[k], s_sumi[k]));
      }
           
    for ( ; i < y.Size(); i++)
      {
        Complex sum = 0;
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), s*sum);
      }
  }



  template <typename FUNC>
  void NgGEMV (Complex s, BareSliceMatrix<Complex,ColMajor> a, SliceVector<Complex> x, SliceVector<Complex> y, FUNC func) NETGEN_NOEXCEPT  
  {
    static Timer t("Complex MatVec, ColMajor, SliceVec"); RegionTimer reg(t);
    t.AddFlops (4*x.Size()*y.Size());
    
    /*
      // dit not improve ...
    typedef void (*pmv)(Complex, BareSliceMatrix<TM,ColMajor>, SliceVector<Complex>, SliceVector<Complex>, FUNC) NETGEN_NOEXCEPT;
    static pmv dispatch_matvec[12];
    static auto init_matvec = [] ()
    {
      Iterate<std::size(dispatch_matvec)> ([&] (auto i)
      { dispatch_matvec[i] = &NgGEMV_Short<TM,i,FUNC>; });
      return 1;
    }();

    if (x.Size() < std::size(dispatch_matvec))
      {
        (*dispatch_matvec[x.Size()])(s, a, x, y, func);
        return;
      }
    */
    /*
    for (size_t i = 0; i < y.Size(); i++)
      {
        Complex sum = 0;
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), s*sum);
      }
    */

    /*
    size_t i = 0;
    for ( ; i+4 <= y.Size(); i+=4)
      {
        Complex sum0 = 0;
        Complex sum1 = 0;
        Complex sum2 = 0;
        Complex sum3 = 0;
        
        for (size_t j = 0; j < x.Size(); j++)
          {
            sum0 += a(i  ,j) * x(j);
            sum1 += a(i+1,j) * x(j);
            sum2 += a(i+2,j) * x(j);
            sum3 += a(i+3,j) * x(j);
          }
        
        func(y(i  ), s*sum0);
        func(y(i+1), s*sum1);
        func(y(i+2), s*sum2);
        func(y(i+3), s*sum3);
      }
    */

    size_t i = 0;
    for ( ; i+8 <= y.Size(); i+=8)
      {
        Complex sum0 = 0;
        Complex sum1 = 0;
        Complex sum2 = 0;
        Complex sum3 = 0;
        Complex sum4 = 0;
        Complex sum5 = 0;
        Complex sum6 = 0;
        Complex sum7 = 0;
        
        for (size_t j = 0; j < x.Size(); j++)
          {
            Complex xj = x(j);
            sum0 += a(i  ,j) * xj;
            sum1 += a(i+1,j) * xj;
            sum2 += a(i+2,j) * xj;
            sum3 += a(i+3,j) * xj;
            sum4 += a(i+4,j) * xj;
            sum5 += a(i+5,j) * xj;
            sum6 += a(i+6,j) * xj;
            sum7 += a(i+7,j) * xj;
          }
        
        func(y(i  ), s*sum0);
        func(y(i+1), s*sum1);
        func(y(i+2), s*sum2);
        func(y(i+3), s*sum3);
        func(y(i+4), s*sum4);
        func(y(i+5), s*sum5);
        func(y(i+6), s*sum6);
        func(y(i+7), s*sum7);
      }

    for ( ; i+2 <= y.Size(); i+=2)
      {
        Complex sum0 = 0;
        Complex sum1 = 0;

        for (size_t j = 0; j < x.Size(); j++)
          {
            Complex xj = x(j);
            sum0 += a(i  ,j) * xj;
            sum1 += a(i+1,j) * xj;
          }
        
        func(y(i  ), s*sum0);
        func(y(i+1), s*sum1);
      }

    for ( ; i < y.Size(); i++)
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

