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
    // static Timer t("Complex MatVec"); RegionTimer reg(t);
    
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
    // static Timer t("Complex MatVec, ColMajor"); RegionTimer reg(t);
    
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
    // static Timer t("Complex MatVec, RowMajor, SliceVec"); RegionTimer reg(t);
    
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

        auto ysub = y.Range(i, i+SW);
        for (size_t k = 0; k < SW; k++)
          func(ysub(k), Complex(s_sumr[k], s_sumi[k]));
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
    NgGEMV_Short<2> (s, a, x, y,[](Complex & y, Complex sum) { y=sum; });
  }
  
  template <typename FUNC>
  void NgGEMV (Complex s, BareSliceMatrix<double,ColMajor> a, SliceVector<Complex> x, SliceVector<Complex> y, FUNC func) NETGEN_NOEXCEPT  
  {
    // static Timer t("Complex MatVec, ColMajor, SliceVec"); RegionTimer reg(t); t.AddFlops (2*x.Size()*y.Size());

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

        auto ysub = y.Range(i, i+MSW);
        Iterate<MSW> ([&] (auto k) {
          func(ysub(k), Complex(s_sumr[k], s_sumi[k]));
        });
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

        auto ysub = y.Range(i, i+SW);        
        for (size_t k = 0; k < SW; k++)
          func(ysub(k), Complex(s_sumr[k], s_sumi[k]));
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
    // static Timer t("Complex MatVec, ColMajor, SliceVec"); RegionTimer reg(t);
    // t.AddFlops (4*x.Size()*y.Size());
    
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





  /* *********************** GEMM Complex ******************************** */


  template <bool ADD, bool POS>
  void Func (Complex & y, Complex sum)
  {
    if constexpr (ADD)
      {
        if constexpr (POS)
          y += sum;
        else
          y -= sum;
      }
    else
      {
        if constexpr (POS)
          y = sum;
        else
          y = -sum;
      }
  }


  template <size_t AW, bool ADD, bool POS, ORDERING OA, ORDERING OB>  
  void NgGEMMFuncSmall (size_t ah, size_t bw,
                        BareSliceMatrix<Complex,OA> a, BareSliceMatrix<Complex,OB> b, BareSliceMatrix<Complex,RowMajor> c)
  {
    // static Timer t("NgGEMMC - small"); RegionTimer reg(t);
    size_t i = 0;
    for( ; i+2 <= bw; i+= 2)
      {
        Vec<AW,SIMD<double,4>> bi;
        for (size_t k = 0; k < AW; k++)
          {
            SIMD<double,2> bi0((double*)(void*)b.Addr(k,i));
            SIMD<double,2> bi1((double*)(void*)b.Addr(k,i+1));
            bi(k) = SIMD<double,4>(bi0,bi1);
          }

        for (size_t j = 0; j < ah; j++)
          {
            SIMD<double,4> sum(0.0);
            for (size_t k = 0; k < AW; k++)
              {
                SIMD<double,2> a0((double*)(void*)a.Addr(j,k));
                SIMD<double,4> aj(a0,a0);
                FMAComplex (aj,bi(k),sum);
              }
            Func<ADD,POS>(c(j,i), Complex(sum[0], sum[1]));
            Func<ADD,POS>(c(j,i+1), Complex(sum[2], sum[3]));
          }
      }
    for( ; i+1 <= bw; i+= 1)
      {
        Vec<AW,SIMD<double,2>> bi;
        for (size_t k = 0; k < AW; k++)
          {
            SIMD<double,2> bi0((double*)(void*)b.Addr(k,i));
            bi(k) = bi0;
          }

        for (size_t j = 0; j < ah; j++)
          {
            SIMD<double,2> sum(0.0);
            for (size_t k = 0; k < AW; k++)
              {
                SIMD<double,2> aj((double*)(void*)a.Addr(j,k));
                FMAComplex (aj,bi(k),sum);
              }
            Func<ADD,POS>(c(j,i), Complex(sum[0], sum[1]));
          }
      }
  }

  
  template <> NGS_DLL_HEADER pmatmatc<RowMajor,RowMajor> dispatch_matmatc<false,false,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<RowMajor,RowMajor> dispatch_matmatc<false,true,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<RowMajor,RowMajor> dispatch_matmatc<true,false,RowMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<RowMajor,RowMajor> dispatch_matmatc<true,true,RowMajor,RowMajor>[9];

  template <> NGS_DLL_HEADER pmatmatc<ColMajor,RowMajor> dispatch_matmatc<false,false, ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<ColMajor,RowMajor> dispatch_matmatc<false,true, ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<ColMajor,RowMajor> dispatch_matmatc<true,false, ColMajor,RowMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<ColMajor,RowMajor> dispatch_matmatc<true,true, ColMajor,RowMajor>[9];

  template <> NGS_DLL_HEADER pmatmatc<RowMajor,ColMajor> dispatch_matmatc<false,false,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<RowMajor,ColMajor> dispatch_matmatc<false,true,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<RowMajor,ColMajor> dispatch_matmatc<true,false,RowMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<RowMajor,ColMajor> dispatch_matmatc<true,true,RowMajor,ColMajor>[9];

  template <> NGS_DLL_HEADER pmatmatc<ColMajor,ColMajor> dispatch_matmatc<false,false, ColMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<ColMajor,ColMajor> dispatch_matmatc<false,true, ColMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<ColMajor,ColMajor> dispatch_matmatc<true,false, ColMajor,ColMajor>[9];
  template <> NGS_DLL_HEADER pmatmatc<ColMajor,ColMajor> dispatch_matmatc<true,true, ColMajor,ColMajor>[9];


  template <bool ADD, bool POS, ORDERING OA, ORDERING OB>  
  void NgGEMMFunc (size_t ah, size_t aw, size_t bw,
                   BareSliceMatrix<Complex,OA> a, BareSliceMatrix<Complex,OB> b, BareSliceMatrix<Complex,RowMajor> c);

  template <size_t AW, bool ADD, bool POS, ORDERING OA, ORDERING OB>  
  void NgGEMMFuncAW (size_t ah, size_t bw,
                     BareSliceMatrix<Complex,OA> a, BareSliceMatrix<Complex,OB> b, BareSliceMatrix<Complex,RowMajor> c)
  {
    NgGEMMFunc<ADD,POS,OA,OB> (ah, AW, bw, a, b, c);
  }
  
  template <bool ADD, bool POS, ORDERING OA, ORDERING OB>
  void InitMatMatC()
  {
    Iterate<std::size(dispatch_matmatc<ADD,POS,OA,OB>)> ([&] (auto i)
    {
#ifdef __aarch64__
      dispatch_matmatc<ADD,POS,OA,OB>[i] = &NgGEMMFuncSmall<i,ADD,POS,OA,OB>;
#else
      // x86: FMAComplex shuffles make the small kernels slow
      if constexpr (i.value == 0)
        dispatch_matmatc<ADD,POS,OA,OB>[i] = &NgGEMMFuncSmall<i,ADD,POS,OA,OB>;
      else
        dispatch_matmatc<ADD,POS,OA,OB>[i] = &NgGEMMFuncAW<i,ADD,POS,OA,OB>;
#endif
    });    
  }

  template <ORDERING OA, ORDERING OB>
  void InitMatMatC2()
  {
    InitMatMatC<false,false,OA,OB>();
    InitMatMatC<false,true,OA,OB>();
    InitMatMatC<true,false,OA,OB>();
    InitMatMatC<true,true,OA,OB>();
  }
  
  auto init_matmatc = [] ()
  {
    InitMatMatC2<RowMajor,RowMajor>();
    InitMatMatC2<RowMajor,ColMajor>();
    InitMatMatC2<ColMajor,RowMajor>();
    InitMatMatC2<ColMajor,ColMajor>();
    return 1;
  }();
  

  // complex GEMM: B block packed into panels of GW simds (re/im interleaved), GH x GW register tiles
  constexpr size_t SWC = SIMD<double>::Size();
  constexpr size_t CW = SWC/2;   // complex numbers per simd
#if defined(__aarch64__)
  constexpr size_t GH = 4, GW = 4;
#elif defined(__AVX512F__)
  constexpr size_t GH = 4, GW = 3;
#else
  constexpr size_t GH = 3, GW = 2;
#endif

  // panel p holds B(0:kb, p*GW*CW + ...) as bb[p*kb*GW + k*w + v], last panel narrower, zero padded
  template <ORDERING OB>
  void PackBC (size_t kb, size_t nb, BareSliceMatrix<Complex,OB> b, SIMD<double> * bb)
  {
    for (size_t j = 0; j < nb; j += GW*CW, bb += kb*GW)
      {
        size_t nj = min(GW*CW, nb-j);
        size_t pw = (nj+CW-1)/CW*CW;
        Complex * pbb = (Complex*)(void*)bb;
        if constexpr (OB == RowMajor)
          for (size_t k = 0; k < kb; k++)
            for (size_t q = 0; q < nj; q++)
              pbb[k*pw+q] = b(k,j+q);
        else
          for (size_t q = 0; q < nj; q++)
            for (size_t k = 0; k < kb; k++)
              pbb[k*pw+q] = b(k,j+q);
        for (size_t k = 0; k < kb; k++)
          for (size_t q = nj; q < pw; q++)
            pbb[k*pw+q] = 0.0;
      }
  }

  template <bool ADD, bool POS>
  INLINE void StoreC (SIMD<double> res, double * pc, size_t nvalid)
  {
    if (nvalid >= CW)
      {
        if constexpr (ADD)
          res = POS ? SIMD<double>(pc) + res : SIMD<double>(pc) - res;
        else if constexpr (!POS)
          res = -res;
        res.Store(pc);
      }
    else
      {
        SIMD<mask64> mask(int(2*nvalid));
        if constexpr (ADD)
          res = POS ? SIMD<double>(pc,mask) + res : SIMD<double>(pc,mask) - res;
        else if constexpr (!POS)
          res = -res;
        res.Store(pc, mask);
      }
  }

  // C(0:H, 0:ncols) op= A(0:H, 0:kb) * panel,  A(h,k) = pa[h*dah+k*dak]
  template <size_t H, size_t W, bool ADD, bool POS>
  INLINE void NgGEMMKernelC (size_t kb, const Complex * pa, size_t dah, size_t dak,
                             const SIMD<double> * pb, Complex * pc, size_t dc, size_t ncols)
  {
#ifdef __aarch64__
    // vcmla, no shuffles
    SIMD<double> sum[H][W];
    for (size_t h = 0; h < H; h++)
      for (size_t w = 0; w < W; w++)
        sum[h][w] = 0.0;
    for (size_t k = 0; k < kb; k++, pb += W, pa += dak)
      {
        SIMD<double> bk[W];
        for (size_t w = 0; w < W; w++) bk[w] = pb[w];
        for (size_t h = 0; h < H; h++)
          {
            SIMD<double> ah((double*)(void*)(pa+h*dah));
            for (size_t w = 0; w < W; w++)
              FMAComplex (ah, bk[w], sum[h][w]);
          }
      }
#else
    // sr += are*b, si += aim*b, combine once: c = sr + i*si
    SIMD<double> sr[H][W], si[H][W];
    for (size_t h = 0; h < H; h++)
      for (size_t w = 0; w < W; w++)
        sr[h][w] = si[h][w] = 0.0;
    for (size_t k = 0; k < kb; k++, pb += W, pa += dak)
      {
        SIMD<double> bk[W];
        for (size_t w = 0; w < W; w++) bk[w] = pb[w];
        for (size_t h = 0; h < H; h++)
          {
            SIMD<double> are(pa[h*dah].real()), aim(pa[h*dah].imag());
            for (size_t w = 0; w < W; w++)
              {
                sr[h][w] = FMA(are, bk[w], sr[h][w]);
                si[h][w] = FMA(aim, bk[w], si[h][w]);
              }
          }
      }
    double sgnmem[SWC];
    for (size_t q = 0; q < SWC; q++) sgnmem[q] = (q%2) ? 1.0 : -1.0;
    SIMD<double> sgn(&sgnmem[0]);
    SIMD<double> sum[H][W];
    for (size_t h = 0; h < H; h++)
      for (size_t w = 0; w < W; w++)
        sum[h][w] = FMA(sgn, SwapPairs(si[h][w]), sr[h][w]);
#endif
    for (size_t h = 0; h < H; h++)
      for (size_t w = 0; w < W; w++)
        StoreC<ADD,POS> (sum[h][w], (double*)(void*)(pc+h*dc+w*CW), ncols-w*CW);
  }

  template <size_t H, bool ADD, bool POS, ORDERING OA>
  INLINE void NgGEMMPanelsC (size_t kb, size_t nb, BareSliceMatrix<Complex,OA> a, size_t i,
                             Complex * mema, const SIMD<double> * bb, Complex * pc, size_t dc)
  {
    // col-major A: kb strided columns thrash L1, copy the tile k-major if B block is wide enough
    const Complex * pa = a.Addr(i,0);
    size_t dah = (OA == RowMajor) ? a.Dist() : 1;
    size_t dak = (OA == RowMajor) ? 1 : a.Dist();
    if (OA == ColMajor && nb >= 32)
      {
        for (size_t k = 0; k < kb; k++)
          for (size_t h = 0; h < H; h++)
            mema[k*H+h] = a(i+h,k);
        pa = mema;
        dah = 1;
        dak = H;
      }

    size_t j = 0;
    for ( ; j+GW*CW <= nb; j += GW*CW, bb += kb*GW)
      NgGEMMKernelC<H,GW,ADD,POS> (kb, pa, dah, dak, bb, pc+j, dc, GW*CW);
    if (j < nb)
      {
        size_t wl = (nb-j+CW-1)/CW;
        Iterate<GW> ([&] (auto w) {
          if (w.value+1 == wl)
            NgGEMMKernelC<H,w.value+1,ADD,POS> (kb, pa, dah, dak, bb, pc+j, dc, nb-j);
        });
      }
  }

  template <bool ADD, bool POS, ORDERING OA>
  void NgGEMMRowsC (size_t ah, size_t kb, size_t nb, BareSliceMatrix<Complex,OA> a,
                    Complex * mema, const SIMD<double> * bb, BareSliceMatrix<Complex,RowMajor> c)
  {
    size_t i = 0;
    for ( ; i+GH <= ah; i += GH)
      NgGEMMPanelsC<GH,ADD,POS> (kb, nb, a, i, mema, bb, c.Addr(i,0), c.Dist());
    if (i < ah)
      Iterate<GH-1> ([&] (auto r) {
        if (r.value+1 == ah-i)
          NgGEMMPanelsC<r.value+1,ADD,POS> (kb, nb, a, i, mema, bb, c.Addr(i,0), c.Dist());
      });
  }

  template <bool ADD, bool POS, ORDERING OA, ORDERING OB>  
  void NgGEMMFunc (size_t ah, size_t aw, size_t bw,
                   BareSliceMatrix<Complex,OA> a, BareSliceMatrix<Complex,OB> b, BareSliceMatrix<Complex,RowMajor> c)
  {
    // packed B block ~200 KB; wider for col-major A to amortize the A copy
    constexpr size_t KB = (OA == RowMajor) ? 128 : 64;
    constexpr size_t NB = ((OA == RowMajor) ? 96 : 192) / (GW*CW) * (GW*CW);
    alignas(64) SIMD<double> bb[KB*NB/CW];
    alignas(64) Complex mema[(OA == RowMajor) ? 1 : GH*KB];

    for (size_t j = 0; j < bw; j += NB)
      {
        size_t nb = min(NB, bw-j);
        for (size_t k = 0; k < aw; k += KB)
          {
            size_t kb = min(KB, aw-k);
            PackBC<OB> (kb, nb, b.Rows(k,k+kb).Cols(j,j+nb), &bb[0]);
            if (k == 0)
              NgGEMMRowsC<ADD,POS,OA> (ah, kb, nb, a.Cols(k,k+kb), mema, &bb[0], c.Cols(j,j+nb));
            else
              NgGEMMRowsC<true,POS,OA> (ah, kb, nb, a.Cols(k,k+kb), mema, &bb[0], c.Cols(j,j+nb));
          }
      }
  }
  
  template <>
  NGS_DLL_HEADER void NgGEMMBare<false,true>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,RowMajor> a, BareSliceMatrix<Complex,RowMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<false,true> (ah, aw, bw, a, b, c); }
  
  template <>
  NGS_DLL_HEADER void NgGEMMBare<false,true>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,ColMajor> a, BareSliceMatrix<Complex,RowMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<false,true> (ah, aw, bw, a, b, c); }
  
  template <>
  NGS_DLL_HEADER void NgGEMMBare<false,true>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,RowMajor> a, BareSliceMatrix<Complex,ColMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<false,true> (ah, aw, bw, a, b, c); }

  template <>
  NGS_DLL_HEADER void NgGEMMBare<false,true>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,ColMajor> a, BareSliceMatrix<Complex,ColMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<false,true> (ah, aw, bw, a, b, c); }


  template <>
  NGS_DLL_HEADER void NgGEMMBare<true,true>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,RowMajor> a, BareSliceMatrix<Complex,RowMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<true,true> (ah, aw, bw, a, b, c); }
  
  template <>
  NGS_DLL_HEADER void NgGEMMBare<true,true>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,ColMajor> a, BareSliceMatrix<Complex,RowMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<true,true> (ah, aw, bw, a, b, c); }
  
  template <>
  NGS_DLL_HEADER void NgGEMMBare<true,true>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,RowMajor> a, BareSliceMatrix<Complex,ColMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<true,true> (ah, aw, bw, a, b, c); }

  template <>
  NGS_DLL_HEADER void NgGEMMBare<true,true>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,ColMajor> a, BareSliceMatrix<Complex,ColMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<true,true> (ah, aw, bw, a, b, c); }




  template <>
  NGS_DLL_HEADER void NgGEMMBare<true,false>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,RowMajor> a, BareSliceMatrix<Complex,RowMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<true,false> (ah, aw, bw, a, b, c); }
  
  template <>
  NGS_DLL_HEADER void NgGEMMBare<true,false>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,ColMajor> a, BareSliceMatrix<Complex,RowMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<true,false> (ah, aw, bw, a, b, c); }
  
  template <>
  NGS_DLL_HEADER void NgGEMMBare<true,false>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,RowMajor> a, BareSliceMatrix<Complex,ColMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<true,false> (ah, aw, bw, a, b, c); }

  template <>
  NGS_DLL_HEADER void NgGEMMBare<true,false>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,ColMajor> a, BareSliceMatrix<Complex,ColMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<true,false> (ah, aw, bw, a, b, c); }


  

  template <>
  NGS_DLL_HEADER void NgGEMMBare<false,false>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,RowMajor> a, BareSliceMatrix<Complex,RowMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<false,false> (ah, aw, bw, a, b, c); }
  
  template <>
  NGS_DLL_HEADER void NgGEMMBare<false,false>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,ColMajor> a, BareSliceMatrix<Complex,RowMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<false,false> (ah, aw, bw, a, b, c); }
  
  template <>
  NGS_DLL_HEADER void NgGEMMBare<false,false>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,RowMajor> a, BareSliceMatrix<Complex,ColMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<false,false> (ah, aw, bw, a, b, c); }

  template <>
  NGS_DLL_HEADER void NgGEMMBare<false,false>(size_t ah, size_t aw, size_t bw,
                              BareSliceMatrix<Complex,ColMajor> a, BareSliceMatrix<Complex,ColMajor> b, BareSliceMatrix<Complex,RowMajor> c)
  { NgGEMMFunc<false,false> (ah, aw, bw, a, b, c); }



}

