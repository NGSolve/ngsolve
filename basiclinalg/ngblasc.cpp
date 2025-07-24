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
  void InitMatMatC()
  {
    Iterate<std::size(dispatch_matmatc<ADD,POS,OA,OB>)> ([&] (auto i)
    { dispatch_matmatc<ADD,POS,OA,OB>[i] = &NgGEMMFuncSmall<i,ADD,POS,OA,OB>; });    
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
  

  template <bool ADD, bool POS, ORDERING OA, ORDERING OB>  
  void NgGEMMFunc2 (size_t ah, size_t aw, size_t bw,
                    BareSliceMatrix<Complex,OA> a, BareSliceMatrix<Complex,OB> b, BareSliceMatrix<Complex,RowMajor> c)
  {
    size_t i = 0;
    for ( ; i+4 <= ah; i+=4)
      {
        size_t j = 0;
        for ( ; j+4 <= bw; j+=4)
          {
            SIMD<double,4> sum00 = 0;
            SIMD<double,4> sum10 = 0;
            SIMD<double,4> sum20 = 0;
            SIMD<double,4> sum30 = 0;
            SIMD<double,4> sum01 = 0;
            SIMD<double,4> sum11 = 0;
            SIMD<double,4> sum21 = 0;
            SIMD<double,4> sum31 = 0;
            for (size_t k = 0; k < aw; k++)
              {
                SIMD<double,2> a0_((double*)(void*)a.Addr(i+0,k));
                SIMD<double,2> a1_((double*)(void*)a.Addr(i+1,k));
                SIMD<double,2> a2_((double*)(void*)a.Addr(i+2,k));
                SIMD<double,2> a3_((double*)(void*)a.Addr(i+3,k));
                SIMD<double,2> b0_((double*)(void*)b.Addr(k,j+0));
                SIMD<double,2> b1_((double*)(void*)b.Addr(k,j+1));
                SIMD<double,2> b2_((double*)(void*)b.Addr(k,j+2));
                SIMD<double,2> b3_((double*)(void*)b.Addr(k,j+3));
                SIMD<double,4> a0(a0_, a0_);
                SIMD<double,4> a1(a1_, a1_);
                SIMD<double,4> a2(a2_, a2_);
                SIMD<double,4> a3(a3_, a3_);
                SIMD<double,4> b0(b0_, b1_);
                SIMD<double,4> b1(b2_, b3_);
                
                FMAComplex (a0,b0,sum00);
                FMAComplex (a1,b0,sum10);
                FMAComplex (a2,b0,sum20);
                FMAComplex (a3,b0,sum30);
                FMAComplex (a0,b1,sum01);
                FMAComplex (a1,b1,sum11);
                FMAComplex (a2,b1,sum21);
                FMAComplex (a3,b1,sum31);
              }
            Func<ADD,POS>(c(i+0,j  ), Complex(sum00[0], sum00[1]));
            Func<ADD,POS>(c(i+0,j+1), Complex(sum00[2], sum00[3]));  
            Func<ADD,POS>(c(i+0,j+2), Complex(sum01[0], sum01[1]));  
            Func<ADD,POS>(c(i+0,j+3), Complex(sum01[2], sum01[3]));  
            Func<ADD,POS>(c(i+1,j  ), Complex(sum10[0], sum10[1]));  
            Func<ADD,POS>(c(i+1,j+1), Complex(sum10[2], sum10[3]));  
            Func<ADD,POS>(c(i+1,j+2), Complex(sum11[0], sum11[1]));  
            Func<ADD,POS>(c(i+1,j+3), Complex(sum11[2], sum11[3]));

            Func<ADD,POS>(c(i+2,j  ), Complex(sum20[0], sum20[1]));
            Func<ADD,POS>(c(i+2,j+1), Complex(sum20[2], sum20[3]));  
            Func<ADD,POS>(c(i+2,j+2), Complex(sum21[0], sum21[1]));  
            Func<ADD,POS>(c(i+2,j+3), Complex(sum21[2], sum21[3]));  
            Func<ADD,POS>(c(i+3,j  ), Complex(sum30[0], sum30[1]));  
            Func<ADD,POS>(c(i+3,j+1), Complex(sum30[2], sum30[3]));  
            Func<ADD,POS>(c(i+3,j+2), Complex(sum31[0], sum31[1]));  
            Func<ADD,POS>(c(i+3,j+3), Complex(sum31[2], sum31[3]));
          }
        for ( ; j < bw; j++)
          {
            Complex sum0 = 0;
            Complex sum1 = 0;
            Complex sum2 = 0;
            Complex sum3 = 0;
            for (size_t k = 0; k < aw; k++)
              {
                sum0 += a(i+0,k)*b(k,j);
                sum1 += a(i+1,k)*b(k,j);
                sum2 += a(i+2,k)*b(k,j);
                sum3 += a(i+3,k)*b(k,j);
              }
            Func<ADD,POS>(c(i+0,j), sum0);  
            Func<ADD,POS>(c(i+1,j), sum1);          
            Func<ADD,POS>(c(i+2,j), sum2);  
            Func<ADD,POS>(c(i+3,j), sum3);          
          }
      }
    for ( ; i < ah; i++)
      for (size_t j = 0; j < bw; j++)
        {
          Complex sum = 0;
          for (size_t k = 0; k < aw; k++)
            sum += a(i,k)*b(k,j);
          Func<ADD,POS>(c(i,j), sum);          
        }

    

    /*
    for (size_t i = 0; i < ah; i++)
      for (size_t j = 0; j < bw; j++)
        {
          Complex sum = 0;
          for (size_t k = 0; k < aw; k++)
            sum += a(i,k)*b(k,j);
          Func<ADD,POS>(c(i,j), sum);          
        }
    */
  }

  template <bool ADD, bool POS, ORDERING OA, ORDERING OB>  
  void NgGEMMFunc (size_t ah, size_t aw, size_t bw,
                   BareSliceMatrix<Complex,OA> a, BareSliceMatrix<Complex,OB> b, BareSliceMatrix<Complex,RowMajor> c)
  {
    constexpr size_t BH=48;
    constexpr size_t BW=48;

    for (size_t i = 0; i < ah; i+= BH)
      {
        size_t imax = min(ah, i+BH);
        size_t jmax = min(aw, BW);
        NgGEMMFunc2<ADD,POS,OA,OB> (imax-i, jmax, bw, a.Rows(i,imax).Cols(0,jmax), b.Rows(0,jmax), c.Rows(i,imax));
        for (size_t j = BW; j < aw; j+= BW)
          {
            size_t jmax = min(aw, j+BW);
            NgGEMMFunc2<true,POS,OA,OB> (imax-i, jmax-j, bw, a.Rows(i,imax).Cols(j,jmax), b.Rows(j,jmax), c.Rows(i,imax));            
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

