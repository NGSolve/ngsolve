#ifndef FILE_NGBLAS
#define FILE_NGBLAS

// optimized matrix kernels

#if defined(__clang__ ) && defined(NETGEN_ARCH_AMD64)
#define REGCALL __regcall
#else
#define REGCALL
#endif

/*
namespace ngcore
{
  template <int S>
  INLINE auto Range (IC<S> s)
  {
    return IntRange(0,s);
  }
}
*/


namespace ngbla
{

  // ***************************** vector operations **************************

  
  extern NGS_DLL_HEADER void SetVector (double val, FlatVector<double> vec) NETGEN_NOEXCEPT;
  extern NGS_DLL_HEADER void SetVector (Complex val, FlatVector<Complex> vec) NETGEN_NOEXCEPT;

  extern NGS_DLL_HEADER void SetVector (double val, BareSliceVector<double> dest, size_t size) NETGEN_NOEXCEPT;
  
  INLINE void SetVector (double val, SliceVector<double> vec) NETGEN_NOEXCEPT
  {
    SetVector (val, vec, vec.Size());
  }
  
  extern NGS_DLL_HEADER void SetVector (Complex val, SliceVector<Complex> vec) NETGEN_NOEXCEPT;
  

  
  
  template <typename T1, typename T2, typename T1S, typename T2S>
  void CopyVector (LinearVector<T1,T1S> src, LinearVector<T2,T2S> dest) NETGEN_NOEXCEPT
  {
    auto cs = CombinedSize(src.Size(), dest.Size());
    for (size_t i : Range(cs))
      dest[i] = src[i];
  }

  template <typename T1, typename T2>
  void CopyVector (BareSliceVector<T1> src, BareSliceVector<T2> dest, size_t size) NETGEN_NOEXCEPT
  {
    for (size_t i : Range(size))
      dest[i] = src[i];
  }

  extern NGS_DLL_HEADER void CopyVector (BareVector<double> src, FlatVector<double> dest) NETGEN_NOEXCEPT;
  extern NGS_DLL_HEADER void CopyVector (BareSliceVector<double> src, BareSliceVector<double> dest, size_t size) NETGEN_NOEXCEPT;
  extern NGS_DLL_HEADER void CopyVector (BareVector<Complex> src, FlatVector<Complex> dest) NETGEN_NOEXCEPT;
  extern NGS_DLL_HEADER void CopyVector (BareSliceVector<Complex> src, BareSliceVector<Complex> dest, size_t size) NETGEN_NOEXCEPT;


  template <typename T0, typename T1, typename T2>
  void CopyVector (T0 alpha, BareVector<T1> src, FlatVector<T2> dest) NETGEN_NOEXCEPT
  {
    for (size_t i : Range(dest))
      dest[i] = alpha * src[i];
  }

  template <typename T0, typename T1, typename T2>
  void CopyVector (T0 alpha, BareSliceVector<T1> src, SliceVector<T2> dest) NETGEN_NOEXCEPT
  {
    for (size_t i : Range(dest))
      dest[i] = alpha * src[i];
  }

  extern NGS_DLL_HEADER void CopyVector (double alpha, BareVector<double> src, FlatVector<double> dest) NETGEN_NOEXCEPT;
  extern NGS_DLL_HEADER void CopyVector (double alpha, BareSliceVector<double> src, SliceVector<double> dest) NETGEN_NOEXCEPT;
  extern NGS_DLL_HEADER void CopyVector (Complex alpha, BareVector<Complex> src, FlatVector<Complex> dest) NETGEN_NOEXCEPT;
  extern NGS_DLL_HEADER void CopyVector (Complex alpha, BareSliceVector<Complex> src, SliceVector<Complex> dest) NETGEN_NOEXCEPT;


  
  template <typename T0, typename T1, typename T2>
  void AddVector (T0 alpha, BareVector<const T1> src, FlatVector<T2> dest) NETGEN_NOEXCEPT
  {
    for (size_t i : Range(dest))
      dest[i] += alpha*src[i];
  }

  template <typename T0, typename T1, typename T2>
  void AddVector (T0 alpha, BareSliceVector<const T1> src, SliceVector<T2> dest) NETGEN_NOEXCEPT
  {
    for (size_t i : Range(dest))
      dest[i] += alpha*src[i];
  }

  extern NGS_DLL_HEADER void AddVector (double alpha, BareVector<const double> src, FlatVector<double> dest) NETGEN_NOEXCEPT;
  extern NGS_DLL_HEADER void AddVector (double alpha, BareSliceVector<const double> src, BareSliceVector<double> dest, size_t size) NETGEN_NOEXCEPT;
  inline void AddVector (double alpha, BareSliceVector<const double> src, SliceVector<double> dest)
  {
    AddVector (alpha, src, dest, dest.Size());
  }



  // ************************ matrix and matrix-vector ops ****************
  
  template <typename TA, typename TB>
  void TransposeMatrix(SliceMatrix<TA> a, SliceMatrix<TB> b)
  {
    b = Trans(a);
  }

  void TransposeMatrix(SliceMatrix<> a, SliceMatrix<> b);

  

  typedef void (*pmult_matvec)(BareSliceMatrix<>, FlatVector<>, FlatVector<>) NETGEN_NOEXCEPT;
  extern NGS_DLL_HEADER pmult_matvec dispatch_matvec[26];
  
  inline void MultMatVec (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y) NETGEN_NOEXCEPT
  {
    size_t dsx = min(x.Size(), std::size(dispatch_matvec)-1);
    (*dispatch_matvec[dsx])  (a, x, y);    
  }


  typedef void (*pmultadd_matvec)(double s, BareSliceMatrix<>, FlatVector<>, FlatVector<>);
  extern NGS_DLL_HEADER pmultadd_matvec dispatch_addmatvec[25];
  
  inline void MultAddMatVec (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    size_t dsx = min(x.Size(), std::size(dispatch_addmatvec)-1);    
    (*dispatch_addmatvec[dsx])  (s, a, x, y);    
  }


  // typedef void (*pmult_mattransvec)(BareSliceMatrix<>, FlatVector<>, FlatVector<>);
  extern NGS_DLL_HEADER pmult_matvec dispatch_mattransvec[13];
  inline void MultMatTransVec (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    size_t dsx = min(x.Size(), std::size(dispatch_mattransvec)-1);    
    (*dispatch_mattransvec[dsx])  (a, x, y);    
  }

  // typedef void (*pmultadd_mattransvec)(double s, BareSliceMatrix<>, FlatVector<>, FlatVector<>);
  extern NGS_DLL_HEADER pmultadd_matvec dispatch_addmattransvec[13];
  inline void MultAddMatTransVec (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    size_t dsx = min(x.Size(), std::size(dispatch_addmattransvec)-1);    
    (*dispatch_addmattransvec[dsx])  (s, a, x, y);    
  }
  

  inline void MultAddMatVec (double s, BareSliceMatrix<double, ColMajor> a, FlatVector<> x, FlatVector<> y)
  {
    MultAddMatTransVec (s, Trans(a), x, y);
  }


  

  extern NGS_DLL_HEADER void MultAddMatTransVecIndirect_intern
  (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y, FlatArray<int> ind);
  typedef void (*pmultadd_mattransvecind)(double s, BareSliceMatrix<>, FlatVector<>, FlatVector<>, FlatArray<int>);
  extern NGS_DLL_HEADER pmultadd_mattransvecind dispatch_addmattransvecI[25];
  
  inline void MultAddMatTransVecIndirect (double s, BareSliceMatrix<> a,
                                          FlatVector<> x, FlatVector<> y, FlatArray<int> ind)
  {
    size_t sy = y.Size();
    if (sy <= 24)
      (*dispatch_addmattransvecI[sy])  (s, a, x, y, ind);
    else
      MultAddMatTransVecIndirect_intern (s, a, x, y, ind);
  }




  
    
  template <typename TA, typename TB, typename TC>
  INLINE void MultMatMat(SliceMatrix<TA> a, SliceMatrix<TB> b, SliceMatrix<TC> c)
  {
    c = a * b;
  }
  extern NGS_DLL_HEADER void REGCALL MultMatMat_intern (size_t ha, size_t wa, size_t wb,
                                                        BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c);

  typedef void REGCALL (*pmultABW)(size_t, size_t, size_t, BareSliceMatrix<>, BareSliceMatrix<>, BareSliceMatrix<>);
  
  extern NGS_DLL_HEADER pmultABW dispatch_multAB[14];
  inline void MultMatMat (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    if (a.Height() == 0 || b.Width() == 0) return;
    size_t wa = std::min(a.Width(), std::size(dispatch_multAB)-1);
    (*dispatch_multAB[wa])  (a.Height(), a.Width(), b.Width(), a, b, c);
  }
  
  extern NGS_DLL_HEADER pmultABW dispatch_minusmultAB[14];
  inline void MinusMultAB (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    if (a.Height() == 0 || b.Width() == 0) return;
    /*
    size_t wa = a.Width();
    if (wa >= std::size(dispatch_minusmultAB))
      wa = std::size(dispatch_minusmultAB)-1;
    */
    size_t wa = std::min(a.Width(), std::size(dispatch_minusmultAB)-1);    
    (*dispatch_minusmultAB[wa])  (a.Height(), a.Width(), b.Width(), a, b, c);
  }
  
  extern NGS_DLL_HEADER pmultABW dispatch_addAB[14];
  inline void AddAB (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    if (a.Height() == 0 || b.Width() == 0) return;
    // size_t wa = a.Width();
    // if (wa >= std::size(dispatch_addAB))
    // wa = std::size(dispatch_addAB)-1;
    size_t wa = std::min(a.Width(), std::size(dispatch_addAB)-1);    
    (*dispatch_addAB[wa])  (a.Height(), a.Width(), b.Width(), a, b, c);
  }

  extern NGS_DLL_HEADER pmultABW dispatch_subAB[14];
  inline void SubAB (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    if (a.Height() == 0 || b.Width() == 0) return;
    // size_t wa = a.Width();
    // if (wa >= std::size(dispatch_subAB))
    // wa = std::size(dispatch_subAB)-1;
    size_t wa = std::min(a.Width(), std::size(dispatch_subAB)-1);        
    (*dispatch_subAB[wa])  (a.Height(), a.Width(), b.Width(), a, b, c);
  }

  
  extern NGS_DLL_HEADER void MultMatMat_intern (size_t ha, size_t wa, size_t wb,
                                 BareSliceMatrix<> a, BareSliceMatrix<SIMD<double>> b, BareSliceMatrix<SIMD<double>> c);
  
  inline void MultMatMat (SliceMatrix<> a, SliceMatrix<SIMD<double>> b, SliceMatrix<SIMD<double>> c)
  {
    MultMatMat_intern (a.Height(), a.Width(), b.Width(), a, b, c);
  }



  

  template <bool ADD, bool POS>
  struct NGS_DLL_HEADER dispatch_atb { static pmultABW ptrs[14]; };
  
  template <bool ADD, bool POS>
  inline void MatMat_AtB (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    if (a.Height() == 0 || b.Width() == 0) return;
    /*
    size_t wa = a.Width();
    if (wa >= std::size(dispatch_atb<ADD,POS>::ptrs))
      wa = std::size(dispatch_atb<ADD,POS>::ptrs)-1;
    */
    size_t wa = std::min(a.Width(), std::size(dispatch_atb<ADD,POS>::ptrs)-1);            
    (*dispatch_atb<ADD,POS>::ptrs[wa])  (a.Height(), a.Width(), b.Width(), a, b, c);
  }
  
  inline void MultAtB (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  { MatMat_AtB<false,true> (a, b, c); }
  


  
  
  //extern NGS_DLL_HEADER void MultABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);

  typedef void REGCALL (*pfunc_abt)(size_t, size_t, BareSliceMatrix<>, BareSliceMatrix<>, BareSliceMatrix<>);
  extern NGS_DLL_HEADER pfunc_abt dispatch_abt[25];
  extern NGS_DLL_HEADER pfunc_abt dispatch_addabt[25];
  
  extern NGS_DLL_HEADER void MultABt_intern (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);
  extern NGS_DLL_HEADER void AddABt_intern (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);  

  inline void MultABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    size_t wa = a.Width();
    if (wa <= 24)
      (*dispatch_abt[wa])  (a.Height(), b.Height(), a, b, c);
    else
      MultABt_intern (a,b,c);
  }

  inline void AddABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    size_t wa = a.Width();
    if (wa <= 24)
      (*dispatch_addabt[wa])  (a.Height(), b.Height(), a, b, c);
    else
      AddABt_intern (a,b,c);
  }

  
  extern NGS_DLL_HEADER void MinusMultABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);  
  extern NGS_DLL_HEADER void SubABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);

  extern NGS_DLL_HEADER void AddABt (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<double>> b, BareSliceMatrix<double> c);  
  extern NGS_DLL_HEADER void SubABt (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<double>> b, BareSliceMatrix<double> c);


  //  copied from symbolicintegrator, needs some rework 
  extern NGS_DLL_HEADER void AddABtSym (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);    
  extern NGS_DLL_HEADER void AddABtSym (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<double>> b, BareSliceMatrix<double> c);
  
  extern NGS_DLL_HEADER void AddABt (FlatMatrix<SIMD<Complex>> a, FlatMatrix<SIMD<Complex>> b, SliceMatrix<Complex> c);
  extern NGS_DLL_HEADER void AddABtSym (FlatMatrix<SIMD<Complex>> a, FlatMatrix<SIMD<Complex>> b, SliceMatrix<Complex> c);
  extern NGS_DLL_HEADER void AddABt (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<Complex>> b, SliceMatrix<Complex> c);
  extern NGS_DLL_HEADER void AddABt (SliceMatrix<SIMD<Complex>> a, SliceMatrix<SIMD<double>> b, SliceMatrix<Complex> c);

  extern NGS_DLL_HEADER void AddABtSym (FlatMatrix<SIMD<double>> a,
                         FlatMatrix<SIMD<Complex>> b,
                         SliceMatrix<Complex> c);
  extern NGS_DLL_HEADER void AddABt (FlatMatrix<SIMD<double>> a,
                      FlatMatrix<SIMD<double>> b,
                      SliceMatrix<Complex> c);

  extern NGS_DLL_HEADER void AddABtSym (FlatMatrix<SIMD<double>> a,
                         FlatMatrix<SIMD<double>> b,
                         SliceMatrix<Complex> c);
  
  extern NGS_DLL_HEADER void AddABt (SliceMatrix<double> a,
                      SliceMatrix<double> b,
                      SliceMatrix<Complex> c);

  extern NGS_DLL_HEADER void AddABtSym (SliceMatrix<double> a,
                         SliceMatrix<double> b,
                         SliceMatrix<Complex> c);
  


  extern NGS_DLL_HEADER void AddABt (SliceMatrix<double> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c);
  extern NGS_DLL_HEADER void AddABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c);
  extern NGS_DLL_HEADER void SubABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c);  

  extern NGS_DLL_HEADER void AddABtSym (SliceMatrix<double> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c);
  extern NGS_DLL_HEADER void AddABtSym (SliceMatrix<Complex> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c);


  extern NGS_DLL_HEADER
  void ScaleCols (SliceMatrix<double,RowMajor> a, BareSliceVector<double> diag);
  extern NGS_DLL_HEADER
  void ScaleCols (SliceMatrix<double,ColMajor> a, BareSliceVector<double> diag);

  template <ORDERING ord>
  INLINE void ScaleRows (SliceMatrix<double,ord> a, BareSliceVector<double> diag)
  {
    ScaleCols (Trans(a), diag);
  }

  

  // for Cholesky and SparseCholesky
  extern NGS_DLL_HEADER
  void SubADBt (SliceMatrix<double> a,
                SliceVector<double> diag,
                SliceMatrix<double> b, SliceMatrix<double> c);

  extern NGS_DLL_HEADER
  void SubAtDB (SliceMatrix<double> a,
                SliceVector<double> diag,
                SliceMatrix<double> b, SliceMatrix<double> c);

  extern NGS_DLL_HEADER
  void SubAtDB (SliceMatrix<Complex> a,
                SliceVector<Complex> diag,
                SliceMatrix<Complex> b, SliceMatrix<Complex> c);

  template <typename T>
  void SubADBt (SliceMatrix<T,ColMajor> a,
                SliceVector<T> diag,
                SliceMatrix<T,ColMajor> b, SliceMatrix<T,ColMajor> c)
  {
    SubAtDB (Trans(b), diag, Trans(a), Trans(c));
  }  





  // MultiVector operations:

  // ip(i,j) = InnerProduct(x_i, y_j)
  extern NGS_DLL_HEADER  
  void PairwiseInnerProduct (size_t n, FlatArray<double*> x, FlatArray<double*> y, BareSliceMatrix<double> ip);

  extern NGS_DLL_HEADER
  void PairwiseInnerProduct (size_t n, FlatArray<Complex*> x, FlatArray<Complex*> y, BareSliceMatrix<Complex> ip, bool conj);


  // x_i += sum_j a(i,j) y_j
  extern NGS_DLL_HEADER  
  void MultiVectorAdd (size_t n, FlatArray<double*> x, FlatArray<double*> y, BareSliceMatrix<double> a);

  extern NGS_DLL_HEADER
  void MultiVectorAdd (size_t n, FlatArray<Complex*> x, FlatArray<Complex*> y, BareSliceMatrix<Complex> a);







  
  template <bool ADD, bool POS, ORDERING orda, ORDERING ordb>
  void NgGEMM (SliceMatrix<double,orda> a, SliceMatrix<double, ordb> b, SliceMatrix<double> c);

  template <bool ADD, bool POS, ORDERING orda, ORDERING ordb>
  void NgGEMM (SliceMatrix<double,orda> a, SliceMatrix<double, ordb> b, SliceMatrix<double,ColMajor> c);
  
  

  

  // ADD/POS 
  // f   f    C = -A*B
  // f   t    C = A*B
  // t   f    C -= A*B
  // t   t    C += A*B
  
  template <bool ADD, bool POS, ORDERING orda, ORDERING ordb>
  inline void NgGEMM (SliceMatrix<double,orda> a, SliceMatrix<double, ordb> b, SliceMatrix<double> c)
  {
    // static Timer t("generic MM, add/pos/ord="+ToString(ADD)+ToString(POS)+ToString(orda)+ToString(ordb));
    // RegionTimer r(t);
    
    // static Timer t("NgGEMM unresolved" + ToString(ADD) + ToString(POS) + ToString(orda) + ToString(ordb));
    // RegionTimer reg(t);
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Height());
    
    if (!ADD)
      {
        if (!POS)
          c = -1*a*b;
        else
          c = 1*a*b;
      }
    else
      {
        if (!POS)
          c -= 1*a*b;
        else
          c += 1*a*b;
      }
  }
  
  template <> INLINE void NgGEMM<false,true> (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  MultMatMat");
    // RegionTimer reg(t);
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Width());

    MultMatMat (a,b,c);
  }


  //  C  ???  A * B
  
  template <> INLINE void NgGEMM<true,true> (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    // nstatic Timer t("NgGEMM  AddAB");
    // RegionTimer reg(t);
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Width());

    AddAB (a,b,c);
  }

  template <> INLINE void NgGEMM<true,false> (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  SubAB");
    // RegionTimer reg(t);
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Width());

    SubAB (a,b,c);
  }

  template <> INLINE void NgGEMM<false,false> (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  MinusAB");
    // RegionTimer reg(t);
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Width());

    MinusMultAB (a, b, c);
  }

  //  C  ???  A * Bt
  
  template <> INLINE void NgGEMM<false,false> (SliceMatrix<> a, SliceMatrix<double,ColMajor> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  MinusABt");
    // RegionTimer reg(t);
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Height());

    MinusMultABt (a, Trans(b), c);
  }

  template <> INLINE void NgGEMM<false,true> (SliceMatrix<> a, SliceMatrix<double,ColMajor> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  MultABt");
    // RegionTimer reg(t);
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Height());

    MultABt (a, Trans(b), c);
  }
  
  template <> INLINE void NgGEMM<true,false> (SliceMatrix<> a, SliceMatrix<double,ColMajor> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  SubABt");
    // RegionTimer reg(t);
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Height());

    SubABt (a, Trans(b), BareSliceMatrix<>(c));
  }
  
  template <> INLINE void NgGEMM<true,true> (SliceMatrix<> a, SliceMatrix<double,ColMajor> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  AddABt");
    // RegionTimer reg(t);
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Height());

    AddABt (a, Trans(b), c);
  }
  
  
  //  C  ???  At * B
  
  template <> INLINE void NgGEMM<false,true> (SliceMatrix<double,ColMajor> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    // MultAtB (Trans(a), b, c);
    MatMat_AtB<false, true> (Trans(a), b, c);
  }
  template <> INLINE void NgGEMM<true,true> (SliceMatrix<double,ColMajor> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    MatMat_AtB<true, true> (Trans(a), b, c);
  }
  template <> INLINE void NgGEMM<true,false> (SliceMatrix<double,ColMajor> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    MatMat_AtB<true, false> (Trans(a), b, c);
  }
  template <> INLINE void NgGEMM<false,false> (SliceMatrix<double,ColMajor> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    MatMat_AtB<false, false> (Trans(a), b, c);
  }

  
  template <bool ADD, bool POS, ORDERING orda, ORDERING ordb>
  INLINE void NgGEMM (SliceMatrix<double,orda> a, SliceMatrix<double, ordb> b, SliceMatrix<double,ColMajor> c)
  {
    NgGEMM<ADD,POS> (Trans(b), Trans(a), Trans(c));
  }


  template <typename TM, typename FUNC, typename TX, typename TY>
  void NgGEMV_fallback (BareSliceMatrix<TM,RowMajor> a, FlatVector<TX> x, FlatVector<TY> y,
                       FUNC func) NETGEN_NOEXCEPT  
  {
    for (size_t i = 0; i < y.Size(); i++)
      {
        TY sum{0.0};
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), sum);
      }
  }  
  template <typename TM, typename FUNC, typename TX, typename TY>
  void NgGEMV_fallback (BareSliceMatrix<TM,ColMajor> a, FlatVector<TX> x, FlatVector<TY> y,
                        FUNC func) NETGEN_NOEXCEPT  
  {
    for (size_t i = 0; i < y.Size(); i++)
      {
        TY sum{0.0};
        for (size_t j = 0; j < x.Size(); j++)
          sum += a(i,j) * x(j);
        func(y(i), sum);
      }
  }

  template <typename TM, typename TVX, typename TVY>
  extern void TestFunc (TM m, TVX x, TVY y);

  template <bool ADD, bool POS, typename TM, ORDERING ORD, typename TX, typename TY>
  INLINE void NgGEMV (BareSliceMatrix<TM,ORD> a, FlatVector<const TX> x, FlatVector<TY> y)
  {
    if (!ADD)
      {
        if (!POS)
          // y = -1*a*x;
          NgGEMV_fallback(a, x, y, [](auto & y, auto sum) { y=-sum; });
        else
          // y = 1*a*x;
          NgGEMV_fallback(a, x, y, [](auto & y, auto sum) { y=sum; });      
      }
    else
      {
        if (!POS)
          // y -= 1*a*x;
          NgGEMV_fallback(a, x, y, [](auto & y, auto sum) { y-=sum; });
        else
          // y += 1*a*x;
          NgGEMV_fallback(a, x, y, [](auto & y, auto sum) { y+=sum; });
      }
  }

  
  
  // template <bool ADD, ORDERING ord>
  // void NgGEMV (double s, SliceMatrix<double,ord> a, BareSliceVector<double> x, BareSliceVector<double> y) NETGEN_NOEXCEPT;

  template <bool ADD, ORDERING ord>
  extern NGS_DLL_HEADER  
  void NgGEMV (double s, BareSliceMatrix<double,ord> a, FlatVector<const double> x, FlatVector<double> y) NETGEN_NOEXCEPT;

  template <bool ADD, ORDERING ord>
  extern NGS_DLL_HEADER
  void NgGEMV (const Complex s, BareSliceMatrix<Complex,ord> a, FlatVector<const Complex> x, FlatVector<Complex> y) NETGEN_NOEXCEPT;
  
  template <bool ADD, ORDERING ord>
  extern NGS_DLL_HEADER  
  void NgGEMV (Complex s, BareSliceMatrix<Complex,ord> a, FlatVector<const double> x, FlatVector<Complex> y) NETGEN_NOEXCEPT;
  template <bool ADD, ORDERING ord>
  extern NGS_DLL_HEADER  
  void NgGEMV (Complex s, BareSliceMatrix<double,ord> a, FlatVector<const Complex> x, FlatVector<Complex> y) NETGEN_NOEXCEPT;
  
  /*
  template <bool ADD, ORDERING ord>
  extern NGS_DLL_HEADER  
  void NgGEMV (double s, BareSliceMatrix<double,ord> a, SliceVector<double> x, SliceVector<double> y) NETGEN_NOEXCEPT;
  */
  
  template <bool ADD, ORDERING ord>
  extern NGS_DLL_HEADER  
  void NgGEMV (double s, BareSliceMatrix<double,ord> a, BareSliceVector<double> x, size_t sx,
               BareSliceVector<double> y, size_t sy) NETGEN_NOEXCEPT;
  
  template <bool ADD, ORDERING ord>
  INLINE void NgGEMV (double s, BareSliceMatrix<double,ord> a, SliceVector<double> x, SliceVector<double> y) NETGEN_NOEXCEPT
  {
    NgGEMV<ADD,ord> (s, a, x, x.Size(), y, y.Size());
  }

  
  template <bool ADD, ORDERING ord>
  extern NGS_DLL_HEADER  
  void NgGEMV (Complex s, BareSliceMatrix<double,ord> a, SliceVector<Complex> x, SliceVector<Complex> y) NETGEN_NOEXCEPT;

  template <bool ADD, ORDERING ord>
  extern NGS_DLL_HEADER  
  void NgGEMV (Complex s, BareSliceMatrix<Complex,ord> a, SliceVector<Complex> x, SliceVector<Complex> y) NETGEN_NOEXCEPT;


  
  template <> INLINE void NgGEMV<false,true> (BareSliceMatrix<double,RowMajor> a, FlatVector<const double> x, FlatVector<double> y)
  {
    MultMatVec (a,x.RemoveConst(),y);
  }
  
  template <> INLINE void NgGEMV<false,true> (BareSliceMatrix<double,ColMajor> a, FlatVector<const double> x, FlatVector<> y)
  {
    MultMatTransVec (Trans(a),x.RemoveConst(),y);
  }
  

  template <> INLINE void NgGEMV<true,true> (BareSliceMatrix<> a, FlatVector<const double> x, FlatVector<> y)
  {
    MultAddMatVec (1,a,x.RemoveConst(),y);
  }

  template <> INLINE void NgGEMV<true,true> (BareSliceMatrix<double,ColMajor> a, FlatVector<const double> x, FlatVector<> y)
  {
    MultAddMatTransVec (1,Trans(a),x.RemoveConst(),y);
  }
  
  template <> INLINE void NgGEMV<true,false> (BareSliceMatrix<> a, FlatVector<const double> x, FlatVector<> y)
  {
    MultAddMatVec (-1,a,x.RemoveConst(),y);
  }

  template <> INLINE void NgGEMV<true,false> (BareSliceMatrix<double,ColMajor> a, FlatVector<const double> x, FlatVector<> y)
  {
    MultAddMatTransVec (-1,Trans(a),x.RemoveConst(),y);
  }








  // bla dispatsches
  
  // vector-vector

  /*
  template <typename OP, typename T, typename TB>
  class assign_trait<OP, T, TB, 
                     enable_if_t<std::is_same_v<OP,typename MatExpr<T>::As> == true &&
                                 IsConvertibleToFlatVector<TB>()&&
                                 IsConvertibleToFlatVector<T>(), int>>
  {
  public:
    static inline T & Assign (MatExpr<T> & self, const Expr<TB> & v)
    {
      CopyVector(make_BareVector(v.Spec()), make_FlatVector(self.Spec()));
      return self.Spec();
    }
  };
  */

  /*
  template <typename OP, typename T, typename TS, typename TB, typename TBS>
  class assign_trait<OP, LinearVector<T,TS>, LinearVector<TB,TBS>, 
                     enable_if_t<std::is_same_v<OP,typename MatExpr<LinearVector<T,TS>>::As> == true, int>>
  {
  public:
    static inline auto & Assign (MatExpr<LinearVector<T,TS>> & self, const Expr<LinearVector<TB,TBS>> & v)
    {
      auto cs = CombinedSize (self.Spec().Size(), v.Spec().Size());
      CopyVector(BareVector<TB>(v.Spec()), self.Spec().Range(0, cs));
      return self.Spec();
    }
  };
  */

  /*
  template <typename OP, typename T, typename TB>
  class assign_trait<OP, T, TB, 
                     enable_if_t<std::is_same_v<OP,typename MatExpr<T>::As> == true &&
                                 ! (IsConvertibleToFlatVector<TB>() && IsConvertibleToFlatVector<T>()) &&
                                 IsConvertibleToSliceVector<TB>()&&
                                 IsConvertibleToSliceVector<T>(), int>>
  {
  public:
    static inline T & Assign (MatExpr<T> & self, const Expr<TB> & v)
    {
      CopyVector(make_BareSliceVector(v.Spec()), make_SliceVector(self.Spec()));
      return self.Spec();
    }
  };
  */

  /*
  template <typename OP, typename T, typename TS, typename TD, typename TB, typename TBS, typename TBD>
  class assign_trait<OP, VectorView<T,TS,TD>, VectorView<TB,TBS,TBD>, 
                     enable_if_t<std::is_same_v<OP,typename MatExpr<VectorView<T,TS,TD>>::As> == true, int>>
  {
    typedef VectorView<T,TS,TD> TVec;
    typedef VectorView<TB,TBS,TBD> TVecB;
  public:
    static inline auto & Assign (MatExpr<TVec> & self, const Expr<TVecB> & v)
    {
      auto cs = CombinedSize (self.Spec().Size(), v.Spec().Size());      
      CopyVector(BareSliceVector<TB>(v.Spec()), SliceVector<T>(self.Spec().Range(0,cs)));
      return self.Spec();
    }
  };
  */

  
  // x = y
  template <typename OP, typename T, typename ...Args, typename TB, typename ...BArgs>
  class assign_trait<OP, VectorView<T,Args...>, VectorView<TB, BArgs...>, 
                     enable_if_t<std::is_same_v<OP,typename MatExpr<VectorView<T,Args...>>::As> == true, int>>
  {
    typedef VectorView<T,Args...> TVec;
    typedef VectorView<TB, BArgs...> TVecB;
  public:
    static inline auto & Assign (MatExpr<TVec> & self, const Expr<TVecB> & v)
    {
      auto cs = CombinedSize (self.Spec().Size(), v.Spec().Size());

      if constexpr (is_IC<decltype(cs)>())
        {
          Vec<cs,typename remove_const<TB>::type> tmp;
          for (size_t i = 0; i<cs; i++)
            tmp[i] = v.Spec()[i];
          for (size_t i = 0; i<cs; i++)
            self.Spec()[i] = tmp[i];
          /*
            // does not allow auto-vectorization
          for (size_t i = 0; i<cs; i++)
            self.Spec()[i] = v.Spec()[i];
          */
        }
      else if constexpr (TVec::IsLinear() && TVecB::IsLinear())
        CopyVector(BareVector<TB>(v.Spec()), FlatVector<T>(self.Spec().Range(0,cs)));
      else
        CopyVector(BareSliceVector<TB>(v.Spec()), BareSliceVector<T>(self.Spec()), cs);
      return self.Spec();
    }
  };
  

  // x += s*y
  template <typename OP, typename T, typename TS, typename TD, typename TB, typename TBS, typename TBD, typename TSCAL>
  class assign_trait<OP, VectorView<T,TS,TD>, ScaleExpr<VectorView<TB,TBS,TBD>,TSCAL>,
                     enable_if_t<OP::IsAdd(), int>>
  {
    typedef VectorView<T,TS,TD> TVec;
    typedef VectorView<TB,TBS,TBD> TVecB;
  public:
    static inline auto & Assign (MatExpr<TVec> & self, const Expr<ScaleExpr<TVecB,TSCAL>> & v)
    {
      auto cs = CombinedSize (self.Spec().Size(), v.Spec().A().Size());
      auto s = v.View().S();
      if constexpr (!OP::IsPos()) s = -s;
      if constexpr (is_IC<decltype(cs)>())
        {
          Vec<cs,typename remove_const<TB>::type> tmp;
          for (size_t i = 0; i<cs; i++)
            tmp[i] = v.Spec()[i];
          for (size_t i = 0; i<cs; i++)
            self.Spec()[i] += s * tmp[i];
        }
      else
        if constexpr (TVec::IsLinear() && TVecB::IsLinear())
          AddVector(s, BareVector<const TB>(v.View().A()), FlatVector<T>(self.Spec().Range(0,cs)));
        else
          AddVector(s, BareSliceVector<const TB>(v.View().A()), SliceVector<T>(self.Spec().Range(0,cs)));
      return self.Spec();
    }
  };
  



  // matrix-vector
  // x OP= M*y
  template <typename OP, typename T, typename TS, typename TD, typename TA, typename TB, typename TBS, typename TBD>
  class assign_trait<OP, VectorView<T,TS,TD>, MultExpr<TA,VectorView<TB,TBS,TBD>>,
                     enable_if_t<IsConvertibleToBareSliceMatrix<TA>(), int>>
    
  {
    typedef VectorView<T,TS,TD> TVec;
    typedef VectorView<TB,TBS,TBD> TVecB;
  public:
    static inline auto & Assign (MatExpr<TVec> & self, const Expr<MultExpr<TA,VectorView<TB,TBS,TBD>>> & prod)
    {
      auto h = CombinedSize(get<0>(self.Spec().Shape()), get<0>(prod.View().A().Shape()));
      auto w = CombinedSize(get<0>(prod.View().B().Shape()), get<1>(prod.View().A().Shape()));

      constexpr bool ADD = OP::IsAdd();
      constexpr bool POS = OP::IsPos();
      
      if constexpr (TVec::IsLinear() && TVecB::IsLinear())
        NgGEMV<ADD,POS> (make_BareSliceMatrix(prod.View().A()).RemoveConst(),
                         FlatVector<const TB>(prod.View().B().Range(0,w)),
                         FlatVector<T>(self.Spec().Range(0,h)));
      else
        NgGEMV<ADD> (POS ? 1.0 : -1.0, make_BareSliceMatrix(prod.View().A()),
                     SliceVector<TB>(prod.View().B().Range(0,w)),
                     SliceVector<T>(self.Spec().Range(0,h)));
      return self.Spec();
    }
  };
  
  // x OP= (s*M)*y
  template <typename OP, typename T, typename TS, typename TD, typename TA, typename TB, typename TBS, typename TBD, typename TC>
  class assign_trait<OP, VectorView<T,TS,TD>, MultExpr<ScaleExpr<TA,TC>,VectorView<TB,TBS,TBD>>,
                     enable_if_t<IsConvertibleToBareSliceMatrix<TA>(), int>>
    
  {
    typedef VectorView<T,TS,TD> TVec;
    typedef VectorView<TB,TBS,TBD> TVecB;
  public:
    static inline auto & Assign (MatExpr<TVec> & self, const Expr<MultExpr<ScaleExpr<TA,TC>,VectorView<TB,TBS,TBD>>> & prod)
    {
      auto h = CombinedSize(get<0>(self.Spec().Shape()), get<0>(prod.View().A().Shape()));
      auto w = CombinedSize(get<0>(prod.View().B().Shape()), get<1>(prod.View().A().Shape()));

      constexpr bool ADD = OP::IsAdd();
      double POS = OP::IsPos() ? 1.0 : -1.0;
      
      if constexpr (TVec::IsLinear() && TVecB::IsLinear())
        NgGEMV<ADD> (POS*prod.View().A().S(), make_BareSliceMatrix(prod.View().A().A()).RemoveConst(),
                     FlatVector<const TB>(prod.View().B().Range(0,w)),
                     FlatVector<T>(self.Spec().Range(0,h)));
      else
        NgGEMV<ADD> (POS*prod.View().A().S(), make_BareSliceMatrix(prod.View().A().A()),
                     SliceVector<TB>(prod.View().B().Range(0,w)),
                     SliceVector<T>(self.Spec().Range(0,h)));
      return self.Spec();
    }
  };
  
  

#ifdef OLDMatVec
  
  template <typename OP, typename T, typename TA, typename TB>
  class assign_trait<OP, T, MultExpr<TA,TB>,
                     enable_if_t<IsConvertibleToSliceMatrix<TA,double>() &&
                                 is_convertible<TB,FlatVector<double>>::value &&
                                 is_convertible<T,FlatVector<double>>::value, int>>
  {
  public:
    static inline  T & Assign (MatExpr<T> & self, const Expr<MultExpr<TA, TB>> & prod)
    {
      auto h = CombinedSize(get<0>(self.Spec().Shape()), get<0>(prod.View().A().Shape()));
      auto w = CombinedSize(get<0>(prod.View().B().Shape()), get<1>(prod.View().A().Shape()));
      
      constexpr bool ADD = std::is_same<OP,typename MatExpr<T>::AsAdd>::value || std::is_same<OP,typename MatExpr<T>::AsSub>::value;
      constexpr bool POS = std::is_same<OP,typename MatExpr<T>::As>::value || std::is_same<OP,typename MatExpr<T>::AsAdd>::value;
      NgGEMV<ADD,POS> (make_SliceMatrix(prod.View().A()),
                       make_FlatVector(prod.View().B().Range(0,w)),
                       make_FlatVector(self.Spec().Range(0,h)));
      return self.Spec();
    }
  };

  
  template <typename OP, typename T, typename TA, typename TB>
  class assign_trait<OP, T, MultExpr<TA,TB>,
                     enable_if_t<IsConvertibleToSliceMatrix<TA,Complex>() &&
                                 IsConvertibleToFlatVector<TB>() &&
                                 IsConvertibleToFlatVector<T>(), int>>
  {
  public:    
    static inline  T & Assign (MatExpr<T> & self, const Expr<MultExpr<TA, TB>> & prod)
    {
      constexpr bool ADD = std::is_same<OP,typename MatExpr<T>::AsAdd>::value || std::is_same<OP,typename MatExpr<T>::AsSub>::value;
      constexpr double POS = (std::is_same<OP,typename MatExpr<T>::As>::value || std::is_same<OP,typename MatExpr<T>::AsAdd>::value) ? 1 : -1;
      NgGEMV<ADD> (POS, BareSliceMatrix(prod.View().A()),
                   make_FlatVector(prod.View().B()),
                   make_FlatVector(self.Spec()));
      return self.Spec();
    }
  };

  template <typename OP, typename T, typename TA, typename TB, typename TC>
  class assign_trait<OP, T, MultExpr<ScaleExpr<TA,TC>,TB>, 
                     enable_if_t<IsConvertibleToSliceMatrix<TA,Complex>() && 
                                 IsConvertibleToFlatVector<TB>() &&
                                 IsConvertibleToFlatVector<T>(), int>>
  {
  public:    
    static inline T & Assign (MatExpr<T> & self, const Expr<MultExpr<ScaleExpr<TA,TC>, TB>> & prod)
    {
      constexpr bool ADD = std::is_same<OP,typename MatExpr<T>::AsAdd>::value || std::is_same<OP,typename MatExpr<T>::AsSub>::value;
      constexpr double POS = (std::is_same<OP,typename MatExpr<T>::As>::value || std::is_same<OP,typename MatExpr<T>::AsAdd>::value) ? 1 : -1;      
      NgGEMV<ADD> (POS*prod.View().A().S(), BareSliceMatrix(prod.View().A().A()),
                   make_FlatVector(prod.View().B()),
                   make_FlatVector(self.Spec()));
      return self.Spec();
    }
  };
  
  template <typename OP, typename T, typename TA, typename TB>
  class assign_trait<OP, T, MultExpr<TA,TB>,
                     enable_if_t< ( (is_same_v<typename T::TELEM,double> ==true)||(is_same_v<typename T::TELEM,Complex> ==true) )&& 
                                  IsConvertibleToBareSliceVector<T>() &&
                                  IsConvertibleToBareSliceVector<TB>() &&
                                  (!IsConvertibleToFlatVector<TB>()||!IsConvertibleToFlatVector<T>()) && 
                                  IsConvertibleToBareSliceMatrix<TA>(),int>>
  {
  public:
    static inline T & Assign (MatExpr<T> & self, const Expr<MultExpr<TA,TB>> & prod)
    {
      auto h = CombinedSize(get<0>(self.Spec().Shape()), get<0>(prod.View().A().Shape()));
      auto w = CombinedSize(get<0>(prod.View().B().Shape()), get<1>(prod.View().A().Shape()));
      
      constexpr bool ADD = std::is_same<OP,typename MatExpr<T>::AsAdd>::value || std::is_same<OP,typename MatExpr<T>::AsSub>::value;
      constexpr double POS = (std::is_same<OP,typename MatExpr<T>::As>::value || std::is_same<OP,typename MatExpr<T>::AsAdd>::value) ? 1 : -1;
      NgGEMV<ADD> (POS, BareSliceMatrix(prod.View().A()),
                   make_BareSliceVector(prod.View().B()).Range(0, w).RemoveConst(),
                   make_BareSliceVector(self.Spec()).Range(0, h));
      return self.Spec();
    }
  };
  
  

  template <typename OP, typename T, typename TA, typename TB, typename TC>
  class assign_trait<OP, T, MultExpr<ScaleExpr<TA,TC>,TB>,
                     enable_if_t< ( (is_same_v<typename T::TELEM,double> ==true)||(is_same_v<typename T::TELEM,Complex> ==true) )&& 
                                  IsConvertibleToBareSliceVector<T>() &&
                                  IsConvertibleToBareSliceVector<TB>() &&
                                  (!IsConvertibleToFlatVector<TB>()||!IsConvertibleToFlatVector<T>()) && 
                                  IsConvertibleToBareSliceMatrix<TA>(),int>>
  {
  public:
    static inline T & Assign (MatExpr<T> & self, const Expr<MultExpr<ScaleExpr<TA,TC>, TB>> & prod)
    {
      auto h = CombinedSize(get<0>(self.Spec().Shape()), get<0>(prod.View().A().Shape()));
      auto w = CombinedSize(get<0>(prod.View().B().Shape()), get<1>(prod.View().A().Shape()));
      
      constexpr bool ADD = std::is_same<OP,typename MatExpr<T>::AsAdd>::value || std::is_same<OP,typename MatExpr<T>::AsSub>::value;
      constexpr double POS = std::is_same<OP,typename MatExpr<T>::As>::value || std::is_same<OP,typename MatExpr<T>::AsAdd>::value ? 1 : -1;
      NgGEMV<ADD> (POS*prod.View().A().S(), BareSliceMatrix(prod.View().A().A()),
                   make_BareSliceVector(prod.View().B()).Range(0, w),
                   make_BareSliceVector(self.Spec()).Range(0, h));
      return self.Spec();
    }
  };

#endif
  


  
  template <typename OP, typename T, typename TA, typename TB>
  class assign_trait<OP, T, MultExpr<TA, TB>,
                     enable_if_t<IsConvertibleToSliceMatrix<TA>() &&
                                 IsConvertibleToSliceMatrix<TB>() &&
                                 IsConvertibleToSliceMatrix<T,double>(), int>>
  {
  public:
    static inline T & Assign (MatExpr<T> & self, const Expr<MultExpr<TA, TB>> & prod) 
    {
      constexpr bool ADD = std::is_same<OP,typename MatExpr<T>::AsAdd>::value || std::is_same<OP,typename MatExpr<T>::AsSub>::value;
      constexpr bool POS = std::is_same<OP,typename MatExpr<T>::As>::value || std::is_same<OP,typename MatExpr<T>::AsAdd>::value;

      size_t n = CombinedSize(prod.View().A().Height(), self.Spec().Height());
      size_t m = CombinedSize(prod.View().B().Width(), self.Spec().Width());
      size_t k = CombinedSize(prod.View().A().Width(), prod.View().B().Height());
      
      NgGEMM<ADD,POS> (make_BareSliceMatrix(prod.View().A()).AddSize(n,k).RemoveConst(),
                       make_BareSliceMatrix(prod.View().B()).AddSize(k,m).RemoveConst(),
                       make_BareSliceMatrix(self.Spec()).AddSize(n,m));
      return self.Spec();
    }
  };

  
  template <typename OP, typename T, typename TA, typename TB>
  class assign_trait<OP, T, MultExpr<MinusExpr<TA>, TB>,
                     enable_if_t<IsConvertibleToSliceMatrix<TA,double>() &&
                                 IsConvertibleToSliceMatrix<TB,double>() &&
                                 IsConvertibleToSliceMatrix<T, double>(), int>>
  {
  public:
    static inline T & Assign (MatExpr<T> & self, const Expr<MultExpr<MinusExpr<TA>, TB>> & prod) 
    {
      constexpr bool ADD = std::is_same<OP,typename MatExpr<T>::AsAdd>::value || std::is_same<OP,typename MatExpr<T>::AsSub>::value;
      constexpr bool POS = std::is_same<OP,typename MatExpr<T>::As>::value || std::is_same<OP,typename MatExpr<T>::AsAdd>::value;
      
      NgGEMM<ADD,!POS> (make_SliceMatrix(prod.View().A().A()),
                        make_SliceMatrix(prod.View().B()),
                        make_SliceMatrix(self.Spec()));
      return self.Spec();
    }
  };

  // rank 1 update
  template <typename OP, typename T, typename TA, typename TB>
  class assign_trait<OP, T, MultExpr<TA, TransExpr<TB>>,
                     enable_if_t<IsConvertibleToSliceMatrix<T,double>() &&
                                 is_convertible<TA,FlatVector<double>>() && 
                                 is_convertible<TB,FlatVector<double>>(), int>>
  {
  public:
    static inline T & Assign (MatExpr<T> & self, const Expr<MultExpr<TA, TransExpr<TB>>> & prod)
    {
      constexpr bool ADD = std::is_same<OP,typename MatExpr<T>::AsAdd>::value || std::is_same<OP,typename MatExpr<T>::AsSub>::value;
      constexpr bool POS = std::is_same<OP,typename MatExpr<T>::As>::value || std::is_same<OP,typename MatExpr<T>::AsAdd>::value;
      
      auto veca = prod.Spec().A();
      FlatMatrix<> mata(veca.Height(), 1, veca.Data());
      auto vecb = prod.Spec().B().A();
      FlatMatrix<> matb(1, vecb.Height(), vecb.Data());
      
      NgGEMM<ADD,POS> (make_SliceMatrix(mata),
                       make_SliceMatrix(matb),
                       make_SliceMatrix(self.Spec()));
      return self.Spec();
    }
  };
  



  



  extern NGS_DLL_HEADER  
  double MatKernelMaskedScalAB (size_t n,
				double * pa, size_t da,
				double * pb, size_t db,
				const BitArray & ba);

  extern string GetTimingHelpString();
  extern list<tuple<string,double>> Timing (int what, size_t n, size_t m, size_t k,
                                            bool lapack, bool doubleprec, size_t maxits);
}


#endif


