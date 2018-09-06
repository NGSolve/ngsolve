#ifndef FILE_NGBLAS
#define FILE_NGBLAS

// optimized matrix kernels

#ifdef __clang__
#define REGCALL __regcall
#else
#define REGCALL
#endif

namespace ngbla
{

  extern NGS_DLL_HEADER void MultMatVec_intern (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y);
  typedef void (*pmult_matvec)(BareSliceMatrix<>, FlatVector<>, FlatVector<>);
  extern NGS_DLL_HEADER pmult_matvec dispatch_matvec[25];
  INLINE void MultMatVec (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    size_t sx = x.Size();
    if (sx <= 24)
      (*dispatch_matvec[sx])  (a, x, y);
    else
      MultMatVec_intern (a, x, y);
  }

  
  extern NGS_DLL_HEADER void MultMatTransVec_intern (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y);
  typedef void (*pmult_mattransvec)(BareSliceMatrix<>, FlatVector<>, FlatVector<>);
  extern NGS_DLL_HEADER pmult_mattransvec dispatch_mattransvec[13];
  
  INLINE void MultMatTransVec (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    size_t sx = x.Size();
    if (sx <= 12)
      (*dispatch_mattransvec[sx])  (a, x, y);
    else
      MultMatTransVec_intern (a, x, y);
  }


    
  template <typename TA, typename TB, typename TC>
  void MultMatMat(SliceMatrix<TA> a, SliceMatrix<TB> b, SliceMatrix<TC> c)
  {
    c = a * b;
  }
  extern NGS_DLL_HEADER void REGCALL MultMatMat_intern (size_t ha, size_t wa, size_t wb,
                                                        BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c);

  typedef void REGCALL (*pmultAB)(size_t, size_t, BareSliceMatrix<>, BareSliceMatrix<>, BareSliceMatrix<>);
  extern NGS_DLL_HEADER pmultAB dispatch_multAB[13];
  
  inline void MultMatMat (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    size_t wa = a.Width();
    if (wa <= 12)
      (*dispatch_multAB[wa])  (a.Height(), b.Width(), a, b, c);
    else
      MultMatMat_intern (a.Height(), a.Width(), b.Width(), a, b, c);    
  }


  
  extern NGS_DLL_HEADER void MultMatMat_intern (size_t ha, size_t wa, size_t wb,
                                 BareSliceMatrix<> a, BareSliceMatrix<SIMD<double>> b, BareSliceMatrix<SIMD<double>> c);
  
  inline void MultMatMat (SliceMatrix<> a, SliceMatrix<SIMD<double>> b, SliceMatrix<SIMD<double>> c)
  {
    MultMatMat_intern (a.Height(), a.Width(), b.Width(), a, b, c);
  }


  extern NGS_DLL_HEADER void REGCALL MinusMultAB_intern (size_t ha, size_t wa, size_t wb,
                                                         BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c);
  inline void MinusMultAB (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    MinusMultAB_intern (a.Height(), a.Width(), b.Width(), a, b, c);
  }

  
  extern NGS_DLL_HEADER void REGCALL AddAB_intern (size_t ha, size_t wa, size_t wb,
                                    BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c);
  inline void AddAB (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    AddAB_intern (a.Height(), a.Width(), b.Width(), a, b, c);
  }

  extern NGS_DLL_HEADER void REGCALL SubAB_intern (size_t ha, size_t wa, size_t wb,
                                                   BareSliceMatrix<> a, BareSliceMatrix<> b, BareSliceMatrix<> c);
  inline void SubAB (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    SubAB_intern (a.Height(), a.Width(), b.Width(), a, b, c);
  }



  
  
  typedef void REGCALL (*pfunc_atb)(size_t, size_t, BareSliceMatrix<>, BareSliceMatrix<>, BareSliceMatrix<>);
  extern NGS_DLL_HEADER pfunc_atb dispatch_atb[13];

  extern NGS_DLL_HEADER void MultAtB_intern (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);
    
  inline void MultAtB (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  {
    size_t wa = a.Width();
    if (wa <= 12)
      (*dispatch_atb[wa])  (a.Height(), b.Width(), a, b, c);
    else
      MultAtB_intern (a,b,c);
  }



  
  extern void MultABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);
  extern void MinusMultABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);  
  extern void AddABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);  
  extern void SubABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);

  extern void AddABt (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<double>> b, BareSliceMatrix<double> c);  
  extern void SubABt (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<double>> b, BareSliceMatrix<double> c);


  //  copied from symbolicintegrator, needs some rework 
  extern void AddABtSym (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);    
  extern void AddABtSym (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<double>> b, BareSliceMatrix<double> c);
  
  extern void AddABt (FlatMatrix<SIMD<Complex>> a, FlatMatrix<SIMD<Complex>> b, SliceMatrix<Complex> c);
  extern void AddABtSym (FlatMatrix<SIMD<Complex>> a, FlatMatrix<SIMD<Complex>> b, SliceMatrix<Complex> c);
  extern void AddABt (SliceMatrix<SIMD<double>> a, SliceMatrix<SIMD<Complex>> b, SliceMatrix<Complex> c);
  extern void AddABtSym (FlatMatrix<SIMD<double>> a,
                         FlatMatrix<SIMD<Complex>> b,
                         SliceMatrix<Complex> c);
  extern void AddABt (FlatMatrix<SIMD<double>> a,
                      FlatMatrix<SIMD<double>> b,
                      SliceMatrix<Complex> c);

  extern void AddABtSym (FlatMatrix<SIMD<double>> a,
                         FlatMatrix<SIMD<double>> b,
                         SliceMatrix<Complex> c);
  
  extern void AddABt (SliceMatrix<double> a,
                      SliceMatrix<double> b,
                      SliceMatrix<Complex> c);

  extern void AddABtSym (SliceMatrix<double> a,
                         SliceMatrix<double> b,
                         SliceMatrix<Complex> c);
  


  // for Cholesky and SparseCholesky
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


  // ADD/POS 
  // f   f    C = -A*B
  // f   t    C = A*B
  // t   f    C -= A*B
  // t   t    C += A*B

  template <bool ADD, bool POS, ORDERING orda, ORDERING ordb>
  void NgGEMM (SliceMatrix<double,orda> a, SliceMatrix<double, ordb> b, SliceMatrix<double> c)
  {
    // static Timer t("NgGEMM unresolved" + ToString(ADD) + ToString(POS) + ToString(orda) + ToString(ordb));
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
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
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Width());

    MultMatMat (a,b,c);
  }


  //  C  ???  A * B
  
  template <> INLINE void NgGEMM<true,true> (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    // nstatic Timer t("NgGEMM  AddAB");
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Width());

    AddAB (a,b,c);
  }

  template <> INLINE void NgGEMM<true,false> (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  SubAB");
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Width());

    SubAB (a,b,c);
  }

  template <> INLINE void NgGEMM<false,false> (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  MinusAB");
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Width());

    MinusMultAB (a, b, c);
  }

  //  C  ???  A * Bt
  
  template <> INLINE void NgGEMM<false,false> (SliceMatrix<> a, SliceMatrix<double,ColMajor> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  MinusABt");
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Height());

    MinusMultABt (a, Trans(b), c);
  }

  template <> INLINE void NgGEMM<false,true> (SliceMatrix<> a, SliceMatrix<double,ColMajor> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  MultABt");
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Height());

    MultABt (a, Trans(b), c);
  }
  
  template <> INLINE void NgGEMM<true,false> (SliceMatrix<> a, SliceMatrix<double,ColMajor> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  SubABt");
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Height());

    SubABt (a, Trans(b), BareSliceMatrix<>(c));
  }
  
  template <> INLINE void NgGEMM<true,true> (SliceMatrix<> a, SliceMatrix<double,ColMajor> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  AddABt");
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Height());

    AddABt (a, Trans(b), c);
  }
  
  
  //  C  ???  At * B
  
  template <> INLINE void NgGEMM<false,true> (SliceMatrix<double,ColMajor> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    // static Timer t("NgGEMM  MinusABt");
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width()*b.Height());

    MultAtB (Trans(a), b, c);
  }
  
  
  template <bool ADD, bool POS, ORDERING orda, ORDERING ordb>
  void NgGEMM (SliceMatrix<double,orda> a, SliceMatrix<double, ordb> b, SliceMatrix<double,ColMajor> c)
  {
    NgGEMM<ADD,POS> (Trans(b), Trans(a), Trans(c));
  }



  template <bool ADD, bool POS, ORDERING ord>
  void NgGEMV (SliceMatrix<double,ord> a, FlatVector<double> x, FlatVector<double> y)
  {
    // static Timer t("NgGEMV unresolved" + ToString(ADD) + ToString(POS) + ToString(ord));
    // ThreadRegionTimer reg(t, TaskManager::GetThreadId());
    // NgProfiler::AddThreadFlops (t, TaskManager::GetThreadId(), a.Height()*a.Width());
    
    if (!ADD)
      {
        if (!POS)
          y = -1*a*x;
        else
          y = 1*a*x;
      }
    else
      {
        if (!POS)
          y -= 1*a*x;
        else
          y += 1*a*x;
      }
  }

  template <> INLINE void NgGEMV<false,true> (SliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    MultMatVec (a,x,y);
  }
  
  template <> INLINE void NgGEMV<false,true> (SliceMatrix<double,ColMajor> a, FlatVector<> x, FlatVector<> y)
  {
    MultMatTransVec (Trans(a),x,y);
  }
  
  
  extern list<tuple<string,double>> Timing (int what, size_t n, size_t m, size_t k);

}


#endif
