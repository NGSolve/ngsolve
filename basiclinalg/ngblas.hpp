#ifndef FILE_NGBLAS
#define FILE_NGBLAS

// optimized matrix kernels

#if defined(__clang__ ) && defined(NETGEN_ARCH_AMD64)
#define REGCALL __regcall
#else
#define REGCALL
#endif

namespace ngbla
{

  extern NGS_DLL_HEADER void MultMatVec_intern (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y);
  typedef void (*pmult_matvec)(BareSliceMatrix<>, FlatVector<>, FlatVector<>);
  extern NGS_DLL_HEADER pmult_matvec dispatch_matvec[26];
  INLINE void MultMatVec (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    size_t dsx = x.Size();
    if (dsx >= std::size(dispatch_matvec))
      dsx = std::size(dispatch_matvec)-1;
    (*dispatch_matvec[dsx])  (a, x, y);    
  }


  extern NGS_DLL_HEADER void MultAddMatVec_intern (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y);
  typedef void (*pmultadd_matvec)(double s, BareSliceMatrix<>, FlatVector<>, FlatVector<>);
  extern NGS_DLL_HEADER pmultadd_matvec dispatch_addmatvec[25];
  INLINE void MultAddMatVec (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    size_t sx = x.Size();
    if (sx <= 24)
      (*dispatch_addmatvec[sx])  (s, a, x, y);
    else
      MultAddMatVec_intern (s, a, x, y);
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

  extern NGS_DLL_HEADER void MultAddMatTransVec_intern (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y);
  typedef void (*pmultadd_mattransvec)(double s, BareSliceMatrix<>, FlatVector<>, FlatVector<>);
  extern NGS_DLL_HEADER pmultadd_mattransvec dispatch_addmattransvec[13];
  
  INLINE void MultAddMatTransVec (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    size_t sx = x.Size();
    if (sx <= 12)
      (*dispatch_addmattransvec[sx])  (s, a, x, y);
    else
      MultAddMatTransVec_intern (s, a, x, y);
  }




  

  extern NGS_DLL_HEADER void MultAddMatTransVecIndirect_intern
  (double s, BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y, FlatArray<int> ind);
  typedef void (*pmultadd_mattransvecind)(double s, BareSliceMatrix<>, FlatVector<>, FlatVector<>, FlatArray<int>);
  extern NGS_DLL_HEADER pmultadd_mattransvecind dispatch_addmattransvecI[25];
  
  INLINE void MultAddMatTransVecIndirect (double s, BareSliceMatrix<> a,
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

  typedef void REGCALL (*pmultAB)(size_t, size_t, BareSliceMatrix<>, BareSliceMatrix<>, BareSliceMatrix<>);
  typedef void REGCALL (*pmultABW)(size_t, size_t, size_t, BareSliceMatrix<>, BareSliceMatrix<>, BareSliceMatrix<>);
  
  extern NGS_DLL_HEADER pmultABW dispatch_multAB[14];
  inline void MultMatMat (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    if (a.Height() == 0 || b.Width() == 0) return;
    size_t wa = a.Width();
    if (wa >= std::size(dispatch_multAB))
      wa = std::size(dispatch_multAB)-1;
    (*dispatch_multAB[wa])  (a.Height(), a.Width(), b.Width(), a, b, c);
  }
  
  extern NGS_DLL_HEADER pmultABW dispatch_minusmultAB[14];
  inline void MinusMultAB (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    if (a.Height() == 0 || b.Width() == 0) return;
    size_t wa = a.Width();
    if (wa >= std::size(dispatch_minusmultAB))
      wa = std::size(dispatch_minusmultAB)-1;
    (*dispatch_minusmultAB[wa])  (a.Height(), a.Width(), b.Width(), a, b, c);
  }
  
  extern NGS_DLL_HEADER pmultABW dispatch_addAB[14];
  inline void AddAB (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    if (a.Height() == 0 || b.Width() == 0) return;
    size_t wa = a.Width();
    if (wa >= std::size(dispatch_addAB))
      wa = std::size(dispatch_addAB)-1;
    (*dispatch_addAB[wa])  (a.Height(), a.Width(), b.Width(), a, b, c);
  }

  extern NGS_DLL_HEADER pmultABW dispatch_subAB[14];
  inline void SubAB (SliceMatrix<> a, SliceMatrix<> b, SliceMatrix<> c)
  {
    if (a.Height() == 0 || b.Width() == 0) return;
    size_t wa = a.Width();
    if (wa >= std::size(dispatch_subAB))
      wa = std::size(dispatch_subAB)-1;
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
    size_t wa = a.Width();
    if (wa >= std::size(dispatch_atb<ADD,POS>::ptrs))
      wa = std::size(dispatch_atb<ADD,POS>::ptrs)-1;
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





  

  // ADD/POS 
  // f   f    C = -A*B
  // f   t    C = A*B
  // t   f    C -= A*B
  // t   t    C += A*B

  template <bool A, bool P, ORDERING oa, ORDERING ob>
  class trait__
  {
  public:
    typedef double TELEM;
  };
  
  
  template <bool ADD, bool POS, ORDERING orda, ORDERING ordb>
  INLINE void NgGEMM (SliceMatrix<double,orda> a, SliceMatrix<double, ordb> b, SliceMatrix<double> c)
  {
    // static Timer t("generic MM, add/pos/ord="+ToString(ADD)+ToString(POS)+ToString(orda)+ToString(ordb));
    // RegionTimer r(t);
    // typename trait__<ADD,POS,orda,ordb>::TELEM x;  // to get a warning
    
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


  template <bool A, bool P, ORDERING oa>
  class vtrait__
  {
  public:
    typedef double TELEM;
  };

  template <bool ADD, bool POS, ORDERING ord>
  void NgGEMV (SliceMatrix<double,ord> a, FlatVector<double> x, FlatVector<double> y)
  {
    typename vtrait__<ADD,POS,ord>::TELEM hv;  // to get a warning
    
    // static Timer t("generic MV, add/pos/ord="+ToString(ADD)+ToString(POS)+ToString(ord));
    // RegionTimer r(t);
    // cout << "generic nggemv , add = " << ADD << ", pos = " << POS << endl;
    // static Timer t("NgGEMV unresolved" + ToString(ADD) + ToString(POS) + ToString(ord));
    // RegionTimer reg(t);
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
  

  template <> INLINE void NgGEMV<true,true> (SliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    MultAddMatVec (1,a,x,y);
  }

  template <> INLINE void NgGEMV<true,true> (SliceMatrix<double,ColMajor> a, FlatVector<> x, FlatVector<> y)
  {
    MultAddMatTransVec (1,Trans(a),x,y);
  }
  
  template <> INLINE void NgGEMV<true,false> (SliceMatrix<> a, FlatVector<> x, FlatVector<> y)
  {
    MultAddMatVec (-1,a,x,y);
  }

  template <> INLINE void NgGEMV<true,false> (SliceMatrix<double,ColMajor> a, FlatVector<> x, FlatVector<> y)
  {
    MultAddMatTransVec (-1,Trans(a),x,y);
  }

  
  extern list<tuple<string,double>> Timing (int what, size_t n, size_t m, size_t k, bool lapack, size_t maxits);


  double MatKernelMaskedScalAB (size_t n,
				double * pa, size_t da,
				double * pb, size_t db,
				const BitArray & ba);
}


#endif
