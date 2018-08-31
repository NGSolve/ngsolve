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

  extern void REGCALL MultMatVec (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y);
  extern void REGCALL MultMatTransVec (BareSliceMatrix<> a, FlatVector<> x, FlatVector<> y);

    
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

  
  
  extern list<tuple<string,double>> Timing (int what, size_t n, size_t m, size_t k);

}


#endif
