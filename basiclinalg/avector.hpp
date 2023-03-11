#ifndef FILE_AVECTOR
#define FILE_AVECTOR



namespace ngbla
{



#if defined(__AVX__) && !defined(__AVX512F__)

  extern void AddABt (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);
  extern void AddABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);

  //extern void AddABtSym (AFlatMatrix<double> a, AFlatMatrix<double> b, BareSliceMatrix<double> c);
  // extern void AddABtSym (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c);
  extern void AddABtSym (SliceMatrix<double> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);
  extern void AddABtSym (SliceMatrix<Complex> a, SliceMatrix<Complex> b, SliceMatrix<Complex> c);

  /*
  extern void SubAtB (SliceMatrix<double> a, SliceMatrix<double> b, SliceMatrix<double> c);
  inline void SubABt (SliceMatrix<double,ColMajor> a, SliceMatrix<double,ColMajor> b, SliceMatrix<double,ColMajor> c)
  {
    SubAtB (Trans(b), Trans(a), Trans(c));
  }
  */
  
  // extern void MultMatDiagMat(AFlatMatrixD a, AFlatVectorD diag, AFlatMatrixD c);


#else // __AVX__

  // INLINE void AddABt (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  // { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }
  
  // INLINE void AddABtSym (SliceMatrix<double> a, SliceMatrix<double> b, BareSliceMatrix<double> c)
  // { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }

  INLINE void AddABt (SliceMatrix<double> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c)
  { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }
  
  INLINE void AddABtSym (SliceMatrix<double> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c)
  { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }

  INLINE void AddABt (SliceMatrix<Complex> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c)
  { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }
  
  INLINE void AddABtSym (SliceMatrix<Complex> a, SliceMatrix<Complex> b, BareSliceMatrix<Complex> c)
  { c.AddSize(a.Height(), b.Height()) += a * Trans(b) | Lapack; }

#endif // __AVX__


  /*
  template <typename TA, typename TB, typename TC, ORDERING ORD>
  INLINE void SubABt (const TA & a, const TB & b, SliceMatrix<TC,ORD> c)
  {
    c -= a * Trans(b) | Lapack;
    // LapackMultAdd (a, Trans(b), 1.0, c, 1.0);
  }
  */
  
}
#endif // FILE_AVECTOR
