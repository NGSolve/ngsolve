#ifndef FILE_TRIANGULAR
#define FILE_TRIANGULAR

/****************************************************************************/
/* File:   triangular.hpp                                                   */
/* Author: Joachim Schoeberl                                                */
/* Date:   Nov 2020                                                         */
/****************************************************************************/

namespace ngbla
{


  enum TRIG_SIDE { LowerLeft, UpperRight };
  enum TRIG_NORMAL { NonNormalized = 0, Normalized = 1 };
  constexpr TRIG_SIDE operator! (TRIG_SIDE s)
  { return (s == LowerLeft) ? UpperRight : LowerLeft; }


  // Solve T X = Y
  // rhs is overwritten by solution
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized,
            typename TT, typename TX, ORDERING OT, ORDERING OX>
  void TriangularSolve (BareSliceMatrix<TT,OT> T, SliceMatrix<TX,OX> X)
  {
    size_t n = X.Height();
  
    if (n == 0) return;
    if (n == 1)
      {
        if (NORM == NonNormalized)
          X.Row(0) *= 1.0/T(0.0);
        return;
      }

    if (n < 8)
      {
        if constexpr (SIDE == UpperRight)
            for (size_t i = n; i-- > 0; )
              {
                for (size_t j = i+1; j < n; j++)
                  X.Row(i) -= T(i,j) * X.Row(j);
                if (NORM==NonNormalized)
                  X.Row(i) *= 1.0/T(i,i);
              }
        else
          for (size_t i = 0; i < n; i++)
            {
              for (size_t j = 0; j < i; j++)
                X.Row(i) -= T(i,j) * X.Row(j);
              if (NORM==NonNormalized)
                X.Row(i) *= 1.0/T(i,i);
            }
        return;
      }



    
    if (X.Width() > 256)
      {
        size_t m = X.Width();
        IntRange r1(0,m/2), r2(m/2,m);      
        TriangularSolve<SIDE,NORM> (T, X.Cols(r1));
        TriangularSolve<SIDE,NORM> (T, X.Cols(r2));
        return;
      }

  
    IntRange r1(0,n/2), r2(n/2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    auto T12 = T.Rows(r1).Cols(r2).AddSize(r1.Size(), r2.Size());
    auto T21 = T.Rows(r2).Cols(r1).AddSize(r2.Size(), r1.Size());
    auto T22 = T.Rows(r2).Cols(r2);
    auto X1 = X.Rows(r1);
    auto X2 = X.Rows(r2);

    if (SIDE == LowerLeft)
      {
        TriangularSolve<SIDE,NORM> (T11.Bare(), X1);
        X2 -= T21 * X1;
        TriangularSolve<SIDE,NORM> (T22.Bare(), X2);
      }
    else
      {
        TriangularSolve<SIDE,NORM> (T22.Bare(), X2);
        X1 -= T12 * X2;
        TriangularSolve<SIDE,NORM> (T11.Bare(), X1);
      }
  }

  extern NGS_DLL_HEADER void TriangularSolveLL (BareSliceMatrix<double> T, SliceMatrix<double> X);
  template <> inline void TriangularSolve<LowerLeft,NonNormalized> (BareSliceMatrix<double> T, SliceMatrix<double> X)
  {
    TriangularSolveLL(T,X);
  }
  extern NGS_DLL_HEADER void TriangularSolveLLN (BareSliceMatrix<double> T, SliceMatrix<double> X);
  template <> inline void TriangularSolve<LowerLeft,Normalized> (BareSliceMatrix<double> T, SliceMatrix<double> X)
  {
    TriangularSolveLLN(T,X);
  }

  extern NGS_DLL_HEADER void TriangularSolveUR (BareSliceMatrix<double> T, SliceMatrix<double> X);
  template <> inline void TriangularSolve<UpperRight,NonNormalized> (BareSliceMatrix<double> T, SliceMatrix<double> X)
  {
    TriangularSolveUR(T,X);
  }
  extern NGS_DLL_HEADER void TriangularSolveURN (BareSliceMatrix<double> T, SliceMatrix<double> X);
  template <> inline void TriangularSolve<UpperRight,Normalized> (BareSliceMatrix<double> T, SliceMatrix<double> X)
  {
    TriangularSolveURN(T,X);
  }

  

  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, typename TX,
            typename enable_if<IsConvertibleToBareSliceMatrix<TT>(),int>::type = 0,
            typename enable_if<IsConvertibleToSliceMatrix<TX>(),int>::type = 0>
  void TriangularSolve (const TT & T, TX && X)
  {
    TriangularSolve<SIDE,NORM> (make_BareSliceMatrix(T), make_SliceMatrix(X));
  }

  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, 
            typename enable_if<IsConvertibleToSliceMatrix<TT>(),int>::type = 0>
  void TriangularSolve (const TT & T, FlatVector<> x)
  {
    TriangularSolve<SIDE,NORM> (make_BareSliceMatrix(T), SliceMatrix<>(x.Size(),1,1,&x(0)));
  }





  // Y = T X
  // input X is overwritten 
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized,
            typename TT, typename TX, ORDERING OT, ORDERING OX>
  void TriangularMult2 (BareSliceMatrix<TT,OT> T, SliceMatrix<TX,OX> X)
  {
    // static Timer t("TriangularMult generic, rec"); RegionTimer r(t);
    
    size_t n = X.Height();
  
    if (n == 0) return;
    if (n == 1)
      {
        if (NORM == NonNormalized)
          X.Row(0) *= T(0.0);
        return;
      }

    if (n < 8) {
      if constexpr (SIDE == UpperRight) { 
          for (size_t i = 0; i < n; i++)
            {
              if (NORM==NonNormalized)
                X.Row(i) *= T(i,i);
              for (size_t j = i+1; j < n; j++)
                X.Row(i) += T(i,j) * X.Row(j);
            }
        }
      else {
        for (size_t i = n; i > 0; i--)
          {
            size_t im = i-1;
            if (NORM==NonNormalized)
              X.Row(im) *= T(im,im);
            for (size_t j = 0; j < im; j++)
              X.Row(im) += T(im,j) * X.Row(j);
          }
      }
      return;
    }

    size_t n2 = n/2;
    // if (n2 > 4) n2 = n2 - n2%4;
    IntRange r1(0,n2), r2(n2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    auto T12 = T.Rows(r1).Cols(r2).AddSize(r1.Size(), r2.Size());;
    auto T21 = T.Rows(r2).Cols(r1).AddSize(r2.Size(), r1.Size());
    auto T22 = T.Rows(r2).Cols(r2);
    auto X1 = X.Rows(r1);
    auto X2 = X.Rows(r2);

    if (SIDE == LowerLeft)
      {
        TriangularMult2<SIDE,NORM> (T22.Bare(), X2);
        X2 += T21 * X1;
        TriangularMult2<SIDE,NORM> (T11.Bare(), X1);
      }
    else
      {
        TriangularMult2<SIDE,NORM> (T11.Bare(), X1);
        X1 += T12 * X2;
        TriangularMult2<SIDE,NORM> (T22.Bare(), X2);
      }
  }

  // Y = T X
  // input X is overwritten 
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized,
            typename TT, typename TX, ORDERING OT, ORDERING OX>
  void TriangularMult (BareSliceMatrix<TT,OT> T, SliceMatrix<TX,OX> X)
  {
    static Timer t("TriangularMult generic"); RegionTimer r(t);
    size_t i = 0;
    constexpr size_t bw = 256;
    for ( ; i+bw <= X.Width(); i += bw)
      TriangularMult2<SIDE,NORM> (T, X.Cols(i,i+bw));
    if (i < X.Width())
      TriangularMult2<SIDE,NORM> (T, X.Cols(i,X.Width()));      
  }

  
  extern NGS_DLL_HEADER void TriangularMultLL (BareSliceMatrix<double> T, SliceMatrix<double> X);
  extern NGS_DLL_HEADER void TriangularMultLLN (BareSliceMatrix<double> T, SliceMatrix<double> X);
  extern NGS_DLL_HEADER void TriangularMultLL (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X);
  extern NGS_DLL_HEADER void TriangularMultLLN (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X);

  template <> inline void TriangularMult<LowerLeft,NonNormalized> (BareSliceMatrix<double> T, SliceMatrix<double> X)
  {
    TriangularMultLL(T,X);
  }

  template <> inline void TriangularMult<LowerLeft,Normalized> (BareSliceMatrix<double> T, SliceMatrix<double> X)
  {
    TriangularMultLLN(T,X);
  }

  template <> inline void TriangularMult<LowerLeft,NonNormalized> (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X)
  {
    TriangularMultLL(T,X);
  }

  template <> inline void TriangularMult<LowerLeft,Normalized> (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X)
  {
    TriangularMultLLN(T,X);
  }


  
  extern NGS_DLL_HEADER void TriangularMultUR (BareSliceMatrix<double> T, SliceMatrix<double> X);
  extern NGS_DLL_HEADER void TriangularMultURN (BareSliceMatrix<double> T, SliceMatrix<double> X);
  extern NGS_DLL_HEADER void TriangularMultUR (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X);
  extern NGS_DLL_HEADER void TriangularMultURN (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X);
  
  template <> inline void TriangularMult<UpperRight,NonNormalized> (BareSliceMatrix<double> T, SliceMatrix<double> X)
  {
    TriangularMultUR(T,X);
  }

  template <> inline void TriangularMult<UpperRight,Normalized> (BareSliceMatrix<double> T, SliceMatrix<double> X)
  {
    TriangularMultURN(T,X);
  }

  template <> inline void TriangularMult<UpperRight,NonNormalized> (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X)
  {
    TriangularMultUR(T,X);
  }
  template <> inline void TriangularMult<UpperRight,Normalized> (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X)
  {
    TriangularMultURN(T,X);
  }

  
  
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, typename TX,
            typename enable_if<IsConvertibleToSliceMatrix<TT>(),int>::type = 0,
            typename enable_if<IsConvertibleToSliceMatrix<TX>(),int>::type = 0>
  void TriangularMult (const TT & T, TX & X)
  {
    TriangularMult<SIDE,NORM> (make_BareSliceMatrix(T), make_SliceMatrix(X));
  }
  
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, 
            typename enable_if<IsConvertibleToSliceMatrix<TT>(),int>::type = 0>
  void TriangularMult (const TT & T, FlatVector<> x)
  {
    TriangularMult<SIDE,NORM> (make_BareSliceMatrix(T), SliceMatrix<>(x.Size(),1,1,&x(0)));
  }


  // X = X * T
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, typename TX,
            typename enable_if<IsConvertibleToSliceMatrix<TT>(),int>::type = 0,
            typename enable_if<IsConvertibleToSliceMatrix<TX>(),int>::type = 0>
  void MultTriangular (TX & X, const TT & T)
  {
    TriangularMult<!SIDE,NORM> (Trans(make_BareSliceMatrix(T)), Trans(make_SliceMatrix(X)));
  }


  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized,
            typename TT, typename TX, ORDERING OT, ORDERING OX>
  void MultTriangular (SliceMatrix<TX,OX> X, SliceMatrix<TT,OT> T)
  {
    TriangularMult<!SIDE,NORM> (Trans(T), Trans(X));
  }

  extern NGS_DLL_HEADER void MultTriangularLLN (SliceMatrix<double> X, BareSliceMatrix<double> T);
  template <> inline void MultTriangular<LowerLeft,Normalized> (SliceMatrix<double> X, SliceMatrix<double> T)
  {
    MultTriangularLLN(X,T);
  }
  extern NGS_DLL_HEADER void MultTriangularUR (SliceMatrix<double> X, BareSliceMatrix<double> T);
  template <> inline void MultTriangular<UpperRight> (SliceMatrix<double> X, SliceMatrix<double> T)
  {
    MultTriangularUR(X,T);
  }



  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, ORDERING OT, ORDERING OXY>
  extern NGS_DLL_HEADER void GeneralizedTriangularMult_SM (SliceMatrix<double, OT> T,
                                                           SliceMatrix<double, OXY> X,
                                                           SliceMatrix<double, OXY> Y);
    
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, typename TX, typename TY>
  void GeneralizedTriangularMult (const TT & T,
                                  const TX & X,
                                  const TY & Y)
  {
    GeneralizedTriangularMult_SM<SIDE,NORM> (make_SliceMatrix(T), make_SliceMatrix(X), make_SliceMatrix(Y));
  }
  


  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, ORDERING OT, ORDERING OXY>
  extern NGS_DLL_HEADER void GeneralizedTriangularSub_SM (SliceMatrix<double, OT> T,
                                                           SliceMatrix<double, OXY> X,
                                                           SliceMatrix<double, OXY> Y);

    
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, typename TX, typename TY>
  void GeneralizedTriangularSub (const TT & T,
                                 const TX & X,
                                 const TY & Y)
  {
    GeneralizedTriangularSub_SM<SIDE,NORM> (make_SliceMatrix(T), make_SliceMatrix(X), make_SliceMatrix(Y));
  }
  

  


  

  /*
  static Timer triginv_loops ("TriangularInvert loops");
  static Timer triginv_L1 ("TriangularInvert L1");
  static Timer triginv_L2 ("TriangularInvert L2");
  static Timer triginv_R1 ("TriangularInvert R1");
  static Timer triginv_R2 ("TriangularInvert r2");
  */
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized,
            typename TT, ORDERING TO>
  void TriangularInvert (SliceMatrix<TT,TO> T)
  {
    size_t n = T.Height();
  
    if (n == 0) return;
    if (n == 1)
      {
        if (NORM == NonNormalized)
          T(0,0) = 1.0/T(0,0);
        return;
      }

    if (n < 16)
      {
        // RegionTimer reg(triginv_loops);
        for (size_t j = 0; j < n; j++)
          {
            auto rowj = T.Row(j);
            TT invdiag = 1.0;
            if constexpr (NORM == NonNormalized) {
                invdiag = 1.0/T(j,j);  
                
                if constexpr (SIDE==LowerLeft) {
                    rowj.Range(0,j) *= invdiag;
                  }
                else
                  rowj.Range(j+1,n) *= invdiag;
                rowj(j) = invdiag;
              }
            
            if constexpr (SIDE == UpperRight) {
                for (size_t k = 0; k < j; k++)
                  {
                    auto rowk = T.Row(k);
                    TT help = rowk(j);
                    rowk.Range(j+1,n) -= help * rowj.Range(j+1,n);
                    rowk(j) = -help*invdiag;
                  }
              }
            else
              for (size_t k = j+1; k < n; k++)
                {
                  auto rowk = T.Row(k);                  
                  TT help = rowk(j);
                  rowk.Range(j) -= help * rowj.Range(j);
                  rowk(j) = -help*invdiag;
                }
          }
        
        return;
      }

    size_t n2 = n/2;
    IntRange r1(0,n2), r2(n2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    auto T12 = T.Rows(r1).Cols(r2);
    auto T21 = T.Rows(r2).Cols(r1);
    auto T22 = T.Rows(r2).Cols(r2);

    TriangularInvert<SIDE,NORM> (T11);
    TriangularInvert<SIDE,NORM> (T22);
    
    if (SIDE == LowerLeft)
      {
        T21 *= -1;
        // triginv_L1.Start();
        TriangularMult<SIDE,NORM> (T22, T21);
        // triginv_L1.Stop();
        // triginv_L2.Start();        
        // TriangularMult<!SIDE,NORM> (Trans(T11), Trans(T21));
        MultTriangular<SIDE,NORM> (T21, T11);
        // triginv_L2.Stop();        
      }
    else
      {
        T12 *= -1;
        // triginv_R1.Start();        
        TriangularMult<SIDE,NORM> (T11, T12);
        // triginv_R1.Stop();        
        // triginv_R2.Start();        
        // TriangularMult<!SIDE,NORM> (Trans(T22), Trans(T12));
        MultTriangular<SIDE,NORM> (T12, T22);
        // triginv_R2.Stop();                
      }
  }

  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, 
            typename enable_if<IsConvertibleToSliceMatrix<TT>(),int>::type = 0>
  void TriangularInvert (const TT & T)
  {
    TriangularInvert<SIDE,NORM> (make_SliceMatrix(T));
  }




  
  extern NGS_DLL_HEADER void CalcLU (SliceMatrix<double> A, FlatArray<int> p);
  extern NGS_DLL_HEADER void InverseFromLU (SliceMatrix<double> A, FlatArray<int> p);
  extern NGS_DLL_HEADER void SolveFromLU (SliceMatrix<double> A, FlatArray<int> p, SliceMatrix<double,ColMajor> X);
  extern NGS_DLL_HEADER void SolveTransFromLU (SliceMatrix<double> A, FlatArray<int> p, SliceMatrix<double,ColMajor> X);
  

  
  
}
#endif
