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
        TriangularSolve<SIDE,NORM> (T11, X1);
        X2 -= T21 * X1;
        TriangularSolve<SIDE,NORM> (T22, X2);
      }
    else
      {
        TriangularSolve<SIDE,NORM> (T22, X2);
        X1 -= T12 * X2;
        TriangularSolve<SIDE,NORM> (T11, X1);
      }
  }



  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, typename TX,
            typename enable_if<IsConvertibleToSliceMatrix<TT>(),int>::type = 0,
            typename enable_if<IsConvertibleToSliceMatrix<TX>(),int>::type = 0>
  void TriangularSolve (const TT & T, TX & X)
  {
    TriangularSolve<SIDE,NORM> (BareSliceMatrix(make_SliceMatrix(T)), make_SliceMatrix(X));
  }

  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, 
            typename enable_if<IsConvertibleToSliceMatrix<TT>(),int>::type = 0>
  void TriangularSolve (const TT & T, FlatVector<> x)
  {
    TriangularSolve<SIDE,NORM> (BareSliceMatrix(make_SliceMatrix(T)), SliceMatrix<>(x.Size(),1,1,&x(0)));
  }





  // Y = T X
  // input X is overwritten 
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized,
            typename TT, typename TX, ORDERING OT, ORDERING OX>
  void TriangularMult (SliceMatrix<TT,OT> T, SliceMatrix<TX,OX> X)
  {
    size_t n = T.Height();
  
    if (n == 0) return;
    if (n == 1)
      {
        if (NORM == NonNormalized)
          X.Row(0) *= T(0.0);
        return;
      }

    if (X.Width() > 256)
      {
        size_t m = X.Width();
        IntRange r1(0,m/2), r2(m/2,m);      
        TriangularMult<SIDE,NORM> (T, X.Cols(r1));
        TriangularMult<SIDE,NORM> (T, X.Cols(r2));
        return;
      }

  
    IntRange r1(0,n/2), r2(n/2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    auto T12 = T.Rows(r1).Cols(r2);
    auto T21 = T.Rows(r2).Cols(r1);
    auto T22 = T.Rows(r2).Cols(r2);
    auto X1 = X.Rows(r1);
    auto X2 = X.Rows(r2);

    if (SIDE == LowerLeft)
      {
        TriangularMult<SIDE,NORM> (T22, X2);
        X2 += T21 * X1;
        TriangularMult<SIDE,NORM> (T11, X1);
      }
    else
      {
        TriangularMult<SIDE,NORM> (T11, X1);
        X1 += T12 * X2;
        TriangularMult<SIDE,NORM> (T22, X2);
      }
  }

  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, 
            typename enable_if<IsConvertibleToSliceMatrix<TT>(),int>::type = 0>
  void TriangularMult (const TT & T, FlatVector<> x)
  {
    TriangularMult<SIDE,NORM> (make_SliceMatrix(T), SliceMatrix<>(x.Size(),1,1,&x(0)));
  }











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


    IntRange r1(0,n/2), r2(n/2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    auto T12 = T.Rows(r1).Cols(r2);
    auto T21 = T.Rows(r2).Cols(r1);
    auto T22 = T.Rows(r2).Cols(r2);

    if (SIDE == LowerLeft)
      {
        TriangularInvert<SIDE,NORM> (T22);
        T21 *= -1;
        TriangularMult<SIDE,NORM> (T22, T21);
        TriangularSolve<!SIDE,NORM> (Trans(T11), Trans(T21));
        TriangularInvert<SIDE,NORM> (T11);
      }
    else
      {
        TriangularInvert<SIDE,NORM> (T11);
        T12 *= -1;
        TriangularMult<SIDE,NORM> (T11, T12);
        TriangularSolve<!SIDE,NORM> (Trans(T22), Trans(T12));
        TriangularInvert<SIDE,NORM> (T22);
      }
  }

  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, typename TT, 
            typename enable_if<IsConvertibleToSliceMatrix<TT>(),int>::type = 0>
  void TriangularInvert (const TT & T)
  {
    TriangularInvert<SIDE,NORM> (make_SliceMatrix(T));
  }




  
}
#endif
