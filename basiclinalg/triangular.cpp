#include <bla.hpp>


namespace ngbla
{
#include "matkernel.hpp"


  // these kernels should go into generate_mat_kernels ...
  // T * X      T .. triangular, X rowwise
  // L/R, Norm/NonNorm,  Solve + Mult
  //
  // not so nice: T * X  ... X colwise    needed ???
  

  template <size_t H, OPERATION OP, ORDERING ORD> 
  INLINE void MatKernel2AddAB (size_t hb, size_t wb, double * pa, size_t da, double * pb, size_t db, double * pc, size_t dc)
  {
    if constexpr (H != 0)
      {
    // static Timer t("matkernel2addab"+ToString(H)+"x"+ToString(hb)+"x"+ToString(wb));
    // t.AddFlops(H*hb*wb);
    
    // cout << "addab, ORD = " << ToString(ORD) << ", hb = " << hb << ", wb = " << wb << ", distb = " << db << endl;
    if constexpr (ORD == RowMajor)
      {
        // RegionTimer r(t);
        constexpr size_t SW = SIMD<double>::Size();
        size_t l = 0, lb = 0;
        for ( ; l+3*SW <= wb; l += 3*SW, lb += 3*SW)
          MatKernelMultAB<H,3,OP> (hb, pa, da, pb+lb, db, pc+l, dc);
        for ( ; l+SW <= wb; l += SW, lb += SW)
          MatKernelMultAB<H,1,OP> (hb, pa, da, pb+lb, db, pc+l, dc);
        if (l < wb)
          MatKernelMultABMask<H,OP>(hb, SIMD<mask64>(wb-l), pa, da, pb+lb, db, pc+l, dc);
      }
    else
      {
        // static Timer tcopy("copya");
        // tcopy.Start();
        alignas (64) double mema[96][H];
        if (hb <= 96)
          {
            for (int i = 0; i < hb; i++)
              for (int j = 0; j < H; j++)
                mema[i][j] = pa[i*da+j];
            pa = &mema[0][0];
            da = H;
          }
        // tcopy.Stop();
        // RegionTimer r(t);        
        MatKernelAtB_SmallWA2<H,OP>(hb, wb, pa, da, pb, db, pc, dc);
      }
      }
  }

  
  template <int SW>
  inline void CopyMatrixIn (size_t h, size_t w,
                            double * ps, size_t dists,
                            SIMD<double,SW> * pd, size_t distd)
  {
    if (w % SW == 0)
      {
        for (size_t i = 0; i < h; i++, pd += distd, ps += dists)
          {
            auto ps2 = ps;
            auto pd2 = pd;
            size_t js = 0; 
            for ( ; js+SW <= w; js+=SW, pd2++, ps2+=SW)
              *pd2 = SIMD<double>(ps2);
          }
      }
    else
      {
        SIMD<mask64> mask(w % SW);
        for (size_t i = 0; i < h; i++, pd += distd, ps += dists)
          {
            auto ps2 = ps;
            auto pd2 = pd;
            
            size_t js = 0; 
            for ( ; js+SW <= w; js+=SW, pd2++, ps2+=SW)
              *pd2 = SIMD<double>(ps2);
            SIMD<double>(ps2, mask).Store((double*) (pd2), mask);
          }
      }
  }




  /* ***************************  TriangularMult - Normalized **************************** */


  template <TRIG_NORMAL NORM, ORDERING ORD>
  void TriangularMultLL3 (BareSliceMatrix<double,ORD> bL, SliceMatrix<double> X)
  {
    size_t n = X.Height();
    size_t i = n;
    auto L = bL.AddSize(n,n);
    
    constexpr size_t HA = 4;
    
    size_t remainder = n % HA;
    if (remainder > 0)
      Switch<HA> (remainder, [L,X,i] (auto r)
                  {
                    KernelTriangularMult<LowerLeft,NORM,ORD,r.value> (X.Width(), &L(i-r.value,i-r.value), L.Dist(), &X(i-r.value,0), X.Dist());
                    if (i > r)
                      MatKernel2AddAB<r.value,ADD,ORD> (i-r.value, X.Width(),
                                                        &L(i-r.value,0), L.Dist(),
                                                        &X(0,0), X.Dist(),
                                                        &X(i-r.value,0), X.Dist());
                  });
    
    i -= remainder;
    for ( ; i >= HA; i -= HA)
      {
        KernelTriangularMult<LowerLeft,NORM,ORD,4> (X.Width(), &L(i-HA,i-HA), L.Dist(), &X(i-HA,0), X.Dist());
        if (i > HA)
          MatKernel2AddAB<HA,ADD,ORD> (i-HA, X.Width(), &L(i-HA,0), L.Dist(), &X(0,0), X.Dist(), &X(i-HA,0), X.Dist());
      }
  }

  template <TRIG_NORMAL NORM, ORDERING ORD>
  void TriangularMultLL2 (BareSliceMatrix<double,ORD> L, SliceMatrix<double> X)  
  {
    size_t i = 0;
    constexpr size_t bw = 256;
    for ( ; i+bw <= X.Width(); i += bw)
      TriangularMultLL3<NORM> (L, X.Cols(i,i+bw));
    if (i < X.Width())
      TriangularMultLL3<NORM> (L, X.Cols(i,X.Width()));      
  }

  template <TRIG_NORMAL NORM, ORDERING ORD>  
  void TriangularMultLL1 (BareSliceMatrix<double,ORD> T, SliceMatrix<double> X)
  {
    size_t n = X.Height();
    if (n < 128)
      {
        TriangularMultLL2<NORM> (T,X);
        return;
      }
    size_t n2 = n/2;
    n2 -= n2 % (3*SIMD<double>::Size());
    // IntRange r1(0,n/2), r2(n/2,n);
    IntRange r1(0,n2), r2(n2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    auto T21 = T.Rows(r2).Cols(r1).AddSize(r2.Size(), r1.Size());
    auto T22 = T.Rows(r2).Cols(r2);
    auto X1 = X.Rows(r1);
    auto X2 = X.Rows(r2);

    TriangularMultLL1<NORM> (T22, X2);
    X2 += T21 * X1;
    TriangularMultLL1<NORM> (T11, X1);
  }

  void TriangularMultLL (BareSliceMatrix<double> T, SliceMatrix<double> X)
  { TriangularMultLL1<NonNormalized> (T,X); }
  void TriangularMultLLN (BareSliceMatrix<double> T, SliceMatrix<double> X)
  { TriangularMultLL1<Normalized> (T,X); }
  void TriangularMultLL (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X)
  { TriangularMultLL1<NonNormalized> (T,X); }
  void TriangularMultLLN (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X)
  { TriangularMultLL1<Normalized> (T,X); }

  
  
  template <TRIG_NORMAL NORM, ORDERING ORD>
  void TriangularMultUR3 (BareSliceMatrix<double,ORD> bU, SliceMatrix<double> X)  
  {
    size_t n = X.Height();
    auto U = bU.AddSize(n,n);

    /*
    for (size_t i = 0; i < n; i++)
      {
        X.Row(i) *= U(i,i);
        for (size_t j = i+1; j < n; j++)
          X.Row(i) += U(i,j) * X.Row(j);
      }
    return;
    */

    constexpr size_t HA = 4;

    size_t i = 0;
    for ( ; i+HA <= n; i += HA)
      {
        KernelTriangularMult<UpperRight,NORM,ORD,4> (X.Width(), &U(i,i), U.Dist(), &X(i,0), X.Dist());
        if (i+HA < n)
          MatKernel2AddAB<HA,ADD,ORD> (n-i-HA, X.Width(), &U(i,i+HA), U.Dist(), &X(i+HA,0), X.Dist(), &X(i,0), X.Dist());
      }

    size_t remainder = n % HA;
    if (remainder > 0)
      Switch<HA> (remainder, [U,X,i] (auto r)
                  {
                    KernelTriangularMult<UpperRight,NORM,ORD,r.value> (X.Width(), &U(i,i), U.Dist(), &X(i,0), X.Dist());
                  });
  }





  
  template <TRIG_NORMAL NORM, ORDERING ORD>
  void TriangularMultUR2 (BareSliceMatrix<double,ORD> L, SliceMatrix<double> X)  
  {
    size_t i = 0;
    constexpr size_t bw = 192;
    for ( ; i+bw <= X.Width(); i += bw)
      TriangularMultUR3<NORM> (L, X.Cols(i,i+bw));
    if (i < X.Width())
      TriangularMultUR3<NORM> (L, X.Cols(i,X.Width()));      
  }

  template <TRIG_NORMAL NORM, ORDERING ORD>  
  void TriangularMultUR1 (BareSliceMatrix<double,ORD> T, SliceMatrix<double> X)
  {
    size_t n = X.Height();
    if (n < 128)
      {
        TriangularMultUR2<NORM> (T,X);
        return;
      }

    size_t n2 = n/2;
    n2 -= n2 % (3*SIMD<double>::Size());    
    IntRange r1(0,n2), r2(n2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    auto T12 = T.Rows(r1).Cols(r2).AddSize(r1.Size(), r2.Size());
    auto T22 = T.Rows(r2).Cols(r2);
    auto X1 = X.Rows(r1);
    auto X2 = X.Rows(r2);

    TriangularMultUR1<NORM> (T11, X1);
    X1 += T12 * X2;
    TriangularMultUR1<NORM> (T22, X2);
  }
  
  void TriangularMultUR (BareSliceMatrix<double> T, SliceMatrix<double> X)
  { TriangularMultUR1<NonNormalized> (T,X); }
  void TriangularMultURN (BareSliceMatrix<double> T, SliceMatrix<double> X)
  { TriangularMultUR1<Normalized> (T,X); }
  void TriangularMultUR (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X)
  { TriangularMultUR1<NonNormalized> (T,X); }
  void TriangularMultURN (BareSliceMatrix<double,ColMajor> T, SliceMatrix<double> X)
  { TriangularMultUR1<Normalized> (T,X); }




  /* ***************************  TriangularSolve **************************** */


  template <TRIG_NORMAL NORM>
  void TriangularSolveLL3 (BareSliceMatrix<double> bL, SliceMatrix<double> X)
  {
    size_t n = X.Height();
    auto L = bL.AddSize(n,n);
    
    constexpr size_t HA = 4;

    size_t i = 0;
    for ( ; i+HA <= n; i += HA)
      {
        if (i > 0)
          MatKernel2AddAB<HA,SUB,RowMajor> (i, X.Width(), &L(i,0), L.Dist(), &X(0,0), X.Dist(), &X(i,0), X.Dist());
        KernelTriangularSolve<LowerLeft,NORM,RowMajor,4> (X.Width(), &L(i,i), L.Dist(), &X(i,0), X.Dist());
      }
    
    size_t remainder = n % HA;
    if (remainder > 0)
      {
        Switch<HA> (remainder, [L,X,i] (auto r)
                    {
                      if constexpr (r.value > 0) {
                          if (i > 0)
                            MatKernel2AddAB<r.value,SUB,RowMajor> (i, X.Width(),
                                                                   &L(i,0), L.Dist(),
                                                                   &X(0,0), X.Dist(),
                                                                   &X(i,0), X.Dist());
                          
                          KernelTriangularSolve<LowerLeft,NORM,RowMajor,r.value> (X.Width(), &L(i,i), L.Dist(), &X(i,0), X.Dist());
                        }
                    });
      }
  }

  template <TRIG_NORMAL NORM>  
  void TriangularSolveLL2 (BareSliceMatrix<double> L, SliceMatrix<double> X)  
  {
    size_t i = 0;
    constexpr size_t bw = 128;
    for ( ; i+bw <= X.Width(); i += bw)
      TriangularSolveLL3<NORM> (L, X.Cols(i,i+bw));
    if (i < X.Width())
      TriangularSolveLL3<NORM> (L, X.Cols(i,X.Width()));      
  }
  
  template <TRIG_NORMAL NORM>
  void TriangularSolveLL1 (BareSliceMatrix<double> T, SliceMatrix<double> X)
  {
    size_t n = X.Height();
    if (n < 128)
      {
        TriangularSolveLL2<NORM> (T,X);
        return;
      }
    IntRange r1(0,n/2), r2(n/2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    // auto T12 = T.Rows(r1).Cols(r2).AddSize(r1.Size(), r2.Size());
    auto T21 = T.Rows(r2).Cols(r1).AddSize(r2.Size(), r1.Size());
    auto T22 = T.Rows(r2).Cols(r2);
    auto X1 = X.Rows(r1);
    auto X2 = X.Rows(r2);

    TriangularSolveLL1<NORM> (T11, X1);
    X2 -= T21 * X1;
    TriangularSolveLL1<NORM> (T22, X2);
  }

  void TriangularSolveLL (BareSliceMatrix<double> L, SliceMatrix<double> X)
  {
    TriangularSolveLL1<NonNormalized> (L,X);
  }
  void TriangularSolveLLN (BareSliceMatrix<double> L, SliceMatrix<double> X)
  {
    TriangularSolveLL1<Normalized> (L,X);
  }




  



  template <TRIG_NORMAL NORM>
  void TriangularSolveUR3 (BareSliceMatrix<double> bU, SliceMatrix<double> X)
  {
    size_t n = X.Height();
    auto U = bU.AddSize(n,n);
    
    constexpr size_t HA = 4;

    size_t i = n;
    for ( ; i >= HA; i -= HA)
      {
        if (i < n)
          MatKernel2AddAB<HA,SUB,RowMajor> (n-i, X.Width(), &U(i-HA,i), U.Dist(), &X(i,0), X.Dist(), &X(i-HA,0), X.Dist());
        KernelTriangularSolve<UpperRight,NORM,RowMajor,4> (X.Width(), &U(i-HA,i-HA), U.Dist(), &X(i-HA,0), X.Dist());
      }
    
    size_t remainder = i;
    if (remainder > 0)
      {
        Switch<HA> (remainder, [U,X,i,n] (auto r)
                    {
                      if constexpr (r.value > 0) {
                          if (i < n)
                            MatKernel2AddAB<r.value,SUB,RowMajor> (n-i, X.Width(),
                                                                   &U(0,i), U.Dist(),
                                                                   &X(i,0), X.Dist(),
                                                                   &X(0,0), X.Dist());
                          
                          KernelTriangularSolve<UpperRight,NORM,RowMajor,r.value> (X.Width(), &U(0,0), U.Dist(), &X(0,0), X.Dist());
                        }
                    });
      }
  }

  template <TRIG_NORMAL NORM>  
  void TriangularSolveUR2 (BareSliceMatrix<double> L, SliceMatrix<double> X)  
  {
    size_t i = 0;
    constexpr size_t bw = 128;
    for ( ; i+bw <= X.Width(); i += bw)
      TriangularSolveUR3<NORM> (L, X.Cols(i,i+bw));
    if (i < X.Width())
      TriangularSolveUR3<NORM> (L, X.Cols(i,X.Width()));      
  }
  
  template <TRIG_NORMAL NORM>
  void TriangularSolveUR1 (BareSliceMatrix<double> T, SliceMatrix<double> X)
  {
    size_t n = X.Height();
    if (n < 128)
      {
        TriangularSolveUR2<NORM> (T,X);
        return;
      }
    IntRange r1(0,n/2), r2(n/2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    auto T12 = T.Rows(r1).Cols(r2).AddSize(r1.Size(), r2.Size());
    // auto T21 = T.Rows(r2).Cols(r1).AddSize(r2.Size(), r1.Size());
    auto T22 = T.Rows(r2).Cols(r2);
    auto X1 = X.Rows(r1);
    auto X2 = X.Rows(r2);

    TriangularSolveUR1<NORM> (T22, X2);
    X1 -= T12 * X2;
    TriangularSolveUR1<NORM> (T11, X1);
  }

  void TriangularSolveUR (BareSliceMatrix<double> U, SliceMatrix<double> X)
  {
    TriangularSolveUR1<NonNormalized> (U,X);
  }
  void TriangularSolveURN (BareSliceMatrix<double> U, SliceMatrix<double> X)
  {
    TriangularSolveUR1<Normalized> (U,X);
  }






  /* ****************************** MultTriangular ************************/ 

  /* **************** LowerLeft - Normalized *************************** */
  

  void MultTriangularLLN3 (SliceMatrix<double> X, BareSliceMatrix<double> T)
  {
    size_t n = X.Width();
    size_t m = X.Height();
    /*
    for (size_t i = 0; i < n; i++)
      for (size_t j = i+1; j < n; j++)
        X.Col(i) += T(j,i) * X.Col(j);
    */
    /*
    for (size_t i = 0; i < n; i++)
      X.Col(i) += X.Cols(i+1,n) * T.Col(i).Range(i+1,n); 
    */
    constexpr size_t BS = 2*SIMD<double>::Size();    
    alignas (SIMD<double>) double memb[BS*128];
    size_t i = 0;
    for ( ; i+BS <= n; i+= BS)
      {
        auto Xrest = X.Cols(i,n);
        // copy b
        for (size_t j = 0; j < n-i; j++)
          for (size_t k = 0; k < BS; k++)
            memb[BS*j+k] = T(i+j,i+k);
        for (int i = 0; i < BS; i++)
          for (int j = i; j < BS; j++)
            memb[BS*i+j] = 0;
        
        FlatMatrix<> subb(n-i, BS, &memb[0]);
        /*
        for (size_t j = 0; j < m; j++)
          {
            Vec<BS> c = Trans(subb) * Xrest.Row(j);
            Xrest.Row(j).Range(0,BS) += c;
          }
        */
        size_t j = 0;
        for ( ; j+4 <= m; j+=4)
          MatKernelMultAB<4,2,ADD> (n-i, &Xrest(j,0), X.Dist(), &memb[0], BS, &Xrest(j,0), X.Dist());
        for ( ; j+1 <= m; j+=1)
          MatKernelMultAB<1,2,ADD> (n-i, &Xrest(j,0), X.Dist(), &memb[0], BS, &Xrest(j,0), X.Dist());
      }
    for ( ; i < n; i++)
      X.Col(i) += X.Cols(i+1,n) * T.Col(i).Range(i+1,n); 
  }
  
  void MultTriangularLLN2 (SliceMatrix<double> X, BareSliceMatrix<double> T)
  {
    size_t n = X.Width();
    if (n <= 128)
      {
        MultTriangularLLN3 (X,T);
        return;
      }
    IntRange r1(0,n/2), r2(n/2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    // auto T12 = T.Rows(r1).Cols(r2).AddSize(r1.Size(), r2.Size());
    auto T21 = T.Rows(r2).Cols(r1).AddSize(r2.Size(), r1.Size());
    auto T22 = T.Rows(r2).Cols(r2);
    auto X1 = X.Cols(r1);
    auto X2 = X.Cols(r2);

    MultTriangularLLN2 (X1,T11);
    X1 += X2 * T21;
    MultTriangularLLN2 (X2,T22);
  }
  
  void MultTriangularLLN (SliceMatrix<double> X, BareSliceMatrix<double> L)  
  {
    size_t i = 0;
    constexpr size_t bw = 128;
    for ( ; i+bw <= X.Height(); i += bw)
      MultTriangularLLN2 (X.Rows(i,i+bw), L);
    if (i < X.Height())
      MultTriangularLLN2 (X.Rows(i,X.Height()), L);      
  }






  /* **************** UpperRight - NonNormalized *************************** */

  

  void MultTriangularUR3 (SliceMatrix<double> X, BareSliceMatrix<double> T)
  {
    size_t n = X.Width();
    size_t m = X.Height();

    /*
    for (size_t i = n; i-->0; )
      {
        X.Col(i) *= T(i,i);
        for (size_t j = 0; j < i; j++)
          X.Col(i) += T(j,i) * X.Col(j);
      }
    */

    constexpr size_t BS = 2*SIMD<double>::Size();
    alignas (SIMD<double>) double memb[BS*128];
    size_t i = n;

    for ( ; i >= BS; i -= BS)
      {
        // copy b
        for (size_t j = 0; j < i; j++)
          for (size_t k = 0; k < BS; k++)
            memb[BS*j+k] = T(j,i+k-BS);
        for (int k = 0; k < BS; k++)
          for (int j = 0; j < k; j++)
            memb[BS*(i-BS+k)+j] = 0;

        /*
        FlatMatrix<> subb(i, BS, &memb[0]);
        for (size_t j = 0; j < m; j++)
          {
            Vec<BS> c = Trans(subb) * X.Row(j).Range(0,i);
            X.Row(j).Range(i-BS,i) = c;
          }
        */
        size_t j = 0;
        for ( ; j+4 <= m; j+=4)
          MatKernelMultAB<4,2,SET> (i, &X(j,0), X.Dist(), &memb[0], BS, &X(j,i-BS), X.Dist());
        for ( ; j+1 <= m; j+=1)
          MatKernelMultAB<1,2,SET> (i, &X(j,0), X.Dist(), &memb[0], BS, &X(j,i-BS), X.Dist());
      }

    for ( ; i >= 1; i--)
      {
        X.Col(i-1) *= T(i-1, i-1);
        if (i > 1)
          X.Col(i-1) += X.Cols(0, i-1) * T.Col(i-1).Range(0, i-1);
      }
  }
  
  void MultTriangularUR2 (SliceMatrix<double> X, BareSliceMatrix<double> T)
  {
    size_t n = X.Width();
    if (n <= 128)
      {
        MultTriangularUR3 (X,T);
        return;
      }
    IntRange r1(0,n/2), r2(n/2,n);
    auto T11 = T.Rows(r1).Cols(r1);
    auto T12 = T.Rows(r1).Cols(r2).AddSize(r1.Size(), r2.Size());
    // auto T21 = T.Rows(r2).Cols(r1).AddSize(r2.Size(), r1.Size());
    auto T22 = T.Rows(r2).Cols(r2);
    auto X1 = X.Cols(r1);
    auto X2 = X.Cols(r2);

    MultTriangularUR2 (X2,T22);
    X2 += X1 * T12;
    MultTriangularUR2 (X1,T11);
  }

  void MultTriangularUR (SliceMatrix<double> X, BareSliceMatrix<double> U)  
  {
    // MultTriangularUR3 (X, L);
    // return;
    
    size_t i = 0;
    constexpr size_t bw = 128;
    for ( ; i+bw <= X.Height(); i += bw)
      MultTriangularUR2 (X.Rows(i,i+bw), U);
    if (i < X.Height())
      MultTriangularUR2 (X.Rows(i,X.Height()), U);      
  }






  /* **************************** Generalized Triangular ************************** */

  template <TRIG_NORMAL NORM, ORDERING OT, ORDERING OXY>
  void GeneralizedTriangularMultLL (SliceMatrix<double, OT> T,
                                    SliceMatrix<double, OXY> X,
                                    SliceMatrix<double, OXY> Y)
  {
    auto [Y1,Y2] = Y.SplitRows(X.Height());
    auto [T1,T2] = T.SplitRows(X.Height());
    
    Y1 = X;
    TriangularMult<LowerLeft, NORM> (T1, Y1);
    Y2 = T2 * X;
  }

  template <TRIG_NORMAL NORM, ORDERING OT, ORDERING OXY>
  void GeneralizedTriangularMultUR (SliceMatrix<double, OT> T,
                                    SliceMatrix<double, OXY> X,
                                    SliceMatrix<double, OXY> Y)
  {
    // static Timer t ("TriangularMultUR"+ToString(OT)); 
    // static Timer tcopy ("TriangularMultUR, copy");
    
    if constexpr (OXY == RowMajor) {
        if (T.Width() <= 96)   
          {
            // RegionTimer reg(t);
            constexpr size_t WX = 96;
            alignas (64) double memx[96*96];
            
            size_t m = T.Height();
            size_t n = T.Width();
            // t.AddFlops ( (n*m-m*m/2) * X.Width());
            
            if (n < m) throw Exception ("generictrig UR should be wide");
            for (size_t j = 0; j < X.Width(); j += WX)
              {
                IntRange cols(j, min(j+WX, X.Width()));
                size_t roundup = cols.Size() + SIMD<double>::Size()-1;
                roundup -= roundup & (SIMD<double>::Size()-1);
                SliceMatrix<double,OXY> hx(X.Height(), cols.Size(), roundup, &memx[0]);
                bool copyx = X.Dist() >= 256;
                // tcopy.Start();
                if (copyx)
                  CopyMatrixIn (X.Height(), cols.Size(), &X(cols.First(),0), X.Dist(),
                                (SIMD<double>*) hx.Data(), hx.Dist() / SIMD<double>::Size());
                  // hx = X.Cols(cols);
                // tcopy.Stop();

                auto Xc = copyx ? hx : X.Cols(cols);
                auto Yc = Y.Cols(cols);

                constexpr size_t HA = 4;

                size_t i = 0;
                for ( ; i + HA <= m; i += HA)
                  {
                    KernelTriangularMultXY<UpperRight,NORM,OT,SET,HA> (Xc.Width(), &T(i,i), T.Dist(), &Xc(i,0), Xc.Dist(), &Yc(i,0), Yc.Dist());
                    if (i+HA < n)
                      MatKernel2AddAB<HA,ADD,OT> (n-i-HA, Xc.Width(), &T(i,i+HA), T.Dist(), &Xc(i+HA,0), Xc.Dist(), &Yc(i,0), Yc.Dist());
                  }

                size_t remainder = m % HA;
                Switch<HA> (remainder, [i,n,T,Xc,Yc] (auto r)
                            {
                              if constexpr (r.value > 0) {
                                  KernelTriangularMultXY<UpperRight,NORM,OT,SET,r.value> (Xc.Width(), &T(i,i), T.Dist(), &Xc(i,0), Xc.Dist(), &Yc(i,0), Yc.Dist());
                                  if (i+r.value < n)
                                    MatKernel2AddAB<r.value,ADD,OT> (n-i-r.value, Xc.Width(), &T(i,i+r.value), T.Dist(), &Xc(i+r.value,0), Xc.Dist(), &Yc(i,0), Yc.Dist());
                                }
                            });
              }
            return;
          }
      }
    
    
    auto [X1,X2] = X.SplitRows(T.Height());
    auto [T1,T2] = T.SplitCols(T.Height());
    
    // static Timer tc("generalizedtrig, copy trig, OXY="+ToString(OXY));
    // tc.Start();
    Y = X1;
    // tc.Stop();
    // static Timer t("generalizedtrig, trig part, side = UpperRight, normalized = "+ToString(NORM)+" OT = "+ToString(OT));
    // t.Start();
    TriangularMult<UpperRight, NORM> (T1, Y);
    // t.Stop();
    Y += T2 * X2;
  }

  template <>
  void GeneralizedTriangularMult_SM<LowerLeft,Normalized> (SliceMatrix<double,ColMajor> T,
                                                          SliceMatrix<double,RowMajor> X,
                                                          SliceMatrix<double,RowMajor> Y)
  { GeneralizedTriangularMultLL<Normalized> (T,X,Y); }

  template <>
  void GeneralizedTriangularMult_SM<LowerLeft,Normalized> (SliceMatrix<double,RowMajor> T,
                                                          SliceMatrix<double,RowMajor> X,
                                                          SliceMatrix<double,RowMajor> Y)
  { GeneralizedTriangularMultLL<Normalized> (T,X,Y); }

  template <>
  void GeneralizedTriangularMult_SM<LowerLeft,Normalized> (SliceMatrix<double,ColMajor> T,
                                                          SliceMatrix<double,ColMajor> X,
                                                          SliceMatrix<double,ColMajor> Y)
  { GeneralizedTriangularMultLL<Normalized> (T,X,Y); }

  template <>
  void GeneralizedTriangularMult_SM<LowerLeft,Normalized> (SliceMatrix<double,RowMajor> T,
                                                          SliceMatrix<double,ColMajor> X,
                                                          SliceMatrix<double,ColMajor> Y)
  { GeneralizedTriangularMultLL<Normalized> (T,X,Y); }

  template <>
  void GeneralizedTriangularMult_SM<UpperRight,Normalized> (SliceMatrix<double,ColMajor> T,
                                                          SliceMatrix<double,RowMajor> X,
                                                          SliceMatrix<double,RowMajor> Y)
  { GeneralizedTriangularMultUR<Normalized> (T,X,Y); }

  template <>
  void GeneralizedTriangularMult_SM<UpperRight,Normalized> (SliceMatrix<double,RowMajor> T,
                                                          SliceMatrix<double,RowMajor> X,
                                                          SliceMatrix<double,RowMajor> Y)
  { GeneralizedTriangularMultUR<Normalized> (T,X,Y); }

  template <>
  void GeneralizedTriangularMult_SM<UpperRight,Normalized> (SliceMatrix<double,ColMajor> T,
                                                          SliceMatrix<double,ColMajor> X,
                                                          SliceMatrix<double,ColMajor> Y)
  { GeneralizedTriangularMultUR<Normalized> (T,X,Y); }

  template <>
  void GeneralizedTriangularMult_SM<UpperRight,Normalized> (SliceMatrix<double,RowMajor> T,
                                                          SliceMatrix<double,ColMajor> X,
                                                          SliceMatrix<double,ColMajor> Y)
  { GeneralizedTriangularMultUR<Normalized> (T,X,Y); }


  

  /*
  template <TRIG_SIDE SIDE, TRIG_NORMAL NORM=NonNormalized, ORDERING OT, ORDERING OXY>
  void GeneralizedTriangularMult_SM (SliceMatrix<double, OT> T,
                                     SliceMatrix<double, OXY> X,
                                     SliceMatrix<double, OYY> Y)
  {
    if constexpr (SIDE == LowerLeft) {
        
        auto [Y1,Y2] = Y.SplitRows(X.Height());
        auto [T1,T2] = T.SplitRows(X.Height());
        
        Y1 = X;
        TriangularMult<SIDE, NORM> (T1, Y1);
        Y2 = T2 * X;
      }
    else {

      auto [X1,X2] = X.SplitRows(T.Height());
      auto [T1,T2] = T.SplitCols(T.Height());

      static Timer tc("generalizedtrig, copy trig, OX="+ToString(OX)+" OY="+ToString(OY));
      tc.Start();
      Y = X1;
      tc.Stop();
      static Timer t("generalizedtrig, trig part, SIDE = "+ToString(SIDE)+" normalized = "+ToString(NORM)+" OT = "+ToString(OT));
      t.Start();
      TriangularMult<SIDE, NORM> (T1, Y);
      t.Stop();
      Y += T2 * X2;
    }
  }
  */






  
  // T  .... lower left
  // Y -= T * X
  template <TRIG_NORMAL NORM, ORDERING OT, ORDERING OXY>
  void GeneralizedTriangularSubLL (SliceMatrix<double, OT> T,
                                   SliceMatrix<double, OXY> X,
                                   SliceMatrix<double, OXY> Y)
  {
    // static Timer t("GeneraliedTrigSubLL,OT = " + ToString(OT) + " OXY = " + ToString(OXY)); RegionTimer reg(t);
    // static Timer trect("rect");
    auto [T1,T2] = T.SplitRows(T.Width());
    auto [Y1,Y2] = Y.SplitRows(T.Width());

    if constexpr (OXY == ColMajor)
                   {
                     static Timer ttrig("trig,LL,generic"); RegionTimer reg(ttrig);
                     Matrix<double,OXY> tmpX = X;
                     TriangularMult<LowerLeft,NORM> (T1, tmpX);
                     Y1 -= tmpX;
                   }
    else
      {
        // static Timer ttrig("trig,LL"); RegionTimer reg(ttrig);
        // ttrig.AddFlops(sqr(X.Height())/2 * X.Width());
        /*
        Matrix<double,OXY> tmpX = X;
        TriangularMult<LowerLeft,NORM> (T1, tmpX);
        Y1 -= tmpX;
        */
        constexpr size_t WX = 96;

        for (size_t j = 0; j < X.Width(); j += WX)
          {
            IntRange cols(j, min(j+WX, X.Width()));

            auto Xc = X.Cols(cols);
            auto Yc = Y.Cols(cols);
            
            size_t n = T.Width();
            constexpr size_t HA = 4;
            size_t remainder = n % HA;
        
            Switch<HA> (remainder, [T,Xc,Yc] (auto r)
                        {
                          KernelTriangularMultXY<LowerLeft,NORM,OT,SUB,r.value> (Xc.Width(), &T(0,0), T.Dist(),
                                                                                 &Xc(0,0), Xc.Dist(), &Yc(0,0), Yc.Dist());
                        });
            
            for (size_t i = remainder; i < n; i += HA)
              {
                KernelTriangularMultXY<LowerLeft,NORM,OT,SUB,HA> (Xc.Width(), &T(i,i), T.Dist(), &Xc(i,0), Xc.Dist(), &Yc(i,0), Yc.Dist());
                if (i > 0)
                  MatKernel2AddAB<HA,SUB,OT> (i, Xc.Width(), &T(i,0), T.Dist(), &Xc(0,0), Xc.Dist(), &Yc(i,0), Yc.Dist());
              }
          }
      }

    {
      // RegionTimer regrect(trect);
      // trect.AddFlops(T2.Height()*T2.Width()*X.Width());
      Y2 -= T2 * X;
    }
  }
  
  template <>
  void GeneralizedTriangularSub_SM<LowerLeft,Normalized> (SliceMatrix<double,ColMajor> T,
                                                          SliceMatrix<double,RowMajor> X,
                                                          SliceMatrix<double,RowMajor> Y)
  { GeneralizedTriangularSubLL<Normalized> (T,X,Y); }

  template <>
  void GeneralizedTriangularSub_SM<LowerLeft,Normalized> (SliceMatrix<double,RowMajor> T,
                                                          SliceMatrix<double,RowMajor> X,
                                                          SliceMatrix<double,RowMajor> Y)
  { GeneralizedTriangularSubLL<Normalized> (T,X,Y); }

  template <>
  void GeneralizedTriangularSub_SM<LowerLeft,Normalized> (SliceMatrix<double,ColMajor> T,
                                                          SliceMatrix<double,ColMajor> X,
                                                          SliceMatrix<double,ColMajor> Y)
  { GeneralizedTriangularSubLL<Normalized> (T,X,Y); }

  template <>
  void GeneralizedTriangularSub_SM<LowerLeft,Normalized> (SliceMatrix<double,RowMajor> T,
                                                          SliceMatrix<double,ColMajor> X,
                                                          SliceMatrix<double,ColMajor> Y)
  { GeneralizedTriangularSubLL<Normalized> (T,X,Y); }
  
  
}


