/****************************************************************************/
/* File:   LUdecomposition.cpp                                              */
/* Author: Joachim Schoeberl                                                */
/* Date:   Nov 2020                                                         */
/****************************************************************************/

#include <core/simd.hpp>
#include <ngs_stdcpp_include.hpp>
#include <matrix.hpp>
#include <vector.hpp>
#include <triangular.hpp>
#include <ngblas.hpp>


namespace ngbla
{

  /*
    void CalcLU (SliceMatrix<double> a, FlatArray<int> p)
    {
    size_t n = a.Height();
  
    for (size_t i = 0; i < n; i++) p[i] = i;
  
    for (size_t i = 0; i < n; i++)
    {
    size_t imax = i;
    for (size_t j = i+1; j < n; j++)
    if (fabs(a(j,i)) > fabs(a(imax,i)))
    imax = j;
      
    if (imax != i)
    {
    Swap (p[i], p[imax]);
    for (int j = 0; j < n; j++)
    Swap (a(i,j), a(imax, j));
    }
      
    a.Row(i).Range(i,n) -= Trans (a.Rows(0,i).Cols(i,n)) * a.Row(i).Range(0,i);
    a.Col(i).Range(i+1,n) -= a.Cols(0,i).Rows(i+1,n) * a.Col(i).Range(0,i);
    a.Col(i).Range(i+1,n) *= 1.0/a(i,i);
    }
    }
  */



  void SwapVectors (FlatVector<> a, BareVector<> b)
  {
    size_t n = a.Size();
    /*
      for (size_t i = 0; i < n; i++)
      Swap (a(i), b(i));
    */
  
    size_t i = 0;
    for ( ; i+SIMD<double>::Size() <= n; i+=SIMD<double>::Size())
      {
        auto va = SIMD<double>(&a(i));
        auto vb = SIMD<double>(&b(i));
        va.Store(&b(i));
        vb.Store(&a(i));
      }

    if (i==n) return;

    // handle SIMD rest
    Switch<SIMD<double>::Size()> ( (n-i), [&] (auto r)
                                   {
                                     Vec<r.value> ha = a.Range(i, i+r.value);
                                     Vec<r.value> hb = b.Range(i, i+r.value);
                                     a.Range(i,i+r.value) = hb;
                                     b.Range(i,i+r.value) = ha;
                                   });

  }
  

  
  void CalcLU1 (SliceMatrix<double> a, FlatArray<int> p)
  {
    size_t n = a.Height();

    static Timer t("CalcLU"); RegionTimer reg(t);
    t.AddFlops (n*n*n/3);

    // static Timer tmm1("CalcLU - mat mat 1"); 
    // static Timer tmm2("CalcLU - mat mat 2"); 
    // static Timer tmm3("CalcLU - mat mat 3"); 

  
    for (size_t i = 0; i < n; i++) p[i] = i;
  
    size_t bs = 48;
    size_t bs2 = 8;
    for (size_t i1 = 0; i1 < n; i1+=bs)
      {
        size_t end1 = min(n, i1+bs);

        for (size_t i2 = i1; i2 < end1; i2 += bs2)
          {
            size_t end2 = min(n, i2+bs2);
            
            for (size_t i = i2; i < end2; i++)
              {
                size_t imax = i;
                double valmax = fabs(a(i,i));
              
                for (size_t j = i+1; j < n; j++)
                  if (double valj = fabs(a(j,i)) > valmax)
                    {
                      valmax = valj;
                      imax = j;
                    }
              
                if (imax != i)
                  {
                    Swap (p[i], p[imax]);
                    SwapVectors (a.Row(i), a.Row(imax));
                  }
              
                if (i+1 < n)
                  {
                    a.Col(i).Range(i+1,n) *= 1.0/a(i,i);
                    // RegionTimer rmm3(tmm3);                                  
                    a.Rows(i+1,n).Cols(i+1,end2) -= a.Rows(i+1,n).Cols(i,i+1) * a.Rows(i,i+1).Cols(i+1,end2);
                  }
              }

            // RegionTimer rmm2(tmm2);                
            if (end2 < end1)
              {
                TriangularSolve<LowerLeft,Normalized> (a.Rows(i2, end2).Cols(i2, end2), a.Rows(i2, end2).Cols(end2, end1));
                a.Rows(end2,n).Cols(end2, end1) -= a.Cols(i2,end2).Rows(end2,n) * a.Rows(i2,end2).Cols(end2, end1);
              }
          }

        // RegionTimer rmm1(tmm1);      
        if (end1 < n)
          {
            TriangularSolve<LowerLeft,Normalized> (a.Rows(i1, end1).Cols(i1, end1), a.Rows(i1, end1).Cols(end1, n));
            a.Rows(end1,n).Cols(end1, n) -= a.Cols(i1,end1).Rows(end1,n) * a.Rows(i1,end1).Cols(end1, n);
          }
      }
  }


  /*
  static Timer calcLUSolveL("CalcLU - SolveL");
  static Timer calcLUMatMat("CalcLU - MatMat");
  static Timer calcLUSimple("CalcLU - simple");  
  static Timer calcLUSimple_simd("CalcLU - simple SIMD");  
  static Timer calcLUSimple2("CalcLU - simple2x");  
  static Timer calcLUSearch("CalcLU - search");  
  static Timer calcLUSwap("CalcLU - swap");  
  */
  void CalcLURec (SliceMatrix<double> a, FlatArray<int> p, IntRange r)
  {
    size_t n = a.Height();
    
    /*
    if (r.Size() == 0) return;
    if (r.Size() == 1)
      {
        size_t i = r.First();
        
        size_t imax = i;
        double valmax = fabs(a(i,i));
      
        for (size_t j = i+1; j < n; j++)
          if (double valj = fabs(a(j,i)) > valmax)
            {
              valmax = valj;
              imax = j;
            }
        
        if (imax != i)
          {
            Swap (p[i], p[imax]);
            SwapVectors (a.Row(i), a.Row(imax));
          }
        
        if (i+1 < n)
          a.Col(i).Range(i+1,n) *= 1.0/a(i,i);
        return;
      }
    */
    constexpr size_t bs = 4;
    if (r.Size() <= bs)
      {
        // RegionTimer reg(calcLUSimple);
        for (auto i : r)
          {
            size_t imax = i;
            double valmax = fabs(a(i,i));

            {
              // RegionTimer reg(calcLUSearch);
            for (size_t j = i+1; j < n; j++)
              {
                double valj = fabs(a(j,i));
                if (valj > valmax)
                  {
                    valmax = valj;
                    imax = j;
                  }
              }
            }
            // singular if the pivot is negligible relative to its own row
            // (row-scaling invariant); also catches an exact 0 and NaN
            double rowmax = 0;
            for (size_t j = i+1; j < r.Next(); j++)
              rowmax = max(rowmax, fabs(a(imax,j)));
            if (!(valmax > 1e-20*rowmax))
              throw Exception ("CalcLU: matrix singular");
            if (imax != i)
              {
                // RegionTimer reg(calcLUSwap);                
                Swap (p[i], p[imax]);
                SwapVectors (a.Row(i), a.Row(imax));
              }


            if (r.Size() == 4)
              {
                // RegionTimer reg(calcLUSimple_simd);                                
                size_t rest = i-r.First();
                double * ptr = &a(i, r.First());
                
                double invaii = 1.0/ptr[rest];

                /*
                double scale[4] = { 0, 0, 0, 0 };
                
                // for (size_t k = 0; k < rest; k++)
                // scale[k] = 0.0;
                scale[rest] = invaii-1;
                for (size_t k = rest+1; k < 4; k++)
                  scale[k] = -a(i,r.First()+k)*invaii;
                SIMD<double,4> scale1(&scale[0]);
                */

                // pivot row right of the pivot, zero elsewhere; the multiplier
                // a(j,i)/a(i,i) is computed directly and stored (folding the
                // scaling into the update as a(j,i)*(1/a(i,i)-1) cancels for
                // large pivots and loses eps*|a(i,i)| of accuracy)
                SIMD<mask64,4> m2(rest+1);
                SIMD<double,4> urow = If (m2, SIMD<double,4>(0.0), SIMD<double,4>(ptr));
                ptr += a.Dist();
                for (size_t j = i+1; j < n; j++, ptr+=a.Dist())
                  {
                    SIMD<double,4> row1(ptr);
                    double fac = ptr[rest]*invaii;
                    row1 -= fac*urow;
                    row1.Store(ptr);
                    ptr[rest] = fac;
                  }
              }
            /*
            else if (r.Size() == 8)
              {
                double scale[8];
                
                size_t rest = i-r.First();
                for (size_t k = 0; k < rest; k++)
                  scale[k] = 0.0;
                double invaii = 1.0/a(i,i);
                scale[rest] = invaii-1;
                for (size_t k = rest+1; k < 8; k++)
                  scale[k] = -a(i,r.First()+k)*invaii;

                SIMD<double,4> scale1(&scale[0]), scale2(&scale[4]);

                double * ptr = &a(i+1, r.First());
                for (size_t j = i+1; j < n; j++, ptr+=a.Dist())
                  {
                    SIMD<double,4> row1(ptr);
                    SIMD<double,4> row2(ptr+4);
                    double fac = ptr[rest];
                    row1 += fac*scale1;
                    row2 += fac*scale2;
                    row1.Store(ptr);
                    row2.Store(ptr+4);
                  }
              }
            */
            else
              {
                if (i+1 < n)
                  a.Col(i).Range(i+1,n) *= 1.0/a(i,i);
                if (i+1 < r.Next())
                  {
                    // RegionTimer reg(calcLUSimple2);                
                    a.Rows(i+1,n).Cols(i+1,r.Next()) -= a.Rows(i+1,n).Cols(i,i+1) * a.Rows(i,i+1).Cols(i+1,r.Next());
                    /*
                      double mem[bs];
                      FlatMatrix<> row(1, r.Size(), &mem[0]);
                      row.Row(0) = a.Row(i).Range(r);
                      row.Row(0).Range(i-r.First()+1) = 0.0;
                      a.Rows(i+1,n).Cols(r) -= a.Rows(i+1,n).Cols(i,i+1) * row;
                    */
                  }
              }
          }

        return;
      }
    
    size_t half = r.Size()/2;
    if (half > bs)
      half = half - (half % bs);
    size_t mid = r.First() + half;
    IntRange r1(r.First(), mid);
    IntRange r2(mid, r.Next());
    
    CalcLURec (a, p, r1);
    
    {
      // RegionTimer r(calcLUSolveL);
      TriangularSolve<LowerLeft,Normalized> (a.Rows(r1).Cols(r1), a.Rows(r1).Cols(r2));
    }
    
    {
      // RegionTimer r(calcLUMatMat);
      a.Rows(mid,n).Cols(r2) -= a.Rows(mid,n).Cols(r1) * a.Rows(r1).Cols(r2);
    }
    CalcLURec (a, p, r2);
  }



  void CalcLU (SliceMatrix<double> a, FlatArray<int> p)
  {
    size_t n = a.Height();
    
    // static Timer t("CalcLU - rec"); RegionTimer reg(t);
    // t.AddFlops (n*n*n/3);
    
    for (size_t i = 0; i < n; i++) p[i] = i;
    
    CalcLURec (a, p, IntRange(n));
  }


  inline double CAbs1 (Complex z) { return fabs(z.real())+fabs(z.imag()); }

  // as the double version
  void CalcLURec (SliceMatrix<Complex> a, FlatArray<int> p, IntRange r)
  {
    size_t n = a.Height();
    constexpr size_t bs = 4;
#ifndef __aarch64__
    // x86 only: no gain over the scalar leaf on arm
    if (r.Size() == bs)
      {
        // simd leaf: the 4 complex entries of a row in one SIMD<double,8>,
        // the pivot search for the next column is fused into the update sweep
        size_t f = r.First();
        size_t dist = 2*a.Dist();
        double sgnmem[8] = { -1, 1, -1, 1, -1, 1, -1, 1 };
        SIMD<double,8> sgn(&sgnmem[0]);

        size_t imax = f;
        double valmax = CAbs1(a(f,f));
        for (size_t j = f+1; j < n; j++)
          if (double valj = CAbs1(a(j,f)); valj > valmax)
            {
              valmax = valj;
              imax = j;
            }

        for (size_t i = f; i < f+bs; i++)
          {
            double rowmax = 0;
            for (size_t j = i+1; j < r.Next(); j++)
              rowmax = max(rowmax, CAbs1(a(imax,j)));
            if (!(CAbs1(a(imax,i)) > 1e-20*rowmax))
              throw Exception ("CalcLU: matrix singular");
            if (imax != i)
              {
                Swap (p[i], p[imax]);
                for (size_t k = 0; k < a.Width(); k++)
                  Swap (a(i,k), a(imax,k));
              }

            size_t rest = i-f;
            Complex invaii = 1.0/a(i,i);
            double * ptr = (double*)(void*)&a(i,f);
            SIMD<double,8> urow = If (SIMD<mask64,8>(int64_t(2*rest+2)), SIMD<double,8>(0.0), SIMD<double,8>(ptr));
            SIMD<double,8> urows = sgn * SwapPairs(urow);
            bool search = i+1 < f+bs;
            double vmax = -1;
            ptr += dist;
            for (size_t j = i+1; j < n; j++, ptr += dist)
              {
                Complex fac = Complex(ptr[2*rest], ptr[2*rest+1]) * invaii;
                SIMD<double,8> row(ptr);
                row = FNMA (SIMD<double,8>(fac.imag()), urows, FNMA (SIMD<double,8>(fac.real()), urow, row));
                row.Store(ptr);
                ptr[2*rest] = fac.real();
                ptr[2*rest+1] = fac.imag();
                if (search)
                  if (double v = fabs(ptr[2*rest+2])+fabs(ptr[2*rest+3]); v > vmax)
                    {
                      vmax = v;
                      imax = j;
                    }
              }
          }
        return;
      }
#endif
    if (r.Size() <= bs)
      {
        for (auto i : r)
          {
            size_t imax = i;
            double valmax = CAbs1(a(i,i));
            for (size_t j = i+1; j < n; j++)
              {
                double valj = CAbs1(a(j,i));
                if (valj > valmax)
                  {
                    valmax = valj;
                    imax = j;
                  }
              }
            double rowmax = 0;
            for (size_t j = i+1; j < r.Next(); j++)
              rowmax = max(rowmax, CAbs1(a(imax,j)));
            if (!(valmax > 1e-20*rowmax))
              throw Exception ("CalcLU: matrix singular");
            if (imax != i)
              {
                Swap (p[i], p[imax]);
                for (size_t k = 0; k < a.Width(); k++)
                  Swap (a(i,k), a(imax,k));
              }

            Complex invaii = 1.0/a(i,i);
            for (size_t j = i+1; j < n; j++)
              {
                Complex fac = a(j,i)*invaii;
                a(j,i) = fac;
                for (size_t k = i+1; k < r.Next(); k++)
                  a(j,k) -= fac*a(i,k);
              }
          }
        return;
      }

    size_t half = r.Size()/2;
    if (half > bs)
      half = half - (half % bs);
    size_t mid = r.First() + half;
    IntRange r1(r.First(), mid);
    IntRange r2(mid, r.Next());

    CalcLURec (a, p, r1);
    TriangularSolve<LowerLeft,Normalized> (a.Rows(r1).Cols(r1), a.Rows(r1).Cols(r2));
    a.Rows(mid,n).Cols(r2) -= a.Rows(mid,n).Cols(r1) * a.Rows(r1).Cols(r2);
    CalcLURec (a, p, r2);
  }

  void CalcLU (SliceMatrix<Complex> a, FlatArray<int> p)
  {
    for (size_t i = 0; i < a.Height(); i++) p[i] = i;
    CalcLURec (a, p, IntRange(a.Height()));
  }



  


  // U .. upper right,
  // L .. lower left, normalized
  // static Timer tmulul1 ("MultUL - matmat");
  // static Timer tmulul2 ("MultUL - trigmultR");
  // static Timer tmulul3 ("MultUL - trigmultL");
  void MultUL (SliceMatrix<> A)
  {
    size_t n = A.Height();  
    if (n <= 1) return;

    if (n <= 8)
      {
        for (size_t i = 0; i < n; i++)
          {
            auto rowi = A.Row(i);
            for (size_t j = 0; j < i; j++)
              {
                double sum = 0;
                for (size_t k = i; k < n; k++)
                  sum += rowi(k) * A(k,j);
                rowi(j) = sum;
              }
            for (size_t j = i; j < n; j++)
              {
                double sum = rowi(j);
                for (size_t k = j+1; k < n; k++)
                  sum += rowi(k) * A(k,j);
                rowi(j) = sum;
              }
          }
        return;
      }

    
    IntRange r1(0,n/2), r2(n/2,n);
    auto A11 = A.Rows(r1).Cols(r1);
    auto A12 = A.Rows(r1).Cols(r2);
    auto A21 = A.Rows(r2).Cols(r1);
    auto A22 = A.Rows(r2).Cols(r2);

    MultUL (A11);
    // tmulul1.Start();
    A11 += A12 * A21;
    // tmulul1.Stop();
    // tmulul1.AddFlops (r1.Size()*r1.Size()*r2.Size());
    // tmulul2.Start();
    TriangularMult<UpperRight> (A22, A21);
    // tmulul2.Stop();
    // tmulul3.Start();
    // TriangularMult<UpperRight,Normalized> (Trans(A22), Trans(A12));
    MultTriangular<LowerLeft,Normalized> (A12, A22);
    // tmulul3.Stop();
    MultUL (A22);
  }



  // X := X * inv(L), L unit lower triangular m x m; recursion over columns keeps X row-wise
  template <typename T>
  void SolveXLUnit (SliceMatrix<T> X, BareSliceMatrix<T> L, size_t m)
  {
    if (m <= 4)
      {
        for (size_t i = 0; i < X.Height(); i++)
          {
            T * x = &X(i,0);
            for (size_t j = m-1; j-- > 0; )
              {
                T sum = x[j];
                for (size_t k = j+1; k < m; k++)
                  sum -= x[k] * L(k,j);
                x[j] = sum;
              }
          }
        return;
      }
    size_t m1 = m/2;
    SolveXLUnit<T> (X.Cols(m1,m), L.Rows(m1,m).Cols(m1,m), m-m1);
    X.Cols(0,m1) -= X.Cols(m1,m) * L.Rows(m1,m).Cols(0,m1).AddSize(m-m1,m1);
    SolveXLUnit<T> (X.Cols(0,m1), L, m1);
  }

  // X := X * U, U upper triangular m x m
  template <typename T>
  void MultXU (SliceMatrix<T> X, BareSliceMatrix<T> U, size_t m)
  {
    if (m <= 4)
      {
        for (size_t i = 0; i < X.Height(); i++)
          {
            T * x = &X(i,0);
            for (size_t j = m; j-- > 0; )
              {
                T sum = x[j] * U(j,j);
                for (size_t k = 0; k < j; k++)
                  sum += x[k] * U(k,j);
                x[j] = sum;
              }
          }
        return;
      }
    size_t m1 = m/2;
    MultXU<T> (X.Cols(m1,m), U.Rows(m1,m).Cols(m1,m), m-m1);
    X.Cols(m1,m) += X.Cols(0,m1) * U.Rows(0,m1).Cols(m1,m).AddSize(m1,m-m1);
    MultXU<T> (X.Cols(0,m1), U, m1);
  }

  // inv(U) in place; the generic MultTriangular works on transposed (strided) rows
  template <typename T>
  void TriangularInvertUR (SliceMatrix<T> U)
  {
    size_t n = U.Height();
    if (n < 16)
      {
        TriangularInvert<UpperRight> (U);
        return;
      }
    size_t n1 = n/2;
    auto U11 = U.Rows(0,n1).Cols(0,n1);
    auto U12 = U.Rows(0,n1).Cols(n1,n);
    auto U22 = U.Rows(n1,n).Cols(n1,n);
    TriangularInvertUR (U11);
    TriangularInvertUR (U22);
    U12 *= -1;
    TriangularMult<UpperRight> (U11, U12);
    MultXU<T> (U12, U22, n-n1);
  }

  template <typename T>
  void T_InverseFromLU (SliceMatrix<T> A, FlatArray<int> p)
  {
    size_t n = A.Height();

    /*
    // testing: lapack-version
    Matrix tmp = Trans(A);
    ArrayMem<integer,100> ipiv(A.Height());
    for (int i = 0; i < n; i++)
      ipiv[i] = i+1;
    integer lda = tmp.Dist();

    integer lwork = 32 * A.Height();
    Array<double> work(lwork);
    
    integer info;
    integer ni = A.Height();
    dgetri(&ni, &tmp(0,0), &lda, &ipiv[0], work.Data(), &lwork, &info);

    A = Trans(tmp);
    if (info != 0)
      cout << "info = " << info << endl;
    */

    // as lapack dgetri: invert U in place, then solve X L = inv(U) by
    // column blocks from the right. Inverting L explicitly and forming
    // inv(U)*inv(L) has the same cost but a worse residual for
    // ill-conditioned matrices (Higham, Accuracy and Stability, 14.3)
    if (is_same_v<T,Complex> && n >= 64)
      TriangularInvertUR (A);
    else
      TriangularInvert<UpperRight> (A);

    constexpr size_t bs = 32;
    ArrayMem<T, bs*128> mem(bs*n);
    for (size_t j1 = n; j1 > 0; )
      {
        size_t j0 = (j1 > bs) ? j1-bs : 0;
        IntRange J(j0, j1);
        size_t nb = j1-j0;
        // W = L(J, J) and L(j1:n, J), the strictly lower part of the columns J
        FlatMatrix<T> W(n-j0, nb, mem.Data());
        W = 0.0;
        for (size_t i = j0; i < n; i++)
          for (size_t j = j0; j < min(i, j1); j++)
            {
              W(i-j0, j-j0) = A(i,j);
              A(i,j) = 0.0;
            }
        for (size_t j = j0; j < j1; j++) W(j-j0, j-j0) = 1.0;
        // X(:,J) -= X(:, j1:n) * L(j1:n, J)
        if (j1 < n)
          A.Cols(J) -= A.Cols(j1, n) * W.Rows(nb, n-j0);
        // X(:,J) = X(:,J) * inv(L(J,J))
        SolveXLUnit<T> (A.Cols(J), W.Rows(0,nb), nb);
        j1 = j0;
      }

    // RegionTimer rperm(tperm);
    VectorMem<100,T> row(n);
    for (size_t i = 0; i < n; i++)
      {
        auto rowi = A.Row(i);
        for (size_t j = 0; j < n; j++)
          row(p[j]) = rowi(j);
        rowi = row;
      }
  }

  void InverseFromLU (SliceMatrix<double> A, FlatArray<int> p)
  {
    T_InverseFromLU (A, p);
  }

  void InverseFromLU (SliceMatrix<Complex> A, FlatArray<int> p)
  {
    T_InverseFromLU (A, p);
  }


  void SolveFromLU (SliceMatrix<double> A, FlatArray<int> p, SliceMatrix<double,ColMajor> X)
  {
    size_t n = X.Height();
    VectorMem<100,double> hv(n);
    for (size_t i = 0; i < X.Width(); i++)
      {
        auto coli = X.Col(i);
        hv = coli;
        for (size_t j = 0; j < n; j++)
          coli(j) = hv(p[j]);
      }
    
    TriangularSolve<LowerLeft,Normalized> (A, X);
    TriangularSolve<UpperRight> (A, X);
  }


  template <typename T>
  void T_SolveTransFromLU (SliceMatrix<T> A, FlatArray<int> p, SliceMatrix<T,ColMajor> X)
  {
    TriangularSolve<LowerLeft> (Trans(A).Bare(), X);
    TriangularSolve<UpperRight,Normalized> (Trans(A).Bare(), X);

    size_t n = X.Height();

    VectorMem<100,T> hv(n);
    for (size_t i = 0; i < X.Width(); i++)
      {
        auto coli = X.Col(i);
        hv = coli;
        for (size_t j = 0; j < n; j++)
          coli(p[j]) = hv(j);
      }
    
  }

  void SolveTransFromLU (SliceMatrix<double> A, FlatArray<int> p, SliceMatrix<double,ColMajor> X)
  {
    T_SolveTransFromLU (A, p, X);
  }

  void SolveTransFromLU (SliceMatrix<Complex> A, FlatArray<int> p, SliceMatrix<Complex,ColMajor> X)
  {
    T_SolveTransFromLU (A, p, X);
  }


  
}
