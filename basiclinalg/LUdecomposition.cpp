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

                SIMD<double,4> scale1 = 0.0;
                SIMD<mask64,4> m1(rest);
                scale1 = If (m1, scale1, SIMD<double,4>(invaii-1));
                SIMD<mask64,4> m2(rest+1);
                scale1 = If (m2, scale1, -SIMD<double,4>(ptr)*invaii);
                ptr += a.Dist();
                for (size_t j = i+1; j < n; j++, ptr+=a.Dist())
                  {
                    SIMD<double,4> row1(ptr);
                    double fac = ptr[rest];
                    row1 += fac*scale1;
                    row1.Store(ptr);
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



  void InverseFromLU (SliceMatrix<double> A, FlatArray<int> p)
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

    // static Timer t("InverseFromLU"); RegionTimer reg(t);
    // t.AddFlops (2*n*n*n/3);

    // static Timer tl("InvertL"); 
    // static Timer tu("InvertU"); 
    // static Timer tperm("permutation");
    
    // tl.Start();
    TriangularInvert<LowerLeft,Normalized> (A);
    // tl.Stop();
    // tl.AddFlops (n*n*n/6);
    // tu.Start();
    TriangularInvert<UpperRight> (A);
    // tu.Stop();
    // tu.AddFlops (n*n*n/6);    
    MultUL (A);

    // RegionTimer rperm(tperm);
    VectorMem<100> row(n);
    for (size_t i = 0; i < n; i++)
      {
        auto rowi = A.Row(i);
        for (size_t j = 0; j < n; j++)
          row(p[j]) = rowi(j);
        rowi = row;
      }
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


  void SolveTransFromLU (SliceMatrix<double> A, FlatArray<int> p, SliceMatrix<double,ColMajor> X)
  {
    TriangularSolve<LowerLeft> (Trans(A), X);
    TriangularSolve<UpperRight,Normalized> (Trans(A), X);

    size_t n = X.Height();

    VectorMem<100,double> hv(n);
    for (size_t i = 0; i < X.Width(); i++)
      {
        auto coli = X.Col(i);
        hv = coli;
        for (size_t j = 0; j < n; j++)
          coli(p[j]) = hv(j);
      }
    
  }


  
}
