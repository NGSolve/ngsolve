/****************************************************************************/
/* File:   LUdecomposition.cpp                                              */
/* Author: Joachim Schoeberl                                                */
/* Date:   Nov 2020                                                         */
/****************************************************************************/


#include <bla.hpp>

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







  void CalcLURec (SliceMatrix<double> a, FlatArray<int> p, IntRange r)
  {
    size_t n = a.Height();
    
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
    
    size_t mid = r.First() + r.Size()/2;
    IntRange r1(r.First(), mid);
    IntRange r2(mid, r.Next());
    
    CalcLURec (a, p, r1);
    
    TriangularSolve<LowerLeft,Normalized> (a.Rows(r1).Cols(r1), a.Rows(r1).Cols(r2));
    a.Rows(mid,n).Cols(r2) -= a.Rows(mid,n).Cols(r1) * a.Rows(r1).Cols(r2);
    
    CalcLURec (a, p, r2);
  }



  void CalcLU (SliceMatrix<double> a, FlatArray<int> p)
  {
    size_t n = a.Height();
    
    static Timer t("CalcLU - rec"); RegionTimer reg(t);
    t.AddFlops (n*n*n/3);
    
    for (size_t i = 0; i < n; i++) p[i] = i;
    
    CalcLURec (a, p, IntRange(n));
  }
  



  


  // U .. upper right,
  // L .. lower left, normalized
  static Timer tmulul1 ("MultUL - matmat");
  static Timer tmulul2 ("MultUL - trigmultR");
  static Timer tmulul3 ("MultUL - trigmultL");
  void MultUL (SliceMatrix<> A)
  {
    size_t n = A.Height();  
    if (n <= 1) return;
  
    IntRange r1(0,n/2), r2(n/2,n);
    auto A11 = A.Rows(r1).Cols(r1);
    auto A12 = A.Rows(r1).Cols(r2);
    auto A21 = A.Rows(r2).Cols(r1);
    auto A22 = A.Rows(r2).Cols(r2);

    MultUL (A11);
    tmulul1.Start();
    A11 += A12 * A21;
    tmulul1.Stop();
    tmulul1.AddFlops (r1.Size()*r1.Size()*r2.Size());
    tmulul2.Start();
    TriangularMult<UpperRight> (A22, A21);
    tmulul2.Stop();
    tmulul3.Start();
    TriangularMult<UpperRight,Normalized> (Trans(A22), Trans(A12));
    tmulul3.Stop();
    MultUL (A22);
  }



  void InverseFromLU (SliceMatrix<double> A, FlatArray<int> p)
  {
    size_t n = A.Height();    
    static Timer t("InverseFromLU"); RegionTimer reg(t);
    t.AddFlops (2*n*n*n/3);

    static Timer tl("InvertL"); 
    static Timer tu("InvertU"); 
    static Timer tperm("permutation");
    
    tl.Start();
    TriangularInvert<LowerLeft,Normalized> (A);
    tl.Stop();
    tl.AddFlops (n*n*n/6);
    tu.Start();
    TriangularInvert<UpperRight> (A);
    tu.Stop();
    tu.AddFlops (n*n*n/6);    
    MultUL (A);

    RegionTimer rperm(tperm);
    VectorMem<100> row(n);
    for (size_t i = 0; i < n; i++)
      {
        auto rowi = A.Row(i);
        for (size_t j = 0; j < n; j++)
          row(p[j]) = rowi(j);
        rowi = row;
      }
  }





}
