/*
  The singular value decomposition: Anatomy of optimizing an algorithm for extreme scale
  J Dongarra, M Gates, A Haidar, J Kurzak, P Luszczek, S Tomov, ...
  SIAM review 60 (4), 808-865 (2017)

  L. Hogben (ed): Handbook of Linear Algebra, Chap 45
  

  * use the two-step bidiagonalization using Householder reflections from both sides
  * bulge chasing (by B. Lang), an  O(N*N*bs) algorithm for diagonalization
  * bisection for singular values
  * twisted factorization (RRR relatively robust representation)
  
  TODO: 
  MRRR ... multiple RRR
  degegenerate cases (with intermediate zeros)
*/


#include <bla.hpp>

namespace ngbla
{

  /*
    int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, 
    doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
    ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
    integer *info);
  */

  void LapackSVD (SliceMatrix<> A,
                  SliceMatrix<double, ColMajor> U,
                  SliceMatrix<double, ColMajor> V)
  {
    static Timer t("LapackSVD"); RegionTimer reg(t);
    ngbla::integer m = A.Width(), n = A.Height();
    // Matrix<> U(m), V(n);
    Vector<> S(min(n,m));
    Array<double> work(n*m+100);
    ngbla::integer info;
    char jobu = 'A', jobv = 'A';
    ngbla::integer lda = A.Dist(), ldu = U.Dist(), ldv = V.Dist();
    ngbla::integer lwork = work.Size();

    dgesvd_ ( &jobu, &jobv, &m, &n, A.Data(), &lda,
              S.Data(),
              U.Data(), &ldu, V.Data(), &ldv,
              work.Data(), &lwork, 
              &info);
    // cout << "info = " << info << endl;
    // if (n <= 100)
    // cout << "S = " << S << endl;
    A.Diag(0) = S;
  }

  void LapackSVD (SliceMatrix<double,ColMajor> A,
                  SliceMatrix<double, ColMajor> U,
                  SliceMatrix<double, ColMajor> V)
  {
    return LapackSVD(Trans(A), V, U);
  }

  

  template <ORDERING OH, ORDERING OM>  
  void ApplyBandHouseholderReflections (int bs, SliceMatrix<double,OH> H, SliceMatrix<double,OM> M)
  {
    int n = H.Width();
    for (int i = 0; i+bs < H.Width(); i++)
      H.Col(i).Range(i+bs, H.Width()) = 0;

    for (int i = 0; i < H.Width(); i += bs)
      {
        int next = min(i+bs, n);
        MultiHouseholderReflection Hv(Trans(H.Cols(i,next).Rows(i,min(i+2*bs,n))));
        Hv.Mult(M.Rows(i,min(i+2*bs,n)));
      }
  }

  template <ORDERING ORD>
  void ReduceBiDiagonal (SliceMatrix<double,ORD> A, 
                         SliceMatrix<double, ColMajor> U1, SliceMatrix<double, ColMajor> V1)
  {
    static Timer ttrig("householder-triangular"); RegionTimer reg(ttrig);
    static Timer tuv("householder-triangular, transform uv"); 
  
    size_t m = A.Height();
    size_t n = A.Width();
    size_t minnm = min(n,m);

    U1 = Identity(m);
    V1 = Identity(n);

    /*
    // no blocking at all
    Vector<> v2(max(n,m));
    for (int i = 0; i < n; i++)
    {
    auto Arest = A.Rows(i,m).Cols(i,n);
    auto vrest = v2.Range(i,m);
    CalcHouseholderVector(Arest.Col(0), vrest);
    HouseholderReflection (vrest).Mult(Arest);
    HouseholderReflection (vrest).Mult(Trans(U1).Rows(i,m));
        
    if (i+1 < n)
    {
    auto Arest = A.Rows(i,m).Cols(i+1,n);
    auto vrest = v2.Range(i+1,n);
    CalcHouseholderVector(Arest.Row(0), vrest);
    HouseholderReflection (vrest).Mult(Trans(Arest));
    HouseholderReflection (vrest).Mult(Trans(V1).Rows(i+1,n));
    }
    }
    */

    /*
    // no blocking of A, but blocking of U and V
    constexpr size_t bs = 16;
    Matrix<> mvc(bs, m);
    Matrix<> mvr(bs, n);
    for (int i1 = 0; i1 < min(n,m); i1 += bs)
    {
    int bs1 = min( min(n-i1, m-i1), bs);
    int bs2 = 0;
    mvc = 0.0;
    mvr = 0.0;
      
    for (int i2 = 0; i2 < bs1; i2++)
    {
    int i = i1+i2;
    auto Arest = A.Rows(i,m).Cols(i,n);
    CalcHouseholderVector(Arest.Col(0), mvc.Row(i2).Range(i,m));
    HouseholderReflection (mvc.Row(i2).Range(i,m)).Mult(Arest);
          
    if (i+1 < n)
    {
    bs2 = i2+1;
    auto Arest = A.Rows(i,m).Cols(i+1,n);
    CalcHouseholderVector (Arest.Row(0), mvr.Row(i2).Range(i+1,n));
    HouseholderReflection (mvr.Row(i2).Range(i+1,n)).Mult(Trans(Arest));
    }
    }

    MultiHouseholderReflection (mvc.Rows(0,bs1).Cols(i1,m)).Mult (Trans(U1).Rows(i1,m));
    if (bs2 > 0)
    MultiHouseholderReflection (mvr.Rows(0,bs2).Cols(i1+1,n)).Mult (Trans(V1).Rows(i1+1,n));
    }
    // cout << "A bidiag = " << endl << Truncate (A) << endl;
    */

    VectorMem<100> hv(max(n,m));
    for (size_t i = 0; i < minnm; i++)
      {
        auto Acol = A.Cols(i,i+1).Rows(i,m);  // a matrix
        double signed_norm = CalcHouseholderVectorInPlace (Acol.Col(0));
        hv.Range(i,m) = Acol.Col(0);
        Acol(0,0) = signed_norm;
        
        // MultiHouseholderReflection H(Trans(Acol));   
        HouseholderReflection H(hv.Range(i,m));     // more efficient
        H.Mult (A.Rows(i,m).Cols(i+1,n));
        
        if (i+1 < n)
          {
            auto Arow = A.Rows(i,i+1).Cols(i+1,n);  // a matrix
            double signed_norm = CalcHouseholderVectorInPlace (Arow.Row(0));
            hv.Range(i+1,n) = Arow.Row(0);            
            Arow(0,0) = signed_norm;
            
            // MultiHouseholderReflection H(Arow);  
            HouseholderReflection H(hv.Range(i+1,n));   // more efficient
            H.Mult (Trans (A.Rows(i+1,m).Cols(i+1,n)));
          }
      }

    RegionTimer ruv(tuv);    

    ApplyHouseholderReflections (A, Trans(U1));
    if (n > 1)
      ApplyHouseholderReflections (Trans(A.Cols(1,n).Rows(min(n-1,m))), Trans(V1).Rows(1,n));

    for (size_t i = 0; i < minnm; i++)
      A.Col(i).Range(i+1, m) = 0.0;
    for (size_t i = 1; i < minnm; i++)
      A.Col(i).Range(0, i-1) = 0.0;
  }



  template <ORDERING OA, ORDERING ORD>
  void ReduceBiDiagonalBlocked (SliceMatrix<double, OA> A, 
                                SliceMatrix<double, ORD> U1, SliceMatrix<double, ORD> V1)
  {
    /*
      Block-factorization. 

      A parallel algorithm for reducing symmetric banded matrices to tridiagonal form
      Bruno Lang - SIAM Journal on Scientific Computing, 1993 - SIAM
    */
    
    // cout << "ReducedBiDiagnoalBlocked" << endl;
    static Timer ttrig("householder-block-triangular");
    static Timer tUV("householder-block-triangular, transform U1,V1");    
    static Timer tbulgechasing("bulge chasing");
    ttrig.Start();
  
    size_t m = A.Height();
    size_t n = A.Width();
    Matrix A_orig = A; // for testint only
    constexpr size_t bs = 24;  // orig: 32

    
    for (size_t i1 = 0; i1 < min(n,m); i1 += bs)
      {
        size_t bs1 = min( min(n-i1, m-i1), bs);
        for (size_t i2 = 0; i2 < bs1; i2++)
          {
            size_t i = i1+i2;

            auto Apanel = A.Rows(i,m).Cols(i,i1+bs1);            
            double signed_norm = CalcHouseholderVectorInPlace (Apanel.Col(0));
            Apanel(0,0) = signed_norm;
            
            MultiHouseholderReflection H(Trans(Apanel.Cols(0,1)));
            H.Mult (Apanel.Cols(1, Apanel.Width()));
          }

        if (i1+bs1 < n)
          {
            MultiHouseholderReflection H(Trans(A.Cols(i1, i1+bs1).Rows(i1, m)));
            H.Mult (A.Rows(i1,m).Cols(i1+bs1,n));
          }
        
        for (int i2 = 0; i2 < bs1; i2++)
          {
            size_t i = i1 + i2;
            if (i+bs1 < n)
              {
                auto Apanel = A.Rows(i,i1+bs1).Cols(i+bs1, n);
                double signed_norm = CalcHouseholderVectorInPlace (Apanel.Row(0));
                Apanel(0,0) = signed_norm;
                
                MultiHouseholderReflection H(Apanel.Rows(0,1));
                H.Mult (Trans(Apanel.Rows(1, Apanel.Height())));
              }
          }
        
        size_t bs2 = min(bs1, n-i1-bs1);
        if (bs2 > 0)
          if (i1+bs1 < m)
            {
              MultiHouseholderReflection H(A.Rows(i1, i1+bs2).Cols(i1+bs1, n));
              H.Mult (Trans(A.Rows(i1+bs1,m).Cols(i1+bs1,n)));
            }
      }
    
    // ttrig.Stop();
    // tUV.Start();

    U1 = Identity(m);
    V1 = Identity(n);

    ApplyHouseholderReflections (A, Trans(U1));
    if (n > bs)
      ApplyHouseholderReflections (Trans(A.Cols(bs,n).Rows(min(n-bs,m))), Trans(V1).Rows(bs,n));      

    size_t minnm = min(n,m);

    for (size_t i = 0; i < minnm; i++)
      A.Col(i).Range(i+1, m) = 0.0;
    for (size_t i = bs; i < minnm; i++)
      A.Col(i).Range(0, i-bs) = 0.0;
    
    // tUV.Stop();

    
    
    /*
      Matrix<> Abandmem(minnm, minnm);
      Abandmem = 0.0;
      // simulate a band-matrix by reducing dist:
      SliceMatrix<> Aband(minnm, minnm, min(3*bs, minnm), Abandmem.Data());
      for (int i = 0; i < min(n, bs+1); i++)
      Aband.Diag(i) = A.Rows(minnm).Cols(minnm).Diag(i);
    */
    // auto Aband = A.Rows(minnm).Cols(minnm);
    
    Matrix Aband = A.Rows(minnm).Cols(minnm);
    
    
    /*
      if (n <= 20)
      {
      cout << "householder, blocked = " << endl << Truncate(A) << endl;
      cout << "U^T A V = " << endl << Truncate( Matrix (Matrix(Trans(U1)*A_orig) * V1)) << endl;
      }
    */

    
    // now do the bulge chasing
    tbulgechasing.Start();
    /*
      Vector<> tmp(bs);
      for (int i = 0; i < n-1; i++)
      {
      // first vector
      IntRange cols(i+1, min(n, i+bs+1));
      auto x = A.Row(i).Range(cols);
      CalcHouseholderVector (x, tmp.Range(x.Size()));

      IntRange rows(i, min(m,i+bs+1));
      HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(A.Rows(rows).Cols(cols)));
      HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(V1.Cols(cols)));      
      
      for (int i1 = i+1; i1 < n; i1 += bs)
      {
      {
      IntRange rows(i1, min(n, i1+bs));
      auto x = A.Col(i1).Range(rows);
      CalcHouseholderVector (x, tmp.Range(x.Size()));
            
      IntRange cols(i1, min(n, i1+2*bs));
      HouseholderReflection (tmp.Range(x.Size())).Mult (A.Rows(rows).Cols(cols));
      HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(U1.Cols(rows)));                  
      }
          
      if (i1+bs < n)
      {
      IntRange cols(i1+bs, min(n, i1+2*bs));
      auto x = A.Row(i1).Range(cols);
      CalcHouseholderVector (x, tmp.Range(x.Size()));
              
      IntRange rows(i1, min(n, i1+2*bs));
      HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(A.Rows(rows).Cols(cols)));
      HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(V1.Cols(cols)));                    
      }
      }
      }
    */



    /*
      TODO: store the bulge-chasing reflections in triangular matrix
      write banded - Multihouseholder
    */
    

    /*
      constexpr size_t bcbs = 16;   // bulge chasing block size
      Vector<> tmp(bs);
      int maxblocks = n/bs+1;
      Matrix<> reflect_bcU(bcbs * maxblocks, bs+bcbs-1);
      Matrix<> reflect_bcV(bcbs * maxblocks, bs+bcbs-1);
      Array<IntRange> ranges;
                                     
      for (int bi = 0; bi < n-1; bi += bcbs)
      {
      int bsi = min(bcbs, n-1-bi);
      reflect_bcU = 0.0;
      reflect_bcV = 0.0;
      
      for (int ii = 0; ii < bsi; ii++)
      {
      int i = bi + ii;
      // first vector
      IntRange cols(i+1, min(n, i+bs+1));
      auto x = Aband.Row(i).Range(cols);
      CalcHouseholderVector (x, tmp.Range(x.Size()));
      reflect_bcV.Row(0*bcbs+ii).Range(ii, ii+cols.Size()) = tmp.Range(x.Size());
                
      IntRange rows(i, min(m,i+bs+1));
      HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(Aband.Rows(rows).Cols(cols)));
      // HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(V1.Cols(cols)));      
      int nr = 0;
      for (int i1 = i+1; i1 < n; i1 += bs, nr++)
      {
      {
      IntRange rows(i1, min(n, i1+bs));
      auto x = Aband.Col(i1).Range(rows);
      CalcHouseholderVector (x, tmp.Range(x.Size()));
      reflect_bcU.Row(nr*bcbs+ii).Range(ii, ii+rows.Size()) = tmp.Range(x.Size());
                
      IntRange cols(i1, min(n, i1+2*bs));
      HouseholderReflection (tmp.Range(x.Size())).Mult (Aband.Rows(rows).Cols(cols));
      // HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(U1.Cols(rows)));
      }

      if (i1+bs < n)
      {
      IntRange cols(i1+bs, min(n, i1+2*bs));
      auto x = Aband.Row(i1).Range(cols);
      CalcHouseholderVector (x, tmp.Range(x.Size()));
      reflect_bcV.Row((nr+1)*bcbs+ii).Range(ii, ii+cols.Size()) = tmp.Range(x.Size());
                  
      IntRange rows(i1, min(n, i1+2*bs));
      HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(Aband.Rows(rows).Cols(cols)));
      // HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(V1.Cols(cols)));                    
      }
      }
      }

      ranges.SetSize(0);
      for (int i1 = bi+1; i1 < n; i1 += bs)
      ranges.Append (IntRange(i1, min(n,i1+bs+bcbs-1)));
      for (int i = ranges.Size()-1; i >= 0; i--)
      {
      int bcbsi = min(bcbs, ranges[i].Size());
      MultiHouseholderReflection H(reflect_bcU.Rows(i*bcbs, i*bcbs+bcbsi).Cols(ranges[i].Size()));
      H.Mult (Trans(U1.Cols(ranges[i])));
      }
      for (int i = ranges.Size()-1; i >= 0; i--)
      {
      int bcbsi = min(bcbs, ranges[i].Size());          
      MultiHouseholderReflection (reflect_bcV.Rows(i*bcbs, i*bcbs+bcbsi).Cols(ranges[i].Size())).Mult (Trans(V1.Cols(ranges[i])));
      }
      }
    */



    // constexpr size_t bcbs = 16;   // bulge chasing block size
    Vector<> tmp(bs);
    // int maxblocks = n/bs+1;
    // Matrix<> reflect_bcU(bcbs * maxblocks, bs+bcbs-1);
    // Matrix<> reflect_bcV(bcbs * maxblocks, bs+bcbs-1);
    // Array<IntRange> ranges;

    Array<tuple<int,int,int>> refu;
    Array<tuple<int,int,int>> refv;

    static Timer tb1 ("bulge chasing - transform band");
    tb1.Start();

    for (size_t i = 0; i+1 < n; i++)
      {
        // first vector
        IntRange cols(i+1, min(n, i+bs+1));
        auto x = Aband.Row(i).Range(cols);
        
        CalcHouseholderVector (x, tmp.Range(x.Size()));
        refv.Append ( tuple(i,cols.First(),cols.Size() ));
                  
        IntRange rows(i, min(m,i+bs+1));
        HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(Aband.Rows(rows).Cols(cols)));
        Aband.Row(i).Range(cols.First()+1, cols.Next()) = tmp.Range(1,x.Size());
            
        int nr = 0;
        for (int i1 = i+1; i1 < n; i1 += bs, nr++)
          {
            {
              IntRange rows(i1, min(n, i1+bs));
              auto x = Aband.Col(i1).Range(rows);
              CalcHouseholderVector (x, tmp.Range(x.Size()));

              Aband.Col(i).Range(rows) = tmp.Range(x.Size());
              refu.Append ( tuple(i1,i,rows.Size() ));
              IntRange cols(i1, min(n, i1+2*bs));
              HouseholderReflection (tmp.Range(x.Size())).Mult (Aband.Rows(rows).Cols(cols));
            }
            
            if (i1+bs < n)
              {
                IntRange cols(i1+bs, min(n, i1+2*bs));
                auto x = Aband.Row(i1).Range(cols);
                CalcHouseholderVector (x, tmp.Range(x.Size()));
                
                Aband.Row(i).Range(cols) = tmp.Range(x.Size());
                refv.Append ( tuple(i,cols.First(),cols.Size() ));
                
                IntRange rows(i1, min(n, i1+2*bs));
                HouseholderReflection (tmp.Range(x.Size())).Mult (Trans(Aband.Rows(rows).Cols(cols)));
              }
          }
      }
    
    tb1.Stop();
    static Timer tb2 ("bulge chasing - mult U");
    static Timer tb3 ("bulge chasing - mult V");
    tb2.Start();
    /*
      for (auto [i,j,size] : refu)
      {
      tmp.Range(size) = Aband.Col(j).Range(i,i+size);
      HouseholderReflection H(tmp.Range(size));
      H.Mult (Trans(U1.Cols(i,i+size)));
      }
    */
    
    // for (auto [i,j,size] : refu)
    // cout << "i,j = " << i << "," << j << ", size = " << size << endl;
    
    int startrow = 1;
    while (startrow+bs < n) startrow += bs;
    while (startrow > 0)
      {
        IntRange rest(startrow, n);
        auto sub = Aband.Rows(rest).Cols(0, n-startrow);
        ApplyBandHouseholderReflections (bs, sub, Trans(U1.Cols(rest)));
        startrow -= bs;
      }
    
    tb2.Stop();
    tb3.Start();

    /*
      for (auto [i,j,size] : refv)
      {
      tmp.Range(size) = Aband.Row(i).Range(j,j+size);
      tmp(0) = 1;
      HouseholderReflection H(tmp.Range(size));
      H.Mult (Trans(V1.Cols(j,j+size)));
      }
    */
    startrow = 1;
    while (startrow+bs < n) startrow += bs;
    while (startrow > 0)
      {
        IntRange rest(startrow, n);
        auto sub = Aband.Cols(rest).Rows(0, n-startrow);
        ApplyBandHouseholderReflections (bs, Trans(sub), Trans(V1.Cols(rest)));
        startrow -= bs;
      }

    
    tb3.Stop();
    tbulgechasing.Stop();
    
    A.Rows(minnm).Cols(minnm) = Aband;
    
    // if (n < 20)
    // cout << "after bulge-chasing: " << endl << Truncate(A) << endl;
  }





  /* ***************** Bi-diagonal SVD by bisection ************ */

  /*
    Othogonal eigenvectors and relative gaps
    I.S. Dhillon, B.N. Parlett
    SIAM J. Matrix Anal. Appl, 25(3), pp 858-899
  */

  /*
    int NegCount (SliceMatrix<> B, double mu)
    {
    double t = -mu;
    int negcnt = 0;
    int n = B.Width();
    for (int k = 0; k < n-1; k++)
    {
    double d = sqr(B(k,k)) + t;
    if (d < 0) negcnt++;
    t = t * sqr (B(k,k+1))/d - mu;
    }
    double d = sqr(B(n-1, n-1)) + t;
    if (d < 0) negcnt++;
    return negcnt;
    }
  */

  int NegCount (SliceMatrix<> B, double mu)
  {
    double t = -mu;
    int negcnt = 0;
    size_t n = B.Width();
    auto diag = B.Diag(0);
    auto super_diag = B.Diag(1);
    // cout << "Bptr = " << B.Data() << "diag.pr = " << diag.Data() << ", supe_diag.data = " << super_diag.Data() << endl;
    for (size_t k = 0; k < n-1; k++)
      {
        double d = sqr(diag(k)) + t;
        if (d < 0) negcnt++;
        // t = t * sqr (B(k,k+1))/d - mu;
        t = t * sqr (super_diag(k))/d - mu;
      }
    double d = sqr(diag(n-1)) + t;
    if (d < 0) negcnt++;
    return negcnt;
  }


  SIMD<int64_t> NegCount (SliceMatrix<> B, SIMD<double> mu)
  {
    SIMD<double> t = -mu;
    SIMD<int64_t> negcnt = 0;
    size_t n = B.Width();
    auto diag = B.Diag(0);
    auto super_diag = B.Diag(1);

    for (size_t k = 0; k < n-1; k++)
      {
        SIMD<double> d = sqr(diag(k)) + t;
        negcnt += If (d < SIMD<double>(0.0), SIMD<int64_t>(1), SIMD<int64_t>(0));
        t = t * sqr (super_diag(k))/d - mu;
      }
    SIMD<double> d = sqr(diag(n-1)) + t;
    negcnt += If (d < SIMD<double>(0.0), SIMD<int64_t>(1), SIMD<int64_t>(0));  
    return negcnt;
  }



  // stationary quotient-difference with shift
  // L D L^t - mu I  =  Lp Dp Lp^t
  void stqds (FlatVector<> D, FlatVector<> L, double mu,
              FlatVector<> Dp, FlatVector<> Lp)
  {
    size_t n = D.Size();
    Dp(0) = D(0)-mu;
    for (size_t i = 0; i < n-1; i++)
      {
        Lp(i) = D(i)*L(i) / Dp(i);
        Dp(i+1) = D(i)*L(i)*L(i) + D(i+1)-Lp(i)*D(i)*L(i) - mu;
      }
  }

  // differential form of stationary quotient-difference with shift
  // L D L^t - mu I  =  Lp Dp Lp^t
  void dstqds (FlatVector<> D, FlatVector<> L, double mu,
               FlatVector<> Dp, FlatVector<> Lp, FlatVector<> gamma)
  {
    size_t n = D.Size();
    double s = -mu;
    for (size_t i = 0; i < n-1; i++)
      {
        gamma(i) = s;
        Dp(i) = s+D(i);
        Lp(i) = D(i)*L(i)/Dp(i);
        s = Lp(i)*L(i)*s - mu;
      }
    Dp(n-1) = s + D(n-1);
    gamma(n-1) = s;
  }



  // differential norm of the progressive qd transform
  // L D L^t - mu I  =  Um Dm Um^t
  void dqds (FlatVector<> D, FlatVector<> L, double mu,
             FlatVector<> Dm, FlatVector<> Um, FlatVector<> gamma)
  {
    size_t n = D.Size();
    double p = D(n-1) - mu;
    for (int i = n-2; i >= 0; i--)
      {
        gamma(i+1) += p + mu;
        Dm(i+1) = D(i)*L(i)*L(i) + p;
        double t = D(i)/Dm(i+1);
        Um(i) = L(i)*t;
        p = p*t-mu;
      }
    Dm(0) = p;
    gamma(0) += p+mu;
  }


  void CalcSVDBiDiagonal1 (SliceMatrix<> B, 
                           FlatMatrix<double,ColMajor> U, FlatMatrix<double,ColMajor> V)
  {
    static Timer tb("tbisect");
    static Timer tbi("tsolvebi");
    static Timer tub("tsmallUV");
    double upperbound = 0;
  
    size_t n = B.Width();
    Vector<> singvals(n);

    // int m = B.Height();
    for (int k = 0; k < B.Width(); k++)
      upperbound = max2(upperbound, fabs(B(k,k)));
    for (int k = 0; k < B.Width()-1; k++)
      upperbound = max2(upperbound, fabs(B(k+1,k)));
    upperbound *= 3; // ????
    // cout << "upperbound = " << upperbound << endl;

    tb.Start();
    /*
      for (int i = 0; i < n; i++)
      {
      double l = 0, u = upperbound;
      for (int bisectit = 0; bisectit < 50; bisectit++)  // 1/2**50 < floating pnt prec
      {
      double mid = l + 0.5 * (u-l);
      int negcnt = NegCount (B, mid*mid);
      u = (negcnt > i) ? mid : u;
      l = (negcnt > i) ? l : mid;
      }
      singvals[i] = l;      
      }
    */

    for (size_t i = 0; i < n; i += SIMD<double>::Size())
      {
        SIMD<double> l = 0.0, u = upperbound;
        SIMD<int64_t> evnr = SIMD<int64_t>(n-i-1)-SIMD<int64_t>::FirstInt(); // start with the largest
        for (int bisectit = 0; bisectit < 50; bisectit++)  // 1/2^50 < floating pnt prec
          {
            SIMD<double> mid = 0.5 * (u+l); 
            SIMD<int64_t> negcnt = NegCount (B, mid*mid);
            u = If(negcnt > evnr, mid, u);
            l = If(negcnt > evnr, l, mid);
          }
        SIMD<mask64> mask(n-i);
        l.Store (&singvals[i], mask);
      }


  
    tb.Stop();
    tb.AddFlops (n * 50 * n);
    // B^T B = L D L^T
    Vector D(n), L(n-1);
    for (int i = 0; i < n; i++)
      D(i) = sqr(B(i,i));
    for (int i = 0; i < n-1; i++)
      L(i) = B(i,i+1) / B(i,i);

    // B = U Sigma V
    // Matrix U(m,m);
    // Matrix V(n,n);
    U = 0.0;

    Vector gamma(n);
    Vector Lp(n-1), Dp(n);
    Vector Um(n-1), Dm(n);
    Vector z(n);   // the eigenvector
      
    for (int lamnr = 0; lamnr < n; lamnr++)
      {
        // cout << "compute sing " << lamnr << ": " << singvals[lamnr]
        // << " ^ 2 = " << sqr(singvals[lamnr]) << endl;
        tbi.Start();
        double lam = sqr(singvals[lamnr]);

        /*
          stqds (D, L, lam, Dp, Lp);
          cout << "Dp = " << Dp << endl;
          cout << "Lp = " << Lp << endl;
        */
        dstqds (D, L, lam, Dp, Lp, gamma);
        // cout << "Dp = " << Dp << endl;
        // cout << "Lp = " << Lp << endl;

        /*
          Matrix Dmat(n), Lmat(n);
          Dmat = 0; Lmat = 0;
          for (int i = 0; i < n; i++)
          {
          Lmat(i,i) = 1;
          Dmat(i,i) = Dp(i);
          }
          for (int i = 0; i < n-1; i++)
          Lmat(i+1,i) = Lp(i);
      
          Matrix LDLt = Lmat * Dmat * Trans(Lmat);
          Vector lams(n);
          Matrix evecs(n,n);
          CalcEigenSystem (LDLt, lams, evecs);
          cout << "check from LpDpLp^t" << endl << lams << endl;
        */
      
        dqds (D, L, lam, Dm, Um, gamma);
        // cout << "gamma = " << gamma << endl;

        /*
        // simple version for gamma:
        Vector gamma1(n);
        for (int i = 0; i < n-1; i++)
        gamma1(i) = Dp(i) - sqr(D(i)*L(i)) / Dm(i+1);
        gamma1(n-1) = Dp(n-1);
        cout << "gamma =?= " << gamma1 << endl;
        */

      
        int r = 0;
        for (int i = 0; i < n; i++)
          if (fabs(gamma(i)) < gamma(r))
            r = i;

        z(r) = 1;
        for (int i = r-1; i >= 0; i--)
          {
            if (Dp(i) != 0)
              z(i) = -Lp(i)*z(i+1);
            else
              z(i) = -(D(i+1)*L(i+1)/(D(i)*L(i))) * z(i+2);
          }
        for (int i = r; i < n-1; i++)
          {
            if (Dm(i+1) != 0)
              z(i+1) = -Um(i)*z(i);
            else
              z(i+1) = -(D(i-1)*L(i-1)/(D(i)*L(i))) * z(i-1);
          }
        // cout << "z = " << z << endl;
        z /= L2Norm(z);
        tbi.Stop();
        tub.Start();
        V.Col(lamnr) = z;
        // U.Col(lamnr) = B * z;
        U.Col(lamnr) = 0.0;
        for (int i = 0; i < n; i++)
          U.Col(lamnr)(i) = B(i,i) * z(i);
        for (int i = 0; i < n-1; i++)
          U.Col(lamnr)(i) += B(i,i+1) * z(i+1);
        U.Col(lamnr) /= singvals[lamnr];
        tub.Stop();
        // cout << "|Ui| = " << L2Norm(U.Col(lamnr)) << endl;
      }

    // cout << "V^T * V = " << Trans(V)*V << endl;
    // cout << "U^T * U = " << Trans(U)*U << endl;
  }


  /* ***************** Bi-diagonal SVD by divide-and-conquer ************ */

  /*
    matrix is (n+1)*n
    alpha = m.Diag(0)
    beta = m.Diag(-1)
  */


  void CalcRecLapack (FlatVector<> alpha, FlatVector<> beta, FlatVector<> sigma,
                      SliceMatrix<double,ColMajor> U, SliceMatrix<double,ColMajor> V)
  {
    size_t n = alpha.Size();
    if (n == 0)
      {
        U(0,0) = 1;
        return;
      }
    
    if (U.Height() != n+1 || U.Width() != n+1) throw Exception ("U has wrong size");
    if (V.Height() != n || V.Width() != n) throw Exception ("V has wrong size");
    Matrix M(n+1, n);
    M = 0;
    M.Diag(0) = alpha;
    M.Diag(-1) = beta;
    Matrix Morig = M;
    Matrix<double,ColMajor> hU = U;
    cout << "lapack SVD of matrix " << endl << M << endl;
    LapackSVD (M, V, hU);
    U = Trans(hU);
    sigma = M.Diag(0);

    cout << "check lapack:" << endl;
    cout << "alpha, beta = " << alpha << ", " << beta << endl;
    cout << Truncate(Trans(U) * Morig * V) << endl;
    cout << "U = " << endl << U << endl;
    cout << "V = " << endl << V << endl;
  }





  class GivensRotation
  {
    double c, s;
    int i1, i2;
  public:
    GivensRotation (double ac, double as, int ai1, int ai2)
      : c(ac), s(as), i1(ai1), i2(ai2) { ; }
  
    void Apply (SliceVector<> v, bool trans = false)
    {
      double v1 = v(i1), v2 = v(i2);
      if (trans)
        {
          v(i1) = c * v1 + s * v2;
          v(i2) = -s * v1 + c * v2;
        }
      else
        {
          v(i1) = c * v1 - s * v2;
          v(i2) = s * v1 + c * v2;
        }
    }

    template <ORDERING OM>
    void Apply (SliceMatrix<double,OM> mat, bool trans = false)
    {
      for (int i = 0; i < mat.Width(); i++)
        Apply (mat.Col(i), trans);
    }
  };




  
  /*
    find solution within (lower, upper) of
    1 + \sum_i  vs(i)**2 / (lams(i)^2 - lam^2) = 0
  */
  // returns t-lower and upper-t
  // the smaller one is more accurate !
  tuple<double,double> SolveSecular (FlatVector<> lams,
                                     FlatVector<> vs, double lower, double upper)
  {
    /*
      cout << "solve sec1ular, lams = " << endl << lams << endl
      << "vs = " << vs << endl
      << "l/u = " << lower << "/" << upper << endl;
    */
    // static Timer timer("SolveSecular"); RegionTimer r(timer);
    // static Timer tloop("SolveSecular loop ");
    
    /*
      Rank-One Modification of the Symmetric Eigenproblem
      Bunch, Nielsen, Sorensen
      Numer Math 31, 31-48, 1978

      with modifications by Gu + Eisenstat for SVD
    */

    // cout << "solvesecular called, lower = " << lower << ", upper = " << upper << endl;
    ArrayMem<double,100> vpsi, vphi;
    ArrayMem<double,100> lampsi_shift, lamphi_shift;
    double lower2 = sqr(lower);
    double upper2 = sqr(upper);
    double mid = 0.5 * (upper2+lower2);
  
    for (int i = 0; i < lams.Size(); i++)
      if (vs[i] != 0)        // only terms with non-zero vs count
        {
          if (sqr(lams[i]) < mid)
            {
              vpsi.Append(vs[i]);
              lampsi_shift.Append (lams[i]-lower);
            }
          else
            {
              vphi.Append(vs[i]);
              lamphi_shift.Append (lams[i]-upper);
            }
        }

    /*
      cout << "vpsi = " << vpsi << endl;
      cout << "lampsi_shift = " << lampsi_shift << endl;
      cout << "vphi = " << vphi << endl;
      cout << "vphi_shift = " << lamphi_shift << endl;
    */
    // RegionTimer rloop(tloop);

    double delta = upper-lower;
    double t = delta/2;

    // double delta2 = upper2-lower2;
    double Delta = delta-t;

    int maxits = 40;
    for (int i = 0; i < maxits; i++)
      {
        // cout << "i = " << i << ", t = " << t << ", Delta = " << Delta << endl;
      
        double t2 = t*t+2*lower*t;  // sqr(t+lower)-sqr(lower);
        double Delta2 = Delta*(2*upper-Delta);    // sqr(Delta)-sqr(upper-Delta)   
        // double Delta2 = delta2-t2;
        // cout << "Delta2 = " << Delta2 << " =?= " << delta2-t2 << endl;
      
        // Delta = delta - t;
        // cout << "Delta = " << Delta << " =?= " << delta-t << " i = " << i << ", t = " << t << endl;
        if (isnan(t2) || isinf(t2))
          {
            cout << "t = " << t2 << endl;
            cout << "vpsi = " << vpsi << endl;
            cout << "vphi = " << vphi << endl;
            cout << "lampsi_shift = " << lampsi_shift << endl;
            cout << "lamphi_shift = " << lamphi_shift << endl;
            throw Exception ("secular failed");
          }

      
        double phi = 0, psi=0, phiprime=0, psiprime=0;
        double calt = 1;

        for (int j = 0; j < vpsi.Size(); j++)
          {
            // double num = ( (lampsi[j]-lower-t) * (lampsi[j]+lower+t) );
            double num = ( (lampsi_shift[j]-t) * (lampsi_shift[j]+2*lower+t) );
            psi += sqr(vpsi[j]) / num;
            psiprime += sqr(vpsi[j]) / sqr( num );
          }
      
        for (int j = 0; j < vphi.Size(); j++)
          {
            // double num = (lamphi[j]-lower-t)*(lamphi[j]+lower+t);
            double num = (lamphi_shift[j]+Delta)*(lamphi_shift[j]+2*upper-Delta);
            phi += sqr(vphi[j]) / num;
            phiprime += sqr(vphi[j]) / sqr(num);
            calt += sqr(vphi[j])  * lamphi_shift[j]*(lamphi_shift[j]+2*upper) / sqr(num);
          }
      
        double w = 1 + phi + psi;
      
        // cout << "i = " << i << ", psi = " << psi << ", phi = " << phi << ", w = " << w << endl;

        if (i == 0 && w > 0)  // need a legal starting value with w < 0
          {
            t /= 2;
            Delta = delta-t;
            i--;
            continue;
          }

        
        if (vphi.Size() > 0)
          {
            // double c1 = 1 + phi - Delta2 * phiprime;
            double c = calt;
            double a = (Delta2*(1+phi)+sqr(psi)/psiprime)/c + psi/psiprime;
            double b = (Delta2 * w * psi) / (psiprime * c);
          
            double disc = a*a-4*b;
            if (disc < 0) disc = 0;
            double incr = 2*b/(a+sqrt(disc));

            if (Delta2-incr < 1e-6 * Delta2) incr *= 0.999999;
            if (t2+incr < 1e-6 * Delta2) incr *= 0.999999;

            t2 += incr;
            Delta2 -= incr;

            t = t2 / (sqrt(t2+lower2)+lower);
          
            double discu = upper2-Delta2;
            if (discu < 0) discu = 0;
            Delta = Delta2 / (sqrt(discu)+upper);

            /*
            if (fabs (incr) < 1e-15 * min(t2, Delta2))
              {
                // cout << "break due to small inc" << endl;
                // chosen, but not necessary ..
                if (i < maxits-3) i = maxits-3;
              }
            */
          }
        else
          {
            double incr = (1+psi)/psiprime*psi;
            t2 += incr;
            t = t2 / (sqrt(t2+lower2)+lower);
            Delta = delta-t;
          }
      
        if (fabs (w) < 1e-14 * (1.0+fabs(phi)+fabs(psi)))
          if (i < maxits-2) i = maxits-2;
        
        // if (i > 20 && i < 20)
        // cout << "many its, i = " << i << endl;
        
        if (t < Delta)
          Delta = delta - t;
        else
          t = delta - Delta;
      }

    return { t, Delta };
  }







  
  







  
  /*
    Matrix is (n+1) * n
    alpha diagonal, beta .. sub-diagonal
    sigma(n) singular values
    U   (n+1) * (n+1)  .... first n cols are singular vectors
    V   n * n


    first n columns form the reduced SVD
    if alpha(0) == 0, then first row of first cols is 0
    i.e. U(1:n+1, 0:n) is U-factor of reduced matrix
  */
  
  void CalcSVDBiDiagonalRec (FlatVector<> alpha, FlatVector<> beta, FlatVector<> sigma,
                             SliceMatrix<double,ColMajor> U, SliceMatrix<double,ColMajor> V)
  {
    /*
    static Timer t("CalcSVDBiDiagonalRec"); RegionTimer r(t);
    static Timer t1("CalcSVDBiDiagonalRec 1"); 
    static Timer t2("CalcSVDBiDiagonalRec 2");
    static Timer t2b("CalcSVDBiDiagonalRec 2b");
    static Timer tWQ("CalcSVDBiDiagonalRec Build WQ");
    static Timer tzD("CalcSVDBiDiagonalRec Build zD");
    static Timer tdef("deflation");
    static Timer tMUV("MU,MV");
    static Timer tsort("sort");
    */

    
    /*
      CalcRecLapack (alpha, beta, sigma, U, V);
      cout << "lapack sigma = " << sigma << endl;
      return;
    */

    // cout << "CalcSVDBiDiagonalRec called, n = " << alpha.Size() << endl;
    // cout << "CalcSVDBiDiagonalRec, alpha " << endl << alpha << endl << "beta = " << endl << beta << endl;

    size_t n = alpha.Size();
    if (n == 0)
      {
        U(0,0) = 1;
        return;
      }

    // t1.Start();

    double normmat = L2Norm(alpha) + L2Norm(beta) + 1e-40;

    /*
    Matrix<> B(n+1, n);  // for testing
    B = 0;
    B.Diag(0) = alpha;
    B.Diag(-1) = beta;
    */
    
    size_t k = n/2;
    double alphak = alpha(k);
    double betak = beta(k);

    // borrow memory ...
    auto U1 = U.Rows(k+1).Cols(k+1);
    auto U2 = U.Rows(k+1,n+1).Cols(k+1,n+1);
    auto W1 = V.Rows(k).Cols(k);
    auto W2 = V.Rows(k,n-1).Cols(k,n-1);
    
    VectorMem<100> D(n);
    D(0) = 0;

    // t1.Stop();
    CalcSVDBiDiagonalRec (alpha.Range(0,k), beta.Range(0,k), D.Range(1,k+1), U1, W1);
    CalcSVDBiDiagonalRec (alpha.Range(k+1,n), beta.Range(k+1,n), D.Range(k+1,n), U2, W2);

    // tWQ.Start();

    /*
    if (isnan(L2Norm(D)) || isinf(L2Norm(D)))
      {
        cout << "D has nans" << endl;
        cout << "alpha = " << alpha << endl;
        cout << "beta = " << beta << endl;
        cout << "D = " << D << endl;
        cout << "U1 = " << U1 << endl;
        cout << "W1 = " << W1 << endl;
        cout << "U2 = " << U2 << endl;
        cout << "W2 = " << W2 << endl;
        throw Exception("D has nans");
      }
    */
    
    /*
      cout << "own recursive: D = " << D << endl;
      cout << "U1,V1 = " << endl << U1 << endl << W1 << endl;
      cout << "U2,V2 = " << endl << U2 << endl << W2 << endl;
    */
    // CalcRecLapack (alpha.Range(0,k), beta.Range(0,k), D.Range(1,k+1), U1, W1);
    // CalcRecLapack (alpha.Range(k+1,n), beta.Range(k+1,n), D.Range(k+1,n), U2, W2);
    /*
    cout << "Lapack recursive: D = " << D << endl;
    cout << "U1,V1 = " << endl << U1 << endl << W1 << endl;
    cout << "U2,V2 = " << endl << U2 << endl << W2 << endl;
    */
    
    
    auto Q1 = U1.Cols(0,k);
    auto q1 = U1.Col(k);
    auto Q2 = U2.Cols(0,n-k-1);
    auto q2 = U2.Col(n-k-1);

    double lam1 = q1(k);
    auto l1 = Q1.Row(k);
    double phi2 = q2(0);
    auto f2 = Q2.Row(0);

    double r0 = hypot(alphak*lam1, betak*phi2);

    double c0, s0;
    if (r0 != 0)
      { c0 = alphak*lam1/r0; s0 = betak*phi2/r0; }
    else
      { c0 = 0; s0 = 1; }

      
    // Matrix<> Q(n+1,n), W(n,n);
    ArrayMem<double,512> memQW((2*n+1)*n);
    FlatMatrix<> Q(n+1, n, memQW.Data());
    FlatMatrix<> W(n, n, memQW.Data()+(n+1)*n);
    
    Q = 0.0;
    auto [Qtop,Qbot] = Q.SplitRows(k+1);
    Qtop.Col(0) = c0*q1;
    Qbot.Col(0) = s0*q2;
    Qtop.Cols(1,1+k) = Q1;
    Qbot.Cols(1+k,n) = Q2;
    W = 0.0;
    W.Rows(0,k).Cols(1,k+1) = W1;
    W(k,0) = 1;
    W.Rows(k+1,n).Cols(k+1,n) = W2;
    VectorMem<100> q(n+1);
    q.Range(0,k+1) = -s0*q1;
    q.Range(k+1,n+1) = c0*q2;
    
    // now, B = Q * M * Trans(W)
    // with M = diag + first col
    /*
      cout << "arrow head before deflation" << endl;
      Matrix M = Matrix(Trans(Q)*B) * W;
      cout << "Q^T B W = " << Truncate(M) << endl;
      cout << "q^T B = " << Truncate( Trans(B)*q ) << endl;
    */

    // tWQ.Stop();
    // tzD.Start();
    
    VectorMem<100> z(n);
    z(0) = r0;
    z.Range(1,k+1) = alphak*l1;
    z.Range(k+1,n) = betak*f2;
    
    // cout << "z = " << endl << z << endl;
    // cout << "D = " << endl << D << endl;
    
    ArrayMem<int,100> index(n);
    for (int i : Range(n)) index[i] = i;
    QuickSortI (FlatArray(D.Size(), D.Data()), index.Range(1,n));

    VectorMem<100> zs(n), Ds(n), Ds2(n);
    for (size_t i = 0; i < n; i++)
      {
        zs(i) = z[index[i]];
        Ds(i) = D[index[i]];
        Ds2(i) = sqr(Ds(i));
      }

    /*
      cout << "before deflation: " << endl;
      cout << "zs = " << endl << zs << endl;
      cout << "Ds = " << endl << Ds << endl;
    */

    // tzD.Stop();
    // tdef.Start();
    
    // deflation  (Gu + Eisenstat, Chap 4.1)
    // if Ds are the same, rotate up
    // BitArray sameDs(n);
    // sameDs.Clear();
    ArrayMem<bool,100> sameDs(n);
    sameDs = false;
    for (int i = n-2; i >= 0; i--)
      {
        if ( (Ds(i+1) - Ds(i)) < 1e-14*normmat)
          {
            sameDs[i+1] = true; // .SetBit(i+1);
            double c = zs(i), s = -zs(i+1);
            double r = hypot (c, s);
            if (r > 1e-20 * normmat)
              {
                c /= r; s /= r;
                GivensRotation G(c, s, index[i], index[i+1]);
                zs(i) = r;
                zs(i+1) = 0;
                G.Apply (SliceMatrix(Trans(Q)));
                if (i > 0)
                  G.Apply (SliceMatrix(Trans(W))); 
              }
            else
              {
                zs(i) = 0;
                zs(i+1) = 0;
              }
          }
      }
    for (size_t i = 1; i < n; i++)
      if (sameDs[i])
        Ds(i) = Ds(i-1);

    ArrayMem<bool,100> deflated(n);
    deflated = false;
    if (fabs(zs(0)) < 1e-14 * normmat) zs(0) = 1e-14 * normmat; 
    for (int i = 1; i < n; i++)
      if (fabs(zs(i)) < 1e-14 * normmat)
        {
          zs(i) = 0;
          deflated[i] = true;
        }

    /*
      cout << "arrowhead, after deflation:" << endl;
      Matrix M2 = Matrix(Trans(Q)*B) * W;
      cout << "Q^T B W = " << Truncate(M2) << endl;
      cout << "q^T B = " << Truncate( Trans(B)*q ) << endl;
    */

    
    VectorMem<100> Ds_next(n);
    double dsn = Ds(n-1)+L2Norm(zs);
    for (int i = n-1; i >= 0; i--)
      {
        Ds_next(i) = dsn;
        if (!deflated[i])
          dsn = Ds(i);
      }

    // tdef.Stop();
    // RegionTimer r2b(t2b);
    
    /*
      cout << "before Secular:" << endl;
      cout << "deflated: " << endl << deflated << endl;
      cout << "zs = " << endl << zs << endl;
      cout << "Ds = " << endl << Ds << endl;
      cout << "Ds_next = " << endl << Ds_next << endl;
    */


    // static Timer tsec("Solve all secular");
    // tsec.Start();
    
    VectorMem<100> omega(n), omega2diff(n), omega2ref(n);

    omega2diff = 0;
    omega2ref = -1;
    // Vector<> uvec(n), vvec(n);
    // cout << "Ds = " << Ds << endl;
    // cout << "deflated = " << deflated << endl;


    for (int i = 0; i < n; i++)
      {
        if (!deflated[i])
          {
            double l = Ds(i);
            double u = Ds_next(i);

            auto [t,Delta] = SolveSecular (Ds, zs, l, u);

            double t2 = t*(t+2*l);
            double Delta2 = Delta*(2*u-Delta);
          
            if (t < Delta)
              {
                omega2ref(i) = Ds2(i);
                omega2diff(i) = t2;
                omega(i) = Ds(i)+t;
              }
            else
              {
                omega2ref(i) = sqr(Ds_next(i));
                omega2diff(i) = -Delta2;
                omega(i) = Ds_next(i)-Delta;  
              }
          
            // omega(i) = sqrt(omega2ref(i)+omega2diff(i));
          }
        else
          {
            omega(i) = Ds(i);
          }
      }

    // tsec.Stop();


    // static Timer tzmod("modified z");
    // tzmod.Start();
    
    VectorMem<100> zsmod(n);
    zsmod = 0;
    size_t last_nondeflated = 0;
    for (size_t i = 0; i < n; i++)
      if (!deflated[i]) last_nondeflated = i;
    for (size_t i = 0; i < n; i++)
      if (!deflated[i])
        {
          double prod = 1;
          for (size_t j = 0; j < i; j++)
            if (!deflated[j] && (j != last_nondeflated))
              prod *= (omega2ref(j)-Ds2(i)+omega2diff(j)) / (Ds2(j)-Ds2(i));
          for (size_t j = i; j < n-1; j++)
            if (!deflated[j] && (j != last_nondeflated))
              prod *= (omega2ref(j)-Ds2(i)+omega2diff(j)) / (sqr(Ds_next(j))-Ds2(i));
          prod *= omega2ref(last_nondeflated)-Ds2(i)+omega2diff(last_nondeflated);
          zsmod(i) = sqrt(prod);
        }

    for (size_t i = 0; i < n; i++)
      if (zs(i) < 0) zsmod(i)*=-1;
    
    // tzmod.Stop();
    // tMUV.Start();

    
    // Matrix<double,ColMajor> MU(n,n), MV(n,n);   // M = MU * Msigma * MV
    ArrayMem<double,512> memMUV(2*n*n);
    FlatMatrix<double,ColMajor> MU(n,n,memMUV.Data());
    FlatMatrix<double,ColMajor> MV(n,n,memMUV.Data()+n*n);
    
    MU = 0.0;
    MV = 0.0;

    for (size_t i = 0; i < n; i++)
      {
        auto colu = MU.Col(n-i-1);
        auto colv = MV.Col(n-i-1);
        
        if (zs(i) != 0)
          {
            for (size_t j = 0; j < n; j++)
              if (zs(j) != 0)
                {
                  double ui = zsmod(j) / (Ds2(j)-omega2ref(i)-omega2diff(i));
                  colu(index[j]) = ui;
                  colv(index[j]) = Ds(j)*ui;
                }
            colv(index[0]) = -1;           
            colu /= L2Norm(colu);
            colv /= L2Norm(colv);
          }
        else
          {
            colu = 0;
            colu(index[i]) = 1;
            colv = 0;
            colv(index[i]) = 1;
          }
      }

    
    /*
      // testing
    if (double erru = L2Norm(Matrix(Trans(MU)*MU)-Identity(n)); erru > 1e-12)
      cout << "Secular, orthogonal U error = " << erru << endl;
    if (double errv = L2Norm(Matrix(Trans(MV)*MV)-Identity(n)); errv > 1e-12)
      cout << "Secular, orthogonal V error = " << errv << endl;
    */
    
    
    // cout << "Ds2 = " << Ds2 << endl;
    // cout << "omega = " << omega << endl;
    /*
    if (isnan(L2Norm(omega)) || isinf(L2Norm(omega)))
      {
        cout << "omega is illegal" << endl;
        cout << "normmat = " << normmat << endl;
        cout << "omega = " << endl << omega;
        cout << "zs = " << endl << zs << endl;
        cout << "D^2s = " << endl << Ds2 << endl;
        cout << "D_next = " << endl << Ds_next << endl;
        cout << "deflated = " << endl << deflated << endl;
      }
    */
    
    for (size_t i = 0; i < n; i++)
      sigma(n-i-1) = omega(i);

    // tMUV.Stop();
    // tsort.Start();
    
    // we still need sorting, deflated and non-deflated may overlap
    for (int i : Range(n)) index[i] = i;    
    QuickSortI (FlatArray(sigma.Size(), sigma.Data()), index, [](double x, double y) { return y < x; });
    
    VectorMem<100> tmp_sigma(n);
    tmp_sigma = sigma;
    auto tmp_MU = U.Rows(n).Cols(n);  // borrow mem
    auto tmp_MV = V.Rows(n).Cols(n);
    tmp_MU = MU;
    tmp_MV = MV;
    for (int i = 0; i < n; i++)
      {
        sigma(i) = tmp_sigma(index[i]);
        MU.Col(i) = tmp_MU.Col(index[i]);
        MV.Col(i) = tmp_MV.Col(index[i]);
      }
    // tsort.Stop();
    /*
      Matrix Sigma(n);
      Sigma = 0.0; Sigma.Diag() = sigma;
      cout << "shold be same arrow-head:" << endl
      << Truncate( Matrix ( Matrix(MU*Sigma) * Trans(MV) ) ) << endl;
    */

    // static Timer tmatmult("matmult");
    // tmatmult.Start();
    U.Cols(n) = Q * MU;
    U.Col(n) = q; 
    V = W * MV;
    // tmatmult.Stop();

    /*
    if (isnan (L2Norm(sigma)))
      throw Exception ("sigma is none");
    */
    
      // testing

    /*
    double errMU = L2Norm(Trans(MU)*MU-Identity(n));
    double errMV = L2Norm(Trans(MV)*MV-Identity(n));
    
    if (max(errMU, errMV) > 1e-10)
      {
        cout << "MU orthogonality err = " << errMU << endl;
        cout << "MV orthogonality err = " << errMV << endl;

        if (max(errMU, errMV) > 1e-8)
          {
            cout << "MU = " << endl << MU << endl;
            cout << "MV = " << endl << MV << endl;
            cout << "Trans(MU)*MU-I = " << endl << Truncate(Trans(MU)*MU-Identity(n)) << endl;
            cout << "Trans(MV)*MV-I = " << endl << Truncate(Trans(MV)*MV-Identity(n)) << endl;
            cout << "alpha = " << endl << alpha << endl;
            cout << "beta = " << endl << beta << endl;
            // cout << "omega2 = " << omega2 << endl;
            // cout << "omega2ref = " << omega2ref << endl;
            // cout << "omega2diff = " << omega2diff << endl;
            cout << "Ds = " << Ds << endl;
            cout << "zs = " << zs << endl;
            cout << "omega = " << omega << endl;
            cout << "degen = " << deflated << endl;
            if (max(errMU,errMV) > 1e-3)
              throw Exception ("bad orthogonality");
          }
      }
    */


    
    /*
      // testing 
    Matrix Sigma = Matrix(Trans(U)*B)*V;
    Sigma.Diag(0) = 0.0;
    if (double errsigma = L2Norm(Sigma); errsigma > 1e-4)
      {
        cout << "in DnC, bidiag, sigma err = " << errsigma << endl;
        if (Sigma.Height() <= 6)
          {
            cout << "Sigma clear diag = " << endl << Sigma << endl;
            cout << "Sigma clear diag = " << endl << Truncate(Sigma, 1e-10) << endl;
            Sigma = Matrix(Trans(U)*B)*V;
            cout << "bidiag, Sigma = " << endl << Truncate(Sigma) << endl;
            cout << "alpha = " << endl << alpha << endl;
            cout << "beta = " << endl << beta << endl;
            cout << "deflation = " << endl << deflated << endl;
            cout << "zs = " << endl << zs << endl;
            cout << "zsmod = " << endl << zsmod << endl;
            cout << "Ds = " << endl << Ds << endl;
          }
        // throw Exception("wrong");
      }
    */
  }


  
  /*
    B ... n*n,   uppper bi-diagonal
    Compute   B = U Sigma V^T
    with Sigma stored as diag B
  */
  template <ORDERING ORD>
  void CalcSVDBiDiagonal (SliceMatrix<double,ORD> B, 
                          SliceMatrix<double,ColMajor> U, SliceMatrix<double,ColMajor> V)
  {
    static Timer t("CalcSVDBiDiagonal"); RegionTimer r(t);

    // Matrix<> Borig = B; // for testing
    size_t n = B.Height();
    Vector<> alpha(n), beta(n), sigma(n);

    alpha(0) = 0;
    alpha.Range(1,n) = B.Diag(1);
    beta = B.Diag(0);

    Matrix<double,ColMajor> U1(n+1, n+1);
    CalcSVDBiDiagonalRec (alpha, beta, sigma, U1, V);
    U = U1.Rows(1,n+1).Cols(0,n);
    B = 0.0;
    B.Diag(0) = sigma;


    /*
      // testing 
    Matrix Sigma = Matrix(Trans(U)*Borig)*V;
    Sigma.Diag(0) = 0.0;
    if (double errsigma = L2Norm(Sigma); errsigma > 1e-4)
      {
        cout << "bidiag, sigma err = " << errsigma << endl;
        if (Sigma.Height() <= 6)
          {
            cout << "Sigma clear diag = " << endl << Sigma << endl;
            cout << "Sigma clear diag = " << endl << Truncate(Sigma, 1e-10) << endl;
            Sigma = Matrix(Trans(U)*Borig)*V;
            cout << "bidiag, Sigma = " << endl << Truncate(Sigma) << endl;
            cout << "alpha = " << endl << alpha << endl;
            cout << "beta = " << endl << beta << endl;
          }
        // throw Exception("wrong");
      }
    */
  }
    


  // A = U * Sigma * V^T
  template <ORDERING ORD>
  void TCalcSVD (SliceMatrix<double,ORD> A, 
                 SliceMatrix<double, ColMajor> U, SliceMatrix<double, ColMajor> V)
  {
    /*
    Matrix Aorig = A;
    Matrix<double,ColMajor> hU = U;
    LapackSVD (A, V, hU);
    U = Trans (hU);
    // cout << "have svd, sigma = " << SliceMatrix(A).Diag(0) << endl;
    // cout << "Trans(U) * A * V = " << Truncate(Trans(U) * Aorig * V) << endl;
    return;
    */
    
    if (A.Width() > A.Height())
      {
        TCalcSVD (Trans(A), V, U);
        return;
      }
    
    // cout << "CalcSVD, h = " << A.Height() << ", w = " << A.Width() << endl;
    if (double norm = L2Norm(A); isnan(norm) || isinf(norm))
      {
        cout << "input matrix norm = " << norm << endl;
        cout << "mat = " << A << endl;
        throw Exception ("called SVD with nan-matrix");
      }
        
    static Timer t("CalcSVD"); RegionTimer reg(t);
    size_t m = A.Height();
    size_t n = A.Width();

    // Matrix<> A_orig = A;   // for testing 
  
    Matrix<double, ColMajor> U1(m);
    Matrix<double, ColMajor> V1(n);

    if (min(n,m) < 100)
      ReduceBiDiagonal (A, U1, V1);
    else
      ReduceBiDiagonalBlocked<ORD,ColMajor> (A, U1, V1);

    // testting
    /*
      if (m < 50)  
      {  
      Matrix B = Matrix(Trans(U1) * A_orig) * V1;   // bi-diagonal
      cout << "U^T * A V" << endl << Truncate(B) << endl;
      }
    */
    
    Matrix<double,ColMajor> UB(m), VB(n);
    UB = Identity(m);
    // cout << "diag A = " << endl << A.Diag(0) << endl;
    // cout << "super-diag A = " << endl << A.Diag(1) << endl;

    CalcSVDBiDiagonal (A.Rows(n), UB.Rows(n).Cols(n), VB);
    
    static Timer tmultUV("CalcSVD, mult U1*UB, V1*VB"); RegionTimer regUV(tmultUV);
    tmultUV.AddFlops (n*n*n+m*m*m);
    U = U1*UB;
    V = V1*VB;


    /* 
       // testing ... 
    double errU = L2Norm (U * Trans(U) - Identity(U.Height()));
    double errV = L2Norm (V * Trans(V) - Identity(V.Height()));
    if (errU + errV > 1e-2)
      {
        cout << "orthogonaliry error = " << errU + errV << endl;
        // cout << "Aorig = " << A_orig << endl;
        // cout << "diag = " << A.Diag(0) << endl;
        // cout << "superd = " << A.Diag(1) << endl;
        cout << "err ortho U: " << errU << endl;
        cout << "err ortho V: " << errV << endl;
        cout << "err ortho U1: " << L2Norm (U1 * Trans(U1) - Identity(U1.Height())) << endl;
        cout << "err ortho V1: " << L2Norm (V1 * Trans(V1) - Identity(V1.Height())) << endl;
        cout << "err ortho UB: " << L2Norm (UB * Trans(UB) - Identity(UB.Height())) << endl;
        cout << "err ortho VB: " << L2Norm (VB * Trans(VB) - Identity(VB.Height())) << endl;
        // cout << "UB = " << endl << UB << endl;
        // cout << "VB = " << endl << VB << endl;
        if (errU+errV > 1e-4)
          throw Exception ("got a bad SVD");
      }
    */

    /*
    // testing ...
    Matrix Sigma = Matrix(Trans(U)*A_orig)*V;
    Sigma.Diag(0) = 0.0;
    if (double errsigma = L2Norm(Sigma); errsigma > 1e-2)
      {
        cout << "sigma err = " << errsigma << endl;
        Sigma = Matrix(Trans(U)*A_orig)*V;
        cout << "Sigma = " << endl << Truncate(Sigma) << endl;
        cout << "A_orig = " << endl << A_orig << endl;
        cout << "diag = " << endl << A.Diag(0) << endl;
        cout << "super-diag = " << endl << A.Diag(1) << endl;
        if (errsigma > 1e-2)
          throw Exception("Sigma not diagonal");
      }
    */
    // cout << "sigma = " << A.Diag(0) << endl;
  }


  void CalcSVD (SliceMatrix<double,RowMajor> A, 
                SliceMatrix<double, ColMajor> U, SliceMatrix<double, ColMajor> V)
  {
    TCalcSVD (A, U, V);
  }
  void CalcSVD (SliceMatrix<double,ColMajor> A, 
                SliceMatrix<double, ColMajor> U, SliceMatrix<double, ColMajor> V)
  {
    TCalcSVD (A, U, V);
  }

  
}
