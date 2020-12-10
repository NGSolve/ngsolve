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

  
  void ReduceBiDiagonal (SliceMatrix<> A, 
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



  template <ORDERING ORD>
  void ReduceBiDiagonalBlocked (SliceMatrix<> A, 
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
    
    ttrig.Stop();

    tUV.Start();

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
    
    tUV.Stop();

    
    
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


  void CalcSVDBiDiagonal (SliceMatrix<> B, 
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





  // A = U * Sigma * V^T
  void CalcSVD (SliceMatrix<> A, 
                SliceMatrix<double, ColMajor> U, SliceMatrix<double, ColMajor> V)
  {
    Matrix Aorig = A;
    Matrix<double,ColMajor> hU = U;
    LapackSVD (A, V, hU);
    U = Trans (hU);
    // cout << "have svd, sigma = " << SliceMatrix(A).Diag(0) << endl;
    // cout << "Trans(U) * A * V = " << Truncate(Trans(U) * Aorig * V) << endl;
    return;
    
    static Timer t("CalcSVD"); RegionTimer reg(t);
    size_t m = A.Height();
    size_t n = A.Width();

    Matrix<> a_orig = A;   // for testing 

  
    Matrix<double, ColMajor> U1(m);
    Matrix<double, ColMajor> V1(n);

    if (min(n,m) < 100)
      ReduceBiDiagonal (A, U1, V1);
    else
      ReduceBiDiagonalBlocked<ColMajor> (A, U1, V1);

    // testting
    /*
    if (m < 50)  
      {  
        Matrix B = Matrix(Trans(U1) * a_orig) * V1;   // bi-diagonal
        cout << "U^T * A V" << endl << Truncate(B) << endl;
      }
    */
    
    Matrix<double,ColMajor> UB(m), VB(n);
    CalcSVDBiDiagonal (A, UB, VB);

    static Timer tmultUV("CalcSVD, mult U1*UB, V1*VB"); RegionTimer regUV(tmultUV);
    tmultUV.AddFlops (n*n*n+m*m*m);
    U = U1*UB;
    V = V1*VB;
  }




#ifdef NONE

  int main()
  {
    cout.precision(15);
    // size_t n = 5, m = 8;    // for now:  n <= m
    size_t n = 1000, m = 1000;    // for now:  n <= m 
    Matrix<> a(m,n);
    for (int i = 0; i < a.Height(); i++)
      for (int j = 0; j < a.Width(); j++)
        a(i,j) = rand() / double(RAND_MAX);
    Matrix a_orig = a;

    if (n <= 1000)
      {
        Matrix alapack = a;
        LapackSVD (alapack);
      }
  
    if (n <= 200)
      {
        Matrix aTa = Trans(a)*a;
        Vector lams(n);
        Matrix evecs(n,n);
        CalcEigenSystem (aTa, lams, evecs);
        for (auto & lami : lams)
          lami = sqrt(lami);
        cout << "sing vals = " << lams;
      }
  
    Matrix<double, ColMajor> U(m,m);
    Matrix<double, ColMajor> V(n,n);
    CalcSVD (a, U, V);
  
    Matrix Sigma = Matrix(Trans(U) * a_orig) * V;   // diagonal

    if (m < 50)
      cout << "u^T * a v" << endl << Truncate(Sigma) << endl;

    NgProfiler::Print(stdout);
  }
#endif

}
