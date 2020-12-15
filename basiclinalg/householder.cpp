#include <bla.hpp>

namespace ngbla
{

  
  // find Householder reflection vector v such that
  // reflection matrix H_v = I - 2 v v^T / (v^T v)
  // leads to H_v x = +/- e_0 ||x||
  void CalcHouseholderVector (SliceVector<> x, FlatVector<> v)
  {
    v = x;
    double norm = L2Norm(x);
    v(0) += (v(0) < 0) ? -norm : norm;
    double v0 = v(0);
    if (v0 != 0)
      v *= 1/v0;
  }


  // find Householder reflection vector v such that
  // reflection matrix H_v = I - 2 v v^T / (v^T v)
  // leads to H_v x = +/- e_0 ||x||
  // scaling such that v(0) = +1   
  // returns (H_v x)(0) = +/- ||x||
  double CalcHouseholderVectorInPlace (SliceVector<> x)
  {
    double signed_norm = L2Norm(x);
    if (x(0) > 0) signed_norm *= -1;
    double v0 = x(0) - signed_norm;
    // cout << "refection in place, norm = " << signed_norm << ", v0 = " << v0 << endl;
    // if (fabs(v0) > 1e-100)
    if (v0 != 0)
      x.Range(1,x.Size()) *= 1/v0;
    else
      {  // with very small numbers it can be Norm(x)=0, but x != 0
        x = 0.0;
        signed_norm = 0.0;
      }
    x(0) = 1;
    return signed_norm;
  }

  
  HouseholderReflection :: HouseholderReflection (FlatVector<> av)
    : v(av)
  {
    factor = L2Norm2(v);
    if (factor != 0) factor = 2/factor;
  }

  
  template <ORDERING ORD>
  void HouseholderReflection :: TMult (SliceMatrix<double,ORD> m2) const  
  {
    // const char * timername = (ORD == ColMajor)
    // ? "Householder, colmajor" : "Householder, rowmajor";    
    // static Timer tcolmajor(timername); RegionTimer reg(tcolmajor);
    // tcolmajor.AddFlops (2*v.Size()*m2.Width());
    
    constexpr size_t bs = 96;
    double mem[bs];

    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        FlatVector<> hv(bsi, &mem[0]);
        auto colsm2 = m2.Cols(i, i+bsi);
        
        hv = Trans(colsm2) * v;
        hv *= factor;
        colsm2 -= v * Trans(hv);
      }
  }

  template void HouseholderReflection :: TMult (SliceMatrix<double,RowMajor> m2) const;
  template void HouseholderReflection :: TMult (SliceMatrix<double,ColMajor> m2) const;

  

  // H = H_{m-1} ... H_1 H_0 = I - V^T T V
  // SliceMatrix<> mv;  // every row one reflection vector
  // Matrix<> T;        // strict lower triangular
  template <ORDERING OMV>
  MultiHouseholderReflection<OMV> :: MultiHouseholderReflection (SliceMatrix<double,OMV> amv)
    : mv(amv), T(amv.Height())
  {
    // static Timer t("multiHouseholder, ctor"); RegionTimer reg(t);
    size_t m = mv.Height();

    //
    // Chiara Puglisi
    // Modification of the Householder method based on the compact WY rep
    // Q = I - Y T Y^T
    //
    // T = mv * Trans(mv);  // inner products (TODO: only strict lower trig is needed)

    auto [mv1,mv2] = mv.SplitCols(m);
    Matrix<> U = mv1; 
    U.Diag() = 1.0;
    for (size_t i = 1; i < m; i++)
      U.Row(i).Range(0,i) = 0;

    T = U * Trans(U);
    T += mv2 * Trans(mv2);

    T.Diag() *= 0.5;
    for (auto & d : T.Diag())
      if (d == 0) d = 1;

    TriangularInvert<LowerLeft> (T);
  }

  
  template <ORDERING OMV> template <ORDERING ORD>
  void MultiHouseholderReflection<OMV> :: TMult (SliceMatrix<double,ORD> m2) const  // Hm-1 * ... * H1 * H0 * m2
  {
    // const char * timername = (ORD == ColMajor)
    // ? "multiHouseholder, colmajor" : "multiHouseholder, rowmajor";
    // static Timer t(timername); RegionTimer reg(t);
    // t.AddFlops (2*mv.Height()*m2.Height()*m2.Width());

    /*
    // naive version
    for (size_t j  = 0; j < mv.Height(); j++)
      HouseholderReflection(mv.Row(j)).Mult(m2);
    return;
    */

    constexpr size_t bs = 96;
    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        auto colsm2 = m2.Cols(i, i+bsi);
        
        // Matrix<> tmp = mv * colsm2;
        Matrix<double,ORD> tmp(mv.Height(), bsi);
        {
          // static Timer t(timername+string("1")); RegionTimer reg(t);
          GeneralizedTriangularMult<UpperRight,Normalized> (mv, colsm2, tmp);
        }
        {
          // static Timer t(timername+string("2")); RegionTimer reg(t);        
          TriangularMult<LowerLeft,NonNormalized> (T, tmp);
        }
        // colsm2 -= Trans(mv)*tmp;
        IntRange r1(0, mv.Height());
        IntRange r2(mv.Height(), colsm2.Height());
        {
          // static Timer t(timername+string("3")); RegionTimer reg(t);        
          colsm2.Rows(r2) -= Trans(mv.Cols(r2)) * tmp;
        }
        {
          // static Timer t(timername+string("4")); RegionTimer reg(t);        
          TriangularMult<LowerLeft,Normalized> (Trans(mv.Cols(r1)), tmp);
        }
        colsm2.Rows(r1) -= tmp;
      }
  }





  template <ORDERING OMV> template <ORDERING ORD>
  void MultiHouseholderReflection<OMV> :: TMultTrans (SliceMatrix<double,ORD> m2) const  // Hm-1 * ... * H1 * H0 * m2
  {
    // const char * timername = (ORD == ColMajor)
    // ? "multiHouseholder trans, colmajor" : "multiHouseholder trans, rowmajor";
    // static Timer t(timername); RegionTimer reg(t);
    // t.AddFlops (2*mv.Height()*m2.Height()*m2.Width());

    constexpr size_t bs = 96;
    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        auto colsm2 = m2.Cols(i, i+bsi);
        
        Matrix<double,ORD> tmp(mv.Height(), bsi);
        GeneralizedTriangularMult<UpperRight,Normalized> (mv, colsm2, tmp);
        
        // TriangularMult<LowerLeft,NonNormalized> (T, tmp);
        TriangularMult<UpperRight,NonNormalized> (Trans(T), tmp);

        IntRange r1(0, mv.Height());
        IntRange r2(mv.Height(), colsm2.Height());

        colsm2.Rows(r2) -= Trans(mv.Cols(r2)) * tmp;
        
        TriangularMult<LowerLeft,Normalized> (Trans(mv.Cols(r1)), tmp);
        colsm2.Rows(r1) -= tmp;
      }
  }





  
  template class MultiHouseholderReflection<RowMajor>;
  template class MultiHouseholderReflection<ColMajor>;
  template void MultiHouseholderReflection<RowMajor> :: TMult (SliceMatrix<double,RowMajor> m2) const;
  template void MultiHouseholderReflection<RowMajor> :: TMult (SliceMatrix<double,ColMajor> m2) const;
  template void MultiHouseholderReflection<ColMajor> :: TMult (SliceMatrix<double,RowMajor> m2) const;
  template void MultiHouseholderReflection<ColMajor> :: TMult (SliceMatrix<double,ColMajor> m2) const;
  template void MultiHouseholderReflection<RowMajor> :: TMultTrans (SliceMatrix<double,RowMajor> m2) const;
  template void MultiHouseholderReflection<RowMajor> :: TMultTrans (SliceMatrix<double,ColMajor> m2) const;
  template void MultiHouseholderReflection<ColMajor> :: TMultTrans (SliceMatrix<double,RowMajor> m2) const;
  template void MultiHouseholderReflection<ColMajor> :: TMultTrans (SliceMatrix<double,ColMajor> m2) const;



  template <ORDERING ORDER>
  void T_QRFactorizationInPlace (SliceMatrix<double, ORDER> A)
  {
    static Timer t("QRFactorization inplace"); RegionTimer reg(t);

    size_t m = A.Height();
    size_t n = A.Width();

    /*
    // non-blocking (for comparison)
    Vector<> v(m);
    for (size_t i = 0; i < min(n,m-1); i++)
      {
        auto x = A.Col(i).Range(i, m);
        CalcHouseholderVector (x, v.Range(i,m));
        HouseholderReflection H(v.Range(i,m));
        H.Mult (A.Rows(i,m).Cols(i,n));
        x.Range(1, x.Size()) = v.Range(i+1,m);
      }
    return;
    */

    // blocking
    constexpr size_t bs = 32;
    size_t k = min(n,m-1);
    for (size_t i1 = 0; i1 < k; i1 += bs)
      {
        size_t bsi = min(bs, k-i1);
        for (size_t i = i1; i < i1+bsi; i++)
          {
            auto Apanel = A.Cols(i,i1+bsi).Rows(i,m);

            double signed_norm = CalcHouseholderVectorInPlace (Apanel.Col(0));
            Apanel(0,0) = signed_norm;

            MultiHouseholderReflection H(Trans(Apanel.Cols(0,1)));
            H.Mult (Apanel.Cols(1, Apanel.Width()));
          }

        MultiHouseholderReflection H(Trans(A.Cols(i1, i1+bsi).Rows(i1, m)));
        H.Mult (A.Rows(i1,m).Cols(i1+bsi,n));
      }
  }

  NGS_DLL_HEADER void QRFactorizationInPlace (SliceMatrix<double,RowMajor> A)
  {
    T_QRFactorizationInPlace (A);
  }
  
  NGS_DLL_HEADER void QRFactorizationInPlace (SliceMatrix<double,ColMajor> A)
  {
    T_QRFactorizationInPlace (A);    
  }

  
  /*
    Factorize A = Q R
    A is m x n, overwritten by R
    Q is either m x m   or m x n
  */
  void QRFactorization (SliceMatrix<> A, SliceMatrix<> Q)
  {
    static Timer t("QRFactorization"); RegionTimer reg(t);
    static Timer tQ("QRFactorization - factor Q");

    QRFactorizationInPlace (A);
    
    size_t m = A.Height();
    size_t n = A.Width();

    size_t minnmA = min(A.Height(), A.Width());
    size_t minnmQ = min(Q.Height(), Q.Width());

    Q = 0.0;
    Q.Rows(minnmQ).Cols(minnmQ).Diag() = 1;
    
    tQ.Start();
    size_t bs = 32;
    size_t k = min(n,m-1);    
    for (size_t i = k; i+bs > bs; i -= bs)
      {
        size_t i1 = (i >= bs) ? i-bs : 0;
        MultiHouseholderReflection H(Trans(A.Cols(i1,i).Rows(i1, m)));
        H.MultTrans (Q.Rows(i1,m).Cols(i1, Q.Width()));
      }
    tQ.Stop();

    // eliminate lower triangle of R
    for (size_t i = 1; i < minnmA; i++)
      A.Row(i).Range(i) = 0.0;
    if (m > n)
      A.Rows(n,m) = 0.0;
  }

}
