#include <bla.hpp>

namespace ngbla
{


  /*
  // find Householder reflection vector v such that
  // reflection matrix H_v = I - 2 v v^T / (v^T v)
  // leads to H_v x = +/- e_0 ||x||
  void CalcHouseholderVector (SliceVector<> x, FlatVector<> v)
  {
    double norm = L2Norm(x);
    v(0) += (v(0) < 0) ? -norm : norm;
    double v0 = v(0);
    if (v0 != 0)
      v *= 1/v0;
  }
  */


  
  // find Householder reflection vector v such that
  // reflection matrix H_v = I - 2 v v^T / (v^T v)
  // leads to H_v x = +/- e_0 ||x||
  // scaling such that v(0) = +1   
  // returns (H_v x)(0) = +/- ||x||
  double CalcHouseholderVectorInPlace (SliceVector<> x)
  {
    double signed_norm = L2Norm(x);

    // have to treat this case,
    // since norm can be zero for very-small nonzero x
    if (signed_norm == 0)  
      { 
        x(0) = 1;
        return 0;
      }
    
    if (x(0) > 0) signed_norm *= -1;
    
    double v0 = x(0) - signed_norm;
    
    if (v0 != 0)
      x.Range(1,x.Size()) *= 1/v0;
    x(0) = 1;
    return signed_norm;
  }

  
  HouseholderReflection :: HouseholderReflection (FlatVector<> av)
    : v(av)
  {
    factor = L2Norm2(v);
    if (factor != 0) factor = 2/factor;
  }

  void HouseholderReflection :: Mult (FlatVector<double> v2) const
  {
    double ip = InnerProduct(v, v2);
    v2 -= (factor*ip) * v;
  }
  
  template <ORDERING ORD>
  void HouseholderReflection :: TMult (SliceMatrix<double,ORD> m2) const  
  {
    // const char * timername = (ORD == ColMajor)
    // ? "Householder, colmajor" : "Householder, rowmajor";    
    // static Timer tcolmajor(timername); RegionTimer reg(tcolmajor);
    // tcolmajor.AddFlops (2*v.Size()*m2.Width());
    
    constexpr size_t bs = 24;
    double mem[bs];

    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        FlatVector<> hv(bsi, &mem[0]);
        auto colsm2 = m2.Cols(i, i+bsi);
        hv = Trans(colsm2) * v;
        hv *= factor;
        // colsm2 -= v * Trans(hv);
        colsm2 -= FlatMatrix<double, ORD> (v.Size(), 1, v.Data()) *
          FlatMatrix<double,ORD> (1, hv.Size(), hv.Data());
      }
  }

  template void HouseholderReflection :: TMult (SliceMatrix<double,RowMajor> m2) const;
  template void HouseholderReflection :: TMult (SliceMatrix<double,ColMajor> m2) const;




  HouseholderReflection1 :: HouseholderReflection1 (SliceVector<> av)
    : v(av)
  {
    factor = 1+L2Norm2(v.Range(1, v.Size()));
    factor = 2/factor;
  }

  
  template <ORDERING ORD>
  void HouseholderReflection1 :: TMult (SliceMatrix<double,ORD> m2) const  
  {
    // const char * timername = (ORD == ColMajor)
    // ? "Householder, colmajor" : "Householder, rowmajor";    
    // static Timer tcolmajor(timername); RegionTimer reg(tcolmajor);
    // tcolmajor.AddFlops (2*v.Size()*m2.Width());
    
    constexpr size_t bs = 96;
    double mem[bs];
    size_t m = v.Size();
    
    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        FlatVector<> hv(bsi, &mem[0]);
        auto colsm2 = m2.Cols(i, i+bsi);

        hv = colsm2.Row(0);
        hv = Trans(colsm2).Cols(1,m) * v.Range(1,m);
        hv *= factor;
        // colsm2 -= v * Trans(hv);
        colsm2.Row(0) -= hv;
        colsm2.Rows(1,m) -= v.Range(1,m) * Trans(hv);
      }
  }

  template void HouseholderReflection1 :: TMult (SliceMatrix<double,RowMajor> m2) const;
  template void HouseholderReflection1 :: TMult (SliceMatrix<double,ColMajor> m2) const;



  

  // H = H_{m-1} ... H_1 H_0 = I - V^T T V
  // SliceMatrix<> mv;  // every row one reflection vector
  // Matrix<> T;        // strict lower triangular

  template <ORDERING OMV>
  void BaseMultiHouseholderReflection<OMV> :: CalcT()
  {
    static Timer t("multiHouseholder, ctor"); RegionTimer reg(t);
    size_t m = mv.Height();

    //
    // Chiara Puglisi
    // Modification of the Householder method based on the compact WY rep
    // Q = I - Y T Y^T
    //
    // T = mv * Trans(mv);  // inner products (TODO: only strict lower trig is needed)

    auto [mv1,mv2] = mv.SplitCols(m);
    T = Trans(mv1);
    T.Diag() = 1.0;
    for (size_t i = 0; i < m; i++)
      T.Row(i).Range(i+1, m) = 0;
    TriangularMult<UpperRight, Normalized> (mv1, T);
      
    T += mv2 * Trans(mv2);

    T.Diag() *= 0.5;
    for (auto & d : T.Diag())
      if (d == 0) d = 1;

    TriangularInvert<LowerLeft> (T);
  }

  
  template <ORDERING OMV> template <ORDERING ORD>
  void BaseMultiHouseholderReflection<OMV> :: TMult (SliceMatrix<double,ORD> m2) const  // Hm-1 * ... * H1 * H0 * m2
  {
    const char * timername = (OMV == ColMajor)
      ? ((ORD == ColMajor) ? "multiHouseholder, H..colmajor, M..colmajor" : "multiHouseholder, H..colmajor, M..rowmajor") 
      : ((ORD == ColMajor) ? "multiHouseholder, H..rowmajor, M..colmajor" : "multiHouseholder, H..rowmajor, M..rowmajor");
    static Timer t(timername); RegionTimer reg(t);
    t.AddFlops (2*mv.Height()*m2.Height()*m2.Width());

    /*
    // naive version
    for (size_t j  = 0; j < mv.Height(); j++)
      HouseholderReflection(mv.Row(j)).Mult(m2);
    return;
    */

    constexpr size_t bs = 96;
    ArrayMem<double,48*96> mem(bs*mv.Height());
    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        auto colsm2 = m2.Cols(i, i+bsi);
        
        FlatMatrix<double,ORD> tmp(mv.Height(), bsi, &mem[0]);
        {
          // static Timer t(timername+string("1")); RegionTimer reg(t);
          GeneralizedTriangularMult<UpperRight,Normalized> (mv, colsm2, tmp);
        }
        {
          // static Timer t(timername+string("2")); RegionTimer reg(t);        
          TriangularMult<LowerLeft,NonNormalized> (T, tmp);
        }

        GeneralizedTriangularSub<LowerLeft,Normalized> (Trans(mv), tmp, colsm2);
      }
  }





  template <ORDERING OMV> template <ORDERING ORD>
  void BaseMultiHouseholderReflection<OMV> :: TMultTrans (SliceMatrix<double,ORD> m2) const  // Hm-1 * ... * H1 * H0 * m2
  {
    const char * timername = (ORD == ColMajor)
      ? "multiHouseholder trans, colmajor" : "multiHouseholder trans, rowmajor";
    static Timer t(timername); RegionTimer reg(t);
    // t.AddFlops (2*mv.Height()*m2.Height()*m2.Width());

    constexpr size_t bs = 96;
    ArrayMem<double,48*96> mem(bs*mv.Height());    
    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        auto colsm2 = m2.Cols(i, i+bsi);
        
        FlatMatrix<double,ORD> tmp(mv.Height(), bsi, &mem[0]);        
        GeneralizedTriangularMult<UpperRight,Normalized> (mv, colsm2, tmp);
        
        TriangularMult<UpperRight,NonNormalized> (Trans(T), tmp);

        GeneralizedTriangularSub<LowerLeft,Normalized> (Trans(mv), tmp, colsm2);
      }
  }





  
  template class MultiHouseholderReflection<RowMajor>;
  template class MultiHouseholderReflection<ColMajor>;
  template void BaseMultiHouseholderReflection<RowMajor> :: TMult (SliceMatrix<double,RowMajor> m2) const;
  template void BaseMultiHouseholderReflection<RowMajor> :: TMult (SliceMatrix<double,ColMajor> m2) const;
  template void BaseMultiHouseholderReflection<ColMajor> :: TMult (SliceMatrix<double,RowMajor> m2) const;
  template void BaseMultiHouseholderReflection<ColMajor> :: TMult (SliceMatrix<double,ColMajor> m2) const;
  template void BaseMultiHouseholderReflection<RowMajor> :: TMultTrans (SliceMatrix<double,RowMajor> m2) const;
  template void BaseMultiHouseholderReflection<RowMajor> :: TMultTrans (SliceMatrix<double,ColMajor> m2) const;
  template void BaseMultiHouseholderReflection<ColMajor> :: TMultTrans (SliceMatrix<double,RowMajor> m2) const;
  template void BaseMultiHouseholderReflection<ColMajor> :: TMultTrans (SliceMatrix<double,ColMajor> m2) const;


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
