#include <bla.hpp>

namespace ngbla
{

  // find Householder reflection vector v such that
  // reflection matrix H_v = I - v v^T
  // leads to H_v x = +/- e_0
  void CalcHouseholderVector (SliceVector<> x, FlatVector<> v)
  {
    v = x;
    double norm = L2Norm(x);
    v(0) += (v(0) < 0) ? -norm : norm;
    double normv = L2Norm(v);
    if (normv != 0.0)
      v *= sqrt(2.0) / normv;
  }


  void HouseholderReflection :: Mult (SliceMatrix<> m2) const
  {
    static Timer t("Householder, rowmajor"); RegionTimer reg(t);
    t.AddFlops (2*v.Size()*m2.Width());
    
    /*
      VectorMem<2000> hv(m2.Width());
      hv = Trans(m2)*v;
      // m2 -= v * Trans(hv);
      m2 -=
      FlatMatrix<> (v.Size(), 1, v.Data()) *
      FlatMatrix<> (1, hv.Size(), hv.Data());
    */
    
    constexpr size_t bs = 96;
    double mem[bs];
    
    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        FlatVector<> hv(bsi, &mem[0]);
        auto colsm2 = m2.Cols(i, i+bsi);
        hv = Trans(colsm2) * v;
        colsm2 -= v * Trans(hv);
      }
  }
  
  void HouseholderReflection :: Mult (SliceMatrix<double,ColMajor> m2) const
  {
    static Timer tcolmajor("Householder, colmajor"); RegionTimer reg(tcolmajor);
    tcolmajor.AddFlops (2*v.Size()*m2.Width());

    /*
      VectorMem<2000> hv(m2.Width());
      hv = Trans(m2)*v;
      // m2 -= v * Trans(hv);
      m2 -=
      FlatMatrix<> (v.Size(), 1, v.Data()) *
      FlatMatrix<> (1, hv.Size(), hv.Data());
    */
    
    constexpr size_t bs = 96;
    double mem[bs];

    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        FlatVector<> hv(bsi, &mem[0]);
        auto colsm2 = m2.Cols(i, i+bsi);
        hv = Trans(colsm2) * v;
        colsm2 -= v * Trans(hv);
      }
  }


  // H = H_{m-1} ... H_1 H_0 = I - V^T T V
  // SliceMatrix<> mv;  // every row one reflection vector
  // Matrix<> T;        // strict lower triangular

  
  MultiHouseholderReflection :: MultiHouseholderReflection (SliceMatrix<> amv) : mv(amv), T(amv.Height())
  {
    static Timer t("multiHouseholder, ctor"); RegionTimer reg(t);
    size_t m = mv.Height();
    Matrix<> ip = mv * Trans(mv);  // inner products (TODO: only strict lower trig is needed)

    /*
    // original version by Schreiber + Van Loan, or Bischof + Van Loan
    T = Identity(m);
    for (int i = 1; i < m; i++)
    T.Row(i).Range(0,i) = -Trans(T.Rows(0,i).Cols(0,i)) * ip.Row(i).Range(0,i);
    return;
    */
    
    //
    // Chiara Puglisi
    // Modification of the Householder method based on the compact WY rep
    // Q = I - Y T Y^T
    // 
    T = Identity(m);
    for (size_t i = 1; i < m; i++)
      T.Row(i).Range(i) = ip.Row(i).Range(i);
    TriangularInvert<LowerLeft, Normalized> (T);
  }

  void MultiHouseholderReflection :: Mult (SliceMatrix<> m2) const  // Hm-1 * ... * H1 * H0 * m2
  {
    static Timer t("multiHouseholder, rowmajor"); RegionTimer reg(t);
    t.AddFlops (2*mv.Height()*m2.Height()*m2.Width());

    /*
    // naive version
    for (size_t j  = 0; j < mv.Height(); j++)
    HouseholderReflection(mv.Row(j)).Mult(m2);
    return;
    */
    
    /*
    // matrix products
    Matrix<> tmp = mv * m2;
    Matrix<> tmp2 = T * tmp;
    m2 -= Trans(mv)*tmp2;  
    */
    
    // matrix products with blocking    
    constexpr size_t bs = 96;
    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        auto colsm2 = m2.Cols(i, i+bsi);   // reuse vertical m2 panel
        Matrix<> tmp = mv * colsm2;
        TriangularMult<LowerLeft,Normalized> (T, tmp);
        colsm2 -= Trans(mv)*tmp;  
      }
  }
  
  void MultiHouseholderReflection :: Mult (SliceMatrix<double,ColMajor> m2) const  // Hm-1 * ... * H1 * H0 * m2
  {
    static Timer t("multiHouseholder, colmajor"); RegionTimer reg(t);
    t.AddFlops (2*mv.Height()*m2.Height()*m2.Width());

    /*
      Matrix<> tmp = mv * m2;
      Matrix<> tmp2 = T * tmp;
      m2 -= Trans(mv)*tmp2;
    */
    constexpr size_t bs = 96;
    for (size_t i = 0; i < m2.Width(); i += bs)
      {
        size_t bsi = min(bs, m2.Width()-i);
        auto colsm2 = m2.Cols(i, i+bsi);
        Matrix<> tmp = mv * colsm2;
        TriangularMult<LowerLeft,Normalized> (T, tmp);        
        colsm2 -= Trans(mv)*tmp;  
      }
  }


  /*
    Factorize A = Q R
    A is m x n, overwritten by R
    Q is either m x m      (TODO: or m x n)
  */
  void QRFactorization (SliceMatrix<> A, SliceMatrix<> Q)
  {
    static Timer t("QRFactorization"); RegionTimer reg(t);
    size_t m = A.Height();
    size_t n = A.Width();
    Q = Identity (m);

    /*
    // non-blocking (for testing)
    Vector<> v(m);
    for (size_t i = 0; i < min(n,m-1); i++)
    {
    auto x = A.Col(i).Range(i, m);
    CalcHouseholderVector (x, v.Range(i,m));
    HouseholderReflection H(v.Range(i,m));
    H.Mult (A.Rows(i,m).Cols(i,n));
    H.Mult (Trans(Q).Rows(i,m));
    }
    */

    // blocking
    constexpr size_t bs = 32;
    Matrix<> mv(bs, m);
    size_t k = min(n,m-1);
    for (size_t i1 = 0; i1 < k; i1 += bs)
      {
        size_t bsi = min(bs, k-i1);
        mv.Rows(0,bsi).Cols(i1,i1+bsi) = 0.0;
        for (size_t i2 = 0; i2 < bsi; i2++)
          {
            auto i = i1+i2;
            auto x = A.Col(i).Range(i, m);
            CalcHouseholderVector (x, mv.Row(i2).Range(i,m));
            HouseholderReflection H(mv.Row(i2).Range(i,m));
            H.Mult (A.Rows(i,m).Cols(i,i1+bsi));
          }
        MultiHouseholderReflection H(mv.Rows(bsi).Cols(i1,m));
        H.Mult (Trans(Q).Rows(i1,m));
        H.Mult (A.Rows(i1,m).Cols(i1+bsi,n));
      }
  }

}
