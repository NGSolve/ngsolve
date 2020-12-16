#ifndef FILE_HOUSEHOLDER
#define FILE_HOUSEHOLDER

/****************************************************************************/
/* File:   householder.hpp                                                  */
/* Author: Joachim Schoeberl                                                */
/* Date:   Nov 2020                                                         */
/****************************************************************************/

namespace ngbla
{



  
  // find Householder reflection vector v such that
  // reflection matrix H_v = I - 2 v v^T / (v^T v)
  // leads to H_v x = +/- e_0 ||x||
  // scaling such that v(0) = 1   
  // returns (H_v x)(0) = +/- ||x||
  extern NGS_DLL_HEADER double CalcHouseholderVectorInPlace (SliceVector<> x);

  
  // find Householder reflection vector v such that
  // reflection matrix H_v = I - 2 v v^T / (v^T v)
  // leads to H_v x = +/- e_0
  // scaling such that v(0) = 1   
  inline void CalcHouseholderVector (SliceVector<> x, FlatVector<> v)
  {
    v = x;
    CalcHouseholderVectorInPlace (v);
  }


  
  class NGS_DLL_HEADER HouseholderReflection
  {
    FlatVector<> v;
    double factor;   // 2 / (v^T v)
  public:
    HouseholderReflection (FlatVector<> av);

    template <ORDERING ORD>
    void TMult (SliceMatrix<double,ORD> m2) const; 
    void Mult (SliceMatrix<double,RowMajor> m2) const { TMult(m2); }
    void Mult (SliceMatrix<double,ColMajor> m2) const { TMult(m2); }
  };

  // H = H_{m-1} ... H_1 H_0 = I - V^T T V
  template <ORDERING OMV>
  class MultiHouseholderReflection
  {
    SliceMatrix<double, OMV> mv;  // every row one reflection vector
    Matrix<> T;        // strict lower triangular
  public:
    MultiHouseholderReflection (SliceMatrix<double, OMV> amv);
    template <ORDERING ORD>
    void TMult (SliceMatrix<double,ORD> m2) const;  // Hm-1 * ... * H1 * H0 * m2
    void Mult (SliceMatrix<double, RowMajor> m2) const { TMult (m2); } 
    void Mult (SliceMatrix<double, ColMajor> m2) const { TMult (m2); }

    template <ORDERING ORD>
    void TMultTrans (SliceMatrix<double,ORD> m2) const; 
    void MultTrans (SliceMatrix<double,RowMajor> m2) const { TMultTrans(m2); }
    void MultTrans (SliceMatrix<double,ColMajor> m2) const { TMultTrans(m2); }
  };


  // H is a generalized, normalized, lower left triangular matrix
  // its columns are the Householder vectors
  template <ORDERING OH, ORDERING OM>
  void ApplyHouseholderReflections (SliceMatrix<double,OH> H, SliceMatrix<double,OM> M)
  {
    size_t bs = 48;
    size_t n = H.Width();
    size_t m = H.Height();

    if (n > m) throw Exception("ApplyHouseholderReflections, n > m");
    if (n == 0) return;
    
    for (size_t i = 0; i < n; i += bs)
      {
        size_t bsi = min(n-i, bs);
        MultiHouseholderReflection Hv(Trans(H.Cols(i,i+bsi).Rows(i, m)));
        Hv.Mult (M.Rows(i,m));
      }
  }

  // H is a generalized, normalized, lower left triangular matrix
  // its columns are the Householder vectors
  // go backwards
  template <ORDERING OH, ORDERING OM>
  void ApplyHouseholderReflectionsTrans (SliceMatrix<double,OH> H, SliceMatrix<double,OM> M)
  {
    size_t bs = 48;
    size_t n = H.Width();
    size_t m = H.Height();

    if (n > m) throw Exception("ApplyHouseholderReflections, n > m");
    if (n == 0) return;
    
    // for (size_t i = 0; i < n; i += bs)
    for (size_t i = n; i+bs > bs; i -= bs)
      {
        size_t first = max(bs, i)-bs;
        size_t last = i;
        MultiHouseholderReflection Hv(Trans(H.Cols(first, last).Rows(first, m)));
        Hv.MultTrans (M.Rows(first,m));
      }
  }

  
  
  /*
    Factorize A = Q R
    A is m x n, overwritten by R
    Q is either m x m  or m x n
  */
  extern NGS_DLL_HEADER void QRFactorization (SliceMatrix<> A, SliceMatrix<> Q);

  extern NGS_DLL_HEADER void QRFactorizationInPlace (SliceMatrix<double,RowMajor> A);
  extern NGS_DLL_HEADER void QRFactorizationInPlace (SliceMatrix<double,ColMajor> A);
  
  extern NGS_DLL_HEADER void CalcSVD (SliceMatrix<double, RowMajor> A, 
                                      SliceMatrix<double, ColMajor> U,
                                      SliceMatrix<double, ColMajor> V);
  
  extern NGS_DLL_HEADER void CalcSVD (SliceMatrix<double, ColMajor> A, 
                                      SliceMatrix<double, ColMajor> U,
                                      SliceMatrix<double, ColMajor> V);
}

#endif
