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
  // reflection matrix H_v = I - v v^T
  // leads to H_v x = +/- e_0
  extern NGS_DLL_HEADER void CalcHouseholderVector (SliceVector<> x, FlatVector<> v);

  class NGS_DLL_HEADER HouseholderReflection
  {
    FlatVector<> v;
  public:
    HouseholderReflection (FlatVector<> av) : v(av) { ; }
  
    void Mult (SliceMatrix<> m2) const;
    void Mult (SliceMatrix<double,ColMajor> m2) const;
  };

  // H = H_{m-1} ... H_1 H_0 = I - V^T T V
  class MultiHouseholderReflection
  {
    SliceMatrix<> mv;  // every row one reflection vector
    Matrix<> T;        // strict lower triangular
  public:
    MultiHouseholderReflection (SliceMatrix<> amv);
    void Mult (SliceMatrix<> m2) const;  // Hm-1 * ... * H1 * H0 * m2
    void Mult (SliceMatrix<double,ColMajor> m2) const;   // Hm-1 * ... * H1 * H0 * m2
  };


  /*
    Factorize A = Q R
    A is m x n, overwritten by R
    Q is either m x m      (TODO: or m x n)
  */
  extern NGS_DLL_HEADER void QRFactorization (SliceMatrix<> A, SliceMatrix<> Q);
}

#endif
