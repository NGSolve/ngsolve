#ifndef CALCINVERSE_HPP
#define CALCINVERSE_HPP


namespace ngbla
{

  /* **************************** Inverse *************************** */


  enum class INVERSE_LIB { INV_NGBLA, INV_NGBLA_LU, INV_LAPACK, INV_NGBLA_QR, INV_CHOOSE };

  /// Calculate inverse. Gauss elimination with row pivoting
  template <class T2>
  extern NGS_DLL_HEADER void CalcInverse (FlatMatrix<T2> inv, INVERSE_LIB il = INVERSE_LIB::INV_CHOOSE);

  extern NGS_DLL_HEADER void CalcInverse (FlatMatrix<double> inv, INVERSE_LIB il = INVERSE_LIB::INV_CHOOSE);

  template <class T, class T2>
  inline void CalcInverse (const FlatMatrix<T> m, FlatMatrix<T2> inv)
  {
    inv = m;
    CalcInverse (inv);
  }

  template <class T, class T2>
  inline void CalcInverse (const FlatMatrix<T> m, Matrix<T2> & inv)
  {
    inv = m;
    CalcInverse (inv);
  }



  /**
     Calculates the inverse of a Matrix.
  */

  template <typename T>
  inline Matrix<T> Inverse (const FlatMatrix<T> & m)
  {
    Matrix<T> inv(m.Height(),m.Height());
    CalcInverse (m, inv);
    return inv;
  }

  
  inline void CalcInverse (double & m)
  {
    m = 1 / m;
  }

  inline void CalcInverse (Complex & m)
  {
    m = 1.0 / m;
  }

  template <int H, int W, typename T>
  inline void CalcInverse (Mat<H,W,T> & m)
  {
    FlatMatrix<T> fm(m);
    CalcInverse (fm);
  }


  INLINE void CalcInverse (const double & m, double & inv)
  {
    inv = 1 / m;
  }

  INLINE void CalcInverse (const Complex & m, Complex & inv)
  {
    inv = 1.0 / m;
  }

  template <int H, int W, typename T, typename TINV>
  inline void CalcInverse (const Mat<H,W,T> & m, TINV & inv)
  {
    FlatMatrix<T> fm(m);
    FlatMatrix<T> finv(inv);
    CalcInverse (fm, finv);
  }

  template <typename T, typename TINV>
  INLINE void CalcInverse (const Mat<0,0,T> & m, TINV & inv)
  {
    ;
  }

  template <typename T, typename TINV>
  INLINE void CalcInverse (const Mat<1,1,T> & m, TINV & inv)
  {
    inv(0,0) = 1.0 / m(0,0);
  }

  template <typename T, typename TINV>
  INLINE void CalcInverse (const Mat<2,2,T> & m, TINV & inv)
  {
    T idet = 1.0 / (m(0,0) * m(1,1) - m(0,1) * m(1,0));
    inv(0,0) = idet * m(1,1);
    inv(0,1) = -idet * m(0,1);
    inv(1,0) = -idet * m(1,0);
    inv(1,1) = idet * m(0,0);
  }

  template <typename T, typename TINV>
  INLINE void CalcInverse (Mat<3,3,T> m, TINV & inv)
  {
    T h0 = m(4)*m(8)-m(5)*m(7);
    T h1 = m(5)*m(6)-m(3)*m(8);
    T h2 = m(3)*m(7)-m(4)*m(6);
    T det = m(0) * h0 + m(1) * h1 + m(2) * h2;
    T idet = 1.0 / det;
    
    inv(0,0) =  idet * h0; 
    inv(0,1) = -idet * (m(1) * m(8) - m(2) * m(7));
    inv(0,2) =  idet * (m(1) * m(5) - m(2) * m(4));
    
    inv(1,0) =  idet * h1; 
    inv(1,1) =  idet * (m(0) * m(8) - m(2) * m(6));
    inv(1,2) = -idet * (m(0) * m(5) - m(2) * m(3));
    
    inv(2,0) =  idet * h2; 
    inv(2,1) = -idet * (m(0) * m(7) - m(1) * m(6));
    inv(2,2) =  idet * (m(0) * m(4) - m(1) * m(3));
    return;
  }


  
  /// Computes the Schur Complement.
  extern NGS_DLL_HEADER void CalcSchurComplement (const FlatMatrix<double> a, 
				   FlatMatrix<double> s,
				   const BitArray & used,
				   LocalHeap & lh);

}



#endif
