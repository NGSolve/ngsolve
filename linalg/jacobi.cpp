 /**************************************************************************/
/* File:   blockjacobi.cpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   14. Aug. 2002                                                  */
/**************************************************************************/


#include <la.hpp>

namespace ngla
{
  using namespace ngla;

  template <class TM, class TV_ROW, class TV_COL>
  JacobiPrecond<TM,TV_ROW,TV_COL> ::
  JacobiPrecond (const SparseMatrix<TM,TV_ROW,TV_COL> & amat, 
		 const BitArray * ainner)
    : mat(amat), inner(ainner)
  { 
    height = mat.Height();
    invdiag = new TM[height];
    for (int i = 0; i < height; i++)
      ngbla::CalcInverse (mat(i,i), invdiag[i]);

    //(*testout) << "inv diag = " << endl;
    //for (int i = 0; i < height; i++)
    //  (*testout) << "inv[" << i << "] = " << invdiag[i] << endl;
  }

  ///
  template <class TM, class TV_ROW, class TV_COL>
  JacobiPrecond<TM,TV_ROW,TV_COL> :: ~JacobiPrecond () 
  {
    // he: this is an array
    delete [] invdiag; 
  }

  ///
  template <class TM, class TV_ROW, class TV_COL>
  void JacobiPrecond<TM,TV_ROW,TV_COL> ::
  MultAdd (TSCAL s, const BaseVector & x, BaseVector & y) const 
  {
    const FlatVector<TV_ROW> fx = 
      dynamic_cast<const T_BaseVector<TV_ROW> &> (x).FV();
    FlatVector<TV_ROW> fy = 
      dynamic_cast<T_BaseVector<TV_ROW> &> (y).FV();

    if (!inner)
      for (int i = 0; i < height; i++)
	fy(i) += s * (invdiag[i] * fx(i));
    else
      for (int i = 0; i < height; i++)
	if (inner->Test(i))
	  fy(i) += s * (invdiag[i] * fx(i));
  }


  ///
  template <class TM, class TV_ROW, class TV_COL>
  BaseVector * JacobiPrecond<TM,TV_ROW,TV_COL> :: CreateVector () const 
  {
    return mat.CreateVector();
  }


  ///
  template <class TM, class TV_ROW, class TV_COL>
  void JacobiPrecond<TM,TV_ROW,TV_COL> ::
  GSSmooth (BaseVector & x, const BaseVector & b) const 
  {
    static int timer = NgProfiler::CreateTimer ("JacobiPrecond::GSSmooth");
    NgProfiler::RegionTimer reg (timer);

    FlatVector<TV_ROW> fx = 
      dynamic_cast<T_BaseVector<TV_ROW> &> (x).FV();
    const FlatVector<TV_ROW> fb = 
      dynamic_cast<const T_BaseVector<TV_ROW> &> (b).FV();

    for (int i = 0; i < height; i++)
      if (!this->inner || this->inner->Test(i))
	{
	  TV_ROW ax = mat.RowTimesVector (i, fx);
	  fx(i) += invdiag[i] * (fb(i) - ax);
	}
  }


  ///
  template <class TM, class TV_ROW, class TV_COL>
  void JacobiPrecond<TM,TV_ROW,TV_COL> ::
  GSSmoothBack (BaseVector & x, const BaseVector & b) const 
  {
    static int timer = NgProfiler::CreateTimer ("JacobiPrecond::GSSmoothBack");
    NgProfiler::RegionTimer reg (timer);

    FlatVector<TV_ROW> fx = 
      dynamic_cast<T_BaseVector<TV_ROW> &> (x).FV();
    const FlatVector<TV_ROW> fb = 
      dynamic_cast<const T_BaseVector<TV_ROW> &> (b).FV();

    for (int i = height-1; i >= 0; i--)
      if (!this->inner || this->inner->Test(i))
	{
	  TV_ROW ax = mat.RowTimesVector (i, fx);
	  fx(i) += invdiag[i] * (fb(i) - ax);
	}
  }

  ///
  template <class TM, class TV_ROW, class TV_COL>
  void JacobiPrecond<TM,TV_ROW,TV_COL> ::
  GSSmoothNumbering (BaseVector & x, const BaseVector & b,
		     const Array<int> & numbering, 
		     int forward) const
  {
    ;
  }




  template <class TM, class TV>
  JacobiPrecondSymmetric<TM,TV> ::
  JacobiPrecondSymmetric (const SparseMatrixSymmetric<TM,TV> & amat, 
			  const BitArray * ainner)
    : JacobiPrecond<TM,TV,TV> (amat, ainner)
  { 
    ;
  }

  ///
  template <class TM, class TV>
  void JacobiPrecondSymmetric<TM,TV> ::
  GSSmooth (BaseVector & x, const BaseVector & b) const 
  {
    static int timer = NgProfiler::CreateTimer ("JacobiPrecondSymmetric::GSSmooth");
    NgProfiler::RegionTimer reg (timer);

    FlatVector<TVX> fx = 
      dynamic_cast<T_BaseVector<TVX> &> (x).FV();
    const FlatVector<TVX> fb = 
      dynamic_cast<const T_BaseVector<TVX> &> (b).FV();

    const SparseMatrixSymmetric<TM,TV> & smat =
      dynamic_cast<const SparseMatrixSymmetric<TM,TV>&> (this->mat);

    // x := b - L^t x
    for (int i = 0; i < this->height; i++)
      if (!this->inner || this->inner->Test(i))
	{
	  smat.AddRowTransToVectorNoDiag (i, -fx(i), fx);
	  fx(i) = fb(i);
	}
      else
	fx(i) = TVX(0);
    
    // x := (L+D)^{-1} x
    for (int i = 0; i < this->height; i++)
      if (!this->inner || this->inner->Test(i))
	{
	  TVX hv = fx(i) - smat.RowTimesVectorNoDiag (i, fx);
	  fx(i) = this->invdiag[i] * hv;
	}
  }



  template <class TM, class TV>
  void JacobiPrecondSymmetric<TM,TV> ::
  GSSmooth (BaseVector & x, const BaseVector & b, BaseVector & y /* , BaseVector & help */) const 
  {
    static int timer = NgProfiler::CreateTimer ("JacobiPrecondSymmetric::GSSmooth-help");
    NgProfiler::RegionTimer reg (timer);

    FlatVector<TVX> fx = 
      dynamic_cast<T_BaseVector<TVX> &> (x).FV();
    FlatVector<TVX> fy = 
      dynamic_cast<T_BaseVector<TVX> &> (y).FV();
    const FlatVector<TVX> fb = 
      dynamic_cast<const T_BaseVector<TVX> &> (b).FV();
    // FlatVector<TVX> fh = 
    // dynamic_cast<const T_BaseVector<TVX> &> (help).FV();

    const SparseMatrixSymmetric<TM,TV> & smat =
      dynamic_cast<const SparseMatrixSymmetric<TM,TV>&> (this->mat);

    // input, y = b - (D+L^t) x
    // (L+D) x_new := b - L^t x  = y + D x
    // D (x_new-x) = b - L x_new
    // y -= (D+L^t) w
    for (int i = 0; i < this->height; i++)
      if (!this->inner || this->inner->Test(i))
	{
	  TVX d = fy(i) - smat.RowTimesVectorNoDiag (i, fx);
	  TVX w = this->invdiag[i] * d;
	  
	  fx(i) += w;
	  smat.AddRowTransToVector (i, -w, fy);
	}
    // else
    // fx(i) = TVX(0);
  }






  ///
  template <class TM, class TV>
  void JacobiPrecondSymmetric<TM,TV> ::
  GSSmoothBack (BaseVector & x, const BaseVector & b) const 
  {
    static int timer = NgProfiler::CreateTimer ("JacobiPrecondSymmetric::GSSmoothBack");
    NgProfiler::RegionTimer reg (timer);

    FlatVector<TVX> fx = 
      dynamic_cast<T_BaseVector<TVX> &> (x).FV();
    const FlatVector<TVX> fb = 
      dynamic_cast<const T_BaseVector<TVX> &> (b).FV();

    const SparseMatrixSymmetric<TM,TV> & smat =
      dynamic_cast<const SparseMatrixSymmetric<TM,TV>&> (this->mat);
    
    for (int i = this->height-1; i >= 0; i--)
      if (!this->inner || this->inner->Test(i))
	{
	  fx(i) = fb(i) - smat.RowTimesVectorNoDiag (i, fx);
	}
      else
	fx(i) = TVX(0);

    
    for (int i = this->height-1; i >= 0; i--)
      if (!this->inner || this->inner->Test(i))
	{
	  TVX val = this->invdiag[i] * fx(i);
	  fx(i) = val;
	  smat.AddRowTransToVectorNoDiag (i, -val, fx);
	}	 
      else
	fx(i) = TVX(0);
  }

  ///
  template <class TM, class TV>
  void JacobiPrecondSymmetric<TM,TV> ::
  GSSmoothNumbering (BaseVector & x, const BaseVector & b,
		     const Array<int> & numbering, 
		     int forward) const
  {
    ;
  }



  template class JacobiPrecond<double>;
  template class JacobiPrecond<Complex>;
  template class JacobiPrecond<double, Complex, Complex>;
#if MAX_SYS_DIM >= 1
  template class JacobiPrecond<Mat<1,1,double> >;
  template class JacobiPrecond<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class JacobiPrecond<Mat<2,2,double> >;
  template class JacobiPrecond<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class JacobiPrecond<Mat<3,3,double> >;
  template class JacobiPrecond<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class JacobiPrecond<Mat<4,4,double> >;
  template class JacobiPrecond<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class JacobiPrecond<Mat<5,5,double> >;
  template class JacobiPrecond<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class JacobiPrecond<Mat<6,6,double> >;
  template class JacobiPrecond<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class JacobiPrecond<Mat<7,7,double> >;
  template class JacobiPrecond<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class JacobiPrecond<Mat<8,8,double> >;
  template class JacobiPrecond<Mat<8,8,Complex> >;
#endif


 
  template class JacobiPrecondSymmetric<double>;
  template class JacobiPrecondSymmetric<Complex>;
  template class JacobiPrecondSymmetric<double,Complex>;
#if MAX_SYS_DIM >= 1
  template class JacobiPrecondSymmetric<Mat<1,1,double> >;
  template class JacobiPrecondSymmetric<Mat<1,1,Complex> >;
#endif
#if MAX_SYS_DIM >= 2
  template class JacobiPrecondSymmetric<Mat<2,2,double> >;
  template class JacobiPrecondSymmetric<Mat<2,2,Complex> >;
#endif
#if MAX_SYS_DIM >= 3
  template class JacobiPrecondSymmetric<Mat<3,3,double> >;
  template class JacobiPrecondSymmetric<Mat<3,3,Complex> >;
#endif
#if MAX_SYS_DIM >= 4
  template class JacobiPrecondSymmetric<Mat<4,4,double> >;
  template class JacobiPrecondSymmetric<Mat<4,4,Complex> >;
#endif
#if MAX_SYS_DIM >= 5
  template class JacobiPrecondSymmetric<Mat<5,5,double> >;
  template class JacobiPrecondSymmetric<Mat<5,5,Complex> >;
#endif
#if MAX_SYS_DIM >= 6
  template class JacobiPrecondSymmetric<Mat<6,6,double> >;
  template class JacobiPrecondSymmetric<Mat<6,6,Complex> >;
#endif
#if MAX_SYS_DIM >= 7
  template class JacobiPrecondSymmetric<Mat<7,7,double> >;
  template class JacobiPrecondSymmetric<Mat<7,7,Complex> >;
#endif
#if MAX_SYS_DIM >= 8
  template class JacobiPrecondSymmetric<Mat<8,8,double> >;
  template class JacobiPrecondSymmetric<Mat<8,8,Complex> >;
#endif

#ifdef CACHEBLOCKSIZE
  template class JacobiPrecond<double, Vec<CACHEBLOCKSIZE>, Vec<CACHEBLOCKSIZE> >;
  template class JacobiPrecondSymmetric<double, Vec<CACHEBLOCKSIZE> >;
#endif


#if MAX_CACHEBLOCKS >= 2
  template class JacobiPrecond<double, Vec<2,double>, Vec<2,double> >;
  template class JacobiPrecondSymmetric<double, Vec<2,double> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class JacobiPrecond<double, Vec<3,double>, Vec<3,double> >;
  template class JacobiPrecondSymmetric<double, Vec<3,double> >;
  template class JacobiPrecond<double, Vec<4,double>, Vec<4,double> >;
  template class JacobiPrecondSymmetric<double, Vec<4,double> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class JacobiPrecond<double, Vec<5,double>, Vec<5,double> >;
  template class JacobiPrecondSymmetric<double, Vec<5,double> >;
  template class JacobiPrecond<double, Vec<6,double>, Vec<6,double> >;
  template class JacobiPrecondSymmetric<double, Vec<6,double> >;
  template class JacobiPrecond<double, Vec<7,double>, Vec<7,double> >;
  template class JacobiPrecondSymmetric<double, Vec<7,double> >;
  template class JacobiPrecond<double, Vec<8,double>, Vec<8,double> >;
  template class JacobiPrecondSymmetric<double, Vec<8,double> >;
  template class JacobiPrecond<double, Vec<9,double>, Vec<9,double> >;
  template class JacobiPrecondSymmetric<double, Vec<9,double> >;
  template class JacobiPrecond<double, Vec<10,double>, Vec<10,double> >;
  template class JacobiPrecondSymmetric<double, Vec<10,double> >;
  template class JacobiPrecond<double, Vec<11,double>, Vec<11,double> >;
  template class JacobiPrecondSymmetric<double, Vec<11,double> >;
  template class JacobiPrecond<double, Vec<12,double>, Vec<12,double> >;
  template class JacobiPrecondSymmetric<double, Vec<12,double> >;
  template class JacobiPrecond<double, Vec<13,double>, Vec<13,double> >;
  template class JacobiPrecondSymmetric<double, Vec<13,double> >;
  template class JacobiPrecond<double, Vec<14,double>, Vec<14,double> >;
  template class JacobiPrecondSymmetric<double, Vec<14,double> >;
  template class JacobiPrecond<double, Vec<15,double>, Vec<15,double> >;
  template class JacobiPrecondSymmetric<double, Vec<15,double> >;
#endif
#if MAX_CACHEBLOCKS >= 2
  template class JacobiPrecond<double, Vec<2,Complex>, Vec<2,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class JacobiPrecond<double, Vec<3,Complex>, Vec<3,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<3,Complex> >;
  template class JacobiPrecond<double, Vec<4,Complex>, Vec<4,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class JacobiPrecond<double, Vec<5,Complex>, Vec<5,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<5,Complex> >;
  template class JacobiPrecond<double, Vec<6,Complex>, Vec<6,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<6,Complex> >;
  template class JacobiPrecond<double, Vec<7,Complex>, Vec<7,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<7,Complex> >;
  template class JacobiPrecond<double, Vec<8,Complex>, Vec<8,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<8,Complex> >;
  template class JacobiPrecond<double, Vec<9,Complex>, Vec<9,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<9,Complex> >;
  template class JacobiPrecond<double, Vec<10,Complex>, Vec<10,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<10,Complex> >;
  template class JacobiPrecond<double, Vec<11,Complex>, Vec<11,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<11,Complex> >;
  template class JacobiPrecond<double, Vec<12,Complex>, Vec<12,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<12,Complex> >;
  template class JacobiPrecond<double, Vec<13,Complex>, Vec<13,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<13,Complex> >;
  template class JacobiPrecond<double, Vec<14,Complex>, Vec<14,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<14,Complex> >;
  template class JacobiPrecond<double, Vec<15,Complex>, Vec<15,Complex> >;
  template class JacobiPrecondSymmetric<double, Vec<15,Complex> >;
#endif

#if MAX_CACHEBLOCKS >= 2
  template class JacobiPrecond<Complex, Vec<2,Complex>, Vec<2,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<2,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 3
  template class JacobiPrecond<Complex, Vec<3,Complex>, Vec<3,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<3,Complex> >;
  template class JacobiPrecond<Complex, Vec<4,Complex>, Vec<4,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<4,Complex> >;
#endif
#if MAX_CACHEBLOCKS >= 5
  template class JacobiPrecond<Complex, Vec<5,Complex>, Vec<5,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<5,Complex> >;
  template class JacobiPrecond<Complex, Vec<6,Complex>, Vec<6,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<6,Complex> >;
  template class JacobiPrecond<Complex, Vec<7,Complex>, Vec<7,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<7,Complex> >;
  template class JacobiPrecond<Complex, Vec<8,Complex>, Vec<8,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<8,Complex> >;
  template class JacobiPrecond<Complex, Vec<9,Complex>, Vec<9,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<9,Complex> >;
  template class JacobiPrecond<Complex, Vec<10,Complex>, Vec<10,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<10,Complex> >;
  template class JacobiPrecond<Complex, Vec<11,Complex>, Vec<11,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<11,Complex> >;
  template class JacobiPrecond<Complex, Vec<12,Complex>, Vec<12,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<12,Complex> >;
  template class JacobiPrecond<Complex, Vec<13,Complex>, Vec<13,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<13,Complex> >;
  template class JacobiPrecond<Complex, Vec<14,Complex>, Vec<14,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<14,Complex> >;
  template class JacobiPrecond<Complex, Vec<15,Complex>, Vec<15,Complex> >;
  template class JacobiPrecondSymmetric<Complex, Vec<15,Complex> >;
#endif
}
