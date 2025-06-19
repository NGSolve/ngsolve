#ifndef FILE_JACOBI
#define FILE_JACOBI

/* *************************************************************************/
/* File:   jacobi.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   06. Oct. 96                                                    */
/* *************************************************************************/


#include "sparsematrix.hpp"

namespace ngla
{

  /**
     Jacobi and Gauss Seidel smoother
     for scalar, block and system matrices
  */
  
  class BaseMSMPrecond : virtual public BaseMatrix
  {
  public:
    virtual void Smooth (BaseVector & x, const BaseVector & b, int steps = 1) const = 0;
    virtual void SmoothBack (BaseVector & x, const BaseVector & b, int steps = 1) const = 0;
  };

  
  class BaseJacobiPrecond : public BaseMSMPrecond
  {
  public:
    virtual void Smooth (BaseVector & x, const BaseVector & b, int steps = 1) const override
    {
      GSSmooth (x, b, steps);
    }
    virtual void SmoothBack (BaseVector & x, const BaseVector & b, int steps = 1) const override
    {
      GSSmoothBack (x, b, steps);
    }
    
    virtual void GSSmooth (BaseVector & x, const BaseVector & b, int steps = 1) const = 0;
    virtual void GSSmooth (BaseVector & x, const BaseVector & b, BaseVector & y) const = 0;
    virtual void GSSmoothBack (BaseVector & x, const BaseVector & b, int steps = 1) const = 0;
  };


  class SymmetricGaussSeidelPrecond : virtual public BaseMatrix
  {
    shared_ptr<BaseJacobiPrecond> jac;
  public:
    SymmetricGaussSeidelPrecond (shared_ptr<BaseSparseMatrix> mat, shared_ptr<BitArray> freedofs);
    int VHeight() const override { return jac->VHeight(); }
    int VWidth() const override { return jac->VHeight(); }

    void Mult (const BaseVector & x, BaseVector & y) const override;
    
    AutoVector CreateRowVector() const override { return jac->CreateRowVector(); }
    AutoVector CreateColVector() const override { return jac->CreateColVector(); }
  };

  
  /// A Jaboci preconditioner for general sparse matrices
  template <class TM, class TV_ROW, class TV_COL>
  class JacobiPrecond : public BaseJacobiPrecond,
			public S_BaseMatrix<typename mat_traits<TM>::TSCAL>
  {
  protected:
    shared_ptr<SparseMatrix<TM,TV_ROW,TV_COL>> mat;
    ///
    shared_ptr<BitArray> inner;
    ///
    int height;
    ///
    Array<TM> invdiag;
  public:
    // typedef typename mat_traits<TM>::TV_ROW TVX;
    typedef typename mat_traits<TM>::TSCAL TSCAL;

    ///
    JacobiPrecond (shared_ptr<SparseMatrix<TM,TV_ROW,TV_COL>> amat, 
		   shared_ptr<BitArray> ainner = nullptr, bool use_par = true);

    int VHeight() const override { return height; }
    int VWidth() const override { return height; }

    FlatArray<TM> GetInverse() const { return invdiag; }
    
    ///
    void MultAdd (TSCAL s, const BaseVector & x, BaseVector & y) const override;

    void MultTransAdd (TSCAL s, const BaseVector & x, BaseVector & y) const override
    { MultAdd (s, x, y); }
    ///
    AutoVector CreateRowVector() const override { return mat->CreateColVector(); }
    AutoVector CreateColVector() const override { return mat->CreateRowVector(); }
    ///
    void GSSmooth (BaseVector & x, const BaseVector & b, int steps) const override;

    /// computes partial residual y
    void GSSmooth (BaseVector & x, const BaseVector & b, BaseVector & y) const override
    {
      GSSmooth (x, b, 1);
    }


    ///
    void GSSmoothBack (BaseVector & x, const BaseVector & b, int steps) const override;

    ///
    virtual void GSSmoothNumbering (BaseVector & x, const BaseVector & b,
				    const Array<int> & numbering, 
				    int forward = 1) const;
  };









  /// A Jaboci preconditioner for symmetric sparse matrices
  template <class TM, class TV>
  class NGS_DLL_HEADER JacobiPrecondSymmetric : public JacobiPrecond<TM,TV,TV>
  {
  public:
    typedef TV TVX;

    ///
    JacobiPrecondSymmetric (shared_ptr<SparseMatrixSymmetric<TM,TV>> amat, 
			    shared_ptr<BitArray> ainner = nullptr, bool use_par = true);

    ///
    virtual void GSSmooth (BaseVector & x, const BaseVector & b, int steps) const override;

    /// computes partial residual y
    virtual void GSSmooth (BaseVector & x, const BaseVector & b, BaseVector & y /* , BaseVector & help */) const override;

    ///
    virtual void GSSmoothBack (BaseVector & x, const BaseVector & b, int steps) const override;
    virtual void GSSmoothBack (BaseVector & x, const BaseVector & b, BaseVector & y) const;

    ///
    virtual void GSSmoothNumbering (BaseVector & x, const BaseVector & b,
				    const Array<int> & numbering, 
				    int forward = 1) const override;
  };

}


#endif
