#ifndef FILE_MULTIVECTOR
#define FILE_MULTIVECTOR

/*********************************************************************/
/* File:   multivector.hpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2010                                                 */
/*********************************************************************/

namespace ngla {

  class MultiVectorExpr
  {
  public:
    virtual ~MultiVectorExpr() { ; }
    virtual void AssignTo (double s, class MultiVector & v) const = 0;
    virtual void AddTo (double s, class MultiVector & v) const = 0;
  };
  
  class MultiVector
  {
    shared_ptr<BaseVector> refvec;
    Array<shared_ptr<BaseVector>> vecs;
  public:
    MultiVector (shared_ptr<BaseVector> v, size_t cnt);
    size_t Size() const { return vecs.Size(); }
    shared_ptr<BaseVector> operator[] (size_t i) const { return vecs[i]; }
    void Expand (size_t nr = 1);
    void Append (shared_ptr<BaseVector> v); 
    void operator= (double v);
    template <typename T>
    void operator= (const MultiVectorExpr & expr)
    { expr.AssignTo(1, *this); }
    void operator+= (const MultiVectorExpr & expr)
    { expr.AddTo(1, *this); }
    void operator-= (const MultiVectorExpr & expr)
    { expr.AddTo(-1, *this); }
  };
  
  void MultAdd (const class BaseMatrix & mat, double s, const MultiVector & x, MultiVector & y);
  
  class MatMultiVecExpr : public MultiVectorExpr
  {
    shared_ptr<BaseMatrix> mat;
    shared_ptr<MultiVector> vec;
  public:
    MatMultiVecExpr (shared_ptr<BaseMatrix> amat, shared_ptr<MultiVector> avec)
      : mat(amat), vec(avec) { ; }
    void AssignTo (double s, MultiVector & res) const override
    {
      res = 0.0;
      MultAdd (*mat, s, *vec, res);
    }
    void AddTo (double s, MultiVector & res) const override
    { MultAdd (*mat, s, *vec, res); }
  };

  
  // y += sum a(i) * x[i]
  void Axpy (const Vector<> & a, const MultiVector & x, BaseVector & y);
  class MultiVecAxpyExpr : public DynamicBaseExpression
  {
    Vector<> a;
    shared_ptr<MultiVector> x;
  public:
    MultiVecAxpyExpr (Vector<> aa, shared_ptr<MultiVector> ax) : a(aa), x(ax) { ; }
    
    void AssignTo (double s, BaseVector & v) const override
    {
      v = 0.0;
      AddTo (s, v);
    }
    void AddTo (double s, BaseVector & v) const override
    {
      Vector<> sa = s*a;
      Axpy(sa, *x, v);
    }
    void AssignTo (Complex s, BaseVector & v2) const override { } // missing
    void AddTo (Complex s, BaseVector & v2) const  override { } // missing
  };

  
  Matrix<> InnerProduct (const MultiVector & x, const MultiVector & y);
  
}
#endif

