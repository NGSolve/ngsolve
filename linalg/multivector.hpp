#ifndef FILE_MULTIVECTOR
#define FILE_MULTIVECTOR

/*********************************************************************/
/* File:   multivector.hpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2010                                                 */
/*********************************************************************/

namespace ngla {

  /* TODO: 
   * (+) Gram-Schmidt
   * (+) default arguments 
   *
  */

  class MultiVector;  

  class MultiVectorExpr
  {
  public:
    virtual ~MultiVectorExpr() { ; }
    virtual void AssignTo (double s, class MultiVector & v) const = 0;
    virtual void AddTo (double s, class MultiVector & v) const = 0;
    virtual void AssignTo (Complex s, class MultiVector & v) const = 0;
    virtual void AddTo (Complex s, class MultiVector & v) const = 0;
  };
  
  class MultiVector
  {
    shared_ptr<BaseVector> refvec;
    Array<shared_ptr<BaseVector>> vecs;
  public:
    MultiVector (shared_ptr<BaseVector> v, size_t cnt): refvec (v)
    {
      Expand (cnt);
    }
    MultiVector (size_t size, size_t cnt, bool is_complex)
    {
      refvec = CreateBaseVector(size, is_complex, 1);
      Expand (cnt);
    }
  
    bool IsComplex() const { return refvec->IsComplex(); }

    size_t Size() const { return vecs.Size(); }
    shared_ptr<BaseVector> operator[] (size_t i) const { return vecs[i]; }
    shared_ptr<BaseVector> RefVec() const { return refvec; }
    void Expand (size_t nr = 1) {
      for ([[maybe_unused]] auto i : Range(nr))
        vecs.Append (refvec->CreateVector());
    }
    void Append (shared_ptr<BaseVector> v)
    {
      vecs.Append (v->CreateVector());
      *vecs.Last() = *v;
    }


    MultiVector & operator= (double val)
    {
      for (auto & vec : vecs)
        *vec = val;
      return *this;
    }

    MultiVector & operator= (Complex val)
    {
      for (auto & vec : vecs) 
        *vec = val;
      return *this;      
    }
      
      MultiVector & operator= (const MultiVectorExpr & expr)
    { 
      expr.AssignTo(1, *this);
      return *this;
    }
    
    void operator+= (const MultiVectorExpr & expr)
    { 
      expr.AddTo(1, *this); 
    }
    
    void operator-= (const MultiVectorExpr & expr)
    { 
      expr.AddTo(-1, *this); 
    }
  };
  

  template <class T>
  void MultAdd (const class BaseMatrix & mat, T s, const MultiVector & x, MultiVector & y);
  
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

    void AssignTo (Complex s, MultiVector & res) const override
    {
      res = 0.0;
      MultAdd (*mat, s, *vec, res);
    }

    void AddTo (Complex s, MultiVector & res) const override
    { MultAdd (*mat, s, *vec, res); }
  };

  
  // y += sum a(i) * x[i]
  template <class T>
  void Axpy (const Vector<T> & a, const MultiVector  & x, BaseVector & y);


  // TODO: do that better!!
  template <typename TSCAL>
  class MultiVecAxpyExpr : public DynamicBaseExpression
  {
    Vector<TSCAL> a;
    shared_ptr<MultiVector> x;

  public:
    MultiVecAxpyExpr (Vector<TSCAL> aa, shared_ptr<MultiVector> ax)
      : a(aa), x(ax) { }

    void AssignTo (double s, BaseVector & v) const override
    {
      v = 0.0;
      AddTo (s, v);
    }
    
    void AddTo (double s, BaseVector & v) const override
    {
      Vector<TSCAL> sa = s*a;
      Axpy (sa, *x, v);
    }

    // these are not called for some reason
    void AssignTo (Complex s, BaseVector & v) const override
    {
      v = 0.0;
      AddTo (s, v);
    }

    void AddTo (Complex s, BaseVector & v) const override
    {
      Vector<Complex> sa = s*a;
      Axpy(sa, *x , v);
    }

  };

  // template <class T>
  // Matrix<T> InnerProduct (const MultiVector & x, const MultiVector & y, bool conjugate);
}
#endif
