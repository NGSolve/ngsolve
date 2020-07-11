#ifndef FILE_MULTIVECTOR
#define FILE_MULTIVECTOR

/*********************************************************************/
/* File:   multivector.hpp                                           */
/* Author: Joachim Schoeberl, Amanda Schoefl                         */
/* Date:   June 2020                                                 */
/*********************************************************************/

namespace ngla {

  /* TODO: 
   * (+) Gram-Schmidt
   * (+) default arguments 
   *
  */

  class MultiVector;  
  class BaseMatrix;

  
  class MultiVectorExpr
  {
  public:
    virtual ~MultiVectorExpr() { ; }
    virtual void AssignTo (double s, class MultiVector & v) const = 0;
    virtual void AddTo (double s, class MultiVector & v) const = 0;
    virtual void AssignTo (Complex s, class MultiVector & v) const = 0;
    virtual void AddTo (Complex s, class MultiVector & v) const = 0;

    virtual size_t Size() const = 0;
    virtual shared_ptr<BaseVector> CreateVector() const = 0;
    virtual void CalcComponent(size_t nr, BaseVector & bv) const = 0;
  };
  
  class MultiVector
  {
    shared_ptr<BaseVector> refvec;
    Array<shared_ptr<BaseVector>> vecs;
  public:
    MultiVector (shared_ptr<BaseVector> v, size_t cnt): refvec (v)
    {
      Extend (cnt);
    }
    MultiVector (size_t size, size_t cnt, bool is_complex)
    {
      refvec = CreateBaseVector(size, is_complex, 1);
      Extend (cnt);
    }
    MultiVector (const MultiVector & v) = default;
    MultiVector (MultiVector && v) = default;
    
    virtual ~MultiVector() { } 
    bool IsComplex() const { return refvec->IsComplex(); }

    size_t Size() const { return vecs.Size(); }
    shared_ptr<BaseVector> operator[] (size_t i) const { return vecs[i]; }
    shared_ptr<BaseVector> RefVec() const { return refvec; }

    MultiVector Range(IntRange r) const;
    
    void Extend (size_t nr = 1) {
      for ([[maybe_unused]] auto i : ngstd::Range(nr))
        vecs.Append (refvec->CreateVector());
    }
    void Append (shared_ptr<BaseVector> v)
    {
      vecs.Append (v->CreateVector());
      *vecs.Last() = *v;
    }

    void AppendOrthogonalize (shared_ptr<BaseVector> v, BaseMatrix * ip);

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

    MultiVector & operator= (const MultiVector & v2)
    { 
      for ([[maybe_unused]] auto i : ngstd::Range(vecs))
        *vecs[i] = *v2.vecs[i];
      return *this;
    }
    
    MultiVector & operator= (const MultiVectorExpr & expr)
    { 
      expr.AssignTo(1, *this);
      return *this;
    }
    
    MultiVector operator+= (const MultiVectorExpr & expr)
    { 
      expr.AddTo(1, *this);
      return *this;
    }
    
    MultiVector & operator-= (const MultiVectorExpr & expr)
    { 
      expr.AddTo(-1, *this);
      return *this;
    }

    void Orthogonalize (BaseMatrix * ipmat);
    
    virtual Matrix<> InnerProductD (const MultiVector & v2) const;
    virtual Matrix<Complex> InnerProductC (const MultiVector & v2, bool conjugate = false) const;    
    virtual Matrix<> InnerProductD (const MultiVectorExpr & v2) const;
    virtual Matrix<Complex> InnerProductC (const MultiVectorExpr & v2, bool conjugate = false) const;    
    virtual Vector<> InnerProductD (const BaseVector & v2) const;
    virtual Vector<Complex> InnerProductC (const BaseVector & v2, bool conjugate = false) const;    
  };
  

  template <class T>
  void MultAdd (const BaseMatrix & mat, T s, const MultiVector & x, MultiVector & y);
  
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

    size_t Size() const override { return vec->Size(); }
    shared_ptr<BaseVector> CreateVector() const override;
    void CalcComponent(size_t nr, BaseVector & bv) const override;
  };

  
  // y += sum a(i) * x[i]
  template <class T>
  void Axpy (const Vector<T> & a, const MultiVector  & x, BaseVector & y);


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
}
#endif
