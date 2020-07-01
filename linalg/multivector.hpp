#ifndef FILE_MULTIVECTOR
#define FILE_MULTIVECTOR

/*********************************************************************/
/* File:   multivector.hpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   June 2010                                                 */
/*********************************************************************/

namespace ngla {

  // TODO: find out how to spread templates over more files

  template <class T>
  class MultiVector;  

  template <class T>
  class MultiVectorExpr
  {
  public:
    virtual ~MultiVectorExpr() { ; }
    virtual void AssignTo (T s, class MultiVector<T> & v) const = 0;
    virtual void AddTo (T s, class MultiVector<T> & v) const = 0;
  };
  
  template <class T>
  class MultiVector
  {
    shared_ptr<BaseVector> refvec;
    Array<shared_ptr<BaseVector>> vecs;
    bool complex;
  public:
    MultiVector (shared_ptr<BaseVector> v, size_t cnt): refvec (v)
    {
      complex = v->IsComplex();
      Expand (cnt);
    }
    MultiVector (size_t size, size_t cnt, bool is_complex) : complex(is_complex)
    {
      refvec = CreateBaseVector(size, is_complex, 1);
      Expand (cnt);
    }
  
    bool IsComplex() const {return complex; }

    size_t Size() const { return vecs.Size(); }
    shared_ptr<BaseVector> operator[] (size_t i) const { return vecs[i]; }
    void Expand (size_t nr = 1) {
      for ([[maybe_unused]] auto i : Range(nr))
        vecs.Append (refvec->CreateVector());
    }
    void Append (shared_ptr<BaseVector> v)
    {
      vecs.Append (v->CreateVector());
      *vecs.Last() = *v;
    }

    void operator= (const T v)
    {
      Complex tmp;
      for (auto vec : vecs) {
        if (complex) *vec = Complex(v);
        else if (typeid(v).name()!= typeid(tmp).name()) *vec = v;
      }
        // TODO: what if attempt of ill usage?
    }
    void operator= (const MultiVectorExpr<T> & expr) { 
      expr.AssignTo(T(1), *this); 
    }
    void operator+= (const MultiVectorExpr<T> & expr){ 
      expr.AddTo(T(1), *this); 
      }
    void operator-= (const MultiVectorExpr<T> & expr){ 
      expr.AddTo(T(-1), *this); 
    }
  };
  
  template <class T>
  void MultAdd (const class BaseMatrix & mat, T s, const MultiVector<T> & x, MultiVector<T> & y)
  {
    for (auto i : Range(x.Size()))
      mat.MultAdd (s, *x[i], *y[i]);
  }
  
  template <class T>
  class MatMultiVecExpr : public MultiVectorExpr<T>
  {
    shared_ptr<BaseMatrix> mat;
    shared_ptr<MultiVector<T>> vec;
  public:
    MatMultiVecExpr (shared_ptr<BaseMatrix> amat, shared_ptr<MultiVector<T>> avec)
      : mat(amat), vec(avec) { ; }
    void AssignTo (T s, MultiVector<T> & res) const override
    {
      res = 0.0;
      MultAdd (*mat, s, *vec, res);
    }
    void AddTo (T s, MultiVector<T> & res) const override
    { MultAdd (*mat, s, *vec, res); }
  };

  
  // y += sum a(i) * x[i]
  template <class T>
  void Axpy (const Vector<T> & a, const MultiVector<T>  & x, BaseVector & y)
  {
    for (auto i : Range(a))
      y += a(i) * *x[i];
  }

  // no template because of overwritten functions of DynamicBaseExpression,
  // TODO: do that better
  class MultiVecAxpyExpr : public DynamicBaseExpression
  {
    Vector<double> a;
    Vector<Complex> ac;

    shared_ptr<MultiVector<double>> x;
    shared_ptr<MultiVector<Complex>> xc;


  public:
    MultiVecAxpyExpr (Vector<double> aa, shared_ptr<MultiVector<double>> ax) : a(aa), x(ax) { ; }
    MultiVecAxpyExpr (Vector<Complex> aa, shared_ptr<MultiVector<Complex>> ax) : ac(aa), xc(ax) { ; }


    void AssignTo (double s, BaseVector & v) const override
    {
      v = 0.0;
      AddTo (s, v);
    }
    void AddTo (double s, BaseVector & v) const override
    {
      Vector<double> sa = s*a;
      Axpy(sa, *x, v);
    }
    void AssignTo (Complex s, BaseVector & v) const override
    {
      v = 0.0;
      AddTo (s, v);
    }

    void AddTo (Complex s, BaseVector & v) const override
    {
      Vector<Complex> sa = s*ac;
      Axpy(sa, *xc , v);
      
    }
    // void AssignTo (Complex s, BaseVector & v2) const override { } // missing
    // void AddTo (Complex s, BaseVector & v2) const  override { } // missing
  };

  template <class T>
  Matrix<T> InnerProduct (const MultiVector<T> & x, const MultiVector<T> & y) {
    Matrix<T> res(x.Size(), y.Size());
    for (auto i : Range(x.Size()))
      for (auto j : Range(y.Size()))
        res(i,j) = S_InnerProduct<T>(*x[i], *y[j]);
    return res;
  }

}
#endif