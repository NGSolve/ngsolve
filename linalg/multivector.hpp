#ifndef FILE_MULTIVECTOR
#define FILE_MULTIVECTOR

/*********************************************************************/
/* File:   multivector.hpp                                           */
/* Author: Joachim Schoeberl, Amanda Schoefl                         */
/* Date:   June 2020                                                 */
/*********************************************************************/

#include "basevector.hpp"

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
    virtual void AssignTo (FlatVector<double> s, class MultiVector & v) const = 0;
    virtual void AddTo (FlatVector<double> s, class MultiVector & v) const = 0;
    virtual void AssignTo (FlatVector<Complex> s, class MultiVector & v) const = 0;
    virtual void AddTo (FlatVector<Complex> s, class MultiVector & v) const = 0;

    virtual size_t Size() const = 0;
    virtual shared_ptr<BaseVector> CreateVector() const = 0;
    virtual void CalcComponent(size_t nr, BaseVector & bv) const = 0;
  };

  
  class MultiVector : public MultiVectorExpr
  {
  protected:
    shared_ptr<BaseVector> refvec;
    Array<shared_ptr<BaseVector>> vecs;
  public:
    MultiVector (const Array<shared_ptr<BaseVector>> & avecs)
      : refvec(avecs.Size()?avecs[0]:nullptr), vecs(avecs) { } 
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

    size_t Size() const override { return vecs.Size(); }
    shared_ptr<BaseVector> operator[] (size_t i) const { return vecs[i]; }
    shared_ptr<BaseVector> RefVec() const { return refvec; }

    virtual unique_ptr<MultiVector> Range(IntRange r) const;
    virtual unique_ptr<MultiVector> VectorRange(IntRange r) const; // range of each component vector
    virtual unique_ptr<MultiVector> SubSet(const Array<int> & indices) const;
    
    void Extend (size_t nr = 1) {
      for ([[maybe_unused]] auto i : ngstd::Range(nr))
        vecs.Append (refvec->CreateVector());
    }
    void Append (shared_ptr<BaseVector> v)
    {
      vecs.Append (v->CreateVector());
      *vecs.Last() = *v;
    }

    void Replace (int i, shared_ptr<BaseVector> v)
    {
      vecs[i] = v;
    }

    void AppendOrthogonalize (shared_ptr<BaseVector> v, BaseMatrix * ip, bool parallel, int iterations);
    template <typename T>
    Vector<T> T_AppendOrthogonalize (shared_ptr<BaseVector> v, BaseMatrix * ip, bool parallel, int iterations);

    MultiVector & operator= (double val)
    {
      SetScalar (val);
      return *this;
    }

    MultiVector & operator= (Complex val)
    {
      SetScalar (val);
      return *this;      
    }

    MultiVector & operator= (const MultiVector & v2);
    MultiVector & operator= (const MultiVectorExpr & expr);
    MultiVector operator+= (const MultiVectorExpr & expr);
    MultiVector & operator-= (const MultiVectorExpr & expr);

    void Orthogonalize (BaseMatrix * ipmat);
    template <class T>
    Matrix<T> T_Orthogonalize (BaseMatrix * ipmat);
    
    virtual void SetScalar (double s);
    virtual void SetScalar (Complex s);
    
    // v2 += vec(i) * v[i]
    virtual void AddTo (FlatVector<double> vec, BaseVector & v2);
    virtual void AddTo (FlatVector<Complex> vec, BaseVector & v2);
    // me[i] += v2[j] mat(j,i)
    virtual void Add (const MultiVector & v2, FlatMatrix<double> mat);
    virtual void Add (const MultiVector & v2, FlatMatrix<Complex> mat);

    template <typename T>
    auto T_InnerProduct (const MultiVector & v2, bool conjugate = false)
    {
      if constexpr (is_same<T, double>())
                     return InnerProductD (v2);
      else
        return InnerProductC (v2, conjugate);
    }
    
    virtual Matrix<> InnerProductD (const MultiVector & v2) const;
    virtual Matrix<Complex> InnerProductC (const MultiVector & v2, bool conjugate = false) const;    
    virtual Matrix<> InnerProductD (const MultiVectorExpr & v2) const;
    virtual Matrix<Complex> InnerProductC (const MultiVectorExpr & v2, bool conjugate = false) const;


    template <typename T>
    auto T_InnerProduct (const BaseVector & v2, bool conjugate = false)
    {
      if constexpr (is_same<T, double>())
                     return InnerProductD (v2);
        else
          return InnerProductC (v2, conjugate);
    }
    
    virtual Vector<> InnerProductD (const BaseVector & v2) const;
    virtual Vector<Complex> InnerProductC (const BaseVector & v2, bool conjugate = false) const;



    virtual void AssignTo (FlatVector<double> s, class MultiVector & v) const override;
    virtual void AddTo (FlatVector<double> s, class MultiVector & v) const override;
    virtual void AssignTo (FlatVector<Complex> s, class MultiVector & v) const override;
    virtual void AddTo (FlatVector<Complex> s, class MultiVector & v) const override;

    virtual shared_ptr<BaseVector> CreateVector() const override;
    virtual void CalcComponent(size_t nr, BaseVector & bv) const override;
  };



  
  class BaseVectorPtrMV : public MultiVector
  {
  public:
    using MultiVector::MultiVector;

    unique_ptr<MultiVector> Range(IntRange r) const override;
    void SetScalar (double s) override;
    void Add (const MultiVector & v2, FlatMatrix<double> mat) override;
    void Add (const MultiVector & v2, FlatMatrix<Complex> mat) override;
    Vector<> InnerProductD (const BaseVector & v2) const override;
    Matrix<> InnerProductD (const MultiVector & v2) const override;
    Matrix<Complex> InnerProductC (const MultiVector & v2, bool conjugate) const override;    
  };
  
    

  

  template <class T>
  void MultAdd (const BaseMatrix & mat, FlatVector<T> s, const MultiVector & x, MultiVector & y);
  
  class MatMultiVecExpr : public MultiVectorExpr
  {
    shared_ptr<BaseMatrix> mat;
    shared_ptr<MultiVector> vec;
  public:
    MatMultiVecExpr (shared_ptr<BaseMatrix> amat, shared_ptr<MultiVector> avec)
      : mat(amat), vec(avec) { ; }

    void AssignTo (FlatVector<double> s, MultiVector & res) const override
    {
      res = 0.0;
      MultAdd (*mat, s, *vec, res);
    }

    void AddTo (FlatVector<double> s, MultiVector & res) const override
    { MultAdd (*mat, s, *vec, res); }

    void AssignTo (FlatVector<Complex> s, MultiVector & res) const override
    {
      res = 0.0;
      MultAdd (*mat, s, *vec, res);
    }

    void AddTo (FlatVector<Complex> s, MultiVector & res) const override
    { MultAdd (*mat, s, *vec, res); }

    size_t Size() const override { return vec->Size(); }
    shared_ptr<BaseVector> CreateVector() const override;
    void CalcComponent(size_t nr, BaseVector & bv) const override;
  };





  template <typename T>
  class MultiVecMatrixExpr : public MultiVectorExpr
  {
    Matrix<T> mat;
    shared_ptr<MultiVector> vec;
  public:
    MultiVecMatrixExpr (Matrix<T> amat, shared_ptr<MultiVector> avec)
      : mat(amat), vec(avec)
    {
      if (vec->Size() != mat.Height())
        throw Exception("Multivector * Matrix don't fit: mv.size = "+ToString(vec->Size())+
                        ", matrix.height = " + ToString(mat.Height()));
      
    }

    void AssignTo (FlatVector<double> s, MultiVector & res) const override
    {
      res = 0.0;
      AddTo (s, res);
    }

    void AddTo (FlatVector<double> s, MultiVector & res) const override
    {
      Matrix<T> hmat = mat;
      for (auto i : Range(mat.Width()))
        hmat.Col(i) *= s(i);
      res.Add (*vec, hmat);
    }

    void AssignTo (FlatVector<Complex> s, MultiVector & res) const override
    {
      res = 0.0;
      AddTo (s, res);      
    }

    void AddTo (FlatVector<Complex> s, MultiVector & res) const override
    {
      Matrix<Complex> hmat = mat;
      for (auto i : Range(mat.Width()))
        hmat.Col(i) *= s(i);
      res.Add (*vec, hmat);
    }

    size_t Size() const override { return mat.Width(); } // vec->Size(); }
    shared_ptr<BaseVector> CreateVector() const override
    {
      return vec->RefVec()->CreateVector();
    }
      
    void CalcComponent(size_t nr, BaseVector & bv) const override
    {
      bv = 0;
      vec->AddTo (Vector<T>(mat.Col(nr)), bv);
    }
  };





  class SumMultiVectorExpr : public MultiVectorExpr
  {
    shared_ptr<MultiVectorExpr> e1;
    shared_ptr<MultiVectorExpr> e2;
  public:
    SumMultiVectorExpr (shared_ptr<MultiVectorExpr> ae1,
                        shared_ptr<MultiVectorExpr> ae2)
      : e1(ae1), e2(ae2) { } 

    void AssignTo (FlatVector<double> s, MultiVector & res) const override
    {
      e1->AssignTo (s, res);
      e2->AddTo (s, res);
    }

    void AddTo (FlatVector<double> s, MultiVector & res) const override
    {
      e1->AddTo (s, res);
      e2->AddTo (s, res);
    }

    void AssignTo (FlatVector<Complex> s, MultiVector & res) const override
    {
      e1->AssignTo (s, res);
      e2->AddTo (s, res);
    }

    void AddTo (FlatVector<Complex> s, MultiVector & res) const override
    {
      e1->AddTo (s, res);
      e2->AddTo (s, res);
    }

    size_t Size() const override { return e1->Size(); }
    shared_ptr<BaseVector> CreateVector() const override
    {
      return e1->CreateVector();
    }
      
    void CalcComponent(size_t nr, BaseVector & bv) const override
    {
      auto tmp = bv.CreateVector();
      e1->CalcComponent(nr, *tmp);
      e2->CalcComponent(nr, bv);
      bv += *tmp;
    }
  };




  template <typename T>
  class ScaledMultiVectorExpr : public MultiVectorExpr
  {
    shared_ptr<MultiVectorExpr> mv;
    Vector<T> scale;
  public:
    ScaledMultiVectorExpr (shared_ptr<MultiVectorExpr> amv,
                           Vector<T> ascale)
      : mv(amv), scale(ascale) { } 

    void AssignTo (FlatVector<double> s, MultiVector & res) const override
    {
      Vector<T> hv = pw_mult(scale, s);
      mv->AssignTo (hv, res);
    }

    void AddTo (FlatVector<double> s, MultiVector & res) const override
    {
      Vector<T> hv = pw_mult(scale, s);
      mv->AddTo (hv, res);
    }

    void AssignTo (FlatVector<Complex> s, MultiVector & res) const override
    {
      Vector<Complex> hv = pw_mult(scale, s);
      mv->AssignTo (hv, res);
    }

    void AddTo (FlatVector<Complex> s, MultiVector & res) const override
    {
      Vector<Complex> hv = pw_mult(scale, s);
      mv->AddTo (hv, res);
    }

    size_t Size() const override { return mv->Size(); }
    shared_ptr<BaseVector> CreateVector() const override
    {
      return mv->CreateVector();
    }
      
    void CalcComponent(size_t nr, BaseVector & bv) const override
    {
      mv->CalcComponent(nr, bv);
      bv *= scale(nr);
    }
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

    AutoVector CreateVector() const override
    { return x->RefVec()->CreateVector(); }    
    
    void AssignTo (double s, BaseVector & v) const override
    {
      v = 0.0;
      AddTo (s, v);
    }
    
    void AddTo (double s, BaseVector & v) const override
    {
      Vector<TSCAL> sa = s*a;
      // Axpy (sa, *x, v);
      x->AddTo(sa, v);
    }

    void AssignTo (Complex s, BaseVector & v) const override
    {
      v = 0.0;
      AddTo (s, v);
    }

    void AddTo (Complex s, BaseVector & v) const override
    {
      Vector<Complex> sa = s*a;
      // Axpy(sa, *x , v);
      x->AddTo(sa, v);      
    }
  };




  shared_ptr<MultiVectorExpr> operator+ (shared_ptr<MultiVectorExpr> e1,
                                         shared_ptr<MultiVectorExpr> e2);
  shared_ptr<MultiVectorExpr> operator- (shared_ptr<MultiVectorExpr> e1);

  shared_ptr<MultiVectorExpr> operator- (shared_ptr<MultiVectorExpr> e1,
                                         shared_ptr<MultiVectorExpr> e2);
                                         
  shared_ptr<MultiVectorExpr> operator* (double s, shared_ptr<MultiVectorExpr> e1);
  
  shared_ptr<MultiVectorExpr> operator* (Complex s, shared_ptr<MultiVectorExpr> e1);
  
}
#endif
