#include <la.hpp>
using namespace ngla;



namespace ngla {
  
  // template <class T>
  // MultiVector<T> :: MultiVector(shared_ptr<BaseVector> v, size_t cnt)
  //   : refvec (v), complex(false)
  // {
  //   Expand (cnt);
  // }

  // template <class T>
  // MultiVector<T> :: MultiVector (size_t size, size_t cnt, bool is_complex)
  //   : complex(is_complex)
  // {
  //   refvec = CreateBaseVector(size, is_complex, 1);
  //   Expand (cnt);
  // }
  
  // template <class T>
  // void MultiVector :: operator= (T v)
  // {
  //   Complex tmp;
  //   for (auto vec : vecs)
  //     if (complex) *vec = Complex(v);
  //     else if (typeid(v).name()!= typeid(tmp).name()) *vec = v;
  //     // TODO: what if attempt of ill usage? 
  // }

  // TODO: cast double inputs to complex for complex vector
  // template <class T>
  // void MultiVector<T> :: operator= (T v)
  // {
    
  //   for (auto vec : vecs)
  //     *vec = v;
  // }

  // void MultiVector :: operator= (Complex v)
  // {
  //   for (auto vec : vecs)
  //     *vec = v;
  // }
  
// void MultiVector :: Expand (size_t nr)
// {
//   for ([[maybe_unused]] auto i : Range(nr))
//     vecs.Append (refvec->CreateVector());
// }

  // void MultiVector :: Append (shared_ptr<BaseVector> v)
  // {
  //   vecs.Append (v->CreateVector());
  //   *vecs.Last() = *v;
  // }


  unique_ptr<MultiVector> MultiVector :: Range(IntRange r) const
  {
    auto mv2 = make_unique<MultiVector>(refvec, 0);
    for (auto i : r)
      mv2->vecs.Append (vecs[i]);
    return mv2;
  }

  

  template <class T>
  void Axpy (const Vector<T> & a, const MultiVector  & x, BaseVector & y)
  {
    for (auto i : Range(a))
      y += a(i) * *x[i];
  }
  
  template void Axpy (const Vector<double> & a, const MultiVector  & x, BaseVector & y);
  template void Axpy (const Vector<Complex> & a, const MultiVector  & x, BaseVector & y);

  template <class T>
  void MultAdd (const class BaseMatrix & mat, T s, const MultiVector & x, MultiVector & y)
  {
    if constexpr (std::is_same<T,double>::value)
      {
        Vector<double> vs(x.Size());
        vs = s;
        mat.MultAdd (vs, x, y);
      }
    else
      for (auto i : Range(x.Size()))
        mat.MultAdd (s, *x[i], *y[i]);
  }
  
  template void MultAdd (const BaseMatrix & mat, double s, const MultiVector & x, MultiVector & y);
  template void MultAdd (const BaseMatrix & mat, Complex s, const MultiVector & x, MultiVector & y);


  shared_ptr<BaseVector> MatMultiVecExpr :: CreateVector() const 
  { return mat->CreateColVector(); }
    
  void MatMultiVecExpr :: CalcComponent(size_t nr, BaseVector & bv) const 
  {
    bv = *mat * *(*vec)[nr];
  }
  


  template <typename T = double>
  inline Matrix<T> InnerProduct (const MultiVector & v1, const MultiVector & v2, bool conjugate = false)
  {
    if constexpr (is_same<T,double>::value)
                   return v1.InnerProductD(v2);
    else
      return v1.InnerProductC(v2, conjugate);
  }


  
  template <class T>
  void T_Orthogonalize (MultiVector & mv, BaseMatrix * ipmat)
  {
    if (ipmat)
      {
        auto tmp = mv.RefVec()->CreateVector();
        for (int i = 0; i < mv.Size(); i++)
          {
            *tmp = *ipmat * *mv[i];
            double norm = sqrt(fabs(InnerProduct<T>(*tmp, *mv[i], true)));
            *mv[i] *= 1.0 / norm;
            for (int j = i+1; j < mv.Size(); j++)
              *mv[j] -= InnerProduct<T>(*tmp, *mv[j], true)/norm * *mv[i];
          }
      }
    else
      {
        if (mv.Size() == 1)
          *mv[0] *= 1.0 / (*mv[0]).L2Norm();
        else
          {
            int half = mv.Size()/2;
            auto mv1 = mv.Range(IntRange(0, half));
            auto mv2 = mv.Range(IntRange(half, mv.Size()));
            mv1->Orthogonalize(ipmat);
            Matrix<T> ip = InnerProduct<T> (*mv1, *mv2, true);
            ip *= -1;
            mv2->Add (*mv1, ip);
            mv2->Orthogonalize(ipmat);
          }
        /*
        for (int i = 0; i < mv.Size(); i++)
          {
            *mv[i] *= 1.0 / (*mv[i]).L2Norm();
            for (int j = i+1; j < mv.Size(); j++)
              *mv[j] -= InnerProduct<T>(*mv[i], *mv[j], true) * *mv[i];
          }
        */
      }
  }
  
  void MultiVector :: Orthogonalize (BaseMatrix * ipmat)
  {
    static Timer t("MultiVector::Orthogonalize");
    RegionTimer reg(t);
    if (IsComplex())
      T_Orthogonalize<Complex> (*this, ipmat);
    else
      T_Orthogonalize<double> (*this, ipmat);
  }



  void MultiVector :: AppendOrthogonalize (shared_ptr<BaseVector> v, BaseMatrix * ipmat)
  {
    vecs.Append (v->CreateVector());
    *vecs.Last() = *v;

    if (IsComplex())    
      T_Orthogonalize<Complex> (*this, ipmat);
    else
      T_Orthogonalize<double> (*this, ipmat);
  }


  void MultiVector :: AddTo (FlatVector<double> vec, BaseVector & v2)
  {
    for (int i = 0; i < vec.Size(); i++)
      v2 += vec(i) * *vecs[i];
  }
  
  void MultiVector :: AddTo (FlatVector<Complex> vec, BaseVector & v2)
  {
    for (int i = 0; i < vec.Size(); i++)
      v2 += vec(i) * *vecs[i];
  }

  // me[i] += v2[j] mat(j,i)
  void MultiVector :: Add (const MultiVector & v2, FlatMatrix<double> mat)
  {
    for (int i = 0; i < mat.Width(); i++)
      for (int j = 0; j < mat.Height(); j++)
        *vecs[i] += mat(j,i) * *v2.vecs[j];
  }
    
  void MultiVector :: Add (const MultiVector & v2, FlatMatrix<Complex> mat)
  {
    for (int i = 0; i < mat.Width(); i++)
      for (int j = 0; j < mat.Height(); j++)
        *vecs[i] += mat(j,i) * *v2.vecs[j];
  }


  
  
  Matrix<> MultiVector ::
  InnerProductD (const MultiVector & y) const
  {
    static Timer t("MultiVector::InnerProductD");
    RegionTimer reg(t);

    Matrix<double> res(Size(), y.Size());
    for (int i = 0; i < Size(); i++)
      for (int j = 0; j < y.Size(); j++)
        res(i,j) = vecs[i]->InnerProductD(*y[j]);
    return res;
  }
  
  Matrix<Complex> MultiVector ::
  InnerProductC (const MultiVector & y, bool conjugate) const
  {
    static Timer t("MultiVector::InnerProductC");
    RegionTimer reg(t);

    Matrix<Complex> res(Size(), y.Size());
    for (int i = 0; i < Size(); i++)
      for (int j = 0; j < y.Size(); j++)
        res(i,j) = vecs[i]->InnerProductC(*y[j], conjugate);
    return res;
  }


  Matrix<> MultiVector ::
  InnerProductD (const MultiVectorExpr & y) const
  {
    static Timer t("MultiVector::InnerProductD");
    RegionTimer reg(t);

    /*
    Matrix<double> res(Size(), y.Size());
    shared_ptr<BaseVector> hy = y.CreateVector();
    for (int j = 0; j < y.Size(); j++)
      {
        y.CalcComponent(j, *hy);
        res.Col(j) = InnerProductD(*hy);
      }
    return res;
    */
    auto mv = y.CreateVector()->CreateMultiVector(y.Size());
    y.AssignTo (1, *mv);
    return InnerProductD (*mv);
  }
  
  Matrix<Complex> MultiVector ::
  InnerProductC (const MultiVectorExpr & y, bool conjugate) const
  {
    static Timer t("MultiVector::InnerProductC");
    RegionTimer reg(t);

    Matrix<Complex> res(Size(), y.Size());
    shared_ptr<BaseVector> hy = y.CreateVector();    
    for (int j = 0; j < y.Size(); j++)
      {
        y.CalcComponent(j, *hy);
        res.Col(j) = InnerProductC(*hy, conjugate);
      }
    return res;
  }


  



  Vector<> MultiVector ::
  InnerProductD (const BaseVector & y) const
  {
    Vector<double> res(Size());
    for (int i = 0; i < Size(); i++)
      res(i) = vecs[i]->InnerProductD(y);
    return res;
  }
  
  Vector<Complex> MultiVector ::
  InnerProductC (const BaseVector & y, bool conjugate) const
  {
    Vector<Complex> res(Size());
    for (int i = 0; i < Size(); i++)
      res(i) = vecs[i]->InnerProductC(y, conjugate);
    return res;
  }



}
