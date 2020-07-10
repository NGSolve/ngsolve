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


  MultiVector MultiVector :: Range(IntRange r) const
  {
    MultiVector mv2(refvec, 0);
    for (auto i : r)
      mv2.vecs.Append (vecs[i]);
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
    // TODO:  exception
    for (auto i : Range(x.Size()))
      mat.MultAdd (s, *x[i], *y[i]);
  }

  template void MultAdd (const BaseMatrix & mat, double s, const MultiVector & x, MultiVector & y);
  template void MultAdd (const BaseMatrix & mat, Complex s, const MultiVector & x, MultiVector & y);


  template <class T>
  void T_Orthogonalize (MultiVector & mv, BaseMatrix * ipmat)
  {
    if (ipmat)
      {
        auto tmp = mv.RefVec()->CreateVector();
        for (int i = 0; i < mv.Size(); i++)
          {
            *tmp = *ipmat * *mv[i];
            T norm = sqrt(InnerProduct(*tmp, *mv[i]));
            *mv[i] *= 1.0 / norm;
            for (int j = i+1; j < mv.Size(); j++)
              *mv[j] -= InnerProduct<T>(*tmp, *mv[j], true)/norm * *mv[i];
          }
      }
    else
      {
        for (int i = 0; i < mv.Size(); i++)
          {
            *mv[i] *= 1.0 / (*mv[i]).L2Norm();
            for (int j = i+1; j < mv.Size(); j++)
              *mv[j] -= InnerProduct<T>(*mv[i], *mv[j], true) * *mv[i];
          }
      }
  }
  
  void MultiVector :: Orthogonalize (BaseMatrix * ipmat)
  {
    if (IsComplex())
      T_Orthogonalize<Complex> (*this, ipmat);
    else
      T_Orthogonalize<double> (*this, ipmat);
  }

  
  Matrix<> MultiVector ::
  InnerProductD (const MultiVector & y) const
  {
    Matrix<double> res(Size(), y.Size());
    for (int i = 0; i < Size(); i++)
      for (int j = 0; j < y.Size(); j++)
        res(i,j) = vecs[i]->InnerProductD(*y[j]);
    return res;
  }
  
  Matrix<Complex> MultiVector ::
  InnerProductC (const MultiVector & y, bool conjugate) const
  {
    Matrix<Complex> res(Size(), y.Size());
    for (int i = 0; i < Size(); i++)
      for (int j = 0; j < y.Size(); j++)
        res(i,j) = vecs[i]->InnerProductC(*y[j], conjugate);
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
