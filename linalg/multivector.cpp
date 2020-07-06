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


  template <class T>
  void Axpy (const Vector<T> & a, const MultiVector  & x, BaseVector & y)
  {
    cout << "type in axpy " << typeid(a(0)).name() << endl;
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

  template void MultAdd (const class BaseMatrix & mat, double s, const MultiVector & x, MultiVector & y);
  template void MultAdd (const class BaseMatrix & mat, Complex s, const MultiVector & x, MultiVector & y);

  template <class T>
  Matrix<T> InnerProduct (const MultiVector & x, const MultiVector & y) {
    Matrix<T> res(x.Size(), y.Size());
    for (auto i : Range(x.Size()))
      for (auto j : Range(y.Size()))
        res(i,j) = S_InnerProduct<T>(*x[i], *y[j]);
    return res;
  }

  template Matrix<Complex> InnerProduct<Complex> (const MultiVector & x, const MultiVector & y);
  template Matrix<double> InnerProduct<double> (const MultiVector & x, const MultiVector & y);


}
