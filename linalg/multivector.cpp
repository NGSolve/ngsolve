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


// template <class T>
// void Axpy (const Vector<T> & a, const MultiVector<T>  & x, BaseVector & y)
// {
//   for (auto i : Range(a))
//     y += a(i) * *x[i];
// }

// Matrix<Complex> InnerProduct (const MultiVector<Complex>  & x, const MultiVector<Complex>  & y)
// {
//   Matrix<T> res(x.Size(), y.Size());
//   for (auto i : Range(x.Size()))
//     for (auto j : Range(y.Size()))
//       res(i,j) = InnerProduct(*x[i], *y[j]);
//   return res;
// }

// template <class T>
// void MultAdd (const BaseMatrix & mat, T s, const MultiVector<T>  & x, MultiVector<T>  & y)
// {
//   for (auto i : Range(x.Size()))
//     mat.MultAdd (s, *x[i], *y[i]);
// }

}
