#include <la.hpp>
using namespace ngla;

namespace ngla {
  
  MultiVector :: MultiVector (shared_ptr<BaseVector> v, size_t cnt)
    : refvec (v)
  {
    Expand (cnt);
  }
  
  void MultiVector :: operator= (double v)
  {
    for (auto vec : vecs)
      *vec = v;
  }
  
void MultiVector :: Expand (size_t nr)
{
  for ([[maybe_unused]] auto i : Range(nr))
    vecs.Append (refvec->CreateVector());
}

  void MultiVector :: Append (shared_ptr<BaseVector> v)
  {
    vecs.Append (v->CreateVector());
    *vecs.Last() = *v;
  }


  
void Axpy (const Vector<> & a, const MultiVector & x, BaseVector & y)
{
  for (auto i : Range(a))
    y += a(i) * *x[i];
}

Matrix<> InnerProduct (const MultiVector & x, const MultiVector & y)
{
  Matrix<> res(x.Size(), y.Size());
  for (auto i : Range(x.Size()))
    for (auto j : Range(y.Size()))
      res(i,j) = InnerProduct(*x[i], *y[j]);
  return res;
}

void MultAdd (const BaseMatrix & mat, double s, const MultiVector & x, MultiVector & y)
{
  for (auto i : Range(x.Size()))
    mat.MultAdd (s, *x[i], *y[i]);
}

}
