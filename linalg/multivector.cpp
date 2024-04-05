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








  MultiVector & MultiVector::operator= (const MultiVector & v2)
  {
    if (Size() != v2.Size())
      throw Exception("MultiVector assignment sizes mismatch, my size = "
                      + ToString(Size()) + " other size = " + ToString(v2.Size()));
    
    for (auto i : ngstd::Range(vecs))
      *vecs[i] = *v2.vecs[i];
    return *this;
  }
  
  MultiVector & MultiVector::operator= (const MultiVectorExpr & expr)
  {
    if (Size() != expr.Size())
      throw Exception("MultiVector assignment sizes mismatch, my size = "
                      + ToString(Size()) + " other size = " + ToString(expr.Size()));
    
    Vector<double> ones(Size());
    ones = 1;
    expr.AssignTo(ones, *this);
    return *this;
  }
  
  MultiVector MultiVector::operator+= (const MultiVectorExpr & expr)
  {
    if (Size() != expr.Size())
      throw Exception("MultiVector assignment-add sizes mismatch, my size = "
                      + ToString(Size()) + " other size = " + ToString(expr.Size()));
    
    Vector<double> ones(Size());
    ones = 1;      
    expr.AddTo(ones, *this);
    return *this;
  }
  
  MultiVector & MultiVector::operator-= (const MultiVectorExpr & expr)
  {
    if (Size() != expr.Size())
      throw Exception("MultiVector assignment-sub sizes mismatch, my size = "
                      + ToString(Size()) + " other size = " + ToString(expr.Size()));

    Vector<double> mones(Size());
    mones = -1;
    expr.AddTo(mones, *this);
    return *this;
  }
  




  
  unique_ptr<MultiVector> MultiVector :: Range(IntRange r) const
  {
    auto mv2 = make_unique<MultiVector>(refvec, 0);
    for (auto i : r)
      mv2->vecs.Append (vecs[i]);
    return mv2;
  }

  unique_ptr<MultiVector> MultiVector :: VectorRange(IntRange r) const
  {
    auto mv2 = make_unique<MultiVector>(refvec->Range(r), 0);
    for (auto v : vecs)
      mv2->vecs.Append (v->Range(r));
    return mv2;
  }

  
  unique_ptr<MultiVector> MultiVector :: SubSet(const Array<int> & indices) const
  {
    auto mv2 = make_unique<MultiVector>(refvec, 0);
    for (auto i : indices)
      mv2->vecs.Append (vecs[i]);
    return mv2;
  }

  

  void MultiVector :: AssignTo (FlatVector<double> s, class MultiVector & v) const
  {
    for (int i = 0; i < s.Size(); i++)
      *v[i] = s[i] * *vecs[i];
  }
  
  void MultiVector :: AddTo (FlatVector<double> s, class MultiVector & v) const 
  {
    for (int i = 0; i < s.Size(); i++)
      *v[i] += s[i] * *vecs[i];
  }
  
  void MultiVector ::AssignTo (FlatVector<Complex> s, class MultiVector & v) const
  {
    for (int i = 0; i < s.Size(); i++)
      *v[i] = s[i] * *vecs[i];
  }
        
  void MultiVector ::AddTo (FlatVector<Complex> s, class MultiVector & v) const
  {
    for (int i = 0; i < s.Size(); i++)
      *v[i] += s[i] * *vecs[i];
  }

  shared_ptr<BaseVector> MultiVector ::CreateVector() const
  {
    return refvec->CreateVector();
  }
    
  void MultiVector ::CalcComponent(size_t nr, BaseVector & bv) const
  {
    bv = *vecs[nr];
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
  void MultAdd (const class BaseMatrix & mat, FlatVector<T> s, const MultiVector & x, MultiVector & y)
  {
    if constexpr (std::is_same<T,double>::value)
      {
        mat.MultAdd (s, x, y);
      }
    else
      for (auto i : Range(x.Size()))
        mat.MultAdd (s(i), *x[i], *y[i]);
  }
  
  template void MultAdd (const BaseMatrix & mat, FlatVector<double> s, const MultiVector & x, MultiVector & y);
  template void MultAdd (const BaseMatrix & mat, FlatVector<Complex> s, const MultiVector & x, MultiVector & y);


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
  Matrix<T> MultiVector::T_Orthogonalize (BaseMatrix * ipmat)
  {
    static Timer t("MultiVector::Orthogonalize");
    RegionTimer reg(t);
    auto & mv = *this;
    Matrix<T> Rfactor(mv.Size());
    Rfactor = T(0.0);
    if (ipmat)
      {
        auto tmp = mv.RefVec()->CreateVector();
        for (int i = 0; i < mv.Size(); i++)
          {
            *tmp = *ipmat * *mv[i];
            double norm = sqrt(fabs(InnerProduct<T>(*tmp, *mv[i], true)));
            Rfactor(i,i) = norm;
            *mv[i] *= 1.0 / norm;
            /*
            for (int j = i+1; j < mv.Size(); j++)
              {
                T rij = InnerProduct<T>(*mv[j], *tmp, true) / norm;
                Rfactor(i,j) = rij;
                *mv[j] -= rij * *mv[i];
              }
            */
            IntRange rest(i+1, mv.Size());
            auto mvrest = mv.Range(rest);
            Rfactor.Row(i).Range(rest) = (1.0/norm) * mvrest->T_InnerProduct<T> (*tmp, true);
            for (auto j : rest)
              *mv[j] -= Rfactor(i,j) * *mv[i];              
          }
      }
    else
      {
        if (mv.Size() == 1)
          {
            double norm = (*mv[0]).L2Norm();
            Rfactor(0,0) = norm;
            *mv[0] *= 1.0 / norm;
          }
        else
          {
            int half = mv.Size()/2;
            auto r1 = IntRange(0,half);
            auto r2 = IntRange(half, mv.Size());
            auto mv1 = mv.Range(r1);
            auto mv2 = mv.Range(r2);
            Rfactor.Rows(r1).Cols(r1) = mv1->T_Orthogonalize<T>(ipmat);
            Matrix<T> ip = Conj(InnerProduct<T> (*mv1, *mv2, true));
            Rfactor.Rows(r2).Cols(r1) = ip;
            ip *= -1;
            mv2->Add (*mv1, ip);
            Rfactor.Rows(r2).Cols(r2) = mv2->T_Orthogonalize<T>(ipmat);
          }
      }
    return Rfactor;
  }
  
  void MultiVector :: Orthogonalize (BaseMatrix * ipmat)
  {
    if (IsComplex())
      this->T_Orthogonalize<Complex> (ipmat);
    else
      this->T_Orthogonalize<double> (ipmat);
  }


  void MultiVector :: AppendOrthogonalize (shared_ptr<BaseVector> v, BaseMatrix * ipmat,
                                           bool parallel, int iterations)
  {
    if (!IsComplex())
      this->T_AppendOrthogonalize<double> (v, ipmat, parallel, iterations);
    else
      this->T_AppendOrthogonalize<Complex> (v, ipmat, parallel, iterations);
  }
  
  template <typename T>
  Vector<T> MultiVector :: T_AppendOrthogonalize (shared_ptr<BaseVector> v, BaseMatrix * ipmat,
                                           bool parallel, int iterations)
  {
    size_t osize = this->Size();
    Vector<T> R(osize+1);
    R = 0;
    auto oldR = R.Range(0, osize);
    auto hv = v->CreateVector();
    hv = *v;

    // if (!IsComplex())
    // {
        if (!ipmat)
          {
            for (int j = 0; j < iterations; j++)
              {
                if (parallel)
                  {
                    // Vector<T> prod = -this->InnerProduct<T> (hv);
                    Vector<T> prod = -Conj(this->T_InnerProduct<T> (hv, true));
                    Axpy (prod, *this, hv);
                    oldR -= prod;
                  }
                else
                  {
                    for (int i = 0; i < osize; i++)
                      {
                        T ip = -Conj(InnerProduct<T> (*(*this)[i], hv, true));
                        hv += ip * *(*this)[i];
                        R(i) -= ip;
                      }
                  }
              }
            
            double norm = hv.L2Norm();
            norm = sqrt(fabs(InnerProduct<T>(hv, hv, true)));
            R(osize) = norm;
            hv /= norm;
          }
        else
          {
            auto hv2 = v->CreateVector();              
            for (int j = 0; j < iterations; j++)          
              {
                if (parallel)
                  {
                    hv2 = *ipmat * hv;
                    Vector<T> prod = -this->T_InnerProduct<T> (hv2);
                    Axpy (prod, *this, hv);
                    oldR -= prod;                    
                  }
                else
                  {
                    for (int i = 0; i < this->Size(); i++)
                      {
                        hv2 = *ipmat * hv;
                        T ip = -::InnerProduct<T> (*(*this)[i], hv2);
                        hv += ip * *(*this)[i];
                        R(i) -= ip;                        
                      }
                  }
              }
            hv2 = *ipmat * hv;
            // double norm = sqrt(InnerProduct(hv, hv2));
            double norm = sqrt(fabs(InnerProduct<T>(hv, hv2)));
            R(osize) = norm;            
            hv /= norm;
          }
        // }
    vecs.Append (hv);
    return R;
  }
  
  void MultiVector :: SetScalar (double s)
  {
    for (auto & vec : vecs) 
      *vec = s;
  }
  
  void MultiVector :: SetScalar (Complex s)
  {
    for (auto & vec : vecs) 
      *vec = s;
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

    Matrix<double> res(y.Size(), Size());
    for (int i = 0; i < Size(); i++)
      for (int j = 0; j < y.Size(); j++)
        res(j,i) = vecs[i]->InnerProductD(*y[j]);
    return res;
  }
  
  Matrix<Complex> MultiVector ::
  InnerProductC (const MultiVector & y, bool conjugate) const
  {
    static Timer t("MultiVector::InnerProductC");
    RegionTimer reg(t);

    Matrix<Complex> res(y.Size(), Size());
    for (int i = 0; i < Size(); i++)
      for (int j = 0; j < y.Size(); j++)
        res(j,i) = vecs[i]->InnerProductC(*y[j], conjugate);
    return res;
  }


  Matrix<> MultiVector ::
  InnerProductD (const MultiVectorExpr & y) const
  {
    static Timer t("InnerProductD (MulitVector,MultiVectorExpr)");
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
    Vector<double> ones(y.Size());
    ones = 1;
    y.AssignTo (ones, *mv);
    return InnerProductD (*mv);
  }
  
  Matrix<Complex> MultiVector ::
  InnerProductC (const MultiVectorExpr & y, bool conjugate) const
  {
    static Timer t("MultiVector::InnerProductC");
    RegionTimer reg(t);

    Matrix<Complex> res(y.Size(), Size());
    shared_ptr<BaseVector> hy = y.CreateVector();    
    for (int j = 0; j < y.Size(); j++)
      {
        y.CalcComponent(j, *hy);
        res.Row(j) = InnerProductC(*hy, conjugate);
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






  shared_ptr<MultiVectorExpr> operator+ (shared_ptr<MultiVectorExpr> e1,
                                                shared_ptr<MultiVectorExpr> e2)
  {
    if (e1->Size() != e2->Size())
      throw Exception("MultiVector+ sizes don't fit: " + ToString(e1->Size()) + " != " + ToString(e2->Size()));
    return make_shared<SumMultiVectorExpr> (e1,e2);
  }
                                         
  shared_ptr<MultiVectorExpr> operator- (shared_ptr<MultiVectorExpr> e1)
  {
    Vector<> mones(e1->Size()); mones = -1;
    return make_shared<ScaledMultiVectorExpr<double>> (e1, mones);
  }

  shared_ptr<MultiVectorExpr> operator- (shared_ptr<MultiVectorExpr> e1,
						shared_ptr<MultiVectorExpr> e2)
  {
    return e1 + (-e2);
  }
                                         
  shared_ptr<MultiVectorExpr> operator* (double s, shared_ptr<MultiVectorExpr> e1)
  {
    Vector<double> scale(e1->Size()); scale = s;
    return make_shared<ScaledMultiVectorExpr<double>> (e1, scale);
  }
  
  shared_ptr<MultiVectorExpr> operator* (Complex s, shared_ptr<MultiVectorExpr> e1)
  {
    Vector<Complex> scale(e1->Size()); scale = s;
    return make_shared<ScaledMultiVectorExpr<Complex>> (e1, scale);
  }  


}
