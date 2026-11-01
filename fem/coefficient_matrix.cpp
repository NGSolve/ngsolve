
#include <core/register_archive.hpp>
#include <ngstd.hpp>
#include <bla.hpp>

#include <../basiclinalg/complex_wrapper.hpp>
#include <coefficient.hpp>
#include <coefficient_impl.hpp>
#include "symbolicintegrator.hpp"
#include <../ngstd/evalfunc.hpp>
#include <algorithm>
#ifdef NGS_PYTHON
#include <core/python_ngcore.hpp> // for shallow archive
#endif // NGS_PYTHON


#include "coefficient_matrix.hpp"

namespace ngfem
{
  MultMatMatCoefficientFunction :: ~MultMatMatCoefficientFunction() { }
  MultDiagMatCoefficientFunction :: ~MultDiagMatCoefficientFunction() { }
  TransposeCoefficientFunction :: ~TransposeCoefficientFunction() { }
  TraceCoefficientFunction :: ~TraceCoefficientFunction() { }    
  IdentityCoefficientFunction :: ~IdentityCoefficientFunction () {  }

  static int TailSize (FlatArray<int> dims, int from)
  {
    int s = 1;
    for (int i = from; i < dims.Size(); i++) s *= dims[i];
    return s;
  }

  static void AppendTail (Array<int> & dims, FlatArray<int> d, int from)
  {
    for (int i = from; i < d.Size(); i++) dims.Append (d[i]);
  }

  // r[h,b,a,t] = c[h,a,b,t], index relabeling only
  class SwapAxesCoefficientFunction : public T_CoefficientFunction<SwapAxesCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> c1;
    int head, n1, n2, tail;
    using BASE = T_CoefficientFunction<SwapAxesCoefficientFunction>;

    int Src (int i) const
    {
      int t = i % tail; i /= tail;
      int a = i % n1; i /= n1;
      int b = i % n2; int h = i / n2;
      return ((h*n1+a)*n2+b)*tail+t;
    }
  public:
    SwapAxesCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                 int ahead, int an1, int an2, int atail)
      : BASE(1, ac1->IsComplex()), c1(ac1), head(ahead), n1(an1), n2(an2), tail(atail)
    {
      if (head*n1*n2*tail != c1->Dimension())
        throw Exception("SwapAxes: dimensions don't fit");
      SetDimensions (Array<int> { head, n2, n1, tail });
      elementwise_constant = c1->ElementwiseConstant();
    }

    auto GetCArgs() const { return tuple { c1, head, n1, n2, tail }; }

    string GetDescription () const override
    { return "swap-axes"; }

    void CalcEquivalenceKey() override
    {
      this->equivalence_key = "swap(" + c1->EquivalenceKey() + "," + ToString(head) + "," + ToString(n1)
        + "," + ToString(n2) + "," + ToString(tail) + ")";
    }

    Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return { c1 }; }

    shared_ptr<CoefficientFunction>
    Transform(CoefficientFunction::T_Transform& transformation) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if(transformation.cache.count(thisptr))
        return transformation.cache[thisptr];
      if(transformation.replace.count(thisptr))
        return transformation.replace[thisptr];
      auto newcf = make_shared<SwapAxesCoefficientFunction>
        (c1->Transform(transformation), head, n1, n2, tail);
      transformation.cache[thisptr] = newcf;
      return newcf;
    }

    void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
    {
      auto dims = Dimensions();
      if (!code_uses_tensors)
        {
          for (int i : Range(Dimension()))
            code.body += "auto &" + Var(index, i, dims).S()
              + " = " + Var(inputs[0], Src(i), c1->Dimensions()).S() + ";\n";
          return;
        }
      code.Declare (index, dims, IsComplex());
      for (int i : Range(Dimension()))
        code.body += Var(index, i, dims).Assign(Var(inputs[0], Src(i), c1->Dimensions()), false);
    }

    void NonZeroPattern (const class ProxyUserData & ud,
                         FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                         FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      for (int i : Range(Dimension()))
        values(i) = input[0](Src(i));
    }

    using BASE::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("SwapAxesCF:: scalar evaluate for tensor called");
    }

    void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> result) const override
    {
      Vector<> va(c1->Dimension());
      c1->Evaluate (ip, va);
      for (int i : Range(Dimension()))
        result(i) = va(Src(i));
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,
                     BareSliceMatrix<T,ORD> values) const
    {
      auto va = input[0];
      size_t np = ir.Size();
      for (int i : Range(Dimension()))
        {
          int si = Src(i);
          for (size_t ip = 0; ip < np; ip++)
            values(i, ip) = va(si, ip);
        }
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir, BareSliceMatrix<T,ORD> values) const
    {
      auto & lh = TLHeap();
      HeapReset hr(lh);
      FlatMatrix<T,ORD> va(c1->Dimension(), mir.Size(), lh);
      c1->Evaluate (mir, va);
      BareSliceMatrix<T,ORD> in[1] = { va };
      T_Evaluate (mir, FlatArray<BareSliceMatrix<T,ORD>>(1, in), values);
    }

    Array<shared_ptr<CoefficientFunction>> InputVJP (shared_ptr<CoefficientFunction> w) const override
    {
      return { SwapAxesCF (w, head, n2, n1, tail, c1->Dimensions()) };
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == this) return dir;
      return SwapAxesCF (c1->Diff(var,dir), head, n1, n2, tail, Dimensions());
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];
      if (this == var)
        return IdentityCF(this->Dimensions());
      Array<int> dimres;
      dimres += Dimensions();
      dimres += var->Dimensions();
      auto res = SwapAxesCF (c1->DiffJacobi(var, cache), head, n1, n2, tail*var->Dimension(), dimres);
      cache[thisptr] = res;
      return res;
    }
  };

  shared_ptr<CoefficientFunction> SwapAxesCF (shared_ptr<CoefficientFunction> c,
                                              int head, int n1, int n2, int tail,
                                              FlatArray<int> dims)
  {
    if (n1 == 1 || n2 == 1)
      return c->Reshape(dims);
    if (c->IsZeroCF())
      return ZeroCF(dims);
    return make_shared<SwapAxesCoefficientFunction> (c, head, n1, n2, tail) -> Reshape(dims);
  }


  // c[i,j,A,B] = sum_k t[i,k,A] m[k,j,B]
  class MultTensorMatCoefficientFunction : public T_CoefficientFunction<MultTensorMatCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> c1, c2;
    int h, inner, w, ra, rb;
    using BASE = T_CoefficientFunction<MultTensorMatCoefficientFunction>;

    int IA (int i, int k, int a) const { return (i*inner+k)*ra+a; }
    int IB (int k, int j, int b) const { return (k*w+j)*rb+b; }
    int IC (int i, int j, int a, int b) const { return ((i*w+j)*ra+a)*rb+b; }
  public:
    MultTensorMatCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                      shared_ptr<CoefficientFunction> ac2)
      : BASE(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
    {
      auto d1 = c1->Dimensions();
      auto d2 = c2->Dimensions();
      if (d1.Size() < 2 || d2.Size() < 2 || d1[1] != d2[0])
        throw Exception("MultTensorMat: dimensions don't fit");
      h = d1[0]; inner = d1[1]; w = d2[1];
      ra = TailSize(d1, 2);
      rb = TailSize(d2, 2);
      Array<int> dims { h, w };
      AppendTail (dims, d1, 2);
      AppendTail (dims, d2, 2);
      SetDimensions (dims);
    }

    auto GetCArgs() const { return tuple { c1, c2 }; }

    string GetDescription () const override
    { return "tensor-matrix multiply"; }

    void CalcEquivalenceKey() override
    {
      this->equivalence_key = "(" + c1->EquivalenceKey() + "@@" + c2->EquivalenceKey() + ")";
    }

    Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return { c1, c2 }; }

    shared_ptr<CoefficientFunction>
    Transform(CoefficientFunction::T_Transform& transformation) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if(transformation.cache.count(thisptr))
        return transformation.cache[thisptr];
      if(transformation.replace.count(thisptr))
        return transformation.replace[thisptr];
      auto newcf = make_shared<MultTensorMatCoefficientFunction>
        (c1->Transform(transformation), c2->Transform(transformation));
      transformation.cache[thisptr] = newcf;
      return newcf;
    }

    void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
    {
      auto dims = Dimensions();
      code.Declare (index, dims, IsComplex());
      auto d1 = c1->Dimensions();
      auto d2 = c2->Dimensions();
      for (int i : Range(h))
        for (int j : Range(w))
          for (int a : Range(ra))
            for (int b : Range(rb))
              {
                CodeExpr s;
                for (int k : Range(inner))
                  s += Var(inputs[0], IA(i,k,a), d1) * Var(inputs[1], IB(k,j,b), d2);
                code.body += Var(index, IC(i,j,a,b), dims).Assign(s, false);
              }
    }

    void NonZeroPattern (const class ProxyUserData & ud,
                         FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                         FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      values = NonZero(false);
      for (int i : Range(h))
        for (int j : Range(w))
          for (int a : Range(ra))
            for (int b : Range(rb))
              for (int k : Range(inner))
                values(IC(i,j,a,b)) += input[0](IA(i,k,a)) * input[1](IB(k,j,b));
    }

    using BASE::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("MultTensorMatCF:: scalar evaluate for tensor called");
    }

    void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> result) const override
    {
      Vector<> va(c1->Dimension()), vb(c2->Dimension());
      c1->Evaluate (ip, va);
      c2->Evaluate (ip, vb);
      result = 0.0;
      for (int i : Range(h))
        for (int j : Range(w))
          for (int a : Range(ra))
            for (int b : Range(rb))
              for (int k : Range(inner))
                result(IC(i,j,a,b)) += va(IA(i,k,a)) * vb(IB(k,j,b));
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,
                     BareSliceMatrix<T,ORD> values) const
    {
      auto va = input[0];
      auto vb = input[1];
      size_t np = ir.Size();
      values.AddSize(Dimension(),np) = T(0.0);
      for (int i : Range(h))
        for (int j : Range(w))
          for (int a : Range(ra))
            for (int b : Range(rb))
              for (int k : Range(inner))
                {
                  auto row_a = va.Row(IA(i,k,a));
                  auto row_b = vb.Row(IB(k,j,b));
                  auto row_c = values.Row(IC(i,j,a,b));
                  for (size_t ip = 0; ip < np; ip++)
                    row_c(ip) += row_a(ip) * row_b(ip);
                }
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir, BareSliceMatrix<T,ORD> values) const
    {
      auto & lh = TLHeap();
      HeapReset hr(lh);
      FlatMatrix<T,ORD> va(c1->Dimension(), mir.Size(), lh);
      FlatMatrix<T,ORD> vb(c2->Dimension(), mir.Size(), lh);
      c1->Evaluate (mir, va);
      c2->Evaluate (mir, vb);
      BareSliceMatrix<T,ORD> in[2] = { va, vb };
      T_Evaluate (mir, FlatArray<BareSliceMatrix<T,ORD>>(2, in), values);
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == this) return dir;
      return MultTensorMatCF (c1->Diff(var,dir), c2) + MultTensorMatCF (c1, c2->Diff(var,dir));
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];
      if (this == var)
        return IdentityCF(this->Dimensions());
      Array<int> dimres;
      dimres += Dimensions();
      dimres += var->Dimensions();
      int nv = var->Dimension();
      auto dc1 = c1->DiffJacobi(var, cache);
      auto dc2 = c2->DiffJacobi(var, cache);
      shared_ptr<CoefficientFunction> res = ZeroCF(dimres);
      if (!dc1->IsZeroCF())   // [h,w,A,V,B] -> [h,w,A,B,V]
        res = SwapAxesCF (MultTensorMatCF (dc1->Reshape(Array<int>{h, inner, ra*nv}), c2),
                          h*w*ra, nv, rb, 1, dimres);
      if (!dc2->IsZeroCF())
        res = res + MultTensorMatCF (c1, dc2->Reshape(Array<int>{inner, w, rb*nv})) -> Reshape(dimres);
      cache[thisptr] = res;
      return res;
    }
  };

  shared_ptr<CoefficientFunction> MultTensorMatCF (shared_ptr<CoefficientFunction> t,
                                                   shared_ptr<CoefficientFunction> m)
  {
    if (t->IsZeroCF() || m->IsZeroCF())
      {
        Array<int> dims { t->Dimensions()[0], m->Dimensions()[1] };
        AppendTail (dims, t->Dimensions(), 2);
        AppendTail (dims, m->Dimensions(), 2);
        return ZeroCF (dims);
      }
    return make_shared<MultTensorMatCoefficientFunction> (t, m);
  }


  // tensor cross product, r[i,I,A,B] = eps_ijk eps_IJK a[j,J,A] b[k,K,B],  Cof(F) = 1/2 F x F
  class TensorCrossCoefficientFunction : public T_CoefficientFunction<TensorCrossCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> c1, c2;
    int ra, rb;
    using BASE = T_CoefficientFunction<TensorCrossCoefficientFunction>;
    using typename BASE::T_DJC;

    template <typename FUNC>
    static void Terms (int i, int ii, FUNC f)
    {
      for (int s1 : { 1, 2 })
        for (int s2 : { 1, 2 })
          f ((i+s1)%3, (i+3-s1)%3, (ii+s2)%3, (ii+3-s2)%3, s1==s2);
    }
    int IA (int j, int jj, int a) const { return (3*j+jj)*ra+a; }
    int IB (int k, int kk, int b) const { return (3*k+kk)*rb+b; }
    int IC (int i, int ii, int a, int b) const { return ((3*i+ii)*ra+a)*rb+b; }
  public:
    TensorCrossCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                    shared_ptr<CoefficientFunction> ac2)
      : BASE(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
    {
      auto d1 = c1->Dimensions();
      auto d2 = c2->Dimensions();
      if (d1.Size() < 2 || d1[0] != 3 || d1[1] != 3 || d2.Size() < 2 || d2[0] != 3 || d2[1] != 3)
        throw Exception("TensorCross: dimensions don't fit");
      ra = TailSize(d1, 2);
      rb = TailSize(d2, 2);
      Array<int> dims { 3, 3 };
      AppendTail (dims, d1, 2);
      AppendTail (dims, d2, 2);
      SetDimensions (dims);
    }

    auto GetCArgs() const { return tuple { c1, c2 }; }

    string GetDescription () const override
    { return "tensor-cross"; }

    void CalcEquivalenceKey() override
    {
      this->equivalence_key = "(" + c1->EquivalenceKey() + "#x#" + c2->EquivalenceKey() + ")";
    }

    Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return { c1, c2 }; }

    shared_ptr<CoefficientFunction>
    Transform(CoefficientFunction::T_Transform& transformation) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if(transformation.cache.count(thisptr))
        return transformation.cache[thisptr];
      if(transformation.replace.count(thisptr))
        return transformation.replace[thisptr];
      auto newcf = make_shared<TensorCrossCoefficientFunction>
        (c1->Transform(transformation), c2->Transform(transformation));
      transformation.cache[thisptr] = newcf;
      return newcf;
    }

    void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
    {
      auto dims = Dimensions();
      code.Declare (index, dims, IsComplex());
      auto d1 = c1->Dimensions();
      auto d2 = c2->Dimensions();
      for (int i : Range(3))
        for (int ii : Range(3))
          for (int a : Range(ra))
            for (int b : Range(rb))
              {
                CodeExpr s;
                Terms (i, ii, [&] (int j, int k, int jj, int kk, bool plus)
                       {
                         auto t = Var(inputs[0], IA(j,jj,a), d1) * Var(inputs[1], IB(k,kk,b), d2);
                         if (plus) s += t; else s -= t;
                       });
                code.body += Var(index, IC(i,ii,a,b), dims).Assign(s, false);
              }
    }

    void NonZeroPattern (const class ProxyUserData & ud,
                         FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                         FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      values = NonZero(false);
      for (int i : Range(3))
        for (int ii : Range(3))
          for (int a : Range(ra))
            for (int b : Range(rb))
              Terms (i, ii, [&] (int j, int k, int jj, int kk, bool plus)
                     { values(IC(i,ii,a,b)) += input[0](IA(j,jj,a)) * input[1](IB(k,kk,b)); });
    }

    using BASE::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("TensorCrossCF:: scalar evaluate for tensor called");
    }

    void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> result) const override
    {
      Vector<> va(c1->Dimension()), vb(c2->Dimension());
      c1->Evaluate (ip, va);
      c2->Evaluate (ip, vb);
      result = 0.0;
      for (int i : Range(3))
        for (int ii : Range(3))
          for (int a : Range(ra))
            for (int b : Range(rb))
              Terms (i, ii, [&] (int j, int k, int jj, int kk, bool plus)
                     {
                       double t = va(IA(j,jj,a)) * vb(IB(k,kk,b));
                       result(IC(i,ii,a,b)) += plus ? t : -t;
                     });
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,
                     BareSliceMatrix<T,ORD> values) const
    {
      auto va = input[0];
      auto vb = input[1];
      size_t np = ir.Size();
      values.AddSize(Dimension(),np) = T(0.0);
      for (int i : Range(3))
        for (int ii : Range(3))
          for (int a : Range(ra))
            for (int b : Range(rb))
              Terms (i, ii, [&] (int j, int k, int jj, int kk, bool plus)
                     {
                       auto row_a = va.Row(IA(j,jj,a));
                       auto row_b = vb.Row(IB(k,kk,b));
                       auto row_c = values.Row(IC(i,ii,a,b));
                       if (plus)
                         for (size_t ip = 0; ip < np; ip++)
                           row_c(ip) += row_a(ip) * row_b(ip);
                       else
                         for (size_t ip = 0; ip < np; ip++)
                           row_c(ip) -= row_a(ip) * row_b(ip);
                     });
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir, BareSliceMatrix<T,ORD> values) const
    {
      auto & lh = TLHeap();
      HeapReset hr(lh);
      FlatMatrix<T,ORD> va(c1->Dimension(), mir.Size(), lh);
      FlatMatrix<T,ORD> vb(c2->Dimension(), mir.Size(), lh);
      c1->Evaluate (mir, va);
      c2->Evaluate (mir, vb);
      BareSliceMatrix<T,ORD> in[2] = { va, vb };
      T_Evaluate (mir, FlatArray<BareSliceMatrix<T,ORD>>(2, in), values);
    }

    Array<shared_ptr<CoefficientFunction>> InputVJP (shared_ptr<CoefficientFunction> w) const override
    {
      if (ra != 1 || rb != 1) return Array<shared_ptr<CoefficientFunction>>();
      return { TensorCrossCF (w, c2), TensorCrossCF (w, c1) };
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == this) return dir;
      return TensorCrossCF (c1->Diff(var,dir), c2) + TensorCrossCF (c1, c2->Diff(var,dir));
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];
      if (this == var)
        return IdentityCF(this->Dimensions());
      Array<int> dimres;
      dimres += Dimensions();
      dimres += var->Dimensions();
      int nv = var->Dimension();
      auto dc1 = c1->DiffJacobi(var, cache);
      auto dc2 = c2->DiffJacobi(var, cache);
      shared_ptr<CoefficientFunction> res = ZeroCF(dimres);
      if (!dc1->IsZeroCF())   // [3,3,A,V,B] -> [3,3,A,B,V]
        res = SwapAxesCF (TensorCrossCF (dc1->Reshape(Array<int>{3, 3, ra*nv}), c2),
                          9*ra, nv, rb, 1, dimres);
      if (!dc2->IsZeroCF())
        res = res + TensorCrossCF (c1, dc2->Reshape(Array<int>{3, 3, rb*nv})) -> Reshape(dimres);
      cache[thisptr] = res;
      return res;
    }
  };

  shared_ptr<CoefficientFunction> TensorCrossCF (shared_ptr<CoefficientFunction> a,
                                                 shared_ptr<CoefficientFunction> b)
  {
    if (a->IsZeroCF() || b->IsZeroCF())
      {
        Array<int> dims { 3, 3 };
        AppendTail (dims, a->Dimensions(), 2);
        AppendTail (dims, b->Dimensions(), 2);
        return ZeroCF (dims);
      }
    return make_shared<TensorCrossCoefficientFunction> (a, b);
  }


  // r[i...,j...] = a[i...] b[j...]
  class TensorProductCoefficientFunction : public T_CoefficientFunction<TensorProductCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> c1, c2;
    int n1, n2;
    using BASE = T_CoefficientFunction<TensorProductCoefficientFunction>;
    using typename BASE::T_DJC;
  public:
    TensorProductCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                      shared_ptr<CoefficientFunction> ac2)
      : BASE(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
    {
      n1 = c1->Dimension();
      n2 = c2->Dimension();
      Array<int> dims;
      dims += c1->Dimensions();
      dims += c2->Dimensions();
      SetDimensions (dims);
    }

    auto GetCArgs() const { return tuple { c1, c2 }; }

    string GetDescription () const override
    { return "tensor-product"; }

    void CalcEquivalenceKey() override
    {
      this->equivalence_key = "(" + c1->EquivalenceKey() + "(x)" + c2->EquivalenceKey() + ")";
    }

    Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return { c1, c2 }; }

    shared_ptr<CoefficientFunction>
    Transform(CoefficientFunction::T_Transform& transformation) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if(transformation.cache.count(thisptr))
        return transformation.cache[thisptr];
      if(transformation.replace.count(thisptr))
        return transformation.replace[thisptr];
      auto newcf = make_shared<TensorProductCoefficientFunction>
        (c1->Transform(transformation), c2->Transform(transformation));
      transformation.cache[thisptr] = newcf;
      return newcf;
    }

    void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
    {
      auto dims = Dimensions();
      code.Declare (index, dims, IsComplex());
      auto d1 = c1->Dimensions();
      auto d2 = c2->Dimensions();
      for (int i : Range(n1))
        for (int j : Range(n2))
          code.body += Var(index, i*n2+j, dims).Assign(Var(inputs[0], i, d1) * Var(inputs[1], j, d2), false);
    }

    void NonZeroPattern (const class ProxyUserData & ud,
                         FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                         FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      for (int i : Range(n1))
        for (int j : Range(n2))
          values(i*n2+j) = input[0](i) * input[1](j);
    }

    using BASE::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      if (Dimension() != 1)
        throw Exception ("TensorProductCF:: scalar evaluate for tensor called");
      return c1->Evaluate(ip) * c2->Evaluate(ip);
    }

    void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> result) const override
    {
      Vector<> va(n1), vb(n2);
      c1->Evaluate (ip, va);
      c2->Evaluate (ip, vb);
      for (int i : Range(n1))
        for (int j : Range(n2))
          result(i*n2+j) = va(i) * vb(j);
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,
                     BareSliceMatrix<T,ORD> values) const
    {
      auto va = input[0];
      auto vb = input[1];
      size_t np = ir.Size();
      for (int i : Range(n1))
        for (int j : Range(n2))
          {
            auto row_a = va.Row(i);
            auto row_b = vb.Row(j);
            auto row_c = values.Row(i*n2+j);
            for (size_t ip = 0; ip < np; ip++)
              row_c(ip) = row_a(ip) * row_b(ip);
          }
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir, BareSliceMatrix<T,ORD> values) const
    {
      auto & lh = TLHeap();
      HeapReset hr(lh);
      FlatMatrix<T,ORD> va(n1, mir.Size(), lh);
      FlatMatrix<T,ORD> vb(n2, mir.Size(), lh);
      c1->Evaluate (mir, va);
      c2->Evaluate (mir, vb);
      BareSliceMatrix<T,ORD> in[2] = { va, vb };
      T_Evaluate (mir, FlatArray<BareSliceMatrix<T,ORD>>(2, in), values);
    }

    Array<shared_ptr<CoefficientFunction>> InputVJP (shared_ptr<CoefficientFunction> w) const override
    {
      auto wm = w->Reshape(n1, n2);
      return { (wm * c2->Reshape(n2)) -> Reshape(c1->Dimensions()),
               (c1->Reshape(1, n1) * wm) -> Reshape(c2->Dimensions()) };
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == this) return dir;
      return TensorProductCF (c1->Diff(var,dir), c2) + TensorProductCF (c1, c2->Diff(var,dir));
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];
      if (this == var)
        return IdentityCF(this->Dimensions());
      Array<int> dimres;
      dimres += Dimensions();
      dimres += var->Dimensions();
      int nv = var->Dimension();
      auto dc1 = c1->DiffJacobi(var, cache);
      auto dc2 = c2->DiffJacobi(var, cache);
      shared_ptr<CoefficientFunction> res = ZeroCF(dimres);
      if (!dc2->IsZeroCF())
        res = TensorProductCF (c1, dc2) -> Reshape(dimres);
      if (!dc1->IsZeroCF())
        res = res + MultTensorMatCF (dc1->Reshape(Array<int>{n1, 1, nv}), c2->Reshape(1, n2)) -> Reshape(dimres);
      cache[thisptr] = res;
      return res;
    }
  };

  shared_ptr<CoefficientFunction> TensorProductCF (shared_ptr<CoefficientFunction> a,
                                                   shared_ptr<CoefficientFunction> b)
  {
    Array<int> dims;
    dims += a->Dimensions();
    dims += b->Dimensions();
    if (a->IsZeroCF() || b->IsZeroCF())
      return ZeroCF (dims);
    if (a->Dimension() == 1 && a->Dimensions().Size() == 0)
      return (a * b) -> Reshape(dims);
    if (b->Dimension() == 1 && b->Dimensions().Size() == 0)
      return (b * a) -> Reshape(dims);
    return make_shared<TensorProductCoefficientFunction> (a, b);
  }


  // r[i,A,B] = eps_ijk a[j,A] b[k,B]
  class CrossTensorCoefficientFunction : public T_CoefficientFunction<CrossTensorCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> c1, c2;
    int ra, rb;
    using BASE = T_CoefficientFunction<CrossTensorCoefficientFunction>;
  public:
    CrossTensorCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                    shared_ptr<CoefficientFunction> ac2)
      : BASE(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
    {
      auto d1 = c1->Dimensions();
      auto d2 = c2->Dimensions();
      if (d1.Size() < 1 || d1[0] != 3 || d2.Size() < 1 || d2[0] != 3)
        throw Exception("CrossTensor: dimensions don't fit");
      ra = TailSize(d1, 1);
      rb = TailSize(d2, 1);
      Array<int> dims { 3 };
      AppendTail (dims, d1, 1);
      AppendTail (dims, d2, 1);
      SetDimensions (dims);
    }

    auto GetCArgs() const { return tuple { c1, c2 }; }

    string GetDescription () const override
    { return "cross-tensor"; }

    void CalcEquivalenceKey() override
    {
      this->equivalence_key = "(" + c1->EquivalenceKey() + "Xt" + c2->EquivalenceKey() + ")";
    }

    Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return { c1, c2 }; }

    shared_ptr<CoefficientFunction>
    Transform(CoefficientFunction::T_Transform& transformation) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if(transformation.cache.count(thisptr))
        return transformation.cache[thisptr];
      if(transformation.replace.count(thisptr))
        return transformation.replace[thisptr];
      auto newcf = make_shared<CrossTensorCoefficientFunction>
        (c1->Transform(transformation), c2->Transform(transformation));
      transformation.cache[thisptr] = newcf;
      return newcf;
    }

    void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
    {
      auto dims = Dimensions();
      code.Declare (index, dims, IsComplex());
      auto d1 = c1->Dimensions();
      auto d2 = c2->Dimensions();
      for (int i : Range(3))
        for (int a : Range(ra))
          for (int b : Range(rb))
            {
              int j = (i+1)%3, k = (i+2)%3;
              code.body += Var(index, (i*ra+a)*rb+b, dims).Assign
                (Var(inputs[0], j*ra+a, d1) * Var(inputs[1], k*rb+b, d2)
                 - Var(inputs[0], k*ra+a, d1) * Var(inputs[1], j*rb+b, d2), false);
            }
    }

    void NonZeroPattern (const class ProxyUserData & ud,
                         FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                         FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      for (int i : Range(3))
        for (int a : Range(ra))
          for (int b : Range(rb))
            {
              int j = (i+1)%3, k = (i+2)%3;
              values((i*ra+a)*rb+b) = input[0](j*ra+a) * input[1](k*rb+b) + input[0](k*ra+a) * input[1](j*rb+b);
            }
    }

    using BASE::Evaluate;
    double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("CrossTensorCF:: scalar evaluate for tensor called");
    }

    void Evaluate (const BaseMappedIntegrationPoint & ip, FlatVector<> result) const override
    {
      Vector<> va(c1->Dimension()), vb(c2->Dimension());
      c1->Evaluate (ip, va);
      c2->Evaluate (ip, vb);
      for (int i : Range(3))
        for (int a : Range(ra))
          for (int b : Range(rb))
            {
              int j = (i+1)%3, k = (i+2)%3;
              result((i*ra+a)*rb+b) = va(j*ra+a)*vb(k*rb+b) - va(k*ra+a)*vb(j*rb+b);
            }
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,
                     BareSliceMatrix<T,ORD> values) const
    {
      auto va = input[0];
      auto vb = input[1];
      size_t np = ir.Size();
      for (int i : Range(3))
        for (int a : Range(ra))
          for (int b : Range(rb))
            {
              int j = (i+1)%3, k = (i+2)%3;
              int ia = j*ra+a, ka = k*ra+a, kb = k*rb+b, jb = j*rb+b, ic = (i*ra+a)*rb+b;
              for (size_t ip = 0; ip < np; ip++)
                values(ic, ip) = va(ia,ip)*vb(kb,ip) - va(ka,ip)*vb(jb,ip);
            }
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir, BareSliceMatrix<T,ORD> values) const
    {
      auto & lh = TLHeap();
      HeapReset hr(lh);
      FlatMatrix<T,ORD> va(c1->Dimension(), mir.Size(), lh);
      FlatMatrix<T,ORD> vb(c2->Dimension(), mir.Size(), lh);
      c1->Evaluate (mir, va);
      c2->Evaluate (mir, vb);
      BareSliceMatrix<T,ORD> in[2] = { va, vb };
      T_Evaluate (mir, FlatArray<BareSliceMatrix<T,ORD>>(2, in), values);
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == this) return dir;
      return CrossTensorCF (c1->Diff(var,dir), c2) + CrossTensorCF (c1, c2->Diff(var,dir));
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];
      if (this == var)
        return IdentityCF(this->Dimensions());
      Array<int> dimres;
      dimres += Dimensions();
      dimres += var->Dimensions();
      int nv = var->Dimension();
      auto dc1 = c1->DiffJacobi(var, cache);
      auto dc2 = c2->DiffJacobi(var, cache);
      shared_ptr<CoefficientFunction> res = ZeroCF(dimres);
      if (!dc1->IsZeroCF())   // [3,A,V,B] -> [3,A,B,V]
        res = SwapAxesCF (CrossTensorCF (dc1->Reshape(3, ra*nv), c2), 3*ra, nv, rb, 1, dimres);
      if (!dc2->IsZeroCF())
        res = res + CrossTensorCF (c1, dc2->Reshape(3, rb*nv)) -> Reshape(dimres);
      cache[thisptr] = res;
      return res;
    }
  };

  shared_ptr<CoefficientFunction> CrossTensorCF (shared_ptr<CoefficientFunction> a,
                                                 shared_ptr<CoefficientFunction> b)
  {
    if (a->Dimension() == 3) a = a->Reshape(3);
    if (b->Dimension() == 3) b = b->Reshape(3);
    if (a->IsZeroCF() || b->IsZeroCF())
      {
        Array<int> dims { 3 };
        AppendTail (dims, a->Dimensions(), 1);
        AppendTail (dims, b->Dimensions(), 1);
        return ZeroCF (dims);
      }
    if (a->Dimensions().Size() == 1 && b->Dimensions().Size() == 1)
      return CrossProduct (a, b);
    return make_shared<CrossTensorCoefficientFunction> (a, b);
  }




  template <int D>
  class InverseCoefficientFunction : public T_CoefficientFunction<InverseCoefficientFunction<D>>
  {
    shared_ptr<CoefficientFunction> c1;
    using BASE = T_CoefficientFunction<InverseCoefficientFunction<D>>;
    using typename BASE::T_DJC;
  public:
    Array<shared_ptr<CoefficientFunction>> InputVJP (shared_ptr<CoefficientFunction> w) const override
    {
      // d A^{-1}: -A^{-T} w A^{-T}
      auto self = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      auto it = self->Transpose();
      return { -1.0 * (it * w * it) };
    }

    InverseCoefficientFunction() = default;
    InverseCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
      : T_CoefficientFunction<InverseCoefficientFunction>(D*D, ac1->IsComplex()), c1(ac1)
    {
      this->SetDimensions (ngstd::IVec<2> (D,D));
    }

    void DoArchive(Archive& ar) override
    {
      BASE::DoArchive(ar);
      ar.Shallow(c1);
    }

    virtual string GetDescription () const override
    { return "inverse"; }

  
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      // auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.res_type+">";
      auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.GetType(this->IsComplex())+">";
      auto mat_var = Var("mat", index);
      auto inv_var = Var("inv", index);
      code.body += mat_var.Declare(mat_type);
      code.body += inv_var.Declare(mat_type);
      for (int j = 0; j < D; j++)
        for (int k = 0; k < D; k++)
          code.body += mat_var(j,k).Assign(Var(inputs[0], j, k), false);

      code.body += inv_var.Assign(mat_var.Func("Inv"), false);

      for (int j = 0; j < D; j++)
        for (int k = 0; k < D; k++)
          code.body += code.AssignNew(index, j*D+k, this->Dimensions(), inv_var(j,k), this->IsComplex());
    }

    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

    /*
      virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
      FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
      {
      nonzero = true;
      nonzero_deriv = true;
      nonzero_dderiv = true;
      }
    */

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      auto v1 = input[0];
      AutoDiffDiff<1,NonZero> sum(false);
      for (int i = 0; i < v1.Size(); i++)
        sum += v1(i);
      values = sum;
      /*
        AutoDiffDiff<1,NonZero> add(true);
        add.DValue(0) = true;
        add.DDValue(0,0) = true;
        values = add;
      */
      /*
        FlatArray<int> hdims = Dimensions();    
        auto in0 = input[0];
        for (size_t j = 0; j < hdims[0]; j++)
        for (size_t k = 0; k < hdims[1]; k++)
        values(j*hdims[1]+k) = in0(k*hdims[0]+j);
      */
    }
    using T_CoefficientFunction<InverseCoefficientFunction<D>>::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("InverseCF:: scalar evaluate for matrix called");
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir,
                     BareSliceMatrix<T,ORD> result) const
    {
      c1->Evaluate (mir, result);
      for (size_t i = 0; i < mir.Size(); i++)
        {
          Mat<D,D,T> hm;
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              hm(j,k) = result(j*D+k, i);
          hm = Inv(hm);
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              result(j*D+k, i) = hm(j,k);
        }
    }  

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      size_t np = ir.Size();
      auto in0 = input[0];

      for (size_t i = 0; i < np; i++)
        {
          Mat<D,D,T> hm;
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              hm(j,k) = in0(j*D+k, i);
          hm = Inv(hm);
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              values(j*D+k, i) = hm(j,k);
        }
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (this == var) return dir;

      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      return (-1) * thisptr * c1->Diff(var,dir) * thisptr;
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];

      if (this == var)
        return IdentityCF(this->Dimensions());

      auto diffc1 = c1->DiffJacobi (var, cache);
      auto inv1 = thisptr;

      Array<int> dimres { D, D };
      dimres += var->Dimensions();
    
      auto prod1 = -inv1 * diffc1->Reshape( D, -1 );
      auto res = MultTensorMatCF (prod1 -> Reshape( Array<int>{D, D, int(var->Dimension())} ), inv1) -> Reshape( Array<int> (dimres) );
      cache[thisptr] = res;
      return res;
    }
  };


  class InverseCoefficientFunctionAnyDim : public CoefficientFunction
  {
    // This is implemented in a separate class because of limitations of msvc in the context of "constexpr if".
    // Otherwise, it could be easily integrated in InverseCoefficientFunction<D>, eg. for D == -1.

    shared_ptr<CoefficientFunction> c1;
    using T_DJC = CoefficientFunction::T_DJC;

  public:
    InverseCoefficientFunctionAnyDim() = default;
    InverseCoefficientFunctionAnyDim(shared_ptr<CoefficientFunction> ac1)
      : CoefficientFunction(ac1->Dimension() * ac1->Dimension(), ac1->IsComplex()), c1(ac1)
    {
      this->SetDimensions(ac1->Dimensions());
    }

    void DoArchive(Archive &ar) override
    {
      CoefficientFunction::DoArchive(ar);
      ar.Shallow(c1);
    }

    virtual Array<shared_ptr<CoefficientFunction>>
    InputCoefficientFunctions() const override
    {
      return Array<shared_ptr<CoefficientFunction>>({c1});
    }

    virtual void NonZeroPattern(const class ProxyUserData &ud,
                                FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      auto v1 = input[0];
      AutoDiffDiff<1,NonZero> sum(false);
      for (auto i : Range(v1))
        sum += v1(i);
      values = sum;
    }

    virtual double Evaluate(const BaseMappedIntegrationPoint &ip) const override
    {
      throw Exception("InverseAnyDimCF:: scalar evaluate for matrix called");
    }

    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<> result) const override
    {
      FlatMatrix<double,ColMajor> mat(Dimension(), 1, &result(0));
      ip.IntegrationRuleFromPoint([this,mat] (const BaseMappedIntegrationRule & ir)
      { this -> T_Evaluate (ir, BareSliceMatrix<double,ColMajor>(mat)); });
    }

    void Evaluate(const BaseMappedIntegrationPoint & ip,
                  FlatVector<Complex> result) const override
    {
      FlatMatrix<Complex,ColMajor> mat(Dimension(), 1, &result(0));
      ip.IntegrationRuleFromPoint([this,mat] (const BaseMappedIntegrationRule & ir)
      { this -> T_Evaluate (ir, BareSliceMatrix<Complex,ColMajor>(mat)); });
    }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<double> values) const override
    { this -> /* template */ T_Evaluate (ir, Trans(values)); }

    virtual void Evaluate (const BaseMappedIntegrationRule & ir,
                           FlatArray<BareSliceMatrix<double,ColMajor>> input,
                           BareSliceMatrix<double,ColMajor> values) const override
    { this -> T_Evaluate (ir, input, values); }


    virtual void Evaluate (const BaseMappedIntegrationRule & ir, BareSliceMatrix<Complex> values) const override
    {
      if (!IsComplex())
        {
          /*
            BareSliceMatrix<double> realvalues(2*values.Dist(), (double*)values.Data(),
            DummySize(values.Height(), values.Width()));
          */
          BareSliceMatrix<double> realvalues(values.Height(), values.Width(), 2*values.Dist(), (double*)values.Data());
          Evaluate (ir, realvalues);
          for (size_t i = 0; i < ir.Size(); i++)
            for (size_t j = Dimension(); j-- > 0; )
              values(i,j) = realvalues(i,j);
          return;
        }
      this -> /* template */ T_Evaluate (ir, Trans(values));
    }


    virtual void Evaluate (const BaseMappedIntegrationRule & mir, 
                           BareSliceMatrix<AutoDiff<1,double>> result) const override
    {
      c1->Evaluate(mir, result);
      
      // A^{-1}' = -A^{-1} A' A^{-1}
      const int D = c1->Dimensions()[0];      
      ArrayMem<double, 1000> mem(4*D*D);
      FlatMatrix<> hm(D, D, mem.Data());
      FlatMatrix<> hmp(D, D, mem.Data()+D*D);      
      FlatMatrix<> h1(D, D, mem.Data()+2*D*D);      
      FlatMatrix<> h2(D, D, mem.Data()+3*D*D);      

      for (auto i : Range(mir))
        {
          for (auto j : Range(D))
            for (auto k : Range(D))
              {
                hm(j, k) = result(i, j*D +k).Value();
                hmp(j, k) = result(i, j*D +k).DValue(0);
              }
          h1 = hmp*hm;
          h2 = hm*h1;
          for (auto j : Range(D))
            for (auto k : Range(D))
              result(i, j*D +k).DValue(0)=-h2(j,k);
        }
    }
    
    
    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate(const MIR &mir, BareSliceMatrix<T, ORD> result) const
    {
      c1->Evaluate(mir, result);

      const int D = c1->Dimensions()[0];
      ArrayMem<T, 1000> mem(D * D);
      FlatMatrix<T> hm(D, D, mem.Data());

      for (auto i : Range(mir))
        {
          for (auto j : Range(D))
            for (auto k : Range(D))
              hm(j, k) = result(j * D + k, i);
          // hm = Inv(hm);
          CalcInverse(hm);
          for (auto j : Range(D))
            for (auto k : Range(D))
              result(j * D + k, i) = hm(j, k);
        }
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate(const MIR &ir, FlatArray<BareSliceMatrix<T, ORD>> input,
                    BareSliceMatrix<T, ORD> values) const
    {
      auto in0 = input[0];

      const int D = c1->Dimensions()[0];
      ArrayMem<T, 1000> mem(D * D);
      FlatMatrix<T> hm(D, D, mem.Data());

      for (auto i : Range(ir))
        {
          for (auto j : Range(D))
            for (auto k : Range(D))
              hm(j, k) = in0(j * D + k, i);
          // hm = Inv(hm);
          CalcInverse(hm);
          for (auto j : Range(D))
            for (auto k : Range(D))
              values(j * D + k, i) = hm(j, k);
        }
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
    {
      if (code.is_simd)
        return CoefficientFunction::GenerateCode(code, inputs, index);

      const auto h = Dimensions()[0];
      const auto h_str = ToString(h);

      const auto mem_var = Var("mem", index);
      stringstream _code;
      _code << "ArrayMem<" << code.res_type << ", 400> " << mem_var.code << "("
            << Dimension() * 2 << ");"
            << "\n";
      code.body += _code.str();

      auto mat_var = Var("mat", index);
      auto inv_var = Var("inv", index);
      stringstream _code2;
      _code2 << "FlatMatrix<" + code.res_type + "> " << mat_var.code << "(" << h_str
             << ", " << mem_var.code << ".Data()"
             << ");"
             << "\n";
      _code2 << "FlatMatrix<" + code.res_type + "> " << inv_var.code << "(" << h_str
             << ", " << mem_var.code << ".Data() + " << Dimension() << ");"
             << "\n";

      code.body += _code2.str();
      for (auto j : Range(h))
        for (auto k : Range(h))
          code.body += mat_var(j, k).Assign(Var(inputs[0], j, k), false);

      code.body += inv_var.Assign(mat_var.Func("Inv"), false);

      for (auto j : Range(h))
        for (auto k : Range(h))
          code.body += code.AssignNew(index, j*h+k, this->Dimensions(), inv_var(j, k), this->IsComplex());
    }

    shared_ptr<CoefficientFunction>
    Diff(const CoefficientFunction *var,
         shared_ptr<CoefficientFunction> dir) const override
    {
      if (this == var)
        return c1->Diff(c1.get(), dir);

      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      return (-1) * thisptr * c1->Diff(var, dir) * thisptr;
    }

    shared_ptr<CoefficientFunction> DiffJacobi(const CoefficientFunction *var,
                                               T_DJC &cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];

      if (this == var)
        return IdentityCF(this->Dimensions());

      auto diffc1 = c1->DiffJacobi (var, cache);
      auto inv1 = thisptr;

      const int D = c1->Dimensions()[0];
      Array<int> dimres { D, D };
      dimres += var->Dimensions();

      auto prod1 = -inv1 * diffc1->Reshape( D, -1 );
      auto res = MultTensorMatCF (prod1 -> Reshape( Array<int>{D, D, int(var->Dimension())} ), inv1) -> Reshape( Array<int> (dimres) );
      cache[thisptr] = res;
      return res;
    }
  };





  template <int D>
  class DeterminantCoefficientFunction : public T_CoefficientFunction<DeterminantCoefficientFunction<D>>
  {
    shared_ptr<CoefficientFunction> c1;
    using BASE = T_CoefficientFunction<DeterminantCoefficientFunction<D>>;
  public:
    Array<shared_ptr<CoefficientFunction>> InputVJP (shared_ptr<CoefficientFunction> w) const override
    {
      return { w * CofactorCF(c1) };
    }

    DeterminantCoefficientFunction() = default;
    DeterminantCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
      : T_CoefficientFunction<DeterminantCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
    {
      auto dims_c1 = c1 -> Dimensions();
      if (dims_c1.Size() != 2)
        throw Exception("Determinant of non-matrix called");
      if (dims_c1[0] != dims_c1[1])
        throw Exception("Determinant of non-square matrix called");
    }

    virtual string GetDescription () const override
    { return "Determinant"; }

  
    void DoArchive(Archive& ar) override
    {
      BASE::DoArchive(ar);
      ar.Shallow(c1);
    }
  
    shared_ptr<CoefficientFunction>
    Transform(CoefficientFunction::T_Transform& transformation) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if(transformation.cache.count(thisptr))
        return transformation.cache[thisptr];
      if(transformation.replace.count(thisptr))
        return transformation.replace[thisptr];
      auto newcf = make_shared<DeterminantCoefficientFunction<D>>
        (c1->Transform(transformation));
      transformation.cache[thisptr] = newcf;
      return newcf;
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      // auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.res_type+">";
      auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.GetType(this->IsComplex())+">";
      auto mat_var = Var("mat", index);
      code.body += mat_var.Declare(mat_type);
      for (int j = 0; j < D; j++)
        for (int k = 0; k < D; k++)
          code.body += mat_var(j,k).Assign(Var(inputs[0], j, k), false);

      // code.Declare (code.res_type, index, this->Dimensions());
      code.Declare (index, this->Dimensions(), this->IsComplex());
    
      code.body += Var(index).Assign(mat_var.Func("Det"), false);
    }

    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

    /*
      virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
      FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
      {
      nonzero = true;
      nonzero_deriv = true;
      nonzero_dderiv = true;
      }
    */
  
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      auto sm = input[0];
      AutoDiffDiff<1,NonZero> res;
      switch (D)
        {
        case 1: 
          res = sm(0);
          break;
        case 2:
          res = sm(0)*sm(3)+sm(1)*sm(2);
          break;
        case 3:
          res = 
            sm(0) * (sm(4) * sm(8) + sm(5) * sm(7)) +
            sm(1) * (sm(5) * sm(6) + sm(3) * sm(8)) +
            sm(2) * (sm(3) * sm(7) + sm(4) * sm(6));
          break;
        default:
          {
            cerr << "general det not implemented" << endl;
          }
        }

      /*
        Mat<D,D,AutoDiffDiff<1,NonZero>> hm;
        for (int j = 0; j < D; j++)
        for (int k = 0; k < D; k++)
        hm(j,k) = in0(j*D+k);
        cout << "Nonzero mat:" << endl << hm << endl;
        values(0) = Det(hm);
        cout << "nonzero det: " << values(0) << endl;
      */
      values(0) = res;
    }
  
    using T_CoefficientFunction<DeterminantCoefficientFunction<D>>::Evaluate;
    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir,
                     BareSliceMatrix<T,ORD> result) const
    {
      STACK_ARRAY(T, hmem, mir.Size()*D*D);
      FlatMatrix<T,ORD> hv(D*D, mir.Size(), &hmem[0]);
      c1->Evaluate (mir, hv);
    
      for (size_t i = 0; i < mir.Size(); i++)
        {
          Mat<D,D,T> hm;
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              hm(j,k) = hv(j*D+k, i);
          result(0,i) = Det(hm);
        }
    }  

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      size_t np = ir.Size();
      auto in0 = input[0];

      for (size_t i = 0; i < np; i++)
        {
          Mat<D,D,T> hm;
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              hm(j,k) = in0(j*D+k, i);
          values(0,i) = Det(hm);        
        }
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (this == var) return dir;
      // return DeterminantCF(c1) * InnerProduct( TransposeCF(InverseCF(c1)), c1->Diff(var,dir) );
      return InnerProduct( CofactorCF(c1), c1->Diff(var,dir) );
    }

    shared_ptr<CoefficientFunction>
    DiffJacobi (const CoefficientFunction * var, typename BASE::T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];

      if (this == var) return make_shared<ConstantCoefficientFunction>(1);
      if (c1.get() == var) return CofactorCF (c1);
      auto input = c1->InputCoefficientFunctions();
      if (input.Size() == 0) return ZeroCF(var->Dimensions());

      auto cof = CofactorCF(c1) -> Reshape( 1, D*D );
      auto diffc1 = c1->DiffJacobi (var, cache) -> Reshape( D*D, var->Dimension() );
      auto prod = cof * diffc1;
      auto res = prod->Reshape(var->Dimensions());

      /*
        auto cof = CofactorCF(c1) -> Reshape( D*D );
        auto diffc1 = c1->DiffJacobi (var, cache) -> Reshape( D*D, var->Dimension() );
        auto res = diffc1->Transpose() * cof;
      */
    
      cache[thisptr] = res;
      return res;
    }

  
  };








  template <int D>
  class CofactorCoefficientFunction : public T_CoefficientFunction<CofactorCoefficientFunction<D>>
  {
    shared_ptr<CoefficientFunction> c1;
    using BASE = T_CoefficientFunction<CofactorCoefficientFunction<D>>;
  public:
    Array<shared_ptr<CoefficientFunction>> InputVJP (shared_ptr<CoefficientFunction> w) const override
    {
      if constexpr (D == 1)
        return { ZeroCF(c1->Dimensions()) };
      else if constexpr (D == 2)
        return { CofactorCF(w) };   // Cof is linear and self-adjoint in 2D
      else if constexpr (D == 3)
        {
          return { TensorCrossCF (w, c1) };
        }
      else
        return Array<shared_ptr<CoefficientFunction>>();
    }

    CofactorCoefficientFunction() = default;
    CofactorCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
      : T_CoefficientFunction<CofactorCoefficientFunction>(D*D, ac1->IsComplex()), c1(ac1)
    {
      this->SetDimensions (ngstd::IVec<2> (D,D));
    }

    void DoArchive(Archive& ar) override
    {
      BASE::DoArchive(ar);
      ar.Shallow(c1);
    }

    virtual string GetDescription () const override
    { return "cofactor"; }
  
    shared_ptr<CoefficientFunction>
    Transform(CoefficientFunction::T_Transform& transformation) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if(transformation.cache.count(thisptr))
        return transformation.cache[thisptr];
      if(transformation.replace.count(thisptr))
        return transformation.replace[thisptr];
      auto newcf = make_shared<CofactorCoefficientFunction<D>>(c1->Transform(transformation));
      transformation.cache[thisptr] = newcf;
      return newcf;
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      // auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.res_type+">";
      auto mat_type = "Mat<"+ToString(D)+","+ToString(D)+","+code.GetType(this->IsComplex())+">";
      auto mat_var = Var("mat", index);
      auto cof_var = Var("cof", index);
      code.body += mat_var.Declare(mat_type);
      code.body += cof_var.Declare(mat_type);
      for (int j = 0; j < D; j++)
        for (int k = 0; k < D; k++)
          code.body += mat_var(j,k).Assign(Var(inputs[0], j, k), false);

      code.body += cof_var.Assign(mat_var.Func("Cof"), false);
    
      // code.Declare (code.res_type, index, this->Dimensions());
      code.Declare (index, this->Dimensions(), this->IsComplex()); 
      for (int j = 0; j < D; j++)
        for (int k = 0; k < D; k++)
          code.body += Var(index, j, k).Assign(cof_var(j,k), false);
    }

    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

    /*
      virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
      FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
      {
      nonzero = true;
      nonzero_deriv = true;
      nonzero_dderiv = true;
      }
    */

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      auto v1 = input[0];

      Mat<D,D,AutoDiffDiff<1,NonZero>> v1mat;
      for (int i = 0; i < D; i++)
        for (int j = 0; j < D; j++)
          v1mat(i,j) = v1(i*D+j);
      auto v2mat = Cof(v1mat);
      for (int i = 0; i < D; i++)
        for (int j = 0; j < D; j++)
          values(i*D+j) = v2mat(i,j);
    }
  
    using T_CoefficientFunction<CofactorCoefficientFunction<D>>::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("CofactorCF:: scalar evaluate for matrix called");
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir,
                     BareSliceMatrix<T,ORD> result) const
    {
      c1->Evaluate (mir, result);
      for (size_t i = 0; i < mir.Size(); i++)
        {
          Mat<D,D,T> hm;
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              hm(j,k) = result(j*D+k, i);
          hm = Cof(hm);
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              result(j*D+k, i) = hm(j,k);
        }
    }  

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      size_t np = ir.Size();
      auto in0 = input[0];

      for (size_t i = 0; i < np; i++)
        {
          Mat<D,D,T> hm;
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              hm(j,k) = in0(j*D+k, i);
          hm = Cof(hm);
          for (int j = 0; j < D; j++)
            for (int k = 0; k < D; k++)
              values(j*D+k, i) = hm(j,k);
        }
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (this == var) return dir;
      if (this->Dimensions()[0] <= 2)
        {
          //Cofactor Matrix linear in 2d (in 1d Cofactor Matrix = 0)
          return CofactorCF(c1->Diff(var,dir));
        }
      else if (this->Dimensions()[0] == 3) //3d
        {
          // Cof(A) = 1/2 A x A
          return TensorCrossCF (c1->Diff(var,dir), c1);
        }
      else
        throw Exception("CofactorCF Diff only implemented for dim <=3");
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, typename BASE::T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];

      if (this == var)
        return IdentityCF(this->Dimensions());

      shared_ptr<CoefficientFunction> res;

      if (this->Dimensions()[0] == 2)
        res = (TraceCF(c1)*IdentityCF(2)-TransposeCF(c1)) -> DiffJacobi (var, cache);
      else if (this->Dimensions()[0] == 3)
        {
          res = TensorCrossCF (c1->DiffJacobi (var, cache), c1);
        }
      else
        res = (DeterminantCF(c1) * InverseCF(c1)->Transpose() ) -> DiffJacobi (var, cache);
      cache[thisptr] = res;
      return res;
    }

  
  };



  class SymmetricCoefficientFunction : public T_CoefficientFunction<SymmetricCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> c1;
    using BASE = T_CoefficientFunction<SymmetricCoefficientFunction>;
  public:
    SymmetricCoefficientFunction() = default;
    SymmetricCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
      : T_CoefficientFunction<SymmetricCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
    {
      auto dims_c1 = c1 -> Dimensions();
      if (dims_c1.Size() != 2)
        throw Exception("Sym of non-matrix called");
      if (dims_c1[0] != dims_c1[1])
        throw Exception("Sym of non-square matrix called");
    
      SetDimensions (ngstd::IVec<2> (dims_c1[0], dims_c1[0]) );
    }

    void DoArchive(Archive& ar) override
    {
      BASE::DoArchive(ar);
      ar.Shallow(c1);
    }
  
    virtual string GetDescription () const override
    { return "symmetric"; }

    shared_ptr<CoefficientFunction>
    Transform(CoefficientFunction::T_Transform& transformation) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if(transformation.cache.count(thisptr))
        return transformation.cache[thisptr];
      if(transformation.replace.count(thisptr))
        return transformation.replace[thisptr];
      auto newcf = SymmetricCF(c1->Transform(transformation));
      transformation.cache[thisptr] = newcf;
      return newcf;
    }


    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      FlatArray<int> hdims = Dimensions();        
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1]))
          code.body += code.AssignNew(index, i*hdims[1]+j, hdims, "0.5*("+Var(inputs[0],i,j).S()+"+"+Var(inputs[0],j,i).S()+")", IsComplex());
    }

    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      int hd = Dimensions()[0];    
      auto in0 = input[0];
      for (int i = 0; i < hd; i++)
        for (int j = 0; j < hd; j++)
          {
            int ii = i*hd+j;
            int jj = j*hd+i;
            values(ii) = in0(ii)+in0(jj);   // logical or 
          }
    }
    using T_CoefficientFunction<SymmetricCoefficientFunction>::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("SymCF:: scalar evaluate for matrix called");
    }
  
    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir,
                     BareSliceMatrix<T,ORD> result) const
    {
      int hd = Dimensions()[0];
      c1->Evaluate (mir, result);
      STACK_ARRAY(T, hmem, hd*hd);
      FlatMatrix<T,ORD> tmp (hd, hd, &hmem[0]);

      for (size_t i = 0; i < mir.Size(); i++)
        {
          for (int j = 0; j < hd; j++)
            for (int k = 0; k < hd; k++)
              tmp(j,k) = result(k*hd+j, i);
          for (int j = 0; j < hd; j++)
            for (int k = 0; k < hd; k++)
              result(j*hd+k, i) = 0.5*(tmp(j,k)+tmp(k,j));
        }
    }  

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      int hd = Dimensions()[0];
      size_t np = ir.Size();
    
      auto in0 = input[0];
      for (size_t j = 0; j < hd; j++)
        for (size_t k = 0; k < hd; k++)
          for (size_t i = 0; i < np; i++)
            values(j*hd+k, i) = 0.5 * (in0(k*hd+j, i)+in0(j*hd+k, i));
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (this == var) return dir;
      return SymmetricCF(c1->Diff(var, dir));
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];

      if (this == var)
        return IdentityCF(this->Dimensions());
  
      auto diffc1 = c1->DiffJacobi (var, cache);
      auto res = 0.5*(diffc1 + diffc1 -> TensorTranspose( 0, 1 ));
      cache[thisptr] = res;
      return res;
    }
  };


  class SkewCoefficientFunction : public T_CoefficientFunction<SkewCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> c1;
    using BASE = T_CoefficientFunction<SkewCoefficientFunction>;
  public:
    SkewCoefficientFunction() = default;
    SkewCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
      : T_CoefficientFunction<SkewCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
    {
      auto dims_c1 = c1 -> Dimensions();
      if (dims_c1.Size() != 2)
        throw Exception("Skew of non-matrix called");
      if (dims_c1[0] != dims_c1[1])
        throw Exception("Skew of non-square matrix called");
    
      SetDimensions (IVec<2> (dims_c1[0], dims_c1[0]) );
    }

    void DoArchive(Archive& ar) override
    {
      BASE::DoArchive(ar);
      ar.Shallow(c1);
    }
  
    virtual string GetDescription () const override
    { return "skew"; }
  
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      FlatArray<int> hdims = Dimensions();        
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1]))
          code.body += code.AssignNew(index, i*hdims[1]+j, hdims, "0.5*("+Var(inputs[0],i,j).S()+"-"+Var(inputs[0],j,i).S()+")", IsComplex());
    }

    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

    /*
      virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
      FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
      {
      int hd = Dimensions()[0];    
      c1->NonZeroPattern (ud, nonzero, nonzero_deriv, nonzero_dderiv);
      for (int i = 0; i < hd; i++)
      for (int j = 0; j < hd; j++)
      {
      int ii = i*hd+j;
      int jj = j*hd+i;
      nonzero(ii) |= nonzero(jj);
      nonzero_deriv(ii) |= nonzero_deriv(jj);
      nonzero_dderiv(ii) |= nonzero_dderiv(jj);
      }
      }
    */

  
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      int hd = Dimensions()[0];    
      auto in0 = input[0];
      for (int i = 0; i < hd; i++)
        for (int j = 0; j < hd; j++)
          {
            int ii = i*hd+j;
            int jj = j*hd+i;
            values(ii) = in0(ii)+in0(jj);   // logical or 
          }
    }
    using T_CoefficientFunction<SkewCoefficientFunction>::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("SkewCF:: scalar evaluate for matrix called");
    }

  
    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir,
                     BareSliceMatrix<T,ORD> result) const
    {
      int hd = Dimensions()[0];
      c1->Evaluate (mir, result);
      STACK_ARRAY(T, hmem, hd*hd);
      FlatMatrix<T,ORD> tmp (hd, hd, &hmem[0]);

      for (size_t i = 0; i < mir.Size(); i++)
        {
          for (int j = 0; j < hd; j++)
            for (int k = 0; k < hd; k++)
              tmp(j,k) = result(j*hd+k, i);
          for (int j = 0; j < hd; j++)
            for (int k = 0; k < hd; k++)
              result(j*hd+k, i) = 0.5*(tmp(j,k)-tmp(k,j));
        }
    }  

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      int hd = Dimensions()[0];
      size_t np = ir.Size();
    
      auto in0 = input[0];
      for (size_t j = 0; j < hd; j++)
        for (size_t k = 0; k < hd; k++)
          for (size_t i = 0; i < np; i++)
            values(j*hd+k, i) = 0.5 * (in0(j*hd+k, i)-in0(k*hd+j, i));
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (this == var) return dir;
      return SkewCF(c1->Diff(var, dir));
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];

      if (this == var)
        return IdentityCF(this->Dimensions());
  
      auto diffc1 = c1->DiffJacobi (var, cache);
      auto res = 0.5*(diffc1 - diffc1 -> TensorTranspose( 0, 1 ));
      cache[thisptr] = res;
      return res;
    }
  };




  

  shared_ptr<CoefficientFunction> IdentityCF (int dim)
  {
     return make_shared<IdentityCoefficientFunction> (dim);
  }

  shared_ptr<CoefficientFunction> IdentityCF (FlatArray<int> dims)
  {

    if (dims.Size() == 0)
      return ConstantCF(1);

    int dim = 1;
    for (auto d : dims)
      dim *= d;
    
    Array<int> tensor_dims;
    tensor_dims.Append(dims);
    tensor_dims.Append(dims);
    return make_shared<IdentityCoefficientFunction>(dim)->Reshape(tensor_dims);
  }


  shared_ptr<CoefficientFunction> InverseCF (shared_ptr<CoefficientFunction> coef)
  {
    auto dims = coef->Dimensions();
    if (dims.Size() != 2) throw Exception("Inverse of non-matrix");
    if (dims[0] != dims[1]) throw Exception("Inverse of non-quadratic matrix");
    switch (dims[0])
      {
      case 1: return make_shared<InverseCoefficientFunction<1>> (coef);
      case 2: return make_shared<InverseCoefficientFunction<2>> (coef);
      case 3: return make_shared<InverseCoefficientFunction<3>> (coef);
      default:
        return make_shared<InverseCoefficientFunctionAnyDim>(coef);
      }
  }



  shared_ptr<CoefficientFunction> DeterminantCF (shared_ptr<CoefficientFunction> coef)
  {
    auto dims = coef->Dimensions();
    if (dims.Size() != 2) throw Exception("Inverse of non-matrix");
    if (dims[0] != dims[1]) throw Exception("Inverse of non-quadratic matrix");

    if (coef->IsZeroCF())
      return ZeroCF(Array<int>());

    if (dynamic_pointer_cast<IdentityCoefficientFunction> (coef) && !coef->IsVariable())
      return make_shared<ConstantCoefficientFunction>(1);


    // common pattern : Det (F^T F) = Det(F)**2, for F square
    if (!coef->IsVariable())
      if (auto mmm = dynamic_pointer_cast<MultMatMatCoefficientFunction> (coef))
        {
          auto AB = mmm->InputCoefficientFunctions();
          if (!AB[0]->IsVariable())
            if (auto trans = dynamic_pointer_cast<TransposeCoefficientFunction> (AB[0]))
              {
                auto At = trans->InputCoefficientFunctions()[0];
                if (At->Dimensions()[0] == At->Dimensions()[1])
                  {
                    auto detF = DeterminantCF (At);
                    return detF*detF;
                  }
              }
        }
    
    switch (dims[0])
      {
      case 1: return make_shared<DeterminantCoefficientFunction<1>> (coef);
      case 2: return make_shared<DeterminantCoefficientFunction<2>> (coef);
      case 3: return make_shared<DeterminantCoefficientFunction<3>> (coef);
      default:
        throw Exception("Determinant of matrix of size "+ToString(dims[0]) + " not available");
      }
  }

  shared_ptr<CoefficientFunction> CofactorCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      return coef;
    
    auto dims = coef->Dimensions();
    if (dims.Size() != 2) throw Exception("Cofactor of non-matrix");
    if (dims[0] != dims[1]) throw Exception("Cofactor of non-quadratic matrix");
    switch (dims[0])
      {
      case 1: return make_shared<CofactorCoefficientFunction<1>> (coef);
      case 2: return make_shared<CofactorCoefficientFunction<2>> (coef);
      case 3: return make_shared<CofactorCoefficientFunction<3>> (coef);
      case 4: return make_shared<CofactorCoefficientFunction<4>> (coef);
      default:
        throw Exception("Cofactor of matrix of size "+ToString(dims[0]) + " not available");
      }
  }


  shared_ptr<CoefficientFunction> TraceCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      return ZeroCF(Array<int>());
    
    return make_shared<TraceCoefficientFunction> (coef);
  }


  shared_ptr<CoefficientFunction> SymmetricCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      return coef;

    // Sym(Grad(u)) -> u.Operator("strain"), if the space provides it
    if (auto proxy = dynamic_pointer_cast<ProxyFunction> (coef))
      if (auto prim = dynamic_pointer_cast<ProxyFunction> (proxy->Primary()))
        if (prim->GetAdditionalEvaluator("strain") &&
            (proxy == prim->Deriv() || proxy == prim->GetAdditionalProxy("Grad")))
          return prim->GetAdditionalProxy("strain");

    return make_shared<SymmetricCoefficientFunction> (coef);
  }

  shared_ptr<CoefficientFunction> SkewCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      return coef;

    return make_shared<SkewCoefficientFunction> (coef);
  }


  shared_ptr<CoefficientFunction>
  MakeMultDiagMatCoefficientFunction (shared_ptr<CoefficientFunction> c1,
                                  shared_ptr<CoefficientFunction> c2)
  {
    return make_shared<MultDiagMatCoefficientFunction> (c1,c2);
  }

  
  
  template class T_CoefficientFunction<MultMatMatCoefficientFunction>;
  template class T_CoefficientFunction<MultDiagMatCoefficientFunction>;  
  template class T_CoefficientFunction<TransposeCoefficientFunction>;
  template class T_CoefficientFunction<TraceCoefficientFunction>;
  template class T_CoefficientFunction<IdentityCoefficientFunction, CoefficientFunctionNoDerivative>;    
  
  static RegisterClassForArchive<TransposeCoefficientFunction, CoefficientFunction> regtransposecf;
  static RegisterClassForArchive<SymmetricCoefficientFunction, CoefficientFunction> regsymmetriccf;
  static RegisterClassForArchive<SkewCoefficientFunction, CoefficientFunction> regskewcf;
  static RegisterClassForArchive<TraceCoefficientFunction, CoefficientFunction> regtracecf;
  static RegisterClassForArchive<IdentityCoefficientFunction, CoefficientFunction> regidentitycf;
  
  static RegisterClassForArchive<MultMatMatCoefficientFunction, CoefficientFunction> regmultmatmatcf;
  static RegisterClassForArchive<MultDiagMatCoefficientFunction, CoefficientFunction> regmultdiagmatcf;  
  static RegisterClassForArchive<SwapAxesCoefficientFunction, CoefficientFunction> regswapaxescf;
  static RegisterClassForArchive<MultTensorMatCoefficientFunction, CoefficientFunction> regmulttensormatcf;
  static RegisterClassForArchive<TensorCrossCoefficientFunction, CoefficientFunction> regtensorcrosscf;
  static RegisterClassForArchive<TensorProductCoefficientFunction, CoefficientFunction> regtensorproductcf;
  static RegisterClassForArchive<CrossTensorCoefficientFunction, CoefficientFunction> regcrosstensorcf;

  static RegisterClassForArchive<InverseCoefficientFunction<1>, CoefficientFunction> reginversecf1;
  static RegisterClassForArchive<InverseCoefficientFunction<2>, CoefficientFunction> reginversecf2;
  static RegisterClassForArchive<InverseCoefficientFunction<3>, CoefficientFunction> reginversecf3;
  static RegisterClassForArchive<InverseCoefficientFunctionAnyDim, CoefficientFunction> reginverseanydimcf;

  static RegisterClassForArchive<DeterminantCoefficientFunction<1>, CoefficientFunction> regdetcf1;
  static RegisterClassForArchive<DeterminantCoefficientFunction<2>, CoefficientFunction> regdetcf2;
  static RegisterClassForArchive<DeterminantCoefficientFunction<3>, CoefficientFunction> regdetcf3;
  static RegisterClassForArchive<CofactorCoefficientFunction<1>, CoefficientFunction> regcof1;
  static RegisterClassForArchive<CofactorCoefficientFunction<2>, CoefficientFunction> regcof2;
  static RegisterClassForArchive<CofactorCoefficientFunction<3>, CoefficientFunction> regcof3;
  
}
