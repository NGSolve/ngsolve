
#include <core/register_archive.hpp>
#include <ngstd.hpp>
#include <bla.hpp>

#include <../basiclinalg/complex_wrapper.hpp>
// #include <fem.hpp>
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
  TransposeCoefficientFunction :: ~TransposeCoefficientFunction() { }
  TraceCoefficientFunction :: ~TraceCoefficientFunction() { }    
  IdentityCoefficientFunction :: ~IdentityCoefficientFunction () {  }




  template <int D>
  class InverseCoefficientFunction : public T_CoefficientFunction<InverseCoefficientFunction<D>>
  {
    shared_ptr<CoefficientFunction> c1;
    using BASE = T_CoefficientFunction<InverseCoefficientFunction<D>>;
    using typename BASE::T_DJC;
  public:
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

  
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      c1->TraverseTree (func);
      func(*this);
    }

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
          code.body += Var(index, j, k).Assign(inv_var(j,k));
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
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      Vector<AutoDiffDiff<1,NonZero>> v1(c1->Dimension());
      c1->NonZeroPattern (ud, v1);
      AutoDiffDiff<1,NonZero> sum(false);
      for (int i = 0; i < v1.Size(); i++)
        sum += v1(i);
      values = sum;
    }
  
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
      auto prod1r = prod1 -> Reshape( Array<int> (dimres) );
      auto trans = prod1r -> TensorTranspose( 0, 1 );
      auto prod2 = inv1->Transpose() * trans -> Reshape( D, -1 );
      auto res = prod2 -> Reshape( Array<int> (dimres) ) -> TensorTranspose( 0, 1);
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

    virtual void TraverseTree(const function<void(CoefficientFunction &)> &func) override
    {
      c1->TraverseTree(func);
      func(*this);
    }

    virtual Array<shared_ptr<CoefficientFunction>>
    InputCoefficientFunctions() const override
    {
      return Array<shared_ptr<CoefficientFunction>>({c1});
    }

    virtual void NonZeroPattern(const class ProxyUserData &ud,
                                FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      Vector<AutoDiffDiff<1,NonZero>> v1(c1->Dimension());
      c1->NonZeroPattern(ud, v1);
      AutoDiffDiff<1,NonZero> sum(false);
      for (auto i : Range(v1))
        sum += v1(i);
      values = sum;
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
          code.body += Var(index, j, k).Assign(inv_var(j, k));
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
      auto prod1r = prod1 -> Reshape( Array<int> (dimres) );
      auto trans = prod1r -> TensorTranspose( 0, 1 );
      auto prod2 = inv1->Transpose() * trans -> Reshape( D, -1 );
      auto res = prod2 -> Reshape( Array<int> (dimres) ) -> TensorTranspose( 0, 1);
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
  
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      c1->TraverseTree (func);
      func(*this);
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
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      /*
        AutoDiffDiff<1,NonZero> add(true);
        add.DValue(0) = true;
        add.DDValue(0,0) = true;
        values = add;
      */
      Vector<AutoDiffDiff<1,NonZero>> in(D*D);
      c1->NonZeroPattern (ud, in);
      Array<FlatVector<AutoDiffDiff<1,NonZero>>> input{1UL};
      input[0].AssignMemory(D*D, &in(0));
      NonZeroPattern (ud, input, values);
    }
  
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
  
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      c1->TraverseTree (func);
      func(*this);
    }

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
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      Vector<AutoDiffDiff<1,NonZero>> v1(c1->Dimension());
      c1->NonZeroPattern (ud, v1);
      /*
        AutoDiffDiff<1,NonZero> sum(false);
        for (int i = 0; i < v1.Size(); i++)
        sum += v1(i);
        values = sum;
      */
      Mat<D,D,AutoDiffDiff<1,NonZero>> v1mat;
      for (int i = 0; i < D; i++)
        for (int j = 0; j < D; j++)
          v1mat(i,j) = v1(i*D+j);
      auto v2mat = Cof(v1mat);
      for (int i = 0; i < D; i++)
        for (int j = 0; j < D; j++)
          values(i*D+j) = v2mat(i,j);
    }
  
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
          //formula follows from Cayley–Hamilton
          //Cof(A) = 0.5*(tr(A)**2 - tr(A**2))I - tr(A)A^T +(AA)^T

          //return (0.5*(TraceCF(c1)*TraceCF(c1) - TraceCF(c1*c1))*IdentityCF(3) - TraceCF(c1)*TransposeCF(c1) + TransposeCF(c1*c1))->Diff(var,dir);
          auto trace_c1 = TraceCF(c1);
          auto diff_c1 = c1->Diff(var,dir);
          auto trace_diff_c1 = TraceCF(diff_c1);
          auto diff_c1_x_c1 = diff_c1 * c1;
          auto c1_x_diff_c1 = c1 * diff_c1;
          return (trace_c1 * trace_diff_c1 - TraceCF(diff_c1_x_c1)) * IdentityCF(3)
            - trace_diff_c1 * TransposeCF(c1)
            - trace_c1 * TransposeCF(diff_c1) + TransposeCF(diff_c1_x_c1 + c1_x_diff_c1);
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
          auto trcf = TraceCF(c1);
          auto prodcf = c1*c1;
          res = (0.5*(trcf*trcf - TraceCF(prodcf))*IdentityCF(3) - trcf*TransposeCF(c1) + TransposeCF(prodcf)) -> DiffJacobi (var, cache);
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
  
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      c1->TraverseTree (func);
      func(*this);
    }

    virtual string GetDescription () const override
    { return "symmetric"; }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      FlatArray<int> hdims = Dimensions();        
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1]))
          code.body += Var(index,i,j).Assign("0.5*("+Var(inputs[0],i,j).S()+"+"+Var(inputs[0],j,i).S()+")");
    }

    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      int hd = Dimensions()[0];
      c1->NonZeroPattern (ud, values);

      for (int i = 0; i < hd; i++)
        for (int j = 0; j < hd; j++)
          {
            int ii = i*hd+j;
            int jj = j*hd+i;
            values(ii) = values(ii)+values(jj);  // logical or
          }
    }
  
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
  
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      c1->TraverseTree (func);
      func(*this);
    }

    virtual string GetDescription () const override
    { return "skew"; }
  
    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      FlatArray<int> hdims = Dimensions();        
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1]))
          code.body += Var(index,i,j).Assign("0.5*("+Var(inputs[0],i,j).S()+"-"+Var(inputs[0],j,i).S()+")");
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
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      int hd = Dimensions()[0];
      c1->NonZeroPattern (ud, values);
      for (int i = 0; i < hd; i++)
        for (int j = 0; j < hd; j++)
          {
            int ii = i*hd+j;
            int jj = j*hd+i;
            values(ii) = values(ii)+values(jj);   // logical or 
          }
    }

  
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

    return make_shared<SymmetricCoefficientFunction> (coef);
  }

  shared_ptr<CoefficientFunction> SkewCF (shared_ptr<CoefficientFunction> coef)
  {
    if (coef->IsZeroCF())
      return coef;

    return make_shared<SkewCoefficientFunction> (coef);
  }

  
  
  template class T_CoefficientFunction<MultMatMatCoefficientFunction>;
  template class T_CoefficientFunction<TransposeCoefficientFunction>;
  template class T_CoefficientFunction<TraceCoefficientFunction>;
  template class T_CoefficientFunction<IdentityCoefficientFunction>;    
  
  static RegisterClassForArchive<TransposeCoefficientFunction, CoefficientFunction> regtransposecf;
  static RegisterClassForArchive<SymmetricCoefficientFunction, CoefficientFunction> regsymmetriccf;
  static RegisterClassForArchive<SkewCoefficientFunction, CoefficientFunction> regskewcf;
  static RegisterClassForArchive<TraceCoefficientFunction, CoefficientFunction> regtracecf;
  static RegisterClassForArchive<IdentityCoefficientFunction, CoefficientFunction> regidentitycf;
  
  static RegisterClassForArchive<MultMatMatCoefficientFunction, CoefficientFunction> regmultmatmatcf;

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
