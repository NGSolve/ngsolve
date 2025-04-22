#ifndef COEFFICIENT_MATRIX
#define COEFFICIENT_MATRIX

namespace ngfem
{


  class IdentityCoefficientFunction : public T_CoefficientFunction<IdentityCoefficientFunction>
  {
    using BASE = T_CoefficientFunction<IdentityCoefficientFunction>;
  public:
    IdentityCoefficientFunction (int dim)
      : T_CoefficientFunction<IdentityCoefficientFunction>(1, false)
    {
      SetDimensions (ngstd::IVec<2> (dim, dim) );
    }

    // For archive
    // IdentityCoefficientFunction() = default;
    virtual ~IdentityCoefficientFunction ();

    auto GetCArgs() const { return tuple { Dimensions()[0] }; }
    void DoArchive(Archive& ar) override { /* BASE::DoArchive(ar); */ }

    virtual string GetDescription () const override
    { return "Identity matrix"; }
  
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      func(*this);
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override
    {
      FlatArray<int> hdims = Dimensions();
      // code.Declare (code.res_type, index, this->Dimensions());
      code.Declare (index, this->Dimensions(), this->IsComplex());

      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1]))
          {
            if (i == j)
              code.body += Var(index,i,j).Assign(string("1.0"), false);
            else
              code.body += Var(index,i,j).Assign(string("0.0"), false);
          }
    }

  
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      int hd = Dimensions()[0];
      values = AutoDiffDiff<1,NonZero>(false);

      for (int i = 0; i < hd; i++)
        values(i*(hd+1)) = AutoDiffDiff<1,NonZero>(true);
    }
  
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      int hd = Dimensions()[0];
      values = AutoDiffDiff<1,NonZero>(false);

      for (int i = 0; i < hd; i++)
        values(i*(hd+1)) = AutoDiffDiff<1,NonZero>(true);
    }

    using T_CoefficientFunction<IdentityCoefficientFunction>::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("IdentityCF:: scalar evaluate for matrix called");
    }
  
    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir,
                     BareSliceMatrix<T,ORD> result) const
    {
      result.AddSize(Dimension(), mir.Size()) = T(0.0);
      int hd = Dimensions()[0];

      for (size_t i = 0; i < mir.Size(); i++)
        for (int j = 0; j < hd; j++)
          result(j*(hd+1), i) = T(1.0);
    }  

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      int hd = Dimensions()[0];
      size_t np = ir.Size();
      values.AddSize(Dimension(), np) = T(0.0);
    
      for (size_t j = 0; j < hd; j++)
        for (size_t i = 0; i < np; i++)
          values(j*(hd+1), i) = T(1.0);
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (this == var) return dir;
      return ZeroCF(this->Dimensions());
    }

    shared_ptr<CoefficientFunction>
    DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      if (this == var)
        return IdentityCF(this->Dimensions());

      int dim = Dimensions()[0];
      Array<int> resdims = { dim, dim };
      resdims += var->Dimensions();
      return ZeroCF( resdims );
    }

  
  };



  
  class MultMatMatCoefficientFunction : public T_CoefficientFunction<MultMatMatCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> c1;
    shared_ptr<CoefficientFunction> c2;
    int inner_dim;
    using BASE = T_CoefficientFunction<MultMatMatCoefficientFunction>;
  public:
    // MultMatMatCoefficientFunction() = default;
    MultMatMatCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                   shared_ptr<CoefficientFunction> ac2)
      : T_CoefficientFunction<MultMatMatCoefficientFunction>(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
    {
      auto dims_c1 = c1 -> Dimensions();
      auto dims_c2 = c2 -> Dimensions();
      if (dims_c1.Size() != 2 || dims_c2.Size() != 2)
        throw Exception("Mult of non-matrices called");
      if (dims_c1[1] != dims_c2[0])
        throw Exception(string("Matrix dimensions don't fit: m1 is ") +
                        ToLiteral(dims_c1[0]) + " x " + ToLiteral(dims_c1[1]) +
                        ", m2 is " + ToLiteral(dims_c2[0]) + " x " + ToLiteral(dims_c2[1]) );
      SetDimensions( ngstd::IVec<2> (dims_c1[0], dims_c2[1]) );
      inner_dim = dims_c1[1];
    }

    // auto GetCArgs() const { return tuple { Shallow(c1), Shallow(c2) }; }
    auto GetCArgs() const { return tuple { c1, c2 }; }
    
    virtual ~MultMatMatCoefficientFunction();
    virtual string GetDescription () const override
    { return "matrix-matrix multiply"; }
  
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      c1->TraverseTree (func);
      c2->TraverseTree (func);
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
      auto newcf = make_shared<MultMatMatCoefficientFunction>
        (c1->Transform(transformation),
         c2->Transform(transformation));
      transformation.cache[thisptr] = newcf;
      return newcf;
    }

    void DoArchive(Archive& ar) override
    {
      /*
      BASE::DoArchive(ar);
      ar.Shallow(c1).Shallow(c2) & inner_dim;
      cout << "archive mamat, c1 = " << c1 << ", c2 = " << c2 << endl;
      */
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      FlatArray<int> hdims = Dimensions();
      // code.Declare (code.res_type, index, hdims);
      code.Declare (index, hdims, IsComplex());

      if (code_uses_tensors)
        {
          code.body += "for (size_t i = 0; i < "+ToString(hdims[0])+"; i++)\n";
          code.body += "for (size_t j = 0; j < "+ToString(hdims[1])+"; j++) { \n";
          code.body += "auto sum = var_" + ToString(inputs[0]) + "(i,0) * var_" + ToString(inputs[1]) + "(0,j); \n";
          code.body += "for (size_t k = 1; k < "+ToString(inner_dim)+"; k++) \n";
          code.body += "sum += var_" + ToString(inputs[0]) + "(i,k) * var_" + ToString(inputs[1]) + "(k,j); \n";
          code.body += "var_" + ToString(index) + "(i,j) = sum; } \n";
        }
      else
        {
          for (int i : Range(hdims[0]))
            for (int j : Range(hdims[1])) {
              CodeExpr s;
              for (int k : Range(inner_dim))
                s += Var(inputs[0], i, k) * Var(inputs[1], k, j);
              code.body += Var(index, i, j).Assign(s, false);
            }
        }
    }

    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ c1, c2 }); }  


    /*
      virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
      FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
      {
      FlatArray<int> hdims = Dimensions();
      Vector<bool> v1(hdims[0]*inner_dim), v2(hdims[1]*inner_dim);
      Vector<bool> d1(hdims[0]*inner_dim), d2(hdims[1]*inner_dim);
      Vector<bool> dd1(hdims[0]*inner_dim), dd2(hdims[1]*inner_dim);
      c1->NonZeroPattern (ud, v1, d1, dd1);
      c2->NonZeroPattern (ud, v2, d2, dd2);
      nonzero = false;
      nonzero_deriv = false;
      nonzero_dderiv = false;
      FlatMatrix<bool> m1(hdims[0], inner_dim, &v1(0));
      FlatMatrix<bool> m2(inner_dim, hdims[1], &v2(0));
      FlatMatrix<bool> md1(hdims[0], inner_dim, &d1(0));
      FlatMatrix<bool> md2(inner_dim, hdims[1], &d2(0));
      FlatMatrix<bool> mdd1(hdims[0], inner_dim, &dd1(0));
      FlatMatrix<bool> mdd2(inner_dim, hdims[1], &dd2(0));

      for (int i = 0; i < hdims[0]; i++)
      for (int j = 0; j < hdims[1]; j++)
      for (int k = 0; k < inner_dim; k++)
      {
      nonzero(i*hdims[1]+j) |= m1(i,k) && m2(k,j);
      nonzero_deriv(i*hdims[1]+j) |= (m1(i,k) && md2(k,j)) || (md1(i,k) && m2(k,j));
      nonzero_dderiv(i*hdims[1]+j) |= (m1(i,k) && mdd2(k,j)) || (md1(i,k) && md2(k,j)) || (mdd1(i,k) && m2(k,j));
      }
      }
    */

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      FlatArray<int> hdims = Dimensions();
      Vector<AutoDiffDiff<1,NonZero>> va(hdims[0]*inner_dim), vb(hdims[1]*inner_dim);
      c1->NonZeroPattern (ud, va);
      c2->NonZeroPattern (ud, vb);
    
      size_t d1 = hdims[1];

      values = NonZero(false);
    
      for (size_t j = 0; j < hdims[0]; j++)
        for (size_t k = 0; k < hdims[1]; k++)
          for (size_t l = 0; l < inner_dim; l++)
            values(j*d1+k) += va(j*inner_dim+l) * vb(l*d1+k);
    }

  
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      auto va = input[0];
      auto vb = input[1];

      FlatArray<int> hdims = Dimensions();    
      size_t d1 = hdims[1];

      values = NonZero(false);
    
      for (size_t j = 0; j < hdims[0]; j++)
        for (size_t k = 0; k < hdims[1]; k++)
          for (size_t l = 0; l < inner_dim; l++)
            values(j*d1+k) += va(j*inner_dim+l) * vb(l*d1+k);
    }

    using T_CoefficientFunction<MultMatMatCoefficientFunction>::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("MultMatMatCF:: scalar evaluate for matrix called");
    }

    virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                           FlatVector<> result) const override
    {
      FlatArray<int> hdims = Dimensions();
      Vector<> va(hdims[0]*inner_dim);
      Vector<> vb(hdims[1]*inner_dim);
      FlatMatrix<> a(hdims[0], inner_dim, &va[0]);
      FlatMatrix<> b(inner_dim, hdims[1], &vb[0]);
    
      c1->Evaluate (ip, va);
      c2->Evaluate (ip, vb);

      FlatMatrix<> c(hdims[0], hdims[1], &result(0));
      c = a*b;
    }  

    virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                           FlatVector<Complex> result) const override
    {
      FlatArray<int> hdims = Dimensions();
      STACK_ARRAY(double,mema,2*hdims[0]*inner_dim);
      STACK_ARRAY(double,memb,2*hdims[1]*inner_dim);
      FlatVector<Complex> va(hdims[0]*inner_dim, reinterpret_cast<Complex*>(&mema[0]));
      FlatVector<Complex> vb(inner_dim*hdims[1], reinterpret_cast<Complex*>(&memb[0]));
    
      c1->Evaluate (ip, va);
      c2->Evaluate (ip, vb);

      FlatMatrix<Complex> a(hdims[0], inner_dim, &va(0));
      FlatMatrix<Complex> b(inner_dim, hdims[1], &vb(0));

      FlatMatrix<Complex> c(hdims[0], hdims[1], &result(0));
      c = a*b;
      //cout << "MultMatMat: complex not implemented" << endl;
    }  

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir, BareSliceMatrix<T,ORD> values) const 
    {
      FlatArray<int> hdims = Dimensions();    
      STACK_ARRAY(T, hmem1, mir.Size()*hdims[0]*inner_dim);
      STACK_ARRAY(T, hmem2, mir.Size()*hdims[1]*inner_dim);
      FlatMatrix<T,ORD> va(hdims[0]*inner_dim, mir.Size(), &hmem1[0]);
      FlatMatrix<T,ORD> vb(hdims[1]*inner_dim, mir.Size(), &hmem2[0]);

      c1->Evaluate (mir, va);
      c2->Evaluate (mir, vb);

      values.AddSize(Dimension(),mir.Size()) = T(0.0);

      size_t d1 = hdims[1];
      size_t mir_size = mir.Size();
      for (size_t j = 0; j < hdims[0]; j++)
        for (size_t k = 0; k < hdims[1]; k++)
          for (size_t l = 0; l < inner_dim; l++)
            {
              auto row_a = va.Row(j*inner_dim+l);
              auto row_b = vb.Row(l*d1+k);
              auto row_c = values.Row(j*d1+k);
              for (size_t i = 0; i < mir_size; i++)
                row_c(i) += row_a(i) * row_b(i);
              // row_c = pw_mult (row_a, row_b);
            }
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      auto va = input[0];
      auto vb = input[1];

      FlatArray<int> hdims = Dimensions();    
      size_t d1 = hdims[1];
      size_t np = ir.Size();

      values.AddSize(Dimension(),np) = T(0.0);
    
      for (size_t j = 0; j < hdims[0]; j++)
        for (size_t k = 0; k < hdims[1]; k++)
          for (size_t l = 0; l < inner_dim; l++)
            {
              auto row_a = va.Row(j*inner_dim+l);
              auto row_b = vb.Row(l*d1+k);
              auto row_c = values.Row(j*d1+k);
              for (size_t i = 0; i < np; i++)
                row_c(i) += row_a(i) * row_b(i);
              // row_c = pw_mult (row_a, row_b);
            }
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == this) return dir;
      return c1->Diff(var,dir)*c2 + c1 * c2->Diff(var,dir);
    }
  
    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {    
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];

      int h = Dimensions()[0];
      int w = Dimensions()[1];
      int dimvar = var->Dimension();

      if (this == var)
        return IdentityCF(this->Dimensions());

      Array<int> dimres{h,w};
      dimres += var->Dimensions();

      auto diffc1 = c1->DiffJacobi (var, cache);
      auto diffc2 = c2->DiffJacobi (var, cache);
    
      auto diffc1_trans = diffc1 -> TensorTranspose( 0, 1 ) -> Reshape( inner_dim, h*dimvar );
      auto prod1 = c2->Transpose() * diffc1_trans;
      Array<int> dimtmp{w, h};
      dimtmp += var->Dimensions();
      auto prod1trans = prod1->Reshape(dimtmp) -> TensorTranspose( 0, 1 );
    
      auto diffc2_trans = diffc2 -> Reshape( Array<int>{inner_dim,w*dimvar} );
      auto prod2 = c1 * diffc2_trans;
      auto prod2trans = prod2->Reshape(dimres);

      auto res = prod1trans + prod2trans;
      cache[thisptr] = res;
      return res;
    }
  };



  /*
    Diagonal-matrix(vector c1) * tensor c2
    plan:
    C = A*B in sense of
    C_ijk = A_ij B_ik
    with i.. single index, and j,k multi-index
   */
  
  class MultDiagMatCoefficientFunction : public T_CoefficientFunction<MultDiagMatCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> c1; 
    shared_ptr<CoefficientFunction> c2;
    size_t numcols;
    using BASE = T_CoefficientFunction<MultDiagMatCoefficientFunction>;
  public:
    MultDiagMatCoefficientFunction (shared_ptr<CoefficientFunction> ac1,
                                    shared_ptr<CoefficientFunction> ac2)
      : T_CoefficientFunction<MultDiagMatCoefficientFunction>(1, ac1->IsComplex()||ac2->IsComplex()), c1(ac1), c2(ac2)
    {
      auto dims_c1 = c1 -> Dimensions();
      auto dims_c2 = c2 -> Dimensions();
      if (dims_c1.Size() != 1)
        throw Exception("MultDiagMat: first argument must be vector");
      if (dims_c2.Size() < 1)
        throw Exception("MultDiagMat: second argument must be tensor");
      if (dims_c1[0] != dims_c2[0])
        throw Exception(string("MultDiagMat dimensions don't fit"));
      
      SetDimensions( dims_c2 );
      numcols = 1;
      for (auto d : dims_c2.Range(1,END))
        numcols *= d;
    }

    auto GetCArgs() const { return tuple { c1, c2 }; }
    
    virtual ~MultDiagMatCoefficientFunction();
    virtual string GetDescription () const override
    { return "diagmat-matrix multiply"; }
  
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      c1->TraverseTree (func);
      c2->TraverseTree (func);
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
      auto newcf = make_shared<MultDiagMatCoefficientFunction>
        (c1->Transform(transformation),
         c2->Transform(transformation));
      transformation.cache[thisptr] = newcf;
      return newcf;
    }

    void DoArchive(Archive& ar) override
    {} 

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      FlatArray<int> hdims = Dimensions();
      code.Declare (index, hdims, IsComplex());
      
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1])) {
          CodeExpr s = Var(inputs[0], i) * Var(inputs[1], i, j);
          code.body += Var(index, i, j).Assign(s, false);
        }
    }

    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ c1, c2 }); }  

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      size_t dim1 = c1->Dimension();
      Vector<AutoDiffDiff<1,NonZero>> va(dim1);
      c1->NonZeroPattern (ud, va);
      c2->NonZeroPattern (ud, values);

      for (size_t i = 0, ii=0; i < dim1; i++)
        for (size_t j = 0; j < numcols; j++, ii++)
          values(ii) *= va(i);
    }

  
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      auto va = input[0];
      auto vb = input[1];

      size_t dim1 = c1->Dimension();
      for (size_t i = 0, ii=0; i < dim1; i++)
        for (size_t j = 0; j < numcols; j++, ii++)
          values(ii) = va(i)*vb(ii);
    }

    using T_CoefficientFunction<MultDiagMatCoefficientFunction>::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("MultDiagMatCF:: scalar evaluate for matrix called");
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir, BareSliceMatrix<T,ORD> values) const 
    {
      size_t dim1 = c1->Dimension();
      STACK_ARRAY(T, hmem1, mir.Size()*dim1);
      FlatMatrix<T,ORD> va(dim1, mir.Size(), &hmem1[0]);

      c1->Evaluate (mir, va);
      c2->Evaluate (mir, values);

      for (size_t i = 0, ii=0; i < dim1; i++)
        for (size_t j = 0; j < numcols; j++, ii++)
          for (size_t k = 0; k < mir.Size(); k++)
            values(ii,k) *= va(i,k);
    }

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      auto va = input[0];
      auto vb = input[1];

      size_t dim1 = c1->Dimension();      
      for (size_t i = 0, ii=0; i < dim1; i++)
        for (size_t j = 0; j < numcols; j++, ii++)
          for (size_t k = 0; k < ir.Size(); k++)
            values(ii,k) = va(i,k) * vb(ii,k);
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (var == this) return dir;
      return make_shared<MultDiagMatCoefficientFunction>(c1->Diff(var,dir), c2) +
        make_shared<MultDiagMatCoefficientFunction>(c1, c2->Diff(var,dir));
    }

    /*
      TODO
    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {    
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];

      int h = Dimensions()[0];
      int w = Dimensions()[1];
      int dimvar = var->Dimension();

      if (this == var)
        return IdentityCF(this->Dimensions());

      Array<int> dimres{h,w};
      dimres += var->Dimensions();

      auto diffc1 = c1->DiffJacobi (var, cache);
      auto diffc2 = c2->DiffJacobi (var, cache);
    
      auto diffc1_trans = diffc1 -> TensorTranspose( 0, 1 ) -> Reshape( inner_dim, h*dimvar );
      auto prod1 = c2->Transpose() * diffc1_trans;
      Array<int> dimtmp{w, h};
      dimtmp += var->Dimensions();
      auto prod1trans = prod1->Reshape(dimtmp) -> TensorTranspose( 0, 1 );
    
      auto diffc2_trans = diffc2 -> Reshape( Array<int>{inner_dim,w*dimvar} );
      auto prod2 = c1 * diffc2_trans;
      auto prod2trans = prod2->Reshape(dimres);

      auto res = prod1trans + prod2trans;
      cache[thisptr] = res;
      return res;
    }
    */
  };




  

  class TransposeCoefficientFunction : public T_CoefficientFunction<TransposeCoefficientFunction>
  {
    shared_ptr<CoefficientFunction> c1;
    using BASE = T_CoefficientFunction<TransposeCoefficientFunction>;
  public:
    // TransposeCoefficientFunction() = default;
    TransposeCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
      : T_CoefficientFunction<TransposeCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
    {
      auto dims_c1 = c1 -> Dimensions();
      if (dims_c1.Size() != 2)
        throw Exception("Transpose of non-matrix called");

      SetDimensions (ngstd::IVec<2> (dims_c1[1], dims_c1[0]) );
    }

    virtual ~TransposeCoefficientFunction();

    shared_ptr<CoefficientFunction>
    Transform(CoefficientFunction::T_Transform& transformation) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if(transformation.cache.count(thisptr))
        return transformation.cache[thisptr];
      if(transformation.replace.count(thisptr))
        return transformation.replace[thisptr];
      auto newcf = make_shared<TransposeCoefficientFunction>
        (c1->Transform(transformation));
      transformation.cache[thisptr] = newcf;
      return newcf;
    }

    auto GetCArgs() const { return tuple { c1 }; }    
    void DoArchive(Archive& ar) override
    {
      /*
      BASE::DoArchive(ar);
      ar.Shallow(c1);
      */
    }

    virtual string GetDescription () const override
    { return "Matrix transpose"; }
  
    virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
    {
      c1->TraverseTree (func);
      func(*this);
    }

    virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
      FlatArray<int> hdims = Dimensions();
      // code.Declare (code.res_type, index, this->Dimensions());
      code.Declare (index, this->Dimensions(), this->IsComplex());
      
      for (int i : Range(hdims[0]))
        for (int j : Range(hdims[1]))
          code.body += Var(index,i,j).Assign( Var(inputs[0],j,i), false );
    }

    virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
    { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

    /*
      virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
      FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
      {
      FlatArray<int> hdims = Dimensions();    
      Vector<bool> v1(hdims[0]*hdims[1]), d1(hdims[0]*hdims[1]), dd1(hdims[0]*hdims[1]);
      c1->NonZeroPattern (ud, v1, d1, dd1);
      {
      FlatMatrix<bool> m1(hdims[1], hdims[0], &v1(0));
      FlatMatrix<bool> m2(hdims[0], hdims[1], &nonzero(0));
      m2 = Trans(m1);
      }
      {
      FlatMatrix<bool> m1(hdims[1], hdims[0], &d1(0));
      FlatMatrix<bool> m2(hdims[0], hdims[1], &nonzero_deriv(0));
      m2 = Trans(m1);
      }
      {
      FlatMatrix<bool> m1(hdims[1], hdims[0], &dd1(0));
      FlatMatrix<bool> m2(hdims[0], hdims[1], &nonzero_dderiv(0));
      m2 = Trans(m1);
      }
      }
    */
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      FlatArray<int> hdims = Dimensions();    
      Vector<AutoDiffDiff<1,NonZero>> v1(hdims[0]*hdims[1]);
      c1->NonZeroPattern (ud, v1);

      for (size_t j = 0; j < hdims[0]; j++)
        for (size_t k = 0; k < hdims[1]; k++)
          values(j*hdims[1]+k) = v1(k*hdims[0]+j);
    }

    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                                 FlatVector<AutoDiffDiff<1,NonZero>> values) const override
    {
      FlatArray<int> hdims = Dimensions();    
      auto in0 = input[0];
      for (size_t j = 0; j < hdims[0]; j++)
        for (size_t k = 0; k < hdims[1]; k++)
          values(j*hdims[1]+k) = in0(k*hdims[0]+j);
    }
    using T_CoefficientFunction<TransposeCoefficientFunction>::Evaluate;
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const override
    {
      throw Exception ("TransposeCF:: scalar evaluate for matrix called");
    }

    virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                           FlatVector<> result) const override
    {
      FlatArray<int> hdims = Dimensions();        
      VectorMem<20> input(result.Size());
      c1->Evaluate (ip, input);    
      FlatMatrix<> reshape1(hdims[1], hdims[0], &input(0));  // source matrix format
      FlatMatrix<> reshape2(hdims[0], hdims[1], &result(0));  // range matrix format
      reshape2 = Trans(reshape1);
    
      /*
        c1->Evaluate (ip, result);
        static Timer t("Transpose - evaluate");
        RegionTimer reg(t);
        FlatMatrix<> reshape(dims[1], dims[0], &result(0));  // source matrix format
        Matrix<> tmp = Trans(reshape);
        FlatMatrix<> reshape2(dims[0], dims[1], &result(0));  // range matrix format
        reshape2 = tmp;
      */
    }  

    virtual void Evaluate (const BaseMappedIntegrationPoint & ip,
                           FlatVector<Complex> result) const override
    {
      FlatArray<int> hdims = Dimensions();        
      STACK_ARRAY(double,meminput,2*hdims[0]*hdims[1]);
      FlatVector<Complex> input(hdims[0]*hdims[1],reinterpret_cast<Complex*>(&meminput[0]));
      c1->Evaluate (ip, input);    
      FlatMatrix<Complex> reshape1(hdims[1], hdims[0], &input(0));  // source matrix format
      FlatMatrix<Complex> reshape2(hdims[0], hdims[1], &result(0));  // range matrix format
      reshape2 = Trans(reshape1);
      //cout << "Transpose: complex not implemented" << endl;
    }  

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & mir,
                     BareSliceMatrix<T,ORD> result) const
    {
      FlatArray<int> hdims = Dimensions();    
      c1->Evaluate (mir, result);
      STACK_ARRAY(T, hmem, hdims[0]*hdims[1]);
      FlatMatrix<T,ORD> tmp (hdims[0], hdims[1], &hmem[0]);

      for (size_t i = 0; i < mir.Size(); i++)
        {
          for (int j = 0; j < hdims[0]; j++)
            for (int k = 0; k < hdims[1]; k++)
              tmp(j,k) = result(k*hdims[0]+j, i);
          for (int j = 0; j < hdims[0]; j++)
            for (int k = 0; k < hdims[1]; k++)
              result(j*hdims[1]+k, i) = tmp(j,k);
        }
    }  

    template <typename MIR, typename T, ORDERING ORD>
    void T_Evaluate (const MIR & ir,
                     FlatArray<BareSliceMatrix<T,ORD>> input,                       
                     BareSliceMatrix<T,ORD> values) const
    {
      FlatArray<int> hdims = Dimensions();
      size_t np = ir.Size();
    
      auto in0 = input[0];
      for (size_t j = 0; j < hdims[0]; j++)
        for (size_t k = 0; k < hdims[1]; k++)
          for (size_t i = 0; i < np; i++)
            values(j*hdims[1]+k, i) = in0(k*hdims[0]+j, i);
    }

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                          shared_ptr<CoefficientFunction> dir) const override
    {
      if (this == var) return dir;
      return TransposeCF (c1->Diff(var, dir));
    }

    shared_ptr<CoefficientFunction> DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
    {
      auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
      if (cache.find(thisptr) != cache.end())
        return cache[thisptr];

      if (this == var)
        return IdentityCF(this->Dimensions());
  
      auto diffc1 = c1->DiffJacobi (var, cache);
      auto res = diffc1 -> TensorTranspose( 0, 1 );
      cache[thisptr] = res;
      return res;
    }
  };




class TraceCoefficientFunction : public T_CoefficientFunction<TraceCoefficientFunction>
{
  shared_ptr<CoefficientFunction> c1;
  using BASE = T_CoefficientFunction<TraceCoefficientFunction>;
public:
  // TraceCoefficientFunction() = default;
  TraceCoefficientFunction (shared_ptr<CoefficientFunction> ac1)
    : T_CoefficientFunction<TraceCoefficientFunction>(1, ac1->IsComplex()), c1(ac1)
  {
    auto dims_c1 = c1 -> Dimensions();
    if (dims_c1.Size() != 2)
      throw Exception("Trace of non-matrix called");
    if (dims_c1[0] != dims_c1[1])
      throw Exception("Trace of non-square matrix called");
  }

  ~TraceCoefficientFunction();
  
  virtual string GetDescription () const override
  { return "trace"; }

  auto GetCArgs() const { return tuple { c1 }; }      
  void DoArchive(Archive& ar) override
  {
    /*
    BASE::DoArchive(ar);
    ar.Shallow(c1);
    */
  }
  
  virtual void TraverseTree (const function<void(CoefficientFunction&)> & func) override
  {
    c1->TraverseTree (func);
    func(*this);
  }

  virtual void GenerateCode(Code &code, FlatArray<int> inputs, int index) const override {
    CodeExpr result;
    int dim1 = c1->Dimensions()[0];
    code.Declare (index, Array<int>(), IsComplex()); 
    for (int i = 0; i < dim1; i++)
      result += Var(inputs[0],i,i);
    code.body += Var(index).Assign(result.S(), false);
  }
  
  virtual Array<shared_ptr<CoefficientFunction>> InputCoefficientFunctions() const override
  { return Array<shared_ptr<CoefficientFunction>>({ c1 } ); }  

  /*
  virtual void NonZeroPattern (const class ProxyUserData & ud, FlatVector<bool> nonzero,
                               FlatVector<bool> nonzero_deriv, FlatVector<bool> nonzero_dderiv) const override
  {
    int dim1 = c1->Dimension();
    Vector<bool> v1(dim1), d1(dim1), dd1(dim1);
    c1->NonZeroPattern (ud, v1, d1, dd1);

    bool v = false, d = false, dd = false;
    for (int i = 0; i < dim1; i++)
      {
        v |= v1(i);
        d |= d1(i);
        dd |= dd1(i);
      }
    nonzero = v;
    nonzero_deriv = d;
    nonzero_dderiv = dd;
  }
  */
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatVector<AutoDiffDiff<1,NonZero>> values) const override
  {
    int dim1 = c1->Dimension();
    int d = c1->Dimensions()[0];
    Vector<AutoDiffDiff<1,NonZero>> v1(dim1);
    c1->NonZeroPattern (ud, v1);
    values(0) = false;
    /*
    for (int i = 0; i < dim1; i++)
      values(0) = values(0)+v1(i);   // logical or
    */
    for (int i = 0; i < d; i++)
      values(0) = values(0)+v1(i*(d+1));   // logical or
  }
  
  virtual void NonZeroPattern (const class ProxyUserData & ud,
                               FlatArray<FlatVector<AutoDiffDiff<1,NonZero>>> input,
                               FlatVector<AutoDiffDiff<1,NonZero>> values) const override
  {
    // int dim1 = c1->Dimension();
    int d = c1->Dimensions()[0];
    auto in0 = input[0];
    values(0) = false;
    /*
    for (int i = 0; i < dim1; i++)
      values(0) = values(0)+in0(i);   // logical or
    */
    for (int i = 0; i < d; i++)
      values(0) = values(0)+in0(i*(d+1));   // logical or
  }

  shared_ptr<CoefficientFunction>
  Transform(CoefficientFunction::T_Transform& transformation) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if(transformation.cache.count(thisptr))
      return transformation.cache[thisptr];
    if(transformation.replace.count(thisptr))
      return transformation.replace[thisptr];
    auto newcf = make_shared<TraceCoefficientFunction>(c1->Transform(transformation));
    transformation.cache[thisptr] = newcf;
    return newcf;
  }

  using T_CoefficientFunction<TraceCoefficientFunction>::Evaluate;
  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & mir,
                   BareSliceMatrix<T,ORD> result) const
  {
    int hd = c1->Dimensions()[0];
    STACK_ARRAY(T, hmem, hd*hd*mir.Size());
    FlatMatrix<T,ORD> m1(hd*hd, mir.Size(), &hmem[0]);
    c1->Evaluate (mir, m1);
    
    for (size_t i = 0; i < mir.Size(); i++)
      {
        T sum{0.0};
        for (int j = 0; j < hd; j++)
          sum += m1(j*(hd+1), i);
        result(0, i) = sum;
      }
  }  

  template <typename MIR, typename T, ORDERING ORD>
  void T_Evaluate (const MIR & ir,
                   FlatArray<BareSliceMatrix<T,ORD>> input,                       
                   BareSliceMatrix<T,ORD> values) const
  {
    int hd = c1->Dimensions()[0];
    size_t np = ir.Size();
    
    auto in0 = input[0];
    for (size_t i = 0; i < np; i++)
      {
        T sum{0.0};
        for (size_t j = 0; j < hd; j++)
          sum += in0(j*(hd+1), i);
        values(0,i) = sum;
      }
  }

  shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var,
                                        shared_ptr<CoefficientFunction> dir) const override
  {
    if (this == var) return dir;
    return TraceCF(c1->Diff(var, dir));
  }

  shared_ptr<CoefficientFunction>
  DiffJacobi (const CoefficientFunction * var, T_DJC & cache) const override
  {
    auto thisptr = const_pointer_cast<CoefficientFunction>(this->shared_from_this());
    if (cache.find(thisptr) != cache.end())
      return cache[thisptr];

    if (this == var) return make_shared<ConstantCoefficientFunction>(1);
    if (c1.get() == var) return IdentityCF (c1->Dimensions()[0]);
    auto input = c1->InputCoefficientFunctions();
    if (input.Size() == 0) return ZeroCF(c1->Dimensions());

    shared_ptr<CoefficientFunction> res;
    if (c1->GetDescription()=="binary operation '-'")
      res = TraceCF(input[0])->DiffJacobi (var, cache) - TraceCF(input[1])->DiffJacobi (var, cache);
    else if (dynamic_pointer_cast<MultMatMatCoefficientFunction>(c1) && !c1->IsComplex())
      {
        auto AB = c1->InputCoefficientFunctions();
        res = InnerProduct(AB[0], AB[1]->Transpose() )->DiffJacobi (var, cache);
      }
    else
      res = MakeTensorTraceCoefficientFunction (c1->DiffJacobi (var, cache), 0, 1);
    cache[thisptr] = res;
    return res;
  }
  
};






  
  extern template class T_CoefficientFunction<MultMatMatCoefficientFunction>;
  extern template class T_CoefficientFunction<MultDiagMatCoefficientFunction>;  
  extern template class T_CoefficientFunction<TransposeCoefficientFunction>;
  extern template class T_CoefficientFunction<TraceCoefficientFunction>;    
  extern template class T_CoefficientFunction<IdentityCoefficientFunction>;      
}


#endif
