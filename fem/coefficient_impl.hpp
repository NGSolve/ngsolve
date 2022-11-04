#ifndef FILE_COEFFICIENT_IMPL
#define FILE_COEFFICIENT_IMPL


namespace ngfem
{

  class IdentityCoefficientFunction : public T_CoefficientFunction<IdentityCoefficientFunction>
  {
    using BASE = T_CoefficientFunction<IdentityCoefficientFunction>;
  public:
    IdentityCoefficientFunction (int dim)
      : T_CoefficientFunction<IdentityCoefficientFunction>(1, false)
    {
      SetDimensions (ngstd::INT<2> (dim, dim) );
    }

    // For archive
    IdentityCoefficientFunction() = default;
    virtual ~IdentityCoefficientFunction ();
    
    void DoArchive(Archive& ar) override
    {
      BASE::DoArchive(ar);
    }

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
                                 FlatVector<AutoDiffDiff<1,bool>> values) const override
    {
      int hd = Dimensions()[0];
      values = AutoDiffDiff<1,bool>(false);

      for (int i = 0; i < hd; i++)
        values(i*(hd+1)) = AutoDiffDiff<1,bool>(true);
    }
  
    virtual void NonZeroPattern (const class ProxyUserData & ud,
                                 FlatArray<FlatVector<AutoDiffDiff<1,bool>>> input,
                                 FlatVector<AutoDiffDiff<1,bool>> values) const override
    {
      int hd = Dimensions()[0];
      values = AutoDiffDiff<1,bool>(false);

      for (int i = 0; i < hd; i++)
        values(i*(hd+1)) = AutoDiffDiff<1,bool>(true);
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


  struct GenericIdentity
  {
    template <typename T> T operator() (T x) const { return x; }
    static string Name() { return  " "; }
    void DoArchive(Archive& ar) {}
  };

}

#endif
