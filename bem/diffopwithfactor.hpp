#ifndef DIFFOPWITHFACTOR_HPP
#define DIFFOPWITHFACTOR_HPP



namespace ngsbem
{
  using namespace ngfem;

  class DifferentialOperatorWithFactor : public DifferentialOperator
  {
    shared_ptr<DifferentialOperator> diffop;
    shared_ptr<CoefficientFunction> factor;

  public:
    DifferentialOperatorWithFactor (shared_ptr<DifferentialOperator> adiffop,
                                    shared_ptr<CoefficientFunction> afactor)
      : DifferentialOperator(afactor->Dimensions()[0], 1, adiffop->VB(), adiffop->DiffOrder()),
        diffop(adiffop), factor(afactor)
    {
      ;
    }

    int DimRef() const override { return diffop->DimRef(); }

    virtual IntRange UsedDofs(const FiniteElement & fel) const override { return diffop->UsedDofs(fel); }

    auto BaseDiffOp() const { return diffop; }

    void CalcMatrix (const FiniteElement & fel,
                     const BaseMappedIntegrationPoint & mip,
                     BareSliceMatrix<double,ColMajor> mat, 
                     LocalHeap & lh) const override
    {
      FlatMatrix<double,ColMajor> hmat(diffop->Dim(), fel.GetNDof(), lh);
      diffop -> CalcMatrix (fel, mip, hmat, lh);

      auto dims = factor->Dimensions();      
      FlatMatrix<double> factorx(dims[0], dims[1], lh);
      factor->Evaluate (mip, factorx.AsVector());

      IntRange used = diffop->UsedDofs(fel);
      mat.Cols(used) = factorx * hmat.Cols(used);
    }
    
    void CalcMatrix (const FiniteElement & fel,
                     const SIMD_BaseMappedIntegrationRule & mir,
                     BareSliceMatrix<SIMD<double>> mat) const override
    {
      // *testout << "CalcMatrix SIMD" << endl;      
      Matrix<SIMD<double>> hmat (fel.GetNDof()*diffop->Dim(), mir.Size());
      // hmat = SIMD<double>(0.0);
      diffop -> CalcMatrix (fel, mir, hmat);

      Matrix<SIMD<double>> fac(factor->Dimension(), mir.Size());
      factor -> Evaluate (mir, fac);

      auto dims = factor -> Dimensions();

      mat.Rows(fel.GetNDof()*dims[1]).Cols(mir.Size()) = SIMD<double>(0.0);
      
      for (size_t i = 0; i < mir.Size(); i++)
        for (size_t j = 0; j < dims[0]; j++)
          for (size_t k = 0; k < dims[1]; k++)
            mat.Col(i).Slice(j,dims[0]) += fac(j*dims[1]+k, i) * hmat.Col(i).Slice(k, dims[1]);
    }
    

    void CalcMatrix (const FiniteElement & fel,
                     const IntegrationPoint & ip,
                     BareSliceMatrix<double,ColMajor> mat,
                     LocalHeap & lh) const override
    {
      diffop -> CalcMatrix(fel, ip, mat, lh);
      /*
      *testout << "calcmatrix mip" << endl
               << mat.Rows(Dim()).Cols(fel.GetNDof()) << endl;
      */
    }

    void CalcTransformationMatrix (const BaseMappedIntegrationPoint & mip,
                                   SliceMatrix<double> trans,
                                   LocalHeap & lh) const override
    {
      HeapReset hr(lh);
      auto dims = factor->Dimensions();
      
      FlatMatrix<double> factorx(dims[0], dims[1], lh);
      factor->Evaluate (mip, factorx.AsVector());

      FlatMatrix<double> basetrans(diffop->Dim(), diffop->DimRef(), lh);
      diffop -> CalcTransformationMatrix(mip, basetrans, lh);

      trans = factorx * basetrans;
      // *testout << "trans = " << trans << endl;
    }


    void Apply (const FiniteElement & fel,
                const SIMD_BaseMappedIntegrationRule & mir,
                BareSliceVector<double> x, 
                BareSliceMatrix<SIMD<double>> flux) const override
    {
      auto dims = factor->Dimensions();

      Matrix<SIMD<double>> tmpflux(dims[1], mir.Size());
      Matrix<SIMD<double>> factorx(dims[0]*dims[1], mir.Size());
      
      diffop -> Apply (fel, mir, x, tmpflux);
      factor -> Evaluate (mir, factorx);
      flux.Rows(0, dims[0]).Cols(0, mir.Size()) = SIMD<double>(0.0);
      for (int i = 0; i < dims[0]; i++)
        for (int j = 0; j < dims[1]; j++)
          flux.Row(i).Range(mir.Size()) += pw_mult(factorx.Row(i*dims[1]+j), tmpflux.Row(j));
    }

    
  };


}

#endif
