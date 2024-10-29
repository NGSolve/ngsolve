#include "globalspace.hpp"

namespace ngcomp
{

  void GlobalSpace::VolDiffOp::CalcMatrix
    (const FiniteElement & bfel,
     const BaseMappedIntegrationPoint & mip,
     BareSliceMatrix<double,ColMajor> mat,
     LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<double> basisvec(basis->Dimension(), lh);
    basis -> Evaluate (mip, basisvec);

    for (int i = 0, ii = 0; i < dim; i++)
      for (int j = 0; j < vecdim; j++)
        mat(j,i) = basisvec(ii++);
  }


  void GlobalSpace::VolDiffOp::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationRule & mir,
         BareSliceVector<double> x, 
         BareSliceMatrix<double> flux,
         LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<double> basisvecs(mir.Size(), basis->Dimension(), lh);
    basis -> Evaluate (mir, basisvecs);
    for (int i = 0; i < mir.Size(); i++)
      flux.Row(i) = Trans(basisvecs.Row(i).AsMatrix(dim, vecdim)) * x;
  }

  void GlobalSpace::VolDiffOp::  
  ApplyTrans (const FiniteElement & fel,
              const BaseMappedIntegrationRule & mir,
              FlatMatrix<double> flux,
              BareSliceVector<double> bx, 
              LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<double> basisvecs(mir.Size(), basis->Dimension(), lh);
    basis -> Evaluate (mir, basisvecs);
    auto x = bx.Range(dim);
    x = 0.0;
    for (int i = 0; i < mir.Size(); i++)
      x += basisvecs.Row(i).AsMatrix(dim, vecdim) * flux.Row(i);
  }

  /*
  void GlobalSpace::VolDiffOp::  
  AddTrans (const FiniteElement & bfel,
            const SIMD_BaseMappedIntegrationRule & bmir,
            BareSliceMatrix<SIMD<double>> flux,
            BareSliceVector<double> x) const
  {
    Matrix<SIMD<double>> basisvecs(basis->Dimension(), bmir.Size());
    basis -> Evaluate (bmir, basisvecs);
    for (size_t i = 0; i < dim; i++)
      for (size_t j = 0; j < vecdim; j++)
        x(i) += InnerProduct (basisvecs.Row(i*vecdim+j), flux.Row(j));
  }
  */
  
  void GlobalSpace::VolDiffOp::CalcMatrix
    (const FiniteElement & bfel,
     const BaseMappedIntegrationPoint & mip,
     BareSliceMatrix<Complex,ColMajor> mat,
     LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatVector<Complex> basisvec(basis->Dimension(), lh);
    basis -> Evaluate (mip, basisvec);

    for (int i = 0, ii = 0; i < dim; i++)
      for (int j = 0; j < vecdim; j++)
        mat(j,i) = basisvec(ii++);
  }

  void GlobalSpace::VolDiffOp::
  Apply (const FiniteElement & fel,
         const BaseMappedIntegrationPoint & mip,
         BareSliceVector<Complex> x, 
         FlatVector<Complex> flux,
         LocalHeap & lh) const
  {
    HeapReset hr(lh);
    FlatMatrix<Complex,ColMajor> mat(Dim(), fel.GetNDof(), lh);
    CalcMatrix (fel, mip, mat, lh);
    flux = mat * x;
  }


  

  GlobalSpace::GlobalSpace(shared_ptr<MeshAccess> ama,
                           const Flags& flags)
    : FESpace(ama, flags)
  {
    if (!flags.NumFlagDefined("order"))
      order = 5;
    
    basis = std::any_cast<shared_ptr<CoefficientFunction>>(flags.GetAnyFlag("basis"));

    dim = CalcDim(basis);
    vecdim = CalcVecDim(basis);
    complex_basis = basis->IsComplex();
    if (complex_basis) iscomplex = true;

    SetNDof(dim);
    evaluator[VOL] = make_shared<VolDiffOp>(basis);
    evaluator[BND] = make_shared<VolDiffOp>(basis,BND);
    evaluator[BBND] = make_shared<VolDiffOp>(basis,BBND);    
    evaluator[BBBND] = make_shared<VolDiffOp>(basis,BBBND);    
  }

  void GlobalSpace :: AddOperator (string name, VorB vb, shared_ptr<CoefficientFunction> dbasis)
  {
    additional_evaluators.Set (name, make_shared<VolDiffOp>(dbasis, vb)); 
  }
  
  FiniteElement & GlobalSpace :: GetFE (ElementId ei, Allocator & lh) const
  {
    return *new (lh) FE(dim, order, ma->GetElType(ei), complex_basis);
  }
  
  void GlobalSpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize0();
    for (int i = 0; i < dim; i++)
      dnums.Append(i);
  }
  
}
