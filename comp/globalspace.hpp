#ifndef NGSOLVE_GLOBALSPACE_HPP
#define NGSOLVE_GLOBALSPACE_HPP

#include "fespace.hpp"

namespace ngcomp
{
  class GlobalSpace : public FESpace
  {
  protected:
    shared_ptr<CoefficientFunction> basis;
    shared_ptr<CoefficientFunction> grad_basis[3];     
    int dim;       // number of basis functions
    int vecdim;    // vector dimension of each basis function
    bool complex_basis;
    
    static int CalcDim(shared_ptr<CoefficientFunction> basis)
    {
      auto dims = basis->Dimensions();
      return (dims.Size() == 0) ? 1 : dims[0]; 
    }

    static int CalcVecDim(shared_ptr<CoefficientFunction> basis)
    {
      auto dims = basis->Dimensions();
      return (dims.Size() <= 1) ? 1 : dims[1];
    }
    
    class FE : public FiniteElement
    {
      ELEMENT_TYPE type;
      bool complex;
    public:
      FE(int adim, ELEMENT_TYPE atype, bool acomplex)
        : FiniteElement(adim, 5), type(atype), complex(acomplex) { };
      ELEMENT_TYPE ElementType() const override { return type; }
      bool ComplexShapes() const override { return complex; }
      
    };
    
    class VolDiffOp : public DifferentialOperator
    {
      shared_ptr<CoefficientFunction> basis;
      int dim, vecdim;
    public:
      VolDiffOp (shared_ptr<CoefficientFunction> abasis, VorB vb = VOL)
        : DifferentialOperator(CalcVecDim(abasis), 1, vb, 0), basis(abasis),
          dim(CalcDim(abasis)), vecdim(CalcVecDim(abasis)) { ; }

      virtual bool SupportsVB (VorB checkvb) const override { return true; }

      void CalcMatrix (const FiniteElement & bfel,
                       const BaseMappedIntegrationPoint & mip,
                       BareSliceMatrix<double,ColMajor> mat,
                       LocalHeap & lh) const override;

      void CalcMatrix (const FiniteElement & bfel,
                       const BaseMappedIntegrationPoint & mip,
                       BareSliceMatrix<Complex,ColMajor> mat,
                       LocalHeap & lh) const override;

      void Apply (const FiniteElement & fel,
                  const BaseMappedIntegrationPoint & mip,
                  BareSliceVector<Complex> x, 
                  FlatVector<Complex> flux,
                  LocalHeap & lh) const override;
      
    };

    
  public:
    GlobalSpace(shared_ptr<MeshAccess> ama, const Flags& flags);

    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.Arg("basis") = "Basis functions.";
      return docu;
    }

    void AddOperator (string name, VorB vb, shared_ptr<CoefficientFunction> dbasis);
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const;
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const;
  };

} // namespace ngcomp

#endif // NGSOLVE_GLOBALSPACE_HPP
