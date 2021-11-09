#ifndef NGSOLVE_GLOBALSPACE_HPP
#define NGSOLVE_GLOBALSPACE_HPP

namespace ngcomp
{
  class GlobalSpace : public FESpace
  {
  protected:
    shared_ptr<CoefficientFunction> basis;
    shared_ptr<CoefficientFunction> grad_basis[3];     
    int dim;       // number of basis functions
    int vecdim;    // vector dimension of each basis function

    class FE : public FiniteElement
    {
      ELEMENT_TYPE type;
    public:
      FE(int adim, ELEMENT_TYPE atype)
        : FiniteElement(adim, 5), type(atype) { };
      ELEMENT_TYPE ElementType() const override { return type; }
    };
    
    class VolDiffOp : public DifferentialOperator
    {
      shared_ptr<CoefficientFunction> basis;
      int dim, vecdim;
    public:
      VolDiffOp (shared_ptr<CoefficientFunction> abasis,
                 int adim, int avecdim)
        : DifferentialOperator(avecdim, 1, VOL, 0), basis(abasis),
          dim(adim), vecdim(avecdim) { ; }

      void CalcMatrix (const FiniteElement & bfel,
                       const BaseMappedIntegrationPoint & mip,
                       SliceMatrix<double,ColMajor> mat,
                       LocalHeap & lh) const override;
    };

    
  public:
    GlobalSpace(shared_ptr<MeshAccess> ama, const Flags& flags);

    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.Arg("basis") = "Basis functions.";
      /*
      docu.Arg("periodic") = "Periodic global  space (in 2d in x and y direction).";
      docu.Arg("periodicu") = "Periodic u-dir (local coordinate system) global  space.";
      docu.Arg("periodicv") = "Periodic v-dir (local coordinate system) global  space.";
      */
      return docu;
    }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const;
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const;
  };

} // namespace ngcomp

#endif // NGSOLVE_GLOBALSPACE_HPP
