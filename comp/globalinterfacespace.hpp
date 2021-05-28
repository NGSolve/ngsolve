#ifndef NGSOLVE_GLOBALINTERFACESPACE_HPP
#define NGSOLVE_GLOBALINTERFACESPACE_HPP

namespace ngcomp
{
  class GlobalInterfaceSpace : public FESpace
  {
  protected:
    shared_ptr<CoefficientFunction> mapping;
    int order;
    bool periodic[2]; // periodic x/ periodic y

    template<typename VOLFE>
    class VolDiffOp : public DifferentialOperator
    {
    public:
      VolDiffOp ()
        : DifferentialOperator(1, 1, VOL, 0) { ; }

      void CalcMatrix (const FiniteElement & bfel,
                       const BaseMappedIntegrationPoint & mip,
                       SliceMatrix<double,ColMajor> mat,
                       LocalHeap & lh) const override;
    };

    template<typename INTERFACEFE>
    class InterfaceDiffOp : public DifferentialOperator
    {
    public:
      InterfaceDiffOp ()
        : DifferentialOperator(1, 1, BND, 0) { ; }

      void CalcMatrix (const FiniteElement & fel,
                       const BaseMappedIntegrationPoint & mip,
                       SliceMatrix<double,ColMajor> mat,
                       LocalHeap & lh) const override;
    };

  public:
    GlobalInterfaceSpace(shared_ptr<MeshAccess> ama, const Flags& flags);
  };

  shared_ptr<GlobalInterfaceSpace> CreateGlobalInterfaceSpace
    (shared_ptr<MeshAccess> ma, shared_ptr<CoefficientFunction> mapping,
     optional<Region> definedon, bool periodic, bool periodicu,
     bool periodicv, int order);
} // namespace ngcomp

#endif // NGSOLVE_GLOBALINTERFACESPACE_HPP
