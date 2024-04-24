#ifndef NGSOLVE_GLOBALINTERFACESPACE_HPP
#define NGSOLVE_GLOBALINTERFACESPACE_HPP


#include "fespace.hpp"

namespace ngcomp
{
  class GlobalInterfaceSpace : public FESpace
  {
  protected:
    shared_ptr<CoefficientFunction> mapping;
    int order;
    bool periodic[2]; // periodic u/ periodic v
    bool polar; // using basis for polar coordinates

    template<typename VOLFE>
    class VolDiffOp : public DifferentialOperator
    {
    public:
      VolDiffOp ()
        : DifferentialOperator(1, 1, VOL, 0) { ; }

      void CalcMatrix (const FiniteElement & bfel,
                       const BaseMappedIntegrationPoint & mip,
                       BareSliceMatrix<double,ColMajor> mat,
                       LocalHeap & lh) const override;
    };

    template<typename VOLFE>
    class ParameterGradDiffOp : public DifferentialOperator
    {
    public:
      ParameterGradDiffOp ()
        : DifferentialOperator(VOLFE::GetDim(), 1, VOL, 0) { ; }

      void CalcMatrix (const FiniteElement & bfel,
                       const BaseMappedIntegrationPoint & mip,
                       BareSliceMatrix<double,ColMajor> mat,
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
                       BareSliceMatrix<double,ColMajor> mat,
                       LocalHeap & lh) const override;
    };

  public:
    GlobalInterfaceSpace(shared_ptr<MeshAccess> ama, const Flags& flags);

    static DocInfo GetDocu()
    {
      auto docu = FESpace::GetDocu();
      docu.Arg("mapping") = "Mapping for global interface space.";
      docu.Arg("periodic") = "Periodic global interface space (in 2d in x and y direction).";
      docu.Arg("periodicu") = "Periodic u-dir (local coordinate system) global interface space.";
      docu.Arg("periodicv") = "Periodic v-dir (local coordinate system) global interface space.";
      docu.Arg("polar") = "Polar mapping (r, phi). Uses stable basis for singularity in origin. Automatically sets first argument (r) not periodic and second (phi) periodic. Mapping must be to [0,1]x[-pi, pi)";
      return docu;
    }
  };

  shared_ptr<GlobalInterfaceSpace> CreateGlobalInterfaceSpace
    (shared_ptr<MeshAccess> ma, shared_ptr<CoefficientFunction> mapping,
     optional<Region> definedon, bool periodic, bool periodicu,
     bool periodicv, int order, bool complex, bool polar, bool autoupdate);
} // namespace ngcomp

#endif // NGSOLVE_GLOBALINTERFACESPACE_HPP
