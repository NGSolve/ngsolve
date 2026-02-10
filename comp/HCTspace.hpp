#ifndef HCTSPACE_HPP
#define HCTSPACE_HPP

#include "fespace.hpp"

namespace ngcomp
{
  class HCT_FESpace : public FESpace
  { 
  protected:
    int order;
    Array<int> first_edge_dof;
    Array<int> first_element_dof;

  public:
    HCT_FESpace(shared_ptr<MeshAccess> ama, const Flags &flags, bool checkflags = false);
    ~HCT_FESpace() override = default;
    static DocInfo GetDocu();
    void Update() override;
    void FinalizeUpdate() override;
    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules() const;

    void GetDofNrs(ElementId ei, Array<int> &dnums) const override;
    FiniteElement &GetFE(ElementId ei, Allocator &alloc) const override;
  };

}

#endif // HCTSPACE_HPP
