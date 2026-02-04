#ifndef JKMSPACE_HPP
#define JKMSPACE_HPP

#include "fespace.hpp"

namespace ngcomp
{
  class JKM_FESpace : public FESpace
  { 
  protected:
    int order;
    Array<int> first_edge_dof;
    Array<int> first_element_dof;

  public:
    JKM_FESpace(shared_ptr<MeshAccess> ama, const Flags &flags, bool checkflags = false);
    ~JKM_FESpace() override = default;
    static DocInfo GetDocu();
    void Update() override;
    void FinalizeUpdate() override;
    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules(int bonus_intorder=2) const;

    void GetDofNrs(ElementId ei, Array<int> &dnums) const override;
    FiniteElement &GetFE(ElementId ei, Allocator &alloc) const override;
  };

}

#endif // JKMSPACE_HPP
