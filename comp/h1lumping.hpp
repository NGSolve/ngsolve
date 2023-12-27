#ifndef FILE_H1LUMPING
#define FILE_H1LUMPING


#include "fespace.hpp"

namespace ngcomp
{

  class H1LumpingFESpace : public FESpace
  {
    size_t nvert, nedge, nface;
  public:
    H1LumpingFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);
    
    string GetClassName () const override { return "h1lumping"; }

    static DocInfo GetDocu();

    void Update() override;
    
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules() const;
  };
}    

#endif
