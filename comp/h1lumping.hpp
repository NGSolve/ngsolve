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
    
    virtual shared_ptr<BaseMatrix> 
    GetMassOperator (shared_ptr<CoefficientFunction> rho,
                     shared_ptr<Region> defon, LocalHeap & lh) const override;


    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    void GetVertexDofNrs (int vnr, Array<int> & dnums) const override;
    void GetEdgeDofNrs (int ednr, Array<int> & dnums) const override;
    void GetFaceDofNrs (int fanr, Array<int> & dnums) const override;
    void GetInnerDofNrs (int elnr, Array<int> & dnums) const override;


    
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules() const;
  };
}    

#endif
