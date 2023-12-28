#ifndef FILE_IRFESPACE
#define FILE_IRFESPACE

/*********************************************************************/
/* File:   irfespace.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   Jan 2021                                                  */
/*********************************************************************/

#include "fespace.hpp"

namespace ngcomp
{

  class IntegrationRuleSpace : public FESpace
  {
    Array<int> firsteldof;
  public:
    IntegrationRuleSpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false);
    void Update() override;

    virtual void UpdateCouplingDofArray() override;    

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    
    virtual void GetDofNrs (ElementId ei, Array<int> & dnums) const override;

    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules() const;
  };


  class IntegrationRuleSpaceSurface : public FESpace
  {
    Array<int> firsteldof;
  public:
    IntegrationRuleSpaceSurface (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false);
    void Update() override;

    virtual void UpdateCouplingDofArray() override;    

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    
    virtual void GetDofNrs (ElementId ei, Array<int> & dnums) const override;

    std::map<ELEMENT_TYPE, IntegrationRule> GetIntegrationRules() const;
  };

}
#endif
