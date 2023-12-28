/**
   High Order Finite Element Space for SurfFacets
*/
/* facetsurffespace.*pp*/
#ifndef FACETSURF_FESPACE_HPP
#define FACETSURF_FESPACE_HPP


#include "fespace.hpp"

namespace ngcomp
{
  class NGS_DLL_HEADER FacetSurfaceFESpace : public FESpace 
  {
  protected:  
    // Level
    // int level;
    // Number of Facets
    int nfa;
    // Number of Elements
    int nel;
  
    Array<int> first_edge_dof;
  
    // relative order to mesh-order
    int rel_order; 
  
    // int ndof;
    // Array<int> ndlevel;
    bool var_order; 
    bool nowirebasket;
    
  public:
    ///
    FacetSurfaceFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~FacetSurfaceFESpace ();
    ///
    virtual string GetClassName () const override
    {
      return "FacetSurfFESpace";
    }
  
    ///
    void Update() override;
    ///
    virtual void UpdateCouplingDofArray() override;

    ///
    // virtual size_t GetNDof () const throw() override;
    ///
    // virtual size_t GetNDofLevel (int level) const override;

    template <ELEMENT_TYPE ET>
    FiniteElement & T_GetFE (int elnr, Allocator & alloc) const;
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;

    using FESpace::GetDofNrs;
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    ///

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;

    IntRange GetEdgeDofs (int nr) const
    { 
      return IntRange (first_edge_dof[nr], first_edge_dof[nr+1]); 
    }
  
    virtual void GetFacetDofNrs (int nr, Array<DofId> & dnums) const
    {
      dnums.SetSize(0);
    }

    ///
    virtual int GetNEdgeDofs (int felnr) const 
    { return (first_edge_dof[felnr+1]-first_edge_dof[felnr] + 1); }
    ///

    virtual int GetFirstEdgeDof(int ednr) const 
    {
      return first_edge_dof[ednr];
    }

    virtual void GetVertexDofNrs ( int nr, Array<DofId> & dnums ) const override
    {
      dnums.SetSize0();
    }

    virtual void GetEdgeDofNrs ( int nr, Array<DofId> & dnums ) const override
    {
      dnums = GetEdgeDofs(nr);
    }

    virtual void GetFaceDofNrs (int nr, Array<DofId> & dnums) const override
    {
      dnums.SetSize(0);
      return;
    }
  
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override
    {
      dnums.SetSize(0);
    }
    
  };

}



#endif
