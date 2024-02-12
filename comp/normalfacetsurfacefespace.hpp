#ifndef NORMAL_FACET_SURFACE_FESPACE_HPP
#define NORMAL_FACET_SURFACE_FESPACE_HPP

/*********************************************************************/
/* File:   normalfacetsurfacefespace.hpp                             */
/* Author: Michael Neunteufel                                        */
/* Date:   2020                                                      */
/*********************************************************************/

#include "fespace.hpp"

namespace ngcomp
{

  class NGS_DLL_HEADER NormalFacetSurfaceFESpace : public FESpace
  {
  protected:
    Array<int> first_facet_dof;

    int rel_order;

    Array<IVec<2> > order_facet;
    Array<bool> fine_facet;

    int ndof;
    Array<int> ndlevel;
    bool var_order;
    bool print;

    
  public:
    NormalFacetSurfaceFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
			bool parseflags = false );

    virtual ~NormalFacetSurfaceFESpace () { ; }
    static DocInfo GetDocu ();

    virtual string GetClassName () const override
    {
      return "NormalFacetSurfaceFESpace";
    }

    void Update() override;
    virtual void UpdateCouplingDofArray() override;

    virtual void SetOrder (NodeId ni, int order) override;
    virtual int GetOrder (NodeId ni) const override;
    
    virtual size_t GetNDof() const throw() override { return ndof; }

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;


    virtual FiniteElement & GetFE(ElementId ei, Allocator & lh) const override;
    
    virtual int GetNFacetDofs (int felnr) const;

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

    auto GetFacetDofs (size_t nr) const
    {
      return Range (first_facet_dof[nr], first_facet_dof[nr+1]);
    }
    
  
    // some utility functions for convenience
    ///
    // virtual void GetVertexNumbers(int elnr, Array<int>& vnums) const;
    ///
    virtual IVec<2> GetFacetOrder(int fnr) const;

    virtual int GetFirstFacetDof(int fanr) const;

    virtual void GetVertexDofNrs (int elnum, Array<DofId> & dnums) const override;
    virtual void GetEdgeDofNrs (int elnum, Array<DofId> & dnums) const override;
    virtual void GetFaceDofNrs (int felnr, Array<DofId> & dnums) const override;
    virtual void GetInnerDofNrs (int felnr, Array<DofId> & dnums) const override;

  };

}

#endif

