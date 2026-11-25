#ifndef TANGENTIAL_FACET_FESPACE_HPP
#define TANGENTIAL_FACET_FESPACE_HPP

/*********************************************************************/
/* File:   tangentialfacetfespace.hpp                                */
/* Author: A. Sinwel, Joachim Schoeberl                              */
/* Date:   2008                                                      */
/*********************************************************************/


#include "fespace.hpp"

namespace ngcomp
{

  class NGS_DLL_HEADER TangentialFacetFESpace : public FESpace
  {
  protected:
    Array<int> first_facet_dof;
    Array<int> first_inner_dof;  // for highest_order_dc

    int rel_order;

    Array<IVec<2> > order_facet;
    Array<bool> fine_facet;

    bool var_order;
    bool print;

    bool highest_order_dc;
    bool hide_highest_order_dc;
    bool all_dofs_together;
    
  public:
    TangentialFacetFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
                            bool parseflags = false );

    virtual ~TangentialFacetFESpace () { ; }
    static DocInfo GetDocu ();

    virtual string GetClassName () const override
    {
      return "TangentialFacetFESpace";
    }

    void Update() override;
    virtual void UpdateCouplingDofArray() override;

    virtual void SetOrder (NodeId ni, int order) override;
    virtual int GetOrder (NodeId ni) const override;
    void SetOrder (const Region & region, int order) override
    { SetRegionOrder (region, max(order, 0)); }
    void SetOrder (const Region & region, IVec<2> order, shared_ptr<CoefficientFunction> direction) override
    { SetRegionOrder (region, max(order[0], 0), max(order[1], 0), direction); }
    void SetOrder (ElementId ei, IVec<3> order) override
    { SetCellOrder (ei, IVec<3> (max(order[0],0), max(order[1],0), max(order[2],0))); }
    IVec<3> GetOrder (ElementId ei) const override
    { bool elementwise; return CellBaseOrder (ei, ma->GetElIndex(ei), elementwise); }

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;


    virtual FiniteElement & GetFE(ElementId ei, Allocator & lh) const override;
    
    virtual void GetFacetDofNrs (int felnr, Array<DofId> & dnums) const;

    virtual int GetNFacetDofs (int felnr) const;

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    ///
    virtual shared_ptr<Array<int>> CreateDirectSolverClusters (const Flags & precflags) const override;
  
    // some utility functions for convenience
    ///
    // virtual void GetVertexNumbers(int elnr, Array<int>& vnums) const;
    ///
    virtual IVec<2> GetFacetOrder(int fnr) const;

    virtual int GetFirstFacetDof(int fanr) const;
    const auto & GetFirstFacetDof() const { return first_facet_dof; }

    IntRange GetFacetDofs (int nr) const
    { 
      return IntRange (first_facet_dof[nr], first_facet_dof[nr+1]); 
    }

    virtual bool UsesHighestOrderDiscontinuous() const {return highest_order_dc;};

    virtual void GetVertexDofNrs (int elnum, Array<DofId> & dnums) const override;
    virtual void GetEdgeDofNrs (int elnum, Array<DofId> & dnums) const override;
    virtual void GetFaceDofNrs (int felnr, Array<DofId> & dnums) const override;
    virtual void GetInnerDofNrs (int felnr, Array<DofId> & dnums) const override;

  };

}

#endif

