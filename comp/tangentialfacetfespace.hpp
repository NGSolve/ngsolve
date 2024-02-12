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
    /// Level
    // int level;
    /// Number of Facets
    // int nfacets;
    /// 
    // int ncfacets;
    ///
    // int nel;

    Array<int> first_facet_dof;
    Array<int> first_inner_dof;  // for highest_order_dc
    // int ndof_lo;

    int rel_order;

    Array<IVec<2> > order_facet;
    Array<bool> fine_facet;

    // int ndof;
    // Array<int> ndlevel;
    bool var_order;
    bool print;

    bool highest_order_dc;
    bool hide_highest_order_dc;
    
  public:
    ///
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

    
    // virtual size_t GetNDof() const throw() override { return ndof; }
    // virtual size_t GetNDofLevel ( int i ) const override { return ndlevel[i]; }

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;

    // virtual int GetNDofLowOrder () const
    // { return ndof_lo; }

    virtual FiniteElement & GetFE(ElementId ei, Allocator & lh) const override;
    
    // virtual const FiniteElement & GetFE ( int elnr, LocalHeap & lh ) const;
    // virtual const FiniteElement & GetSFE ( int selnr, LocalHeap & lh ) const;

    virtual void GetFacetDofNrs (int felnr, Array<DofId> & dnums) const;

    virtual int GetNFacetDofs (int felnr) const;

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

    virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & precflags) const override;
    ///
    virtual shared_ptr<Array<int>> CreateDirectSolverClusters (const Flags & precflags) const override;
  
    // some utility functions for convenience
    ///
    // virtual void GetVertexNumbers(int elnr, Array<int>& vnums) const;
    ///
    virtual IVec<2> GetFacetOrder(int fnr) const;

    virtual int GetFirstFacetDof(int fanr) const;

    virtual bool UsesHighestOrderDiscontinuous() const {return highest_order_dc;};

    virtual void GetVertexDofNrs (int elnum, Array<DofId> & dnums) const override;
    virtual void GetEdgeDofNrs (int elnum, Array<DofId> & dnums) const override;
    virtual void GetFaceDofNrs (int felnr, Array<DofId> & dnums) const override;
    virtual void GetInnerDofNrs (int felnr, Array<DofId> & dnums) const override;

  };

}

#endif

