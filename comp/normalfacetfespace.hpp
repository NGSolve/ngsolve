#ifndef NORMAL_FACET_FESPACE_HPP
#define NORMAL_FACET_FESPACE_HPP

/*********************************************************************/
/* File:   normalfacetfespace.hpp                                    */
/* Author: Joachim Schoeberl                                         */
/* Date:   2019                                                      */
/*********************************************************************/


#include "fespace.hpp"

namespace ngcomp
{

  class NGS_DLL_HEADER NormalFacetFESpace : public FESpace
  {
  protected:
    int level;
    Array<int> first_facet_dof;
    Array<int> first_inner_dof;  // for highest_order_dc

    int rel_order;

    Array<IVec<2> > order_facet;
    Array<bool> fine_facet;

    int ndof;
    Array<int> ndlevel;
    bool var_order;
    bool print;

    bool highest_order_dc;
    bool hide_highest_order_dc;
    
  public:
    NormalFacetFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
			bool parseflags = false );

    virtual ~NormalFacetFESpace () { ; }
    static DocInfo GetDocu ();

    virtual string GetClassName () const override
    {
      return "NormalFacetFESpace";
    }

    void Update() override;
    virtual void UpdateCouplingDofArray() override;

    virtual void SetOrder (NodeId ni, int order) override;
    virtual int GetOrder (NodeId ni) const override;
    
    virtual size_t GetNDof() const throw() override { return ndof; }

    virtual size_t GetNDofLevel ( int i ) const override { return ndlevel[i]; }

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;


    virtual FiniteElement & GetFE(ElementId ei, Allocator & lh) const override;
    
    virtual void GetFacetDofNrs (int felnr, Array<DofId> & dnums) const;

    virtual int GetNFacetDofs (int felnr) const;

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

    auto GetFacetDofs (size_t nr) const
    {
      return Range (first_facet_dof[nr], first_facet_dof[nr+1]);
    }
    

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

