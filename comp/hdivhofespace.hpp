#ifndef FILE_HDIVHOFESPACE
#define FILE_HDIVHOFESPACE

/*********************************************************************/
/* File:   hdivhofespace.hpp                                         */
/* Author: START                                                     */
/* Date:   Feb. 2003, update (SZ) Jan. 2007                          */
/*********************************************************************/

#include "fespace.hpp"

namespace ngcomp
{

  /**
     HDiv High Order Finite Element Space
  */

  class NGS_DLL_HEADER HDivHighOrderFESpace : public FESpace
  {
  protected:
    // Level
    int level;
    // Number Dofs 
    // int ndof;
    // order of curl-fields 
    int curl_order; 

    Array<DofId> first_facet_dof;
    Array<DofId> first_inner_dof;

    /// relative order to mesh-order
    int rel_order; 
    // curl-order relative to mesh order 
    int rel_curl_order; 
    // space is of variable order 
    bool var_order;
    // space is continuous/discontinuous 
    bool discont; 
    
    Array<IVec<3> > order_inner;
    Array<IVec<3> > order_inner_curl;
    Array<IVec<2> > order_facet; 
    Array<bool> fine_facet; 
    Array<bool> boundary_facet; 
 
    // Array<int> ndlevel;
    int uniform_order_inner; 
    int uniform_order_facet; 

    // high order divergence free
    bool ho_div_free; 
    ///
    bool fixed_order;

    bool RT = false; 

    // highest-order facet dofs discontinuous ?
    bool highest_order_dc;
    // hide all used dofs
    bool hide_all_dofs;
    // discontinuous facet pairs:
    Array<IVec<2>> dc_pairs;
  public:
    HDivHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
                          bool parseflags=false);
    ///
    virtual ~HDivHighOrderFESpace ();
    static DocInfo GetDocu ();

    virtual void UpdateDofTables() override;  

    virtual void UpdateCouplingDofArray() override;

    virtual void SetOrder (NodeId ni, int order) override;
    virtual int GetOrder (NodeId ni) const override;
    using FESpace::GetOrder;

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override
    {
      static VorB nodes[] = { VOL, BND };
      return FlatArray<VorB> (2, &nodes[0]); 
    }
    
    virtual string GetClassName () const override
    {
      return "HDivHighOrderFESpace";
    }

    /// averaging over high-order dc dofs
    void Average (BaseVector & vec) const;
    
    ///
    void Update() override;
    ///
    // virtual size_t GetNDof () const throw() override;
    ///
    // virtual size_t GetNDofLevel (int level) const override;
    ///
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;
    ///
    // virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const override;
    ///
    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetFE (int elnr, bool onlyhdiv, Allocator & lh) const;
    ///
    virtual const FiniteElement & GetHODivFE (int elnr, LocalHeap & lh) const;
    ///
    // virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const override; // 2D: array =0.;
    ///
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    ///
    virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & precflags) const override;
    /// 
    virtual shared_ptr<Array<int>> CreateDirectSolverClusters (const Flags & precflags) const override;
    /// 
    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
    /// 
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
    /// 
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
    /// 
    virtual void GetFacetDofNrs(int fanr, Array<DofId> & dnums) const
    { 
      if (ma->GetDimension() == 2) GetEdgeDofNrs(fanr,dnums); 
      else if (ma->GetDimension() == 3) GetFaceDofNrs(fanr,dnums); 
    } 
    ///
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override; 
    /// 
    void GetFacetOrder (Array<IVec<2> > & of, Array<bool> & ff) const 
    {of = order_facet; ff = fine_facet;};

    /// 
    int GetNElemDofs(int elnr) const 
    {
      if(discont) return(first_inner_dof[elnr+1] - first_inner_dof[elnr]); 
      else 
	{ 
	  Array<DofId> dnums; 
	  this->GetDofNrs(elnr,dnums); 
	  return(dnums.Size()); 
	} 
    }

    int GetFirstInnerDof(int elnr) const { return(first_inner_dof[elnr]);}; 
    // virtual int LowOrderDof() const { if(discont) return(0); else return(1);} 


    const Array<IVec<2>> & GetDCPairs () const { return dc_pairs; }

    virtual bool VarOrder() const override { return var_order; } 
    virtual int GetRelOrder() const override { return rel_order; } 


    auto GetFacetDofs (size_t nr) const
    {
      return Range (first_facet_dof[nr], first_facet_dof[nr+1]);
    }

    auto GetElementDofs (size_t nr) const
    {
      return Range (first_inner_dof[nr], first_inner_dof[nr+1]);
    }

  };

}

#endif





