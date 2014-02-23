#ifndef FILE_H1HOFESPACE
#define FILE_H1HOFESPACE

/*********************************************************************/
/* File:   h1hofespace.hpp                                           */
/* Author: Start                                                     */
/* Date:   10. Feb. 2003                                             */
/*********************************************************************/

namespace ngcomp
{

  /**
     High Order Finite Element Space
  */


  class NGS_DLL_HEADER H1HighOrderFESpace : public FESpace
  {
  protected:
    int level;

    Array<int> first_edge_dof;
    Array<int> first_face_dof;
    Array<int> first_element_dof;

    // typedef short TORDER;
    typedef unsigned char TORDER;
    
    /// relative order to mesh-order
    int rel_order; 
    bool var_order; 
    bool fixed_order;
    bool wb_loedge;
    Array<TORDER> order_edge;
    Array<INT<2,TORDER> > order_face;
    Array<INT<3,TORDER> > order_inner;
    Array<bool> used_vertex; 
    Array<bool> used_edge; 
    Array<bool> used_face; 

    int ndof;
    int uniform_order_inner;
    int uniform_order_face;
    int uniform_order_edge;
    int uniform_order_quad;
    int uniform_order_trig;
    Array<INT<3> > dom_order_min; 
    Array<INT<3> > dom_order_max;
    int smoother; 
  
    Array<int> ndlevel;

    bool level_adapted_order; 
    bool nodalp2;
  public:

    H1HighOrderFESpace (const MeshAccess & ama, const Flags & flags, bool checkflags=false);
    ///
    virtual ~H1HighOrderFESpace ();

    virtual string GetClassName () const
    {
      return "H1HighOrderFESpace";
    }

    ///
    virtual void Update(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual void DoArchive (Archive & archive);

    ///
    virtual int GetNDof () const;
    ///
    virtual int GetNDofLevel (int alevel) const;
    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    ///
    template <ELEMENT_TYPE ET>
    const FiniteElement & T_GetFE (int elnr, LocalHeap & lh) const;
    ///
    virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const;
    ///
    template <ELEMENT_TYPE ET>
    const FiniteElement & T_GetSFE (int elnr, LocalHeap & lh) const;
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;

    virtual void GetDofRanges (ElementId ei, Array<IntRange> & dranges) const;

    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const override;
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
  
    virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const; 
    // virtual void CreateSmoothingBlocks2 (SmoothingBlocksCreator & sbc, const Flags & precflags) const; 
    ///
    virtual Array<int> * CreateDirectSolverClusters (const Flags & flags) const;

    void UpdateDofTables ();
    ///
    virtual void UpdateCouplingDofArray();    
    
    void SetEdgeOrder (int enr, int eo) { order_edge[enr] = eo; }
    void SetFaceOrder (int fnr, INT<2> fo) { order_face[fnr] = fo; }
    void SetElementOrder (int elnr, INT<3> elo) { order_inner[elnr] = elo; }

    /// get relative (to mesh) order of finite elements
    virtual int GetRelOrder() const { return rel_order; }
    virtual bool VarOrder() const { return var_order; }

    IntRange GetEdgeDofs (int nr) const
    {
      return IntRange (first_edge_dof[nr], first_edge_dof[nr+1]);
    }

    IntRange GetFaceDofs (int nr) const
    {
      return IntRange (first_face_dof[nr], first_face_dof[nr+1]);
    }

    IntRange GetElementDofs (int nr) const
    {
      return IntRange (first_element_dof[nr], first_element_dof[nr+1]);
    }

  };

}

#endif

