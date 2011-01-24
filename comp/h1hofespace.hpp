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

    bool print; 

    Array<int> first_edge_dof;
    Array<int> first_face_dof;
    Array<int> first_element_dof;


    /// relative order to mesh-order
    int rel_order; 
    bool var_order; 
    bool fixed_order;
    bool wb_loedge;
    Array<int> order_edge;
    Array<INT<2> > order_face;
    Array<INT<3> > order_inner;
    Array<bool> fine_edge; 
    Array<bool> fine_face; 

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


    Array<INT<2> > defined_on_one_side_of_bounding_curve;

  public:

    H1HighOrderFESpace (const MeshAccess & ama, const Flags & flags, bool checkflags=false);
    ///
    virtual ~H1HighOrderFESpace ();

    // static FESpace * Create (const MeshAccess & ma, const Flags & flags);

    virtual string GetClassName () const
    {
      return "H1HighOrderFESpace";
    }

    ///
    virtual void Update(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    ///
    virtual int GetNDof () const;
    ///
    virtual int GetNDofLevel (int alevel) const;
    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    ///
    virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const;
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;

    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const;
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;

    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
  
    virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const; 
    // virtual void CreateSmoothingBlocks2 (SmoothingBlocksCreator & sbc, const Flags & precflags) const; 
    ///
    virtual Array<int> * CreateDirectSolverClusters (const Flags & flags) const;

    ///
    int GetFirstFaceDof(int i) const {return(first_face_dof[i]);} ;  
    ///
    int GetFirstEdgeDof(int i) const {return(first_edge_dof[i]);} ; 
    ///
    int GetFirstElementDof(int i) const {return(first_element_dof[i]);} ; 

    void UpdateDofTables ();
    ///
    void UpdateCouplingDofArray();    
    
    void SetEdgeOrder (int enr, int eo) { order_edge[enr] = eo; }
    void SetFaceOrder (int fnr, int fo) { order_face[fnr] = INT<2> (fo, fo); }
    void SetFaceOrder (int fnr, int ox, int oy) { order_face[fnr] = INT<2> (ox, oy); }
    void SetElementOrder (int elnr, int elo) 
    { order_inner[elnr] = INT<3> (elo, elo, elo); }
    void SetElementOrder (int elnr, int ox, int oy, int oz) 
    { order_inner[elnr] = INT<3> (ox, oy, oz); }

    // int GetAugmented() const { return augmented; }

    /// get relative (to mesh) order of finite elements
    virtual int GetRelOrder() const { return rel_order; }
    virtual bool VarOrder() const { return var_order; }

    void RestrictToOneSideOfBoundingCurve(int index1, int index2);
    void DeleteOneSideRestrictions(void);

#ifdef PARALLEL
    virtual void UpdateParallelDofs_loproc();
#endif
    
  protected:
    IntRange GetEdgeDofs (int nr) const
    {
      return IntRange (first_edge_dof[nr], 
                       first_edge_dof[nr+1]);
    }

    IntRange GetFaceDofs (int nr) const
    {
      return IntRange (first_face_dof[nr], 
                       first_face_dof[nr+1]);
    }

    IntRange GetElementDofs (int nr) const
    {
      return IntRange (first_element_dof[nr], 
                       first_element_dof[nr+1]);
    }

  };

}

#endif

