#ifndef FILE_L2HOFESPACE
#define FILE_L2HOFESPACE

/*********************************************************************/
/* File:   l2hofespace.hpp                                           */
/* Author: Start                                                     */
/* Date:   23.Feb. 2003                                              */
/*********************************************************************/

namespace ngcomp
{
  
  /**
     High Order Finite Element Space for L2 (element by element)
  */

  class NGS_DLL_HEADER L2HighOrderFESpace : public FESpace
  {
  protected:
    // Number of Elements
    int nel;
    // Degrees of Freedom 
    int ndof;
    // Levels 
    Array<int> ndlevel;
    // if order is relative to mesh order 
    bool var_order;
    // variable order is set to mesh_order + rel_order 
    int rel_order;
    // order of elements 
    Array<INT<3> > order_inner;
    // table of first element dofnumber 
    Array<int> first_element_dof;
    bool all_dofs_together;
  public:

    L2HighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~L2HighOrderFESpace ();
    // Create if order=0 ElementFESpace Constructor, else L2HOFE 
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags);
    // Creates also for order=0 a L2HighOrderFESpace 
    static shared_ptr<FESpace> CreateHO (shared_ptr<MeshAccess> ma, const Flags & flags) 
    {
      return make_shared<L2HighOrderFESpace> (ma, flags);
    }  
  
    virtual string GetClassName () const
    {
      return "L2HighOrderFESpace";
    }

    bool AllDofsTogether(){return all_dofs_together;};
    ///
    virtual void Update(LocalHeap & lh);
    /// 
    virtual void UpdateDofTables();
    ///
    virtual void UpdateCouplingDofArray();    
    ///
    virtual int GetNDof () const throw();
    ///
    virtual int GetNDofLevel (int level) const;
    ///

    using FESpace::GetFE;
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    template <ELEMENT_TYPE ET>
    const FiniteElement & T_GetFE (int elnr, LocalHeap & lh) const;

    ///
    virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const;
    ///
    ///
    virtual const FiniteElement & GetFacetFE (int fnr, LocalHeap & lh) const;

    virtual void GetDofRanges (ElementId ei, Array<IntRange> & dranges) const;

    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
    ///
    virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;
    /// 
 
    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const;
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;


    IntRange GetElementDofs (int nr) const
    {
      return IntRange (first_element_dof[nr], 
                       first_element_dof[nr+1]);
    }

    /*
    int GetFirstInnerDof(int elnr) const  
    { return (first_element_dof[elnr]); }

  
    virtual int GetNElemDofs(int elnr) const 
    { return(first_element_dof[elnr+1] - first_element_dof[elnr]+1); }
    */
  };












  class NGS_DLL_HEADER L2SurfaceHighOrderFESpace : public FESpace
  {
  protected:
  
    // Level
    int level;

    // Number of Elements
    int nel;

    Array<int> first_element_dof;
    int ndof;


    // if order is relative to mesh order 
    bool var_order;
    // variable order is set to mesh_order + rel_order 
    int rel_order;
    // order of elements 
    Array<INT<3> > order_cell;

  public:

    L2SurfaceHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~L2SurfaceHighOrderFESpace ();

    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags);

    virtual string GetClassName () const
    {
      return "L2SurfaceHighOrderFESpace";
    }

    ///
    virtual void Update(LocalHeap & lh);
    /// 
    //virtual void UpdateDofTables();
    ///
    virtual int GetNDof () const throw();
    ///
    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    ///
    virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const;
    ///
    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    ///
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;
  
    virtual Table<int> * CreateSmoothingBlocks (const Flags & precflags) const;
    /// 
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const;
    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const;
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const;
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const;

    virtual bool VarOrder() const { return var_order; } 
    virtual int GetRelOrder() const { return rel_order; }   

  };

}

#endif

