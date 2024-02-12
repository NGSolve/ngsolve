#ifndef FILE_L2HOFESPACE
#define FILE_L2HOFESPACE

/*********************************************************************/
/* File:   l2hofespace.hpp                                           */
/* Author: Start                                                     */
/* Date:   23.Feb. 2003                                              */
/*********************************************************************/

#include "fespace.hpp"

namespace ngcomp
{
  
  /**
     High Order Finite Element Space for L2 (element by element)
  */

  class NGS_DLL_HEADER L2HighOrderFESpace : public FESpace
  {
  protected:
    /*
      ----> baseclass
    // Number of Elements
    int nel;
    // Degrees of Freedom 
    int ndof;
    // Levels 
    Array<int> ndlevel;
    */
      
    // if order is relative to mesh order 
    bool var_order;
    // variable order is set to mesh_order + rel_order 
    int rel_order;
    // order of elements 
    Array<IVec<3> > order_inner;
    // table of first element dofnumber 
    Array<DofId> first_element_dof;
    bool all_dofs_together;
    // set all used dofs to hidden_dofs
    bool hide_all_dofs;
    COUPLING_TYPE lowest_order_ct;
    bool tensorproduct;
  public:

    L2HighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~L2HighOrderFESpace ();

    static DocInfo GetDocu ();
    
    // Create if order=0 ElementFESpace Constructor, else L2HOFE 
    static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags);
    // Creates also for order=0 a L2HighOrderFESpace 
    static shared_ptr<FESpace> CreateHO (shared_ptr<MeshAccess> ma, const Flags & flags) 
    {
      return make_shared<L2HighOrderFESpace> (ma, flags);
    }  
  
    virtual string GetClassName () const override
    {
      return "L2HighOrderFESpace";
    }

    bool AllDofsTogether(){return all_dofs_together;};
    ///
    void Update() override;
    /// 
    virtual void UpdateDofTables() override;
    ///
    virtual void UpdateCouplingDofArray() override;    
    ///
    // virtual size_t GetNDof () const throw() override;
    ///
    // virtual size_t GetNDofLevel (int level) const override;
    ///
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    using FESpace::GetFE;
    // virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const override;
    // ///
    // virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const override;
    ///
    virtual const FiniteElement & GetFacetFE (int fnr, LocalHeap & lh) const;

    virtual void GetDofRanges (ElementId ei, Array<IntRange> & dranges) const;

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

    virtual void SetOrder (NodeId ni, int order) override;
    virtual int GetOrder (NodeId ni) const override;
    using FESpace::GetOrder;
    
    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;

    ///
    virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & precflags) const override;
    /// 
 
    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override;


    auto GetElementDofs (size_t nr) const
    {
      return Range (first_element_dof[nr], first_element_dof[nr+1]);
    }

    virtual shared_ptr<BaseMatrix> GetMassOperator (shared_ptr<CoefficientFunction> rho,
                                                    shared_ptr<Region> defon,
                                                    LocalHeap & lh) const override;

    virtual shared_ptr<BaseMatrix> CreateMassOperator (shared_ptr<CoefficientFunction> rho,
                                                       shared_ptr<Region> defon,
                                                       bool inverse,
                                                       LocalHeap & lh) const override;
    
    virtual void SolveM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                         LocalHeap & lh) const override;
    virtual void ApplyM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                         LocalHeap & lh) const override;

    virtual shared_ptr<BaseMatrix> GetTraceOperator (shared_ptr<FESpace> tracespace, bool avg) const override;
    virtual void GetTrace (const FESpace & tracespace, const BaseVector & in, BaseVector & out, bool avg,
                           LocalHeap & lh) const override;
    virtual void GetTraceTrans (const FESpace & tracespace, const BaseVector & in, BaseVector & out, bool avg,
                                LocalHeap & lh) const override;

  protected:

    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetFE (int elnr, Allocator & alloc) const;
  };












  class NGS_DLL_HEADER L2SurfaceHighOrderFESpace : public FESpace
  {
  protected:
    // Level
    // int level;

    // Number of Elements
    // int nel;

    Array<int> first_element_dof;
    // int ndof;


    // if order is relative to mesh order 
    bool var_order;
    // variable order is set to mesh_order + rel_order 
    int rel_order;
    // order of elements 
    Array<IVec<3> > order_inner;

    bool lowest_order_wb;
    bool discontinuous;
    bool dual_mapping; // u(x) = 1/measure * u(hatx)
  public:

    L2SurfaceHighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool parseflags=false);
    ///
    virtual ~L2SurfaceHighOrderFESpace ();

    static DocInfo GetDocu ();

    // static shared_ptr<FESpace> Create (shared_ptr<MeshAccess> ma, const Flags & flags);

    virtual string GetClassName () const override
    {
      return "L2SurfaceHighOrderFESpace";
    }

    ///
    void Update() override;
    /// 
    virtual void UpdateCouplingDofArray() override;    
    //virtual void UpdateDofTables() override;
    ///
    // virtual size_t GetNDof () const throw() override;

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    ///
    // virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const override;
    // ///
    // virtual const FiniteElement & GetSFE (int elnr, LocalHeap & lh) const override;
    ///
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
  
    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;

    virtual shared_ptr<BaseMatrix> GetMassOperator (shared_ptr<CoefficientFunction> rho,
                                                    shared_ptr<Region> defon,
                                                    LocalHeap & lh) const override;
    
    virtual void SolveM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                         LocalHeap & lh) const override;
    virtual void ApplyM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                         LocalHeap & lh) const override;
    virtual shared_ptr<Table<int>> CreateSmoothingBlocks (const Flags & precflags) const override;
    /// 
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override;
    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;

    virtual bool VarOrder() const override { return var_order; } 
    virtual int GetRelOrder() const override { return rel_order; }

    virtual void SetOrder (NodeId ni, int order) override;
    virtual int GetOrder (NodeId ni) const override;
    using FESpace::GetOrder;
    
    auto GetElementDofs (size_t nr) const
    {
      return Range (first_element_dof[nr], first_element_dof[nr+1]);
    }

  };


  class VectorL2FESpace : public CompoundFESpace
  {
    bool piola = false;
    bool covariant = false;
  public:
    VectorL2FESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags = false);

    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    static DocInfo GetDocu ();
    
    void GetDofNrs (ElementId ei, Array<int> & dnums) const override;

    virtual shared_ptr<BaseMatrix> GetMassOperator (shared_ptr<CoefficientFunction> rho,
                                                    shared_ptr<Region> defon,
                                                    LocalHeap & lh) const override;

    virtual shared_ptr<BaseMatrix> CreateMassOperator (shared_ptr<CoefficientFunction> rho,
                                                       shared_ptr<Region> defon,
                                                       bool inverse,
                                                       LocalHeap & lh) const override;

    template <int DIM>
    shared_ptr<BaseMatrix> CreateMassOperator_Dim (shared_ptr<CoefficientFunction> rho,
                                                   shared_ptr<Region> defon,
                                                   bool inverse,
                                                   LocalHeap & lh) const;
    
    virtual void SolveM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                         LocalHeap & lh) const override;
    virtual void ApplyM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                         LocalHeap & lh) const override;

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;

    template <int DIM>
    void SolveM_Dim (CoefficientFunction * rho, BaseVector & vec, Region * region,
                     LocalHeap & lh) const;

    /*
    template <int DIM>
    void SolveMPiola (CoefficientFunction * rho, BaseVector & vec,
                      LocalHeap & lh) const;
    template <int DIM>
    void SolveMCovariant (CoefficientFunction * rho, BaseVector & vec,
                          LocalHeap & lh) const;
    */
    
    template <int DIM>
    void ApplyM_Dim (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                     LocalHeap & lh) const;

    template <int DIM>
    void ApplyMPiola (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                      LocalHeap & lh) const;
    template <int DIM>
    void ApplyMCovariant (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                          LocalHeap & lh) const;

    virtual string GetClassName () const override
    {
      return "VectorL2FESpace";
    }

  };







  class TangentialSurfaceL2FESpace : public CompoundFESpace
  {
    bool piola;
  public:
    TangentialSurfaceL2FESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags = false);

    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    static DocInfo GetDocu ();
    
    void GetDofNrs (ElementId ei, Array<int> & dnums) const override;

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;

    virtual shared_ptr<BaseMatrix> GetMassOperator (shared_ptr<CoefficientFunction> rho,
                                                    shared_ptr<Region> defon,
                                                    LocalHeap & lh) const override;
    
    virtual void SolveM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                         LocalHeap & lh) const override;
    virtual void ApplyM (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                         LocalHeap & lh) const override;


    template <int DIM>
    void SolveM_Dim (CoefficientFunction * rho, BaseVector & vec, Region * region,
                     LocalHeap & lh) const;

    template <int DIM>
    void ApplyM_Dim (CoefficientFunction * rho, BaseVector & vec, Region * definedon,
                     LocalHeap & lh) const;

    virtual string GetClassName () const override
    {
      return "TangentialSurfaceL2FESpace";
    }

  };



  
}

#endif

