#ifndef FILE_H1HOFESPACE
#define FILE_H1HOFESPACE

/*********************************************************************/
/* File:   h1hofespace.hpp                                           */
/* Author: Start                                                     */
/* Date:   10. Feb. 2003                                             */
/*********************************************************************/

#include "fespace.hpp"

namespace ngcomp
{

  /**
     High Order Finite Element Space
  */

  class NGS_DLL_HEADER H1HighOrderFESpace : public FESpace
  {
  protected:
    int level;

    Array<DofId> first_edge_dof;
    Array<DofId> first_face_dof;
    Array<DofId> first_element_dof;

    typedef unsigned char TORDER;

    bool fixed_order;
    bool wb_loedge, wb_edge;
    Array<TORDER> order_edge;
    Array<IVec<2,TORDER>> order_face;
    Array<IVec<3,TORDER>> order_inner;
    Array<bool> used_vertex; 
    Array<bool> used_edge; 
    Array<bool> used_face; 

    int uniform_order_quad;
    int uniform_order_trig;
    Array<IVec<3>> dom_order_min; 
    Array<IVec<3>> dom_order_max;
  
    bool level_adapted_order; 
    bool nodalp2;
    bool nodal;
    bool highest_order_dc;
    bool test_ho_prolongation;

    // D^-1 of full-order trigs, quads, tets, per vertex-order class
    mutable std::atomic<double*> dual_inverse[3][64] {};
    void ClearDualInverse ();
  public:

    H1HighOrderFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false);
    ///
    virtual ~H1HighOrderFESpace ();

    static DocInfo GetDocu ();

    
    virtual string GetClassName () const override
    {
      return "H1HighOrderFESpace";
    }

    ///
    void Update() override;

    virtual void DoArchive (Archive & archive) override;
    Array<MemoryUsage> GetMemoryUsage () const override;

    ///
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    ///
    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetFE (int elnr, Allocator & alloc) const;
    ///
    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetSFE (int elnr, Allocator & alloc) const;
    ///
    template <ELEMENT_TYPE ET>
      FiniteElement & T_GetCD2FE(int elnr, Allocator & alloc) const;
    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;

    virtual void GetDofRanges (ElementId ei, Array<IntRange> & dranges) const;

    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override;
    ///
    virtual shared_ptr<Table<int>> CreateSpecificSmoothingBlocks (const Flags & precflags) const override; 
    // virtual void CreateSmoothingBlocks2 (SmoothingBlocksCreator & sbc, const Flags & precflags) const; 
    ///
    virtual shared_ptr<Array<int>> CreateDirectSolverClusters (const Flags & flags) const override;

    virtual void UpdateDofTables () override;
    virtual void UpdateCouplingDofArray() override;    

    virtual void SetOrder (NodeId ni, int order) override;
    virtual int GetOrder (NodeId ni) const override;
    using FESpace::GetOrder;
    void SetOrder (const Region & region, int order) override { SetRegionOrder (region, max(order, 1)); }
    void SetOrder (const Region & region, IVec<2> order, shared_ptr<CoefficientFunction> direction) override
    { SetRegionOrder (region, max(order[0], 1), max(order[1], 1), direction); }
    void SetOrder (ElementId ei, IVec<3> order) override
    { SetCellOrder (ei, IVec<3> (max(order[0],1), max(order[1],1), max(order[2],1))); }
    IVec<3> GetOrder (ElementId ei) const override
    { bool elementwise; return CellBaseOrder (ei, ma->GetElIndex(ei), elementwise); }
    ORDER_POLICY GetOrderPolicy () const override { return DerivedOrderPolicy(); }
    
    void SetEdgeOrder (int enr, int eo) { order_edge[enr] = eo; }
    void SetFaceOrder (int fnr, IVec<2> fo) { order_face[fnr] = fo; }
    void SetElementOrder (int elnr, IVec<3> elo) { order_inner[elnr] = elo; }

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;

    auto GetEdgeDofs (size_t nr) const
    {
      return Range (first_edge_dof[nr], first_edge_dof[nr+1]);
    }

    auto GetFaceDofs (size_t nr) const
    {
      return Range (first_face_dof[nr], first_face_dof[nr+1]);
    }

    auto GetElementDofs (size_t nr) const
    {
      return Range (first_element_dof[nr], first_element_dof[nr+1]);
    }

  };


  class VectorH1FESpace : public CompoundFESpace
  {
    bool interleaved;
  public:
    VectorH1FESpace (shared_ptr<MeshAccess> ama, const Flags & flags, 
                     bool checkflags = false);

    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    static DocInfo GetDocu ();

    virtual void SetOrder (ELEMENT_TYPE et, TORDER order) override;

    virtual void SetOrder (NodeId ni, int order) override;
    void SetOrder (const Region & region, int order) override;
    void SetOrder (const Region & region, IVec<2> order, shared_ptr<CoefficientFunction> direction) override;
    void SetOrder (ElementId ei, IVec<3> order) override;
    IVec<3> GetOrder (ElementId ei) const override
    { bool elementwise; return CellBaseOrder (ei, ma->GetElIndex(ei), elementwise); }
    // single node orders are in the components
    const OrderSpecification & GetOrderSpecification () const override
    { return spaces[0]->GetOrderSpecification(); }
    void SetOrderSpecification (const OrderSpecification & spec, bool creator) override
    {
      FESpace::SetOrderSpecification (spec, creator);
      for (auto & spc : spaces) spc->SetOrderSpecification (spec, false);
    }
    virtual int GetOrder (NodeId ni) const override;
    using FESpace::GetOrder;

    virtual void FinalizeUpdate() override;

    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    void GetDofNrs (NodeId ni, Array<DofId> & dnums) const override;

    virtual string GetClassName () const override
    {
      return "VectorH1FESpace";
    }
  };

}

#endif

