/*********************************************************************/
/* File:   hdivdivfespace.h                                          */
/* Author: Joachim (cpp), Astrid                                     */
/* Date:   2017                                                      */
/*********************************************************************/

namespace ngcomp
{

  typedef size_t index_edge;

  class HDivDivFESpace : public FESpace
  {
    size_t ndof;
    Array<int> first_facet_dof;
    Array<int> first_element_dof;
    Array<INT<2,int> > order_facet;
    Array<INT<3,int> > order_inner;
    
    // add divdiv-free inner bubbles
    bool plus;
    bool discontinuous;
    int uniform_order_facet;
    int uniform_order_inner;

    Array<bool> fine_facet;

    bool quadfullpol;  // use quad elements with full polynomial space, but divdiv larger than L^{k+1}*
    bool algebraic_mapping;  // shapes are mapped algebraically from the reference element,
    // if false, then shapes are derivatives of mapped potentials
  public:
    HDivDivFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false);
    static DocInfo GetDocu ();

    virtual string GetClassName () const override
    {
      return "HDivDivFESpace";
    }

    void Update() override;

    virtual size_t GetNDof () const throw() override { return ndof; }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    virtual void SetOrder (NodeId ni, int order) override;
    virtual int GetOrder (NodeId ni) const override;

    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const override
    {
      dnums.SetSize0();
    }
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const override;

    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const override;
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const override;

    void GetDofNrs (ElementId ei, Array<int> & dnums) const override;
    
    virtual void UpdateCouplingDofArray() override;

    virtual SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const override;
  };

}
