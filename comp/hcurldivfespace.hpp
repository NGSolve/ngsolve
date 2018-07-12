#pragma once
/*********************************************************************/
/* File:   hcurldivfespace.h                                           */
/* Author: Philip                                                  */
/* Date:   2017/2018                                             */
/*********************************************************************/

namespace ngcomp
{


  typedef size_t index_edge;

    class HCurlDivFESpace : public FESpace
  {
    size_t ndof;
    Array<int> first_facet_dof;
    Array<int> first_element_dof;
    Array<int> order_facet;
    Array<int> order_inner;
    Array<int> order_trace;

    bool hiddeneldofs;
    
    // add curldiv-free inner bubbles
    bool discontinuous;
    //bool withtrace;
    
    int uniform_order_facet;
    int uniform_order_inner;
    int uniform_order_trace;

  public:
    HCurlDivFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false);

    virtual string GetClassName () const override
    {
      return "HCurlDiv FESpace";
    }

    virtual void Update(LocalHeap & lh) override;

    virtual size_t GetNDof () const throw() override { return ndof; }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;
    

    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const override
    {
      dnums.SetSize0();
    }
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const override;

    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const override;
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const override;

    void GetDofNrs (ElementId ei, Array<int> & dnums) const override;
    
    virtual void UpdateCouplingDofArray();

    virtual SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const override;

    
  };

}
