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

    Array<bool> fine_facet;

    //bool hiddeneldofs;
    bool alllocaldofs;
    bool discontinuous;

    // add curl of Nedelec bubbles
    //bool curlbubbles;
    // GG bubbles
    bool GGbubbles;    
    
    int uniform_order_facet;
    int uniform_order_inner;
    int uniform_order_trace;

  public:
    HCurlDivFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false);

    virtual string GetClassName () const override
    {
      return "HCurlDiv FESpace";
    }

    static DocInfo GetDocu ();

    void Update() override;

    virtual size_t GetNDof () const throw() override { return ndof; }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;
    

    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const override
    {
      dnums.SetSize0();
    }
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const override;

    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const override;
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const override;
    
    virtual void GetLoDofNrs (int elnr, Array<int> & dnums) const;

    virtual void GetFacetDofNrs(int fanr, Array<DofId> & dnums) const
    { 
      if (ma->GetDimension() == 2) GetEdgeDofNrs(fanr,dnums); 
      else if (ma->GetDimension() == 3) GetFaceDofNrs(fanr,dnums); 
    } 

    void GetDofNrs (ElementId ei, Array<int> & dnums) const override;
    
    virtual void UpdateCouplingDofArray() override;

    virtual SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const override;

    
  };

}
