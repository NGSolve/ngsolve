#pragma once
/*********************************************************************/
/* File:   wrapperfespace.hpp                                        */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   12. July 2018                                             */
/*********************************************************************/


#include "fespace.hpp"

namespace ngcomp
{

  class NGS_DLL_HEADER CompressedFESpace : public FESpace
  {
  protected:
    shared_ptr<FESpace> space;
    Array<DofId> comp2all;
    Array<DofId> all2comp;
    shared_ptr<BitArray> active_dofs = nullptr;

    //TODO: flag to hide HIDDEN_DOFs

  public:
    CompressedFESpace (shared_ptr<FESpace> bfes);
    virtual ~CompressedFESpace () {};
    void Update() override;
    shared_ptr<FESpace> GetBaseSpace() const { return space; }

    void WrapDofs(Array<DofId> & dnums) const
    {
      /*
      // (JS): this is buggy: dnums[i] may be -1
      for (int i : Range(dnums.Size()))
        if (all2comp[dnums[i]] != -1)
          dnums[i] = all2comp[dnums[i]];
        else 
          dnums[i] = -1;
      */
      for (DofId & d : dnums)
        if (IsRegularDof (d))
          d = all2comp[d];
    }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    
    FlatArray<VorB> GetDualShapeNodes (VorB vb) const override;

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    virtual void GetDofNrs (NodeId ni, Array<DofId> & dnums) const override;

    virtual SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const override
    { return space->GetAdditionalEvaluators (); }

    /// get dof-nrs of the element of certain coupling type
    virtual void GetElementDofsOfType (ElementId ei, Array<DofId> & dnums, COUPLING_TYPE ctype) const override;

    virtual void SetActiveDofs(shared_ptr<BitArray> actdof) 
    {
      active_dofs = actdof;
    }

    shared_ptr<BitArray> GetActiveDofs() const { return active_dofs; }

    // a name for our new fe-space
    virtual string GetClassName () const override
    {
      return "CompressedFESpace(" + space->GetClassName() + ")";
    }

    /// update element coloring
    void FinalizeUpdate() override
    {
      space->FinalizeUpdate();
      FESpace::FinalizeUpdate();
    }

    ProxyNode MakeProxyFunction (bool testfunction,
                                 const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const override;

    virtual bool DefinedOn (ElementId id) const
    {
      return space->DefinedOn(id);
    }

    virtual void GetVertexDofNrs (int vnr, Array<DofId> & dnums) const override;
    /// get dofs on edge enr
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
    /// get dofs on face fnr
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
    /// get dofs on element (=cell) elnr
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override;

    virtual void VTransformMR (ElementId ei,
			       SliceMatrix<double> mat, TRANSFORM_TYPE tt) const override
    { space-> VTransformMR(ei, mat, tt); }
    virtual void VTransformMC (ElementId ei,
                               SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const override
    { space->VTransformMC (ei, mat, tt); }
    virtual void VTransformVR (ElementId ei,
                               SliceVector<double> vec, TRANSFORM_TYPE tt) const override
    { space->VTransformVR(ei, vec, tt); }
    virtual void VTransformVC (ElementId ei, 
                               SliceVector<Complex> vec, TRANSFORM_TYPE tt) const override
    { space->VTransformVC(ei, vec, tt); }    


  };

}
