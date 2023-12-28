#ifndef FILE_HIDDEN_
#define FILE_HIDDEN_

/*********************************************************************/
/* File:   hidden.hpp                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   Feb 2021                                                  */
/*********************************************************************/


#include "fespace.hpp"

namespace ngcomp
{

 // Hide all dofs of a space 

  class HiddenFESpace : public FESpace
  {
  protected:
    shared_ptr<FESpace> space;
    
  public:
    HiddenFESpace (shared_ptr<FESpace> space, const Flags & flags);
    
    virtual ~HiddenFESpace () { ; }
    void Update () override;
    
    void FinalizeUpdate() override
    {
      space->FinalizeUpdate();
      FESpace::FinalizeUpdate();
    }

    virtual string GetClassName() const override { return "Hidden" + space->GetClassName(); }
    shared_ptr<FESpace> GetBaseSpace() const { return space; }
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    virtual void GetDofNrs (NodeId ni, Array<DofId> & dnums) const override;

    virtual SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const override
    { return space->GetAdditionalEvaluators (); }

    [[deprecated("Use GetDofNrs(NODE_TYPE(NT_VERTEX,nr) instead")]]    
    virtual void GetVertexDofNrs (int vnr,  Array<DofId> & dnums) const override;

    [[deprecated("Use GetDofNrs(NODE_TYPE(NT_EDGE,nr) instead")]]        
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override;
    
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override;
    
    virtual void GetInnerDofNrs (int elnr, Array<DofId> & dnums) const override
    { space->GetInnerDofNrs(elnr, dnums); }

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


#endif // FILE_HIDDEN_
