#ifndef FILE_DISCONTINUOUS_
#define FILE_DISCONTINUOUS_

/*********************************************************************/
/* File:   discontinuous.hpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   Jun 2019                                                  */
/*********************************************************************/

#include "fespace.hpp"

namespace ngcomp
{

 // A discontinuous wrapper class for fespaces 

  class DiscontinuousFESpace : public FESpace
  {
  protected:
    Array<DofId> first_element_dof; 
    shared_ptr<FESpace> space;
    VorB vb;
    
  public:
    DiscontinuousFESpace (shared_ptr<FESpace> space, const Flags & flags);
    
    virtual ~DiscontinuousFESpace () { ; }
    void Update () override;
    
    void FinalizeUpdate() override
    {
      space->FinalizeUpdate();
      FESpace::FinalizeUpdate();
    }

    virtual FlatArray<VorB> GetDualShapeNodes (VorB vb) const override
    {
      return space->GetDualShapeNodes(vb);
    }


    virtual string GetClassName() const override { return "Discontinuous" + space->GetClassName(); }
    shared_ptr<FESpace> GetBaseSpace() const { return space; }
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    // virtual size_t GetNDof () const override { return space->GetNDof(); }
    // virtual size_t GetNDofLevel (int level) const override { return space->GetNDofLevel(level); }

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


#endif // FILE_DISCONTINUOUS_
