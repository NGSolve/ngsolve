#ifndef FILE_PERIODIC_
#define FILE_PERIODIC_

/*********************************************************************/
/* File:   periodic.hpp                                              */
/* Author: Christopher Lackner                                       */
/* Date:   Feb. 2017                                                 */
/*********************************************************************/

namespace ngcomp
{

 // A periodic wrapper class for fespaces 

  class PeriodicFESpace : public FESpace
  {
  protected:
    Array<int> dofmap; // mapping of dofs
    shared_ptr<FESpace> space;
    shared_ptr<Array<int>> used_idnrs;
    
  public:
    PeriodicFESpace (shared_ptr<FESpace> space, const Flags & flags, shared_ptr<Array<int>> aused_idnrs);
    
    virtual ~PeriodicFESpace () { ; }
    virtual void Update (LocalHeap & lh) override;
    
    virtual void FinalizeUpdate (LocalHeap & lh) override {
      space->FinalizeUpdate(lh);
      FESpace::FinalizeUpdate(lh);
    }

    virtual string GetClassName() const override { return "Periodic" + space->GetClassName(); }
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    virtual size_t GetNDof () const override { return space->GetNDof(); }
    virtual size_t GetNDofLevel (int level) const override { return space->GetNDofLevel(level); }

    virtual void GetDofNrs(ElementId ei, Array<int> & dnums) const override;

    virtual SymbolTable<shared_ptr<DifferentialOperator>> GetAdditionalEvaluators () const override
    { return space->GetAdditionalEvaluators (); }

    virtual void GetVertexDofNrs (int vnr,  Array<DofId> & dnums) const override
    { space->GetVertexDofNrs(vnr, dnums); }
    virtual void GetEdgeDofNrs (int ednr, Array<DofId> & dnums) const override
    { space->GetEdgeDofNrs (ednr, dnums); }
    virtual void GetFaceDofNrs (int fanr, Array<DofId> & dnums) const override
    { space->GetFaceDofNrs(fanr, dnums); }
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

  protected:
    // overload in quasiperiodic space
    virtual void DofMapped(size_t from, size_t to, size_t idnr) { ; }
  };

  class QuasiPeriodicFESpace : public PeriodicFESpace
  {
    shared_ptr<Array<Complex>> factors;
    Array<Complex> dof_factors;

  public:
    QuasiPeriodicFESpace (shared_ptr<FESpace> fespace, const Flags & flag, shared_ptr<Array<int>> aused_idnrs, shared_ptr<Array<Complex>> afactors);

    virtual void Update (LocalHeap & lh) override;

    virtual void VTransformMR (ElementId ei, SliceMatrix<double> mat, TRANSFORM_TYPE tt) const override;
    virtual void VTransformMC (ElementId ei, SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const override;
    virtual void VTransformVR (ElementId ei, SliceVector<double> vec, TRANSFORM_TYPE tt) const override;
    virtual void VTransformVC (ElementId ei, SliceVector<Complex> vec, TRANSFORM_TYPE tt) const override;


    
  protected:
    virtual void DofMapped(size_t from, size_t to, size_t idnr) override;
  };
  

}


#endif // FILE_PERIODIC_
