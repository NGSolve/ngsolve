#ifndef FILE_PERIODIC_
#define FILE_PERIODIC_

/*********************************************************************/
/* File:   periodic.hpp                                              */
/* Author: Christopher Lackner                                       */
/* Date:   Feb. 2017                                                 */
/*********************************************************************/

#include <set>
#include "fespace.hpp"

namespace ngcomp
{

 // A periodic wrapper class for fespaces 

  class NGS_DLL_HEADER PeriodicFESpace : public FESpace
  {
  protected:
    Array<int> dofmap; // mapping of dofs
    Array<int> vertex_map; // mapping of vertices
    shared_ptr<FESpace> space;
    shared_ptr<Array<int>> used_idnrs;
    
  public:
    PeriodicFESpace (shared_ptr<FESpace> space, const Flags & flags, shared_ptr<Array<int>> aused_idnrs);
    
    virtual ~PeriodicFESpace () { ; }
    void Update() override;
    void DoArchive (Archive & archive) override;
    auto GetCArgs () { return std::make_tuple(Shallow(space), GetFlags(), used_idnrs); }
    
    void FinalizeUpdate() override
    {
      space->FinalizeUpdate();
      FESpace::FinalizeUpdate();
    }

    ProxyNode MakeProxyFunction (bool testfunction,
                                 const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const override;
    
    shared_ptr<Array<int>> GetUsedIdnrs() const { return used_idnrs; }
    virtual string GetClassName() const override { return "Periodic" + space->GetClassName(); }
    shared_ptr<FESpace> GetBaseSpace() const { return space; }
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    virtual size_t GetNDof () const override { return space->GetNDof(); }
    virtual size_t GetNDofLevel (int level) const override { return space->GetNDofLevel(level); }

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

    auto & GetVertexMap() const { return vertex_map; }
  protected:
    // overload in quasiperiodic space
    virtual void DofMapped(size_t from, size_t to, size_t idnr) { ; }
  };

  template<typename TSCAL>
  class QuasiPeriodicFESpace : public PeriodicFESpace
  {
    shared_ptr<Array<TSCAL>> factors;
    Array<TSCAL> dof_factors;
    Array<set<size_t>> master_dofs;

  public:
    QuasiPeriodicFESpace (shared_ptr<FESpace> fespace, const Flags & flag, shared_ptr<Array<int>> aused_idnrs, shared_ptr<Array<TSCAL>> afactors);

    void Update() override;
    void DoArchive (Archive & archive) override {
      PeriodicFESpace::DoArchive(archive);
      archive & factors & dof_factors & master_dofs;
    }
    auto GetCArgs () {
      return std::make_tuple(Shallow(GetBaseSpace()), GetFlags(), used_idnrs, factors);
    }

    shared_ptr<Array<TSCAL>> GetFactors() const { return factors; }

    virtual void VTransformMR (ElementId ei, SliceMatrix<double> mat, TRANSFORM_TYPE tt) const override;
    virtual void VTransformMC (ElementId ei, SliceMatrix<Complex> mat, TRANSFORM_TYPE tt) const override;
    virtual void VTransformVR (ElementId ei, SliceVector<double> vec, TRANSFORM_TYPE tt) const override;
    virtual void VTransformVC (ElementId ei, SliceVector<Complex> vec, TRANSFORM_TYPE tt) const override;


    
  protected:
    virtual void DofMapped(size_t from, size_t to, size_t idnr) override;
  };
  

}


#endif // FILE_PERIODIC_
