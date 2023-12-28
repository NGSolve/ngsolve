#ifndef FILE_REORDERED_
#define FILE_REORDERED_

/*********************************************************************/
/* File:   reordered.hpp                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   Jun 2019                                                  */
/*********************************************************************/

#include "fespace.hpp"

namespace ngcomp
{

 // A reordered wrapper class for fespaces 

  class ReorderedFESpace : public FESpace
  {
  protected:
    Array<DofId> dofmap;
    shared_ptr<FESpace> space;
    shared_ptr<Table<DofId>> clusters;
    
  public:
    ReorderedFESpace (shared_ptr<FESpace> space, const Flags & flags);
    
    virtual ~ReorderedFESpace () { ; }

    virtual void Update() override;

    virtual void FinalizeUpdate() override;

    ProxyNode MakeProxyFunction (bool testfunction,
                                 const function<shared_ptr<ProxyFunction>(shared_ptr<ProxyFunction>)> & addblock) const override
    {
      return GetBaseSpace()->MakeProxyFunction (testfunction, addblock);
    }
    
    virtual string GetClassName() const override { return "Reordered" + space->GetClassName(); }
    shared_ptr<FESpace> GetBaseSpace() const { return space; }
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    virtual void GetDofNrs (NodeId ni, Array<DofId> & dnums) const override;

    auto GetClusters() const { return clusters; }

    
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


#endif // FILE_REORDERED_
