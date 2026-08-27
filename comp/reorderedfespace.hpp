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

  // element orderings for locality: geometric (Morton curve of centroids) and
  // topological (reverse Cuthill-McKee on the vertex-sharing element graph).
  // LocalityElementOrder is the common entry point used for first-touch dof
  // numbering and locality-aware element batching -- keep all users on it so
  // dof numbering and batching stay consistent.
  NGS_DLL_HEADER Array<int> MortonElementOrder (const shared_ptr<MeshAccess> & ma);
  NGS_DLL_HEADER Array<int> RCMElementOrder (const shared_ptr<MeshAccess> & ma);
  NGS_DLL_HEADER Array<int> LocalityElementOrder (const shared_ptr<MeshAccess> & ma);

 // A reordered wrapper class for fespaces

  class ReorderedFESpace : public FESpace
  {
  protected:
    Array<DofId> dofmap;
    Array<int> elorder;   // space-filling-curve element order used for first-touch numbering
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
      // build the proxy (tree) via the base space (keeps compound structure),
      // but bind every proxy to this wrapper space, so that assembly uses the
      // reordered dof numbering (same pattern as PeriodicFESpace)
      return GetBaseSpace()->MakeProxyFunction
        (testfunction,
         [&] (shared_ptr<ProxyFunction> proxy)
         {
           shared_ptr<FESpace> fes = dynamic_pointer_cast<FESpace> (const_cast<ReorderedFESpace*>(this)->shared_from_this());
           proxy->SetFESpace(fes);
           return addblock (proxy);
         });
    }
    
    virtual string GetClassName() const override { return "Reordered" + space->GetClassName(); }
    shared_ptr<FESpace> GetBaseSpace() const { return space; }
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    virtual void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    virtual void GetDofNrs (NodeId ni, Array<DofId> & dnums) const override;

    auto GetClusters() const { return clusters; }

    // element processing order (Morton curve); process elements in this order
    // to profit from the first-touch dof numbering
    FlatArray<int> GetElementOrder() const { return elorder; }

    
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
