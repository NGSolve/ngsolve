#ifndef FILE_PERIODIC_
#define FILE_PERIODIC_

namespace ngcomp
{

  class PeriodicFESpace : public FESpace
  {
    Array<int> dofmap; // mapping of dofs
    shared_ptr<FESpace> space;
    
  public:
    PeriodicFESpace (shared_ptr<FESpace> space, const Flags & flags);
    
    virtual ~PeriodicFESpace () { ; }
    virtual void Update (LocalHeap & lh) override;
    virtual void FinalizeUpdate (LocalHeap & lh) override { space->FinalizeUpdate(); }
    
    
    virtual FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    virtual size_t GetNDof () const { return space->GetNDof(); }

    virtual void GetDofNrs(ElementId ei, Array<int> & dnums) const override;
    
  private:
    void GetPeriodicNodeIds(Array<std::tuple<NodeId,NodeId>> & node_ids,int idnr) const;

  };
}


#endif // FILE_PERIODIC_
