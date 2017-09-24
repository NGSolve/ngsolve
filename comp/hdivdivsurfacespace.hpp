#ifndef FILE_HDIVSYMSURFACESPACE
#define FILE_HDIVSYMSURFACESPACE


/*********************************************************************/
/* File:   hdivdivsurfacespace.hpp                                   */
/* Author: Michael Neunteufel, Joachim Schöberl, Astrid Pechstein    */
/* Date:   April 2016                                                */
/*********************************************************************/

namespace ngcomp
{
  
  class HDivDivSurfaceSpace : public FESpace
  {
  private:
    int nfa;
    int nel;

    int order;
    Array<int> first_face_dof;
    Array<int> first_element_dof;
    int ndof;

    Array<bool> fine_face; 

    int noncontinuous;
  public:

    HDivDivSurfaceSpace(shared_ptr<MeshAccess> ama, const Flags & aflags, bool parseflags = false);

    virtual ~HDivDivSurfaceSpace ();

    virtual string GetClassName () const
    {
      return "HDivSymSurfaceSpace";
    }

    virtual void Update(LocalHeap & lh);
  
    virtual size_t GetNDof () const
    {
      return ndof;
    }

    int GetFirstFaceDof (const int f) const { return first_face_dof[f]; }
    int GetFirstCellDof (const int c) const { return first_element_dof[c]; }


    // new style
    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const
    {
      if (ei.IsVolume())
        return GetFE_old(ei.Nr(), dynamic_cast<LocalHeap&> (lh));
      else if (ei.VB() == BND)
        return GetSFE_old(ei.Nr(), dynamic_cast<LocalHeap&> (lh));
      else if (ei.VB() == BBND)
        return * new (lh) DummyFE<ET_SEGM>();
      else
        return * new (lh) DummyFE<ET_POINT>();
    }
  
    virtual FiniteElement & GetFE_old (int elnr, LocalHeap & lh) const;
    virtual FiniteElement & GetSFE_old (int elnr, LocalHeap & lh) const;

    virtual void GetDofNrs(ElementId id, Array<int> & dnums) const
    {
      if (id.IsVolume()) GetDofNrs(id.Nr(), dnums);
      else if (id.IsBoundary()) GetSDofNrs(id.Nr(), dnums);
      else if (id.VB() == BBND) GetEdgeDofNrs(ma->GetElEdges(id)[0], dnums);
      else dnums.SetSize0();
    }
  
    virtual void GetDofNrs(int elnr, Array<int> & dnums) const;
    virtual void GetSDofNrs (int elnr, Array<int> & dnums) const;
 
    virtual void GetWireBasketDofNrs (int elnr, Array<int> & dnums) const
    { dnums.SetSize(0); }
    virtual void GetVertexDofNrs (int vnr, Array<int> & dnums) const
    { dnums.SetSize(0); }
    virtual void GetEdgeDofNrs (int ednr, Array<int> & dnums) const
    { dnums = IntRange(first_face_dof[ednr],first_face_dof[ednr+1]); }
    virtual void GetFaceDofNrs (int fanr, Array<int> & dnums) const
    { dnums.SetSize(0); }
    virtual void GetInnerDofNrs (int elnr, Array<int> & dnums) const
    { GetDofNrs (elnr, dnums); }
    virtual std::shared_ptr<ngmg::Prolongation> GetProlongation() const
    { return 0; }
  };

}

#endif
