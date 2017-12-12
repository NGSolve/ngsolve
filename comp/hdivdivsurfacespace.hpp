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

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const;
    virtual void GetDofNrs(ElementId id, Array<int> & dnums) const;

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
