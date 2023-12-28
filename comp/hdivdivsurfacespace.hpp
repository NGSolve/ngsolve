#ifndef FILE_HDIVSYMSURFACESPACE
#define FILE_HDIVSYMSURFACESPACE


/*********************************************************************/
/* File:   hdivdivsurfacespace.hpp                                   */
/* Author: Michael Neunteufel, Joachim Schöberl, Astrid Pechstein    */
/* Date:   April 2016                                                */
/*********************************************************************/


#include "fespace.hpp"

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
    static DocInfo GetDocu ();
    

    string GetClassName () const override
    {
      return "HDivSymSurfaceSpace";
    }

    void Update() override;

    void UpdateCouplingDofArray() override;
  
    size_t GetNDof () const override
    {
      return ndof;
    }

    int GetFirstFaceDof (const int f) const { return first_face_dof[f]; }
    int GetFirstCellDof (const int c) const { return first_element_dof[c]; }

    FiniteElement & GetFE (ElementId ei, Allocator & lh) const override;
    void GetDofNrs(ElementId id, Array<int> & dnums) const override;

    virtual void GetWireBasketDofNrs (int elnr, Array<int> & dnums) const
    { dnums.SetSize(0); }
    void GetVertexDofNrs (int vnr, Array<int> & dnums) const override
    { dnums.SetSize(0); }
    void GetEdgeDofNrs (int ednr, Array<int> & dnums) const override
    { dnums = IntRange(first_face_dof[ednr],first_face_dof[ednr+1]); }
    void GetFaceDofNrs (int fanr, Array<int> & dnums) const override
    { dnums.SetSize(0); }
    void GetInnerDofNrs (int elnr, Array<int> & dnums) const override
    { GetDofNrs (elnr, dnums); }
    std::shared_ptr<ngmg::Prolongation> GetProlongation() const override
    { return 0; }

  };

}

#endif
