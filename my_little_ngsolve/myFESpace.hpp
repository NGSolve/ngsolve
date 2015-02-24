#ifndef FILE_MYFESPACE_HPP
#define FILE_MYFESPACE_HPP

/*********************************************************************/
/* File:   myFESpace.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*

  My own FESpace for linear and quadratic triangular elements An
  fe-space is the connection between the local reference element, and
  the global mesh.

*/


namespace ngcomp
{

  class MyFESpace : public FESpace
  {
    bool secondorder;
    int ndof, nvert;

    FiniteElement * reference_element;
    FiniteElement * reference_surface_element;
    
  public:
    /*
      constructor. 
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
    */
    MyFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

    // destructor
    virtual ~MyFESpace ();

    // a name for our new fe-space
    virtual string GetClassName () const
    {
      return "MyFESpace";
    }


    virtual void Update(LocalHeap & lh);
    virtual int GetNDof () const { return ndof; }

    virtual void GetDofNrs (int elnr, Array<int> & dnums) const;
    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const;

    virtual const FiniteElement & GetFE (int elnr, LocalHeap & lh) const;
    virtual const FiniteElement & GetSFE (int selnr, LocalHeap & lh) const;
  };

}    

#endif
