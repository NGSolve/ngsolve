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
    size_t nvert;
    
  public:
    /*
      constructor. 
      Arguments are the access to the mesh data structure,
      and the flags from the define command in the pde-file
      or the kwargs in the Python constructor.
    */
    MyFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

    // a name for our new fe-space
    string GetClassName () const override { return "MyFESpace"; }

    static DocInfo GetDocu();

    void Update() override;
    
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    // some new functionality our space should have in Python
    size_t GetNVert() { return nvert; }
  };

}    

void ExportMyFESpace(py::module m);


#endif
