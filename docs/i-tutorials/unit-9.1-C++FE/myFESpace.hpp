#ifndef FILE_MYFESPACE_HPP
#define FILE_MYFESPACE_HPP


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
      Arguments are the MeshAccess view of the mesh data structure,
      and the kwargs from the Python constructor converted to C++ Flags.
    */
    MyFESpace (shared_ptr<MeshAccess> ama, const Flags & flags);

    // a name for our new fe-space
    string GetClassName () const override { return "MyFESpace"; }

    // documentation
    static DocInfo GetDocu();

    // organzize the FESpace, called after every mesh update
    void Update() override;
    
    // dof-numbers for element-id ei
    void GetDofNrs (ElementId ei, Array<DofId> & dnums) const override;
    
    // generate FiniteElement for element-id ei
    FiniteElement & GetFE (ElementId ei, Allocator & alloc) const override;

    // some new functionality our space should have in Python
    size_t GetNVert() { return nvert; }
  };

}    

void ExportMyFESpace(py::module m);


#endif
