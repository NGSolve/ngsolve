#ifndef FILE_BASEGEOM
#define FILE_BASEGEOM

/**************************************************************************/
/* File:   basegeom.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   23. Aug. 09                                                    */
/**************************************************************************/


struct Tcl_Interp;

namespace netgen
{

  class DLL_HEADER NetgenGeometry
  {
  public:
    virtual ~NetgenGeometry () { ; }

    virtual int GenerateMesh (shared_ptr<Mesh> & mesh, MeshingParameters & mparam, 
			      int perfstepsstart, int perfstepsend);

    virtual const Refinement & GetRefinement () const;

    virtual void Save (string filename) const;
    virtual void SaveToMeshFile (ostream & /* ost */) const { ; }
  };





  class DLL_HEADER GeometryRegister
  {
  public:
    virtual ~GeometryRegister();
    virtual NetgenGeometry * Load (string filename) const = 0;
    virtual NetgenGeometry * LoadFromMeshFile (istream & /* ist */) const { return NULL; }
    virtual class VisualScene * GetVisualScene (const NetgenGeometry * /* geom */) const
    { return NULL; }
    virtual void SetParameters (Tcl_Interp * /* interp */) { ; }
  };

  class DLL_HEADER GeometryRegisterArray : public Array<GeometryRegister*>
  {
  public:
    virtual ~GeometryRegisterArray()
    {
      for (int i = 0; i < Size(); i++)
        delete (*this)[i];
    }
  };

  // extern DLL_HEADER Array<GeometryRegister*> geometryregister; 
  extern DLL_HEADER GeometryRegisterArray geometryregister; 
}



#endif
