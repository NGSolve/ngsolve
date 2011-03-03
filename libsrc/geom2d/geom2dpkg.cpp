#include <meshing.hpp>
#include <geometry2d.hpp>
#include <visual.hpp>

#include "vsgeom2d.hpp"

// extern "C" int Ng_CSG_Init (Tcl_Interp * interp);

namespace netgen
{
	

  extern DLL_HEADER NetgenGeometry * ng_geometry;
  static VisualSceneGeometry2d vsgeom2d;



  class SplineGeometryRegister : public GeometryRegister
  {
  public:
    virtual NetgenGeometry * Load (string filename) const;
    virtual VisualScene * GetVisualScene (const NetgenGeometry * geom) const;
  };


  NetgenGeometry *  SplineGeometryRegister :: Load (string filename) const
  {
    const char * cfilename = filename.c_str();
    if (strcmp (&cfilename[strlen(cfilename)-4], "in2d") == 0)
      {
	PrintMessage (1, "Load 2D-Spline geometry file ", cfilename);
	

	ifstream infile(cfilename);

	SplineGeometry2d * hgeom = new SplineGeometry2d();
	hgeom -> Load (cfilename);
	return hgeom;
      }
    
    return NULL;
  }



  VisualScene * SplineGeometryRegister :: GetVisualScene (const NetgenGeometry * geom) const
  {
    SplineGeometry2d * geometry = dynamic_cast<SplineGeometry2d*> (ng_geometry);
    if (geometry)
      {
	vsgeom2d.SetGeometry (geometry);
	return &vsgeom2d;
      }
    return NULL;
  }


}


using namespace netgen;
#ifdef WIN32
extern "C" __declspec(dllexport) int Ng_geom2d_Init (Tcl_Interp * interp);
#else
extern "C" int Ng_geom2d_Init (Tcl_Interp * interp);
#endif

int Ng_geom2d_Init (Tcl_Interp * interp)
{
  geometryregister.Append (new SplineGeometryRegister);
  return TCL_OK;
}
