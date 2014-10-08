#include <meshing.hpp>
#include <geometry2d.hpp>
#include <visual.hpp>
#include <inctcl.hpp>

#include "vsgeom2d.hpp"

// extern "C" int Ng_CSG_Init (Tcl_Interp * interp);

namespace netgen
{
	

  // extern DLL_HEADER NetgenGeometry * ng_geometry;
  static VisualSceneGeometry2d vsgeom2d;





  class SplineGeometryVisRegister : public GeometryRegister
  {
  public:
    virtual NetgenGeometry * Load (string filename) const { return NULL; }
    virtual VisualScene * GetVisualScene (const NetgenGeometry * geom) const;
  };


  VisualScene * SplineGeometryVisRegister :: GetVisualScene (const NetgenGeometry * geom) const
  {
    const SplineGeometry2d * geometry = dynamic_cast<const SplineGeometry2d*> (geom);
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
  geometryregister.Append (new SplineGeometryVisRegister);
  return TCL_OK;
}
