#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

  DLL_HEADER Array<GeometryRegister*> geometryregister;

  GeometryRegister :: ~GeometryRegister()
  { ; }




  
  int NetgenGeometry :: GenerateMesh (Mesh*& mesh, MeshingParameters & mparam,
				      int perfstepsstart, int perfstepsend)
  {
    if (!mesh) return 1;

    if (perfstepsstart <= MESHCONST_MESHVOLUME)
      {
	multithread.task = "Volume meshing";
	
	MESHING3_RESULT res =
	  MeshVolume (mparam, *mesh);
	
	if (res != MESHING3_OK) return 1;
	
	if (multithread.terminate) return 0;
	
	RemoveIllegalElements (*mesh);
	if (multithread.terminate) return 0;

	MeshQuality3d (*mesh);
      }

    
    if (multithread.terminate || perfstepsend <= MESHCONST_MESHVOLUME)
      return 0;


    if (perfstepsstart <= MESHCONST_OPTVOLUME)
      {
	multithread.task = "Volume optimization";
	
	OptimizeVolume (mparam, *mesh);
	if (multithread.terminate) return 0;
      }
    
    return 0;
  }    
  

  const Refinement & NetgenGeometry :: GetRefinement () const
  {
    return *new Refinement;;
  }


  void NetgenGeometry :: Save (string filename) const
  {
    throw NgException("Cannot save geometry - no geometry available");
  }

}
