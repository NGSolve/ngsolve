#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

  DLL_HEADER void Optimize2d (Mesh & mesh, MeshingParameters & mp)
  {
    static int timer = NgProfiler::CreateTimer ("optimize2d");
    NgProfiler::RegionTimer reg(timer);

    mesh.CalcSurfacesOfNode();

    const char * optstr = mp.optimize2d.c_str();
    int optsteps = mp.optsteps2d;

    for (int i = 1; i <= optsteps; i++)
      for (size_t j = 1; j <= strlen(optstr); j++)
	{
	  if (multithread.terminate) break;
	  switch (optstr[j-1])
	    {
	    case 's': 
	      {  // topological swap
		MeshOptimize2d meshopt;
                meshopt.SetMetricWeight (mp.elsizeweight);
		meshopt.EdgeSwapping (mesh, 0);
		break;
	      }
	    case 'S': 
	      {  // metric swap
		MeshOptimize2d meshopt;
                meshopt.SetMetricWeight (mp.elsizeweight);
		meshopt.EdgeSwapping (mesh, 1);
		break;
	      }
	    case 'm': 
	      {
		MeshOptimize2d meshopt;
                meshopt.SetMetricWeight (mp.elsizeweight);
		meshopt.ImproveMesh(mesh, mp);
		break;
	      }
	    case 'c': 
	      {
		MeshOptimize2d meshopt;
                meshopt.SetMetricWeight (mp.elsizeweight);
		meshopt.CombineImprove(mesh);
		break;
	      }
	    default:
	      cerr << "Optimization code " << optstr[j-1] << " not defined" << endl;
	    }  
	}
  }

}
