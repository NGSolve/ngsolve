#ifdef PARALLEL

#include "dlfcn.h"


// #include <mystdlib.h>

#include <meshing.hpp>

// #include <visual.hpp>


#include <meshing.hpp>

void (*NGS_ParallelRun) (const string & message) = NULL;


namespace netgen
{
#include "../interface/writeuser.hpp"
  extern string ngdir;
}

void Parallel_Exit();


namespace netgen {
  extern AutoPtr<Mesh>  mesh;
  // extern VisualSceneMesh vsmesh;
  extern DLL_HEADER MeshingParameters mparam;
}

using namespace netgen;
using netgen::RegisterUserFormats;

namespace netgen
{
  // int id, ntasks;
  MPI_Comm mesh_comm;
}


void ParallelRun()
{   
  string message;
  MPI_Status status;
      

  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  while ( true )
    {
      message = MyMPI_RecvCmd();

      if ( message.compare(0, 3, "ngs") == 0 ) 
        {
	  if (NGS_ParallelRun == NULL)
	    {
	      static int timer = NgProfiler::CreateTimer ("load shared library ngsolve");
	      NgProfiler::RegionTimer reg (timer);
  

	      void * handle = dlopen ("libngsolve.so", RTLD_NOW | RTLD_GLOBAL);
	      if (!handle)
		{
		  cerr << "cannot load shared library libngsolve.so" << endl;
		  exit(1);
		}
	      
	      NGS_ParallelRun = (void (*) (const string & message))  dlsym (handle, "NGS_ParallelRun");
	      
	      if (!NGS_ParallelRun)
		{
		  cerr << "cannot bind function NGS_ParallelRun" << endl;
		  exit(1);
		}
	    }
          (*NGS_ParallelRun) (message);
        }
      else if ( message == "mesh" )
	{
	  VT_USER_START ("Mesh::ReceiveParallelMesh");
	  mesh.Reset( new netgen::Mesh);
	  mesh->SendRecvMesh();
	  VT_USER_END ("Mesh::ReceiveParallelMesh");
	}

      else if ( message == "visualize" )
	{
	  cout << "parallel message visualize depreciated" << endl;
	}
      
      else if ( message == "bcastparthread" )
	{
	  MyMPI_Bcast (mparam.parthread);
	}

      else if ( message ==  "end" )
	{
	  break;
	}
      
      else
	{
	  PrintMessage ( 1, "received unidentified message '" + message + "'\n");
	  break;
	}
      
    }
}

#endif
