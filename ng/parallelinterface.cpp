#ifdef PARALLEL
#include <mystdlib.h>

#include <meshing.hpp>

#include "../libsrc/include/parallelinterface.hpp"

namespace netgen
{
  extern AutoPtr<Mesh> mesh;
  
  
  using namespace netgen;
  
  
//   int NgPar_Glob2Loc_SurfEl ( int globnum ) 
//   { 
//     return mesh->GetParallelTopology().Glob2Loc_SurfEl(globnum+1) - 1; 
//   }
  
//   int NgPar_Glob2Loc_VolEl ( int globnum ) 
//   { 
//     return mesh->GetParallelTopology().Glob2Loc_VolEl(globnum+1) - 1; 
//   }
  
//   int NgPar_Glob2Loc_Segm (   int globnum )   
//   { 
//     return mesh->GetParallelTopology().Glob2Loc_Segm(globnum+1) - 1; 
//   }
  
//   int NgPar_Glob2Loc_Vert ( int globnum )  
//   { 
//     return mesh->GetParallelTopology().Glob2Loc_Vert(globnum+1) -1;
//   }
  
  
  int NgPar_GetLoc2Glob_VolEl ( int locnum )
  { 
    return mesh -> GetParallelTopology().GetLoc2Glob_VolEl ( locnum+1) -1; 
  }

  // gibt anzahl an distant pnums zurueck
  // * pnums entspricht ARRAY<int[2] >
  int NgPar_GetDistantNodeNums ( int nodetype, int locnum, int * distnums )
  {
    int size;
    switch ( nodetype )
      {
      case 0:
	size = mesh->GetParallelTopology().GetDistantPNums( locnum+1, distnums ); 
	break;
      case 1:
	size = mesh->GetParallelTopology().GetDistantEdgeNums( locnum+1, distnums ); 
	break;
      case 2:
	size = mesh->GetParallelTopology().GetDistantFaceNums( locnum+1, distnums );
	break;
      case 3:
	size = mesh->GetParallelTopology().GetDistantElNums( locnum+1, distnums );
	break;
      default:
	cerr << "NgPar_GetDistantNodeNums() Unknown nodetype " << nodetype << endl;
	size = -1;
      }
    // 0 - based 
    for ( int i = 0; i < size; i++ )
      distnums[2*i+1]--;

    return size;
  }

  int NgPar_GetNDistantNodeNums ( int nodetype, int locnum )
  {
    switch ( nodetype )
      {
      case 0:
	return mesh->GetParallelTopology().GetNDistantPNums( locnum+1 );
      case 1:
	return mesh->GetParallelTopology().GetNDistantEdgeNums( locnum+1 );
      case 2:
	return mesh->GetParallelTopology().GetNDistantFaceNums( locnum+1 );
      case 3:
	return mesh->GetParallelTopology().GetNDistantElNums( locnum+1 );
      }
    return -1;
  }

  int NgPar_GetDistantPNum ( int proc, int locpnum )  
  { 
    return mesh->GetParallelTopology().GetDistantPNum( proc, locpnum+1) - 1; 
  }

  int NgPar_GetDistantEdgeNum ( int proc, int locpnum )  
  { 
    return mesh->GetParallelTopology().GetDistantEdgeNum( proc, locpnum+1) - 1; 
  }

  int NgPar_GetDistantFaceNum ( int proc, int locpnum )  
  { 
    return mesh->GetParallelTopology().GetDistantFaceNum (proc, locpnum+1 ) - 1; 
  }

  int NgPar_GetDistantElNum ( int proc, int locelnum )
  { 
    return mesh->GetParallelTopology().GetDistantElNum (proc, locelnum+1 ) - 1; 
  }

  bool NgPar_IsExchangeFace ( int fnr ) 
  { 
    return (mesh->GetParallelTopology().GetNDistantFaceNums( fnr+1 ) > 0);
    // return mesh->GetParallelTopology().IsExchangeFace ( fnr+1 ); 
  }

  bool NgPar_IsExchangeVert ( int vnum )
  { 
    return (mesh->GetParallelTopology().GetNDistantPNums( vnum+1 ) > 0);
    // return mesh->GetParallelTopology().IsExchangeVert ( vnum+1 ); 
  }

  bool NgPar_IsExchangeEdge ( int ednum )  
  { 
    return (mesh->GetParallelTopology().GetNDistantEdgeNums( ednum+1 ) > 0);
    // return mesh->GetParallelTopology().IsExchangeEdge ( ednum+1 ); 
  }

  bool NgPar_IsExchangeElement ( int elnum )  
  { 
    return (mesh->GetParallelTopology().GetNDistantElNums( elnum+1 ) > 0);
    // return mesh->GetParallelTopology().IsExchangeElement ( elnum+1 ); 
  }


  void NgPar_PrintParallelMeshTopology ()
  { 
    mesh -> GetParallelTopology().Print (); 
  }
 
  bool NgPar_IsElementInPartition ( const int elnum, const int dest )
  {
    return mesh -> GetParallelTopology().IsElementInPartition ( elnum+1, dest ); 
  }



  bool NgPar_IsGhostFace ( const int facenum )
  { 
    return mesh -> GetParallelTopology().IsGhostFace ( facenum+1); 
  }

  bool NgPar_IsGhostEdge ( const int edgenum )
  { 
    return mesh -> GetParallelTopology().IsGhostEdge ( edgenum+1); 
  }

}

#endif
