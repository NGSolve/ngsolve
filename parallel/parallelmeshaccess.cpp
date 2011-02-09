#ifdef PARALLEL

#include <parallelngs.hpp>
// #include <parallel.hpp>
// #include "../parallel/parallelinterface.hpp"

namespace netgen {
#include <parallelinterface.hpp>
}
namespace ngparallel
{

  using namespace ngparallel;
  using namespace netgen;

  ParallelMeshAccess :: ParallelMeshAccess ( const MeshAccess & ama )
    : ma( ama )
  {
   ;
  }

  int ParallelMeshAccess ::GetDistantNodeNum (int proc, NODE_TYPE nt, int locnum) const
  {
    switch (nt)
      {
      case NT_VERTEX: return netgen::NgPar_GetDistantPNum ( proc, locnum ); 
      case NT_EDGE: return netgen::NgPar_GetDistantEdgeNum ( proc, locnum ); 
      case NT_FACE: return netgen::NgPar_GetDistantFaceNum ( proc, locnum ); 
      case NT_CELL: return netgen::NgPar_GetDistantElNum ( proc, locnum ); 
      }
    return -1;
  }


  int ParallelMeshAccess :: GetDistantNodeNums ( NODE_TYPE nt, int locnum, 
                                                 ngstd::Array<int[2]> & distnums ) const
  {
    distnums.SetSize( netgen::NgPar_GetNDistantNodeNums(nt, locnum) );
    int nr = netgen::NgPar_GetDistantNodeNums ( nt, locnum, &distnums[0][0] );
    return nr;
  }

  int  ParallelMeshAccess ::GetDistantPNum ( int proc, int locpnum ) const
  {
    return netgen::NgPar_GetDistantPNum ( proc, locpnum ); 
  }

  int  ParallelMeshAccess ::GetDistantEdgeNum ( int proc, int locedgenum ) const
  { return netgen::NgPar_GetDistantEdgeNum ( proc, locedgenum ); }

  int  ParallelMeshAccess ::GetDistantFaceNum ( int proc, int locfacenum ) const
  { return netgen::NgPar_GetDistantFaceNum ( proc, locfacenum ); }

  int  ParallelMeshAccess ::GetDistantElNum ( int proc, int locelnum ) const
  { return netgen::NgPar_GetDistantElNum ( proc, locelnum ); }

  bool ParallelMeshAccess :: IsExchangeNode ( NODE_TYPE nt, int nr ) const
  {
    switch (nt)
      {
      case NT_VERTEX: return netgen::NgPar_IsExchangeVert ( nr ); 
      case NT_EDGE:   return netgen::NgPar_IsExchangeEdge ( nr ); 
      case NT_FACE:   return netgen::NgPar_IsExchangeFace ( nr ); 
      case NT_CELL:   return netgen::NgPar_IsExchangeElement ( nr ); 
      }
    return 0;
  }

  bool ParallelMeshAccess ::IsExchangeFace ( int fnr ) const 
  { 
    return netgen::NgPar_IsExchangeFace ( fnr ); 
  }

  bool ParallelMeshAccess ::IsExchangeVert (int vnum ) const 
  {
    return netgen::NgPar_IsExchangeVert ( vnum ); 
  }

  bool ParallelMeshAccess ::IsExchangeEdge (int ednum ) const
  { 
    return netgen::NgPar_IsExchangeEdge ( ednum ); 
  }

  bool ParallelMeshAccess ::IsExchangeElement (int elnum ) const
  { 
    return netgen::NgPar_IsExchangeElement ( elnum ); 
  }

  void ParallelMeshAccess :: PrintParallelMeshTopology () const
  { 
    netgen::NgPar_PrintParallelMeshTopology (); 
  }




  bool ParallelMeshAccess :: IsElementInPartition ( int elnum, int dest ) const
  { 
    return netgen::NgPar_IsElementInPartition ( elnum, dest ); 
  }

  int ParallelMeshAccess :: GetLoc2Glob_VolEl ( int locnum )
  { 
    return NgPar_GetLoc2Glob_VolEl ( locnum ); 
  }



  bool ParallelMeshAccess :: IsGhostFace ( int facenum ) const
  { 
    return NgPar_IsGhostFace ( facenum );
  }

  bool ParallelMeshAccess :: IsGhostEdge ( int edgenum ) const
  {
    return NgPar_IsGhostEdge ( edgenum ); 
  }
}


#endif
