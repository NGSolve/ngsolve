#ifdef PARALLEL
 
#include <parallelngs.hpp>
// #include <parallel.hpp>
// #include "../parallel/parallelinterface.hpp"


#include <parallelinterface.hpp>


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
      case NT_VERTEX: return NgPar_GetDistantPNum ( proc, locnum ); 
      case NT_EDGE: return NgPar_GetDistantEdgeNum ( proc, locnum ); 
      case NT_FACE: return NgPar_GetDistantFaceNum ( proc, locnum ); 
      case NT_CELL: return NgPar_GetDistantElNum ( proc, locnum ); 
      }
    return -1;
  }


  int ParallelMeshAccess :: GetDistantNodeNums ( NODE_TYPE nt, int locnum, 
                                                 ngstd::Array<int[2]> & distnums ) const
  {
    distnums.SetSize( NgPar_GetNDistantNodeNums(nt, locnum) );
    int nr = NgPar_GetDistantNodeNums ( nt, locnum, &distnums[0][0] );
    return nr;
  }

  int  ParallelMeshAccess ::GetDistantPNum ( int proc, int locpnum ) const
  {
    return NgPar_GetDistantPNum ( proc, locpnum ); 
  }

  int  ParallelMeshAccess ::GetDistantEdgeNum ( int proc, int locedgenum ) const
  { return NgPar_GetDistantEdgeNum ( proc, locedgenum ); }

  int  ParallelMeshAccess ::GetDistantFaceNum ( int proc, int locfacenum ) const
  { return NgPar_GetDistantFaceNum ( proc, locfacenum ); }

  int  ParallelMeshAccess ::GetDistantElNum ( int proc, int locelnum ) const
  { return NgPar_GetDistantElNum ( proc, locelnum ); }

  bool ParallelMeshAccess :: IsExchangeNode ( NODE_TYPE nt, int nr ) const
  {
    switch (nt)
      {
      case NT_VERTEX: return NgPar_IsExchangeVert ( nr ); 
      case NT_EDGE:   return NgPar_IsExchangeEdge ( nr ); 
      case NT_FACE:   return NgPar_IsExchangeFace ( nr ); 
      case NT_CELL:   return NgPar_IsExchangeElement ( nr ); 
      }
    return 0;
  }

  bool ParallelMeshAccess ::IsExchangeFace ( int fnr ) const 
  { 
    return NgPar_IsExchangeFace ( fnr ); 
  }

  bool ParallelMeshAccess ::IsExchangeVert (int vnum ) const 
  {
    return NgPar_IsExchangeVert ( vnum ); 
  }

  bool ParallelMeshAccess ::IsExchangeEdge (int ednum ) const
  { 
    return NgPar_IsExchangeEdge ( ednum ); 
  }

  bool ParallelMeshAccess ::IsExchangeElement (int elnum ) const
  { 
    return NgPar_IsExchangeElement ( elnum ); 
  }

  void ParallelMeshAccess :: PrintParallelMeshTopology () const
  { 
    NgPar_PrintParallelMeshTopology (); 
  }




  bool ParallelMeshAccess :: IsElementInPartition ( int elnum, int dest ) const
  { 
    return NgPar_IsElementInPartition ( elnum, dest ); 
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
