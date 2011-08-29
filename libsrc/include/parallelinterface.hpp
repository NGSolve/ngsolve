gibt's nicht mehr


#ifndef FILE_PARALLELINTERFACE
#define FILE_PARALLELINTERFACE

#ifdef PARALLEL

#ifdef __cplusplus
extern "C" {
#endif

  // this interface is 0-base  !!


  int NgPar_GetLoc2Glob_VolEl ( int locnum );

  // int NgPar_GetDistantNodeNums ( int nt, int locnum, int * procs, int * distnum);

  // number on distant processor 

  // gibt anzahl an distant pnums zurueck
  // * pnums entspricht ARRAY<int[2] >
  int NgPar_GetDistantNodeNums ( int nodetype, int locnum, int * pnums );
  int NgPar_GetNDistantNodeNums ( int nodetype, int locnum );

  int NgPar_GetDistantPNum ( int proc, int locnum ) ;
  int NgPar_GetDistantEdgeNum ( int proc, int locnum ) ;
  int NgPar_GetDistantFaceNum ( int proc, int locnum ) ;
  int NgPar_GetDistantElNum ( int proc, int locnum );

  bool NgPar_IsExchangeFace ( int fnr ) ;
  bool NgPar_IsExchangeVert ( int vnum );
  bool NgPar_IsExchangeEdge ( int ednum );
  bool NgPar_IsExchangeElement ( int elnum );

  void NgPar_PrintParallelMeshTopology ();
  bool NgPar_IsElementInPartition ( int elnum, int dest );

  bool NgPar_IsGhostFace ( int facenum );
  bool NgPar_IsGhostEdge ( int edgenum );

#ifdef __cplusplus
}
#endif

#endif

#endif
