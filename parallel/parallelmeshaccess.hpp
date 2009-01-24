#ifndef FILE_PARALLELTOP
#define FILE_PARALLELTOP

#ifdef PARALLEL


class ParallelMeshAccess
{
private :

  const MeshAccess & ma;

public : 

  ParallelMeshAccess ( const MeshAccess & ama );
  ~ParallelMeshAccess ()
  {
    ;
  }

//   int Glob2Loc_SurfEl ( int globnum );
//   int Glob2Loc_VolEl ( int globnum );
//   int Glob2Loc_Segm ( int globcnum );
//   int Glob2Loc_Vert ( int globnum );

  int GetLoc2Glob_VolEl ( int locnum );

  // returns number of distant processes which share this node
  int GetDistantNodeNums ( NODE_TYPE nt, int locnum, Array<int[2]> & distnums ) const;
  // ist schöner als proc und dofnr in einem flachen intarray

  int GetDistantNodeNum ( int proc, NODE_TYPE nt, int locnum) const;
  int GetDistantPNum ( const int & proc, const int & locpnum ) const;
  int GetDistantEdgeNum ( const int & proc, const int & locedgenum ) const;
  int GetDistantFaceNum ( const int & proc, const int & locfacenum ) const;
  int GetDistantElNum ( const int & proc, const int & locelnum ) const;

  bool IsExchangeNode ( NODE_TYPE nt, int nr ) const;
  bool IsExchangeFace ( const int fnr ) const;
  bool IsExchangeVert ( const int vnum ) const;
  bool IsExchangeEdge ( const int ednum ) const;
  bool IsExchangeElement ( const int elnum ) const;


  void PrintParallelMeshTopology () const;

  bool IsElementInPartition ( const int elnum, const int dest ) const;

  bool IsGhostFace ( const int facenum ) const;

  bool IsGhostEdge ( const int edgenum ) const;
};





#endif //PARALLEL

#endif
