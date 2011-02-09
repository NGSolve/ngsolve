#ifndef FILE_PARALLELTOP
#define FILE_PARALLELTOP

#ifdef PARALLEL


class ParallelMeshAccess
{
private :

  const MeshAccess & ma;

public : 

  ParallelMeshAccess (const MeshAccess & ama);
  ~ParallelMeshAccess () { ; }


  int GetLoc2Glob_VolEl (int locnum);

  // returns number of distant processes which share this node
  int GetDistantNodeNums (NODE_TYPE nt, int locnum, Array<int[2]> & distnums) const;

  int GetDistantNodeNum (int proc, NODE_TYPE nt, int locnum) const;
  int GetDistantPNum (int proc, int locpnum) const;
  int GetDistantEdgeNum (int proc, int locedgenum) const;
  int GetDistantFaceNum (int proc, int locfacenum) const;
  int GetDistantElNum (int proc, int locelnum) const;

  bool IsExchangeNode (NODE_TYPE nt, int nr) const;
  bool IsExchangeFace (int fnr) const;
  bool IsExchangeVert (int vnum) const;
  bool IsExchangeEdge (int ednum) const;
  bool IsExchangeElement (int elnum) const;


  void PrintParallelMeshTopology () const;
  bool IsElementInPartition (int elnum, int dest) const;
  bool IsGhostFace (int facenum) const;
  bool IsGhostEdge (int edgenum) const;
};





#endif //PARALLEL

#endif
