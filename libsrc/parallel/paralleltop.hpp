#ifndef FILE_PARALLELTOP
#define FILE_PARALLELTOP

#include <meshing.hpp>

extern int ntasks;

class ParallelMeshTopology
{
  const Mesh & mesh;

  // number of local elements, vertices, points (?), edges, faces
  int ne, nv, np, ned, nfa;

  // number of local segments and surface elements
  int nseg, nsurfel;

  // number of global elements, vertices, ???,  faces
  int neglob, nvglob;
  int nparel;
  int nfaglob;

  /**
     mapping from local to distant vertex number
     each row of the table corresponds to one vertex
     each row contains a list of pairs (procnr, dist_vnum)
   */
  TABLE<int,PointIndex::BASE> loc2distvert;
  TABLE<int,0> loc2distedge, loc2distface, loc2distel;
  TABLE<int,0> loc2distsegm, loc2distsurfel;

  BitArray * isexchangeface, * isexchangevert, * isexchangeedge, * isexchangeel;

  BitArray isghostedge, isghostface;

  bool coarseupdate;
  int overlap;

public:

  ParallelMeshTopology (const Mesh & amesh);
  ~ParallelMeshTopology ();

  /// set number of local vertices, reset sizes of loc2dist_vert, isexchangevert...
  void SetNV ( int anv );
  void SetNE ( int ane );
  void SetNSE ( int anse );
  void SetNSegm ( int anseg );

  void Reset ();

  void SetLoc2Glob_Vert   ( int locnum, int globnum ) { loc2distvert[locnum][0] = globnum; }
  void SetLoc2Glob_VolEl  ( int locnum, int globnum ) { loc2distel[locnum-1][0] = globnum; }
  void SetLoc2Glob_SurfEl ( int locnum, int globnum ) { loc2distsurfel[locnum-1][0] = globnum; }
  void SetLoc2Glob_Segm   ( int locnum, int globnum ) { loc2distsegm[locnum-1][0] = globnum; }

  int GetLoc2Glob_Vert  ( int locnum ) const { return loc2distvert[locnum][0]; }
  int GetLoc2Glob_VolEl ( int locnum ) const { return loc2distel[locnum-1][0]; }

  void GetVertNeighbours ( int vnum, Array<int> & dests ) const;

  int Glob2Loc_SurfEl ( int globnum );
  int Glob2Loc_VolEl ( int globnum );
  int Glob2Loc_Segm ( int globnum );
  int Glob2Loc_Vert ( int globnum );

  int GetNDistantPNums ( int locpnum ) const        { return loc2distvert[locpnum].Size() / 2 + 1; } 
  int GetNDistantFaceNums ( int locfacenum ) const  { return loc2distface[locfacenum-1].Size() / 2 + 1; } 
  int GetNDistantEdgeNums ( int locedgenum ) const  { return loc2distedge[locedgenum-1].Size() / 2 + 1; }
  int GetNDistantElNums ( int locelnum ) const      { return loc2distel[locelnum-1].Size() / 2 + 1; }

  int GetDistantPNum ( int proc, int locpnum ) const;
  int GetDistantEdgeNum ( int proc, int locedgenum ) const;
  int GetDistantFaceNum ( int proc, int locedgenum ) const;
  int GetDistantElNum ( int proc, int locelnum ) const;

  int GetDistantPNums ( int locpnum, int * distpnums ) const;
  int GetDistantEdgeNums ( int locedgenum, int * distedgenums ) const;
  int GetDistantFaceNums ( int locedgenum, int * distfacenums ) const;
  int GetDistantElNums ( int locelnum, int * distfacenums ) const;

  void Print() const;

  void SetExchangeFace ( int fnr )      { isexchangeface->Set( (fnr-1)*(ntasks+1) ); }
  void SetExchangeVert ( int vnum )     { isexchangevert->Set( (vnum-1)*(ntasks+1) ); }
  void SetExchangeElement ( int elnum ) { isexchangeel->Set((elnum-1)*(ntasks+1) ); }
  void SetExchangeEdge ( int ednum )    { isexchangeedge->Set ((ednum-1)*(ntasks+1)); }


  // fuer diese Fkts  kann man ja O(N) - bitarrays lassen
  bool IsExchangeVert ( PointIndex vnum ) const
  {
    return loc2distvert[vnum].Size() > 1;
    // return isexchangevert->Test((vnum-1)*(ntasks+1));
  }

  bool IsExchangeEdge ( int ednum ) const
  {
    /*
    if ( isexchangeedge->Test( (ednum-1)*(ntasks+1) )  !=
	 ( loc2distedge[ednum-1].Size() > 1)  )
      {
	cerr << "isexedge differs, id = " << id << ", ednum = " << ednum << endl;
      }
    */
    // return loc2distedge[ednum-1].Size() > 1;
    int i = (ednum-1)*(ntasks+1);
    return isexchangeedge->Test(i);
  }

  bool IsExchangeFace ( int fnum ) const
  {
    /*
    if ( isexchangeface->Test( (fnum-1)*(ntasks+1) )  !=
	 ( loc2distface[fnum-1].Size() > 1)  )
      {
	cerr << "it differs, id = " << id << ", fnum = " << fnum << endl;
      }
    */
    // nur hier funktioniert's so nicht ???? (JS)
    // return loc2distface[fnum-1].Size() > 1;
    return isexchangeface->Test( (fnum-1)*(ntasks+1) );
  }

  bool IsExchangeElement ( int elnum ) const
  {
    // return loc2distel[elnum-1].Size() > 1;
    return isexchangeel->Test((elnum-1)*(ntasks+1));
  }

  bool IsExchangeSEl ( int selnum ) const
  {
    return loc2distsurfel[selnum-1].Size() > 1;
  }


  void  SetExchangeFace (int dest, int fnr )
  {
    // isexchangeface->Set((fnr-1)*(ntasks+1) + (dest+1));
  }

  void  SetExchangeVert (int dest, int vnum )
  {
    // isexchangevert->Set((vnum-1)*(ntasks+1) + (dest+1) );
  }

  void  SetExchangeElement (int dest, int vnum )
  {
    isexchangeel->Set( (vnum-1)*(ntasks+1) + (dest+1) );
  }

  void  SetExchangeEdge (int dest, int ednum )
  {
    // isexchangeedge->Set ( (ednum-1)*(ntasks+1) + (dest+1) );
  }


  bool IsExchangeVert (int dest, int vnum ) const
  {
    FlatArray<int> exchange = loc2distvert[vnum];
    for (int i = 1; i < exchange.Size(); i += 2)
      if (exchange[i] == dest) return true;
    return false;
    // return isexchangevert->Test((vnum-1)*(ntasks+1) + (dest+1) );
  }

  bool IsExchangeEdge (int dest, int ednum ) const
  {
    FlatArray<int> exchange = loc2distedge[ednum-1];
    for (int i = 1; i < exchange.Size(); i += 2)
      if (exchange[i] == dest) return true;
    return false;

    // return isexchangeedge->Test((ednum-1)*(ntasks+1) + (dest+1));
  }

  bool IsExchangeFace (int dest, int fnum ) const
  {
    FlatArray<int> exchange = loc2distface[fnum-1];
    for (int i = 1; i < exchange.Size(); i += 2)
      if (exchange[i] == dest) return true;
    return false;

    // return isexchangeface->Test((fnum-1)*(ntasks+1) + (dest+1) );
  }

  bool IsExchangeElement (int dest, int elnum ) const
  {
    // das geht auch nicht
    /*
    FlatArray<int> exchange = loc2distel[elnum-1];
    for (int i = 1; i < exchange.Size(); i += 2)
      if (exchange[i] == dest) return true;
    return false;
    */
    return isexchangeel->Test((elnum-1)*(ntasks+1) + (dest+1));
  }

  int Overlap() const
  { return overlap; }

  int IncreaseOverlap ()
  { overlap++; return overlap; }

  void Update();

  void UpdateCoarseGrid();
  void UpdateCoarseGridOverlap();
  void UpdateRefinement ();
  void UpdateTopology ();
  void UpdateExchangeElements();

  void UpdateCoarseGridGlobal();

  bool DoCoarseUpdate() const 
  { return !coarseupdate; }


  void SetDistantFaceNum ( int dest, int locnum, int distnum );
  void SetDistantPNum ( int dest, int locnum, int distnum );
  void SetDistantEdgeNum ( int dest, int locnum, int distnum );
  void SetDistantEl ( int dest, int locnum, int distnum );
  void SetDistantSurfEl ( int dest, int locnum, int distnum );
  void SetDistantSegm ( int dest, int locnum, int distnum );

  bool IsGhostEl ( int elnum ) const   { return mesh.VolumeElement(elnum).IsGhost(); }
  bool IsGhostFace ( int fanum ) const { return isghostface.Test(fanum-1); }
  bool IsGhostEdge ( int ednum ) const { return isghostedge.Test(ednum-1); }
  bool IsGhostVert ( int pnum ) const  { return mesh.Point(pnum).IsGhost(); }

//    inline void SetGhostVert ( const int pnum )
//    { isghostvert.Set(pnum-1); }

  void SetGhostEdge ( int ednum )
  { isghostedge.Set(ednum-1); }

  void ClearGhostEdge ( int ednum )
  { isghostedge.Clear(ednum-1); }

  void SetGhostFace ( int fanum )
  { isghostface.Set(fanum-1); }

  void ClearGhostFace ( int fanum )
  { isghostface.Clear(fanum-1); }

  bool IsElementInPartition ( int elnum, int dest ) const
  { 
    return IsExchangeElement ( dest, elnum); 
  }
  
  void SetElementInPartition ( int elnum, int dest )
  { 
    SetExchangeElement ( dest, elnum );
  }

  void SetNVGlob ( int anvglob )   { nvglob = anvglob; }
  void SetNEGlob ( int aneglob )   { neglob = aneglob; }

  int GetNVGlob ()  { return nvglob; }
  int GetNEGlob ()  { return neglob; }
};
 






#endif
