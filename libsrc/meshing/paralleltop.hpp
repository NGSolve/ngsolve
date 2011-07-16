#ifndef FILE_PARALLELTOP
#define FILE_PARALLELTOP

namespace netgen
{

  //  extern int ntasks;

  class ParallelMeshTopology
  {
    const Mesh & mesh;

    // number of local elements, vertices, points (?), edges, faces
    int ne, nv, np, ned, nfa;

    // number of local segments and surface elements
    int nseg, nsurfel;

    // number of global elements, vertices, ???,  faces
    int neglob, nseglob, nvglob;
    int nparel;
    int nfaglob;
    int nedglob;

    /**
       mapping from local to distant vertex number
       each row of the table corresponds to one vertex
       each row contains a list of pairs (procnr, dist_vnum)
    */
    TABLE<int,PointIndex::BASE> loc2distvert;
    TABLE<int,0> loc2distedge, loc2distface;
    TABLE<int,0> loc2distel, loc2distsegm, loc2distsurfel;

    bool coarseupdate;

  public:

    ParallelMeshTopology (const Mesh & amesh);
    ~ParallelMeshTopology ();

    /// set number of local vertices, reset sizes of loc2dist_vert, isexchangevert...
    void SetNV (int anv);
    void SetNE (int ane);
    void SetNSE (int anse);
    void SetNSegm (int anseg);

    void SetNVGlob ( int anvglob )   { nvglob = anvglob; }
    void SetNEGlob ( int aneglob )   { neglob = aneglob; }
    void SetNSEGlob ( int anseglob )   { nseglob = anseglob; }

    int GetNVGlob ()  { return nvglob; }
    int GetNEGlob ()  { return neglob; }


    void Reset ();

    void SetLoc2Glob_Vert   ( int locnum, int globnum ) { loc2distvert[locnum][0] = globnum; }
    void SetLoc2Glob_VolEl  ( int locnum, int globnum ) { loc2distel[locnum-1][0] = globnum; }
    void SetLoc2Glob_SurfEl ( int locnum, int globnum ) { loc2distsurfel[locnum-1][0] = globnum; }
    void SetLoc2Glob_Segm   ( int locnum, int globnum ) { loc2distsegm[locnum-1][0] = globnum; }

    int GetLoc2Glob_Vert  ( int locnum ) const { return loc2distvert[locnum][0]; }
    int GetLoc2Glob_VolEl ( int locnum ) const { return loc2distel[locnum-1][0]; }
    int GetLoc2Glob_SurfEl ( int locnum ) const { return loc2distsurfel[locnum-1][0]; }

    void GetVertNeighbours ( int vnum, Array<int> & dests ) const;


    int GetNDistantPNums ( int locpnum ) const       
    { return loc2distvert[locpnum].Size() / 2 + 1; } 

    int GetNDistantFaceNums ( int locfacenum ) const 
    { return loc2distface[locfacenum-1].Size() / 2 + 1; } 

    int GetNDistantEdgeNums ( int locedgenum ) const  
    { return loc2distedge[locedgenum-1].Size() / 2 + 1; }

    int GetNDistantElNums ( int locelnum ) const      
    { return loc2distel[locelnum-1].Size() / 2 + 1; }

    int GetDistantPNum ( int proc, int locpnum ) const;
    int GetDistantEdgeNum ( int proc, int locedgenum ) const;
    int GetDistantFaceNum ( int proc, int locedgenum ) const;
    int GetDistantElNum ( int proc, int locelnum ) const;

    int GetDistantPNums ( int locpnum, int * distpnums ) const;
    int GetDistantEdgeNums ( int locedgenum, int * distedgenums ) const;
    int GetDistantFaceNums ( int locedgenum, int * distfacenums ) const;
    int GetDistantElNums ( int locelnum, int * distfacenums ) const;

    void Print() const;


    bool IsExchangeVert ( PointIndex vnum ) const  { return loc2distvert[vnum].Size() > 1; }
    bool IsExchangeEdge ( int ednum ) const  { return loc2distedge[ednum-1].Size() > 1; }
    bool IsExchangeFace ( int fnum ) const   { return loc2distface[fnum-1].Size() > 1; }
    bool IsExchangeElement ( int elnum ) const   { return false; }


    bool IsExchangeSEl ( int selnum ) const { return loc2distsurfel[selnum-1].Size() > 1; }


    bool IsExchangeVert (int dest, int vnum ) const
    {
      FlatArray<int> exchange = loc2distvert[vnum];
      for (int i = 1; i < exchange.Size(); i += 2)
	if (exchange[i] == dest) return true;
      return false;
    }

    bool IsExchangeEdge (int dest, int ednum ) const
    {
      FlatArray<int> exchange = loc2distedge[ednum-1];
      for (int i = 1; i < exchange.Size(); i += 2)
	if (exchange[i] == dest) return true;
      return false;
    }

    bool IsExchangeFace (int dest, int fnum ) const
    {
      FlatArray<int> exchange = loc2distface[fnum-1];
      for (int i = 1; i < exchange.Size(); i += 2)
	if (exchange[i] == dest) return true;
      return false;
    }

    bool IsExchangeElement (int dest, int elnum ) const  { return false; }

    void Update();

    void UpdateCoarseGrid();
    void UpdateRefinement ();
    void UpdateTopology ();
    void UpdateExchangeElements();

    void UpdateCoarseGridGlobal();

    bool DoCoarseUpdate() const { return !coarseupdate; }

    void SetDistantFaceNum ( int dest, int locnum, int distnum );
    void SetDistantPNum ( int dest, int locnum, int distnum );
    void SetDistantEdgeNum ( int dest, int locnum, int distnum );
    void SetDistantEl ( int dest, int locnum, int distnum );
    void SetDistantSurfEl ( int dest, int locnum, int distnum );
    void SetDistantSegm ( int dest, int locnum, int distnum );

    bool IsGhostEl ( int elnum ) const   { return mesh.VolumeElement(elnum).IsGhost(); }
  };
 

}




#endif
