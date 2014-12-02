#ifndef FILE_PARALLELTOP
#define FILE_PARALLELTOP

namespace netgen
{


  class ParallelMeshTopology
  {
    const Mesh & mesh;

    /**
       mapping from local to distant vertex number
       each row of the table corresponds to one vertex
       each row contains a list of pairs (procnr, dist_vnum)
    */

    TABLE<int> loc2distvert, loc2distedge, loc2distface;

    Array<int> glob_vert, glob_edge, glob_face;
    Array<int> glob_el, glob_surfel, glob_segm;

    bool is_updated;

  public:

    ParallelMeshTopology (const Mesh & amesh);
    ~ParallelMeshTopology ();

    void Reset ();
    void Print() const;

    void UpdateCoarseGrid();
    void UpdateCoarseGridGlobal();
    // bool DoCoarseUpdate() const { return !coarseupdate; }



    /// set number of local vertices, reset sizes of loc2dist_vert, isexchangevert...
    void SetNV (int anv);
    void SetNE (int ane);
    void SetNSE (int anse);
    void SetNSegm (int anseg);


    void SetLoc2Glob_Vert   (int locnum, int globnum) { glob_vert[locnum-1] = globnum; }
    void SetLoc2Glob_Edge   (int locnum, int globnum) { glob_edge[locnum-1] = globnum; }
    void SetLoc2Glob_Face   (int locnum, int globnum) { glob_face[locnum-1] = globnum; }
    void SetLoc2Glob_VolEl  (int locnum, int globnum) { glob_el[locnum-1] = globnum; }
    void SetLoc2Glob_SurfEl (int locnum, int globnum) { glob_surfel[locnum-1] = globnum; }
    void SetLoc2Glob_Segm   (int locnum, int globnum) { glob_segm[locnum-1] = globnum; }

    int GetGlobalPNum    (int locnum) const { return glob_vert[locnum-1]; }
    int GetGlobalEdgeNum (int locnum) const { return glob_edge[locnum-1]; }
    int GetGlobalFaceNum (int locnum) const { return glob_face[locnum-1]; }
    int GetGlobalElNum   (int locnum) const { return glob_el[locnum-1]; }
    int GetGlobalSElNum  (int locnum) const { return glob_surfel[locnum-1]; }


    void SetDistantFaceNum (int dest, int locnum);
    void SetDistantPNum    (int dest, int locnum);
    void SetDistantEdgeNum (int dest, int locnum);

    int GetNDistantPNums (int locpnum) const { return loc2distvert[locpnum-1].Size(); } 
    int GetNDistantFaceNums (int locfacenum) const { return loc2distface[locfacenum-1].Size(); } 
    int GetNDistantEdgeNums ( int locedgenum) const { return loc2distedge[locedgenum-1].Size(); }
    
    void GetDistantPNums (int locpnum, int * distpnums ) const
    {
      for (int i = 0; i < loc2distvert[locpnum-1].Size(); i++ )
	distpnums[i] = loc2distvert[locpnum-1][i];
    } 
    
    void GetDistantFaceNums (int locfacenum, int * distfacenums ) const
    {
      for ( int i = 0; i < loc2distface[locfacenum-1].Size(); i++ )
	distfacenums[i] = loc2distface[locfacenum-1][i];
    } 

    void GetDistantFaceNums (int locfacenum, Array<int> & distfacenums ) const
    {
      distfacenums = loc2distface[locfacenum-1];
    }
    
    void GetDistantEdgeNums (int locedgenum, int * distedgenums ) const
    {
      for (int i = 0; i < loc2distedge[locedgenum-1].Size(); i++ )
	distedgenums[i] = loc2distedge[locedgenum-1][i];
    } 

    void GetDistantEdgeNums (int locedgenum, Array<int> & distedgenums ) const
    {
      distedgenums = loc2distedge[locedgenum-1];
    } 

    FlatArray<int> GetDistantPNums (int locnum) const { return loc2distvert[locnum]; }
    FlatArray<int> GetDistantFaceNums (int locnum) const { return loc2distface[locnum]; }
    FlatArray<int> GetDistantEdgeNums (int locnum) const { return loc2distedge[locnum]; }

    bool IsExchangeVert (int dest, int vnum) const
    {
      return loc2distvert[vnum-1].Contains (dest);
    }
  };
 

}




#endif
