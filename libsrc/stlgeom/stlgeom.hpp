#ifndef FILE_STLGEOM
#define FILE_STLGEOM

/**************************************************************************/
/* File:   stlgeom.hpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Author2: Johannes Gerstmayr                                            */
/* Date:   26. Jul. 99                                                    */
/**************************************************************************/

/**
   STL Geometry


   Terminology:
   
   Point ... coordinates of STL triangles
   Triangle  (short Trig)  STL triangle
   TopEdge .... edge in topology, boundary of STL triangles (many)
   Edge .... Edges which will occur in the mesh (confirmed edges, less)
*/


#include <meshing.hpp>


namespace netgen
{
  inline int IsInArray(int n, const Array<int>& ia)
  {
    return ia.Contains(n); 
  }

  inline bool AddIfNotExists(Array<int>& list, int x)
  {
    if (list.Contains(x)) return false;
    list.Append(x);
    return true;
  }
  
  extern DLL_HEADER MeshingParameters mparam;
  


#include "stltopology.hpp"
#include "stltool.hpp"
#include "stlline.hpp"
 






  class STLEdgeDataList
  {
    Array<int> storedstatus;
    STLTopology & geom;
  public:
  
    STLEdgeDataList(STLTopology & ageom);
    ~STLEdgeDataList();

    void Store ();
    void Restore ();

    void SetSize(int /* size */) { };
    void Clear() { };
    int Size() const { return geom.GetNTE(); }
    const STLTopEdge & Get(int i) const { return geom.GetTopEdge(i); }
    STLTopEdge & Elem(int i) { return geom.GetTopEdge(i); }

    int GetNEPP(int pn) const {return geom.NTopEdgesPerPoint(pn); }
    int GetEdgePP(int pn, int vi) const {return geom.TopEdgePerPoint(pn, vi);};

    //void AddEdgePP(int pn, int vn) { } ;

    void ResetAll();
    void ChangeStatus(int status1, int status2);

    int GetEdgeNum(int np1, int np2) const
    { return geom.GetTopEdgeNum (np1, np2); }

    int GetNConfEdges() const;

    void Write(ofstream& of) const;
    void Read(ifstream& ifs);

    void BuildLineWithEdge(int ep1, int ep2, Array<twoint>& line);
    void BuildClusterWithEdge(int ep1, int ep2, Array<twoint>& line);

    int GetNEPPStat(int p, int status) const;
    int GetNConfCandEPP(int p) const;
  };






  class STLGeometry : public STLTopology, public NetgenGeometry
  {
    // edges to be meshed:
    Array<STLEdge> edges;
    //edges per point
    TABLE<int> edgesperpoint;

    // line: a connection of edges
    Array<STLLine*> lines;
    Array<int> lineendpoints; //per geometrypoint, 1 = is endpoint; 0 = no endpoint,

    Array<Vec3d> normals; //normals belong to points!

    Array<twoint> externaledges;

    int undoexternaledges;
    Array<twoint> storedexternaledges;

    STLEdgeDataList * edgedata;
    //  STLEdgeDataList edgedata_store;
    int calcedgedataanglesnew;

    int edgedatastored;



    int facecnt; 
    //meshpoint is only set, if an edge is at this point!!!

    Array<int> vicinity; //is one, if a triangle belongs to vicinity (eg. of selecttrig)
    Array<int> markedtrigs; //is one, if a triangle belongs to marked triangles (calcdirtystrigs)
    Array<Point3d> markedsegs; //every pointpair is a segment!!!  
    Array<twoint> selectedmultiedge;


    //spiralpoints:
    Array<int> spiralpoints;
    //
    Array<STLChart*> atlas;
    //marks all already charted trigs with chartnumber
    Array<int> chartmark; 
    //outerchartspertrig, ascending sorted
    TABLE<int> outerchartspertrig;


    //for meshing and project:
    Array<int> meshcharttrigs; //per trig: 1=belong to chart, 0 not
    int meshchart;

    Array<int> ha_points;  // help array, np long, filled with 0 


    // sharp geometric edges not declared as edges
    // (not considered for spiral check)
    INDEX_2_HASHTABLE<int> * smoothedges;


    //transformation:
    Vec<3> meshtrignv;
    Vec<3> ex, ey, ez;
    Point<3> p1;

    mutable class RefinementSTLGeometry * ref; 
    
  public:
    int edgesfound;
    int surfacemeshed;
    int surfaceoptimized;
    int volumemeshed;

    int trigsconverted; //when STLTriangles exist -> 1

    //for selecting nodes
    //int selecttrig, nodeofseltrig;

    //only for testing;
    Array<STLLine*> meshlines;
    Array<Point3d> meshpoints;

    double area;
  public:
    STLGeometry();
    virtual ~STLGeometry();


    void Clear();

    virtual void Save (string filename) const;


	DLL_HEADER void STLInfo(double* data);
    //stldoctor:
	DLL_HEADER void SmoothNormals();
	DLL_HEADER void MarkNonSmoothNormals();

	DLL_HEADER void CalcEdgeData();
	DLL_HEADER void CalcEdgeDataAngles();

    const STLEdgeDataList& EdgeDataList() const {return *edgedata;}

	DLL_HEADER void UndoEdgeChange();
	DLL_HEADER void StoreEdgeData();
	DLL_HEADER void RestoreEdgeData();

    //void ClearSelectedMultiEdge() {selectedmultiedge.SetSize(0);}
    //void AddSelectedMultiEdge(twoint ep) {selectedmultiedge.Append(ep);}
    //int SelectedMultiEdgeSize() {return selectedmultiedge.Size();}
    const Array<twoint>& SelectedMultiEdge() {return selectedmultiedge;}
    twoint GetNearestSelectedDefinedEdge();
    void BuildSelectedMultiEdge(twoint ep);
    void BuildSelectedEdge(twoint ep);
    void BuildSelectedCluster(twoint ep);

	DLL_HEADER void ImportEdges();
	DLL_HEADER void AddEdges(const Array<Point<3> >& eps);
	DLL_HEADER void ExportEdges();
	DLL_HEADER void LoadEdgeData(const char* file);
	DLL_HEADER void SaveEdgeData(const char* file);
    //  void SetEdgeAtSelected(int mode);
  

	DLL_HEADER void STLDoctorConfirmEdge();
	DLL_HEADER void STLDoctorCandidateEdge();
	DLL_HEADER void STLDoctorExcludeEdge();
	DLL_HEADER void STLDoctorUndefinedEdge();

	DLL_HEADER void STLDoctorSetAllUndefinedEdges();
	DLL_HEADER void STLDoctorEraseCandidateEdges();
	DLL_HEADER void STLDoctorConfirmCandidateEdges();
	DLL_HEADER void STLDoctorConfirmedToCandidateEdges();

	DLL_HEADER void STLDoctorDirtyEdgesToCandidates();
	DLL_HEADER void STLDoctorLongLinesToCandidates();

	DLL_HEADER void UndoExternalEdges();
	DLL_HEADER void StoreExternalEdges();
	DLL_HEADER void RestoreExternalEdges();

	DLL_HEADER void ImportExternalEdges(const char * filename);  // Flame edges, JS
    //  void LoadExternalEdges();

	DLL_HEADER void BuildExternalEdgesFromEdges();
	DLL_HEADER void SaveExternalEdges();
	DLL_HEADER void AddExternalEdgeAtSelected();
	DLL_HEADER void AddClosedLinesToExternalEdges();
	DLL_HEADER void AddLongLinesToExternalEdges();
	DLL_HEADER void AddAllNotSingleLinesToExternalEdges();
	DLL_HEADER void STLDoctorBuildEdges();
	DLL_HEADER void AddExternalEdgesFromGeomLine();
	DLL_HEADER void DeleteDirtyExternalEdges();
	DLL_HEADER void DeleteExternalEdgeAtSelected();
	DLL_HEADER void DeleteExternalEdgeInVicinity();
    void AddExternalEdge(int p1, int p2);
    void DeleteExternalEdge(int p1, int p2);
    int IsExternalEdge(int p1, int p2);
    int NOExternalEdges() const {return externaledges.Size();}
    twoint GetExternalEdge(int i) const {return externaledges.Get(i);}

	DLL_HEADER void DestroyDirtyTrigs();
	DLL_HEADER void CalcNormalsFromGeometry();
	DLL_HEADER void MoveSelectedPointToMiddle();
	DLL_HEADER void NeighbourAnglesOfSelectedTrig();
	DLL_HEADER void PrintSelectInfo();
	DLL_HEADER void ShowSelectedTrigChartnum();
	DLL_HEADER void ShowSelectedTrigCoords();
	DLL_HEADER void SmoothGeometry ();


	DLL_HEADER void LoadMarkedTrigs();
	DLL_HEADER void SaveMarkedTrigs();
	void ClearMarkedSegs() {markedsegs.SetSize(0);}
    void AddMarkedSeg(const Point<3> & ap1, const Point<3> & ap2) 
    {
      markedsegs.Append(ap1);markedsegs.Append(ap2);
    }

    void GetMarkedSeg(int i, Point<3> & ap1, Point<3> & ap2) 
    {
      ap1=markedsegs.Get(i*2-1); 
      ap2=markedsegs.Get(i*2);
    }
    int GetNMarkedSegs() {return markedsegs.Size()/2;}
	DLL_HEADER void CalcVicinity(int starttrig);
	DLL_HEADER void GetVicinity(int starttrig, int size, Array<int>& vic);

	DLL_HEADER int Vicinity(int trig) const;

	DLL_HEADER void InitMarkedTrigs();
	DLL_HEADER void MarkDirtyTrigs();
	DLL_HEADER void SmoothDirtyTrigs();
	DLL_HEADER void GeomSmoothRevertedTrigs();
	DLL_HEADER void MarkRevertedTrigs();
	DLL_HEADER double CalcTrigBadness(int i);
	DLL_HEADER int IsMarkedTrig(int trig) const;
	DLL_HEADER void SetMarkedTrig(int trig, int num);
	DLL_HEADER void MarkTopErrorTrigs ();

    //Selected triangle
	DLL_HEADER void SetSelectTrig(int trig);
	DLL_HEADER int GetSelectTrig() const;
	DLL_HEADER void SetNodeOfSelTrig(int n);
	DLL_HEADER int GetNodeOfSelTrig() const;


    int AddNormal(const Vec3d& n) {return normals.Append(n);}
    const Vec3d & GetNormal(int nr) const {return normals.Get(nr);}
    void SetNormal(int nr, const Vec3d& n) {normals.Elem(nr) = n;}

    int AddEdge(const STLEdge& v) {return edges.Append(v);}
    int AddEdge(int p1, int p2);

    STLEdge GetEdge(int nr) {return edges.Get(nr);}
    int GetNE() {return edges.Size();}

    double Area();

    double GetAngle(int t1, int t2);
    double GetGeomAngle(int t1, int t2);
    //if triangles t1 and t2 touch, return 1 and in p1, p2 the touching points
    //int TrigsTouch(int t1, int t2, int& p1, int& p2);


  
    ///

    ///ReadTriangle->STLTriangle, initialise some important variables, always after load!!!
    virtual void InitSTLGeometry (const Array<STLReadTriangle> & readtrigs);
    virtual void TopologyChanged(); //do some things, if topology changed!
    int CheckGeometryOverlapping();

    //get NO edges per point
    int GetEPPSize() const {return edgesperpoint.Size();};
    int GetNEPP(int pn) 
    {
      if (edgesperpoint.Size() == 0) {BuildEdgesPerPoint();}
      return edgesperpoint.EntrySize(pn);
    };
    int GetEdgePP(int pn, int vi)
    {
      if (edgesperpoint.Size() == 0) {BuildEdgesPerPoint();}
      return edgesperpoint.Get(pn,vi);
    };
    void AddEdgePP(int pn, int vn) {edgesperpoint.Add1(pn,vn);};
    //von 2 punkten ermitteln, ob sie eine Kante sind
    int IsEdge(int p1, int p2);
    int IsEdgeNum(int p1, int p2);

    ///Build EdgeSegments
    void ClearEdges();
    void BuildEdges();
    void BuildEdgesPerPoint();
    void UseExternalEdges();


    void FindEdgesFromAngles();
    void CalcFaceNums();
    int GetNOBodys();
    int GetNOFaces() {return facecnt;}
    void LinkEdges();

    void AddConeAndSpiralEdges();
    void AddFaceEdges(); //each face should have at least one starting edge (outherwise it won't be meshed)

    void GetDirtyChartTrigs(int chartnum, STLChart& chart, const Array<int>& outercharttrigs, 
			    Array<int>& chartpointchecked, Array<int>& dirtytrigs);

    void ClearSpiralPoints();
    void SetSpiralPoint(int pn) {spiralpoints.Elem(pn) = 1;};
    int GetSpiralPoint(int pn) const {return spiralpoints.Get(pn);};

    void GetSortedTrianglesAroundPoint(int p, int starttrig, Array<int>& trigs);

    // smooth edges: sharp geometric edges not declared as edges
    void BuildSmoothEdges ();
    int IsSmoothEdge (int pi1, int pi2) const;


    //make charts with regions of a max. angle
    void MakeAtlas(class Mesh & mesh);

    //outerchartspertrig, sorted!
    int GetOCPTSize() const {return outerchartspertrig.Size();};
    int GetNOCPT(int tn) const {return outerchartspertrig.EntrySize(tn);};
    int GetOCPT(int tn, int vi) const {return outerchartspertrig.Get(tn,vi);};
    void SetOCPT(int tn, int vi, int ocn) {outerchartspertrig.Set(tn,vi,ocn);};
    void AddOCPT(int tn, int ocn) {outerchartspertrig.Add1(tn, ocn);};
    int TrigIsInOC(int tn, int ocn) const;
 
    //get chart number of a trig or 0 if unmarked
    int GetChartNr(int i) const;
    int GetMarker(int i) const 
    { return chartmark.Get(i); }
    void SetMarker(int nr, int m);
    int GetNOCharts() const;
    //get a chart from atlas
    const STLChart& GetChart(int nr) const;
    STLChart& GetChart(int nr) {return *(atlas.Get(nr));};
    int AtlasMade() const;
  
    void GetInnerChartLimes(Array<twoint>& limes, int chartnum);

    //FOR MESHING
    int GetMeshChartNr () { return meshchart; }
    void GetMeshChartBoundary (Array<Point2d > & points,
			       Array<Point3d > & points3d,
			       Array<INDEX_2> & lines, double h);


    Point<3> PointBetween(const Point<3> & p1, int t1, const Point<3> & p2, int t2);

    //select triangles in meshcharttrigs of actual (defined by trig) whole chart
    void PrepareSurfaceMeshing();
    //
    void DefineTangentialPlane(const Point<3> & ap1, const Point<3> & ap2, int trig);
    //
    void SelectChartOfTriangle (int trignum);
    //
    void SelectChartOfPoint (const Point<3> & p);
    //
    const Vec<3> & GetChartNormalVector () const { return meshtrignv; }

    // list of trigs
    void ToPlane (const Point<3> & locpoint, int * trigs, Point<2> & plainpoint, 
		  double h, int& zone, int checkchart);
    //return 0, wenn alles OK, 1 sonst
    int FromPlane (const Point<2> & plainpoint, Point<3> & locpoint, double h);
  
    //get nearest point in actual chart and return any triangle where it lies on
    int ProjectNearest(Point<3> & p3d) const;
    //project point with normal nv from last define tangential plane

    int LastTrig() const;
    int Project(Point<3> & p3d) const;
    int ProjectOnWholeSurface (Point<3> & p3d) const;

    int GetNLines() const {return lines.Size();}
    int AddLine(STLLine* line) {return lines.Append(line);}
    STLLine* GetLine(int nr) const {return lines.Get(nr);}
    int GetLineP(int lnr, int pnr) const {return lines.Get(lnr)->PNum(pnr);}
    int GetLineNP(int nr) const {return lines.Get(nr)->NP();}

    void SetLineEndPoint(int pn);
    int IsLineEndPoint(int pn);
    int LineEndPointsSet() const {return lineendpoints.Size() == GetNP();}
    void ClearLineEndPoints();

	DLL_HEADER void RestrictLocalH(class Mesh & mesh, double gh);
    void RestrictLocalHCurv(class Mesh & mesh, double gh);
    void RestrictHChartDistOneChart(int chartnum, Array<int>& acttrigs, class Mesh & mesh, 
				    double gh, double fact, double minh);

    friend class MeshingSTLSurface;


    virtual int GenerateMesh (shared_ptr<Mesh> & mesh, MeshingParameters & mparam,
			      int perfstepsstart, int perfstepsend);
    
    virtual const Refinement & GetRefinement () const;
  };
 

#include "meshstlsurface.hpp"



  extern int STLMeshingDummy (STLGeometry* stlgeometry, shared_ptr<Mesh> & mesh, MeshingParameters & mparam,
			      int perfstepsstart, int perfstepsend);


}
#endif
