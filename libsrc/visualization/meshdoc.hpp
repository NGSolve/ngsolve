namespace netgen

{

class VisualSceneMeshDoctor : public VisualScene
{
  int filledlist;
  int outlinelist;
  int edgelist;

  int selelement, locpi;
  int selpoint, selpoint2;

  // for edgemarking:
  Array<int> edgedist;
  int markedgedist;
  

public:
  DLL_HEADER VisualSceneMeshDoctor ();
  DLL_HEADER virtual ~VisualSceneMeshDoctor ();

  DLL_HEADER virtual void BuildScene (int zoomall = 0);
  DLL_HEADER virtual void DrawScene ();
  DLL_HEADER virtual void MouseDblClick (int px, int py);

  DLL_HEADER void SetMarkEdgeDist (int dist);
  DLL_HEADER void ClickElement (int elnr);
  DLL_HEADER void UpdateTables ();
  DLL_HEADER int IsSegmentMarked (int segnr) const;
};

class MeshDoctorParameters 
{
public:
  int active;
};


DLL_HEADER extern MeshDoctorParameters meshdoctor;

}
