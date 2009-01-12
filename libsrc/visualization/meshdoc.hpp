
class VisualSceneMeshDoctor : public VisualScene
{
  int filledlist;
  int outlinelist;
  int edgelist;

  int selelement, locpi;
  int selpoint, selpoint2;

  // for edgemarking:
  ARRAY<int> edgedist;
  int markedgedist;
  

public:
  VisualSceneMeshDoctor ();
  virtual ~VisualSceneMeshDoctor ();

  virtual void BuildScene (int zoomall = 0);
  virtual void DrawScene ();
  virtual void MouseDblClick (int px, int py);

  void SetMarkEdgeDist (int dist);
  void ClickElement (int elnr);
  void UpdateTables ();
  int IsSegmentMarked (int segnr) const;
};

class MeshDoctorParameters 
{
public:
  int active;
};


extern MeshDoctorParameters meshdoctor;
