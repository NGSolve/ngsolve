#ifndef FILE_MVDRAW
#define FILE_MVDRAW


namespace netgen
{


  extern void InitDrawMesh ();
  extern void DrawMesh ();
  extern void MouseMove(int oldx, int oldy,
                        int newx, int newy,
                        char mode);

  extern void Render ();




  class VisualScene
  {
  protected:
    static Point3d center;
    static double rad;

    static float lookatmat[16];
    static float transmat[16];
    static float rotmat[16];
    static float centermat[16];
    static float transformationmat[16];

    GLdouble clipplane[4];

    int changeval;
    static GLdouble backcolor;

    static int selface;
    static int selelement;
    static int selpoint;
    static int selpoint2;
    static int locpi;
    static int seledge;

    static int selecttimestamp;

  public:

    // static GLubyte * colortexture;
    static GLuint coltexname;
    static int ntexcols;
    // static bool linear_colors;
    int invcolor;


  public:
    VisualScene ();
    virtual ~VisualScene();

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
  
    void CalcTransformationMatrices();
    void StandardRotation (const char * dir);
    void ArbitraryRotation (const Array<double> & alpha, const Array<Vec3d> & vec);
    void ArbitraryRotation (const double alpha, const Vec3d & vec);

    void MouseMove(int oldx, int oldy,
                   int newx, int newy,
                   char mode);

    void LookAt (const Point<3> & cam, const Point<3> & obj,
                 const Point<3> & camup);

    void SetClippingPlane ();

    virtual void MouseDblClick (int px, int py);

    void SetLight ();
    static void SetBackGroundColor (double col)
    { backcolor = col; }

    void CreateTexture (int ncols, int linear, int typ = GL_DECAL);
    void DrawColorBar (double minval, double maxval, int logscale = 0, bool linear = 1);
    void DrawCoordinateCross ();
    void DrawNetgenLogo ();
    void SetOpenGlColor(double val, double valmin, double valmax, int logscale = 0);


#ifdef PARALLELGL
    void InitParallelGL ();
    void Broadcast ();
#endif 
  };


  extern void MyOpenGLText (const char * text);



  class VisualSceneGeometry : public VisualScene
  {
    Array<int> trilists;
    int selsurf;
  public:
    VisualSceneGeometry ();
    virtual ~VisualSceneGeometry ();

    virtual void SelectSurface (int aselsurf);
    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
  };



  class VisualSceneSTLGeometry : public VisualScene
  {
    Array<int> trilists;
  
  public:
    VisualSceneSTLGeometry ();
    virtual ~VisualSceneSTLGeometry ();

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
  };


  class VisualSceneGeometry2d : public VisualScene
  {
  public:
    VisualSceneGeometry2d ();
    virtual ~VisualSceneGeometry2d ();

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
  };


#ifdef OCCGEOMETRY
  class VisualSceneOCCGeometry : public VisualScene
  {
    Array<int> trilists;
    Array<int> linelists;
    int selsurf;
  public:
    VisualSceneOCCGeometry ();
    virtual ~VisualSceneOCCGeometry ();

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
    virtual void MouseDblClick (int px, int py);
  };
#endif










#ifdef STEP
  class VisualSceneSTEPGeometry : public VisualScene
  {
    Array<int> gllists;
  
  public:
    VisualSceneSTEPGeometry ();
    virtual ~VisualSceneSTEPGeometry ();

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
  };
#endif


  class VisualSceneSTLMeshing : public VisualScene
  {
    Array<int> trilists;
    int selecttrig, nodeofseltrig;

  public:
    VisualSceneSTLMeshing ();
    virtual ~VisualSceneSTLMeshing ();

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
    virtual void MouseDblClick (int px, int py);

    int seltria;
  };




  class VisualSceneSurfaceMeshing : public VisualScene
  {
  public:
    VisualSceneSurfaceMeshing ();
    virtual ~VisualSceneSurfaceMeshing ();

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
  };







  class VisualSceneMesh : public VisualScene
  {
    int filledlist;
    int linelist;
    int edgelist;
    int pointnumberlist;

    int tetlist;
    int prismlist;
    int pyramidlist;
    int hexlist;

    int badellist;
    int identifiedlist;
    int domainsurflist;

    int vstimestamp;//, selecttimestamp;
    int filledtimestamp;
    int linetimestamp;
    int edgetimestamp;
    int pointnumbertimestamp;

    int tettimestamp;
    int prismtimestamp;
    int pyramidtimestamp;
    int hextimestamp;

    int badeltimestamp;
    int identifiedtimestamp;
    int domainsurftimestamp;


#ifdef PARALLELGL
    Array<int> par_linelists;
    Array<int> par_filledlists;
#endif


    NgLock *lock;

    //  int selface, selelement;
    //  int selpoint, selpoint2, locpi;
    //  int seledge;

    double minh, maxh; // for meshsize coloring

  public:
    VisualSceneMesh ();
    virtual ~VisualSceneMesh ();

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();
    virtual void MouseDblClick (int px, int py);

    int SelectedFace () const
    { return selface; }
    void SetSelectedFace (int asf);
    //    { selface = asf; selecttimestamp = GetTimeStamp(); }

    int SelectedEdge () const
    { return seledge; }
    int SelectedElement () const
    { return selelement; }
    int SelectedPoint () const
    { return selpoint; }
    void BuildFilledList();
    // private:
    void BuildLineList();
    void BuildEdgeList();
    void BuildPointNumberList();

    void BuildTetList();
    void BuildPrismList();
    void BuildPyramidList();
    void BuildHexList();

    void BuildBadelList();
    void BuildIdentifiedList();
    void BuildDomainSurfList();
  };







  class VisualSceneSpecPoints : public VisualScene
  {
  public:
    VisualSceneSpecPoints ();
    virtual ~VisualSceneSpecPoints ();

    virtual void BuildScene (int zoomall = 0);
    virtual void DrawScene ();

    double len;
  };

  // extern struct Tcl_Interp * hinterp;


  extern void AddVisualizationScene (const string & name, 
                                     VisualScene * vs);


  void MouseDblClickSelect (const int px, const int py,
                            const GLdouble * clipplane, const GLdouble backcolor,
                            const float * transformationmat,
                            const Point3d & center,
                            const double rad,
                            const int displaylist,
                            int & selelement, int & selface, int & seledge, int & selpoint,
                            int & selpoint2, int & locpi);


}


#endif

