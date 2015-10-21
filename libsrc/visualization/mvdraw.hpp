#ifndef FILE_MVDRAW
#define FILE_MVDRAW


namespace netgen
{


  class VisualScene
  {
  protected:
    static DLL_HEADER Point3d center;
    static DLL_HEADER double rad;

    static double lookatmat[16];
    static double transmat[16];
    static double rotmat[16];
    static double centermat[16];
	
    static DLL_HEADER double transformationmat[16];

    GLdouble clipplane[4];

    int changeval;
    static DLL_HEADER GLdouble backcolor;

    static int DLL_HEADER selface;
    static int selelement;
    static int DLL_HEADER selpoint;
    static int selpoint2;
    static int locpi;
    static int DLL_HEADER seledge;

    static int selecttimestamp;

  public:
    static int viewport[4];
    static GLuint coltexname;
    static int ntexcols;
    int invcolor;


  public:
    DLL_HEADER VisualScene ();
    DLL_HEADER virtual ~VisualScene();

    DLL_HEADER virtual void BuildScene (int zoomall = 0);
    DLL_HEADER virtual void DrawScene ();
  
    DLL_HEADER void CalcTransformationMatrices();
    DLL_HEADER void StandardRotation (const char * dir);
    DLL_HEADER void ArbitraryRotation (const Array<double> & alpha, const Array<Vec3d> & vec);
    DLL_HEADER void ArbitraryRotation (const double alpha, const Vec3d & vec);

    DLL_HEADER void MouseMove(int oldx, int oldy,
                   int newx, int newy,
                   char mode);

    DLL_HEADER void LookAt (const Point<3> & cam, const Point<3> & obj,
                 const Point<3> & camup);

    DLL_HEADER void SetClippingPlane ();

    DLL_HEADER virtual void MouseDblClick (int px, int py);

    DLL_HEADER void SetLight ();
    static void SetBackGroundColor (double col)
    { backcolor = col; }

    DLL_HEADER void CreateTexture (int ncols, int linear, double alpha, int typ);
    DLL_HEADER void DrawColorBar (double minval, double maxval, int logscale = 0, bool linear = 1);
    DLL_HEADER void DrawCoordinateCross ();
    DLL_HEADER void DrawNetgenLogo ();
    DLL_HEADER void SetOpenGlColor(double val, double valmin, double valmax, int logscale = 0);
	

#ifdef PARALLELGL
    DLL_HEADER void InitParallelGL ();
    DLL_HEADER void Broadcast ();
#endif 
  };


  DLL_HEADER extern void MyOpenGLText (const char * text);
  DLL_HEADER extern void Set_OpenGLText_Callback ( void (*fun) (const char * text) );










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

    MouseEventHandler * user_me_handler;

    NgLock *lock;

    //  int selface, selelement;
    //  int selpoint, selpoint2, locpi;
    //  int seledge;

    double minh, maxh; // for meshsize coloring

    // weak_ptr<Mesh> wp_mesh;

  public:
    DLL_HEADER VisualSceneMesh ();
	DLL_HEADER virtual ~VisualSceneMesh ();

	DLL_HEADER virtual void BuildScene (int zoomall = 0);
	DLL_HEADER virtual void DrawScene ();
	DLL_HEADER virtual void MouseDblClick (int px, int py);

    // void SetMesh (shared_ptr<Mesh> mesh) { wp_mesh = mesh; }
    // shared_ptr<Mesh> GetMesh () { return shared_ptr<Mesh>(wp_mesh); }
    shared_ptr<Mesh> GetMesh () const { return shared_ptr<Mesh>(global_mesh); }

    void SetMouseEventHandler (MouseEventHandler * handler)
    { user_me_handler = handler; }


	DLL_HEADER int SelectedFace () const
    { return selface; }
	DLL_HEADER void SetSelectedFace (int asf);
    //    { selface = asf; selecttimestamp = GetTimeStamp(); }

	DLL_HEADER int SelectedEdge () const
    { return seledge; }
	DLL_HEADER int SelectedElement () const
    { return selelement; }
	DLL_HEADER int SelectedPoint () const
    { return selpoint; }
    void BuildFilledList (bool names);
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

  DLL_HEADER extern VisualSceneMesh vsmesh;
  

  class DLL_HEADER VisualSceneSpecPoints : public VisualScene
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
                            const double * transformationmat,
                            const Point3d & center,
                            const double rad,
                            const int displaylist,
                            int & selelement, int & selface, int & seledge, int & selpoint,
                            int & selpoint2, int & locpi);


}


#endif

