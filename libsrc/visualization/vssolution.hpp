#ifndef FILE_VSSOLUTION
#define FILE_VSSOLUTION



#ifndef SMALLLIB
#ifndef NOTCL

extern 
void ImportSolution (const char * filename);

extern int Ng_Vis_Set (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[]);
#endif
#endif

class FieldLineCalc;

class VisualSceneSolution : public VisualScene
{
  friend class FieldLineCalc;
  
  class ClipPlaneTrig
  {
  public:
    struct ps 
    {
      int pnr, locpnr;
    };
    ps points[3];
    ElementIndex elnr;
  };

  class ClipPlanePoint
  {
  public:
    ElementIndex elnr;
    Point<3> lami;
    Point<3> p;
  };


  int surfellist;
  int linelist;
  int clipplanelist_scal;
  int clipplanelist_vec;
  int isolinelist;
  int clipplane_isolinelist;
  int surface_vector_list;
  // int cone_list;
  int isosurface_list;

  int pointcurvelist;

  bool draw_fieldlines;
  bool drawpointcurves;
  bool draw_isosurface;
  int num_fieldlines;
  bool fieldlines_randomstart;
  int fieldlineslist;
  int num_fieldlineslists;
  int fieldlines_startarea;
  Array<double> fieldlines_startarea_parameter;
  int fieldlines_startface;
  string fieldlines_filename;
  double fieldlines_reltolerance;
  int fieldlines_rktype;
  double fieldlines_rellength;
  double fieldlines_relthickness;
  int fieldlines_vecfunction;
  bool fieldlines_fixedphase;
  float fieldlines_phase;
  int fieldlines_maxpoints;


  int surfeltimestamp, clipplanetimestamp, solutiontimestamp;
  int surfellinetimestamp;
  int fieldlinestimestamp, surface_vector_timestamp;
  int pointcurve_timestamp;
  int isosurface_timestamp;
  int subdivision_timestamp;
  int timetimestamp;
  double minval, maxval;

  NgLock *lock;


#ifdef PARALLELGL
  Array<int> par_linelists;
  Array<int> par_surfellists;
#endif


public:

  enum EvalFunc { 
    FUNC_ABS = 1, 
    FUNC_ABS_TENSOR = 2,
    FUNC_MISES = 3, 
    FUNC_MAIN = 4
  };
  int evalfunc;

  enum SolType
    { 
      SOL_NODAL = 1, 
      SOL_ELEMENT = 2, 
      SOL_SURFACE_ELEMENT = 3, 
      SOL_NONCONTINUOUS = 4, 
      SOL_SURFACE_NONCONTINUOUS = 5,
      SOL_VIRTUALFUNCTION = 6,
      SOL_MARKED_ELEMENTS = 10,
      SOL_ELEMENT_ORDER = 11,
    };

  class SolData
  {
  public:
    SolData ();
    ~SolData ();
    
    char * name;
    double * data;
    int components;
    int dist;
    int order;
    bool iscomplex;
    bool draw_volume;
    bool draw_surface;
    SolType soltype;
    SolutionData * solclass;

    // internal variables:
    int size;
  };

  


  Array<SolData*> soldata;
  



  int usetexture;    // 0..no, 1..1D texture (standard), 2..2D-texture (complex)
  int clipsolution;  // 0..no, 1..scal, 2..vec
  int scalfunction, scalcomp, vecfunction;
  int gridsize;
  double xoffset, yoffset;

  int autoscale, logscale;
  double mminval, mmaxval;
  int numisolines;
  int subdivisions;

  bool showclipsolution;
  bool showsurfacesolution;
  bool lineartexture;
  int numtexturecols;

  int multidimcomponent;

  // bool fieldlineplot;
  double time;

  int deform;
  double scaledeform;
  bool imag_part;

private:
  void BuildFieldLinesFromFile(Array<Point3d> & startpoints);
  void BuildFieldLinesFromFace(Array<Point3d> & startpoints);
  void BuildFieldLinesFromBox(Array<Point3d> & startpoints);
  void BuildFieldLinesFromLine(Array<Point3d> & startpoints);

public:
  VisualSceneSolution ();
  virtual ~VisualSceneSolution ();

  virtual void BuildScene (int zoomall = 0);
  virtual void DrawScene ();
  virtual void MouseDblClick (int px, int py);

  void BuildFieldLinesPlot ();

  void AddSolutionData (SolData * soldata);
  void ClearSolutionData ();
  void UpdateSolutionTimeStamp ();
  SolData * GetSolData (int i);
  int GetNSolData () { return soldata.Size(); }

  void SaveSolutionData (const char * filename);


  static void RealVec3d (const double * values, Vec3d & v, 
			 bool iscomplex, bool imag);
  static void RealVec3d (const double * values, Vec3d & v, 
			 bool iscomplex, double phaser, double phasei);


  void SetSubdivision (int sd)
  {
    subdivisions = sd;
    subdivision_timestamp = solutiontimestamp = NextTimeStamp();
  }

  void GetMinMax (int funcnr, int comp, double & minv, double & maxv) const;

private:
  void GetClippingPlaneTrigs (Array<ClipPlaneTrig> & trigs, Array<ClipPlanePoint> & pts);
  void GetClippingPlaneGrid (Array<ClipPlanePoint> & pts);
  void DrawCone (const Point<3> & p1, const Point<3> & p2, double r);
  void DrawCylinder (const Point<3> & p1, const Point<3> & p2, double r);


  // Get Function Value, local coordinates lam1, lam2, lam3, 
  bool GetValue (const SolData * data, ElementIndex elnr, 
		   double lam1, double lam2, double lam3,
		   int comp, double & val) const;

  bool GetValue (const SolData * data, ElementIndex elnr,
		 const double xref[], const double x[], const double dxdxref[], 
		 int comp, double & val) const;

  bool GetValueComplex (const SolData * data, ElementIndex elnr, 
			double lam1, double lam2, double lam3,
			int comp, complex<double> & val) const;

  bool GetValues (const SolData * data, ElementIndex elnr, 
		  double lam1, double lam2, double lam3,
		  double * values) const;

  bool GetValues (const SolData * data, ElementIndex elnr, 
		  const double xref[], const double x[], const double dxdxref[], 
		  double * values) const;

  bool GetMultiValues (const SolData * data, ElementIndex elnr, int npt,
		       const double * xref, int sxref,
		       const double * x, int sx,
		       const double * dxdxref, int sdxdxref,
		       double * val, int sval) const;


  bool GetSurfValue (const SolData * data, SurfaceElementIndex elnr,
		     double lam1, double lam2, 
		     int comp, double & val) const;

  bool GetSurfValue (const SolData * data, SurfaceElementIndex elnr,
		     const double xref[], const double x[], const double dxdxref[], 
		     int comp, double & val) const;

  
  bool GetSurfValueComplex (const SolData * data, SurfaceElementIndex elnr,
			    double lam1, double lam2, 
			    int comp, complex<double> & val) const;

  bool GetSurfValues (const SolData * data, SurfaceElementIndex elnr,
		      double lam1, double lam2, 
		      double * values) const;

  bool GetSurfValues (const SolData * data, SurfaceElementIndex elnr,
		      const double xref[], const double x[], const double dxdxref[], 
		      double * values) const;

  bool GetMultiSurfValues (const SolData * data, SurfaceElementIndex elnr, int npt,
                           const double * xref, int sxref,
                           const double * x, int sx,
                           const double * dxdxref, int sdxdxref,
                           double * val, int sval) const;
  
  double ExtractValue (const SolData * data, int comp, double * values) const;
  complex<double> ExtractValueComplex (const SolData * data, int comp, double * values) const;


  Vec<3> GetDeformation (ElementIndex elnr, const Point<3> & p) const;
  Vec<3> GetSurfDeformation (SurfaceElementIndex selnr, double lam1, double lam2) const;

  void GetPointDeformation (int pnum, Point<3> & p, SurfaceElementIndex elnr = -1) const;

public:
  /// draw elements (build lists)
  void DrawSurfaceElements ();
  void DrawSurfaceElementLines ();
  void DrawSurfaceVectors ();
  void DrawTrigSurfaceVectors(const Array< Point<3> > & lp, const Point<3> & pmin, const Point<3> & pmax,
			      const int sei, const SolData * vsol);
  void DrawIsoSurface(const SolData * sol, const SolData * grad, int comp);
  
  void DrawIsoLines (const Point<3> & p1, 
		     const Point<3> & p2, 
		     const Point<3> & p3,
		     double val1, double val2, double val3);

  // draw isolines between lines (p1,p2) and (p3,p4)
  void DrawIsoLines2 (const Point<3> & p1, 
		      const Point<3> & p2, 
		      const Point<3> & p3,
		      const Point<3> & p4,
		      double val1, double val2, double val3, double val4);


  void DrawClipPlaneTrigs (); // const SolData * sol, int comp);
		  
  void SetOpenGlColor(double val);

  // 0 .. non, 1 .. scalar, 2 .. complex
  void SetTextureMode (int texturemode) const;

#ifndef SMALLLIB  
#ifndef NOTCL

  friend int Ng_Vis_Set (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[]);

#endif
#endif


#ifdef PARALLELGL
  void Broadcast ();
#endif


};




class RKStepper
{
private:
  Array<double> c,b;
  TABLE<double> *a;
  int steps;
  int order;
  
  double tolerance;
  
  Array<Vec3d> K;
  
  int stepcount;
  
  double h;
  double startt;
  double startt_bak;
  Point3d startval;
  Point3d startval_bak;
  
  bool adaptive;
  int adrun;
  Point3d valh;
  
  int notrestarted;

public:
  
  ~RKStepper();
    
  RKStepper(int type = 0);

  void SetTolerance(const double tol){tolerance = tol;}
        
  void StartNextValCalc(const Point3d & astartval, const double astartt, const double ah, const bool aadaptive = false);

  bool GetNextData(Point3d & val, double & t, double & ah);

  bool FeedNextF(const Vec3d & f);
};





class FieldLineCalc
{
private:
  const Mesh & mesh;
  
  VisualSceneSolution & vss;
  
  const VisualSceneSolution::SolData * vsol;

  RKStepper stepper;

  double maxlength;

  int maxpoints;
  
  int direction;
  
  Point3d pmin, pmax;
  double rad;
  double phaser, phasei;
  
  double critical_value;

  bool randomized;

  double thickness;

public:
  FieldLineCalc(const Mesh & amesh, VisualSceneSolution & avss, const VisualSceneSolution::SolData * solution, 
		const double rel_length, const int amaxpoints = -1, 
		const double rel_thickness = -1, const double rel_tolerance = -1, const int rk_type = 0, const int adirection = 0);

  void SetPhase(const double real, const double imag) { phaser = real; phasei = imag; }
  
  void SetCriticalValue(const double val) { critical_value = val; }

  void Randomized(void) { randomized = true; }
  void NotRandomized(void) { randomized = false; }

  void Calc(const Point3d & startpoint, Array<Point3d> & points, Array<double> & vals, Array<bool> & drawelems, Array<int> & dirstart);
  
  void GenerateFieldLines(Array<Point3d> & potential_startpoints, const int numlines, const int gllist,
			  const double minval, const double maxval, const int logscale, double phaser, double phasei);
};




extern VisualSceneSolution vssolution;




#endif

