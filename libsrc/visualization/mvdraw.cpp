#include <mystdlib.h>
#include <myadt.hpp>
#include <meshing.hpp>
#include <csg.hpp>
#include <geometry2d.hpp>
#include <stlgeom.hpp>

#include <visual.hpp>
// #include <parallel.hpp>




#ifndef WIN32
#define GLX_GLXEXT_LEGACY

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>  /* for XA_RGB_DEFAULT_MAP atom */
// #include <GL/glx.h>    // for parallel GL ???
#endif





namespace netgen
{
  Point3d VisualScene :: center;
  double VisualScene :: rad;
  GLdouble VisualScene :: backcolor;

  /*
#if TOGL_MAJOR_VERSION!=2
  GLuint VisualScene :: fontbase = 0;
#else
  Tcl_Obj * VisualScene :: fontbase = NULL;
  Togl * VisualScene :: globtogl;
#endif
  */

  // texture for color decoding
  // GLubyte * VisualScene :: colortexture = NULL;
  GLuint VisualScene :: coltexname = 1;
  int VisualScene :: ntexcols = -1;


  float VisualScene :: lookatmat[16];
  float VisualScene :: transmat[16];
  float VisualScene :: rotmat[16];
  float VisualScene :: centermat[16];
  float VisualScene :: transformationmat[16];

  int VisualScene :: selface;
  int VisualScene :: selelement;
  int VisualScene :: selpoint;
  int VisualScene :: selpoint2;
  int VisualScene :: locpi;
  int VisualScene :: seledge;

  int VisualScene :: selecttimestamp;


  VisualizationParameters :: VisualizationParameters()
  {
    lightamb = 0.3;
    lightdiff = 0.7;
    lightspec = 1;
    shininess = 50;
    transp = 0.3;
    locviewer = 0;
    showstltrias = 0;
    centerpoint = 0;
    usedispllists = 1;
    strcpy (selectvisual, "cross");

    use_center_coords = false;
  };
  VisualizationParameters vispar;



  double dist = 0;
  // double dist = 6;
  // vorher: pnear = 2;
  double pnear = 0.1;
  double pfar = 10;



  extern STLGeometry * stlgeometry;
  extern AutoPtr<SplineGeometry2d> geometry2d;
  extern AutoPtr<Mesh> mesh;
  extern Array<SpecialPoint> specpoints;
  extern Array<Box<3> > boxes;

  VisualScene :: VisualScene ()
  {
    changeval = -1;
    backcolor = 0;
  }


  VisualScene :: ~VisualScene()
  {
    ;
  }


  void Render ()
  {
    multithread.redraw = 1;
  }


  void VisualScene :: BuildScene (int zoomall)
  {
    center = Point3d (0,0,0);
    rad = 1;

    CalcTransformationMatrices();

    glEnable(GL_DEPTH_TEST);
    glDisable (GL_DITHER);
  
    GLfloat ambvals[] = { 0.4f, 0.4f, 0.4f, 1.0f };
    GLfloat diffvals[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat specvals[] =  { 0.7f, 0.7f, 0.7f, 1.0f };
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambvals);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffvals);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specvals);
  
    GLfloat light_position[] = { 1, 3, 3, 0 };
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  
    glLightModeli (GL_LIGHT_MODEL_TWO_SIDE, 0);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
  }


  void VisualScene :: DrawScene ()
  {
    if (changeval == -1)
      BuildScene();
    changeval = 0;

    glClearColor(backcolor, backcolor, backcolor, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable (GL_COLOR_MATERIAL);
    glColor3f (1.0f, 1.0f, 1.0f);
    glLineWidth (1.0f);

    DrawCoordinateCross ();
    DrawNetgenLogo ();
    glFinish();  
  }


  void VisualScene :: CalcTransformationMatrices()
  {
    // prepare model view matrix
  
    glPushMatrix();

    glLoadIdentity();
    gluLookAt (0, 0, 6, 0, 0, 0, 0, 1, 0);
    glGetFloatv (GL_MODELVIEW_MATRIX, lookatmat);

    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, -dist);
    glGetFloatv (GL_MODELVIEW_MATRIX, transmat);
  
    glLoadIdentity();
    glGetFloatv (GL_MODELVIEW_MATRIX, rotmat);

    glScalef (1/rad, 1/rad, 1/rad);
    glTranslated (-center.X(), -center.Y(), -center.Z());
    glGetFloatv (GL_MODELVIEW_MATRIX, centermat);

    glLoadIdentity();
    glMultMatrixf (lookatmat);
    glMultMatrixf (transmat);
    glMultMatrixf (rotmat);
    glMultMatrixf (centermat);
    glGetFloatv (GL_MODELVIEW_MATRIX, transformationmat);

    glPopMatrix();
  }


  void VisualScene :: ArbitraryRotation (const Array<double> & alpha, const Array<Vec3d> & vec)
  {
    glPushMatrix();

    glLoadIdentity();

    for(int i=0; i<alpha.Size() && i<vec.Size(); i++)
      {
	glRotatef(alpha[i], vec[i].X(), vec[i].Y(), vec[i].Z());
      }

    glGetFloatv (GL_MODELVIEW_MATRIX, rotmat);

    glLoadIdentity();
    glMultMatrixf (lookatmat);
    glMultMatrixf (transmat);
    glMultMatrixf (rotmat);
    glMultMatrixf (centermat);
    glGetFloatv (GL_MODELVIEW_MATRIX, transformationmat);
  
    glPopMatrix();
  } 



  void VisualScene :: ArbitraryRotation (const double alpha, const Vec3d & vec)
  {
    Array<double> a(1); a[0] = alpha;
    Array<Vec3d> v(1); v[0] = vec;

    ArbitraryRotation(a,v);
  } 

  void VisualScene :: StandardRotation (const char * dir)
  {
    glPushMatrix();

    glLoadIdentity();
  
    if (strcmp (dir, "xy") == 0)
      ;
    else if (strcmp (dir, "yx") == 0)
      glRotatef(180.0, 1.0f, 1.0f, 0.0f);    
    else if (strcmp (dir, "xz") == 0)
      glRotatef(-90.0, 1.0f, 0.0f, 0.0f);    
    else if (strcmp (dir, "zx") == 0)
      {
	glRotatef(180.0, 1.0f, 1.0f, 0.0f);    
	glRotatef(-90.0, 1.0f, 0.0f, 0.0f);    
      }
    else if (strcmp (dir, "yz") == 0)
      {
	glRotatef(-90.0, 0.0f, 0.0f, 1.0f);    
	glRotatef(-90.0, 0.0f, 1.0f, 0.0f);    
      }
    else if (strcmp (dir, "zy") == 0)
      glRotatef(90.0, 0.0f, 1.0f, 0.0f);    


    glGetFloatv (GL_MODELVIEW_MATRIX, rotmat);

    glLoadIdentity();
    glMultMatrixf (lookatmat);
    glMultMatrixf (transmat);
    glMultMatrixf (rotmat);
    glMultMatrixf (centermat);
    glGetFloatv (GL_MODELVIEW_MATRIX, transformationmat);
  
    glPopMatrix();
  }

  void VisualScene :: MouseMove(int oldx, int oldy,
				int newx, int newy,
				char mode)
  {
    int deltax = newx - oldx;
    int deltay = newy - oldy;
  
    glPushMatrix();
    glLoadIdentity ();
  
    switch (mode)
      {
      case 'r':
	{	
	  glRotatef(float(deltax)/2, 0.0f, 1.0f, 0.0f);
	  glRotatef(float(deltay)/2, 1.0f, 0.0f, 0.0f);
	  glMultMatrixf (rotmat);
	  glGetFloatv (GL_MODELVIEW_MATRIX, rotmat);
	  break;
	}
      case 'm':
	{
	  GLdouble projmat[16], modelviewmat[16];
	  GLint viewport[4];
	  glGetDoublev (GL_PROJECTION_MATRIX, projmat);
	  glGetDoublev (GL_MODELVIEW_MATRIX, modelviewmat);
	  glGetIntegerv (GL_VIEWPORT, viewport);
	
	  // vorher pvz1/2 = 0
	  GLdouble pvx1 = 0, pvy1 = 0, pvz1 = 0.99; //  0.95;
	  GLdouble pvx2 = deltax, pvy2 = -deltay, pvz2 = 0.99; // 0.95;

	  GLdouble px1, py1, pz1;
	  GLdouble px2, py2, pz2;
	
	  gluUnProject (pvx1, pvy1, pvz1, 
			modelviewmat, projmat, viewport,
			&px1, &py1, &pz1);
	  gluUnProject (pvx2, pvy2, pvz2, 
			modelviewmat, projmat, viewport,
			&px2, &py2, &pz2);
	  /*
	    gluUnProject (oldx, oldy, 1, 
	    modelviewmat, projmat, viewport,
	    &px1, &py1, &pz1);
	    gluUnProject (newx, newy, 1, 
	    modelviewmat, projmat, viewport,
	    &px2, &py2, &pz2);
	  */

	  /*	
	    cout << "pv1 = " << pvx1 << ", " << pvy1 << ", " << pvz1 << endl;
	    cout << "p1 = " << px1 << ", " << py1 << ", " << pz1 << endl;
	  */

	  glTranslated (px2-px1, py2-py1, pz2-pz1);
	
	  glMultMatrixf (transmat);
	  glGetFloatv (GL_MODELVIEW_MATRIX, transmat);
	  break;
	}
      case 'z':
	{
	  // glTranslatef(0.0f, 0.0f, -dist);

	  // cout << "deltay = " << deltay << endl;
	  // cout << "float_bug = " << (float(deltay)/100) << endl;   gives wrong result with icc 9.0.021
	  glScaled (exp (double (-deltay)/100), 
		    exp (double (-deltay)/100), 
		    exp (double (-deltay)/100));
	  // glTranslatef(0.0f, 0.0f, dist);
	  glMultMatrixf (transmat);
	  glGetFloatv (GL_MODELVIEW_MATRIX, transmat);
	  break;
	}
      }

    glLoadIdentity();
    glMultMatrixf (lookatmat);
    glMultMatrixf (transmat);
    glMultMatrixf (rotmat);
    glMultMatrixf (centermat);
    glGetFloatv (GL_MODELVIEW_MATRIX, transformationmat);
  
    glPopMatrix();
  }


  void VisualScene :: LookAt (const Point<3> & cam, const Point<3> & obj,
			      const Point<3> & camup)
  {
    glPushMatrix();
    glLoadIdentity ();
    gluLookAt (cam(0), cam(1), cam(2), 
	       obj(0), obj(1), obj(2),
	       camup(0), camup(1), camup(2));
    glMultMatrixf (centermat);
    glGetFloatv (GL_MODELVIEW_MATRIX, transformationmat);
    glPopMatrix();
  }

  
  void VisualScene :: SetClippingPlane ()
  {
    if (vispar.clipenable)
      {
	Vec3d n = vispar.clipnormal;
	n /= (n.Length()+1e-10);
	clipplane[0] = n.X();
	clipplane[1] = n.Y();
	clipplane[2] = n.Z();
	clipplane[3] = -(Vec3d(center) * n) + rad * vispar.clipdist;

	glClipPlane(GL_CLIP_PLANE0, clipplane);
	glEnable(GL_CLIP_PLANE0);
      }
    else
      glDisable (GL_CLIP_PLANE0);
  }




  void VisualScene :: MouseDblClick (int /* px */, int /* py */)
  {
    ;
  }



  void VisualScene :: SetLight()
  {
    GLfloat vals[3];
    double lightamb = vispar.lightamb;
    vals[0] = vals[1] = vals[2] = lightamb;
    glLightfv(GL_LIGHT0, GL_AMBIENT, vals);

    double lightdiff = vispar.lightdiff;
    vals[0] = vals[1] = vals[2] = lightdiff;
    glLightfv(GL_LIGHT0, GL_DIFFUSE, vals);

    double lightspec = vispar.lightspec;
    vals[0] = vals[1] = vals[2] = lightspec;
    glLightfv(GL_LIGHT0, GL_SPECULAR, vals);

    glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, vispar.shininess);
    glLightModeli (GL_LIGHT_MODEL_LOCAL_VIEWER, vispar.locviewer);

    float mat_spec_col[] = { 1, 1, 1, 1 };
    glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_spec_col);

    glEnable (GL_LIGHTING);
    glEnable (GL_LIGHT0);
  }




  void VisualScene :: SetOpenGlColor(double val, double valmin, double valmax,
				     int logscale)
  {
    double value;

    if (!logscale)
      value = (val - valmin) / (valmax - valmin);
    else
      {
	if (valmax <= 0) valmax = 1;
	if (valmin <= 0) valmin = 1e-4 * valmax;
	value = (log(fabs(val)) - log(valmin)) / (log(valmax) - log(valmin));
      }

    if (!invcolor)
      value = 1 - value;

    glTexCoord1f ( 0.998 * value + 0.001);
    // glTexCoord1f ( val ); 

    glTexCoord2f ( 0.998 * value + 0.001, 1.5);
    // glTexCoord1f ( value ); 

    if (value > 1) value = 1;
    if (value < 0) value = 0;

    value *= 4;

    static const double colp[][3] =
      {
	{ 1, 0, 0 },
	{ 1, 1, 0 },
	{ 0, 1, 0 },
	{ 0, 1, 1 },
	{ 0, 0, 1 },
	//	{ 1, 0, 1 },
	//	{ 1, 0, 0 },
      };
  
    int i = int(value);
    double r = value - i;

    GLdouble col[3];
    for (int j = 0; j < 3; j++)
      col[j] = (1-r) * colp[i][j] + r * colp[i+1][j];
  
    glColor3d (col[0], col[1], col[2]);
  }


  /*
  void VisualScene :: CreateTexture (int ncols, int linear, int typ)
  {

    static const double colp[][3] =
      {
	{ 1, 0, 0 },
	{ 1, 1, 0 },
	{ 0, 1, 0 },
	{ 0, 1, 1 },
	{ 0, 0, 1 },
      };
  


    if (ntexcols != 1024)
      {
	ntexcols = 1024;
	
	// glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
	glGenTextures (1, &coltexname);
	glBindTexture (GL_TEXTURE_1D, coltexname);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
	
	for (int level = 0; level <= 11; level++)
	  {
	    ncols = 2048 >> level;
	    cout << "ncols = " << ncols << endl;

	    colortexture = new GLubyte[4*ncols+12];

	    for (int i = 0; i < ncols; i++)
	      {
		double value = 4.0 * i / (ncols-1);
		
		int iv = int(value);
		double r = value - iv;
		
		GLdouble col[3];
		for (int j = 0; j < 3; j++)
		  col[j] = (1-r) * colp[iv][j] + r * colp[iv+1][j];
		
		colortexture[4*i] = GLubyte (255 * col[0]);
		colortexture[4*i+1] = GLubyte (255 * col[1]);
		colortexture[4*i+2] = GLubyte (255 * col[2]);
		colortexture[4*i+3] = GLubyte(255);
		
		if (ncols > 20)
		  if ( i % (ncols / 10) == 0)
		    {
		      colortexture[4*i] = GLubyte (0);
		      colortexture[4*i+1] = GLubyte (0);
		      colortexture[4*i+2] = GLubyte (0);
		      colortexture[4*i+4] = GLubyte (0);
		      colortexture[4*i+5] = GLubyte (0);
		      colortexture[4*i+6] = GLubyte (0);
		    }
	      }
	    
	    glTexImage1D (GL_TEXTURE_1D, level, 4, ncols, 0, GL_RGBA, GL_UNSIGNED_BYTE, colortexture);
	  }
      }
    
    glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, typ);  // DECAL or MODULATE
    glBindTexture (GL_TEXTURE_1D, coltexname);
  }
  */





  void VisualScene :: CreateTexture (int ncols, int linear, int typ)
  {
    if (linear) ncols = 32;
    else   ncols = 8;


    if (ntexcols != ncols) 
      {
	ntexcols = ncols;
      
	GLubyte colortexture[4*32];

	const double colp[][3] =
	  {
	    { 1, 0, 0 },
	    { 1, 1, 0 },
	    { 0, 1, 0 },
	    { 0, 1, 1 },
	    { 0, 0, 1 },
	  };
  
	for (int i = 0; i < ncols; i++)
	  {
	    double value = 4.0 * i / (ncols-1);

	    int iv = int(value);
	    double r = value - iv;

	    GLdouble col[3];

	    if(r > 1e-3)
	      for (int j = 0; j < 3; j++)
		col[j] = (1.-r) * colp[iv][j] + r * colp[iv+1][j];
	    else
	      for (int j = 0; j < 3; j++)
		col[j] = colp[iv][j];

	    colortexture[4*i] = GLubyte (255 * col[0]);
	    colortexture[4*i+1] = GLubyte (255 * col[1]);
	    colortexture[4*i+2] = GLubyte (255 * col[2]);
	    colortexture[4*i+3] = GLubyte(255);
	  }

	// glPixelStorei (GL_UNPACK_ALIGNMENT, 1);

     	glTexImage1D (GL_TEXTURE_1D, 0, 4, ncols, 0, GL_RGBA, GL_UNSIGNED_BYTE, colortexture);
	glTexImage2D (GL_TEXTURE_2D, 0, 4, ncols, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, colortexture);

	glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, typ);  // DECAL or MODULATE
	
	GLfloat bcol[] = { 1, 1, 1, 1.0 };
	glTexParameterfv (GL_TEXTURE_1D, GL_TEXTURE_BORDER_COLOR, bcol);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);

	glTexParameterfv (GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, bcol);
	glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	
	if (linear)
	  {
	    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

	    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	  }
	else
	  {
	    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

	    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	  }
      }
  }
  




  /*
  void VisualScene :: CreateTexture (int ncols, int linear, int typ)
  {
    if (ncols < 2) ncols = 2;

    if (linear) ncols = 32;
    else   ncols = 8;

    if (ntexcols != ncols) 
      {
	if (colortexture) 
	  {
	    glDeleteTextures (1, &coltexname);
	    delete colortexture;
	  }
      
	ntexcols = ncols;
      
	colortexture = new GLubyte[4*ncols+12];

	const double colp[][3] =
	  {
	    { 1, 0, 0 },
	    { 1, 1, 0 },
	    { 0, 1, 0 },
	    { 0, 1, 1 },
	    { 0, 0, 1 },
	  };
  
	for (int i = 0; i < ncols; i++)
	  {
	    double value = 4.0 * i / (ncols-1);

	    int iv = int(value);
	    double r = value - iv;

	    GLdouble col[3];
	    for (int j = 0; j < 3; j++)
	      col[j] = (1-r) * colp[iv][j] + r * colp[iv+1][j];
	  
	    colortexture[4*i+4] = GLubyte (255 * col[0]);
	    colortexture[4*i+5] = GLubyte (255 * col[1]);
	    colortexture[4*i+6] = GLubyte (255 * col[2]);
	    colortexture[4*i+7] = GLubyte(255);
	  }
	for (int j = 0; j < 4; j++)
	  {
	    colortexture[j] = colortexture[4+j];
	    colortexture[ncols*4+4+j] = colortexture[ncols*4+j];
	  }

	//      glPixelStorei (GL_UNPACK_ALIGNMENT, 1);

	glGenTextures (1, &coltexname);
	glBindTexture (GL_TEXTURE_1D, coltexname);

	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

	glTexImage1D (GL_TEXTURE_1D, 0, 4, ncols, 0, GL_RGBA, GL_UNSIGNED_BYTE, colortexture+4);
	int bcol[] = { 0, 0, -1, -1 };
	glTexParameteriv (GL_TEXTURE_1D, GL_TEXTURE_BORDER_COLOR, bcol);
      }

    glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, typ);  // DECAL or MODULATE

    glBindTexture (GL_TEXTURE_1D, coltexname);
    if (linear)
      {
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      }
    else
      {
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      }
  }
  */





  void VisualScene :: DrawColorBar (double minval, double maxval, int logscale, bool linear)
  {
    if (!vispar.drawcolorbar) return;

    CreateTexture (8, linear, GL_DECAL);

    if (logscale && maxval <= 0) maxval = 1;
    if (logscale && minval <= 0) minval = 1e-4 * maxval;

    double minx = -1;
    double maxx = 1;
    double miny = 0.75;
    double maxy = 0.8;

    glDisable (GL_LIGHTING);
    glEnable (GL_COLOR_MATERIAL);
    glEnable (GL_TEXTURE_1D);
    glNormal3d (0, 0, 1);
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    
    glDisable (GL_DEPTH_TEST);
    glBegin (GL_QUAD_STRIP);

    for (double x = minx; x <= maxx; x += (maxx - minx) / 50)
      {
	SetOpenGlColor (x, minx, maxx);
	glVertex3d (x, miny, -5);
	glVertex3d (x, maxy, -5);
      }
    glEnd();

    glDisable (GL_TEXTURE_1D);
    
    glEnable (GL_COLOR_MATERIAL);
    GLfloat textcol[3] = { 1 - backcolor, 1 - backcolor, 1 - backcolor };
    glColor3fv (textcol);
    
    glPushAttrib (GL_LIST_BIT);
    // glListBase (fontbase);

    char buf[20];
    for (int i = 0; i <= 4; i++)
      {
	double x = minx + i * (maxx-minx) / 4;
	glRasterPos3d (x, 0.7,-5);
      
	double val;
	if (logscale)
	  val = minval * pow (maxval / minval, i / 4.0);
	else
	  val = minval + i * (maxval-minval) / 4;

	sprintf (buf, "%8.3e", val);
	// glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
	MyOpenGLText (buf);
      }

    glPopAttrib ();
    glEnable (GL_DEPTH_TEST);
  }


  void VisualScene :: DrawCoordinateCross ()
  {
    if (!vispar.drawcoordinatecross) return;

    glDisable (GL_DEPTH_TEST);
    glMatrixMode (GL_PROJECTION); 
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode (GL_MODELVIEW); 
    glPushMatrix();
    glLoadIdentity();

    GLint viewport[4];
    glGetIntegerv (GL_VIEWPORT, viewport);

    glTranslatef (-1, -1, 0.0);
    glScalef (40.0 / viewport[2], 40.0 / viewport[3], 1);
    glTranslatef (2.0, 2.0, 0.0);
    glMultMatrixf (rotmat);

    glEnable (GL_COLOR_MATERIAL);
    glDisable (GL_LIGHTING);

    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

    GLfloat textcol[3] = { 1 - backcolor,
			   1 - backcolor,
			   1 - backcolor };
    glColor3fv (textcol);

    glLineWidth (1.0f);

    double len = 1;

    glBegin(GL_LINES);
    glVertex3d (0, 0, 0);
    glVertex3d (len, 0, 0);
    glVertex3d (0.0f, 0.0f, 0.0f);
    glVertex3d (0.0f, len, 0.0f);
    glVertex3d (0.0f, 0.0f, 0.0f);
    glVertex3d (0.0f, 0.0f, len);
    glEnd ();

    glPushAttrib (GL_LIST_BIT);
    // glListBase (fontbase);

    char buf[20];

    glRasterPos3d (len, 0.0f, 0.0f);
    sprintf (buf, "x");
    // glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
    MyOpenGLText (buf);
    glRasterPos3d (0.0f, len, 0.0f);
    sprintf (buf, "y");
    // glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
    MyOpenGLText (buf);
    glRasterPos3d (0.0f, 0.0f, len);
    sprintf (buf, "z");
    // glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
    MyOpenGLText (buf);

    glPopAttrib ();

    glEnable (GL_LIGHTING);

    glMatrixMode (GL_PROJECTION); 
    glPopMatrix();
    glMatrixMode (GL_MODELVIEW); 
    glPopMatrix();
    glEnable (GL_DEPTH_TEST);
  }


  void VisualScene :: DrawNetgenLogo ()
  {
    if (!vispar.drawnetgenlogo) return;

    glDisable (GL_DEPTH_TEST);
    glMatrixMode (GL_PROJECTION); 
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode (GL_MODELVIEW); 
    glPushMatrix();
    glLoadIdentity();

    GLint viewport[4];
    glGetIntegerv (GL_VIEWPORT, viewport);

    glTranslatef (1, -1, 0.0);
    glScalef (40.0 / viewport[2], 40.0 / viewport[3], 1);
    glTranslatef (-7.0, 2.0, 0.0);

    glDisable (GL_CLIP_PLANE0);
    glDisable (GL_LIGHTING);

    glEnable (GL_COLOR_MATERIAL);
    GLfloat textcol[3] = { 1 - backcolor,
			   1 - backcolor,
			   1 - backcolor };
    glColor3fv (textcol);
    glLineWidth (1.0f);

    glPushAttrib (GL_LIST_BIT);
    // glListBase (fontbase);

    char buf[] = "Netgen " PACKAGE_VERSION;

    glRasterPos3d (0.0f, 0.0f, 0.0f);
    // glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
    MyOpenGLText (buf);

    glPopAttrib ();

    glEnable (GL_LIGHTING);
    glMatrixMode (GL_PROJECTION); 
    glPopMatrix();
    glMatrixMode (GL_MODELVIEW); 
    glPopMatrix();
    glEnable (GL_DEPTH_TEST);
  }







#ifdef PARALLELGL
  void VisualScene :: InitParallelGL ()
  {
    static int init = 0;

    if (!init)
      {
	init = 1;

	if (id == 0)
	  {
	    string displname;
	    
	    Display * dpy = glXGetCurrentDisplay();
	    GLXDrawable drawable = glXGetCurrentDrawable();
	    GLXContext ctx = glXGetCurrentContext();
	    GLXContextID xid = glXGetContextIDEXT (ctx);
	    
	    displname = XDisplayName (0);
	    
	    cout << "Init Parallel GL" << endl;
	    cout << "DisplayName = " << displname << endl;
	    cout << "current display = " << dpy << endl;
	    cout << "current drawable = " << drawable << endl;                  
	    cout << "current context = " << ctx << endl;                  
	    
	    cout << "contextid = " << xid << endl;
	    cout << "isdirect = " << glXIsDirect ( dpy, ctx ) << endl;                  
	    cout << "extensionstring = " << glXQueryExtensionsString( dpy, 0 ) << endl;

	    Array<MPI_Request> request(ntasks);
	    MPI_Status status;
	    for ( int dest = 1; dest < ntasks; dest++ )
	      {
		MyMPI_Send ("redraw", dest);
		MyMPI_Send ("init", dest);
		
		MyMPI_Send (displname, dest);
		MyMPI_Send (int (drawable), dest);
		MyMPI_Send (int (xid), dest);

		int hi;
		MPI_Irecv( &hi, 1, MPI_INT, dest, MPI_ANY_TAG, MPI_COMM_WORLD, &request[dest]);
		// MyMPI_IRecv (hi, dest, request[dest]);
	      } 
	    for ( int dest = 1; dest < ntasks; dest++ )
	      {
		MPI_Wait(&request[dest], &status);
	      }
	  }
      }
  }


  void VisualScene :: Broadcast ()
  {
    if (ntasks == 1) return;

    if (id == 0)
      {
	for (int dest = 1; dest < ntasks; dest++)
	  {
	    MyMPI_Send ("redraw", dest);
	    MyMPI_Send ("broadcast", dest);
	  }
      }

    MyMPI_Bcast (selface);

    vssolution.Broadcast ();
  }
#endif 






  /* *********************** Draw 2D Geometry **************** */


  VisualSceneGeometry2d :: VisualSceneGeometry2d ()
    : VisualScene()
  {
    ;
  }

  VisualSceneGeometry2d :: ~VisualSceneGeometry2d ()
  {
    ;
  }



  void VisualSceneGeometry2d :: DrawScene ()
  {
    if (changeval != geometry2d->GetSplines().Size())
      BuildScene();
    changeval = geometry2d->GetSplines().Size();

  
    glClearColor(backcolor, backcolor, backcolor, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    SetLight();

    //  glEnable (GL_LIGHT0);
    glDisable (GL_LIGHTING);
    glPushMatrix();
    glMultMatrixf (transformationmat);

    //  SetClippingPlane ();

    glShadeModel (GL_SMOOTH);
    glEnable (GL_COLOR_MATERIAL);
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
  
    //  float mat_col[] = { 0, 0, 1, 1 };
    //  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);
    glColor3f (0, 0, 1);
  

    Array<Point<2> > points, otherpoints;

    for (int i = 1; i <= geometry2d->GetSplines().Size(); i++)
      {
	geometry2d->GetSplines().Get(i)->GetPoints (20, points);
      
	glBegin (GL_LINE_STRIP);
	for (int j = 0; j < points.Size(); j++)
	  glVertex3f (points[j](0), points[j](1), 0);
	glEnd(); 
      }

    glColor3f (1, 0, 0);

    for (int i = 1; i <= geometry2d->GetSplines().Size(); i++)
      {
	int other = geometry2d->GetSplines().Get(i)->copyfrom;
	if (other != -1)
	  {
	    geometry2d->GetSplines().Get(i)->GetPoints (6, points);
	    geometry2d->GetSplines().Get(other)->GetPoints (6, otherpoints);
	    glBegin (GL_LINES);
	    for (int j = 1; j < 5; j++)
	      {
		glVertex3f (points[j](0), points[j](1), 0);
		glVertex3f (otherpoints[j](0), otherpoints[j](1), 0);
	      }
	    glEnd ();
	  }
      }



    glPopMatrix();
  
    DrawCoordinateCross ();
    DrawNetgenLogo ();

    glFinish();  
  }


  void VisualSceneGeometry2d :: BuildScene (int zoomall)
  {
    Box<2> bbox;

    geometry2d->GetBoundingBox (bbox);
  
    Point<2> c = Center (bbox.PMin(), bbox.PMax());

    center = Point3d (c(0), c(1), 0);
    rad = Dist (bbox.PMin(), bbox.PMax()) / 2;

    CalcTransformationMatrices();
  }










  /* *********************** Draw STL Geometry **************** */


  VisualSceneSTLGeometry :: VisualSceneSTLGeometry ()
    : VisualScene()
  {
    ;
  }

  VisualSceneSTLGeometry :: ~VisualSceneSTLGeometry ()
  {
    ;
  }

  void VisualSceneSTLGeometry :: DrawScene ()
  {
    if (changeval != stlgeometry->GetNT())
      BuildScene();

    changeval = stlgeometry->GetNT();


    glClearColor(backcolor, backcolor, backcolor, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    SetLight();


    glPushMatrix();
    glMultMatrixf (transformationmat);




    glShadeModel (GL_SMOOTH);
    glDisable (GL_COLOR_MATERIAL);
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    double shine = vispar.shininess;
    // double transp = vispar.transp;

    glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, shine);
    glLogicOp (GL_COPY);


    float mat_col[] = { 0.2f, 0.2f, 0.8f, 1.0f};
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);

    glPolygonOffset (1, 1);
    glEnable (GL_POLYGON_OFFSET_FILL);

    glCallList (trilists.Get(1));

    glDisable (GL_POLYGON_OFFSET_FILL);


    int showtrias = vispar.showstltrias;

    if (showtrias)
      {
	float mat_coll[] = { 0.2f, 0.2f, 0.2f, 1.0f };
	glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_coll);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      
	glCallList (trilists.Get(1));
      }

    /*

    glBegin (GL_TRIANGLES);
    for (j = 1; j <= stlgeometry -> GetNT(); j++)
    {
    const STLTriangle & tria = stlgeometry -> GetTriangle(j);
    glNormal3f (tria.normal.X(),
    tria.normal.Y(),
    tria.normal.Z());
		  
    for (k = 0; k < 3; k++)
    {
    glVertex3f (tria.pts[k].X(),
    tria.pts[k].Y(),
    tria.pts[k].Z());
    }
    }    
    glEnd ();
    */  



 
    glPopMatrix();
    glFinish();  
  }


  void VisualSceneSTLGeometry :: BuildScene (int zoomall)
  {
    //  cout << "rebuild stl geometry scene" << endl;

    center = stlgeometry -> GetBoundingBox().Center();
    rad = stlgeometry -> GetBoundingBox().Diam() / 2;


    CalcTransformationMatrices();

    for (int i = 1; i <= trilists.Size(); i++)
      glDeleteLists (trilists.Elem(i), 1);
    trilists.SetSize(0);


    trilists.Append (glGenLists (1));
    glNewList (trilists.Last(), GL_COMPILE);

    glEnable (GL_NORMALIZE);

    glBegin (GL_TRIANGLES);
    for (int j = 1; j <= stlgeometry -> GetNT(); j++)
      {
	const Vec3d & n = stlgeometry->GetTriangle(j).Normal();
	glNormal3f (n.X(), n.Y(), n.Z());
      
	for (int k = 1; k <= 3; k++)
	  {
	    const Point3d & p = 
	      stlgeometry->GetPoint (stlgeometry -> GetTriangle(j).PNum(k));
	    glVertex3f (p.X(),p.Y(), p.Z());
	  }
      }    
    glEnd ();
      
    glEndList ();
  }












  VisualSceneSpecPoints :: VisualSceneSpecPoints ()
    : VisualScene()
  {
    ;
  }

  VisualSceneSpecPoints :: ~VisualSceneSpecPoints ()
  {
    ;
  }


  void VisualSceneSpecPoints :: DrawScene ()
  {
    if (!mesh) 
      {
	VisualScene::DrawScene();
	return;
      }

    if (changeval != specpoints.Size())
      BuildScene();
    changeval = specpoints.Size();



    glClearColor(backcolor, backcolor, backcolor, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable (GL_COLOR_MATERIAL);
    glColor3f (1.0f, 1.0f, 1.0f);
    glLineWidth (1.0f);

    glPushMatrix();
    glMultMatrixf (transformationmat);

    //  glEnable (GL_COLOR);
    //  glDisable (GL_COLOR_MATERIAL);
    if (vispar.drawedtangents)
      {
	glColor3d (1, 0, 0);
	glBegin (GL_LINES);
	for (int i = 1; i <= specpoints.Size(); i++)
	  {
	    const Point3d p1 = specpoints.Get(i).p;
	    const Point3d p2 = specpoints.Get(i).p + len * specpoints.Get(i).v;
	    glVertex3d (p1.X(), p1.Y(), p1.Z());
	    glVertex3d (p2.X(), p2.Y(), p2.Z());
	  }
	glEnd();
      }

    if (vispar.drawededges)
      {
	glColor3d (1, 0, 0);
	glBegin (GL_LINES);
	for (int i = 1; i <= mesh->GetNSeg(); i++)
	  {
	    const Segment & seg = mesh -> LineSegment (i);
	    glVertex3dv ( (*mesh)[seg[0]] );
            glVertex3dv ( (*mesh)[seg[1]] );
	    // glVertex3dv ( &(*mesh)[seg[0]].X() );
	    // glVertex3dv ( &(*mesh)[seg[1]].X() );
	  }
	glEnd();
      }

    glColor3d (1, 0, 0);
    glBegin (GL_LINES);
    for (int i = 0; i < boxes.Size(); i++)
      {
	glVertex3dv ( boxes[i].PMin() );
	glVertex3dv ( boxes[i].PMax() );
      }
    glEnd();



    if (vispar.drawededgenrs)
      {
	glEnable (GL_COLOR_MATERIAL);
	GLfloat textcol[3] = { 1 - backcolor,
			       1 - backcolor,
			       1 - backcolor };
	glColor3fv (textcol);
	glNormal3d (0, 0, 1);
	glPushAttrib (GL_LIST_BIT);
	// glListBase (fontbase);

	char buf[20];
	for (int i = 1; i <= mesh->GetNSeg(); i++)
	  {
	    const Segment & seg = mesh -> LineSegment (i);
	    const Point3d p1 = mesh -> Point (seg[0]);
	    const Point3d p2 = mesh -> Point (seg[1]);

	    const Point3d p = Center (p1, p2);
	    glRasterPos3d (p.X(), p.Y(), p.Z());
	  
	    sprintf (buf, "%d", seg.edgenr);
	    // glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
	    MyOpenGLText (buf);
	  }
      
	glPopAttrib ();
	glDisable (GL_COLOR_MATERIAL);
      }


    if (vispar.drawedpoints)
      {

	glColor3d (0, 0, 1);
	/*
	  glPointSize( 3.0 );

	float range[2];
	glGetFloatv(GL_POINT_SIZE_RANGE, &range[0]);
	cout << "max ptsize = " << range[0] << "-" << range[1] << endl;
      

	glBegin( GL_POINTS );
	for (int i = 1; i <= mesh -> GetNP(); i++)
	  {
	    const Point3d & p = mesh -> Point(i);
	    if (i % 2)
	      glVertex3f( p.X(), p.Y(), p.Z());
	  }
	glEnd();
	*/

	static GLubyte knoedel[] = 
	  {
	    0xfe, 0xfe, 0xfe, 0xfe, 0xfe, 0xfe, 0xfe,
	  };
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	glDisable (GL_COLOR_MATERIAL);
	glDisable (GL_LIGHTING);
	glDisable (GL_CLIP_PLANE0);
      
	for (int i = 1; i <= mesh -> GetNP(); i++)
	  {
	    const Point3d & p = mesh -> Point(i);
	    glRasterPos3d (p.X(), p.Y(), p.Z());
	    glBitmap (7, 7, 3, 3, 0, 0, &knoedel[0]);
	  }
      }

    if (vispar.drawedpointnrs)
      {
	glEnable (GL_COLOR_MATERIAL);
	GLfloat textcol[3] = { 1 - backcolor,
			       1 - backcolor,
			       1 - backcolor };
	glColor3fv (textcol);
	glNormal3d (0, 0, 1);
	glPushAttrib (GL_LIST_BIT);
	// glListBase (fontbase);
      
	char buf[20];
	for (int i = 1; i <= mesh->GetNP(); i++)
	  {
	    const Point3d & p = mesh->Point(i);
	    glRasterPos3d (p.X(), p.Y(), p.Z());
	  
	    sprintf (buf, "%d", i);
	    // glCallLists (GLsizei(strlen (buf)), GL_UNSIGNED_BYTE, buf);
	    MyOpenGLText (buf);
	  }
      
	glPopAttrib ();
	glDisable (GL_COLOR_MATERIAL);
      }




    glPopMatrix();

    if (vispar.drawcoordinatecross)
      DrawCoordinateCross ();
    DrawNetgenLogo ();

    glFinish();  
  }


  void VisualSceneSpecPoints :: BuildScene (int zoomall)
  {
    if (!mesh) 
      {
	VisualScene::BuildScene(zoomall);
	return;
      }
  
    Box3d box;
  
    if (mesh->GetNSeg())
      {
	box.SetPoint (mesh->Point (mesh->LineSegment(1)[0]));
	for (int i = 1; i <= mesh->GetNSeg(); i++)
	  {
	    box.AddPoint (mesh->Point (mesh->LineSegment(i)[0]));
	    box.AddPoint (mesh->Point (mesh->LineSegment(i)[1]));
	  }
      }
    else if (specpoints.Size() >= 2)
      {
	box.SetPoint (specpoints.Get(1).p);
	for (int i = 2; i <= specpoints.Size(); i++)
	  box.AddPoint (specpoints.Get(i).p);
      }
    else
      {
	box = Box3d (Point3d (0,0,0), Point3d (1,1,1));
      }
  
    if (zoomall == 2 && ((vispar.centerpoint >= 1 && vispar.centerpoint <= mesh->GetNP()) ||
			 vispar.use_center_coords))
      {
	if (vispar.use_center_coords)
	  {
	    center.X() = vispar.centerx; center.Y() = vispar.centery; center.Z() = vispar.centerz; 
	  }
	else
	  center = mesh->Point (vispar.centerpoint);
      }
    else
      center = Center (box.PMin(), box.PMax());
        

    rad = 0.5 * Dist (box.PMin(), box.PMax());
  
  
    CalcTransformationMatrices();
  }

}
