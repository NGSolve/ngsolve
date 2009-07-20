#ifndef NOTCL

#include <mystdlib.h>

#include <meshing.hpp>

// #include "incvis.hpp"


#include <visual.hpp>


namespace netgen
{
  // #include "meshdoc.hpp"


MeshDoctorParameters meshdoctor;
VisualSceneMeshDoctor vsmeshdoc;

extern AutoPtr<Mesh> mesh;

  int Ng_MeshDoctor (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
{
  cout << "Mesh Doctor:" << endl;
  int i;
  for (i = 0; i < argc; i++)
    cout << argv[i] << " ";
  cout << endl;

  meshdoctor.active = 
      atoi (Tcl_GetVar (interp, "::meshdoctor.active", 0)); 


  if (argc >= 2)
    {
      if (strcmp (argv[1], "markedgedist") == 0)
	{
	  vsmeshdoc.SetMarkEdgeDist (atoi (argv[2]));
	}

      if (strcmp (argv[1], "deletemarkedsegments") == 0)
	{
	  for (i = 1; i <= mesh->GetNSeg(); i++)
	    if (vsmeshdoc.IsSegmentMarked (i))
	      mesh->DeleteSegment (i);

	  //	  for (i = 1; i <= mesh->GetNSE(); i++)
	  //	    mesh->SurfaceElement(i).SetIndex (1);
	  mesh->Compress();
	}
    }


  vsmeshdoc.UpdateTables ();
  vsmeshdoc.BuildScene();
  return TCL_OK;
}





VisualSceneMeshDoctor :: VisualSceneMeshDoctor ()
  : VisualScene()
{
  filledlist = 0;
  outlinelist = 0;
  edgelist = 0;
  selelement = 0;
  locpi = 1;
  selpoint = 0;
  selpoint2 = 0;
  markedgedist = 1;

  UpdateTables ();
}

VisualSceneMeshDoctor :: ~VisualSceneMeshDoctor ()
{
  ;
}

void VisualSceneMeshDoctor :: DrawScene ()
{
  if (!mesh) return;

  int hchval = mesh->GetNP() + mesh->GetNE() + mesh->GetNSE();
  if (changeval != hchval)
    {
      changeval = hchval;
      BuildScene();
    }


  glClearColor(backcolor, backcolor, backcolor, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnable (GL_COLOR_MATERIAL);
  glColor3f (1.0f, 1.0f, 1.0f);
  glLineWidth (1.0f);

  SetLight();

  glPushMatrix();
  glMultMatrixf (transformationmat);
  
  glInitNames ();
  glPushName (0);
  
  glPolygonOffset (1, 1);
  glEnable (GL_POLYGON_OFFSET_FILL);

  SetClippingPlane ();

  if (vispar.drawfilledtrigs)
    glCallList (filledlist);

  glDisable (GL_POLYGON_OFFSET_FILL);
  
  if (vispar.drawoutline)
    glCallList (outlinelist);
  
  glPolygonOffset (-1, -1);
  glEnable (GL_POLYGON_OFFSET_LINE);

  if (vispar.drawedges)
    glCallList (edgelist);
  

  glDisable (GL_POLYGON_OFFSET_LINE);

  
  
  glPopName();

  if (selpoint > 0 && selpoint <= mesh->GetNP())
    {
      GLfloat matcolblue[] = { 0, 0, 1, 1 };

      glPointSize (10);
      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matcolblue);
      glBegin (GL_POINTS);
      
      const Point3d p = mesh->Point(selpoint);
      glVertex3f (p.X(), p.Y(), p.Z());
      glEnd();
    }

  glDisable(GL_CLIP_PLANE0);


  glPopMatrix();
  glFinish();  
}




void VisualSceneMeshDoctor :: BuildScene (int zoomall)
{
  int i, j;
 
  
  if (zoomall)
    {
      Point3d pmin, pmax;
      mesh->GetBox (pmin, pmax, -1);

      if (vispar.centerpoint)
	center = mesh->Point (vispar.centerpoint);
      else
	center = Center (pmin, pmax);
  
      rad = 0.5 * Dist (pmin, pmax);

      glEnable (GL_NORMALIZE);
  
      CalcTransformationMatrices();
    }




  if (filledlist)
    {
      glDeleteLists (filledlist, 1);
      glDeleteLists (outlinelist, 1);
      glDeleteLists (edgelist, 1);
    }

  
  filledlist = glGenLists (1);
  glNewList (filledlist, GL_COMPILE);

  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
  
  glLineWidth (1.0f);
  
  glDisable (GL_COLOR_MATERIAL);
    
  for (i = 1; i <= mesh->GetNSE(); i++)
    {
      glLoadName (i);

      // copy to be thread-safe
      Element2d el = mesh->SurfaceElement (i);

      int drawel = 1;
      for (j = 1; j <= el.GetNP(); j++)
	{
	  if (!el.PNum(j))
	    drawel = 0;
	}

      if (!drawel)
	continue;

      GLfloat matcol[] = { 0, 1, 0, 1 };
      GLfloat matcolsel[] = { 1, 0, 0, 1 };

      if (i == selelement)
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, matcolsel);
      else
	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, matcol);

      if (el.GetNP() == 3)
	{
	  glBegin (GL_TRIANGLES);
	  
	  const Point3d & lp1 = mesh->Point (el.PNum(1));
	  const Point3d & lp2 = mesh->Point (el.PNum(2));
	  const Point3d & lp3 = mesh->Point (el.PNum(3));
	  Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
	  n /= (n.Length()+1e-12);
	  glNormal3d (n.X(), n.Y(), n.Z());

	  if (!vispar.colormeshsize)
	    {
	      glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
	      glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
	      glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	    }
	  else
	    {
	      double h1 = mesh->GetH (lp1);
	      double h2 = mesh->GetH (lp2);
	      double h3 = mesh->GetH (lp3);
	      
	      SetOpenGlColor  (h1, 0.1, 10);
	      glVertex3d (lp1.X(), lp1.Y(), lp1.Z());

	      SetOpenGlColor  (h2, 0.1, 10);
	      glVertex3d (lp2.X(), lp2.Y(), lp2.Z());

	      SetOpenGlColor  (h3, 0.1, 10);
	      glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	    }	    
	  glEnd();
	}
      else if (el.GetNP() == 4)
	{
	  glBegin (GL_QUADS);
	  
	  const Point3d & lp1 = mesh->Point (el.PNum(1));
	  const Point3d & lp2 = mesh->Point (el.PNum(2));
	  const Point3d & lp3 = mesh->Point (el.PNum(4));
	  const Point3d & lp4 = mesh->Point (el.PNum(3));
	  Vec3d n = Cross (Vec3d (lp1, lp2), 
			   Vec3d (lp1, Center (lp3, lp4)));
	  n /= (n.Length()+1e-12);
	  glNormal3d (n.X(), n.Y(), n.Z()); 
	  glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
	  glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
	  glVertex3d (lp4.X(), lp4.Y(), lp4.Z());
	  glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	  glEnd();
	}
      else if (el.GetNP() == 6)
	{
	  glBegin (GL_TRIANGLES);
	  static int trigs[4][3] = {
	    { 1, 6, 5 },
	    { 2, 4, 6 },
	    { 3, 5, 4 },
	    { 4, 5, 6 } };

	  for (j = 0; j < 4; j++)
	    {
	      const Point3d & lp1 = mesh->Point (el.PNum(trigs[j][0]));
	      const Point3d & lp2 = mesh->Point (el.PNum(trigs[j][1]));
	      const Point3d & lp3 = mesh->Point (el.PNum(trigs[j][2]));
	      Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
	      n /= (n.Length() + 1e-12);
	      glNormal3d (n.X(), n.Y(), n.Z());
	      glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
	      glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
	      glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	    }
	  glEnd();
	}
    }
  glLoadName (0);
  
  glEndList ();

  
  
  outlinelist = glGenLists (1);
  glNewList (outlinelist, GL_COMPILE);

  glLineWidth (1.0f);
  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);

  glColor3f (0.0f, 0.0f, 0.0f);
  glEnable (GL_COLOR_MATERIAL);
  
  for (i = 1; i <= mesh->GetNSE(); i++)
    {
      Element2d el = mesh->SurfaceElement(i);

      int drawel = 1;
      for (j = 1; j <= el.GetNP(); j++)
	{
	  if (!el.PNum(j))
	    drawel = 0;
	}

      if (!drawel)
	continue;


      if (el.GetNP() == 3)
	{
	  glBegin (GL_TRIANGLES);
	  
	  const Point3d & lp1 = mesh->Point (el.PNum(1));
	  const Point3d & lp2 = mesh->Point (el.PNum(2));
	  const Point3d & lp3 = mesh->Point (el.PNum(3));
	  Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
	  n /= (n.Length() + 1e-12);
	  glNormal3d (n.X(), n.Y(), n.Z());
	  glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
	  glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
	  glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	  glEnd();
	}
      else if (el.GetNP() == 4)
	{
	  glBegin (GL_QUADS);
	  
	  const Point3d & lp1 = mesh->Point (el.PNum(1));
	  const Point3d & lp2 = mesh->Point (el.PNum(2));
	  const Point3d & lp3 = mesh->Point (el.PNum(4));
	  const Point3d & lp4 = mesh->Point (el.PNum(3));
	  Vec3d n = Cross (Vec3d (lp1, lp2), 
			   Vec3d (lp1, Center (lp3, lp4)));
	  n /= (n.Length() + 1e-12);
	  glNormal3d (n.X(), n.Y(), n.Z());
	  glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
	  glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
	  glVertex3d (lp4.X(), lp4.Y(), lp4.Z());
	  glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	  glEnd();
	}
      else if (el.GetNP() == 6)
	{
	  glBegin (GL_LINES);
	  
	  const Point3d & lp1 = mesh->Point (el.PNum(1));
	  const Point3d & lp2 = mesh->Point (el.PNum(2));
	  const Point3d & lp3 = mesh->Point (el.PNum(3));
	  const Point3d & lp4 = mesh->Point (el.PNum(4));
	  const Point3d & lp5 = mesh->Point (el.PNum(5));
	  const Point3d & lp6 = mesh->Point (el.PNum(6));

	  Vec3d n = Cross (Vec3d (lp1, lp2), Vec3d (lp1, lp3));
	  n /= (n.Length()+1e-12);
	  glNormal3d (n.X(), n.Y(), n.Z());

	  glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
	  glVertex3d (lp6.X(), lp6.Y(), lp6.Z());
	  glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
	  glVertex3d (lp6.X(), lp6.Y(), lp6.Z());

	  glVertex3d (lp1.X(), lp1.Y(), lp1.Z());
	  glVertex3d (lp5.X(), lp5.Y(), lp5.Z());
	  glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	  glVertex3d (lp5.X(), lp5.Y(), lp5.Z());

	  glVertex3d (lp2.X(), lp2.Y(), lp2.Z());
	  glVertex3d (lp4.X(), lp4.Y(), lp4.Z());
	  glVertex3d (lp3.X(), lp3.Y(), lp3.Z());
	  glVertex3d (lp4.X(), lp4.Y(), lp4.Z());
	  glEnd();
	}
    }
  glLoadName (0);  
  glEndList ();





  edgelist = glGenLists (1);
  glNewList (edgelist, GL_COMPILE);

  glDisable (GL_COLOR_MATERIAL);

  GLfloat matcoledge[] = { 0, 0, 1, 1 };
  GLfloat matcolseledge[] = { 1, 0, 1, 1 };

  glLineWidth (2.0f);

  for (i = 1; i <= mesh->GetNSeg(); i++)
    {
      const Segment & seg = mesh->LineSegment(i);
      const Point3d & p1 = mesh->Point(seg[0]);
      const Point3d & p2 = mesh->Point(seg[1]);

      if (edgedist.Get(seg[0]) <= markedgedist &&
	  edgedist.Get(seg[1]) <= markedgedist)
	{
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, 
			matcolseledge);
	  glLineWidth (4.0f);
	}
      else
	{
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, 
			matcoledge);
	  glLineWidth (2.0f);
	}
      glBegin (GL_LINES);
      glVertex3f (p1.X(), p1.Y(), p1.Z());
      glVertex3f (p2.X(), p2.Y(), p2.Z());
      glEnd(); 
    }

  glLineWidth (1.0f);
  glEndList ();
}




void VisualSceneMeshDoctor :: MouseDblClick (int px, int py)
{
  cout << "dblclick: " << px << " - " << py << endl;
  
  int i, hits;

  // select surface triangle by mouse click
  GLuint selbuf[10000];
  glSelectBuffer (10000, selbuf);


  glRenderMode (GL_SELECT);

  GLint viewport[4];
  glGetIntegerv (GL_VIEWPORT, viewport);

  glMatrixMode (GL_PROJECTION); 
  glPushMatrix();

  GLdouble projmat[16];
  glGetDoublev (GL_PROJECTION_MATRIX, projmat);

  glLoadIdentity(); 
  gluPickMatrix (px, viewport[3] - py, 1, 1, viewport); 
  glMultMatrixd (projmat);
  

  glClearColor(backcolor, backcolor, backcolor, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode (GL_MODELVIEW); 

  glPushMatrix();
  glMultMatrixf (transformationmat);

  glInitNames();
  glPushName (1);

  glPolygonOffset (1, 1);
  glEnable (GL_POLYGON_OFFSET_FILL);

  glCallList (filledlist);

  glDisable (GL_POLYGON_OFFSET_FILL);
  
  glPopName();

  glMatrixMode (GL_PROJECTION); 
  glPopMatrix();

  glMatrixMode (GL_MODELVIEW); 
  glPopMatrix();

  glFlush();  

	
  hits = glRenderMode (GL_RENDER);

  cout << "hits = " << hits << endl;

  int minname = 0;
  GLuint mindepth = 0;
  for (i = 0; i < hits; i++)
    {
      int curname = selbuf[4*i+3];
      GLuint curdepth = selbuf[4*i+1];

      if (curname &&
	  (curdepth < mindepth || !minname))
	{
	  mindepth = curdepth;
	  minname = curname;
	}
    }

  cout << "clicked element: " << minname << endl;

  ClickElement (minname);

  BuildScene ();
}




void VisualSceneMeshDoctor :: SetMarkEdgeDist (int dist)
{
  markedgedist = dist;
  BuildScene();
}

void VisualSceneMeshDoctor :: ClickElement (int elnr)
{
  selelement = elnr;

  int oldlocpi = locpi;
  locpi = locpi % 3 + 1;
  
  if (selelement > 0 && selelement <= mesh->GetNSE())
    {
      selpoint = mesh->SurfaceElement(selelement).PNum(locpi);
      selpoint2 = mesh->SurfaceElement(selelement).PNum(oldlocpi);
      cout << "selpts = " << selpoint << ", " << selpoint2 << endl;
    }

  UpdateTables();
}


void VisualSceneMeshDoctor :: UpdateTables ()
{
  if (!mesh) return;

  edgedist.SetSize(mesh->GetNP());
  int i, changed;

  for (i = 1; i <= mesh->GetNP(); i++)
    edgedist.Elem(i) = 10000;

  for (i = 1; i <= mesh->GetNSeg(); i++)
    {
      const Segment & seg = mesh->LineSegment(i);
      if ( (seg[0] == selpoint && seg[1] == selpoint2) ||
           (seg[1] == selpoint && seg[0] == selpoint2) )
	{
	  edgedist.Elem(selpoint) = 1;
	  edgedist.Elem(selpoint2) = 1;
	}
    }

  do
    {
      changed = 0;

      for (i = 1; i <= mesh->GetNSeg(); i++)
	{
	  const Segment & seg = mesh->LineSegment(i);
	  
	  int edist = min2 (edgedist.Get(seg[0]), edgedist.Get(seg[1]));
	  edist++;

	  if (edgedist.Get(seg[0]) > edist)
	    {
	      edgedist.Elem(seg[0]) = edist;
	      changed = 1;
	    }
	  if (edgedist.Get(seg[1]) > edist)
	    {
	      edgedist.Elem(seg[1]) = edist;
	      changed = 1;
	    }
	}	    
    }
  while (changed);
}

int VisualSceneMeshDoctor :: IsSegmentMarked (int segnr) const
{
  const Segment & seg = mesh->LineSegment(segnr);
  return (edgedist.Get(seg[0]) <= markedgedist &&
	  edgedist.Get(seg[1]) <= markedgedist);
}
}


#endif // NOTCL
