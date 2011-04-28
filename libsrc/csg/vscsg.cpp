#include <mystdlib.h>
#include "incvis.hpp"

#include <myadt.hpp>
#include <meshing.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>

#include <visual.hpp>

#include "vscsg.hpp"

namespace netgen
{



  /* *********************** Draw Geometry **************** */

  extern AutoPtr<Mesh> mesh;
  extern Array<SpecialPoint> specpoints;
  extern Array<Box<3> > boxes;


  extern Array<Point<3> > project1, project2;


  // extern AutoPtr<CSGeometry> geometry;


  VisualSceneGeometry :: VisualSceneGeometry ()
    : VisualScene()
  {
    selsurf = 0;
  }

  VisualSceneGeometry :: ~VisualSceneGeometry ()
  {
    ;
  }

  void VisualSceneGeometry :: SelectSurface (int aselsurf)
  {
    selsurf = aselsurf;
    DrawScene();
  }


  void VisualSceneGeometry :: DrawScene ()
  {
    int i;

    if (changeval != geometry->GetChangeVal())
      BuildScene();
    changeval = geometry->GetChangeVal();

    glClearColor(backcolor, backcolor, backcolor, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
    SetLight();


    glPushMatrix();
    glMultMatrixf (transformationmat);

    SetClippingPlane ();

    glShadeModel (GL_SMOOTH);
    glDisable (GL_COLOR_MATERIAL);
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    /*
      float mat_spec_col[] = { 1, 1, 1, 1 };
      glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_spec_col);
    */

    double shine = vispar.shininess;
    double transp = vispar.transp;


    glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, shine);
    glLogicOp (GL_COPY);
  
    glEnable (GL_NORMALIZE);

    for (i = 0; i < geometry->GetNTopLevelObjects(); i++)
      {
	const TopLevelObject * tlo = geometry -> GetTopLevelObject (i);
	if (tlo->GetVisible() && !tlo->GetTransparent())
	  {
	    float mat_col[] = { tlo->GetRed(), tlo->GetGreen(), tlo->GetBlue(), 1 };
	    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);
	  
	    glCallList (trilists[i]);
	  }
      }


    glPolygonOffset (1, 1);
    glEnable (GL_POLYGON_OFFSET_FILL);

    glLogicOp (GL_NOOP);
    for (i = 0; i < geometry->GetNTopLevelObjects(); i++)
      {
	const TopLevelObject * tlo = geometry -> GetTopLevelObject (i);
	if (tlo->GetVisible() && tlo->GetTransparent())
	  {
	    float mat_col[] = { tlo->GetRed(), tlo->GetGreen(), tlo->GetBlue(), transp };

	    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);
	  
	    glCallList (trilists[i]);
	  }
      }

    glDisable (GL_POLYGON_OFFSET_FILL);

    /*
      cout << "draw " << project1.Size() << " lines " << endl;
      glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      glLineWidth (1.0f);
      glEnable (GL_COLOR_MATERIAL);

      glColor3f (1.0f, 0.0f, 0.0f);

      glBegin (GL_LINES);
      for (int i = 0; i < project1.Size(); i++)
      {
      glVertex3dv (project1[i]);
      glVertex3dv (project2[i]);
      }
      glEnd();
    */


    glPopMatrix();
    glDisable(GL_CLIP_PLANE0);
 


    /*
      glFlush();
  
      int err;
      do
      {
      err = glGetError();
      // cout << "glerr,1 = " << err << endl;
      }
      while (err != GL_NO_ERROR);


      // CreateTexture (0, 1, GL_DECAL);
      CreateTexture (0, 1, GL_MODULATE);
      glEnable (GL_TEXTURE_1D);

      float mat_col[] = { 1.0, 1.0, 1.0 }; 
      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);

      glDisable (GL_BLEND);
      glDisable (GL_COLOR_MATERIAL);
      glEnable (GL_NORMALIZE);

      if (geometry->GetNTopLevelObjects())
      {
      cout << "call list" << endl;
      glCallList (trilists[0]);
      }

      glColor3d (1.0, 1.0, 1.0);
  
      glBegin (GL_TRIANGLES);
      glNormal3f (0, 0, 1);
      SetOpenGlColor  (-1.0, 0, 1, 0);
      glVertex3f (0.0, 0.0, 0.0);
      SetOpenGlColor  (0.5, 0, 1, 0);
      glNormal3f (0, 0, 1);
      glVertex3f (1.0, 0.0, 0.0);
      SetOpenGlColor  (2.0, 0, 1, 0);
      glNormal3f (0, 0, 1);
      glVertex3f (0.0, 1.0, 0.0);

      glEnd ();

      cout << "trig drawn" << endl;

      glDisable (GL_TEXTURE_1D);
      glDisable (GL_COLOR_MATERIAL);
      glFlush();

      cout << "glerr,2 = " << glGetError() << endl;
    */



    DrawCoordinateCross ();
    DrawNetgenLogo ();  

    glFinish();  
  }


  void VisualSceneGeometry :: BuildScene (int zoomall)
  {
    Box<3> box;
    int hasp = 0;
    for (int i = 0; i < geometry->GetNTopLevelObjects(); i++)
      {
	const TriangleApproximation & ta =
	  *geometry->GetTriApprox(i);
	if (!&ta) continue;

	for (int j = 0; j < ta.GetNP(); j++)      
	  {
	    if (hasp)
	      box.Add (ta.GetPoint(j));
	    else
	      {
		hasp = 1;
		box.Set (ta.GetPoint(j));
	      }
	  }
      }
    if (hasp)
      {
	center = box.Center();
	rad = box.Diam() / 2;
      }
    else
      {
	center = Point3d(0,0,0);
	rad = 1;
      }

    CalcTransformationMatrices();

    for (int i = 0; i < trilists.Size(); i++)
      glDeleteLists (trilists[i], 1);
    trilists.SetSize(0);

    for (int i = 0; i < geometry->GetNTopLevelObjects(); i++)
      {
	trilists.Append (glGenLists (1));
	glNewList (trilists.Last(), GL_COMPILE);

	glEnable (GL_NORMALIZE);
	const TriangleApproximation & ta =
	  *geometry->GetTriApprox(i);
	if (&ta) 
	  {
	    glBegin (GL_TRIANGLES);
	    for (int j = 0; j < ta.GetNT(); j++)
	      {
		for (int k = 0; k < 3; k++)
		  {
		    int pi = ta.GetTriangle(j)[k];
		    glNormal3dv (ta.GetNormal (pi));
		    glVertex3dv (ta.GetPoint(pi));
		  }
	      }
	    glEnd ();
	  }
	glEndList ();
      }

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
    int edges[12][2] = 
      { { 0, 1 },
	{ 2, 3 },
	{ 4, 5 },
	{ 6, 7 },
	{ 0, 2 },
	{ 1, 3 },
	{ 4, 6 },
	{ 5, 7 },
	{ 0, 4 },
	{ 1, 5 },
	{ 2, 6 },
	{ 3, 7 } };
    for (int i = 0; i < boxes.Size(); i++)
      {
	for (int j = 0; j < 12; j++)
	  {
	    glVertex3dv ( boxes[i].GetPointNr(edges[j][0]) );
	    glVertex3dv ( boxes[i].GetPointNr(edges[j][1]) );
	  }
	/*
	glVertex3dv ( boxes[i].PMin() );
	glVertex3dv ( boxes[i].PMax() );
	*/
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


