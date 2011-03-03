#include <meshing.hpp>
#include <geometry2d.hpp>
#include <visual.hpp>

#include "vsgeom2d.hpp"

namespace netgen
{


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
	  glVertex3d (points[j](0), points[j](1), 0);
	glEnd(); 
      }

    glColor3f (1, 0, 0);

    for (int i = 1; i <= geometry2d->GetSplines().Size(); i++)
      {
	int other = geometry2d->GetSpline(i-1).copyfrom;
	if (other != -1)
	  {
	    geometry2d->GetSplines().Get(i)->GetPoints (6, points);
	    geometry2d->GetSplines().Get(other)->GetPoints (6, otherpoints);
	    glBegin (GL_LINES);
	    for (int j = 1; j < 5; j++)
	      {
		glVertex3d (points[j](0), points[j](1), 0);
		glVertex3d (otherpoints[j](0), otherpoints[j](1), 0);
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
}
