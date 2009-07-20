#ifndef NOTCL

#include <mystdlib.h>
#include "incvis.hpp"

#include <myadt.hpp>
#include <meshing.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>

#include <visual.hpp>


namespace netgen
{

/* *********************** Draw Geometry **************** */

  extern Array<Point<3> > project1, project2;


extern AutoPtr<CSGeometry> geometry;


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





}


#endif // NOTCL
