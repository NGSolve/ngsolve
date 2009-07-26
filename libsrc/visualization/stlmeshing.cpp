#include <mystdlib.h>
#include <myadt.hpp>

#include <linalg.hpp>
#include <stlgeom.hpp>

#include <meshing.hpp>
#ifndef NOTCL
#include <visual.hpp>
#endif

namespace netgen
{

/*
//mmm
#include "stlgeom/modeller.hpp"
*/

/* *********************** Draw STL Geometry **************** */

extern STLGeometry * stlgeometry;
extern AutoPtr<Mesh> mesh;


#ifdef OPENGL

// #include "../../ngtcltk/mvdraw.hpp"


VisualSceneSTLMeshing :: VisualSceneSTLMeshing ()
  : VisualScene()
{
  selecttrig = 0;
  nodeofseltrig = 1;
  stlgeometry->SetSelectTrig(selecttrig);
  stlgeometry->SetNodeOfSelTrig(nodeofseltrig);
}

VisualSceneSTLMeshing :: ~VisualSceneSTLMeshing ()
{
  ;
}

void VisualSceneSTLMeshing :: DrawScene ()
{
  int i, j, k;

  if (changeval != stlgeometry->GetNT())
    BuildScene();
  changeval = stlgeometry->GetNT();

  int colormeshsize = vispar.colormeshsize;
  
  double hmin = 0.0, hmax = 1.0;

  if (colormeshsize)
    {
      hmax = -1E50;
      hmin = +1E50;
      double ms;

      for (i = 1; i <= stlgeometry->GetNP(); i++)
	{
	  ms = mesh->GetH (stlgeometry->GetPoint(i));
	  hmin = min2(hmin,ms);
	  hmax = max2(hmax,ms);
	}

      //hmax = mparam.maxh;
      //hmin = mesh->GetMinH (stlgeometry->GetBoundingBox().PMin(),
      //			    stlgeometry->GetBoundingBox().PMax());
  
      if (hmin == 0) hmin = 0.1 * hmax;
      //hmax *= 1.1;
    }
  


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

  float mat_spec_col[] = { 1, 1, 1, 1 };
  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_spec_col);

  double shine = vispar.shininess;
  // double transp = vispar.transp;

  glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, shine);
  glLogicOp (GL_COPY);

  float mat_colred[]    = { 0.9f, 0.0f, 0.0f, 1.0f };
  float mat_colgreen[]  = { 0.0f, 0.9f, 0.0f, 1.0f };
  float mat_colblue[]   = { 0.1f, 0.1f, 1.0f, 1.0f };

  float mat_colbluegreen[] = { 0.1f, 0.5f, 0.9f, 1.0f };
  // float mat_colpink[]      = { 1.0f, 0.1f, 0.5f, 1.0f };
  float mat_colviolet[]    = { 1.0f, 0.1f, 1.0f, 1.0f };
  float mat_colbrown[]     = { 0.8f, 0.6f, 0.1f, 1.0f };
  // float mat_colorange[]    = { 0.9f, 0.7f, 0.1f, 1.0f };
  // float mat_colturquis[]   = { 0.0f, 1.0f, 0.8f, 1.0f };

  float mat_colgrey[] = { 0.3f, 0.3f, 0.3f, 1.0f };

  float mat_collred[]   = { 1.0f, 0.5f, 0.5f, 1.0f };
  float mat_collgreen[] = { 0.2f, 1.9f, 0.2f, 1.0f };
  float mat_collbrown[] = { 1.0f, 0.8f, 0.3f, 1.0f };

  float mat_collgrey[] = { 0.8f, 0.8f, 0.8f, 1.0f };
  // float mat_colmgrey[] = { 0.4f, 0.4f, 0.4f, 1.0f };

  float mat_colstlbody[] = { 0.0f, 0.0f, 0.8f, 1.0f };
  float mat_colseltrig[] = { 0.7f, 0.7f, 0.3f, 1.0f };
  float mat_colseledge[] = { 0.7f, 0.7f, 1.0f, 1.0f };

  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colblue);

  float pgoff = 0.5f;

  glPolygonOffset (pgoff*1, pgoff*1);
  glEnable (GL_POLYGON_OFFSET_FILL);

  glEnable (GL_NORMALIZE);

  /*
  {
    //mmm
    //test modeller
    Modeller model;
    
    //MoZylinder z1(Point3d(0,0,0),Vec3d(100,0,0),20,0.01);
    //model.Add(&z1);
    //MoZylinder z2(Point3d(50,50,0),Vec3d(0,-100,0),20,0.01);
    //model.Add(&z2);
    
    MoZylinder z1(Point3d(0,0,0),Vec3d(100,0,0),20,0.01);
    MoZylinder z2(Point3d(50,50,0),Vec3d(0,-100,0),20,0.01);
    MoCombine cb1(&z1,&z2);
    model.Add(&cb1);
    
    Array<MoTriangle> trigs;
    model.GetTriangles(trigs);
    int i, k;
    glBegin (GL_TRIANGLES);
    for (i = 1; i <= trigs.Size(); i++)
      {
	const MoTriangle & tria = trigs.Get(i);
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
    

  }

*/



  
  if (!stlgeometry->trigsconverted)
    {
      glBegin (GL_TRIANGLES);
      for (j = 1; j <= stlgeometry -> GetNT(); j++)
	{
	  /*
	  if (j % 10 == seltria)
	    glMaterialfv (GL_FRONT_AND_BACK, 
			  GL_AMBIENT_AND_DIFFUSE, mat_colred);
	  */

	  const Vec3d & n = stlgeometry->GetTriangle(j).Normal();
	  glNormal3f (n.X(), n.Y(), n.Z());
	  /*
	  const STLReadTriangle & tria = stlgeometry -> GetReadTriangle(j);
	  glNormal3f (tria.normal.X(),
		      tria.normal.Y(),
		      tria.normal.Z());
	  */

	  
	  for (k = 1; k <= 3; k++)
	    {
	      const Point3d & tp = stlgeometry->GetPoint(stlgeometry->GetTriangle(j).PNum(k));
	      glVertex3f (tp.X(), tp.Y(), tp.Z());

	    }
	  /*
	  if (j%10 == seltria)
	    glMaterialfv (GL_FRONT_AND_BACK, 
			  GL_AMBIENT_AND_DIFFUSE, mat_colblue);
	  */
	}    
      glEnd ();
  
      glDisable (GL_POLYGON_OFFSET_FILL);

      int showtrias = vispar.stlshowtrias;

      if (showtrias)
	{
	  float mat_coll[] = { 0.2f, 0.2f, 0.2f, 1.f };
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_coll);
	  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      
	  glEnable (GL_NORMALIZE);
      
	  glBegin (GL_TRIANGLES);
	  for (j = 1; j <= stlgeometry -> GetNT(); j++)
	    {
	      const Vec3d & n = stlgeometry->GetTriangle(j).Normal();
	      glNormal3f (n.X(), n.Y(), n.Z());
	      /*
	      const STLReadTriangle & tria = stlgeometry -> GetReadTriangle(j);
	      glNormal3f (tria.normal.X(),
			  tria.normal.Y(),
			  tria.normal.Z());
	      */  

	      for (k = 1; k <= 3; k++)
		{
		  const Point3d & tp = 
		    stlgeometry->GetPoint(stlgeometry->GetTriangle(j).PNum(k));
		  glVertex3f (tp.X(), tp.Y(), tp.Z());
		  
		}
	      
	      /*
	      for (k = 0; k < 3; k++)
		{
		  glVertex3f (tria.pts[k].X(),
			      tria.pts[k].Y(),
			      tria.pts[k].Z());
		}
	      */
	    }    
	  glEnd ();
	}
    }
  else
    {
      int showfilledtrias = vispar.stlshowfilledtrias;

      //(*mycout) << "in " << showfilledtrias << ", NT=" << stlgeometry -> GetNT() << endl;

      int chartnumber;
      if (vispar.stlshowmarktrias)
	chartnumber = vispar.stlchartnumber + vispar.stlchartnumberoffset;
      else
	chartnumber = stlgeometry->GetMeshChartNr();

      if (showfilledtrias)
	{
	  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	  if (colormeshsize)
	    glEnable (GL_COLOR_MATERIAL);
	  
	  glPolygonOffset (pgoff*4, pgoff*4);
	  glEnable (GL_POLYGON_OFFSET_FILL);
	  glEnable (GL_NORMALIZE);


	  glBegin (GL_TRIANGLES);

	  int selt = stlgeometry -> GetSelectTrig();
	  if (stldoctor.selectmode != 0) 
	    {selt = 0; } //do not show selected triangle!!!!

	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colstlbody);

	  for (j = 1; j <= stlgeometry -> GetNT(); j++)
	    {
	      if (stldoctor.showvicinity && !stlgeometry->Vicinity(j)) {continue;}

	      if (j == selt)
		{
		  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colseltrig);
		}
	      else if (j == selt+1)
		{
		  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colstlbody);
		}
	      
	      const STLTriangle& st = stlgeometry -> GetTriangle(j);

	      const Vec3d & n = stlgeometry->GetTriangle(j).Normal();
	      glNormal3f (n.X(), n.Y(), n.Z());
	  
	      /*
	      const STLReadTriangle& tria = stlgeometry -> GetReadTriangle(j);
	      glNormal3f (tria.normal.X(),
			  tria.normal.Y(),
			  tria.normal.Z());
	      */
	      for (k = 0; k < 3; k++)
		{
		  const Point3d & p = stlgeometry->GetPoint(st[k]);
		  if (colormeshsize)
		    {
		      SetOpenGlColor (mesh->GetH (p), hmin, hmax, 1);
		    }

		  glVertex3f (p.X(), p.Y(), p.Z());
		}
	    } 
   
	  glEnd ();
	}
      
      int foundseltrig = stlgeometry -> GetSelectTrig();
      if (foundseltrig == 0 || foundseltrig > stlgeometry->GetNT() ||
	  (stldoctor.showvicinity && !stlgeometry->Vicinity(foundseltrig)))
	{foundseltrig = 0;}

      if (foundseltrig)
	{

	  glPolygonOffset (pgoff*0, 0);
	  glEnable (GL_POLYGON_OFFSET_FILL);

	  //glDisable (GL_POLYGON_OFFSET_FILL);      
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colseledge);
	  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	  
	  glEnable (GL_NORMALIZE);

	  if (stldoctor.selectmode == 2)
	    {
	      //point
	      const STLTriangle& st = stlgeometry -> GetTriangle(foundseltrig);
	      const Point3d & p1 = stlgeometry->GetPoint(st[0]);
	      const Point3d & p2 = stlgeometry->GetPoint(st[1]);
	      const Point3d & p3 = stlgeometry->GetPoint(st[2]);

	      double cs = (Dist(p1,p2)+Dist(p2,p3)+Dist(p3,p1))/100.;

	      const Point3d & p = stlgeometry->GetPoint(st[nodeofseltrig-1]);
	      
	      glLineWidth (4);
	      glBegin (GL_LINES);
	      glVertex3f(p.X()+cs, p.Y()+cs, p.Z()+cs);
	      glVertex3f(p.X()-cs, p.Y()-cs, p.Z()-cs);
	      
	      glVertex3f(p.X()-cs, p.Y()+cs, p.Z()+cs);
	      glVertex3f(p.X()+cs, p.Y()-cs, p.Z()-cs);

	      glVertex3f(p.X()-cs, p.Y()+cs, p.Z()+cs);
	      glVertex3f(p.X()+cs, p.Y()-cs, p.Z()-cs);
	      
	      glVertex3f(p.X()+cs, p.Y()-cs, p.Z()+cs);
	      glVertex3f(p.X()-cs, p.Y()+cs, p.Z()-cs);
	      
	      glEnd ();	  
	      glLineWidth (1);
	    }
	  else if (stldoctor.selectmode == 1 || 
		   stldoctor.selectmode == 3 || 
		   stldoctor.selectmode == 4)
	    {
	      //multiedge
	      
	      const Array<twoint>& me = stlgeometry->SelectedMultiEdge();
	      if (stlgeometry->GetSelectTrig() > 0 && 
		  stlgeometry->GetSelectTrig() <= stlgeometry->GetNT() &&
		  me.Size())
		{

		  int en = stlgeometry->EdgeDataList().GetEdgeNum(me.Get(1).i1,me.Get(1).i2);
		  int status = stlgeometry->EdgeDataList().Get(en).GetStatus();
		  
		  switch (status)
		    {
		    case ED_CONFIRMED:
		      glMaterialfv (GL_FRONT_AND_BACK, 
				    GL_AMBIENT_AND_DIFFUSE, mat_collgreen);
		      break;
		    case ED_CANDIDATE:
		      glMaterialfv (GL_FRONT_AND_BACK, 
				    GL_AMBIENT_AND_DIFFUSE, mat_collbrown);
		      break;
		    case ED_EXCLUDED:
		      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_collred);
		      break;
		    }

		  glLineWidth (2);
		  glBegin (GL_LINES);
		  for (j = 1; j <= me.Size(); j++)
		    { 
		      Point3d p1 = stlgeometry->GetPoint(me.Get(j).i1);
		      Point3d p2 = stlgeometry->GetPoint(me.Get(j).i2);
		      
		      glVertex3f(p1.X(), p1.Y(), p1.Z());
		      glVertex3f(p2.X(), p2.Y(), p2.Z());
		    }
		  glEnd ();
		  glLineWidth (1);
		}
	    }
	}

      int showmarktrias = vispar.stlshowmarktrias || vispar.stlshowactivechart;
 
      if (stldoctor.showmarkedtrigs)
	{
	  //(*mycout) << "marked" << endl;
	  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE); //GL_LINE
	  glPolygonOffset (pgoff*1, pgoff*1);
	  glEnable (GL_POLYGON_OFFSET_FILL);
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colbluegreen);
	  glEnable (GL_NORMALIZE);

	  glBegin (GL_TRIANGLES);

	  for (j = 1; j <= stlgeometry -> GetNT(); j++)
	    {
	      if (stldoctor.showvicinity && !stlgeometry->Vicinity(j)) 
		{continue;}

	      if (!stlgeometry->IsMarkedTrig(j)) 
		{continue;}
	      
	      const STLTriangle& st = stlgeometry -> GetTriangle(j);

	      const Vec3d & n = stlgeometry->GetTriangle(j).Normal();
	      glNormal3f (n.X(), n.Y(), n.Z());
	      /*
	      const STLReadTriangle& tria = stlgeometry -> GetReadTriangle(j);
	      glNormal3f (tria.normal.X(),
			  tria.normal.Y(),
			  tria.normal.Z());
	      */
	      for (k = 0; k < 3; k++)
		{
		  const Point3d & p = stlgeometry->GetPoint(st[k]);
		  glVertex3f (p.X(), p.Y(), p.Z());
		}
	    }    
	  glEnd ();

	  //show OpenSegments on original geometry
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colviolet);
	  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	  glPolygonOffset (pgoff*1, 1);
      
	  glEnable (GL_NORMALIZE);
      
	  glBegin (GL_LINES);

	  if (stlgeometry->GetNMarkedSegs())
	    {
	      Point<3> p1,p2;	      
	      for (j = 1; j <= stlgeometry -> GetNMarkedSegs(); j++)
		{
		  stlgeometry->GetMarkedSeg(j,p1,p2);
		  glVertex3dv(&p1(0));
		  glVertex3dv(&p2(0));
		}
	    }
	  glEnd ();
	}


      if (stldoctor.showfaces)
	{
	  int facenumber = vispar.stlchartnumber + vispar.stlchartnumberoffset;

	  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	  glPolygonOffset (pgoff*3, 3);
	  glEnable (GL_POLYGON_OFFSET_FILL);
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_collgrey);
	  glEnable (GL_NORMALIZE);

	  glBegin (GL_TRIANGLES);

	  for (j = 1; j <= stlgeometry -> GetNT(); j++)
	    {
	      if (stldoctor.showvicinity && !stlgeometry->Vicinity(j)) 
		{continue;}

	      //(*mycout) << " facenum = " << stlgeometry->GetTriangle(j).GetFaceNum() << " ";
	      if (stlgeometry->GetTriangle(j).GetFaceNum() != facenumber) 
		{continue;}
	      
	      const STLTriangle& st = stlgeometry -> GetTriangle(j);

	      const Vec3d & n = stlgeometry->GetTriangle(j).Normal();
	      glNormal3f (n.X(), n.Y(), n.Z());
	      /*
	      const STLReadTriangle& tria = stlgeometry -> GetReadTriangle(j);
	      glNormal3f (tria.normal.X(),
			  tria.normal.Y(),
			  tria.normal.Z());
	      */
	      for (k = 0; k < 3; k++)
		{
		  Point3d p = stlgeometry->GetPoint(st[k]);
		  glVertex3f (p.X(), p.Y(), p.Z());
		}
	    }    
	  glEnd ();
	}

      if (showmarktrias && stlgeometry->AtlasMade())
	{
	  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	  glPolygonOffset (pgoff*3, 3);
	  glEnable (GL_POLYGON_OFFSET_FILL);

	  glBegin (GL_TRIANGLES);
	  
	  if (chartnumber >= 1 && chartnumber <= stlgeometry->GetNOCharts())
	    {
	      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colbrown);
	      const STLChart& chart = stlgeometry->GetChart(chartnumber);
	      for (j = 1; j <= chart.GetNChartT(); j++)
		{
		  /*
		  if (j == charttrignumber) 
		    {glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colred);}
		  else
		    {glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colbrown);}
		  */
		  const STLTriangle& st = stlgeometry -> GetTriangle(chart.GetChartTrig(j));

		  
		  const Vec3d & n = stlgeometry->GetTriangle(chart.GetChartTrig(j)).Normal();
		  glNormal3f (n.X(), n.Y(), n.Z());
		  /*
		  const STLReadTriangle& tria = stlgeometry -> GetReadTriangle(chart.GetChartTrig(j));
		  glNormal3f (tria.normal.X(),
			      tria.normal.Y(),
			      tria.normal.Z());
		  */
		  for (k = 0; k < 3; k++)
		    {
		      glVertex3f (stlgeometry->GetPoint(st[k])(0),
				  stlgeometry->GetPoint(st[k])(1),
				  stlgeometry->GetPoint(st[k])(2));
		    }
		}
	      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colgreen);
	      
	      for (j = 1; j <= chart.GetNOuterT(); j++)
		{
		  
		  const STLTriangle& st = stlgeometry -> GetTriangle(chart.GetOuterTrig(j));

		  const Vec3d & n = stlgeometry->GetTriangle(chart.GetOuterTrig(j)).Normal();
		  glNormal3f (n.X(), n.Y(), n.Z());


		  /*
		  const STLReadTriangle& tria = stlgeometry -> GetReadTriangle(chart.GetOuterTrig(j));
		  glNormal3f (tria.normal.X(),
			      tria.normal.Y(),
			      tria.normal.Z());
		  */
		  for (k = 0; k < 3; k++)
		    {
		      glVertex3f (stlgeometry->GetPoint(st[k])(0),
				  stlgeometry->GetPoint(st[k])(1),
				  stlgeometry->GetPoint(st[k])(2));
		    }
		}
	    }
	  glEnd ();
	}

      int showtrias = vispar.stlshowtrias;

      if (showtrias)
	{
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colgrey);
	  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	  glPolygonOffset (pgoff*2, 2);
	  glEnable (GL_POLYGON_OFFSET_FILL);
	  glEnable (GL_NORMALIZE);

	  glBegin (GL_TRIANGLES);
	  
	  for (j = 1; j <= stlgeometry -> GetNT(); j++)
	    {	  
	      if (stldoctor.showvicinity && !stlgeometry->Vicinity(j)) {continue;}

	      const STLTriangle& st = stlgeometry -> GetTriangle(j);

	      const Vec3d & n = stlgeometry->GetTriangle(j).Normal();
	      glNormal3f (n.X(), n.Y(), n.Z());
	      /*
	      const STLReadTriangle& tria = stlgeometry -> GetReadTriangle(j);
	      glNormal3f (tria.normal.X(),
			  tria.normal.Y(),
			  tria.normal.Z());
	      */	  
	      for (k = 0; k < 3; k++)
		{
		  glVertex3f (stlgeometry->GetPoint(st[k])(0),
			      stlgeometry->GetPoint(st[k])(1),
			      stlgeometry->GetPoint(st[k])(2));
		}
	    }    
	  glEnd ();
	} 

      int showedges = vispar.stlshowedges;
      
      if (showedges)
	{
	  glPolygonOffset (pgoff*1, 1);
	  glEnable (GL_POLYGON_OFFSET_FILL);
	  //glDisable (GL_POLYGON_OFFSET_FILL);      

	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colgreen);
	  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      
	  glEnable (GL_NORMALIZE);
      
	  glBegin (GL_LINES);

	  /*
	  if (stldoctor.useexternaledges)
	    {
	      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colorange);
	      for (j = 1; j <= stlgeometry -> NOExternalEdges(); j++)
		{
		  twoint v = stlgeometry->GetExternalEdge(j);
		  Point3d p1 = stlgeometry->GetPoint(v.i1);
		  Point3d p2 = stlgeometry->GetPoint(v.i2);
		  
		  Vec3d n1 = stlgeometry->GetNormal(v.i1);
		  Vec3d n2 = stlgeometry->GetNormal(v.i2);
		  
		  glNormal3f(n1.X(), n1.Y(), n1.Z());
		  glVertex3f(p1.X(), p1.Y(), p1.Z());
		  glNormal3f(n2.X(), n2.Y(), n2.Z());
		  glVertex3f(p2.X(), p2.Y(), p2.Z());
		}
	    }
	  */

	  
	  if (!stlgeometry->meshlines.Size() || !stldoctor.drawmeshededges)
	    {
	      /*
	      for (j = 1; j <= stlgeometry -> GetNE(); j++)
		{
		  STLEdge v = stlgeometry->GetEdge(j);
		  Point3d p1 = stlgeometry->GetPoint(v.pts[0]);
		  Point3d p2 = stlgeometry->GetPoint(v.pts[1]);
		  
		  Vec3d n1 = stlgeometry->GetNormal(v.pts[0]);
		  Vec3d n2 = stlgeometry->GetNormal(v.pts[1]);
		  
		  glNormal3f(n1.X(), n1.Y(), n1.Z());
		  glVertex3f(p1.X(), p1.Y(), p1.Z());
		  glNormal3f(n2.X(), n2.Y(), n2.Z());
		  glVertex3f(p2.X(), p2.Y(), p2.Z());
		}
	      */
	      const STLEdgeDataList& ed = stlgeometry->EdgeDataList();
	      for (i = 1; i <= ed.Size(); i++)
		{
		  if (ed.Get(i).GetStatus() != ED_UNDEFINED)
		    {
		      switch (ed.Get(i).GetStatus())
			{
			case ED_CONFIRMED:
			  glMaterialfv (GL_FRONT_AND_BACK, 
					GL_AMBIENT_AND_DIFFUSE, mat_colgreen);
			  break;
			case ED_CANDIDATE:
			  glMaterialfv (GL_FRONT_AND_BACK, 
					GL_AMBIENT_AND_DIFFUSE, mat_colbrown);
			  break;
			case ED_EXCLUDED:
			  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colred);
			  break;
			}

		      if (ed.Get(i).GetStatus() == ED_EXCLUDED && !stldoctor.showexcluded) continue;

		      Point3d p1 = stlgeometry->GetPoint(ed.Get(i).PNum(1));
		      Point3d p2 = stlgeometry->GetPoint(ed.Get(i).PNum(2));
		      glVertex3f(p1.X(), p1.Y(), p1.Z());
		      glVertex3f(p2.X(), p2.Y(), p2.Z());		   
		    }
		}
	    }

	  /*
	  else     
	  if (stlgeometry->meshlines.Size() == 0)
	    {
	      for (j = 1; j <= stlgeometry->GetNLines(); j++)
		{
		  STLLine* line = stlgeometry->GetLine(j);
		  int pn1, pn2;
		  for (int k = 1; k <= line->NP()-1; k++)
		    {
		      pn1 = line->PNum(k);
		      pn2 = line->PNum(k+1);

		      Point3d p1 = stlgeometry->GetPoint(pn1);
		      Point3d p2 = stlgeometry->GetPoint(pn2);
		  
		      Vec3d n1 = stlgeometry->GetNormal(pn1);
		      Vec3d n2 = stlgeometry->GetNormal(pn2);
		  
		      glNormal3f(n1.X(), n1.Y(), n1.Z());
		      glVertex3f(p1.X(), p1.Y(), p1.Z());
		      glNormal3f(n2.X(), n2.Y(), n2.Z());
		      glVertex3f(p2.X(), p2.Y(), p2.Z());
		    }
		}    
	    }
	  */
	    
	  else if (stlgeometry->meshlines.Size() != 0)
	    {
	      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colgreen);
	      for (j = 1; j <= stlgeometry->meshlines.Size(); j++)
		{
		  STLLine* line = stlgeometry->meshlines.Get(j);
		  int pn1, pn2;
		  for (int k = 1; k <= line->NP()-1; k++)
		    {
		      pn1 = line->PNum(k);
		      pn2 = line->PNum(k+1);

		      Point3d p1 = stlgeometry->meshpoints.Get(pn1);
		      Point3d p2 = stlgeometry->meshpoints.Get(pn2);
		  		  
		      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colgreen);
		      glVertex3f(p1.X(), p1.Y(), p1.Z());
		      glVertex3f(p2.X(), p2.Y(), p2.Z());

		      
		      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colred);
		      double cs = 0.02*Dist(p1,p2);
		      glVertex3f(p1.X()+cs, p1.Y()+cs, p1.Z()+cs);
		      glVertex3f(p1.X()-cs, p1.Y()-cs, p1.Z()-cs);
		      glVertex3f(p2.X()+cs, p2.Y()+cs, p2.Z()+cs);
		      glVertex3f(p2.X()-cs, p2.Y()-cs, p2.Z()-cs);

		      glVertex3f(p1.X()-cs, p1.Y()+cs, p1.Z()+cs);
		      glVertex3f(p1.X()+cs, p1.Y()-cs, p1.Z()-cs);
		      glVertex3f(p2.X()-cs, p2.Y()+cs, p2.Z()+cs);
		      glVertex3f(p2.X()+cs, p2.Y()-cs, p2.Z()-cs);
		      
		    }
		}
	    }
	    

	  glEnd ();
	}

      if (stldoctor.showedgecornerpoints && stlgeometry->LineEndPointsSet())
	{
	  glPointSize (5);
	  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_colred);
	  glBegin (GL_POINTS);
	  for (i = 1; i <= stlgeometry->GetNP(); i++)
	    {
	      if (stlgeometry->IsLineEndPoint(i))
		{
		  const Point3d p = stlgeometry->GetPoint(i);
		  glVertex3f (p.X(), p.Y(), p.Z());
		}
	    }
	  glEnd();
	  
	}


    }

 
  glPopMatrix();

  if (vispar.colormeshsize)
    DrawColorBar (hmin, hmax, 1);

  glFinish();  
}


void VisualSceneSTLMeshing :: BuildScene (int zoomall)
{
  if (selecttrig && zoomall == 2)
    center = stlgeometry -> GetPoint ( stlgeometry->GetTriangle(selecttrig).PNum(nodeofseltrig));
  else
    center = stlgeometry -> GetBoundingBox().Center();

  rad = stlgeometry -> GetBoundingBox().Diam() / 2;

  CalcTransformationMatrices();
}



void VisualSceneSTLMeshing :: MouseDblClick (int px, int py)
{
  //  (*mycout) << "dblclick: " << px << " - " << py << endl;
  

  int i, j, k, hits;

  // select surface triangle by mouse click

  GLuint selbuf[10000];
  glSelectBuffer (10000, selbuf);


  glRenderMode (GL_SELECT);

  GLint viewport[4];
  glGetIntegerv (GL_VIEWPORT, viewport);

  /*  
  (*mycout) << "viewport = " << viewport[0] << " " 
       << viewport[1] << " " << viewport[2] << " " << viewport[3] << endl;
  */

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


  glEnable (GL_POLYGON_OFFSET_FILL);
  for (j = 1; j <= stlgeometry -> GetNT(); j++)
    {
      if (stldoctor.showvicinity && !stlgeometry->Vicinity(j)) {continue;}

      const STLTriangle& st = stlgeometry -> GetTriangle(j);
			
      //const STLReadTriangle& tria = stlgeometry -> GetReadTriangle(j);
      //glNormal3f (tria.normal.X(), tria.normal.Y(), tria.normal.Z());
      
      if (stldoctor.selectmode == 0)
	{
	  glLoadName (j);
	  glBegin (GL_TRIANGLES);
	  for (k = 0; k < 3; k++)
	    {
	      Point3d p = stlgeometry->GetPoint(st[k]);
	      glVertex3f (p.X(), p.Y(), p.Z());
	    }
	  glEnd ();
	} 
      else if (stldoctor.selectmode == 1 || stldoctor.selectmode == 3
	        || stldoctor.selectmode == 4)
	{
	  Point3d pm = Center(stlgeometry->GetPoint(st[0]),
			      stlgeometry->GetPoint(st[1]),
			      stlgeometry->GetPoint(st[2]));

	  for (k = 0; k < 3; k++)
	    {
	      glLoadName (j*3+k-2);
	      glBegin (GL_TRIANGLES);

	      Point3d p1 = stlgeometry->GetPoint(st[k]);
	      Point3d p2 = stlgeometry->GetPoint(st[(k+1)%3]);
	      glVertex3f (p1.X(), p1.Y(), p1.Z());
	      glVertex3f (p2.X(), p2.Y(), p2.Z());
	      glVertex3f (pm.X(), pm.Y(), pm.Z());

	      glEnd ();
	    }
	}
      else
	{
	  Point3d pm1 = Center(stlgeometry->GetPoint(st[0]),
			       stlgeometry->GetPoint(st[1]));
	  Point3d pm2 = Center(stlgeometry->GetPoint(st[1]),
			       stlgeometry->GetPoint(st[2]));
	  Point3d pm3 = Center(stlgeometry->GetPoint(st[2]),
			       stlgeometry->GetPoint(st[0]));

	  Point3d p1 = stlgeometry->GetPoint(st[0]);
	  Point3d p2 = stlgeometry->GetPoint(st[1]);
	  Point3d p3 = stlgeometry->GetPoint(st[2]);

	  glLoadName (j*4-3);
	  glBegin (GL_TRIANGLES);
	  glVertex3f (p1.X(), p1.Y(), p1.Z());
	  glVertex3f (pm1.X(), pm1.Y(), pm1.Z());
	  glVertex3f (pm3.X(), pm3.Y(), pm3.Z());
	  glEnd ();

	  glLoadName (j*4-2);
	  glBegin (GL_TRIANGLES);
	  glVertex3f (p2.X(), p2.Y(), p2.Z());
	  glVertex3f (pm2.X(), pm2.Y(), pm2.Z());
	  glVertex3f (pm1.X(), pm1.Y(), pm1.Z());
	  glEnd ();

	  glLoadName (j*4-1);
	  glBegin (GL_TRIANGLES);
	  glVertex3f (p3.X(), p3.Y(), p3.Z());
	  glVertex3f (pm3.X(), pm3.Y(), pm3.Z());
	  glVertex3f (pm2.X(), pm2.Y(), pm2.Z());
	  glEnd ();

	  glLoadName (j*4);
	  glBegin (GL_TRIANGLES);
	  glVertex3f (pm1.X(), pm1.Y(), pm1.Z());
	  glVertex3f (pm2.X(), pm2.Y(), pm2.Z());
	  glVertex3f (pm3.X(), pm3.Y(), pm3.Z());
	  glEnd ();
	}
    }    

  glPopName();

  glMatrixMode (GL_PROJECTION); 
  glPopMatrix();

  glMatrixMode (GL_MODELVIEW); 
  glPopMatrix();

  glFlush();  

	
  hits = glRenderMode (GL_RENDER);

  //  (*mycout) << "hits = " << hits << endl;

  //int minrec = -1;
  int minname = 0;
  GLuint mindepth = 0;
  for (i = 0; i < hits; i++)
    {
      int curname = selbuf[4*i+3];
      GLuint curdepth = selbuf[4*i+1];

      /*      
      (*mycout) << selbuf[4*i] << " " << selbuf[4*i+1] << " " 
	   << selbuf[4*i+2] << " " << selbuf[4*i+3] << endl;
      */
      if (curname &&
	  (curdepth < mindepth || !minname))
	{
	  //minrec = i;
	  mindepth = curdepth;
	  minname = curname;
	}
    }

  if (!minname) {return;}
  
  if (stldoctor.selectmode == 0)
    {
      int oldtrig = selecttrig;
      selecttrig = minname;
      if (selecttrig == oldtrig)
	nodeofseltrig = (nodeofseltrig % 3) + 1;
      else
	nodeofseltrig = 1;

      stlgeometry->SetSelectTrig(selecttrig);
      stlgeometry->SetNodeOfSelTrig(nodeofseltrig);
      stlgeometry->PrintSelectInfo();
      
    }
  else if (stldoctor.selectmode == 1 || stldoctor.selectmode == 3 || stldoctor.selectmode == 4)
    {
      selecttrig = (minname-1) / 3 + 1;
      nodeofseltrig = minname-selecttrig*3+3;

      stlgeometry->SetSelectTrig(selecttrig);
      stlgeometry->SetNodeOfSelTrig(nodeofseltrig);
      stlgeometry->PrintSelectInfo();

      if (stldoctor.selectmode == 1)
	{
	  stlgeometry->BuildSelectedEdge(twoint(stlgeometry->GetTriangle(selecttrig).PNumMod(nodeofseltrig),
						stlgeometry->GetTriangle(selecttrig).PNumMod(nodeofseltrig+1)));
	}
      if (stldoctor.selectmode == 3)
	{
	  stlgeometry->BuildSelectedMultiEdge(twoint(stlgeometry->GetTriangle(selecttrig).PNumMod(nodeofseltrig),
						     stlgeometry->GetTriangle(selecttrig).PNumMod(nodeofseltrig+1)));
	}
      else if (stldoctor.selectmode == 4)
	{
	  stlgeometry->BuildSelectedCluster(twoint(stlgeometry->GetTriangle(selecttrig).PNumMod(nodeofseltrig),
						   stlgeometry->GetTriangle(selecttrig).PNumMod(nodeofseltrig+1)));
	}
 
      switch (stldoctor.edgeselectmode)
	{
	case 1: stlgeometry->STLDoctorUndefinedEdge(); break;
	case 2: stlgeometry->STLDoctorConfirmEdge(); break;
	case 3: stlgeometry->STLDoctorCandidateEdge(); break;
	case 4: stlgeometry->STLDoctorExcludeEdge(); break;
	default: break;
	}
    }
  else if (stldoctor.selectmode == 2)
    {
      selecttrig = (minname-1) / 4 + 1;
      nodeofseltrig = minname-selecttrig*4+4;
      if (nodeofseltrig == 4) {nodeofseltrig = 1;}

      stlgeometry->SetSelectTrig(selecttrig);
      stlgeometry->SetNodeOfSelTrig(nodeofseltrig);
      stlgeometry->PrintSelectInfo();

    }

  if (stldoctor.showtouchedtrigchart && stlgeometry->AtlasMade() && stlgeometry->GetSelectTrig())
    {
      vispar.stlchartnumber =  stlgeometry->GetChartNr(stlgeometry->GetSelectTrig());
      vispar.stlchartnumberoffset = 0;
    }
  
}




VisualSceneSTLMeshing vsstlmeshing;

#endif



 
}
