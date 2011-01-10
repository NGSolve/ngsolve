#ifndef NOTCL

#ifdef OCCGEOMETRY

#include <mystdlib.h>
#include <myadt.hpp>
#include <meshing.hpp>

#include <occgeom.hpp>

#include "TopoDS_Shape.hxx"
#include "TopoDS_Vertex.hxx"
#include "TopExp_Explorer.hxx"
#include "BRep_Tool.hxx"
#include "TopoDS.hxx"
#include "gp_Pnt.hxx"
#include "Geom_Curve.hxx"
#include "Poly_Triangulation.hxx"
#include "Poly_Array1OfTriangle.hxx"
#include "TColgp_Array1OfPnt2d.hxx"
#include "Poly_Triangle.hxx"
#include "Poly_Polygon3D.hxx"
#include "Poly_PolygonOnTriangulation.hxx"

#include <visual.hpp>

#include "vsocc.hpp"

namespace netgen
{
  // extern OCCGeometry * occgeometry;

   /* *********************** Draw OCC Geometry **************** */

   VisualSceneOCCGeometry :: VisualSceneOCCGeometry ()
   : VisualScene()
   {
      trilists.SetSize(0);
      linelists.SetSize(1);

   }

   VisualSceneOCCGeometry :: ~VisualSceneOCCGeometry ()
   {
      ;
   }

   void VisualSceneOCCGeometry :: DrawScene ()
   {
      if ( occgeometry->changed )
      {
         BuildScene();
         occgeometry -> changed = 0;
      }

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
      
      //  glEnable (GL_LIGHTING);

      double shine = vispar.shininess;
      // double transp = vispar.transp;

      glMaterialf (GL_FRONT_AND_BACK, GL_SHININESS, shine);
      glLogicOp (GL_COPY);

      glEnable (GL_NORMALIZE);

      float mat_col[] = {  0.2f, 0.2f, 0.8f, 1.0f};
      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);

      glPolygonOffset (1, 1);
      glEnable (GL_POLYGON_OFFSET_FILL);

      // Philippose - 30/01/2009
      // Added clipping planes to Geometry view
      SetClippingPlane();

      GLfloat matcoledge[] = {  0, 0, 1, 1};
      GLfloat matcolhiedge[] = {  1, 0, 0, 1};

      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matcoledge);
      glLineWidth (1.0f);

      if (vispar.occshowedges) glCallList (linelists.Get(1));
      if (vispar.occshowsurfaces) glCallList (trilists.Get(1));

      glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matcolhiedge);
      glLineWidth (5.0f);

      if (vispar.occshowedges) glCallList (linelists.Get(2));

      for (int i = 1; i <= occgeometry->vmap.Extent(); i++)
      if (occgeometry->vvispar[i-1].IsHighlighted())
      {
         glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matcolhiedge);
         glLineWidth (5.0f);

         glBegin (GL_LINES);

         gp_Pnt p = BRep_Tool::Pnt(TopoDS::Vertex(occgeometry->vmap(i)));
         double d = rad/100;
         glVertex3f (p.X()-d, p.Y(), p.Z());
         glVertex3f (p.X()+d, p.Y(), p.Z());
         glVertex3f (p.X(), p.Y()-d, p.Z());
         glVertex3f (p.X(), p.Y()+d, p.Z());
         glVertex3f (p.X(), p.Y(), p.Z()-d);
         glVertex3f (p.X(), p.Y(), p.Z()+d);
         glEnd();
      }

      glDisable (GL_POLYGON_OFFSET_FILL);

      glPopMatrix();
      //  DrawCoordinateCross ();
      //  DrawNetgenLogo ();
      glFinish();

      glDisable (GL_POLYGON_OFFSET_FILL);
   }

   /*
    void VisualSceneOCCGeometry :: BuildScene (int zoomall)
    {
    int i = 0, j, k;

    TopExp_Explorer ex, ex_edge;

    if (vispar.occvisproblemfaces || (occgeometry -> changed != 2))
    {
    Box<3> bb = occgeometry -> GetBoundingBox();

    center = bb.Center();
    rad = bb.Diam() / 2;



    if (vispar.occvisproblemfaces)
    {
    for (i = 1; i <= occgeometry->fmap.Extent(); i++)
    if (occgeometry->facemeshstatus[i-1] == -1)
    {
    GProp_GProps system;
    BRepGProp::LinearProperties(occgeometry->fmap(i), system);
    gp_Pnt pnt = system.CentreOfMass();
    center = Point<3> (pnt.X(), pnt.Y(), pnt.Z());
    cout << "Setting center to mid of face " << i << " = " << center << endl;
    }
    }


    CalcTransformationMatrices();
    }


    for (i = 1; i <= linelists.Size(); i++)
    glDeleteLists (linelists.Elem(i), 1);
    linelists.SetSize(0);

    linelists.Append (glGenLists (1));
    glNewList (linelists.Last(), GL_COMPILE);

    i = 0;
    for (ex_edge.Init(occgeometry -> shape, TopAbs_EDGE);
    ex_edge.More(); ex_edge.Next())
    {
    if (BRep_Tool::Degenerated(TopoDS::Edge(ex_edge.Current()))) continue;
    i++;


    TopoDS_Edge edge = TopoDS::Edge(ex_edge.Current());

    Handle(Poly_PolygonOnTriangulation) aEdgePoly;
    Handle(Poly_Triangulation) T;
    TopLoc_Location aEdgeLoc;
    BRep_Tool::PolygonOnTriangulation(edge, aEdgePoly, T, aEdgeLoc);

    if(aEdgePoly.IsNull())
    {
    cout << "cannot visualize edge " << i << endl;
    continue;
    }

    glBegin (GL_LINE_STRIP);

    int nbnodes = aEdgePoly -> NbNodes();
    for (j = 1; j <= nbnodes; j++)
    {
    gp_Pnt p = (T -> Nodes())(aEdgePoly->Nodes()(j)).Transformed(aEdgeLoc);
    glVertex3f (p.X(), p.Y(), p.Z());
    }

    glEnd ();


    }

    glEndList ();

    for (i = 1; i <= trilists.Size(); i++)
    glDeleteLists (trilists.Elem(i), 1);
    trilists.SetSize(0);


    trilists.Append (glGenLists (1));
    glNewList (trilists.Last(), GL_COMPILE);

    i = 0;

    TopExp_Explorer exp0, exp1, exp2, exp3;
    int shapenr = 0;
    for (exp0.Init(occgeometry -> shape, TopAbs_SOLID); exp0.More(); exp0.Next())
    {
    shapenr++;

    if (vispar.occshowvolumenr != 0 &&
    vispar.occshowvolumenr != shapenr) continue;

    float mat_col[4];
    mat_col[3] = 1;
    switch (shapenr)
    {
    case 1:
    mat_col[0] = 0.2;
    mat_col[1] = 0.2;
    mat_col[2] = 0.8;
    break;
    case 2:
    mat_col[0] = 0.8;
    mat_col[1] = 0.2;
    mat_col[2] = 0.8;
    break;
    case 3:
    mat_col[0] = 0.2;
    mat_col[1] = 0.8;
    mat_col[2] = 0.8;
    break;
    case 4:
    mat_col[0] = 0.8;
    mat_col[1] = 0.2;
    mat_col[2] = 0.2;
    break;
    case 5:
    mat_col[0] = 0.8;
    mat_col[1] = 0.8;
    mat_col[2] = 0.8;
    break;
    case 6:
    mat_col[0] = 0.6;
    mat_col[1] = 0.6;
    mat_col[2] = 0.6;
    break;
    case 7:
    mat_col[0] = 0.2;
    mat_col[1] = 0.8;
    mat_col[2] = 0.2;
    break;
    case 8:
    mat_col[0] = 0.8;
    mat_col[1] = 0.8;
    mat_col[2] = 0.2;
    break;
    default:
    //	  mat_col[0] = 1-(1.0/double(shapenr));
    //	  mat_col[1] = 0.5;
    mat_col[0] = 0.5+double((shapenr*shapenr*shapenr*shapenr) % 10)/20.0;
    mat_col[1] = 0.5+double(int(shapenr*shapenr*shapenr*shapenr*sin(double(shapenr))) % 10)/20.0;
    mat_col[2] = 0.5+double((shapenr*shapenr*shapenr) % 10)/20.0;
    }

    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);

    for (exp1.Init(exp0.Current(), TopAbs_SHELL); exp1.More(); exp1.Next())
    for (exp2.Init(exp1.Current().Composed(exp0.Current().Orientation()), TopAbs_FACE); exp2.More(); exp2.Next())
    {
    TopoDS_Face face = TopoDS::Face (exp2.Current().Composed(exp1.Current().Orientation()));

    i = occgeometry->fmap.FindIndex(face);

    TopLoc_Location loc;
    Handle(Geom_Surface) surf = BRep_Tool::Surface (face);
    BRepAdaptor_Surface sf(face, Standard_False);
    BRepLProp_SLProps prop(sf, 1, 1e-5);
    Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation (face, loc);

    if (triangulation.IsNull())
    {
    cout << "cannot visualize face " << i << endl;
    continue;
    }

    if (vispar.occvisproblemfaces)
    {
    switch (occgeometry->facemeshstatus[i-1])
    {
    case 0:
    mat_col[0] = 0.2;
    mat_col[1] = 0.2;
    mat_col[2] = 0.8;
    break;
    case 1:
    mat_col[0] = 0.2;
    mat_col[1] = 0.8;
    mat_col[2] = 0.2;
    break;
    case -1:
    mat_col[0] = 0.8;
    mat_col[1] = 0.2;
    mat_col[2] = 0.2;
    break;
    }
    glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);

    }
    glBegin (GL_TRIANGLES);

    int ntriangles = triangulation -> NbTriangles();
    for (j = 1; j <= ntriangles; j++)
    {
    Poly_Triangle triangle = (triangulation -> Triangles())(j);
    for (k = 1; k <= 3; k++)
    {
    gp_Pnt2d uv = (triangulation -> UVNodes())(triangle(k));
    gp_Pnt pnt;
    gp_Vec du, dv;
    prop.SetParameters (uv.X(), uv.Y());
    surf->D0 (uv.X(), uv.Y(), pnt);
    gp_Vec n;

    if (prop.IsNormalDefined())
    n = prop.Normal();
    else
    n = gp_Vec (0,0,0);

    if (face.Orientation() == TopAbs_REVERSED) n *= -1;
    glNormal3f (n.X(), n.Y(), n.Z());
    glVertex3f (pnt.X(), pnt.Y(), pnt.Z());
    }
    }
    glEnd ();

    }
    }


    glEndList ();

    }
    */

   void VisualSceneOCCGeometry :: BuildScene (int zoomall)
   {
     if (occgeometry -> changed == OCCGEOMETRYVISUALIZATIONFULLCHANGE)
       {
         occgeometry -> BuildVisualizationMesh (vispar.occdeflection);

         center = occgeometry -> Center();
         rad = occgeometry -> GetBoundingBox().Diam() / 2;

         if (vispar.occzoomtohighlightedentity)
         {
            bool hilite = false;
            bool hiliteonepoint = false;
            Bnd_Box bb;

            for (int i = 1; i <= occgeometry->fmap.Extent(); i++)
            if (occgeometry->fvispar[i-1].IsHighlighted())
            {
               hilite = true;
               BRepBndLib::Add (occgeometry->fmap(i), bb);
            }

            for (int i = 1; i <= occgeometry->emap.Extent(); i++)
            if (occgeometry->evispar[i-1].IsHighlighted())
            {
               hilite = true;
               BRepBndLib::Add (occgeometry->emap(i), bb);
            }

            for (int i = 1; i <= occgeometry->vmap.Extent(); i++)
            if (occgeometry->vvispar[i-1].IsHighlighted())
            {
               hiliteonepoint = true;
               BRepBndLib::Add (occgeometry->vmap(i), bb);
            }

            if (hilite || hiliteonepoint)
            {
               double x1,y1,z1,x2,y2,z2;
               bb.Get (x1,y1,z1,x2,y2,z2);
               Point<3> p1 = Point<3> (x1,y1,z1);
               Point<3> p2 = Point<3> (x2,y2,z2);
               Box<3> boundingbox(p1,p2);

               center = boundingbox.Center();
               if (hiliteonepoint)
               rad = occgeometry -> GetBoundingBox().Diam() / 100;
               else
               rad = boundingbox.Diam() / 2;
            }
         }

         CalcTransformationMatrices();
      }

      // Clear lists

      for (int i = 1; i <= linelists.Size(); i++)
      glDeleteLists (linelists.Elem(i), 1);
      linelists.SetSize(0);

      for (int i = 1; i <= trilists.Size(); i++)
      glDeleteLists (trilists.Elem(i), 1);
      trilists.SetSize(0);

      // Total wireframe

      linelists.Append (glGenLists (1));
      glNewList (linelists.Last(), GL_COMPILE);

      for (int i = 1; i <= occgeometry->emap.Extent(); i++)
      {
         TopoDS_Edge edge = TopoDS::Edge(occgeometry->emap(i));
         if (BRep_Tool::Degenerated(edge)) continue;
         if (occgeometry->evispar[i-1].IsHighlighted()) continue;

         Handle(Poly_PolygonOnTriangulation) aEdgePoly;
         Handle(Poly_Triangulation) T;
         TopLoc_Location aEdgeLoc;
         BRep_Tool::PolygonOnTriangulation(edge, aEdgePoly, T, aEdgeLoc);

         if(aEdgePoly.IsNull())
         {
            (*testout) << "visualizing edge " << occgeometry->emap.FindIndex (edge)
            << " without using the occ visualization triangulation" << endl;

            double s0, s1;
            Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);

            glBegin (GL_LINE_STRIP);
            for (int i = 0; i<=50; i++)
            {
               gp_Pnt p = c->Value (s0 + i*(s1-s0)/50.0);
               glVertex3f (p.X(),p.Y(),p.Z());
            }
            glEnd ();

            continue;
         }

         int nbnodes = aEdgePoly -> NbNodes();
         glBegin (GL_LINE_STRIP);
         for (int j = 1; j <= nbnodes; j++)
         {
            gp_Pnt p = (T -> Nodes())(aEdgePoly->Nodes()(j)).Transformed(aEdgeLoc);
            glVertex3f (p.X(), p.Y(), p.Z());
         }
         glEnd ();
      }

      glEndList ();

      // Highlighted edge list

      linelists.Append (glGenLists (1));
      glNewList (linelists.Last(), GL_COMPILE);

      for (int i = 1; i <= occgeometry->emap.Extent(); i++)
      if (occgeometry->evispar[i-1].IsHighlighted())
      {
         TopoDS_Edge edge = TopoDS::Edge(occgeometry->emap(i));
         if (BRep_Tool::Degenerated(edge)) continue;

         Handle(Poly_PolygonOnTriangulation) aEdgePoly;
         Handle(Poly_Triangulation) T;
         TopLoc_Location aEdgeLoc;
         BRep_Tool::PolygonOnTriangulation(edge, aEdgePoly, T, aEdgeLoc);

         if(aEdgePoly.IsNull())
         {
            (*testout) << "visualizing edge " << occgeometry->emap.FindIndex (edge)
            << " without using the occ visualization triangulation" << endl;

            double s0, s1;
            Handle(Geom_Curve) c = BRep_Tool::Curve(edge, s0, s1);

            glBegin (GL_LINE_STRIP);
            for (int i = 0; i<=50; i++)
            {
               gp_Pnt p = c->Value (s0 + i*(s1-s0)/50.0);
               glVertex3f (p.X(),p.Y(),p.Z());
            }
            glEnd ();

            continue;
         }

         int nbnodes = aEdgePoly -> NbNodes();
         glBegin (GL_LINE_STRIP);
         for (int j = 1; j <= nbnodes; j++)
         {
            gp_Pnt p = (T -> Nodes())(aEdgePoly->Nodes()(j)).Transformed(aEdgeLoc);
            glVertex3f (p.X(), p.Y(), p.Z());
         }
         glEnd ();
      }

      glEndList ();

      // display faces

      trilists.Append (glGenLists (1));
      glNewList (trilists.Last(), GL_COMPILE);

      for (int i = 1; i <= occgeometry->fmap.Extent(); i++)
      {
         if (!occgeometry->fvispar[i-1].IsVisible()) continue;

         glLoadName (i);
         float mat_col[4];
         mat_col[3] = 1;

         TopoDS_Face face = TopoDS::Face(occgeometry->fmap(i));

         if (!occgeometry->fvispar[i-1].IsHighlighted())
         {
            // Philippose - 30/01/2009
            // OpenCascade XDE Support
            Quantity_Color face_colour;
            // Philippose - 23/02/2009
            // Check to see if colours have been extracted first!!
            // Forum bug-fox (Jean-Yves - 23/02/2009)
            if(!(occgeometry->face_colours.IsNull())
               && (occgeometry->face_colours->GetColor(face,XCAFDoc_ColorSurf,face_colour)))
            {
               mat_col[0] = face_colour.Red();
               mat_col[1] = face_colour.Green();
               mat_col[2] = face_colour.Blue();
            }
            else
            {
               mat_col[0] = 0.0;
               mat_col[1] = 1.0;
               mat_col[2] = 0.0;
            }
         }
         else
         {
            mat_col[0] = 0.8;
            mat_col[1] = 0.2;
            mat_col[2] = 0.2;
         }

         glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_col);

         TopLoc_Location loc;
         Handle(Geom_Surface) surf = BRep_Tool::Surface (face);
         BRepAdaptor_Surface sf(face, Standard_False);
         BRepLProp_SLProps prop(sf, 1, 1e-5);
         Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation (face, loc);

         if (triangulation.IsNull())
         {
            cout << "cannot visualize face " << i << endl;
            occgeometry->fvispar[i-1].SetNotDrawable();
            continue;
         }

         gp_Pnt2d uv;
         gp_Pnt pnt;
         gp_Vec n;

         glBegin (GL_TRIANGLES);

         int ntriangles = triangulation -> NbTriangles();
         for (int j = 1; j <= ntriangles; j++)
         {
            Poly_Triangle triangle = (triangulation -> Triangles())(j);
            gp_Pnt p[3];
            for (int k = 1; k <= 3; k++)
            p[k-1] = (triangulation -> Nodes())(triangle(k)).Transformed(loc);

            for (int k = 1; k <= 3; k++)
            {
               uv = (triangulation -> UVNodes())(triangle(k));
               prop.SetParameters (uv.X(), uv.Y());

               //	      surf->D0 (uv.X(), uv.Y(), pnt);

               if (prop.IsNormalDefined())
               n = prop.Normal();
               else
               {
                  (*testout) << "Visualization of face " << i
                  << ": Normal vector not defined" << endl;
                  //		  n = gp_Vec (0,0,0);
                  gp_Vec a(p[0],p[1]);
                  gp_Vec b(p[0],p[2]);
                  n = b^a;
               }

               if (face.Orientation() == TopAbs_REVERSED) n *= -1;
               glNormal3f (n.X(), n.Y(), n.Z());
               glVertex3f (p[k-1].X(), p[k-1].Y(), p[k-1].Z());
            }
         }
         glEnd ();

      }
      glEndList ();

   }

   void SelectFaceInOCCDialogTree (int facenr);

   void VisualSceneOCCGeometry :: MouseDblClick (int px, int py)
   {
      int hits;

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

      glDisable(GL_CLIP_PLANE0);

      // Philippose - 30/01/2009
      // Enable clipping planes for Selection mode in OCC Geometry
      if (vispar.clipenable)
      {
         Vec<3> n(clipplane[0], clipplane[1], clipplane[2]);
         double len = Abs(n);
         double mu = -clipplane[3] / (len*len);
         Point<3> p (mu * n);
         n /= len;
         Vec<3> t1 = n.GetNormal ();
         Vec<3> t2 = Cross (n, t1);

         double xi1mid = (center - p) * t1;
         double xi2mid = (center - p) * t2;

         glLoadName (0);
         glBegin (GL_QUADS);
         glVertex3dv (p + (xi1mid-rad) * t1 + (xi2mid-rad) * t2);
         glVertex3dv (p + (xi1mid+rad) * t1 + (xi2mid-rad) * t2);
         glVertex3dv (p + (xi1mid+rad) * t1 + (xi2mid+rad) * t2);
         glVertex3dv (p + (xi1mid-rad) * t1 + (xi2mid+rad) * t2);
         glEnd ();
      }

      glCallList (trilists.Get(1));

      glDisable (GL_POLYGON_OFFSET_FILL);

      glPopName();

      glMatrixMode (GL_PROJECTION);
      glPopMatrix();

      glMatrixMode (GL_MODELVIEW);
      glPopMatrix();

      glFlush();

      hits = glRenderMode (GL_RENDER);

      int minname = 0;
      GLuint mindepth = 0;

      // find clippingplane
      GLuint clipdepth = 0; // GLuint(-1);

      for (int i = 0; i < hits; i++)
      {
         int curname = selbuf[4*i+3];
         if (!curname) clipdepth = selbuf[4*i+1];
      }

      for (int i = 0; i < hits; i++)
      {
         int curname = selbuf[4*i+3];
         GLuint curdepth = selbuf[4*i+1];
         if (curname && (curdepth> clipdepth) &&
               (curdepth < mindepth || !minname))
         {
            mindepth = curdepth;
            minname = curname;
         }
      }

      occgeometry->LowLightAll();

      if (minname)
      {
         occgeometry->fvispar[minname-1].Highlight();

         if (vispar.occzoomtohighlightedentity)
         occgeometry->changed = OCCGEOMETRYVISUALIZATIONFULLCHANGE;
         else
         occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
         cout << "Selected face: " << minname << endl;
      }
      else
      {
         occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
      }

      glDisable(GL_CLIP_PLANE0);

      SelectFaceInOCCDialogTree (minname);

      // Philippose - 30/01/2009
      // Set the currently selected face in the array
      // for local face mesh size definition
      occgeometry->SetSelectedFace(minname);

      //  selecttimestamp = NextTimeStamp();
   }

}

#endif

#endif // NOTCL
