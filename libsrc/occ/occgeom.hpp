#ifndef FILE_OCCGEOM
#define FILE_OCCGEOM

/* *************************************************************************/
/* File:   occgeom.hpp                                                     */
/* Author: Robert Gaisbauer                                                */
/* Date:   26. May  03                                                     */
/* *************************************************************************/

#ifdef OCCGEOMETRY

#include <meshing.hpp>

#include "BRep_Tool.hxx"
#include "Geom_Curve.hxx"
#include "Geom2d_Curve.hxx"
#include "Geom_Surface.hxx"
#include "GeomAPI_ProjectPointOnSurf.hxx"
#include "GeomAPI_ProjectPointOnCurve.hxx"
#include "BRepTools.hxx"
#include "TopExp.hxx"
#include "BRepBuilderAPI_MakeVertex.hxx"
#include "BRepBuilderAPI_MakeShell.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "BRepOffsetAPI_Sewing.hxx"
#include "BRepLProp_SLProps.hxx"
#include "BRepAdaptor_Surface.hxx"
#include "Poly_Triangulation.hxx"
#include "Poly_Array1OfTriangle.hxx"
#include "TColgp_Array1OfPnt2d.hxx"
#include "Poly_Triangle.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "Geom_Surface.hxx"
#include "TopExp.hxx"
#include "gp_Pnt.hxx"
#include "TopoDS.hxx"
#include "TopoDS_Solid.hxx"
#include "TopExp_Explorer.hxx"
#include "BRep_Tool.hxx"
#include "Geom_Curve.hxx"
#include "Geom2d_Curve.hxx"
#include "Geom_Surface.hxx"
#include "GeomAPI_ProjectPointOnSurf.hxx"
#include "GeomAPI_ProjectPointOnCurve.hxx"
#include "TopoDS_Wire.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "BRepTools.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp.hxx"
#include "BRepBuilderAPI_MakeVertex.hxx"
#include "BRepBuilderAPI_MakeShell.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "BRepOffsetAPI_Sewing.hxx"
#include "BRepLProp_CLProps.hxx"
#include "BRepLProp_SLProps.hxx"
#include "BRepAdaptor_Surface.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "Poly_Triangulation.hxx"
#include "Poly_Array1OfTriangle.hxx"
#include "TColgp_Array1OfPnt2d.hxx"
#include "Poly_Triangle.hxx"
#include "GProp_GProps.hxx"
#include "BRepGProp.hxx"
#include "IGESControl_Reader.hxx"
#include "STEPControl_Reader.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Face.hxx"
#include "IGESToBRep_Reader.hxx"
#include "Interface_Static.hxx"
#include "GeomAPI_ExtremaCurveCurve.hxx"
#include "Standard_ErrorHandler.hxx"
#include "Standard_Failure.hxx"
#include "ShapeUpgrade_ShellSewing.hxx"
#include "ShapeFix_Shape.hxx"
#include "ShapeFix_Wireframe.hxx"
#include "BRepMesh.hxx"
#include "BRepMesh_IncrementalMesh.hxx"
#include "BRepBndLib.hxx"
#include "Bnd_Box.hxx"
#include "ShapeAnalysis.hxx"
#include "ShapeBuild_ReShape.hxx"
#include "IGESControl_Writer.hxx"
#include "STEPControl_Writer.hxx"
#include "StlAPI_Writer.hxx"
#include "STEPControl_StepModelType.hxx"

namespace netgen
{

#include "../visualization/vispar.hpp"
  //  class VisualizationParameters;
  //  extern VisualizationParameters vispar;


#include "occmeshsurf.hpp"

#define PROJECTION_TOLERANCE 1e-10


#define ENTITYISVISIBLE 1
#define ENTITYISHIGHLIGHTED 2
#define ENTITYISDRAWABLE 4

class EntityVisualizationCode
{
  int code;

public:

  EntityVisualizationCode()
  { code = ENTITYISVISIBLE + !ENTITYISHIGHLIGHTED + ENTITYISDRAWABLE; }

  int IsVisible ()
  { return code & ENTITYISVISIBLE; }

  int IsHighlighted ()
  { return code & ENTITYISHIGHLIGHTED; }

  int IsDrawable ()
  { return code & ENTITYISDRAWABLE; }

  void Show ()
  { code |= ENTITYISVISIBLE; }

  void Hide ()
  { code &= ~ENTITYISVISIBLE; }

  void Highlight ()
  { code |= ENTITYISHIGHLIGHTED; }

  void Lowlight ()
  { code &= ~ENTITYISHIGHLIGHTED; }

  void SetDrawable ()
  { code |= ENTITYISDRAWABLE; }

  void SetNotDrawable ()
  { code &= ~ENTITYISDRAWABLE; }
};



inline double Det3 (double a00, double a01, double a02,
		    double a10, double a11, double a12,
		    double a20, double a21, double a22)
{
  return a00*a11*a22 + a01*a12*a20 + a10*a21*a02 - a20*a11*a02 - a10*a01*a22 - a21*a12*a00;
}



#define OCCGEOMETRYVISUALIZATIONNOCHANGE   0
#define OCCGEOMETRYVISUALIZATIONFULLCHANGE 1
  // == compute transformation matrices and redraw
#define OCCGEOMETRYVISUALIZATIONHALFCHANGE 2
  // == redraw

class OCCGeometry
{
  Point<3> center;

public:
  TopoDS_Shape shape;
  TopTools_IndexedMapOfShape fmap, emap, vmap, somap, shmap, wmap;
  ARRAY<bool> fsingular, esingular, vsingular;
  Box<3> boundingbox;

  int changed; 
  ARRAY<int> facemeshstatus;

  ARRAY<EntityVisualizationCode> fvispar, evispar, vvispar;

  double tolerance;
  bool fixsmalledges;
  bool fixspotstripfaces;
  bool sewfaces;
  bool makesolids;
  bool splitpartitions;


  OCCGeometry()
  {
    somap.Clear();
    shmap.Clear();
    fmap.Clear();
    wmap.Clear();
    emap.Clear();
    vmap.Clear();
  }


  void BuildFMap();

  Box<3> GetBoundingBox()
  { return boundingbox; }

  int NrSolids()
  { return somap.Extent(); }

  void SetCenter()
  { center = boundingbox.Center(); }

  Point<3> Center()
  { return center; }

  void Project (int surfi, Point<3> & p) const;
  bool FastProject (int surfi, Point<3> & ap, double& u, double& v) const;

 
  OCCSurface GetSurface (int surfi)
  {
    cout << "OCCGeometry::GetSurface using PLANESPACE" << endl;
    return OCCSurface (TopoDS::Face(fmap(surfi)), PLANESPACE);
  }
  

  void BuildVisualizationMesh ();

  void RecursiveTopologyTree (const TopoDS_Shape & sh,
			      stringstream & str,
			      TopAbs_ShapeEnum l,
			      bool free,
			      const char * lname);

  void GetTopologyTree (stringstream & str);

  void PrintNrShapes ();

  void CheckIrregularEntities (stringstream & str);

  void SewFaces();

  void MakeSolid();

  void HealGeometry();

  void LowLightAll()
  {
    for (int i = 1; i <= fmap.Extent(); i++)
      fvispar[i-1].Lowlight();
    for (int i = 1; i <= emap.Extent(); i++)
      evispar[i-1].Lowlight();
    for (int i = 1; i <= vmap.Extent(); i++)
      vvispar[i-1].Lowlight();
  }

  void GetUnmeshedFaceInfo (stringstream & str);
  void GetNotDrawableFaces (stringstream & str);
  bool ErrorInSurfaceMeshing ();

  void WriteOCC_STL(char * filename);
};


void PrintContents (OCCGeometry * geom);

OCCGeometry * LoadOCC_IGES (const char * filename);
OCCGeometry * LoadOCC_STEP (const char * filename);
OCCGeometry * LoadOCC_BREP (const char * filename);

}

#endif

#endif
