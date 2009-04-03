#include <mystdlib.h>


#include <meshing.hpp>
#include <csg.hpp>
#include <geometry2d.hpp>
#include <stlgeom.hpp>


#ifdef OCCGEOMETRY
#include <occgeom.hpp>
#endif

#ifdef ACIS
#include <acisgeom.hpp>
#endif

#ifdef SOCKETS
#include "../sockets/sockets.hpp"
#endif

#ifndef NOTCL
#include <visual.hpp>
#endif


#include "nginterface.h"
#include "nginterface_v2.hpp"

// #include <FlexLexer.h>


// #include <mystdlib.h>


namespace netgen
{
#include "writeuser.hpp"

	extern AutoPtr<Mesh> mesh;
#ifndef NOTCL
  extern VisualSceneMesh vsmesh;
  extern Tcl_Interp * tcl_interp;
#endif

  extern AutoPtr<SplineGeometry2d> geometry2d;
  extern AutoPtr<CSGeometry> geometry;
  extern STLGeometry * stlgeometry;

#ifdef OCCGEOMETRY
  extern OCCGeometry * occgeometry;
#endif
#ifdef ACIS
  extern ACISGeometry * acisgeometry;
#endif

#ifdef OPENGL
  extern VisualSceneSolution vssolution;
#endif
  extern CSGeometry * ParseCSG (istream & istr);

#ifdef SOCKETS
  extern AutoPtr<ClientSocket> clientsocket;
  //extern Array< AutoPtr < ServerInfo > > servers;
  extern Array< ServerInfo* > servers;
#endif
}


using namespace netgen;



template <> int Ng_GetNElements<1> ()
{
  return mesh->GetNSeg();
}

template <> int Ng_GetNElements<2> ()
{
  return mesh->GetNSE();
}

template <> int Ng_GetNElements<3> ()
{
  return mesh->GetNE();
}




template <> Ng_Element Ng_GetElement<1> (int nr)
{
  const Segment & el = mesh->LineSegment (SegmentIndex(nr));

  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  ret.npoints = el.GetNP();
  ret.points = (int*)&(el[0]);
  
  return ret;
}

template <>
Ng_Element Ng_GetElement<2> (int nr)
{
  const Element2d & el = mesh->SurfaceElement (SurfaceElementIndex (nr));
  
  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  ret.npoints = el.GetNP();
  ret.points = (int*)&el[0];
  return ret;
}

template <>
Ng_Element Ng_GetElement<3> (int nr)
{
  const Element & el = mesh->VolumeElement (ElementIndex (nr));
  
  Ng_Element ret;
  ret.type = NG_ELEMENT_TYPE(el.GetType());
  ret.npoints = el.GetNP();
  ret.points = (int*)&el[0];
  return ret;
}


