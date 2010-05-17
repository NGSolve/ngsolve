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


namespace netgen
{



  template <> int DLL_HEADER Ng_GetNElements<1> ()
  {
    return mesh->GetNSeg();
  }

  template <> DLL_HEADER int Ng_GetNElements<2> ()
  {
    return mesh->GetNSE();
  }

  template <> DLL_HEADER int Ng_GetNElements<3> ()
  {
    return mesh->GetNE();
  }




  template <> DLL_HEADER Ng_Element Ng_GetElement<1> (int nr)
  {
    const Segment & el = mesh->LineSegment (SegmentIndex(nr));

    Ng_Element ret;
    ret.type = NG_ELEMENT_TYPE(el.GetType());

    ret.points.num = el.GetNP();
    ret.points.ptr = (int*)&(el[0]);

    ret.vertices.num = 2;
    ret.vertices.ptr = (int*)&(el[0]);

    ret.edges.num = 1;
    ret.edges.ptr = mesh->GetTopology().GetSegmentElementEdgesPtr (nr);

    ret.faces.num = 0;
    ret.faces.ptr = NULL;

    return ret;
  }

  template <> DLL_HEADER Ng_Element Ng_GetElement<2> (int nr)
  {
    const Element2d & el = mesh->SurfaceElement (SurfaceElementIndex (nr));
  
    Ng_Element ret;
    ret.type = NG_ELEMENT_TYPE(el.GetType());
    ret.points.num = el.GetNP();
    ret.points.ptr  = (int*)&el[0];

    ret.vertices.num = el.GetNV();
    ret.vertices.ptr = (int*)&(el[0]);

    ret.edges.num = MeshTopology::GetNEdges (el.GetType());
    ret.edges.ptr = mesh->GetTopology().GetSurfaceElementEdgesPtr (nr);

    ret.faces.num = MeshTopology::GetNFaces (el.GetType());
    ret.faces.ptr = mesh->GetTopology().GetSurfaceElementFacesPtr (nr);

    return ret;
  }

  template <> DLL_HEADER Ng_Element Ng_GetElement<3> (int nr)
  {
    const Element & el = mesh->VolumeElement (ElementIndex (nr));
  
    Ng_Element ret;
    ret.type = NG_ELEMENT_TYPE(el.GetType());
    ret.points.num = el.GetNP();
    ret.points.ptr = (int*)&el[0];

    ret.vertices.num = el.GetNV();
    ret.vertices.ptr = (int*)&(el[0]);

    ret.edges.num = MeshTopology::GetNEdges (el.GetType());
    ret.edges.ptr = mesh->GetTopology().GetElementEdgesPtr (nr);

    ret.faces.num = MeshTopology::GetNFaces (el.GetType());
    ret.faces.ptr = mesh->GetTopology().GetElementFacesPtr (nr);

    return ret;
  }


  DLL_HEADER Ng_Point Ng_GetPoint (int nr)
  {
    Ng_Point ret;
    ret.pt = &mesh->Point(nr + PointIndex::BASE)(0);
    return ret;
  }


  template <>
  DLL_HEADER int Ng_GetElementIndex<1> (int nr)
  {
    return (*mesh)[SegmentIndex(nr)].si;
  }
  
  template <>
  DLL_HEADER int Ng_GetElementIndex<2> (int nr)
  {
    int ind = (*mesh)[SurfaceElementIndex(nr)].GetIndex(); 
    return mesh->GetFaceDescriptor(ind).BCProperty();
  }
  
  template <>
  DLL_HEADER int Ng_GetElementIndex<3> (int nr)
  {
    return (*mesh)[ElementIndex(nr)].GetIndex();
  }
  
  
  template <>
  DLL_HEADER void Ng_MultiElementTransformation<3,3> (int elnr, int npts,
                                                      const double * xi, size_t sxi,
                                                      double * x, size_t sx,
                                                      double * dxdxi, size_t sdxdxi)
  {
    mesh->GetCurvedElements().CalcMultiPointElementTransformation (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }
  
  template <>
  DLL_HEADER void Ng_MultiElementTransformation<2,2> (int elnr, int npts,
                                                      const double * xi, size_t sxi,
                                                      double * x, size_t sx,
                                                      double * dxdxi, size_t sdxdxi)
  {
    mesh->GetCurvedElements().CalcMultiPointSurfaceTransformation<2> (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }

  template <>
  DLL_HEADER void Ng_MultiElementTransformation<2,3> (int elnr, int npts,
                                                      const double * xi, size_t sxi,
                                                      double * x, size_t sx,
                                                      double * dxdxi, size_t sdxdxi)
  {
    mesh->GetCurvedElements().CalcMultiPointSurfaceTransformation<3> (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }

  template <>
  DLL_HEADER void Ng_MultiElementTransformation<1,2> (int elnr, int npts,
                                                      const double * xi, size_t sxi,
                                                      double * x, size_t sx,
                                                      double * dxdxi, size_t sdxdxi)
  {
    for (int ip = 0; ip < npts; ip++)
      {
        Point<3> xg;
        Vec<3> dx;

        mesh->GetCurvedElements().CalcSegmentTransformation (xi[ip*sxi], elnr, xg, dx);
      
        if (x)
          for (int i = 0; i < 2; i++)
            x[ip*sx+i] = xg(i);
	  
        if (dxdxi)
          for (int i=0; i<2; i++)
            dxdxi[ip*sdxdxi+i] = dx(i);
      }
  }

  template <>
  DLL_HEADER void Ng_MultiElementTransformation<1,1> (int elnr, int npts,
                                                      const double * xi, size_t sxi,
                                                      double * x, size_t sx,
                                                      double * dxdxi, size_t sdxdxi)
  {
    cout << "1D not supported" << endl;
  }



  template <> DLL_HEADER int Ng_GetNNodes<1> ()
  {
    return mesh->GetTopology().GetNEdges();
  }

  template <> DLL_HEADER int Ng_GetNNodes<2> ()
  {
    return mesh->GetTopology().GetNFaces();
  }


  template <> DLL_HEADER Ng_Node<1> Ng_GetNode<1> (int nr)
  {
    Ng_Node<1> node;
    node.vertices.ptr = mesh->GetTopology().GetEdgeVerticesPtr(nr);
    return node;
  }



  template <> DLL_HEADER Ng_Node<2> Ng_GetNode<2> (int nr)
  {
    Ng_Node<2> node;
    node.vertices.ptr = mesh->GetTopology().GetFaceVerticesPtr(nr);
    node.vertices.nv = (node.vertices.ptr[3] == 0) ? 3 : 4;
    return node;
  }

}
