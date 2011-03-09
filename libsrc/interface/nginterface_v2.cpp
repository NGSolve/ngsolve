#include <mystdlib.h>
#include <meshing.hpp>



#ifdef SOCKETS
#include "../sockets/sockets.hpp"
#endif

#ifndef NOTCL
#include <visual.hpp>
#endif


#include "nginterface.h"
#include "nginterface_v2.hpp"



namespace netgen
{
#include "writeuser.hpp"

  extern AutoPtr<Mesh> mesh;
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
    mesh->GetCurvedElements().CalcMultiPointSegmentTransformation<2> (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
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


int link_it_nginterface_v2;

