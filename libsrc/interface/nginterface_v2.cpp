#include <meshing.hpp>

#ifdef SOCKETS
#include "../sockets/sockets.hpp"
#endif

#include "nginterface.h"
#include "nginterface_v2.hpp"
// #include <visual.hpp>

namespace netgen
{
#include "writeuser.hpp"
  extern shared_ptr<Mesh> mesh;
}


namespace netgen
{
#define NGX_INLINE
#include "nginterface_v2_impl.hpp"

  shared_ptr<Mesh> Ngx_Mesh :: SelectMesh () const
  {
    shared_ptr<Mesh> hmesh = netgen::mesh;

    netgen::mesh = mesh;
    SetGlobalMesh (mesh);

    return hmesh;
  }
  

  Ngx_Mesh :: Ngx_Mesh (shared_ptr<Mesh> amesh) 
  {
    if (amesh)
      mesh = amesh;
    else
      mesh = netgen::mesh;
  }
  
  Ngx_Mesh * LoadMesh (const string & filename)
  {
    netgen::mesh.reset();
    Ng_LoadMesh (filename.c_str());
    return new Ngx_Mesh (netgen::mesh);
  }

  void Ngx_Mesh :: LoadMesh (const string & filename)
  {
    netgen::mesh.reset();
    Ng_LoadMesh (filename.c_str());
    // mesh = move(netgen::mesh);
    mesh = netgen::mesh;
  }

  void Ngx_Mesh :: LoadMesh (istream & ist)
  {
    netgen::mesh = make_shared<Mesh>();
    netgen::mesh -> Load (ist);
    // mesh = move(netgen::mesh);
    mesh = netgen::mesh;
    SetGlobalMesh (mesh);
  }

  void Ngx_Mesh :: SaveMesh (ostream & ost) const
  {
    mesh -> Save (ost);
  }

  void Ngx_Mesh :: DoArchive (ngstd::Archive & archive)
  {

    if (archive.Output())
      {
        stringstream str;
        SaveMesh (str);
        string st = str.str();
        archive & st;
      }
    else
      {
        string st;
        archive & st;
        stringstream str(st);
        LoadMesh (str);
      }
  }

  void Ngx_Mesh :: UpdateTopology ()
  {
    if (mesh)
      mesh -> UpdateTopology();
  }


  /*
  Ngx_Mesh :: Ngx_Mesh (Mesh * amesh)
    : mesh(amesh)
  { ; }
  */

  Ngx_Mesh :: ~Ngx_Mesh ()
  {
    if (mesh == netgen::mesh)
      netgen::mesh = nullptr;
  }

  int Ngx_Mesh :: GetDimension() const
  {
    return mesh -> GetDimension();
  }

  int Ngx_Mesh :: GetNLevels() const
  {
    return mesh -> mglevels;
  }
  
  int Ngx_Mesh :: GetNElements (int dim) const
  {
    switch (dim)
      {
      case 0: return mesh -> pointelements.Size();
      case 1: return mesh -> GetNSeg();
      case 2: return mesh -> GetNSE();
      case 3: return mesh -> GetNE();
      }
    return -1;
  }
  
  int Ngx_Mesh :: GetNNodes (int nt) const
  {
    switch (nt)
      {
      case 0: return mesh -> GetNV();
      case 1: return mesh->GetTopology().GetNEdges();
      case 2: return mesh->GetTopology().GetNFaces();
      case 3: return mesh -> GetNE();
      }
    return -1;
  }

  /*
  Ng_Point Ngx_Mesh :: GetPoint (int nr) const
  {
    return Ng_Point (&mesh->Point(nr + PointIndex::BASE)(0));
  }
  */

  /*
  template <> DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<0> (int nr) const
  {
    const Element0d & el = mesh->pointelements[nr];
    
    Ng_Element ret;
    ret.type = NG_PNT;
    ret.index = el.index;

    ret.points.num = 1;
    ret.points.ptr = (int*)&el.pnum;

    ret.vertices.num = 1;
    ret.vertices.ptr = (int*)&el.pnum;

    ret.edges.num = 0;
    ret.edges.ptr = NULL;

    ret.faces.num = 0;
    ret.faces.ptr = NULL;

    return ret;
  }
  */
  /*
  template <> DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<1> (int nr) const
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

  template <> DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<2> (int nr) const
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

  template <> DLL_HEADER Ng_Element Ngx_Mesh :: GetElement<3> (int nr) const
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
  */

  
  /*
  template <>
  DLL_HEADER int Ngx_Mesh :: GetElementIndex<0> (int nr) const
  {
    return 0;
  }

  template <>
  DLL_HEADER int Ngx_Mesh :: GetElementIndex<1> (int nr) const
  {
    return (*mesh)[SegmentIndex(nr)].si;
  }
  
  template <>
  DLL_HEADER int Ngx_Mesh :: GetElementIndex<2> (int nr) const
  {
    int ind = (*mesh)[SurfaceElementIndex(nr)].GetIndex(); 
    return mesh->GetFaceDescriptor(ind).BCProperty();
  }
  
  template <>
  DLL_HEADER int Ngx_Mesh :: GetElementIndex<3> (int nr) const
  {
    return (*mesh)[ElementIndex(nr)].GetIndex();
  }
  */

















  /*
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

  template <> int DLL_HEADER Ng_GetNElements<0> ()
  {
    return 0;
  }

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

  template <> DLL_HEADER Ng_Element Ng_GetElement<0> (int nr)
  {
    cout << "Netgen does not support 0-D elements" << endl;
	Ng_Element ret;
	return ret;
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
  */







  
  
  template <> DLL_HEADER void Ngx_Mesh :: 
  ElementTransformation<3,3> (int elnr, 
                              const double * xi,
                              double * x, 
                              double * dxdxi) const
  {
    Point<3> xl(xi[0], xi[1], xi[2]);
    Point<3> xg;
    Mat<3,3> dx;
    mesh->GetCurvedElements().CalcElementTransformation (xl, elnr, xg, dx);

    if (x)
      for (int i = 0; i < 3; i++) x[i] = xg(i);
    
    if (dxdxi)
      for (int i=0; i<3; i++)
        {
          dxdxi[3*i] = dx(i,0);
          dxdxi[3*i+1] = dx(i,1);
          dxdxi[3*i+2] = dx(i,2);
        }
  }

  template <> DLL_HEADER void Ngx_Mesh ::
  ElementTransformation<2,3> (int elnr,
                              const double * xi, 
                              double * x, 
                              double * dxdxi) const
  {
    Point<2> xl(xi[0], xi[1]);
    Point<3> xg;
    Mat<3,2> dx;

    mesh->GetCurvedElements().CalcSurfaceTransformation (xl, elnr, xg, dx);
    
    if (x)
      for (int i = 0; i < 3; i++) x[i] = xg(i);
    
    if (dxdxi)
      for (int i=0; i<3; i++)
        {
          dxdxi[2*i] = dx(i,0);
          dxdxi[2*i+1] = dx(i,1);
        }
  }


  template <> DLL_HEADER void Ngx_Mesh ::
  ElementTransformation<2,2> (int elnr,
                                 const double * xi, 
                                 double * x, 
                                 double * dxdxi) const
  {
    Point<2> xl(xi[0], xi[1]);
    Point<3> xg;
    Mat<3,2> dx;

    mesh->GetCurvedElements().CalcSurfaceTransformation (xl, elnr, xg, dx);
    
    if (x)
      for (int i = 0; i < 2; i++) x[i] = xg(i);
    
    if (dxdxi)
      for (int i=0; i<2; i++)
        {
          dxdxi[2*i] = dx(i,0);
          dxdxi[2*i+1] = dx(i,1);
        }
  }




  template <> DLL_HEADER void Ngx_Mesh :: 
  ElementTransformation<1,2> (int elnr,
                              const double * xi,
                              double * x,
                              double * dxdxi) const
  {
    Point<3> xg;
    Vec<3> dx;

    mesh->GetCurvedElements().CalcSegmentTransformation (xi[0], elnr, xg, dx);
    
    if (x)
      for (int i = 0; i < 2; i++) x[i] = xg(i);
    
    if (dxdxi)
      for (int i=0; i < 2; i++)
        dxdxi[i] = dx(i);
  }



  template <> DLL_HEADER void Ngx_Mesh :: 
  ElementTransformation<1,1> (int elnr,
                              const double * xi,
                              double * x,
                              double * dxdxi) const
  {
    Point<3> xg;
    Vec<3> dx;

    mesh->GetCurvedElements().CalcSegmentTransformation (xi[0], elnr, xg, dx);
    
    if (x) x[0] = xg(0);
    if (dxdxi) dxdxi[0] = dx(0);
  }


  template <> DLL_HEADER void Ngx_Mesh :: 
  ElementTransformation<0,1> (int elnr, 
                              const double * xi,
                              double * x,
                              double * dxdxi) const
  {
    PointIndex pnum = mesh->pointelements[elnr].pnum;
    if (x) x[0] = (*mesh)[pnum](0);
    if (dxdxi) dxdxi[0] = 0;
  }









  
  
  template <> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<3,3> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointElementTransformation (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }
  
  template <> DLL_HEADER void Ngx_Mesh ::
  MultiElementTransformation<2,2> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSurfaceTransformation<2> (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }

  template <> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<2,3> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSurfaceTransformation<3> (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }

  template <> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<1,2> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    mesh->GetCurvedElements().CalcMultiPointSegmentTransformation<2> (elnr, npts, xi, sxi, x, sx, dxdxi, sdxdxi);
  }

  template <> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<1,1> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    for (int i = 0; i < npts; i++)
      ElementTransformation<1,1> (elnr, xi + i*sxi, x+i*sx, dxdxi+i*sdxdxi);
  }


  template <> DLL_HEADER void Ngx_Mesh :: 
  MultiElementTransformation<0,1> (int elnr, int npts,
                                   const double * xi, size_t sxi,
                                   double * x, size_t sx,
                                   double * dxdxi, size_t sdxdxi) const
  {
    for (int i = 0; i < npts; i++)
      ElementTransformation<0,1> (elnr, xi + i*sxi, x+i*sx, dxdxi+i*sdxdxi);
  }






  template <> DLL_HEADER int Ngx_Mesh :: GetNNodes<1> ()
  {
    return mesh->GetTopology().GetNEdges();
  }

  template <> DLL_HEADER int Ngx_Mesh :: GetNNodes<2> ()
  {
    return mesh->GetTopology().GetNFaces();
  }

  template <> DLL_HEADER Ng_Node<1> Ngx_Mesh :: GetNode<1> (int nr) const
  {
    Ng_Node<1> node;
    node.vertices.ptr = mesh->GetTopology().GetEdgeVerticesPtr(nr);
    return node;
  }

  template <> DLL_HEADER Ng_Node<2> Ngx_Mesh :: GetNode<2> (int nr) const
  {
    Ng_Node<2> node;
    node.vertices.ptr = mesh->GetTopology().GetFaceVerticesPtr(nr);
    node.vertices.nv = (node.vertices.ptr[3] == 0) ? 3 : 4;
    return node;
  }


  template <>
  DLL_HEADER int Ngx_Mesh :: FindElementOfPoint <1> 
  (double * hp, double * lami,
   bool build_searchtree, 
   int * const indices, int numind) const

  {
    if (mesh->GetDimension() != 1)
      throw NgException("FindElementOfPoint<1> called for multidim mesh");

    Point<3> p(hp[0], 0,0);
    for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
      {
        auto & seg = (*mesh)[si];
        Point<3> p1 = (*mesh)[seg[0]];
        Point<3> p2 = (*mesh)[seg[1]];
        double lam = (p(0)-p1(0)) / (p2(0)-p1(0));
        if (lam >= -1e-10 && lam <= 1+1e-10)
          {
            lami[0] = 1-lam;
            return si;
          }
      }
    return -1;
  }

  template <>
  DLL_HEADER int Ngx_Mesh :: FindElementOfPoint <2> 
  (double * p, double * lami,
   bool build_searchtree, 
   int * const indices, int numind) const

  {
    Array<int> dummy(numind);
    for (int i = 0; i < numind; i++) dummy[i] = indices[i]+1;
    
    double lam3[3];
    Point<3> p2d(p[0], p[1], 0);
    int ind = 
      mesh->GetElementOfPoint(p2d, lam3, &dummy, build_searchtree);
    
    if (ind > 0)
      {
        if(mesh->SurfaceElement(ind).GetType()==QUAD)
          {
            lami[0] = lam3[0];
            lami[1] = lam3[1];
          }
        else 
          {
            lami[0] = 1-lam3[0]-lam3[1];
            lami[1] = lam3[0];
          }
      }
    return ind-1;
  }


  template <>
  DLL_HEADER int Ngx_Mesh :: FindElementOfPoint <3>
  (double * p, double * lami,
   bool build_searchtree, 
   int * const indices, int numind) const

  {
    Array<int> dummy(numind);
    for (int i = 0; i < numind; i++) dummy[i] = indices[i]+1;
    
    Point<3> p3d(p[0], p[1], p[2]);
    int ind = 
      mesh->GetElementOfPoint(p3d, lami, &dummy, build_searchtree);
    return ind-1;
  }


#ifdef PARALLEL
  
  std::tuple<int,int*>  Ngx_Mesh :: GetDistantProcs (int nodetype, int locnum) const
  {
    
    switch (nodetype)
      {
      case 0:
	{
	  FlatArray<int> dn = mesh->GetParallelTopology().GetDistantPNums(locnum);
	  return std::tuple<int,int*>(dn.Size(), &dn[0]);
	}
      case 1:
	{
	  FlatArray<int> dn = mesh->GetParallelTopology().GetDistantEdgeNums(locnum);
	  return std::tuple<int,int*>(dn.Size(), &dn[0]);
	}
      case 2:
	{
	  FlatArray<int> dn = mesh->GetParallelTopology().GetDistantFaceNums(locnum);
	  return std::tuple<int,int*>(dn.Size(), &dn[0]);
	}
      default:
	return std::tuple<int,int*>(0,nullptr);
      }
  }

#endif

}


int link_it_nginterface_v2;

