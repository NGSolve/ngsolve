#ifndef FILE_MESHACCESS
#define FILE_MESHACCESS

/*********************************************************************/
/* File:   meshaccess.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#include <nginterface.h>
#include <nginterface_v2.hpp>

namespace ngfem
{
  class ElementTransformation;
  class IntegrationPoint;
}



namespace ngcomp
{
  using netgen::Ng_Node;
  
  class MeshAccess;
  class Ngs_Element;


  enum VorB { VOL, BND };
  inline void operator++(VorB & vb, int)  { vb = VorB(vb+1); } 
  inline ostream & operator<< (ostream & ost, VorB vb)
  {
    if (vb == VOL) ost << "VOL"; else ost << "BND";
    return ost;
  }

  class ElementId
  {
    VorB vb;
    int nr;
  public:
    ElementId (VorB avb, int anr) : vb(avb), nr(anr) { ; }
    ElementId (int anr) : vb(VOL), nr(anr) { ; }
    int Nr() const { return nr; }
    // operator int () const { return nr; }
    bool IsVolume() const { return vb == VOL; }
    bool IsBoundary() const { return vb == BND; }

    bool operator< (int nr2) { return nr < nr2; }
    ElementId operator++ (int) { return ElementId(vb,nr++); }
    ElementId operator++ () { return ElementId(vb,++nr); }
    ElementId operator*() const { return *this; }
    bool operator!=(ElementId id2) const { return nr != id2.nr || vb != id2.vb; }
  };

  inline ostream & operator<< (ostream & ost, ElementId id)
  {
    return ost << (id.IsVolume() ? "VEl" : "BEl") << ' ' << id.Nr();
  }


  class Ngs_Element : public netgen::Ng_Element, public ElementId
  {
  public:
    Ngs_Element (const netgen::Ng_Element & el, ElementId id) 
      : netgen::Ng_Element(el), ElementId(id) { ; }
    auto Vertices() const -> decltype (ArrayObject(vertices)) { return vertices; }
    auto Edges() const -> decltype (ArrayObject(edges)) { return edges; }
    auto Faces() const -> decltype (ArrayObject(faces)) { return faces; }
    
    /*
      Converts element-type from Netgen to element-types of NGSolve.
      E.g. the Netgen-types NG_TRIG and NG_TRIG6 are merged to NGSolve type ET_TRIG.
    */
    static INLINE ELEMENT_TYPE ConvertElementType (NG_ELEMENT_TYPE type)
    {
      switch (type)
        {
        case NG_PNT:                    return ET_POINT;
        case NG_SEGM: case NG_SEGM3:    return ET_SEGM;
        case NG_TRIG: case NG_TRIG6:    return ET_TRIG;
        case NG_QUAD: case NG_QUAD6:    return ET_QUAD;
        case NG_TET: case NG_TET10:     return ET_TET;
        case NG_PRISM: case NG_PRISM12: return ET_PRISM;
        case NG_PYRAMID:                return ET_PYRAMID;
        case NG_HEX:                    return ET_HEX;
        default:
#ifndef __CUDA_ARCH__
	  throw Exception ("Netgen2NgS type conversion: Unhandled element type");
#else
	  return ET_POINT;
#endif	  
        }
    }

    ELEMENT_TYPE GetType () const 
    { return ConvertElementType (Ng_Element::GetType()); }
  };




  class ElementIterator
  {
    const MeshAccess & ma;
    VorB vb;
    int nr;
  public:
    ElementIterator (const MeshAccess & ama, VorB avb, int anr) 
      : ma(ama), vb(avb), nr(anr) { ; }
    ElementIterator operator++ () { return ElementIterator(ma, vb,++nr); }
    // ElementId operator*() const { return ElementId(vb,nr); }
    // Ngs_Element operator*() const { return ma[ElementId(vb,nr)]; }
    INLINE Ngs_Element operator*() const;
    bool operator!=(ElementIterator id2) const { return nr != id2.nr || vb != id2.vb; }
  };
  
  class ElementRange : public IntRange
  {
    const MeshAccess & ma;
    VorB vb;
    // IntRange r;
  public:
    ElementRange (const MeshAccess & ama, VorB avb, IntRange ar) 
      : IntRange(ar), ma(ama), vb(avb) { ; }
    ElementId First() const { return ElementId(vb, IntRange::First()); }
    ElementIterator begin () const { return ElementIterator(ma, vb,IntRange::First()); }
    ElementIterator end () const { return ElementIterator(ma, vb,IntRange::Next()); }
    ElementId operator[] (int nr) { return ElementId(vb, IntRange::First()+nr); }
  };

  template <VorB VB>
  class TElementIterator
  {
    const MeshAccess & ma;
    int nr;
  public:
    INLINE TElementIterator (const MeshAccess & ama, int anr) : ma(ama), nr(anr) { ; }
    INLINE TElementIterator operator++ () { return TElementIterator(ma, ++nr); }
    // ElementId operator*() const { return ElementId(VB,nr); }
    INLINE Ngs_Element operator*() const; 
    INLINE bool operator!=(TElementIterator id2) const { return nr != id2.nr; }
  };
  
  template <VorB VB>
  class TElementRange
  {
    const MeshAccess & ma;
    IntRange r;
  public:
    INLINE TElementRange (const MeshAccess & ama, IntRange ar) : ma(ama), r(ar) { ; }
    INLINE TElementIterator<VB> begin () const { return TElementIterator<VB>(ma, r.First()); }
    INLINE TElementIterator<VB> end () const { return TElementIterator<VB>(ma, r.Next()); }
  };
    

  /** 
      Access to mesh topology and geometry.

      MeshAccess provides information such as element-to-vertex table,
      elemenet-to-edge table etc. 
      
      It provides also element-geometry (straight or curved elements) 
      via GetTrafo. 

      Internally, MeshAccess calls functions from Netgen.
  */

  class GridFunction;

  class NGS_DLL_HEADER MeshAccess : public BaseStatusHandler
  {
    netgen::Ngx_Mesh mesh;

    /// buffered global quantities:
    /// dimension of the domain. Set to -1 if no mesh is present
    int dim;
  
    /// number of vertex, edge, face, and cell nodes
    int nnodes[4];

    // number of nodes of co-dimension i 
    // these are NC, NF, NE, NV  in 3D,
    // and NF, NE, NV, undef, in 2D
    int nnodes_cd[4];


    /// number of elements of dimension i
    int nelements[4];  
    /// number of elements of co-dimension i
    int nelements_cd[4];
    ///
    int ne_vb[2];  // index with VorB
    /// number of multigrid levels 
    int nlevels;

    /// max domain index
    int ndomains;

    /// max boundary index
    int nboundaries;

    /// for ALE
    shared_ptr<GridFunction> deformation;  

  public:
    /// connects to Netgen - mesh
    MeshAccess (shared_ptr<netgen::Mesh> amesh = NULL);
    /// loads mesh from file
    MeshAccess (string filename)
      : MeshAccess() 
    {
      LoadMesh (filename);
    }
    /// select this mesh in netgen visuaization
    void SelectMesh() const;
    /// not much to do 
    virtual ~MeshAccess ();

    /// the spatial dimension of the mesh
    int GetDimension() const { return dim; }  

    /// number of points. needed for isoparametric nodal elements
    int GetNP() const  { return Ng_GetNP(); }

    /// number of vertices
    int GetNV() const  { return nnodes[0]; }  

    /// number of elements in the domain
    int GetNE() const  { return nelements_cd[0]; }  

    /// number of boundary elements
    int GetNSE() const { return nelements_cd[1]; }  

    /// number of volume or boundary elements
    int GetNE(VorB vb) const { return ne_vb[vb]; } 

    /// number of edges in the whole mesh
    int GetNEdges() const { return nnodes[1]; }     

    /// number of faces in the whole mesh
    int GetNFaces() const { return nnodes[2]; }    


    /// maximal sub-domain (material) index. range is [0, ndomains)
    int GetNDomains () const  { return ndomains; }

    /// maximal boundary condition index. range is [0, nboundaries)
    int GetNBoundaries () const { return nboundaries; }

    /// returns point coordinate
    template <int D>
    void GetPoint (int pi, Vec<D> & p) const
    { 
      auto pt = mesh.GetPoint (pi);
      for (int j = 0; j < D; j++) p(j) = pt[j];
    }

    /// returns point coordinate
    template <int D>
    Vec<D> GetPoint (int pi) const
    { 
      Vec<D> p;
      auto pt = mesh.GetPoint (pi);
      for (int j = 0; j < D; j++) p(j) = pt[j];
      return p;
    }

    ElementRange Elements (VorB vb = VOL) const
    {
      return ElementRange (*this, vb, IntRange (0, GetNE(vb)));
    }

    template <VorB VB>
    TElementRange<VB> Elements () const
    {
      return TElementRange<VB> (*this, IntRange (0, GetNE(VB)));
    }

    /// the geometry type of the element
    ELEMENT_TYPE GetElType (int elnr) const
    { return GetElement(elnr).GetType(); }

    /// the geometry type of the boundary element
    ELEMENT_TYPE GetSElType (int elnr) const
    { return GetSElement(elnr).GetType(); }

    ELEMENT_TYPE GetElType (ElementId ei) const
    { return GetElement(ei).GetType(); }


    /// the sub-domain index of the element
    int GetElIndex (int elnr) const
    { 
      switch (dim)
        {
        case 1: return mesh.GetElementIndex<1>(elnr) - 1;
        case 2: return mesh.GetElementIndex<2>(elnr) - 1;
        case 3: default:
          return mesh.GetElementIndex<3>(elnr) - 1;
        }
    }

    /// the boundary-condition index of the boundary element
    int GetSElIndex (int elnr) const
    { 
      switch (dim)
        {
        case 1: return mesh.GetElementIndex<0>(elnr) - 1;
        case 2: return mesh.GetElementIndex<1>(elnr) - 1;
        case 3: default:
          return mesh.GetElementIndex<2>(elnr) - 1;
        }
    }

    int GetElIndex (ElementId ei) const
    {
      if (ei.IsVolume()) return GetElIndex (ei.Nr());
      else return GetSElIndex (ei.Nr());
    }


    /// change sub-domain index (???) 
    // void SetElIndex (int elnr, int index) const
    // { Ng_SetElementIndex (elnr+1,index+1); }

    /// the material of the element
    string GetElMaterial (int elnr) const
    { return Ng_GetElementMaterial (elnr+1); }

    /// the material of the sub-domain
    string GetDomainMaterial (int domnr) const
    { return Ng_GetDomainMaterial (domnr+1); }
      

    /// the boundary condition name of surface element
    string GetSElBCName (int selnr) const
    { return Ng_GetSurfaceElementBCName (selnr+1); }

    /// the boundary condition name of boundary condition number
    string GetBCNumBCName (int bcnr) const
    { return Ng_GetBCNumBCName (bcnr); }



    /// not sure who needs that
    int GetSElSurfaceNumber (const int elnr) const
    { return Ng_GetSurfaceElementSurfaceNumber (elnr+1)-1; }

    /// not sure who needs that
    int GetSElFDNumber (const int elnr) const
    { return Ng_GetSurfaceElementFDNumber (elnr+1)-1; }



    /// the sub-domain indices next to boundary element. 
    /// returns -1 for void
    void GetSElNeighbouringDomains(const int elnr, int & in, int & out) const
    { 
      Ng_GetSurfaceElementNeighbouringDomains(elnr+1,in,out);
      //in--; out--;
    }
  

    /// update buffered global quantities.
    /// Must be called after every change of the mesh topology
    void UpdateBuffers();



    /*
      Nodes are an abstraction for vertices, edges, faces, and cells
    */

    /// number of elements of dimension dim
    int GetNElements (int dim) const { return nelements[dim]; }

    /// number of nodes of type nt
    int GetNNodes (NODE_TYPE nt) const { return nnodes[nt]; }  


    /// the topology of a domain - element
    /// experimental, not recommended for use yet
    // void GetTopologicElement (int elnr, TopologicElement & topel) const;


    /**
       returns topology of a Netgen - element.  This is the new
       (2008), unified concept. The function returns a direct access
       to the Netgen mesh structure instead of copying point numbers
       etc. The nasty 1-0 conversion is done on the fly.
     */
    Ngs_Element GetElement (int elnr, bool boundary = 0) const
    {
      switch (dim-boundary)
	{
	case 1:	return Ngs_Element (mesh.GetElement<1> (elnr), 
                                    ElementId(boundary ? BND : VOL, elnr));
	case 2: return Ngs_Element (mesh.GetElement<2> (elnr), 
                                    ElementId(boundary ? BND : VOL, elnr));
	case 3:
        default: return Ngs_Element (mesh.GetElement<3> (elnr), 
                                     ElementId(boundary ? BND : VOL, elnr));
	}
    }

    Ngs_Element GetElement (ElementId ei) const
    {
      int hdim = dim;
      if (ei.IsBoundary()) hdim--;
      switch (hdim)
	{
	case 1:	return Ngs_Element (mesh.GetElement<1> (ei.Nr()), ei);
	case 2: return Ngs_Element (mesh.GetElement<2> (ei.Nr()), ei);
	case 3:
        default: return Ngs_Element (mesh.GetElement<3> (ei.Nr()), ei);
	}
    }

    Ngs_Element operator[] (ElementId ei) const    
    {
      return GetElement (ei);
    }


    /**
       returns topology of a Netgen - element.  This is the new
       (2008), unified concept. The function returns a direct access
       to the Netgen mesh structure instead of copying point numbers
       etc. The nasty 1-0 convertion is done on the fly.
     */
    Ngs_Element GetSElement (int elnr) const
    {
      switch (dim)
	{
	case 1:	return Ngs_Element (mesh.GetElement<0> (elnr), ElementId(BND,elnr));
	case 2: return Ngs_Element (mesh.GetElement<1> (elnr), ElementId(BND,elnr));
	case 3: 
        default: return Ngs_Element (mesh.GetElement<2> (elnr), ElementId(BND,elnr));
	}
    }

    /**
       returns element of compile-time fixed dimension
     */
    template <int DIM, VorB vb>
    Ngs_Element GetElement (int elnr) const
    {
      return Ngs_Element (mesh.GetElement<DIM> (elnr), ElementId(vb, elnr));
    }


    void SetDeformation (shared_ptr<GridFunction> def)
    {
      deformation = def;
    }
    
    /**
       returns node topology.
       A facet or edge-node knows its vertices etc.
       The method is not yet fully functional.
     */
    template <int DIM>
    netgen::Ng_Node<DIM> GetNode (int nr) const
    {
      return mesh.GetNode<DIM> (nr);
    }

    /// returns the points of an element.
    /// vertices and possibly edge-midpoints
    void GetElPNums (int elnr, Array<int> & pnums) const
    { pnums = ArrayObject (GetElement(elnr).points); }

    /// returns the points of a boundary element.
    /// vertex and possibly edge-midpoints
    void GetSElPNums (int selnr, Array<int> & pnums) const
    { pnums = ArrayObject (GetSElement(selnr).points); }

    void GetElPNums (ElementId ei, Array<int> & pnums) const
    { pnums = ArrayObject (GetElement(ei).points); }


    /// returns the vertices of an element
    void GetElVertices (int elnr, Array<int> & vnums) const
    { vnums = ArrayObject (GetElement(elnr).vertices); }
    
    ///
    void GetElVertices (ElementId ei, Array<int> & vnums) const
    { vnums = GetElement(ei).Vertices(); }

    /// returns the vertices of a boundary element
    void GetSElVertices (int selnr, Array<int> & vnums) const
    { vnums = GetSElement(selnr).Vertices(); }

    /// returns the edges of an element
    void GetElEdges (int elnr, Array<int> & ednums) const
    { ednums = GetElement(elnr).Edges(); }

    // returns edge numbers and edge orientation of an element. (old style function)
    void GetElEdges (int elnr, Array<int> & ednums, Array<int> & orient) const;

    /// returns the edges of a boundary element
    void GetSElEdges (int selnr, Array<int> & ednums) const
    { ednums = ArrayObject (GetSElement(selnr).edges); }

    // returns edge numbers and edge orientation of an element. (old style function)
    void GetSElEdges (int selnr, Array<int> & ednums, Array<int> & orient) const;

    /// returns the faces of an element
    void GetElFaces (int elnr, Array<int> & fnums) const
    { fnums = GetElement(elnr).Faces(); }

    // returns face numbers and face orientation of an element. (old style function)
    void GetElFaces (int elnr, Array<int> & fnums, Array<int> & orient) const;

    /// returns face number of surface element
    int GetSElFace (int selnr) const;

    // returns face number and orientation of surface element
    void GetSElFace (int selnr, int & fnum, int & orient) const;

    /// returns vertex numbers of face
    void GetFacePNums (int fnr, Array<int> & pnums) const;
    /// returns vertex numbers of edge
    void GetEdgePNums (int enr, int & pn1, int & pn2) const;
    /// returns vertex numbers of edge
    void GetEdgePNums (int enr, Array<int> & pnums) const;
    /// returns all elements connected to an edge
    void GetEdgeElements (int enr, Array<int> & elnums) const;
    /// returns all elements connected to an edge
    void GetEdgeSurfaceElements (int enr, Array<int> & elnums) const;
    /// returns all edges of a face
    void GetFaceEdges (int fnr, Array<int> & edges) const;
    /// returns elements connected to a face
    void GetFaceElements (int fnr, Array<int> & elnums) const;
    /// returns surface elements connected to a face
    void GetFaceSurfaceElements (int fnr, Array<int> & elnums) const;
    /// point numbers of a 1D element
    // void GetSegmentPNums (int snr, Array<int> & pnums) const;
    /// index of 1D element
    // int GetSegmentIndex (int snr) const;
    
    /// number of facets of an element. 
    /// facets are edges (2D) or faces (3D)
    int GetNFacets() const { return nnodes_cd[1]; } 
    /// facets of an element
    void GetElFacets (int elnr, Array<int> & fnums) const;
    /// facet of a surface element
    void GetSElFacets (int selnr, Array<int> & fnums) const;
    /// vertices of a facet
    void GetFacetPNums (int fnr, Array<int> & pnums) const;
    /// geometry type of facet
    ELEMENT_TYPE GetFacetType (int fnr) const;
    /// elements connected to facet
    void GetFacetElements (int fnr, Array<int> & elnums) const
    {
      switch (dim)
        {
        case 1: GetVertexElements (fnr, elnums); break;
        case 2: GetEdgeElements (fnr, elnums); break;
        case 3: GetFaceElements (fnr, elnums); break;
        }
    }
    void GetFacetSurfaceElements (int fnr, Array<int> & elnums) const
    {
      switch (dim)
        {
          // case 1: GetVertexElements (fnr, elnums); break;
        case 2: GetEdgeSurfaceElements (fnr, elnums); break;
        case 3: GetFaceSurfaceElements (fnr, elnums); break;
        }
    }

    // void GetVertexElements (int vnr, Array<int> & elnrs) const;
    /// element order stored in Netgen
    int GetElOrder (int enr) const
    { return Ng_GetElementOrder (enr+1); } 
    /// anisotropic order stored in Netgen
    INT<3> GetElOrders (int enr) const
    { 
      INT<3> eo; 
      Ng_GetElementOrders(enr+1,&eo[0],&eo[1],&eo[2]); 
      return eo; 
    } 
    /// set element order in Netgen
    void SetElOrder (int enr, int order) const
    { Ng_SetElementOrder (enr+1,order); }
    /// set anisotropic element order in Netgen
    void SetElOrders (int enr, int ox, int oy, int oz) const
    { Ng_SetElementOrders (enr+1, ox,oy,oz); }
    /// order of suface element
    int GetSElOrder (int enr) const
    { return Ng_GetSurfaceElementOrder (enr+1); } 
    /// anisotropic order of suface element
    INT<2> GetSElOrders (int enr) const
    { 
      INT<2> eo; 
      Ng_GetSurfaceElementOrders(enr+1,&eo[0],&eo[1]); 
      return eo; 
    }
    /// set surface element order
    void SetSElOrder (int enr, int order) const
    { Ng_SetSurfaceElementOrder (enr+1,order); }
    /// set anisotropic surface element order
    void SetSElOrders (int enr, int ox, int oy) const
    { Ng_SetSurfaceElementOrders (enr+1, ox,oy); }
    

    /// the volume of an element (mid-point rule)
    double ElementVolume (int elnr) const;
    /// the area of a boundary element (mid-point rule)
    double SurfaceElementVolume (int selnr) const;
      

    /// number of multigrid levels
    int GetNLevels() const { return nlevels; }  
    /// the two parent vertices of a vertex. -1 for coarse-grid vertices
    void GetParentNodes (int pi, int * parents) const
    { 
      Ng_GetParentNodes (pi+1, parents);
      parents[0]--; parents[1]--; 
    }
    /// number of parent element on next coarser mesh
    int GetParentElement (int ei) const
    { return Ng_GetParentElement (ei+1)-1; }
    /// number of parent boundary element on next coarser mesh
    int GetParentSElement (int ei) const
    { return Ng_GetParentSElement (ei+1)-1; }
  
    /// representant of vertex for anisotropic meshes
    int GetClusterRepVertex (int pi) const
    { return Ng_GetClusterRepVertex (pi+1)-1; }
    /// representant of edge for anisotropic meshes
    int GetClusterRepEdge (int pi) const
    { return Ng_GetClusterRepEdge (pi+1)-1; }
    /// representant of face for anisotropic meshes
    int GetClusterRepFace (int pi) const
    { return Ng_GetClusterRepFace (pi+1)-1; }
    /// representant of element for anisotropic meshes
    int GetClusterRepElement (int pi) const
    { return Ng_GetClusterRepElement (pi+1)-1; }
    
    
    /// returns the transformation from the reference element to physical element.
    /// Given a point in the refrence element, the ElementTransformation can 
    /// compute the mapped point as well as the Jacobian
    ngfem::ElementTransformation & GetTrafo (int elnr, bool boundary, LocalHeap & lh) const;
    
    ngfem::ElementTransformation & GetTrafo (ElementId ei, LocalHeap & lh) const
    {
      return GetTrafo (ei.Nr(), ei.IsBoundary(), lh);
    }


    // (old style optimization)
    void SetPointSearchStartElement(const int el) const;


    
    int FindElementOfPoint (FlatVector<double> point,
			    IntegrationPoint & ip, 
			    bool build_searchtree,
			    const Array<int> * const indices = NULL) const;
    int FindElementOfPoint (FlatVector<double> point,
			    IntegrationPoint & ip, 
			    bool build_searchtree,
			    int index) const;
    int FindSurfaceElementOfPoint (FlatVector<double> point,
				   IntegrationPoint & ip, 
				   bool build_searchtree,
				   const Array<int> * const indices = NULL) const;
    int FindSurfaceElementOfPoint (FlatVector<double> point,
				   IntegrationPoint & ip, 
				   bool build_searchtree,
				   int index) const;

    /// is element straight or curved ?
    bool IsElementCurved (int elnr) const
    { return bool (Ng_IsElementCurved (elnr+1)); }
    
    
    void GetPeriodicVertices ( Array<ngstd::INT<2> > & pairs) const;
    int GetNPairsPeriodicVertices () const;
    void GetPeriodicVertices (int idnr, Array<ngstd::INT<2> > & pairs) const;
    int GetNPairsPeriodicVertices (int idnr) const;  

    void GetPeriodicEdges ( Array<ngstd::INT<2> > & pairs) const;
    int GetNPairsPeriodicEdges () const;
    void GetPeriodicEdges (int idnr, Array<ngstd::INT<2> > & pairs) const;
    int GetNPairsPeriodicEdges (int idnr) const;  


    virtual void PushStatus (const char * str) const;
    virtual void PopStatus () const;
    virtual void SetThreadPercentage (double percent) const;
    virtual void GetStatus (string & str, double & percent) const;

    virtual void SetTerminate(void) const;
    virtual void UnSetTerminate(void) const;
    virtual bool ShouldTerminate(void) const;
  
    ///// Added by Roman Stainko ....
    void GetVertexElements( int vnr, Array<int>& elems) const;

    void GetVertexSurfaceElements( int vnr, Array<int>& elems) const;


  private:
    Array<bool> higher_integration_order;
  public:
    void SetHigherIntegrationOrder(int elnr);
    void UnSetHigherIntegrationOrder(int elnr);

    void LoadMesh (const string & filename);
    void LoadMesh (istream & str);
    void SaveMesh (ostream & str) const;
    void ArchiveMesh (Archive & archive);
    // void LoadMeshFromString(const string & str);

    // void PrecomputeGeometryData(int intorder);

    void InitPointCurve(double red = 1, double green = 0, double blue = 0) const;
    void AddPointCurvePoint(const Vec<3> & point) const;


    
    template <int DIMS, int DIMR> friend class Ng_ElementTransformation;
    template <int DIMS, int DIMR> friend class Ng_ConstElementTransformation;


    MPI_Comm GetCommunicator () const { return ngs_comm; }

    /**
       Returns the list of other MPI - processes where node is present.
       The ordering coincides for all processes.
    */
    void GetDistantProcs (Node node, Array<int> & procs) const;

    /**
       Returns the global number of the node.
       Currently, this function works only for vertex-nodes.
     */
    int GetGlobalNodeNum (Node node) const;
    
    
    FlatArray<int> GetDistantProcs (Node node) const
    {
#ifdef PARALLEL
      std::tuple<int,int*> tup = mesh.GetDistantProcs(node.GetType(), node.GetNr());
      return FlatArray<int> (std::get<0>(tup), std::get<1>(tup));
#else
      return FlatArray<int>(0,nullptr);
#endif
    }


    /// Reduces MPI - distributed data associated with mesh-nodes
    template <typename T>
    void AllReduceNodalData (NODE_TYPE nt, Array<T> & data, MPI_Op op) const;
  };





  INLINE Ngs_Element ElementIterator :: operator*() const
  {
    return ma[ElementId(vb,nr)]; 
  }
  
  template <VorB VB>
  INLINE Ngs_Element TElementIterator<VB>::operator*() const 
  {
    return ma[ElementId(VB,nr)]; 
  }




  /**
     Controls the progress - output.
     It controls the Netgen - progressbar as well as console progress update.
     In parallel, all processes must enter and call the Done method.
   */
  class ProgressOutput
  {
    shared_ptr<MeshAccess> ma;
    string task;
    int total;
    double prevtime;
    bool is_root;

    int cnt;
    bool done_called;
  public:
    NGS_DLL_HEADER ProgressOutput (shared_ptr<MeshAccess> ama,
		    string atask, int atotal);
    NGS_DLL_HEADER ~ProgressOutput();

    NGS_DLL_HEADER void Update();
    NGS_DLL_HEADER void Update(int nr);
    NGS_DLL_HEADER void Done();
  };



  /// old style, compatibility for a while
  template <typename T>
  void AllReduceNodalData (NODE_TYPE nt, Array<T> & data, MPI_Op op, const MeshAccess & ma)
  {
    ma.AllReduceNodalData (nt, data, op);
  }


#ifdef PARALLEL

  template <typename T>
  void MeshAccess::
  AllReduceNodalData (NODE_TYPE nt, Array<T> & data, MPI_Op op) const
  {
    MPI_Comm comm = GetCommunicator();
    MPI_Datatype type = MyGetMPIType<T>();

    int ntasks = MyMPI_GetNTasks (comm);
    if (ntasks <= 1) return;

    Array<int> cnt(ntasks), distp;
    cnt = 0;
    for (int i = 0; i < GetNNodes(nt); i++)
      {
	GetDistantProcs (Node (nt, i), distp);
	for (int j = 0; j < distp.Size(); j++)
	  cnt[distp[j]]++;
      }

    Table<T> dist_data(cnt), recv_data(cnt);

    cnt = 0;
    for (int i = 0; i < GetNNodes(nt); i++)
      {
	GetDistantProcs (Node (nt, i), distp);
	for (int j = 0; j < distp.Size(); j++)
	  dist_data[distp[j]][cnt[distp[j]]++] = data[i];
      }

    Array<MPI_Request> requests;
    for (int i = 0; i < ntasks; i++)
      if (cnt[i])
	{
	  requests.Append (MyMPI_ISend (dist_data[i], i, MPI_TAG_SOLVE, comm));
	  requests.Append (MyMPI_IRecv (recv_data[i], i, MPI_TAG_SOLVE, comm));
	}
    MyMPI_WaitAll (requests);
    
    cnt = 0;
    for (int i = 0; i < data.Size(); i++)
      {
	GetDistantProcs (Node (nt, i), distp);
	for (int j = 0; j < distp.Size(); j++)
	  MPI_Reduce_local (&recv_data[distp[j]][cnt[distp[j]]++],
			    &data[i], 1, type, op);
      }
  }


#else

  template <typename T>
  void MeshAccess::
  AllReduceNodalData (NODE_TYPE nt, Array<T> & data, MPI_Op op) const { ; }

#endif

}

#endif
