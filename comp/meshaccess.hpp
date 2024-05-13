#ifndef FILE_MESHACCESS
#define FILE_MESHACCESS

/*********************************************************************/
/* File:   meshaccess.hpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#include <nginterface_v2.hpp>
#include <core/ranges.hpp>

#include <elementtopology.hpp>

namespace ngfem
{
  class IntegrationPoint;
  class CoefficientFunction;
  class ElementTransformation;
}



namespace ngcomp
{
  class PML_Transformation;
  
  // using ngcore::INT;
  using netgen::Ng_Node;
  using ngfem::ELEMENT_TYPE;
  
  using namespace ngfem;
  
  class MeshAccess;
  class Ngs_Element;
  class Region;

  class Ngs_Element : public netgen::Ng_Element
  {
    ElementId ei;
  public:
    // static string defaultstring;
    Ngs_Element (const netgen::Ng_Element & el, ElementId id) 
      : netgen::Ng_Element(el), ei(id) { ; }
    AOWrapper<decltype(vertices)> Vertices() const { return vertices; }
    AOWrapper<decltype(points)> Points() const { return points; }
    AOWrapper<decltype(edges)> Edges() const { return edges; }
    AOWrapper<decltype(faces)> Faces() const { return faces; }
    AOWrapper<decltype(facets)> Facets() const { return facets; }
    // string_view GetMaterial() const { return mat ? *mat : defaultstring; }
    string_view GetMaterial() const { return mat; }
    operator ElementId() const { return ei; }
    auto VB() const { return ei.VB(); }
    auto Nr() const { return ei.Nr(); }
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
        case NG_QUAD: case NG_QUAD6: case NG_QUAD8: return ET_QUAD;
        case NG_TET:  case NG_TET10:     return ET_TET;
        case NG_PRISM: case NG_PRISM12: case NG_PRISM15: return ET_PRISM;
        case NG_PYRAMID: case NG_PYRAMID13: return ET_PYRAMID;
        case NG_HEX7:   return ET_HEXAMID;          
        case NG_HEX: case NG_HEX20:     return ET_HEX;
          // default:
        }
      __assume (false);
      return ET_POINT;   // for some compiler
    }

    ELEMENT_TYPE GetType () const 
    { return ConvertElementType (Ng_Element::GetType()); }
  };

  inline ostream & operator<< (ostream & ost, Ngs_Element & el)
  {
    ost << ElementId(el);
    return ost;
  }

  class ElementIterator
  {
    const MeshAccess & ma;
    ElementId ei;
  public:
    ElementIterator (const MeshAccess & ama, ElementId aei) : ma(ama), ei(aei) { ; }
    ElementIterator operator++ () { return ElementIterator(ma, ++ei); }
    INLINE Ngs_Element operator*() const;
    bool operator!=(ElementIterator id2) const { return ei != id2.ei; }
    bool operator==(ElementIterator id2) const { return ei == id2.ei; }
  };

  class ElementRange : public IntRange
  {
    const MeshAccess & ma;
    VorB vb;
  public:
    ElementRange (const MeshAccess & ama, VorB avb, IntRange ar) 
      : IntRange(ar), ma(ama), vb(avb) { ; } 
    ElementId First() const { return ElementId(vb, IntRange::First()); }
    ElementIterator begin () const { return ElementIterator(ma, ElementId(vb,IntRange::First())); }
    ElementIterator end () const { return ElementIterator(ma, ElementId(vb,IntRange::Next())); }
    ElementId operator[] (size_t nr) { return ElementId(vb, IntRange::First()+nr); }
    ElementRange Split(size_t nr, size_t tot) { return ElementRange(ma, vb, IntRange::Split(nr,tot)); }
  };

  template <VorB VB>
  class TElementIterator
  {
    const MeshAccess & ma;
    size_t nr;
  public:
    INLINE TElementIterator (const MeshAccess & ama, size_t anr) : ma(ama), nr(anr) { ; }
    INLINE TElementIterator operator++ () { return TElementIterator(ma, ++nr); }
    INLINE Ngs_Element operator*() const;  // implemented below, after MeshAccess
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

  template <VorB VB, int DIM>
  class DimElementIterator
  {
    const MeshAccess & ma;
    size_t nr;
  public:
    INLINE DimElementIterator (const MeshAccess & ama, size_t anr) : ma(ama), nr(anr) { ; }
    INLINE DimElementIterator operator++ () { return DimElementIterator(ma, ++nr); }
    INLINE Ngs_Element operator*() const;  // implemented below, after MeshAccess
    INLINE bool operator!=(DimElementIterator id2) const { return nr != id2.nr; }
  };
  
  template <VorB VB, int DIM>
  class DimElementRange
  {
    const MeshAccess & ma;
    IntRange r;
  public:
    INLINE DimElementRange (const MeshAccess & ama, IntRange ar) : ma(ama), r(ar) { ; }
    INLINE auto begin () const { return DimElementIterator<VB,DIM>(ma, r.First()); }
    INLINE auto end () const { return DimElementIterator<VB,DIM>(ma, r.Next()); }
  };

  /*
  class NodeIterator
  {
    NodeId ni;
  public:
    NodeIterator (NodeId ani) : ni(ani) { ; }
    NodeIterator operator++ () { return NodeIterator(++ni); }
    INLINE NodeId operator*() const { return ni; }
    bool operator!=(NodeIterator id2) const { return ni != id2.ni; }
    bool operator==(NodeIterator id2) const { return ni == id2.ni; }
  };

  class NodeRange : public IntRange
  {
    NODE_TYPE nt;
  public:
    NodeRange (NODE_TYPE ant, IntRange ar) 
      : IntRange(ar), nt(ant) { ; } 
    NodeIterator begin () const { return NodeIterator(NodeId(nt,IntRange::First())); }
    NodeIterator end () const { return NodeIterator(NodeId(nt,IntRange::Next())); }
    NodeId operator[] (size_t nr) { return NodeId(nt, IntRange::First()+nr); }
  };
  */

  /** 
      Access to mesh topology and geometry.

      MeshAccess provides information such as element-to-vertex table,
      elemenet-to-edge table etc. 
      
      It provides also element-geometry (straight or curved elements) 
      via GetTrafo. 

      Internally, MeshAccess calls functions from Netgen.
  */

  class GridFunction;

  class NGS_DLL_HEADER MeshAccess : public enable_shared_from_this<MeshAccess>
  {
    netgen::Ngx_Mesh mesh;

    /// buffered global quantities:
    /// dimension of the domain. Set to -1 if no mesh is present
    int dim;
  
    /// number of vertex, edge, face, cell, element, facet, and global nodes
    size_t nnodes[7];

    // number of nodes of co-dimension i 
    // these are NC, NF, NE, NV  in 3D,
    // and NF, NE, NV, undef, in 2D
    size_t nnodes_cd[4];


    /// number of elements of dimension i
    size_t nelements[4];  
    /// number of elements of co-dimension i
    size_t nelements_cd[4];
    /// number of multigrid levels 
    int nlevels;

    /// number of regions of co-dimension i
    int nregions[4];

    //ngfem::ElementTransformation & GetTrafoDim (size_t elnr, Allocator & lh) const;
    typedef ngfem::ElementTransformation & (MeshAccess::*pfunc) (size_t elnr, Allocator & lh) const;    
    pfunc trafo_jumptable[4];

    int mesh_timestamp = -1; // timestamp of Netgen-mesh
    size_t timestamp = 0;
    
    /// for ALE
    shared_ptr<GridFunction> deformation;  

    /// pml trafos per sub-domain
    Array<shared_ptr <PML_Transformation>> pml_trafos;
    
    Array<std::tuple<int,int>> identified_facets;

    /// store periodic vertex mapping for each identification number
    // shared ptr because Meshaccess is copy constructible
    shared_ptr<Array<Array<IVec<2>>>> periodic_node_pairs[3] = {make_shared<Array<Array<IVec<2>>>>(),
                                                               make_shared<Array<Array<IVec<2>>>>(),
                                                               make_shared<Array<Array<IVec<2>>>>()};

    DynamicTable<size_t> neighbours[4][4];
    friend class Region;
  public:
    SimpleSignal updateSignal;

    /// for achiving ...
    MeshAccess ();
    /// connects to Netgen - mesh
    MeshAccess (shared_ptr<netgen::Mesh> amesh);
    /// loads new mesh from file
    MeshAccess (string filename, NgMPI_Comm amesh_comm = NgMPI_Comm{});
    /// select this mesh in netgen visuaization
    void SelectMesh() const;
    /// not much to do 
    virtual ~MeshAccess ();

    /// the spatial dimension of the mesh
    int GetDimension() const { return dim; }  

    /// number of points. needed for isoparametric nodal elements
    size_t GetNP() const  { return mesh.GetNP(); }

    /// number of vertices
    size_t GetNV() const  { return nnodes[0]; }  

    /// number of elements in the domain
    size_t GetNE() const  { return nelements_cd[0]; }  

    /// number of boundary elements
    size_t GetNSE() const { return nelements_cd[1]; }

    /// number of elements of co dimension 2
    size_t GetNCD2E() const { return nelements_cd[2]; }

    /// number of volume or boundary elements
    size_t GetNE(VorB vb) const { return nelements_cd[vb]; } 

    /// number of edges in the whole mesh
    size_t GetNEdges() const { return nnodes[1]; }     

    /// number of faces in the whole mesh
    size_t GetNFaces() const { return nnodes[2]; }    

    /// maximal sub-domain (material) index. range is [0, ndomains)
    int GetNDomains () const  { return nregions[VOL]; }

    /// maximal boundary condition index. range is [0, nboundaries)
    int GetNBoundaries () const { return nregions[BND]; }

    [[deprecated("Use GetNRegions (BBND) instead!")]]
    int GetNBBoundaries() const { return nregions[BBND]; }

    int GetNRegions (VorB vb) const { return nregions[vb]; }

    /// returns point coordinate
    template <int D>
    void GetPoint (size_t pi, Vec<D> & p) const
    { 
      auto pt = mesh.GetPoint (pi);
      for (int j = 0; j < D; j++) p(j) = pt[j];
    }

    /// returns point coordinate
    template <int D>
    Vec<D> GetPoint (size_t pi) const
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
    
    template <VorB VB, int DIM>
      auto Elements () const
    {
      return DimElementRange<VB,DIM> (*this, IntRange (0, GetNE(VB)));
    }

    auto Nodes (NODE_TYPE nt) const
    {
      return T_Range<NodeId> (NodeId(nt, 0), NodeId(nt, GetNNodes(nt)));
    }

    template <NODE_TYPE nt>
      auto Nodes () const
    {
      return T_Range<T_NodeId<nt>> (0, GetNNodes(nt));
    }

    auto Vertices() const { return Nodes<NT_VERTEX>(); }
    auto Edges() const { return Nodes<NT_EDGE>(); }
    auto Faces() const { return Nodes<NT_FACE>(); }
    auto Cells() const { return Nodes<NT_CELL>(); }

    template <typename TFUNC>
    void IterateElements (VorB vb, 
                          LocalHeap & clh, 
                          const TFUNC & func) const
    {
      if (task_manager)
        {
          SharedLoop sl(GetNE(vb));

          task_manager -> CreateJob
            ( [&] (const TaskInfo & ti) 
              {
                LocalHeap lh = clh.Split(ti.thread_nr, ti.nthreads);

                for (size_t mynr : sl)
                  {
                    HeapReset hr(lh);
                    ElementId ei(vb, mynr);
                    func (GetElement(ei), lh);
                  }
              } );
        }
      else
        {
          for (auto i : Range(GetNE(vb)))
            {
              HeapReset hr(clh);
              ElementId ei(vb, i);
	      func (GetElement(ei), clh);
            }
        }
    }

    template <typename TFUNC>
    void IterateElements (VorB vb,const TFUNC & func) const
    {
      if (task_manager)
        {
          SharedLoop sl(GetNE(vb));

          task_manager -> CreateJob
            ( [&] (const TaskInfo & ti)
              {
                for (size_t mynr : sl)
                  {
                    ElementId ei(vb, mynr);
                    func (GetElement(ei));
                  }
              } );
        }
      else
        {
          for (auto ei : Elements(vb))
            func (GetElement(ei));
        }
    }

    /// the geometry type of the element
    [[deprecated("Use GetElType with ElementId(VOL, elnr) instead!")]]            
    ELEMENT_TYPE GetElType (int elnr) const
      { return GetElement(ElementId(VOL,elnr)).GetType(); }

    /// the geometry type of the boundary element
    [[deprecated("Use GetElType with ElementId(BND, elnr) instead!")]]                
    ELEMENT_TYPE GetSElType (int elnr) const
      { return GetElement(ElementId(BND,elnr)).GetType(); }

    /// the geometry type of the element    
    ELEMENT_TYPE GetElType (ElementId ei) const
    { return GetElement(ei).GetType(); }


    /// the sub-domain index of the element
    [[deprecated("Use GetElIndex with ElementId(VOL, elnr) instead!")]]                    
    int GetElIndex (int elnr) const
    { 
      return GetElement(ElementId(VOL,elnr)).GetIndex();
    }

    /// the boundary-condition index of the boundary element
    [[deprecated("Use GetElIndex with ElementId(BND, elnr) instead!")]]                        
    int GetSElIndex (int elnr) const
    { 
      return GetElement(ElementId(BND,elnr)).GetIndex();
    }

    [[deprecated("Use GetElIndex with ElementId(BBND, elnr) instead!")]]                            
    int GetCD2ElIndex(int elnr) const
    {
      return GetElement(ElementId(BBND,elnr)).GetIndex();
    }

    int GetElIndex (ElementId ei) const
    {
      return GetElement(ei).GetIndex();
    }

    string_view GetMaterial(ElementId ei) const
    { return GetElement(ei).GetMaterial(); }
    
    // const string & GetMaterial(VorB vb, int region_nr) const
    string_view GetMaterial(VorB vb, int region_nr) const
    {
      switch (vb)
        {
        case VOL:   return mesh.GetMaterialCD<0> (region_nr);
        case BND:   return mesh.GetMaterialCD<1> (region_nr);
        case BBND:  return mesh.GetMaterialCD<2> (region_nr);
        case BBBND: return mesh.GetMaterialCD<3> (region_nr);
        default:    throw Exception("GetMaterial not implemented for " + ToString(vb));
        }
    }

    auto GetMaterials (VorB vb) const
    {
      return ArrayObject (GetNRegions(vb),
                          [this,vb] (size_t i) { return this->GetMaterial(vb, i); });
    }

    
    /// the material of the element
    [[deprecated("Use GetMaterial with ElementId(VOL, elnr) instead!")]]        
    string_view GetElMaterial (int elnr) const
    { return GetMaterial(ElementId(VOL, elnr)); }
      // { return Ng_GetElementMaterial (elnr+1); }

    /// the material of the sub-domain
    [[deprecated("Use GetMaterial(VOL, region_nr) instead!")]]                
    string_view GetDomainMaterial (int domnr) const
      { return GetMaterial(VOL, domnr); }
      // { return Ng_GetDomainMaterial (domnr+1); }
      

    /// the boundary condition name of surface element
    [[deprecated("Use GetMaterial with ElementId(BND, elnr) instead!")]]            
    string_view GetSElBCName (int selnr) const
    { return GetMaterial(ElementId(BND, selnr)); }      
      // { return Ng_GetSurfaceElementBCName (selnr+1); }

    /// the boundary condition name of boundary condition number
    [[deprecated("Use GetMaterial(BND, region_nr) instead!")]]            
    string_view GetBCNumBCName (int bcnr) const
      { return GetMaterial(BND, bcnr); }      
    // { return Ng_GetBCNumBCName (bcnr); }

    [[deprecated("Use GetMaterial(BBND, region_nr) instead!")]]                
    string_view GetCD2NumCD2Name (int cd2nr) const
      { return GetMaterial(BBND, cd2nr); }      
    // { return Ng_GetCD2NumCD2Name (cd2nr); }


    /// not sure who needs that
    int GetSElSurfaceNumber (const int elnr) const
    { return mesh.GetSurfaceElementSurfaceNumber (elnr+1)-1; }

    /// not sure who needs that
    int GetSElFDNumber (const int elnr) const
    { return mesh.GetSurfaceElementFDNumber (elnr+1)-1; }



    /// the sub-domain indices next to boundary element. 
    /// returns -1 for void
    void GetSElNeighbouringDomains(const int elnr, int & in, int & out) const;
  

    /// update buffered global quantities.
    /// Must be called after every change of the mesh topology
    void UpdateBuffers();



    /*
      Nodes are an abstraction for vertices, edges, faces, and cells
    */

    /// number of elements of dimension dim
    size_t GetNElements (int dim) const { return nelements[dim]; }

    /// number of nodes of type nt
    size_t GetNNodes (NODE_TYPE nt) const { return nnodes[nt]; }  


    /// the topology of a domain - element
    /// experimental, not recommended for use yet
    // void GetTopologicElement (int elnr, TopologicElement & topel) const;


    /**
       returns topology of a Netgen - element.  This is the new
       (2008), unified concept. The function returns a direct access
       to the Netgen mesh structure instead of copying point numbers
       etc. The nasty 1-0 conversion is done on the fly.
     */
    [[deprecated("Use GetElement (ElementId ei) instead!")]]        
    INLINE Ngs_Element GetElement (int elnr, bool boundary = 0) const
    {
      switch (dim-boundary)
	{
        case 0:	return Ngs_Element (mesh.GetElement<0> (elnr), 
                                    ElementId(boundary ? BND : VOL, elnr));
	case 1:	return Ngs_Element (mesh.GetElement<1> (elnr), 
                                    ElementId(boundary ? BND : VOL, elnr));
	case 2: return Ngs_Element (mesh.GetElement<2> (elnr), 
                                    ElementId(boundary ? BND : VOL, elnr));
	case 3:
        default: return Ngs_Element (mesh.GetElement<3> (elnr), 
                                     ElementId(boundary ? BND : VOL, elnr));
	}
    }
    
    INLINE Ngs_Element GetElement (ElementId ei) const
    {
      int hdim = dim - int(ei.VB());
      switch (hdim)
	{
        case 0:	return Ngs_Element (mesh.GetElement<0> (ei.Nr()), ei);
	case 1:	return Ngs_Element (mesh.GetElement<1> (ei.Nr()), ei);
	case 2: return Ngs_Element (mesh.GetElement<2> (ei.Nr()), ei);
	case 3:
        default: return Ngs_Element (mesh.GetElement<3> (ei.Nr()), ei);
	}
    }

    template <VorB VB, int DIM>
      INLINE Ngs_Element GetElement (T_ElementId<VB,DIM> ei) const
    {
      constexpr int HDIM = DIM - int(VB);
      return Ngs_Element (mesh.GetElement<HDIM> (ei.Nr()), ei);
    }

    
    INLINE Ngs_Element operator[] (ElementId ei) const    
    {
      return GetElement (ei);
    }


    /**
       returns topology of a Netgen - element.  This is the new
       (2008), unified concept. The function returns a direct access
       to the Netgen mesh structure instead of copying point numbers
       etc. The nasty 1-0 conversion is done on the fly.
     */
    [[deprecated("Use GetElement (ElementId (BND,nr)) instead!")]]            
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
    
    [[deprecated("Use GetElement (ElementId (BBND,nr)) instead!")]]            
    Ngs_Element GetCD2Element(int elnr) const
    {
      switch(dim)
	{
	case 1: throw Exception("No CoDim 2 Element for dimension 1");
	case 2: return Ngs_Element(mesh.GetElement<0>(elnr),ElementId(BBND,elnr));
	case 3:
	default: return Ngs_Element(mesh.GetElement<1>(elnr),ElementId(BBND,elnr));
	}
    }

    /**
       returns element of compile-time fixed dimension
     */
    template <int DIM, VorB vb>
      inline Ngs_Element GetElement (size_t elnr) const
    {
      return Ngs_Element (mesh.GetElement<DIM> (elnr), ElementId(vb, elnr));
    }

    auto GetTimeStamp() const { return timestamp; }
    
    void SetRefinementFlag (ElementId ei, bool ref)
    {
      switch (dim-int(ei.VB()))
        {
        case 2: mesh.SetRefinementFlag<2>(ei.Nr(), ref); break;
        case 3: mesh.SetRefinementFlag<3>(ei.Nr(), ref); break;
        default: ; 
        }
      /*      
      if (id.IsVolume())
        Ng_SetRefinementFlag (id.Nr()+1, ref);
      else
        Ng_SetSurfaceRefinementFlag (id.Nr()+1, ref);
      */
    }

    void Refine (bool onlyonce);
    void Curve (int order);
    int GetCurveOrder ();

    void EnableTable (string name, bool set = true)
    {
      mesh.EnableTable(name, set);
    }
    
    void HPRefinement (int levels, double factor = 0.125)
    {
      mesh.HPRefinement(levels, factor);
      UpdateBuffers();
    }
    
    void SetDeformation (shared_ptr<GridFunction> def = nullptr);

    const shared_ptr<GridFunction> & GetDeformation () const
    {
      return deformation;
    }

    void SetPML (const shared_ptr<PML_Transformation> & pml_trafo, int _domnr);
    void UnSetPML (int _domnr);

    Array<shared_ptr<PML_Transformation>> & GetPMLTrafos();
    shared_ptr<PML_Transformation> GetPML(int _domnr);

    shared_ptr<netgen::Mesh> GetNetgenMesh () const
    { return mesh.GetMesh(); }

    netgen::Ngx_Mesh & GetNetgenMeshX () { return mesh; }
      
    
    /**
       returns node topology.
       A facet or edge-node knows its vertices etc.
       The method is not yet fully functional.
     */
    template <int DIM>
    netgen::Ng_Node<DIM> GetNode (size_t nr) const
    {
      return mesh.GetNode<DIM> (nr);
    }

    /// returns the points of an element.
    /// vertices and possibly edge-midpoints
    [[deprecated("Use GetElPNums(ElementId) instead!")]]                        
    void GetElPNums (int elnr, Array<int> & pnums) const
    { pnums = ArrayObject (GetElement(ElementId(VOL,elnr)).points); }

    /// returns the points of a boundary element.
    /// vertex and possibly edge-midpoints
    [[deprecated("Use GetElPNums(ElementId) instead!")]]                        
    void GetSElPNums (int selnr, Array<int> & pnums) const
    { pnums = ArrayObject (GetElement(ElementId(BND,selnr)).points); }

    [[deprecated("Use GetElPNums(ElementId) instead!")]]                    
    void GetElPNums (ElementId ei, Array<int> & pnums) const
    { pnums = ArrayObject (GetElement(ei).points); }

    auto GetElPNums (ElementId ei) const
    { return ArrayObject(GetElement(ei).points); }
    
    /// returns the vertices of an element
    [[deprecated("Use GetElVertices(ElementId) instead!")]]                
    void GetElVertices (int elnr, Array<int> & vnums) const
    { vnums = GetElement(ElementId(VOL,elnr)).Vertices(); }
    ///
    [[deprecated("Use GetElVertices(ElementId) -> vnums instead!")]]                
    void GetElVertices (ElementId ei, Array<int> & vnums) const
    { vnums = GetElement(ei).Vertices(); }

    auto GetElVertices (ElementId ei) const
    { return GetElement(ei).Vertices(); }

    /// returns the vertices of a boundary element
    [[deprecated("Use vnums = GetElVertices(ElementId) instead!")]]
    void GetSElVertices (int selnr, Array<int> & vnums) const
    { vnums = GetElement(ElementId(BND,selnr)).Vertices(); }

    [[deprecated("Use enums = GetElEdges(ElementId) instead! ")]]    
    void GetElEdges (ElementId ei, Array<int> & ednums) const
    { ednums = GetElement(ei).Edges(); }

    auto GetElEdges (ElementId ei) const { return GetElement(ei).Edges(); }

    /// returns the edges of an element
    [[deprecated("Use GetElEdges(ElementId) instead!")]]            
    void GetElEdges (int elnr, Array<int> & ednums) const
    { ednums = GetElement(ElementId(VOL,elnr)).Edges(); }

    // returns edge numbers and edge orientation of an element. (old style function)
    // [[deprecated("Use GetElEdges(ElementId) instead!")]]                
    void GetElEdges (int elnr, Array<int> & ednums, Array<int> & orient) const;

    /// returns the edges of a boundary element
    [[deprecated("Use GetElEdges(ElementId) instead!")]]                    
    void GetSElEdges (int selnr, Array<int> & ednums) const
    { ednums = ArrayObject (GetElement(ElementId(BND,selnr)).edges); }

    // returns edge numbers and edge orientation of an element. (old style function)
    // [[deprecated("Use GetElEdges(ElementId) instead, orient is deprecated!")]]
    void GetSElEdges (int selnr, Array<int> & ednums, Array<int> & orient) const;

    /// returns the faces of an element
    // [[deprecated("Use fanums = GetElFaces(ElementId) instead!")]]        
    void GetElFaces (ElementId ei, Array<int> & fnums) const
    { fnums = GetElement(ei).Faces(); }

    auto GetElFaces (ElementId ei) const
    { return GetElement(ei).Faces(); }

    /// returns the faces of an element
    [[deprecated("Use GetElFaces(ElementId) instead!")]]                    
    void GetElFaces (int elnr, Array<int> & fnums) const
    { fnums = GetElement(ElementId(VOL,elnr)).Faces(); }

    // returns face numbers and face orientation of an element. (old style function)
    // [[deprecated("Use GetElFaces(ElementId) instead!")]]                        
    void GetElFaces (int elnr, Array<int> & fnums, Array<int> & orient) const;

    /// returns face number of surface element
    int GetSElFace (int selnr) const;

    // returns face number and orientation of surface element
    // [[deprecated("orientation is deprecated! use GetSElFace(nr) instead")]]
    void GetSElFace (int selnr, int & fnum, int & orient) const;

    /// returns vertex numbers of face
    [[deprecated("Use GetFacePNums(fnr) instead!")]]                        
    void GetFacePNums (int fnr, Array<int> & pnums) const;
    
    /// returns vertex numbers of face
    auto GetFacePNums (size_t fnr) const
    {
      return ArrayObject (mesh.GetNode<2> (fnr).vertices);
    }
    /// returns vertex numbers of edge
    [[deprecated("Use GetEdgePNums(enr) instead!")]]                            
    void GetEdgePNums (int enr, int & pn1, int & pn2) const
    {
      auto edge = mesh.GetNode<1>(enr);
      pn1 = edge.vertices[0];
      pn2 = edge.vertices[1];
    }
    /// returns vertex numbers of edge
    [[deprecated("Use GetEdgePNums(enr) instead!")]]                                
    void GetEdgePNums (int enr, Array<int> & pnums) const;
    /// returns vertex numbers of edge
    /*
    auto GetEdgePNums (size_t enr) const -> decltype(ArrayObject(IVec<2>()))
    {
      int v2[2];
      Ng_GetEdge_Vertices (enr+1, v2);
      return ArrayObject (IVec<2> (v2[0]-1, v2[1]-1));
    }
    */
    auto GetEdgePNums (size_t enr) const
    {
      auto vts = mesh.GetNode<1>(enr).vertices;
      return IVec<2>(vts[0],vts[1]);
    }
    /// returns all elements connected to an edge
    void GetEdgeElements (int enr, Array<int> & elnums) const;
    /// returns all elements connected to an edge
    void GetEdgeSurfaceElements (int enr, Array<int> & elnums) const;
    /// returns all edges of a face
    [[deprecated("Use GetFaceEdges(fnr) -> edges instead!")]]                
    void GetFaceEdges (int fnr, Array<int> & edges) const;
    INLINE auto GetFaceEdges (size_t fnr) const
    { return ArrayObject(mesh.GetFaceEdges(fnr)); }
    /*
    {
      ArrayMem<int,4> f2ed;
      GetFaceEdges (fnr, f2ed);
      f2ed.NothingToDelete(); // dynamic allocation never needed
      return f2ed;
    }
    */
    
    void GetEdgeFaces (int enr, Array<int> & faces) const;
    /// returns elements connected to a face
    void GetFaceElements (int fnr, Array<int> & elnums) const;
    /// returns surface elements connected to a face
    void GetFaceSurfaceElements (int fnr, Array<int> & elnums) const;
    /// point numbers of a 1D element
    // void GetSegmentPNums (int snr, Array<int> & pnums) const;
    /// index of 1D element
    // int GetSegmentIndex (int snr) const;

    [[deprecated("Use sels = GetVertexElements(vnr) instead!")]]                    
    void GetVertexElements (size_t vnr, Array<int> & elems) const;
    auto GetVertexElements (size_t vnr) const 
    { return ArrayObject(mesh.GetNode<0> (vnr).elements); }

    [[deprecated("Use sels = GetVertexSurfaceElements(vnr) instead!")]]                
    void GetVertexSurfaceElements (size_t vnr, Array<int> & elems) const;
    auto GetVertexSurfaceElements (size_t vnr) const 
    { return ArrayObject(mesh.GetNode<0> (vnr).bnd_elements); }

    auto GetVertexElements (size_t vnr, VorB vb) const 
    {
      switch (vb)
        {
        case VOL: return ArrayObject(mesh.GetNode<0> (vnr).elements);
        case BND: return ArrayObject(mesh.GetNode<0> (vnr).bnd_elements);
        default: throw Exception ("GetVertexElements, unhandled vb");
        }
    }

    
    /// number of facets of an element. 
    /// facets are edges (2D) or faces (3D)
    size_t GetNFacets() const { return nnodes_cd[1]; }
    
    /// facets of an element
    auto GetElFacets (ElementId ei) const { return GetElement(ei).Facets(); }

    [[deprecated("Use fanums = GetElFacets(ElementId) instead!")]]            
    void GetElFacets (ElementId ei, Array<int> & fnums) const
      { fnums = GetElFacets(ei); }
    [[deprecated("Use GetElFacets(ElementId) instead!")]]        
    void GetElFacets (int elnr, Array<int> & fnums) const
      { fnums = GetElFacets(ElementId(VOL, elnr)); }      
    /// facet of a surface element
    [[deprecated("Use GetElFacets(ElementId) instead!")]]            
    void GetSElFacets (int selnr, Array<int> & fnums) const
      { fnums = GetElFacets(ElementId(BND, selnr)); }

    /// vertices of a facet
    void GetFacetPNums (int fnr, Array<int> & pnums) const;
    /// geometry type of facet
    ELEMENT_TYPE GetFaceType (int fnr) const
    { return (mesh.GetNode<2>(fnr).vertices.Size() == 3) ? ET_TRIG : ET_QUAD; }
    ELEMENT_TYPE GetFacetType (int fnr) const;    
    /// elements connected to facet
    void GetFacetElements (int fnr, Array<int> & elnums) const
    {
      switch (dim)
        {
        case 1: elnums = GetVertexElements (fnr); break;
        case 2: GetEdgeElements (fnr, elnums); break;
        case 3: GetFaceElements (fnr, elnums); break;
        }
    }
    void GetFacetSurfaceElements (int fnr, Array<int> & elnums) const
    {
      switch (dim)
        {
        case 1: elnums = GetVertexSurfaceElements (fnr); break;
        case 2: GetEdgeSurfaceElements (fnr, elnums); break;
        case 3: GetFaceSurfaceElements (fnr, elnums); break;
        }
    }

    void CalcIdentifiedFacets();
    int GetPeriodicFacet(int fnr) const
    {
      if(get<1>(identified_facets[fnr]) == 2) // check if identification type is periodic
	return get<0>(identified_facets[fnr]);
      else
	return fnr;
    }


    // void GetVertexElements (int vnr, Array<int> & elnrs) const;
    /// element order stored in Netgen
    int GetElOrder (int enr) const
    { return mesh.GetElementOrder (enr+1); } 
    /// anisotropic order stored in Netgen
    IVec<3> GetElOrders (int enr) const
    { 
      IVec<3> eo; 
      mesh.GetElementOrders(enr+1,&eo[0],&eo[1],&eo[2]); 
      return eo; 
    } 
    /// set element order in Netgen
    void SetElOrder (int enr, int order) 
    { mesh.SetElementOrder (enr+1,order); }
    /// set anisotropic element order in Netgen
    void SetElOrders (int enr, int ox, int oy, int oz) 
    { mesh.SetElementOrders (enr+1, ox,oy,oz); }
    /// order of suface element
    int GetSElOrder (int enr) const
    { return mesh.GetSurfaceElementOrder (enr+1); } 
    /// anisotropic order of suface element
    IVec<2> GetSElOrders (int enr) const
    { 
      IVec<2> eo; 
      mesh.GetSurfaceElementOrders(enr+1,&eo[0],&eo[1]); 
      return eo; 
    }
    /// set surface element order
    void SetSElOrder (int enr, int order) 
    { mesh.SetSurfaceElementOrder (enr+1,order); }
    /// set anisotropic surface element order
    void SetSElOrders (int enr, int ox, int oy) 
    { mesh.SetSurfaceElementOrders (enr+1, ox,oy); }
    

    /// the volume of an element (mid-point rule)
    double ElementVolume (int elnr) const;
    /// the area of a boundary element (mid-point rule)
    double SurfaceElementVolume (int selnr) const;
      

    /// number of multigrid levels
    int GetNLevels() const { return nlevels; }  
    auto GetNVLevel(int level) const { return mesh.GetNVLevel(level); }
    int GetHPElementLevel(ElementId ei) const { return mesh.GetHPElementLevel(ei.Nr(),1); }

    /// the two parent vertices of a vertex. -1 for coarse-grid vertices
    void GetParentNodes (int pi, int * parents) const
    { 
      mesh.GetParentNodes (pi, parents);
    }
    IVec<2> GetParentNodes (int pi) const
    {
      IVec<2,int> parents;
      mesh.GetParentNodes (pi, &parents[0]);
      return parents;
    }
    /// number of parent element on next coarser mesh
    [[deprecated("Use GetParentElement(ElementId) instead!")]]                
    int GetParentElement (int ei) const
    { return mesh.GetParentElement (ei); }
    /// number of parent boundary element on next coarser mesh
    [[deprecated("Use GetParentElement(ElementId) instead!")]]                    
    int GetParentSElement (int ei) const
    { return mesh.GetParentSElement (ei); }

    ElementId GetParentElement (ElementId ei) const
    {
      if (ei.VB() == VOL)
        return ElementId(VOL, mesh.GetParentElement(ei.Nr()));
      else if (ei.VB() == BND)
        return ElementId(BND, mesh.GetParentSElement(ei.Nr()));
      else
        throw Exception ("GetParentElement only supported for VOL and BND");
    }

    bool HasParentEdges () const { return mesh.HasParentEdges(); }
    auto GetParentEdges (int enr) const { return mesh.GetParentEdges(enr); }
    auto GetParentFaces (int fnr) const { return mesh.GetParentFaces(fnr); }

    Array<uint64_t> BuildRefinementTree() const;
    void RefineFromTree(const Array<uint64_t> & tree);

    /// representant of vertex for anisotropic meshes
    int GetClusterRepVertex (int pi) const
    { return mesh.GetClusterRepVertex (pi+1)-1; }
    /// representant of edge for anisotropic meshes
    int GetClusterRepEdge (int pi) const
    { return mesh.GetClusterRepEdge (pi+1)-1; }
    /// representant of face for anisotropic meshes
    int GetClusterRepFace (int pi) const
    { return mesh.GetClusterRepFace (pi+1)-1; }
    /// representant of element for anisotropic meshes
    int GetClusterRepElement (int pi) const
    { return mesh.GetClusterRepElement (pi+1)-1; }
    
    
    /// returns the transformation from the reference element to physical element.
    /// Given a point in the reference element, the ElementTransformation can 
    /// compute the mapped point as well as the Jacobian
    [[deprecated("Use GetTrafo with ElementId(vb, elnr) instead!")]]    
      ngfem::ElementTransformation & GetTrafo (int elnr, VorB vb, Allocator & lh) const
      {
        return GetTrafo(ElementId(vb, elnr),lh);
      }
    
    ngfem::ElementTransformation & GetTrafo (ElementId ei, Allocator & lh) const
    {
      auto ptr = trafo_jumptable[ei.VB()];
      if (ptr)
        return (this->*ptr)(ei.Nr(),lh);
      else
        return GetTrafoOld(ei, lh);
    }
    
    ngfem::ElementTransformation & GetTrafoOld (ElementId ei, Allocator & lh) const;
    
    template <int DIM>
      ngfem::ElementTransformation & GetTrafoDim (size_t elnr, Allocator & lh) const;
    template <int DIM>
      ngfem::ElementTransformation & GetSTrafoDim (size_t elnr, Allocator & lh) const;
    template <int DIM>
      ngfem::ElementTransformation & GetCD2TrafoDim (size_t elnr, Allocator & lh) const;
    template <int DIM>
      ngfem::ElementTransformation & GetCD3TrafoDim (size_t elnr, Allocator & lh) const;
    
    template <VorB VB,int DIM>
      ngfem::ElementTransformation & GetTrafo (T_ElementId<VB,DIM> ei, Allocator & lh) const
    {
      switch(VB)
	{
	case VOL:
	  return GetTrafoDim<DIM> (ei.Nr(), lh);
	case BND:
	  return GetSTrafoDim<DIM> (ei.Nr(), lh);
	case BBND:
	  return GetCD2TrafoDim<DIM> (ei.Nr(),lh);
	default:
	  __assume(false);
	}
    }
    

    // (old style optimization)
    [[deprecated("functionality not useful anymore, just remove function call!")]]            
    void SetPointSearchStartElement(const int el) const;


    
    ElementId FindElementOfPoint (FlatVector<double> point,
                                  IntegrationPoint & ip, 
                                  bool build_searchtree,
                                  const Array<int> * const indices = NULL) const;
    ElementId FindElementOfPoint (FlatVector<double> point,
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
    [[deprecated("Use GetElement(id).is_curved instead!")]]        
    bool IsElementCurved (int elnr) const
    { return GetElement(ElementId(VOL,elnr)).is_curved; }
      // { return bool (Ng_IsElementCurved (elnr+1)); }

    [[deprecated("Use GetPeriodicNodes(NT_VERTEX, pairs) instead!")]]
    void GetPeriodicVertices ( Array<IVec<2> > & pairs) const;
    [[deprecated("Use GetNPeriodicNodes(NT_VERTEX) instead!")]]
    int GetNPairsPeriodicVertices () const;
    [[deprecated("Use GetPeriodicNodes(NT_VERTEX, idnr) instead")]]
    void GetPeriodicVertices (int idnr, Array<IVec<2> > & pairs) const;
    [[deprecated("Use GetPeriodicNodes(NT_VERTEX, idnr).Size() instead")]]
    int GetNPairsPeriodicVertices (int idnr) const;  

    [[deprecated("Use GetPeriodicNodes(NT_EDGE, pairs) instead!")]]
    void GetPeriodicEdges ( Array<IVec<2> > & pairs) const;
    [[deprecated("Use GetNPeriodicNodes(NT_EDGE) instead!")]]
    int GetNPairsPeriodicEdges () const;
    [[deprecated("Use GetPeriodicNodes(NT_EDGE, idnr) instead")]]
    void GetPeriodicEdges (int idnr, Array<IVec<2> > & pairs) const;
    [[deprecated("Use GetPeriodicNodes(NT_EDGE, idnr).Size() instead")]]
    int GetNPairsPeriodicEdges (int idnr) const;

    int GetNPeriodicIdentifications() const
    {
      return periodic_node_pairs[NT_VERTEX]->Size();
    }
    // get number of all periodic nodes of nodetype nt
    size_t GetNPeriodicNodes(NODE_TYPE nt) const;
    // write all the node pairs of type nt into array pairs
    void GetPeriodicNodes(NODE_TYPE nt, Array<IVec<2>>& pairs) const;
    // Uses 0 based identification numbers! Returns periodic node pairs of given identifcation number
    const Array<IVec<2>>& GetPeriodicNodes(NODE_TYPE nt, int idnr) const;

    shared_ptr<CoefficientFunction> RegionCF(VorB vb, shared_ptr<CoefficientFunction> default_value,
                                             const Array<pair<variant<string, Region>, shared_ptr<CoefficientFunction>>>& region_values);

    shared_ptr<CoefficientFunction> MaterialCF(shared_ptr<CoefficientFunction> default_value,
                                               const Array<pair<variant<string, Region>, shared_ptr<CoefficientFunction>>>& region_values)
    {
      return RegionCF(VOL, default_value, region_values);
    }

    shared_ptr<CoefficientFunction> BoundaryCF(shared_ptr<CoefficientFunction> default_value,
                                               const Array<pair<variant<string, Region>, shared_ptr<CoefficientFunction>>>& region_values)
    {
      return RegionCF(BND, default_value, region_values);
    }

  private:
    mutable Table<size_t> elements_of_class[4];
    mutable size_t elements_of_class_timestamp[4] = { 0, 0, 0, 0 };
  public:
    const Table<size_t> & GetElementsOfClass(VorB vb = VOL) const; // classification by vertex numbers

  private:
    Array<bool> higher_integration_order;
  public:
    void SetHigherIntegrationOrder(int elnr);
    void UnSetHigherIntegrationOrder(int elnr);

    // void LoadMesh (const string & filename);
    // void LoadMesh (istream & str);
    void SaveMesh (ostream & str) const;
    void DoArchive(Archive& ar);
    // void LoadMeshFromString(const string & str);

    // void PrecomputeGeometryData(int intorder);

    void InitPointCurve(double red = 1, double green = 0, double blue = 0) const;
    void AddPointCurvePoint(const Vec<3> & point) const;


    
    template <int DIMS, int DIMR> friend class Ng_ElementTransformation;
    template <int DIMS, int DIMR> friend class Ng_ConstElementTransformation;


    NgMPI_Comm GetCommunicator () const { return mesh.GetCommunicator(); }

    /**
       Returns the list of other MPI - processes where node is present.
       The ordering coincides for all processes.
    */
    [[deprecated("Use GetDistantProcs (NodeId) instead!")]]                
    void GetDistantProcs (NodeId node, Array<int> & procs) const;

    /**
       Returns the global number of the node.
       Currently, this function works only for vertex-nodes.
     */

    [[deprecated("should not need global numbers")]]                    
    size_t GetGlobalNodeNum (NodeId node) const;

    size_t GetGlobalVertexNum (int locnum) const;
    
    FlatArray<int> GetDistantProcs (NodeId node) const
    {
      return mesh.GetDistantProcs(StdNodeType(node.GetType(), GetDimension()), node.GetNr());
    }

    /// Reduces MPI - distributed data associated with mesh-nodes
    template <typename T>
    void AllReduceNodalData (NODE_TYPE nt, Array<T> & data, NG_MPI_Op op) const;

  private:
    void BuildNeighbours();
  };


  


  INLINE Ngs_Element ElementIterator :: operator*() const { return ma[ei]; }
  
  template <VorB VB>
  INLINE Ngs_Element TElementIterator<VB>::operator*() const 
  {
    return ma[ElementId(VB,nr)]; 
  }

  template <VorB VB, int DIM>
  INLINE Ngs_Element DimElementIterator<VB,DIM>::operator*() const 
  {
    return ma.GetElement(T_ElementId<VB,DIM>(nr));
  }

  
  class Region
  {
    shared_ptr<MeshAccess> mesh;
    VorB vb;
    shared_ptr<BitArray> mask;
  public:
    Region() {}
    NGS_DLL_HEADER Region (const shared_ptr<MeshAccess> & amesh, VorB avb, string pattern);
    NGS_DLL_HEADER Region (const shared_ptr<MeshAccess> & amesh, VorB avb, const char* pattern) : Region(amesh, avb, string(pattern)) {}
    NGS_DLL_HEADER Region (const shared_ptr<MeshAccess> & amesh, VorB avb, const BitArray & amask);
    NGS_DLL_HEADER Region (const shared_ptr<MeshAccess> & amesh, VorB avb, bool all = false);
    Region (const Region &) = default;
    explicit operator VorB () const { return vb; }
    VorB VB() const { return vb; }
    bool IsVolume () const { return vb == VOL; }
    bool IsBoundary () const { return vb == BND; }
    bool IsCoDim2() const { return vb == BBND; }
    const BitArray & Mask() const { return *mask; }
    BitArray& Mask() { return *mask; }
    operator const BitArray & () const { return *mask; }
    shared_ptr<BitArray> MaskPtr() { return mask; }
    const shared_ptr<MeshAccess> & Mesh() const { return mesh; }
    Region operator+ (const Region & r2) const
    {
      return Region (mesh, vb, BitArray(*mask).Or(r2.Mask()));
    }
    Region operator- (const Region & r2) const
    {
      return Region (mesh, vb, BitArray(*mask).And(BitArray(r2.Mask()).Invert()));
    }
    Region operator~ () const
    {
      return Region (mesh, vb, BitArray(*mask).Invert());
    }
    Region operator+ (const string & pattern2) const
    {
      return *this + Region(mesh, vb, pattern2);
    }
    Region operator- (const string & pattern2) const
    {
      return *this - Region(mesh, vb, pattern2);
    }

    Region operator* (const Region& r2) const
    {
      return Region(mesh, vb, BitArray(*mask).And(r2.Mask()));
    }

    Region operator* (const string& pattern) const
    {
      return *this * Region(mesh, vb, pattern);
    }

    bool operator==(const Region& other) const;
    size_t Hash() const;

    // Get adjacent regions of codim other_vb
    Region GetNeighbours(VorB other_vb);
    Region GetBoundaries() { return GetNeighbours(VorB(int(vb)+1)); }

    auto GetElements() const
    {
      return mesh->Elements(vb)
        | filter([&](auto ei) { return mask->Test(mesh->GetElIndex(ei)); });
    }
  };


  shared_ptr<CoefficientFunction>
  MakeBoundaryFromVolumeCoefficientFunction  (shared_ptr<CoefficientFunction> avol_cf);

  shared_ptr<CoefficientFunction>
  MakeTrafoCF(shared_ptr<CoefficientFunction> func,
              shared_ptr<CoefficientFunction> trafo,
              Region region);
  
  /**
     Controls the progress - output.
     It controls the Netgen - progressbar as well as console progress update.
     In parallel, all processes must enter and call the Done method.
   */
  class ProgressOutput
  {
    shared_ptr<MeshAccess> ma;
    NgMPI_Comm comm;
    string task;
    size_t total;
    double prevtime;
    bool is_root;
    bool use_mpi = true;
    bool done_called;

    static atomic<size_t> cnt;
    static thread_local size_t thd_cnt;
    // static thread_local double thd_prev_time;
    static thread_local size_t thd_prev_time;
  public:
    NGS_DLL_HEADER ProgressOutput (shared_ptr<MeshAccess> ama,
                                   string atask, size_t atotal);
    NGS_DLL_HEADER ~ProgressOutput();

    // update thd-local counter, and after some time also atomic node-local cnt
    NGS_DLL_HEADER void Update();
    // transfer thd-local counter to atomic node-local cnt    
    NGS_DLL_HEADER static void SumUpLocal();    
    NGS_DLL_HEADER void Update(size_t nr);
    NGS_DLL_HEADER void Done();
    // NGS_DLL_HEADER void DoneThread();
  };



  /// old style, compatibility for a while
  template <typename T>  [[deprecated("use ma.AllReduceNodalData instead")]]
  void AllReduceNodalData (NODE_TYPE nt, Array<T> & data, NG_MPI_Op op, const MeshAccess & ma)
  {
    ma.AllReduceNodalData (nt, data, op);
  }


#ifdef PARALLEL

  template <typename T>
  void MeshAccess::
  AllReduceNodalData (NODE_TYPE nt, Array<T> & data, NG_MPI_Op op) const
  {
    NgMPI_Comm comm = GetCommunicator();
    if (comm.Size() <= 1) return;

    Array<int> cnt(comm.Size());
    cnt = 0;
    for (auto i : Range(GetNNodes(nt)))
      for (auto p : GetDistantProcs(Node(nt,i)))
        cnt[p]++;
    
    Table<T> dist_data(cnt), recv_data(cnt);
    
    cnt = 0;
    for (auto i : Range(GetNNodes(nt)))
      for (auto p : GetDistantProcs(Node(nt, i)))
        dist_data[p][cnt[p]++] = data[i];

    Array<NG_MPI_Request> requests;
    for (auto i : cnt.Range())
      if (cnt[i])
	{
	  requests.Append (comm.ISend(dist_data[i], i, NG_MPI_TAG_SOLVE));
	  requests.Append (comm.IRecv(recv_data[i], i, NG_MPI_TAG_SOLVE));
	}
    MyMPI_WaitAll (requests);
    
    cnt = 0;
    NG_MPI_Datatype type = GetMPIType<T>();
    for (auto i : Range(GetNNodes(nt)))
      for (auto p : GetDistantProcs(Node(nt, i)))
        NG_MPI_Reduce_local (&recv_data[p][cnt[p]++],
                          &data[i], 1, type, op);
  }


#else

  template <typename T>
  void MeshAccess::
  AllReduceNodalData (NODE_TYPE nt, Array<T> & data, NG_MPI_Op op) const { ; }

#endif

}

#endif
