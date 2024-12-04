#ifndef FILE_TOPOLOGY
#define FILE_TOPOLOGY

/*********************************************************************/
/* File:   topology.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/



#include <bla.hpp>

namespace ngfem
{
  using namespace ngbla;
  
  /*
    Toplogy of reference elements
  */


  /**
     Geometry of element.
     Possible are ET_POINT, ET_SEGM, ET_TRIG, ET_QUAD, ET_TET, ET_PYRAMID, ET_PRISM, ET_HEX
  */
  enum NGS_DLL_HEADER ELEMENT_TYPE 
    { ET_POINT = 0, ET_SEGM = 1,
	ET_TRIG = 10, ET_QUAD = 11, 
	ET_TET = 20, ET_PYRAMID = 21, ET_PRISM = 22, ET_HEXAMID = 23, ET_HEX = 24 };
  
  // #ifndef WIN32
  static constexpr initializer_list<ELEMENT_TYPE> element_types =
    { ET_POINT, ET_SEGM, ET_TRIG, ET_QUAD,
      ET_TET, ET_PYRAMID, ET_PRISM, ET_HEXAMID, ET_HEX };
  /*
#else
  // less efficient, but MSVC doesn't like the constexpr here
  static initializer_list<ELEMENT_TYPE> element_types =
    { ET_POINT, ET_SEGM, ET_TRIG, ET_QUAD,
      ET_TET, ET_PYRAMID, ET_PRISM, ET_HEX };
#endif
  */
  

  /**
     Type of node. 
     vertex nodes are 0 dimensional, 
     edge nodes are 1 dimensional
     face nodes are 2 dimensional
     cell nodes are 3 dimensional
     global nodes are shared by all procs and have no attached dimension

     2D elements have vertex, edge, and face nodes

     3D elements have vertex, edge, face and cell nodes
  */

  enum NODE_TYPE { NT_VERTEX = 0, NT_EDGE = 1, NT_FACE = 2, NT_CELL = 3, NT_ELEMENT = 4, NT_FACET = 5, NT_GLOBAL =  6 };

  INLINE NODE_TYPE StdNodeType (NODE_TYPE nt, int meshdim)
  {
    switch (nt)
      {
      case NT_ELEMENT: return NODE_TYPE(meshdim);
      case NT_FACET: return NODE_TYPE(meshdim-1);
      default: return nt;
      }
  }

  INLINE int CoDimension (NODE_TYPE nt, int meshdim)
  {
    int dim(nt);
    if (dim <= 3)  // V, E, F C
      return meshdim-dim;
    else
      return dim-NT_ELEMENT; 
  }
  // INLINE void operator++(NODE_TYPE & nt, int)  { nt = NODE_TYPE(nt+1); } 

  INLINE constexpr int Dim (ELEMENT_TYPE et)
  {
    return (et == ET_POINT) ? 0 :
      (et == ET_SEGM) ? 1 : 
      (et == ET_TRIG || et == ET_QUAD) ? 2 : 3;
  }



  /// point coordinates
  typedef double POINT3D[3];  

  /// initial point, end point
  typedef int EDGE[2];      

  /// points, last one is -1 for trig
  typedef int FACE[4];      

  /// normal vector
  typedef double NORMAL[3];



  enum VorB : uint8_t { VOL, BND, BBND, BBBND };
  inline VorB operator++(VorB & vb, int)  { VorB vbo = vb; vb = VorB(vb+1); return vbo; }
  inline VorB & operator++(VorB & vb)  { vb = VorB(vb+1); return vb; }   
  inline ostream & operator<< (ostream & ost, VorB vb)
  {
    if (vb == VOL) ost << "VOL";
    else if (vb==BND) ost << "BND";
    else if (vb==BBND) ost << "BBND";
    else ost << "BBBND";
    return ost;
  }

  class ElementId
  {
    typedef size_t int_type;
    VorB vb;
    int_type nr;
  public:    
    ElementId (VorB avb, int_type anr) : vb(avb), nr(anr) { ; }
    ElementId (int_type anr) : vb(VOL), nr(anr) { ; }
    int_type Nr() const { return nr; }
    explicit operator int_type () const { return nr; }
    explicit operator VorB () const { return vb; }
    VorB VB() const { return vb; }
    bool IsVolume() const { return vb == VOL; }
    bool IsBoundary() const { return vb == BND; }
    bool IsInvalid() const { return nr == -1; }
    bool operator< (int_type nr2) { return nr < nr2; }
    ElementId operator++ (int) { return ElementId(vb,nr++); }
    ElementId operator++ () { return ElementId(vb,++nr); }
    ElementId operator*() const { return *this; }
    bool operator!=(const ElementId id2) const { return nr != id2.nr || vb != id2.vb; }
    bool operator==(const ElementId id2) const { return nr == id2.nr && vb == id2.vb; }
  };

  inline ostream & operator<< (ostream & ost, ElementId id)
  {
    const char * name[4] = { "VEl", "BEl", "CD2El", "CD3El" };
    return ost << name[id.VB()] << ' ' << id.Nr();
    // return ost << (id.VB()==VOL ? "VEl " : (id.VB()==BND ? "BEl " : "CD2El ")) << ' ' << id.Nr();
  }

  template <VorB VB,int DIM>
  class T_ElementId
  {
    size_t nr;
  public:
    T_ElementId (size_t anr) : nr(anr) { ; }
    T_ElementId (ElementId ei) : nr(ei.Nr()) { ; }
    operator ElementId() const { return ElementId(VB, nr); }
    size_t Nr() const { return nr; } 
  };

  
  /// Topology and coordinate information of master element:
  class NGS_DLL_HEADER ElementTopology
  {
    ELEMENT_TYPE myet;
  public:
    ElementTopology (ELEMENT_TYPE amyet) : myet(amyet) { ; }
    /// returns name of element type
    static const char * GetElementName (ELEMENT_TYPE et);
    const char * GetElementName () { return GetElementName(myet); }

    /// returns space dimension of element type
    static INLINE int GetSpaceDim (ELEMENT_TYPE et)
    {
      switch (et)
	{
	case ET_POINT: return 0;
	case ET_SEGM: return 1;
	case ET_TRIG: return 2;
	case ET_QUAD: return 2;
	case ET_TET: return 3;
	case ET_PYRAMID: return 3;
	case ET_PRISM: return 3;
	case ET_HEXAMID: return 3;
	case ET_HEX: return 3;
	}
      return 0;
    }
    INLINE int GetSpaceDim () const { return GetSpaceDim(myet); }

    /// returns number of vertices
    static INLINE int GetNVertices (ELEMENT_TYPE et)
    {
      switch (et)
	{
	case ET_POINT: return 1;
	case ET_SEGM: return 2;
	case ET_TRIG: return 3;
	case ET_QUAD: return 4;
	case ET_TET: return 4;
	case ET_PYRAMID: return 5;
	case ET_PRISM: return 6;
	case ET_HEXAMID: return 7;
	case ET_HEX: return 8;
	}
      return 0;
    }
    INLINE int GetNVertices () const { return GetNVertices(myet); } 
      
    /// returns number of edges
    static INLINE int GetNEdges (ELEMENT_TYPE et)
    { 
      switch (et)
	{
	case ET_POINT: return 0;
	case ET_SEGM: return 1;
	case ET_TRIG: return 3;
	case ET_QUAD: return 4;
	case ET_TET: return 6;
	case ET_PYRAMID: return 8;
	case ET_PRISM: return 9;
	case ET_HEXAMID: return 11;
	case ET_HEX: return 12;
	}
      return 0;
    }


    /// returns number of faces
    static INLINE int GetNFaces (ELEMENT_TYPE et)
    {
      switch (et)
	{
	case ET_POINT: return 0;
	case ET_SEGM: return 0;
	case ET_TRIG: return 1;
	case ET_QUAD: return 1;
	case ET_TET: return 4;
	case ET_PYRAMID: return 5;
	case ET_PRISM: return 5;
	case ET_HEXAMID: return 6;
	case ET_HEX: return 6;
	}  
      return 0;
    }


    /// returns face type
    static INLINE ELEMENT_TYPE GetFaceType (ELEMENT_TYPE et, int k)
    {
      switch (et)
	{
	case ET_SEGM: return ET_POINT;
	case ET_TRIG: return ET_TRIG;
	case ET_QUAD: return ET_QUAD;
	case ET_TET: return ET_TRIG;
	case ET_PYRAMID: return (k<4 ? ET_TRIG : ET_QUAD); 
	case ET_PRISM: return (k<2 ? ET_TRIG : ET_QUAD);
	case ET_HEXAMID: return (k==1||k==4) ? ET_TRIG : ET_QUAD;
	case ET_HEX: return ET_QUAD;
	default:
	  return ET_SEGM;
	}  
    }


    static INLINE Vec<4,int> GetNNodes (ELEMENT_TYPE et)
    {
      switch (et)
        {
        case ET_POINT  : return Vec<4,int> (1,0,0,0);
        case ET_SEGM   : return Vec<4,int> (2,1,0,0);
        case ET_TRIG   : return Vec<4,int> (3,3,1,0);
        case ET_QUAD   : return Vec<4,int> (4,4,1,0);
        case ET_TET    : return Vec<4,int> (4,6,4,1);
        case ET_PYRAMID : return Vec<4,int> (5,8,5,1);
        case ET_PRISM   : return Vec<4,int> (6,9,5,1);
        case ET_HEXAMID : return Vec<4,int> (7,11,6,1);          
        case ET_HEX     : return Vec<4,int> (8,12,6,1);
        }
      return 0;
    }

    /// returns number of nodes of type nt
    static int GetNNodes (ELEMENT_TYPE et, NODE_TYPE nt)
    {
      static const int nn_point[] = { 1, 0, 0, 0 };
      static const int nn_segm[] = { 2, 1, 0, 0 };
      static const int nn_trig[] = { 3, 3, 1, 0 };
      static const int nn_quad[] = { 4, 4, 1, 0 };
      static const int nn_tet[] = { 4, 6, 4, 1 };
      static const int nn_pyramid[] = { 5, 8, 5, 1 };
      static const int nn_prism[] = { 6, 9, 5, 1 };
      static const int nn_hexamid[] = { 7, 11, 6, 1 };
      static const int nn_hex[] = { 8, 12, 6, 1 };
      switch (et)
	{
	case ET_POINT: return nn_point[nt];
	case ET_SEGM: return nn_segm[nt];
	case ET_TRIG: return nn_trig[nt];
	case ET_QUAD: return nn_quad[nt];
	case ET_TET: return nn_tet[nt];
	case ET_PYRAMID: return nn_pyramid[nt];
	case ET_PRISM: return nn_prism[nt];
	case ET_HEXAMID: return nn_hexamid[nt];
	case ET_HEX: return nn_hex[nt];
	}  
      return 0;
    }







    /// returns number of facets: == GetNFaces in 3D, GetNEdges in 2D
    static INLINE int GetNFacets (ELEMENT_TYPE et)
    {
      switch (et)
	{
	case ET_POINT: return 0;
	case ET_SEGM: return 2;
	case ET_TRIG: return 3;
	case ET_QUAD: return 4;
	case ET_TET: return 4;
	case ET_PYRAMID: return 5;
	case ET_PRISM: return 5;
	case ET_HEXAMID: return 6;
	case ET_HEX: return 6;          
	}  
      return 0;
    }




    /// returns number of facets: == GetNFaces in 3D, GetNEdges in 2D
    static INLINE ELEMENT_TYPE GetFacetType (ELEMENT_TYPE et, int k)
    {
      switch (et)
	{
	case ET_SEGM: return ET_POINT;
	case ET_TRIG: return ET_SEGM;
	case ET_QUAD: return ET_SEGM;
	case ET_TET: return ET_TRIG;
	case ET_PYRAMID: return (k<4 ? ET_TRIG : ET_QUAD); 
	case ET_PRISM: return (k<2 ? ET_TRIG : ET_QUAD);
	case ET_HEXAMID: return (k==1||k==4) ? ET_TRIG : ET_QUAD;          
	case ET_HEX: return ET_QUAD;
	default:
	  return ET_POINT; // dummy
	}  
    }

    /// returns vertex coordinates (as 3D points)
    static  const POINT3D * GetVertices (ELEMENT_TYPE et);
    INLINE const POINT3D * GetVertices () const
    { return GetVertices(myet); }

    /// returns edges of elements. zero-based pairs of integers
    static const EDGE * GetEdges (ELEMENT_TYPE et)
    {
      static const int segm_edges[1][2] =
	{ { 0, 1 }};
    

      static const int trig_edges[3][2] =
	{ { 2, 0 },
	  { 1, 2 },
	  { 0, 1 }};

      static const int quad_edges[4][2] =
	{ { 0, 1 },
	  { 2, 3 },
	  { 3, 0 },
	  { 1, 2 }};
    
      static const int tet_edges[6][2] =
	{ { 3, 0 },
	  { 3, 1 },
	  { 3, 2 }, 
	  { 0, 1 }, 
	  { 0, 2 },
	  { 1, 2 }};
    
      static const int prism_edges[9][2] =
	{ { 2, 0 },
	  { 0, 1 },
	  { 2, 1 },
	  { 5, 3 },
	  { 3, 4 },
	  { 5, 4 },
	  { 2, 5 },
	  { 0, 3 },
	  { 1, 4 }};

      static const int pyramid_edges[8][2] =
	{ { 0, 1 },
	  { 1, 2 },
	  { 0, 3 },
	  { 3, 2 },
	  { 0, 4 },
	  { 1, 4 },
	  { 2, 4 },
	  { 3, 4 }};

      static const int hexamid_edges[11][2] =
        {
          { 0, 1 },
          { 2, 3 },
          { 3, 0 },
          { 1, 2 },
          { 4, 5 },
          { 6, 4 },
          { 5, 6 },
          { 0, 4 },
          { 1, 5 },
          { 2, 6 },
          { 3, 6 },
        };

      
      static const int hex_edges[12][2] =
	{
	  { 0, 1 },
	  { 2, 3 },
	  { 3, 0 },
	  { 1, 2 },
	  { 4, 5 },
	  { 6, 7 },
	  { 7, 4 },
	  { 5, 6 },
	  { 0, 4 },
	  { 1, 5 },
	  { 2, 6 },
	  { 3, 7 },
	};
    
      switch (et)
	{
        case ET_POINT: return nullptr;
	case ET_SEGM: return segm_edges;
	case ET_TRIG: return trig_edges;
	case ET_QUAD: return quad_edges;
	case ET_TET:  return tet_edges;
	case ET_PYRAMID: return pyramid_edges;
	case ET_PRISM: return prism_edges;
	case ET_HEXAMID: return hexamid_edges;
	case ET_HEX: return hex_edges;
	default:
	  break;
	}
      cerr << "Ng_GetEdges, illegal element type " << et << endl;
      return 0;  
    }

    /// returns faces of elements. zero-based array of 4 integers, last one is -1 for triangles
    static const FACE * GetFaces (ELEMENT_TYPE et)
    {
      static int tet_faces[4][4] =
	{ { 3, 1, 2, -1 },
	  { 3, 2, 0, -1 },
	  { 3, 0, 1, -1 },
	  { 0, 2, 1, -1 } }; // all faces point into interior!
  
      static int prism_faces[5][4] =
	{
	  { 0, 2, 1, -1 },
	  { 3, 4, 5, -1 },
	  { 2, 0, 3, 5 },
	  { 0, 1, 4, 3 },
	  { 1, 2, 5, 4 } 
	};

      static int pyramid_faces[5][4] =
	{
	  { 0, 1, 4, -1 },
	  { 1, 2, 4, -1 },
	  { 2, 3, 4, -1 },
	  { 3, 0, 4, -1 },
	  { 0, 1, 2, 3 } // points into interior!
	};
  
      static int hexamid_faces[6][4] =
	{
          { 0, 3, 2, 1 },
          { 4, 5, 6, -1},
          { 0, 1, 5, 4 },
          { 1, 2, 6, 5 },
          { 2, 3, 6, -1},
          { 3, 0, 4, 6 }
	};
    
      static int hex_faces[6][4] =
	{
	  { 0, 3, 2, 1 },
	  { 4, 5, 6, 7 },
	  { 0, 1, 5, 4 },
	  { 1, 2, 6, 5 },
	  { 2, 3, 7, 6 },
	  { 3, 0, 4, 7 }
	};
      static int trig_faces[1][4] = 
	{
	  { 0, 1, 2, -1 },
	};
      
      static int quad_faces[1][4] = 
	{
	  { 0, 1, 2, 3 },
	};
    
      switch (et)
	{
	case ET_TET: return tet_faces;
	case ET_PRISM: return prism_faces;
	case ET_PYRAMID: return pyramid_faces;
	case ET_HEXAMID: return hexamid_faces;
	case ET_HEX: return hex_faces;          

	case ET_TRIG: return trig_faces;
	case ET_QUAD: return quad_faces;
        
	case ET_SEGM: return nullptr;
        case ET_POINT: return nullptr;          
	default:
	  break;
	}
    
      cerr << "Ng_GetFaces, illegal element type " << et << endl;
      return 0;
    }

    /// return normals on facets (old style)
    static NORMAL * GetNormals(ELEMENT_TYPE et);

    template <int D>
    static FlatVector<Vec<D> > GetNormals(ELEMENT_TYPE et);

  
    /// returns number of edge from vertex v1 to vertex v2
    static int GetEdgeNr (ELEMENT_TYPE et, int v1, int v2);

    /// returns number of face containing vertices v1, v2, v3.  (trig only ?)
    static int GetFaceNr (ELEMENT_TYPE et, int v1, int v2, int v3);
  };



  /**
     A Node of an element.  

     A Node has a node type such such NT_VERTEX or NT_FACE, and a node
     number. The number can be with respect to the local element
     numbering, or can be the global numbering on the mesh.
  */
  class /* NGS_DLL_HEADER */ NodeId
  {
    NODE_TYPE nt;
    size_t nodenr;

  public:
    /// do nothing
    NodeId () { ; }
  
    /// construct node from type and number
    NodeId (NODE_TYPE ant, size_t anodenr)
      : nt(ant), nodenr(anodenr) { ; }

    /// copy constructor
    NodeId (const NodeId & n2)
    { nt = n2.nt; nodenr = n2.nodenr; }

    /// returns type of the node
    NODE_TYPE GetType () const { return nt; }

    /// returns number of the node
    size_t GetNr() const { return nodenr; }

    operator size_t () const { return nodenr; }
    NodeId operator++ (int) { return NodeId(nt,nodenr++); }
    NodeId operator++ () { return NodeId(nt,++nodenr); }
    NodeId operator+ (size_t i) { return NodeId(nt,nodenr+i); }    
    // NodeId operator*() const { return *this; }
    bool operator!=(const NodeId id2) const { return nodenr != id2.nodenr || nt != id2.nt; }
    bool operator==(const NodeId id2) const { return nodenr == id2.nodenr && nt == id2.nt; }
    size_t operator- (NodeId id2) const { return nodenr-id2.nodenr; }
  };
  typedef NodeId Node;

  template <NODE_TYPE nt>
  class T_NodeId
  {
    size_t nodenr;

  public:
    T_NodeId () = default;
    T_NodeId (const T_NodeId & n2) = default;
  
    /// construct node from number
    T_NodeId (size_t anodenr)
      : nodenr(anodenr) { ; }

    /// returns type of the node
    NODE_TYPE GetType () const { return nt; }

    /// returns number of the node
    size_t GetNr() const { return nodenr; }

    operator size_t () const { return nodenr; }
    T_NodeId operator++ (int) { return T_NodeId(nodenr++); }
    T_NodeId operator++ () { return T_NodeId(++nodenr); }
    bool operator!=(const T_NodeId id2) const { return nodenr != id2.nodenr; }
    bool operator==(const T_NodeId id2) const { return nodenr == id2.nodenr; }
    size_t operator- (T_NodeId id2) const { return nodenr-id2.nodenr; }
    operator NodeId () const { return NodeId(nt, nodenr); }
  };

  
  
  inline int CalcNodeId (ELEMENT_TYPE et, const NodeId & node)
  {
    switch (et)
      {
      case ET_TRIG: 
	{
	  static const int nodebase[] = { 0, 3, 6 };
	  return nodebase[node.GetType()] + node.GetNr();
	}

      case ET_TET: 
	{
	  static const int nodebase[] = { 0, 4, 10, 14 };
	  return nodebase[node.GetType()] + node.GetNr();
	}
      default:
	throw Exception (string ("CalcNodeId not implemented for element ") +
			 ElementTopology::GetElementName (et));
      }
  }

  inline NodeId CalcNodeFromId (ELEMENT_TYPE et, int nodeid)
  {
    switch (et)
      {
      case ET_TRIG: 
	{
	  static const int nodetypes[] = { 0, 0, 0, 1, 1, 1, 2 };
	  static const int nodenrs[]   = { 0, 1, 2, 0, 1, 2, 0 };
	  return NodeId (NODE_TYPE(nodetypes[nodeid]), nodenrs[nodeid]);
	}
      case ET_TET: 
	{
	  static const int nodetypes[] = { 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3 };
	  static const int nodenrs[]   = { 0, 1, 2, 3, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 0 };
	  return NodeId (NODE_TYPE(nodetypes[nodeid]), nodenrs[nodeid]);
	}

      default:
	throw Exception (string ("CalcNodeFromId not implemented for element ") +
			 ElementTopology::GetElementName (et));
      }
  }


  NGS_DLL_HEADER ostream & operator<< (ostream & ost, ELEMENT_TYPE et);
  NGS_DLL_HEADER ostream & operator<< (ostream & ost, NODE_TYPE nt);

  NGS_DLL_HEADER ostream & operator<< (ostream & ost, const NodeId & node);

  /**
     A binary representation of selecting V-E-F-C Nodes.
  */


  class NodeSet
  {
    int set;
  public:
    
    NodeSet (NODE_TYPE nt1) 
    {
      set = 1 << nt1;
    }

    NodeSet (NODE_TYPE nt1, NODE_TYPE nt2) 
    {
      set =  (1 << nt1) + (1 << nt2);
    }

    NodeSet (NODE_TYPE nt1, NODE_TYPE nt2, NODE_TYPE nt3) 
    {
      set = (1 << nt1) + (1 << nt2) + (1 << nt3);
    }

    NodeSet (NODE_TYPE nt1, NODE_TYPE nt2, NODE_TYPE nt3, NODE_TYPE nt4) 
    {
      set = (1 << nt1) + (1 << nt2) + (1 << nt3) + (1 << nt4);
    }
    operator int () const { return set; }
  };




  /*
  class TopologicElement
  {
    ELEMENT_TYPE et;
    ArrayMem<Node, 27> nodes;

  public:
    void SetElementType (ELEMENT_TYPE aet) { et = aet; }
    void Clear() { nodes.SetSize(0); }
    void AddNode (const Node & node) { nodes.Append (node); }


    ELEMENT_TYPE GetType () const { return et; }
    int GetNNodes () const { return nodes.Size(); }
    const Node & GetNode (int i) const { return nodes[i]; }

    Node GlobalNode (const Node & locnode)
    { return Node (locnode.GetType(), nodes[CalcNodeId (et, locnode)].GetNr() ); }
  };

  ostream & operator<< (ostream & ost, const TopologicElement & etop);
  */

  template <int D> class DIM_trait { };

  template <> class DIM_trait<0>
  {
  public:
    enum { MAX_VERTEX = 1 };
    enum { MAX_EDGE = 0 };
    enum { MAX_FACE = 0 };
    enum { MAX_CELL = 0 };
  };

  template <> class DIM_trait<1>
  {
  public:
    enum { MAX_VERTEX = 2 };
    enum { MAX_EDGE = 1 };
    enum { MAX_FACE = 0 };
    enum { MAX_CELL = 0 };
  };

  template <> class DIM_trait<2>
  {
  public:
    enum { MAX_VERTEX = 4 };
    enum { MAX_EDGE = 4 };
    enum { MAX_FACE = 1 };
    enum { MAX_CELL = 0 };
  };

  template <> class DIM_trait<3>
  {
  public:
    enum { MAX_VERTEX = 8 };
    enum { MAX_EDGE = 12 };
    enum { MAX_FACE = 6 };
    enum { MAX_CELL = 1 };
  };




  
  template <int ET> class ET_trait { };

  template<> class ET_trait<ET_POINT>
  {
  public:
    enum { DIM = 0 };
    enum { N_VERTEX = 1 };
    enum { N_EDGE = 0 };
    enum { N_FACE = 0 };
    enum { N_CELL = 0 };
    enum { N_FACET = 0 };

    static constexpr ELEMENT_TYPE ElementType() { return ET_POINT; }
    constexpr operator ELEMENT_TYPE() const { return ET_POINT; }
    
    static INLINE int PolDimension (IVec<1> order) { return 1; }
    static INLINE int PolBubbleDimension (IVec<1> order) { return 0; }

    static INLINE ELEMENT_TYPE FaceType(int i) { return ET_POINT; }  // dummy



    static INLINE IVec<2> GetEdge (int /* i */);
    /*
    {
      static const int edges[][2] = 
	{ { 0, 1 } };
      return IVec<2> (edges[0][0], edges[0][1]);
    }
    */

    template <typename TVN>
    static INLINE IVec<2> GetEdgeSort (int i, const TVN & vnums);
    /*
    {
      IVec<2> e = GetEdge (i);
      if (vnums[e[0]] > vnums[e[1]]) swap (e[0], e[1]);
      return e;
    }
    */

    static IVec<4> GetFace (int /* i */ );
    /*
    {
      return IVec<4> (-1, -1, -1, -1);
    }
    */
    template <typename TVN>
    static IVec<4> GetFaceSort (int /*  i */ , const TVN & vnums);
    /*
    {
      return GetFace(0);
    }
    */

    template <typename TVN>
    static INLINE int GetClassNr (const TVN & vnums)
    {
      return 0;
    }

    template <typename TVN>
    static INLINE int GetFacetClassNr (int facet, const TVN & vnums)
    {
      return 0;
    } 

  };





  template<> class ET_trait<ET_SEGM>
  {
  public:
    enum { DIM = 1 };
    enum { N_VERTEX = 2 };
    enum { N_EDGE = 1 };
    enum { N_FACE = 0 };
    enum { N_CELL = 0 };
    enum { N_FACET = 2 };
    
    static constexpr ELEMENT_TYPE ElementType() { return ET_SEGM; }
    constexpr operator ELEMENT_TYPE() const { return ET_SEGM; }
    
    static INLINE int PolDimension (IVec<1> order) { return order[0]+1; }
    static INLINE int PolBubbleDimension (IVec<1> order) { return (order[0] <= 1) ? 0 : order[0]-1; }

    static INLINE ELEMENT_TYPE FaceType(int i) { return ET_TRIG; }

    static INLINE IVec<2> GetEdge (int /* i */)
    {
      const int edges[][2] = 
	{ { 0, 1 } };
      return IVec<2> (edges[0][0], edges[0][1]);
    }

    template <typename TVN>
    static INLINE IVec<2> GetEdgeSort (int i, const TVN & vnums)
    {
      IVec<2> e = GetEdge (i);
      if (vnums[e[0]] > vnums[e[1]]) Swap (e[0], e[1]);
      return e;
    }


    static INLINE IVec<4> GetFace (int /* i */ )
    {
      return IVec<4> (-1, -1, -1, -1);
    }

    template <typename TVN>
    static INLINE IVec<4> GetFaceSort (int /*  i */ , const TVN & vnums)
    {
      return GetFace(0);
    }

    template <typename TVN>
    static INLINE int GetClassNr (const TVN & vnums)
    {
      int classnr = 0;
      int sort[3] = { 0, 1 };
      if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 1; }
      return classnr;
    }

    template <typename TVN>
    static INLINE int GetFacetClassNr (int facet, const TVN & vnums)
    {
      return facet;
    } 
  };

  template<> class ET_trait<ET_TRIG>
  {
  public:
    enum { DIM = 2 };
    enum { N_VERTEX = 3 };
    enum { N_EDGE = 3 };
    enum { N_FACE = 1 };
    enum { N_CELL = 0 };
    enum { N_FACET = 3 };
    
    static constexpr ELEMENT_TYPE ElementType() { return ET_TRIG; }
    constexpr operator ELEMENT_TYPE() const { return ET_TRIG; }
    
    static constexpr int PolDimension (int order) { return (order+1)*(order+2)/2; }
    static INLINE int PolDimension (IVec<2> order) { return (order[0]+1)*(order[0]+2)/2; }
    static INLINE int PolBubbleDimension (IVec<2> order) { return (order[0] <= 2) ? 0 : (order[0]-1)*(order[0]-2)/2; }

    static INLINE ELEMENT_TYPE FaceType(int i) { return ET_TRIG; }
    
    template <typename Tx, typename Tlam>
    static INLINE void CalcLambda (const Tx & x, Tlam & lam)
    { lam[0] = x[0]; lam[1] = x[1], lam[2] = 1-x[0]-x[1]; }

    
    // static constexpr IVec<2> edges[3] = { { 2, 0 }, { 1, 2 }, { 0, 1 } };
    static constexpr array<IVec<2>,3> edges = { IVec<2>{ 2, 0 }, IVec<2>{ 1, 2 }, IVec<2>{ 0, 1 } };
    
    static constexpr IVec<2> GetEdge (int i)
    {
      constexpr IVec<2> edges[] =
        { { 2, 0 },
	  { 1, 2 },
	  { 0, 1 } };
      return edges[i];
    }

    template <typename TVN>
    static INLINE IVec<2> GetEdgeSort (int i, const TVN & vnums)
    {
      IVec<2> e = GetEdge (i);
      if (vnums[e[0]] > vnums[e[1]]) Swap (e[0], e[1]);
      return e;
    }


    static INLINE IVec<4> GetFace (int /* i */ )
    {
      static const int face[] = { 0, 1, 2, -1 };
      return IVec<4> (face[0], face[1], face[2], -1);
    }

    template <typename TVN>
    static INLINE IVec<4> GetFaceSort (int /*  i */ , const TVN & vnums)
    {
      IVec<4> f = GetFace (0);
      if(vnums[f[0]] > vnums[f[1]]) Swap(f[0],f[1]); 
      if(vnums[f[1]] > vnums[f[2]]) Swap(f[1],f[2]);
      if(vnums[f[0]] > vnums[f[1]]) Swap(f[0],f[1]); 	

      return f;
    }


    template <typename TVN>
    static INLINE int GetClassNr (const TVN & vnums)
    {
      int classnr = 0;
      int sort[3] = { 0, 1, 2 };
      if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 1; }
      if (vnums[sort[1]] > vnums[sort[2]]) { Swap (sort[1], sort[2]); classnr += 2; }
      if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 2; } // tricky !
      return classnr;
    }

    template <typename TVN>
    static INLINE int GetFacetClassNr (int facet, const TVN & vnums)
    {
      int sort[3] = { 0, 1, 2 };
      if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); }
      if (vnums[sort[1]] > vnums[sort[2]]) { Swap (sort[1], sort[2]); }
      if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); }
      
      static
	int f2vop[] = { 1, 0, 2 };
      int vop = f2vop[facet];
      for (int i = 0; i < 3; i++)
	if (vop == sort[i])
	  return i;
      return -1;  // not possible

      // return GetClassNr (vnums) * 3 + facet;
    } 
  };

  // ******************** QUAD ***********************************
  
  template<> class ET_trait<ET_QUAD>
  {
  public:
    enum { DIM = 2 };
    enum { N_VERTEX = 4 };
    enum { N_EDGE = 4 };
    enum { N_FACE = 1 };
    enum { N_CELL = 0 };
    enum { N_FACET = 4 };

    static constexpr ELEMENT_TYPE ElementType() { return ET_QUAD; }
    constexpr operator ELEMENT_TYPE() const { return ET_QUAD; }
    
    static constexpr INLINE int PolDimension (int order) { return (order+1)*(order+1); }
    static INLINE int PolDimension (IVec<2> order) { return (order[0]+1)*(order[1]+1); }
    static INLINE int PolBubbleDimension (IVec<2> order) { return (order[0] <= 1 || order[1] <= 1) ? 0 : (order[0]-1)*(order[1]-1); }


    static INLINE ELEMENT_TYPE FaceType(int i) { return ET_QUAD; }




    static INLINE IVec<2> GetEdge (int i)
    {
      static const int edges[][2] = 
	{ { 0, 1 },
	  { 2, 3 },
	  { 3, 0 },
	  { 1, 2 }};
      return IVec<2> (edges[i][0], edges[i][1]);
    }

    template <typename TVN>
    static INLINE IVec<2> GetEdgeSort (int i, const TVN & vnums)
    {
      IVec<2> e = GetEdge (i);
      if (vnums[e[0]] > vnums[e[1]]) swap (e[0], e[1]);
      return e;
    }


    template <typename Tx, typename TVN>
    static INLINE Tx XiEdge (int i, Tx hx[], const TVN & vnums)
    {
      IVec<2> e = GetEdgeSort (i, vnums);
      static const int vi[4][2] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
      Tx xi(0.0);
      for (int j = 0; j < 2; j++)
        {
          int edir = vi[e[1]][j] - vi[e[0]][j];
          if (edir == 1) { xi = 2*hx[j]-1; break; }
          if (edir == -1) { xi = 1-2*hx[j]; break; }
        }
      return xi;
    }

    template <typename Tx>
    static INLINE Tx LamEdge (int i, Tx hx[])
    {
      switch (i)
        {
        case 0: return 1-hx[1];
        case 1: return hx[1];
        case 2: return 1-hx[0];
        case 3: default: return hx[0];
        }
    }
 

    static INLINE IVec<4> GetFace (int /* i */ )
    {
      static const int face[] = 
        { 0, 1, 2, 3 };

      return IVec<4> (face[0], face[1], face[2], face[3]);
    }

    template <typename TVN>
    static INLINE IVec<4> GetFaceSort (int /*  i */ , const TVN & vnums)
    {
      IVec<4> f = GetFace (0);
      
      int fmax = 0;
      for (int j=1; j<4; j++) 
        if (vnums[j] < vnums[fmax]) fmax = j;  
      
      int f1 = (fmax+3)%4;
      int f2 = (fmax+1)%4; 
      int fop = (fmax+2)%4; 
      
      if(vnums[f2]<vnums[f1]) swap(f1,f2);  // fmax > f1 > f2 

      f[0] = fmax;
      f[1] = f1;
      f[2] = fop;
      f[3] = f2;

      return f;
    }

    template <typename Tx, typename TVN>
    static INLINE Vec<2,Tx> XiFace (int /* i */, Tx hx[], const TVN & vnums)
    {
      IVec<4> f = GetFaceSort (0, vnums); 
      static const int vi[4][2] = { { 0, 0 }, { 1, 0 }, { 1, 1 }, { 0, 1 } };
      Tx xi(0), eta(0);
      for (int j = 0; j < 2; j++)
        {
          int edir = vi[f[0]][j] - vi[f[1]][j];
          if (edir == 1) { xi = 2*hx[j]-1; break; }
          if (edir == -1) { xi = 1-2*hx[j]; break; }
        }
      for (int j = 0; j < 2; j++)
        {
          int edir = vi[f[0]][j] - vi[f[3]][j];
          if (edir == 1) { eta = 2*hx[j]-1; break; }
          if (edir == -1) { eta = 1-2*hx[j]; break; }
        }
      return Vec<2,Tx> (xi, eta);
    }

    template <typename TVN>
    static INLINE int GetClassNr (const TVN & vnums)
    {
      int classnr = 0;

      int sort[4] = { 0, 1, 2, 3 };
      if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 1; }
      if (vnums[sort[2]] > vnums[sort[3]]) { Swap (sort[2], sort[3]); classnr += 2; }
      if (vnums[sort[0]] > vnums[sort[2]]) { Swap (sort[0], sort[2]); classnr += 4; }
      if (vnums[sort[1]] > vnums[sort[3]]) { Swap (sort[1], sort[3]); classnr += 8; }
      if (vnums[sort[1]] > vnums[sort[2]]) { Swap (sort[1], sort[2]); classnr += 16; }

      return classnr;
    }

    
    template <typename TVN>
    static INLINE int GetFacetClassNr (int facet, const TVN & vnums)
    {
      return GetClassNr (vnums) * 4 + facet;
    }
  };


  template<> class ET_trait<ET_TET>
  {
  public:
    enum { DIM = 3 };
    enum { N_VERTEX = 4 };
    enum { N_EDGE = 6 };
    enum { N_FACE = 4 };
    enum { N_CELL = 1 };
    enum { N_FACET = 4 };

    static constexpr ELEMENT_TYPE ElementType() { return ET_TET; }
    constexpr operator ELEMENT_TYPE() const { return ET_TET; }
    
    static constexpr INLINE int PolDimension (int p) { return (p+1)*(p+2)*(p+3)/6; }
    static INLINE int PolDimension (IVec<3> p) { return (p[0]+1)*(p[0]+2)*(p[0]+3)/6; }
    static INLINE int PolBubbleDimension (IVec<3> p) { return (p[0] <= 3) ? 0 : (p[0]-1)*(p[0]-2)*(p[0]-3)/6;  }

    static INLINE ELEMENT_TYPE FaceType(int i) { return ET_TRIG; }

    static INLINE IVec<2> GetEdge (int i)
    {
      static const int edges[6][2] = 
	{ { 3, 0 },
	  { 3, 1 },
	  { 3, 2 }, 
	  { 0, 1 }, 
	  { 0, 2 },
	  { 1, 2 }};
      return IVec<2> (edges[i][0], edges[i][1]);
    }

    template <typename TVN>
    static INLINE IVec<2> GetEdgeSort (int i, const TVN & vnums)
    {
      IVec<2> e = GetEdge (i);
      if (vnums[e[0]] > vnums[e[1]]) swap (e[0], e[1]);
      return e;
    }


    static INLINE IVec<4> GetFace (int i )
    {
      static const int faces[][4]  =
	{ { 3, 1, 2, -1 },
	  { 3, 2, 0, -1 },
	  { 3, 0, 1, -1 },
	  { 0, 2, 1, -1 } }; 

      return IVec<4> (faces[i][0], faces[i][1], faces[i][2], -1);
    }

    template <typename TVN>
    static INLINE IVec<4> GetFaceSort (int i, const TVN & vnums)
    {
      IVec<4> f = GetFace (i);
      if(vnums[f[0]] > vnums[f[1]]) swap(f[0],f[1]); 
      if(vnums[f[1]] > vnums[f[2]]) swap(f[1],f[2]);
      if(vnums[f[0]] > vnums[f[1]]) swap(f[0],f[1]); 	

      return f;
    }


    template <typename TVN>
    static INLINE int GetClassNr (const TVN & vnums)
    {
      int classnr = 0;

      int sort[4] = { 0, 1, 2, 3 };
      if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 1; }
      if (vnums[sort[2]] > vnums[sort[3]]) { Swap (sort[2], sort[3]); classnr += 2; }
      if (vnums[sort[0]] > vnums[sort[2]]) { Swap (sort[0], sort[2]); classnr += 4; }
      if (vnums[sort[1]] > vnums[sort[3]]) { Swap (sort[1], sort[3]); classnr += 8; }
      if (vnums[sort[1]] > vnums[sort[2]]) { Swap (sort[1], sort[2]); classnr += 16; }

      return classnr;
    }

    template <typename TVN>
    static INLINE int GetFacetClassNr (int facet, const TVN & vnums)
    {
      int sort[4] = { 0, 1, 2, 3 };
      if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); }
      if (vnums[sort[2]] > vnums[sort[3]]) { Swap (sort[2], sort[3]); }
      if (vnums[sort[0]] > vnums[sort[2]]) { Swap (sort[0], sort[2]); }
      if (vnums[sort[1]] > vnums[sort[3]]) { Swap (sort[1], sort[3]); }
      if (vnums[sort[1]] > vnums[sort[2]]) { Swap (sort[1], sort[2]); }

      for (int i = 0; i < 4; i++)
	if (facet == sort[i])
	  return i;
      return -1;  // not possible
      
      // return GetClassNr (vnums) * 4 + facet;
    }


  };

  template<> class ET_trait<ET_PRISM>
  {
  public:
    enum { DIM = 3 };
    enum { N_VERTEX = 6 };
    enum { N_EDGE = 9 };
    enum { N_FACE = 5 };
    enum { N_CELL = 1 };
    enum { N_FACET = 5 };

    static constexpr ELEMENT_TYPE ElementType() { return ET_PRISM; }
    constexpr operator ELEMENT_TYPE() const { return ET_PRISM; }
    
    static INLINE int PolDimension (IVec<3> order) { return (order[0]+1)*(order[0]+2)*(order[2]+1)/2; }
    static INLINE int PolBubbleDimension (IVec<3> p) { return (p[0] <= 2) ? 0 : (p[0]-1)*(p[0]-2)*(p[2]-1)/2; }

    static INLINE ELEMENT_TYPE FaceType(int i) { return (i < 2) ? ET_TRIG : ET_QUAD; }


    static INLINE IVec<2> GetEdge (int i)
    {
#ifndef __CUDA_ARCH__
      static 
#endif
	const int edges[][2] = 
	{ { 2, 0 },
	  { 0, 1 },
	  { 2, 1 },
	  { 5, 3 },
	  { 3, 4 },
	  { 5, 4 },
	  { 2, 5 },
	  { 0, 3 },
	  { 1, 4 }};
      return IVec<2> (edges[i][0], edges[i][1]);
    }

    template <typename TVN>
    static INLINE IVec<2> GetEdgeSort (int i, const TVN & vnums)
    {
      IVec<2> e = GetEdge (i);
      if (vnums[e[0]] > vnums[e[1]]) swap (e[0], e[1]);
      return e;
    }


    static INLINE IVec<4> GetFace (int i )
    {
#ifndef __CUDA_ARCH__
      static 
#endif
	const int faces[][4]  =
	{
	  { 0, 2, 1, -1 },
	  { 3, 4, 5, -1 },
	  { 2, 0, 3, 5 },
	  { 0, 1, 4, 3 },
	  { 1, 2, 5, 4 } 
	};

      return IVec<4> (faces[i][0], faces[i][1], faces[i][2], faces[i][3]);
    }

    template <typename TVN>
    static INLINE IVec<4> GetFaceSort (int i, const TVN & vnums)
    {
      IVec<4> f = GetFace (i);
      if (i < 2)
	{
	  if(vnums[f[0]] > vnums[f[1]]) swap(f[0],f[1]); 
	  if(vnums[f[1]] > vnums[f[2]]) swap(f[1],f[2]);
	  if(vnums[f[0]] > vnums[f[1]]) swap(f[0],f[1]); 	
	  return f;
	}
      else
	{
	  int fmax = 0;
	  for (int j=1; j<4; j++) 
	    if (vnums[f[j]] < vnums[f[fmax]]) fmax = j;  
	  
	  int f1 = (fmax+3)%4;
	  int f2 = (fmax+1)%4; 
	  int fop = (fmax+2)%4; 
	  
	  if(vnums[f[f2]]<vnums[f[f1]]) swap(f1,f2);  // fmax > f1 > f2 
	  
	  return IVec<4> (f[fmax], f[f1], f[fop], f[f2]);
	}
    }

    template <typename TVN>
    static INLINE int GetClassNr (const TVN & vnums)
    {
      int classnr = 0;
      int sort[3] = { 0, 1, 2 };
      if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 1; }
      if (vnums[sort[1]] > vnums[sort[2]]) { Swap (sort[1], sort[2]); classnr += 2; }
      if (vnums[sort[0]] > vnums[sort[1]]) { Swap (sort[0], sort[1]); classnr += 2; } // tricky !
      return classnr;
    }

    template <typename TVN>
    static INLINE int GetFacetClassNr (int facet, const TVN & vnums)
    {
      return 0;
    }

  };

  template<> class ET_trait<ET_PYRAMID>
  {
  public:
    enum { DIM = 3 };
    enum { N_VERTEX = 5 };
    enum { N_EDGE = 8 };
    enum { N_FACE = 5 };
    enum { N_CELL = 1 };
    enum { N_FACET = 5 };

    static constexpr ELEMENT_TYPE ElementType() { return ET_PYRAMID; }
    constexpr operator ELEMENT_TYPE() const { return ET_PYRAMID; }
    
    static INLINE ELEMENT_TYPE FaceType(int i) { return (i < 4) ? ET_TRIG : ET_QUAD; }

    static INLINE int PolDimension (IVec<3> order) { return (order[0]+2)*(order[0]+1)*(2*order[0]+3) / 6; }
    static INLINE int PolBubbleDimension (IVec<3> p) { return (p[0] <= 2) ? 0 :  (p[0]-1)*(p[0]-2)*(2*p[0]-3)/6; }

    static INLINE IVec<2> GetEdge (int i)
    {
#ifndef __CUDA_ARCH__
      static 
#endif
	const int edges[][2] = 
	{ { 0, 1 },
	  { 1, 2 },
	  { 0, 3 },
	  { 3, 2 },
	  { 0, 4 },
	  { 1, 4 },
	  { 2, 4 },
	  { 3, 4 }};

      return IVec<2> (edges[i][0], edges[i][1]);
    }

    template <typename TVN>
    static INLINE IVec<2> GetEdgeSort (int i, const TVN & vnums)
    {
      IVec<2> e = GetEdge (i);
      if (vnums[e[0]] > vnums[e[1]]) swap (e[0], e[1]);
      return e;
    }



    static INLINE IVec<4> GetFace (int i )
    {
#ifndef __CUDA_ARCH__
      static 
#endif
	const int faces[][4]  =
	{
	  { 0, 1, 4, -1 },
	  { 1, 2, 4, -1 },
	  { 2, 3, 4, -1 },
	  { 3, 0, 4, -1 },
	  { 0, 3, 2, 1 } 
	};

      return IVec<4> (faces[i][0], faces[i][1], faces[i][2], faces[i][3]);
    }

    template <typename TVN>
    static INLINE IVec<4> GetFaceSort (int i, const TVN & vnums)
    {
      IVec<4> f = GetFace (i);
      if (i < 4)
	{
	  if(vnums[f[0]] > vnums[f[1]]) swap(f[0],f[1]); 
	  if(vnums[f[1]] > vnums[f[2]]) swap(f[1],f[2]);
	  if(vnums[f[0]] > vnums[f[1]]) swap(f[0],f[1]); 	
	  return f;
	}
      else
	{
	  int fmax = 0;
	  for (int j=1; j<4; j++) 
	    if (vnums[f[j]] < vnums[f[fmax]]) fmax = j;  
	  
	  int f1 = (fmax+3)%4;
	  int f2 = (fmax+1)%4; 
	  int fop = (fmax+2)%4; 
	  
	  if(vnums[f[f2]]<vnums[f[f1]]) swap(f1,f2);  // fmax > f1 > f2 
	  
	  return IVec<4> (f[fmax], f[f1], f[fop], f[f2]);
	}
    }

    template <typename TVN>
    static INLINE int GetClassNr (const TVN & vnums)
    {
      return 0;
    }

    template <typename TVN>
    static INLINE int GetFacetClassNr (int facet, const TVN & vnums)
    {
      return 0;
    }

  };



  template<> class ET_trait<ET_HEXAMID>
  {
  public:
    enum { DIM = 3 };
    enum { N_VERTEX = 7 };
    enum { N_EDGE = 11 };
    enum { N_FACE = 6 };
    enum { N_CELL = 1 };
    enum { N_FACET = 6 };

    static constexpr ELEMENT_TYPE ElementType() { return ET_HEXAMID; }
    constexpr operator ELEMENT_TYPE() const { return ET_HEXAMID; }
    
    static INLINE ELEMENT_TYPE FaceType(int i) { return (i==1||i==4) ? ET_TRIG : ET_QUAD; }

    // TODO
    static INLINE int PolDimension (IVec<3> order) {
      int p = order[0];
      return (p+1)*(p+1)*(p+1) - (p-1) - p*(p-1);
    }
    // TODO
    static INLINE int PolBubbleDimension (IVec<3> p) { return (p[0] < 2) ? 0 :  (p[0]-1)*(p[0]-1)*(p[0]-1); }

    static INLINE IVec<2> GetEdge (int i)
    {
#ifndef __CUDA_ARCH__
      static 
#endif
	const int edges[][2] = 
	{
          { 0, 1 },
          { 2, 3 },
          { 3, 0 },
          { 1, 2 },
          { 4, 5 },
          { 6, 4 },
          { 5, 6 },
          { 0, 4 },
          { 1, 5 },
          { 2, 6 },
          { 3, 6 },
        };
      return IVec<2> (edges[i][0], edges[i][1]);
    }

    template <typename TVN>
    static INLINE IVec<2> GetEdgeSort (int i, const TVN & vnums)
    {
      IVec<2> e = GetEdge (i);
      if (vnums[e[0]] > vnums[e[1]]) swap (e[0], e[1]);
      return e;
    }



    static INLINE IVec<4> GetFace (int i )
    {
#ifndef __CUDA_ARCH__
      static 
#endif
	const int faces[][4]  =
	{
          { 0, 3, 2, 1 },
          { 4, 5, 6, -1},
          { 0, 1, 5, 4 },
          { 1, 2, 6, 5 },
          { 2, 3, 6, -1},
          { 3, 0, 4, 6 }
	};

      return IVec<4> (faces[i][0], faces[i][1], faces[i][2], faces[i][3]);
    }

    // TODO
    template <typename TVN>
    static INLINE IVec<4> GetFaceSort (int i, const TVN & vnums)
    {
      IVec<4> f = GetFace (i);
      if (f[3] < 0)
	{
	  if(vnums[f[0]] > vnums[f[1]]) swap(f[0],f[1]); 
	  if(vnums[f[1]] > vnums[f[2]]) swap(f[1],f[2]);
	  if(vnums[f[0]] > vnums[f[1]]) swap(f[0],f[1]); 	
	  return f;
	}
      else
	{
	  int fmax = 0;
	  for (int j=1; j<4; j++) 
	    if (vnums[f[j]] < vnums[f[fmax]]) fmax = j;  
	  
	  int f1 = (fmax+3)%4;
	  int f2 = (fmax+1)%4; 
	  int fop = (fmax+2)%4; 
	  
	  if(vnums[f[f2]]<vnums[f[f1]]) swap(f1,f2);  // fmax > f1 > f2 
	  
	  return IVec<4> (f[fmax], f[f1], f[fop], f[f2]);
	}
    }

    template <typename TVN>
    static INLINE int GetClassNr (const TVN & vnums)
    {
      return 0;
    }

    template <typename TVN>
    static INLINE int GetFacetClassNr (int facet, const TVN & vnums)
    {
      return 0;
    }

  };



  

  template<> class ET_trait<ET_HEX>
  {
  public:
    enum { DIM = 3 };
    enum { N_VERTEX = 8 };
    enum { N_EDGE = 12 };
    enum { N_FACE = 6 };
    enum { N_CELL = 1 };
    enum { N_FACET = 6 };

    static constexpr ELEMENT_TYPE ElementType() { return ET_HEX; }
    constexpr operator ELEMENT_TYPE() const { return ET_HEX; }
    
    static ELEMENT_TYPE FaceType(int i) { return ET_QUAD; }

    static inline int PolDimension (IVec<3> order) { return (order[0]+1)*(order[1]+1)*(order[2]+1); }
    static inline int PolBubbleDimension (IVec<3> p) { return (p[0] <= 1) ? 0 : (p[0]-1)*(p[1]-1)*(p[2]-1); }

    static IVec<2> GetEdge (int i)
    {
      static const int edges[][2] = 
	{
	  { 0, 1 },
	  { 2, 3 },
	  { 3, 0 },
	  { 1, 2 },
	  { 4, 5 },
	  { 6, 7 },
	  { 7, 4 },
	  { 5, 6 },
	  { 0, 4 },
	  { 1, 5 },
	  { 2, 6 },
	  { 3, 7 },
	};
      return IVec<2> (edges[i][0], edges[i][1]);
    }

    template <typename TVN>
    static IVec<2> GetEdgeSort (int i, const TVN & vnums)
    {
      IVec<2> e = GetEdge (i);
      if (vnums[e[0]] > vnums[e[1]]) swap (e[0], e[1]);
      return e;
    }


    static IVec<4> GetFace (int i )
    {
      static const int faces[][4]  =
	{
	  { 0, 3, 2, 1 },
	  { 4, 5, 6, 7 },
	  { 0, 1, 5, 4 },
	  { 1, 2, 6, 5 },
	  { 2, 3, 7, 6 },
	  { 3, 0, 4, 7 }
	};

      return IVec<4> (faces[i][0], faces[i][1], faces[i][2], faces[i][3]);
    }


    template <typename TVN>
    static IVec<4> GetFaceSort (int i, const TVN & vnums)
    {
      IVec<4> f = GetFace (i);

      int fmax = 0;
      for (int j=1; j<4; j++) 
	if (vnums[f[j]] < vnums[f[fmax]]) fmax = j;  
      
      int f1 = (fmax+3)%4;
      int f2 = (fmax+1)%4; 
      int fop = (fmax+2)%4; 
      
      if(vnums[f[f2]]<vnums[f[f1]]) swap(f1,f2);  // fmax > f1 > f2 
      
      return IVec<4> (f[fmax], f[f1], f[fop], f[f2]);
    }


    template <typename TVN>
    static int GetClassNr (const TVN & vnums)
    {
      int verts[8];
      for (int j = 0; j < 8; j++)
        verts[j] = vnums[j];

      int classnr = 0;
      for (int j = 0; j < 7; j++)
        {
          int maxk = 0;
          for (int k = 0; k < 8-j; k++)
            if (verts[k] > verts[maxk]) maxk = k;
          // compress
          for (int k = maxk; k < 7-j; k++)
            verts[k] = verts[k+1];
          classnr = maxk + (8-j) * classnr;
        }      
      return classnr;
    }

    template <typename TVN>
    static int GetFacetClassNr (int facet, const TVN & vnums)
    {
      return 0;
    }

  };

  
  template <typename T>
  INLINE int PolBubbleDimension (ELEMENT_TYPE ET, IVec<1,T> order)
  {
    switch (ET)
      {
      case ET_SEGM: return ET_trait<ET_SEGM>::PolBubbleDimension (order);
      default: return 0;
      }
  }

  template <typename T>
  inline int PolBubbleDimension (ELEMENT_TYPE ET, IVec<2,T> order)
  {
    switch (ET)
      {
      case ET_TRIG: return ET_trait<ET_TRIG>::PolBubbleDimension (order);
      case ET_QUAD: return ET_trait<ET_QUAD>::PolBubbleDimension (order);
      default: return 0;
      }
  }

  template <typename T>
  inline int PolBubbleDimension (ELEMENT_TYPE ET, IVec<3,T> order)
  {
    switch (ET)
      {
      case ET_TET: return ET_trait<ET_TET>::PolBubbleDimension (order);
      case ET_PRISM: return ET_trait<ET_PRISM>::PolBubbleDimension (order);
      case ET_PYRAMID: return ET_trait<ET_PYRAMID>::PolBubbleDimension (order);
      case ET_HEX: return ET_trait<ET_HEX>::PolBubbleDimension (order);
      default: return 0;
      }
  }

  template <typename FUNC>
  decltype(auto) SwitchET (ELEMENT_TYPE et, FUNC f)
  {
    switch (et)
      {
      case ET_POINT:   return f(ET_trait<ET_POINT>()); 
      case ET_SEGM:    return f(ET_trait<ET_SEGM>()); 
      case ET_TRIG:    return f(ET_trait<ET_TRIG>()); 
      case ET_QUAD:    return f(ET_trait<ET_QUAD>()); 
      case ET_TET:     return f(ET_trait<ET_TET>()); 
      case ET_PRISM:   return f(ET_trait<ET_PRISM>()); 
      case ET_PYRAMID: return f(ET_trait<ET_PYRAMID>()); 
      case ET_HEXAMID: return f(ET_trait<ET_HEXAMID>()); 
      case ET_HEX:     return f(ET_trait<ET_HEX>()); 
      default:
        unreachable();
      }
  }

  template<ELEMENT_TYPE ET1, typename FUNC>
  decltype(auto) SwitchET (ELEMENT_TYPE et, FUNC f)
  {
    if (et != ET1)
      throw Exception("Element type not defined!");
    return f(ET_trait<ET1>());
  }
  
  template<ELEMENT_TYPE ET1, ELEMENT_TYPE ET2, ELEMENT_TYPE ... ET_REST, typename FUNC>
  decltype(auto) SwitchET (ELEMENT_TYPE et, FUNC f)
  {
    if (et==ET1)
      return f(ET_trait<ET1>());
    else
      return SwitchET<ET2,ET_REST...>(et,f);
  }
}

#endif
