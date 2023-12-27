/*********************************************************************/
/* File:   topology.cpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   Toplogy of reference elements
*/


#include "elementtopology.hpp"

namespace ngfem
{


  ostream & operator<< (ostream & ost, ELEMENT_TYPE et)
  {
    ost << ElementTopology::GetElementName (et);
    return ost;
  }

  ostream & operator<< (ostream & ost, NODE_TYPE nt)
  {
    switch (nt)
      {
      case NT_VERTEX: ost << "Vertex"; break;
      case NT_EDGE: ost << "Edge"; break;
      case NT_FACE: ost << "Face"; break;
      case NT_CELL: ost << "Cell"; break;
      case NT_ELEMENT: ost << "Element"; break;
      case NT_FACET: ost << "Facet"; break;
      case NT_GLOBAL: ost << "Global"; break;
      }
    return ost;
  }


  const char * ElementTopology :: GetElementName (ELEMENT_TYPE et)
  {
    switch (et)
      {
      case ET_POINT: return "Point";
      case ET_SEGM: return "Segm";
      case ET_TRIG: return "Trig";
      case ET_QUAD: return "Quad";
      case ET_TET:  return "Tet";
      case ET_PYRAMID: return "Pyramid";
      case ET_PRISM: return "Prism";
      case ET_HEXAMID:  return "Hexamid";
      case ET_HEX:  return "Hex";
      }
    throw Exception("illegal element type");
  }

  const POINT3D * ElementTopology :: GetVertices (ELEMENT_TYPE et)
  {
    static double segm_points [][3] = 
      { { 1 },
	{ 0 } };
    
    static double trig_points [][3] = 
      { { 1, 0 },
	{ 0, 1 },
	{ 0, 0 } };

    static double quad_points [][3] = 
      { { 0, 0 },
	{ 1, 0 },
	{ 1, 1 },
	{ 0, 1 } };

    static double tet_points [][3] = 
      { { 1, 0, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1 },
	{ 0, 0, 0 } };

    static double pyramid_points [][3] =
      {
	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{ 1, 1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1-1e-12 },
      };    
  
    static double prism_points[][3] = 
      {
	{ 1, 0, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 0 },
	{ 1, 0, 1 },
	{ 0, 1, 1 },
	{ 0, 0, 1 }
      };

    static double hexamid_points[][3] = 
      { 
	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{ 1, 1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1 },
	{ 1, 0, 1 },
	{ 0, 1, 1 }
      };

    static double hex_points[][3] = 
      { 
	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{ 1, 1, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1 },
	{ 1, 0, 1 },
	{ 1, 1, 1 },
	{ 0, 1, 1 }
      };

    static double point_points[][3] = 
      { 
	{ 0, 0, 0 },
      };

    
    switch (et)
      {
      case ET_POINT: return point_points;
      case ET_SEGM: return segm_points;
      case ET_TRIG: return trig_points;
      case ET_QUAD: return quad_points;
      case ET_TET:  return tet_points;
      case ET_PYRAMID: return pyramid_points;
      case ET_PRISM: return prism_points;
      case ET_HEXAMID: return hexamid_points;
      case ET_HEX: return hex_points;        
      default:
        break;
      }

    stringstream str;
    str << "Ng_GetVertices, illegal element type " << et << "\n";
    throw Exception (str.str());
  }



#ifdef XXX
  const EDGE * ElementTopology::GetEdges (ELEMENT_TYPE et)
  {
    static int segm_edges[1][2] =
      { { 0, 1 }};
    

    static int trig_edges[3][2] =
      { { 2, 0 },
	{ 1, 2 },
	{ 0, 1 }};

    /*
    static int trig_edges[3][2] =
      { { 1, 2 },
	{ 2, 0 },
	{ 0, 1 }};
    */

    static int quad_edges[4][2] =
      { { 0, 1 },
	{ 2, 3 },
	{ 3, 0 },
	{ 1, 2 }};
    
    static int tet_edges[6][2] =
      { { 3, 0 },
	{ 3, 1 },
	{ 3, 2 }, 
	{ 0, 1 }, 
	{ 0, 2 },
	{ 1, 2 }};
    
    static int prism_edges[9][2] =
      { { 2, 0 },
	{ 0, 1 },
	{ 2, 1 },
	{ 5, 3 },
	{ 3, 4 },
	{ 5, 4 },
	{ 2, 5 },
	{ 0, 3 },
	{ 1, 4 }};

    static int pyramid_edges[8][2] =
      { { 0, 1 },
	{ 1, 2 },
	{ 0, 3 },
	{ 3, 2 },
	{ 0, 4 },
	{ 1, 4 },
	{ 2, 4 },
	{ 3, 4 }};

  static int hex_edges[12][2] =
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
      case ET_SEGM: return segm_edges;
      case ET_TRIG: return trig_edges;
      case ET_QUAD: return quad_edges;
      case ET_TET:  return tet_edges;
      case ET_PYRAMID: return pyramid_edges;
      case ET_PRISM: return prism_edges;
      case ET_HEX: return hex_edges;
      default:
	break;
      }
    cerr << "Ng_GetEdges, illegal element type " << et << endl;
    return 0;  
  }


  const FACE * ElementTopology::GetFaces (ELEMENT_TYPE et)
  {
    static int tet_faces[4][4] =
      { { 3, 1, 2, -1 },
	{ 3, 2, 0, -1 },
	{ 3, 0, 1, -1 },
        { 0, 2, 1, -1 } }; // all faces point into interior!
  
    /*
      // orig
    static int prism_faces[5][4] =
      {
	{ 0, 1, 2, -1 },
	{ 3, 5, 4, -1 },
	{ 2, 0, 3, 5 },
	{ 0, 1, 4, 3 },
	{ 1, 2, 5, 4 } 
      };
    */

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
      case ET_HEX: return hex_faces;

      case ET_TRIG: return trig_faces;
      case ET_QUAD: return quad_faces;

      case ET_SEGM:
      default:
	break;
      }
    
    cerr << "Ng_GetFaces, illegal element type " << et << endl;
    return 0;
  }
#endif  








  
  NORMAL * ElementTopology::GetNormals(ELEMENT_TYPE et)
  {
    // length of normal: 
    //  in 2d:  length of edge (sum of intweights = 1)
    //  in 3d:  area of parallelogramm (sum of intweights = 1(quad), 1/2 trig)
    // this is property is preserved by piola transformation: n = det(F) * F^(-T) n_ref

    static double segm_normals [][3] = 
    { { 1, 0, 0 },
      { -1, 0, 0 } };

    static double trig_normals [][3] = 
    { { 0, -1, 0 },
    { -1, 0, 0 },
    { 1, 1, 0 } };

    static double quad_normals [][3] = 
    { { 0, -1, 0 },
    { 0, 1, 0 },
    { -1, 0, 0 },
    { 1, 0, 0 } };

    static double tet_normals [][3] = 
    { { -1, 0, 0 },
    { 0, -1, 0 },
    { 0, 0, -1 },
    { 1, 1, 1 } };

    static double pyramid_normals [][3] =
    {
      { 0, -1, 0 },
      { 1, 0, 1 },
      { 0, 1, 1 },
      { -1, 0, 0 },
      { 0, 0, -1 }, 
    };    
  
    static double prism_normals[][3] = 
    {
      { 0, 0, -1 },
      { 0, 0, 1 },
      { 0, -1, 0 },
      { 1, 1, 0 },
      { -1, 0, 0 },
    };

    static double hex_normals[][3] = 
    { 
      { 0, 0, -1 },
      { 0, 0, 1 },
      { 0, -1, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { -1, 0, 0 },
    };

    
    switch (et)
    {
      case ET_SEGM: return segm_normals;
      case ET_TRIG: return trig_normals;
      case ET_QUAD: return quad_normals;
      case ET_TET:  return tet_normals;
      case ET_PYRAMID: return pyramid_normals;
      case ET_PRISM: return prism_normals;
      case ET_HEX: return hex_normals;
      default:
        break;
    }
  
    stringstream str;
    str << "Ng_GetNormals, illegal element type " << et << "\n";
    throw Exception (str.str());
  }

  template <int D>
  FlatVector< Vec<D> > ElementTopology::GetNormals(ELEMENT_TYPE et)
  {
    // length of normal: 
    //  in 2d:  length of edge (sum of intweights = 1)
    //  in 3d:  area of parallelogramm (sum of intweights = 1(quad), 1/2 trig)
    // this is property is preserved by piola transformation: n = det(F) * F^(-T) n_ref

    static Vec<0> point_normals [] = 
      { { } };  // only to avoid 0-size array
    
    static Vec<1> segm_normals [] = 
    { { 1 },
      { -1 } };
    
    static Vec<2> trig_normals [] = 
    { { 0, -1 },
    { -1, 0 },
    { 1, 1 } };

    static Vec<2> quad_normals [] = 
    { { 0, -1 },
    { 0, 1 },
    { -1, 0 },
    { 1, 0 } };

    static Vec<3> tet_normals [] = 
    { { -1, 0, 0 },
    { 0, -1, 0 },
    { 0, 0, -1 },
    { 1, 1, 1 } };

    static Vec<3> pyramid_normals []=
    {
      { 0, -1, 0 },
      { 1, 0, 1 },
      { 0, 1, 1 },
      { -1, 0, 0 },
      { 0, 0, -1 }, 
    };    
  
    static Vec<3> prism_normals[] = 
    {
      { 0, 0, -1 },
      { 0, 0, 1 },
      { 0, -1, 0 },
      { 1, 1, 0 },
      { -1, 0, 0 },
    };

    static Vec<3> hex_normals[] = 
    { 
      { 0, 0, -1 },
      { 0, 0, 1 },
      { 0, -1, 0 },
      { 1, 0, 0 },
      { 0, 1, 0 },
      { -1, 0, 0 },
    };


    if constexpr (D==0)
      switch (et)
        {
        case ET_POINT: return FlatVector< Vec<D> > (0, point_normals);
        default: ;
        }
        
    if constexpr (D==1)
      switch (et)
        {
        case ET_SEGM: return FlatVector< Vec<D> > (2, segm_normals);
        default: ;
        }
    
    if constexpr (D==2)
      switch (et)
        {
        case ET_TRIG: return FlatVector< Vec<D> > (3, trig_normals);
        case ET_QUAD: return FlatVector< Vec<D> > (4, quad_normals);
        default: ;
        }
    
    if constexpr (D==3)
      switch (et)
        {
        case ET_TET:  return FlatVector< Vec<D> > (4, tet_normals);
        case ET_PYRAMID: return FlatVector< Vec<D> > (5, pyramid_normals);
        case ET_PRISM: return FlatVector< Vec<D> > (5, prism_normals);
        case ET_HEX: return FlatVector< Vec<D> > (6, hex_normals);
        default: ;
        }

    stringstream str;
    str << "Ng_GetNormals, illegal element type " << et << "\n";
    throw Exception (str.str());
  }
  

  template NGS_DLL_HEADER FlatVector< Vec<0> > ElementTopology::GetNormals<0> (ELEMENT_TYPE et);  
  template NGS_DLL_HEADER FlatVector< Vec<1> > ElementTopology::GetNormals<1> (ELEMENT_TYPE et);
  template NGS_DLL_HEADER FlatVector< Vec<2> > ElementTopology::GetNormals<2> (ELEMENT_TYPE et);
  template NGS_DLL_HEADER FlatVector< Vec<3> > ElementTopology::GetNormals<3> (ELEMENT_TYPE et);

  
  int ElementTopology :: GetEdgeNr (ELEMENT_TYPE et, int v1, int v2)
  {
    const EDGE * edges = GetEdges (et);
    int nedges = GetNEdges (et);

    for (int i = 0; i < nedges; i++)
      {
	if (edges[i][0] == v1 && edges[i][1] == v2) return i;
	if (edges[i][1] == v1 && edges[i][0] == v2) return i;
      }

    stringstream str;
    str << "no element edge, eltype = " << et << ", nedges = " << nedges << ", v1,2 = " << v1 << ", " << v2 << endl;
    throw Exception (str.str());
  }
  


  int ElementTopology :: GetFaceNr (ELEMENT_TYPE et, int v1, int v2, int v3)
  {
    const FACE * faces = GetFaces (et);
    int nfaces = GetNFaces (et);

    for (int i = 0; i < nfaces; i++)
      {
	if (faces[i][0] == v1 && faces[i][1] == v2 && faces[i][2] == v3) return i;
	if (faces[i][0] == v1 && faces[i][1] == v3 && faces[i][2] == v2) return i;
	if (faces[i][0] == v2 && faces[i][1] == v1 && faces[i][2] == v3) return i;
	if (faces[i][0] == v2 && faces[i][1] == v3 && faces[i][2] == v1) return i;
	if (faces[i][0] == v3 && faces[i][1] == v1 && faces[i][2] == v2) return i;
	if (faces[i][0] == v3 && faces[i][1] == v2 && faces[i][2] == v1) return i;
      }

    stringstream str;
    str << "no element face, eltype = " << et << ", nfaces = " << nfaces << ", v1,2,3 = " << v1 << ", " << v2 << ", " << v3 << endl;
    throw Exception (str.str());
  }

  




  ostream & operator<< (ostream & ost, const NodeId & node)
  {
    switch (node.GetType())
      {
      case 0: ost << "V"; break;
      case 1: ost << "E"; break;
      case 2: ost << "F"; break;
      case 3: ost << "C"; break;
      default: ost << "undef"; break;
      }
    ost << node.GetNr();
    return ost;
  }

  /*
  ostream & operator<< (ostream & ost, const TopologicElement & etop)
  {
    ost << ElementTopology::GetElementName(etop.GetType()) << endl;
    cout << "nd = " << etop.GetNNodes() << endl;
    for (int i = 0; i < etop.GetNNodes(); i++)
      ost << etop.GetNode(i) << endl;
    return ost;
  }
  */


  


}
