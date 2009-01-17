/*********************************************************************/
/* File:   meshaccess.cpp                                            */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   Access to fe mesh
*/


#include <comp.hpp>
#include <parallelngs.hpp>


namespace ngcomp

{
  using namespace ngcomp;

  MeshAccess :: MeshAccess ()
  {
    UpdateBuffers();
  }


  MeshAccess :: ~MeshAccess ()
  {
    ;
  }

  void MeshAccess :: UpdateBuffers()
  {
    dim = Ng_GetDimension();

    if (dim == -1)
      for (int i = 0; i < 4; i++)  
        {
          nnodes[i] = 0;
          nelements[i] = 0;
        }
    else
      for (int i = 0; i < 4; i++)  
        {
          nnodes[i] = Ng_GetNNodes(i);
          nelements[i] = Ng_GetNElements(i);
        }

    for (int i = 0; i <= dim; i++)
      {
        nnodes_cd[i] = nnodes[dim-i];
        nelements_cd[i] = nelements[dim-i];
      }

    nlevels = Ng_GetNLevels(); 
  }

  /*
  int MeshAccess :: GetNNodes (NODE_TYPE nt) const
  {
    // return Ng_GetNNodes (nt);

    switch (nt)
      {
      case NT_VERTEX: return Ng_GetNV();
      case NT_EDGE: return Ng_GetNEdges();
      case NT_FACE: return Ng_GetNFaces();
      case NT_CELL: return Ng_GetNE();
      }
  }
  */


  void MeshAccess::GetTopologicElement (int elnr, TopologicElement & topel) const
  {
    int help[54];
    int nnodes;

    nnodes = Ng_GetElementClosureNodes (dim, elnr, 
                                        NodeSet (NT_VERTEX, NT_EDGE, NT_FACE, NT_CELL), help);
        
    topel.SetElementType (GetElType(elnr));
    topel.Clear();
    for (int i = 0; i < nnodes; i++)
      topel.AddNode (Node (NODE_TYPE (help[2*i]), help[2*i+1]));
    
    /*    
    // 2D case
    ArrayMem<int, 12> nums;
    
    topel.SetElementType (GetElType(elnr));
    topel.Clear();

    GetElVertices (elnr, nums);
    for (int j = 0; j < nums.Size(); j++)
      topel.AddNode (Node (NT_VERTEX, nums[j]));

    GetElEdges (elnr, nums);
    for (int j = 0; j < nums.Size(); j++)
      topel.AddNode (Node (NT_EDGE, nums[j]));

    if (GetDimension() == 3)
      {
        GetElFaces (elnr, nums);
        for (int j = 0; j < nums.Size(); j++)
          topel.AddNode (Node (NT_FACE, nums[j]));
        
        topel.AddNode (Node (NT_CELL, elnr));
      }
    else
      topel.AddNode (Node (NT_FACE, elnr));
    */
  }






  int MeshAccess :: GetNDomains () const
  {
    int maxind = -1;
    int ne = GetNE(); 
    for (int i = 0; i < ne; i++)
      { 
	int ind = GetElIndex(i);
	if (ind > maxind) maxind = ind;
      }
    return maxind+1;
  }


  int MeshAccess :: GetNBoundaries () const
  {
    int maxind = -1;
    int nse = GetNSE(); 
    for (int i = 0; i < nse; i++)
      { 
	int ind = GetSElIndex(i);
	if (ind > maxind) maxind = ind;
      }
    return maxind+1;
  }


  
  ELEMENT_TYPE MeshAccess :: GetElType (int elnr) const
{
  elnr++;
  switch (Ng_GetElementType (elnr))
    {
    case NG_TRIG:
    case NG_TRIG6:
      {
	return ET_TRIG;
	break;
      }

    case NG_QUAD:
    case NG_QUAD6:
      {
	return ET_QUAD;
	break;
      }

    case NG_TET:
    case NG_TET10:
      {
	return ET_TET;
	break;
      }
    case NG_PRISM:
    case NG_PRISM12:
      {
	return ET_PRISM;
	break;
      }
    case NG_PYRAMID:
      {
	return ET_PYRAMID;
	break;
      }
    case NG_HEX:
      {
	return ET_HEX;
	break;
      }

    default:
      {
	cerr << "MeshAccess::GetElType undefined type "
	     << int(Ng_GetElementType (elnr)) << endl;
      }
    }
  return ET_TET;
}


ELEMENT_TYPE MeshAccess :: GetSElType (int elnr) const
{
  int pi[NG_SURFACE_ELEMENT_MAXPOINTS];
  switch (Ng_GetSurfaceElement (elnr+1, pi))
    {
    case NG_SEGM:
    case NG_SEGM3:
      {
	return ET_SEGM;
	break;
      }

    case NG_TRIG:
    case NG_TRIG6:
      {
	return ET_TRIG;
	break;
      }

    case NG_QUAD:
    case NG_QUAD6:
      {
	return ET_QUAD;
	break;
      }
    default:
      {
	throw Exception ("GetSElType: Unhandled element type");
      }
    }


  stringstream str;
  str << "MeshAccess::GetSElType undefined type "
      << int(Ng_GetSurfaceElement (elnr, pi)) << endl;
  throw Exception (str.str());
}



// von astrid
int MeshAccess :: GetElNV ( int elnr ) const 
{
  int nv;
  NG_ELEMENT_TYPE typ = Ng_GetElementType (elnr+1);

  switch (typ)
    {
    case NG_TRIG: nv = 3; break;
    case NG_TRIG6: nv = 6; break;
    case NG_QUAD: nv = 4; break;
    case NG_QUAD6: nv = 6; break;

    case NG_TET: nv = 4; break;
    case NG_TET10: nv = 10; break;
    case NG_PYRAMID: nv = 5; break;
    case NG_PRISM: nv = 6; break;
    case NG_PRISM12: nv = 12; break;
    case NG_HEX: nv = 8; break;

    default:
      {
	stringstream err;
	err << "MeshAccess::GetElNV elnr = " << elnr 
            << " has undefined type " << int(typ) << endl;
	throw Exception (err.str());
      }
    }


  return nv;


} 




void MeshAccess :: GetElPNums (int elnr, ARRAY<int> & pnums) const
{
  pnums.SetSize (NG_ELEMENT_MAXPOINTS);
  NG_ELEMENT_TYPE typ = Ng_GetElement (elnr+1, &pnums[0]);
  
  switch (typ)
    {
    case NG_TRIG: pnums.SetSize(3); break;
    case NG_TRIG6: pnums.SetSize(6); break;
    case NG_QUAD: pnums.SetSize(4); break;
    case NG_QUAD6: pnums.SetSize(6); break;

    case NG_TET: pnums.SetSize(4); break;
    case NG_TET10: pnums.SetSize(10); break;
    case NG_PYRAMID: pnums.SetSize(5); break;
    case NG_PRISM: pnums.SetSize(6); break;
    case NG_PRISM12: pnums.SetSize(12); break;
    case NG_HEX: pnums.SetSize(8); break;

    default:
      {
	stringstream err;
	err << "MeshAccess::GetElPNums elnr = " << elnr 
            << " has undefined type " << int(typ) << endl;
	throw Exception (err.str());
      }
    }

  for (int i = 0; i < pnums.Size(); i++)
    pnums[i]--;
}


void MeshAccess :: GetSegmentPNums (int snr, ARRAY<int> & pnums) const
{
  pnums.SetSize(3);
  int np;
  NG_ELEMENT_TYPE typ = Ng_GetSegment (snr+1, &pnums[0], &np);
  pnums.SetSize(np);

  for (int i = 0; i < np; i++)
    pnums[i]--;
}


int MeshAccess :: GetSegmentIndex (int snr) const
{
  return Ng_GetSegmentIndex(snr+1);
}




void MeshAccess :: GetElVertices (int elnr, ARRAY<int> & vnums) const
{
  vnums.SetSize (NG_ELEMENT_MAXPOINTS);
  NG_ELEMENT_TYPE typ = Ng_GetElement (elnr+1, &vnums[0]);
  
  switch (typ)
    {
    case NG_TRIG: vnums.SetSize(3); break;
    case NG_TRIG6: vnums.SetSize(3); break;
    case NG_QUAD: vnums.SetSize(4); break;
    case NG_QUAD6: vnums.SetSize(4); break;

    case NG_TET: vnums.SetSize(4); break;
    case NG_TET10: vnums.SetSize(4); break;
    case NG_PYRAMID: vnums.SetSize(5); break;
    case NG_PRISM: vnums.SetSize(6); break;
    case NG_PRISM12: vnums.SetSize(6); break;
    case NG_HEX: vnums.SetSize(8); break;

    default:
      {
	stringstream err;
	err << "MeshAccess::GetElVertices undefined type "
	     << int(typ) << endl;
	throw Exception (err.str());
      }
    }

  for (int i = 0; i < vnums.Size(); i++)
    vnums[i]--;
}

  
void MeshAccess :: 
GetElEdges (int elnr, ARRAY<int> & ednums, ARRAY<int> & orient) const
{
  ednums.SetSize (12);
  orient.SetSize (12);
  int ned = 
    Ng_GetElement_Edges (elnr+1, &ednums[0], &orient[0]);
  ednums.SetSize (ned);
  orient.SetSize (ned);
  for (int i = 0; i < ned; i++)
    ednums[i]--;
}

void MeshAccess :: 
GetElEdges (int elnr, ARRAY<int> & ednums) const
{
  ednums.SetSize (12);
  int ned = 
    Ng_GetElement_Edges (elnr+1, &ednums[0], 0);
  ednums.SetSize (ned);
  for (int i = 0; i < ned; i++)
    ednums[i]--;
}

  
void MeshAccess :: 
GetSElEdges (int selnr, ARRAY<int> & ednums) const
{
  ednums.SetSize (4);
  int ned = 
    Ng_GetSurfaceElement_Edges (selnr+1, &ednums[0], 0);
  ednums.SetSize (ned);
  for (int i = 0; i < ned; i++)
    ednums[i]--;
}

void MeshAccess :: 
GetSElEdges (int selnr, ARRAY<int> & ednums, ARRAY<int> & orient) const
{
  ednums.SetSize (4);
  orient.SetSize (4);
  int ned = 
    Ng_GetSurfaceElement_Edges (selnr+1, &ednums[0], &orient[0]);
  ednums.SetSize (ned);
  orient.SetSize (ned);

  for (int i = 0; i < ned; i++)
    ednums[i]--;
}



void MeshAccess :: GetEdgeElements (int enr, ARRAY<int> & elnums) const
{
  // old routine: Slow !!!
  /* 
  int nel = GetNE();
  int edges[8];
  elnums.SetSize (0);
  for (int i = 0; i < nel; i++)
    {
      int ned = Ng_GetElement_Edges (i+1, edges, 0);
      for (int j = 0; j < ned; j++)
	if (edges[j]-1 == enr)
	  elnums.Append (i);
    }
  */ 
  
  // fast (he)
  elnums.SetSize(0);
  ARRAY<int> pnums, velems0, velems1; 
  GetEdgePNums(enr, pnums);
  GetVertexElements(pnums[0], velems0);
  GetVertexElements(pnums[1], velems1);
  
  //   now compare
  for (int i=0; i<velems0.Size(); i++) {
    for (int j=0; j<velems1.Size(); j++) {
      if (velems0[i] == velems1[j]) {
        //  hope that this is a fine edge !! 
        // (run several tests and compared with old routine. seems ok)
        elnums.Append(velems0[i]);
        continue;
      }
    }
  }
}
  
  
void MeshAccess :: 
GetElFaces (int elnr, ARRAY<int> & fnums) const
{
  fnums.SetSize (6);
  int nfa = 
    Ng_GetElement_Faces (elnr+1, &fnums[0], 0);
  fnums.SetSize (nfa);

  for (int i = 0; i < nfa; i++)
    fnums[i]--;
}

void MeshAccess :: 
GetElFaces (int elnr, ARRAY<int> & fnums, ARRAY<int> & orient) const
{
  fnums.SetSize (6);
  orient.SetSize (6);
  int nfa = 
    Ng_GetElement_Faces (elnr+1, &fnums[0], &orient[0]);
  fnums.SetSize (nfa);
  orient.SetSize (nfa);

  for (int i = 0; i < nfa; i++)
    fnums[i]--;
}


int MeshAccess :: 
GetSElFace (int selnr) const
{
  return Ng_GetSurfaceElement_Face (selnr+1, 0)-1;
}
  
void MeshAccess :: 
GetSElFace (int selnr, int & fnum, int & orient) const
{
  fnum = Ng_GetSurfaceElement_Face (selnr+1, &orient);
  fnum--;
}

// string MeshAccess :: GetSElBCName (const int elnr) const
//     {
// 	char auxstr[256];
// 	string retval;
// 	Ng_GetSurfaceElementBCName (elnr+1,auxstr);
// 	retval = auxstr;
// 	return retval;
//     }

// string MeshAccess :: GetBCNumBCName (const int bcnr) const
//     {
// 	char auxstr[256];
// 	string retval;
// 	Ng_GetBCNumBCName(bcnr,auxstr);
// 	retval = auxstr;
// 	return retval;
//     }
  
void MeshAccess :: GetFacePNums (int fnr, ARRAY<int> & pnums) const
{
  pnums.SetSize(4);
  int nv = Ng_GetFace_Vertices (fnr+1, &pnums[0]);
  pnums.SetSize(nv);
  for (int i = 0; i < nv; i++) 
    pnums[i]--;
}

 
void MeshAccess :: GetFaceEdges (int fnr, ARRAY<int> & edges) const
{
  edges.SetSize(4);
  int ned = Ng_GetFace_Edges (fnr+1, &edges[0]);
  edges.SetSize(ned);
  for (int i = 0; i < ned; i++) 
    edges[i]--;
}
 

void MeshAccess :: GetFaceElements (int fnr, ARRAY<int> & elnums) const
{
  int nel = GetNE();
  int faces[8];
  elnums.SetSize (0);
  for (int i = 0; i < nel; i++)
    {
      int nfa = Ng_GetElement_Faces (i+1, faces, 0);
      for (int j = 0; j < nfa; j++)
	if (faces[j]-1 == fnr)
	  elnums.Append (i);
    }
}


void MeshAccess :: GetEdgePNums (int fnr, int & pn1, int & pn2) const
{
  int v2[2];
  Ng_GetEdge_Vertices (fnr+1, v2);
  pn1 = v2[0]-1;
  pn2 = v2[1]-1;
}

void MeshAccess :: GetEdgePNums (int fnr, ARRAY<int> & pnums) const
{
  pnums.SetSize(2);
  Ng_GetEdge_Vertices (fnr+1, &pnums[0]);
  pnums[0] -= 1;
  pnums[1] -= 1;
}

  
void MeshAccess :: GetSElPNums (int selnr, ARRAY<int> & pnums) const
{
  pnums.SetSize (NG_SURFACE_ELEMENT_MAXPOINTS);
  NG_ELEMENT_TYPE typ = Ng_GetSurfaceElement (selnr+1, &pnums[0]);
  switch (typ)
    {
    case NG_SEGM: pnums.SetSize(2); break;
    case NG_SEGM3: pnums.SetSize(3); break;

    case NG_TRIG: pnums.SetSize(3); break;
    case NG_TRIG6: pnums.SetSize(6); break;
    case NG_QUAD: pnums.SetSize(4); break;
    case NG_QUAD6: pnums.SetSize(6); break;

    default:
      {
	cerr << "MeshAccess::GetSElPNums undefined type "
	     << int(Ng_GetSurfaceElement (selnr, &pnums[0])) << endl;
      }

    }
  for (int i = 0; i < pnums.Size(); i++) pnums[i]--;
}


void MeshAccess :: GetSElVertices (int selnr, ARRAY<int> & vnums) const
{
  vnums.SetSize (NG_SURFACE_ELEMENT_MAXPOINTS);
  NG_ELEMENT_TYPE typ = Ng_GetSurfaceElement (selnr+1, &vnums[0]);
  switch (typ)
    {
    case NG_SEGM: vnums.SetSize(2); break;
    case NG_SEGM3: vnums.SetSize(2); break;

    case NG_TRIG: vnums.SetSize(3); break;
    case NG_TRIG6: vnums.SetSize(3); break;
    case NG_QUAD: vnums.SetSize(4); break;
    case NG_QUAD6: vnums.SetSize(4); break;

    default:
      {
	cerr << "MeshAccess::GetSElPNums undefined type "
	     << int(Ng_GetSurfaceElement (selnr+1, &vnums[0])) << endl;
      }

    }
  for (int i = 0; i < vnums.Size(); i++) vnums[i]--;
}



// he: some utility for Facets
void MeshAccess :: GetElFacets (int elnr, ARRAY<int> & fnums) const
{
  if (GetDimension() == 2) GetElEdges(elnr, fnums);
  else GetElFaces(elnr, fnums);
} 

void MeshAccess :: GetElFacets (int elnr, ARRAY<int> & fnums, ARRAY<int> & orient) const
{
  if (GetDimension() == 2) GetElEdges(elnr, fnums, orient);
  else GetElFaces(elnr, fnums, orient);
} 
      
void MeshAccess :: GetSElFacets (int selnr, ARRAY<int> & fnums) const
{
  if (GetDimension() == 2) 
    GetSElEdges(selnr, fnums);
  else
  {
    fnums.SetSize(1);
    fnums[0] = GetSElFace(selnr);
  }
}
  
void MeshAccess :: GetSElFacet (int selnr, ARRAY<int> & fnums, ARRAY<int> & orient) const
{
  if (GetDimension() == 2) 
    GetSElEdges(selnr, fnums, orient);
  else
  {
    fnums.SetSize(1);
    orient.SetSize(1);
    GetSElFace(selnr, fnums[0], orient[0]);
  }
}
  
void MeshAccess :: GetFacetPNums (int fnr, ARRAY<int> & pnums) const
{
  if (GetDimension() == 2)
    GetEdgePNums(fnr, pnums);
  else
    GetFacePNums(fnr, pnums);
}







void MeshAccess::PushStatus (const char * str) const
{ 
  Ng_PushStatus (str);
}
void MeshAccess::PopStatus () const
{ 
  Ng_PopStatus ();
}
void MeshAccess::SetThreadPercentage (double percent) const
{ 
  Ng_SetThreadPercentage (percent); 
}
void MeshAccess :: GetStatus (string & str, double & percent) const
{
  char * s;
  Ng_GetStatus(&s,percent);
  str = s;
}



void MeshAccess :: SetTerminate(void) const
{
  Ng_SetTerminate();
}
void MeshAccess :: UnSetTerminate(void) const
{
  Ng_UnSetTerminate();
}
bool MeshAccess :: ShouldTerminate(void) const
{
  return (Ng_ShouldTerminate() != 0);
}


void MeshAccess :: GetVertexElements (int vnr, ARRAY<int> & elnrs) const
{
  int nel = Ng_GetNVertexElements (vnr+1);
  elnrs.SetSize (nel);
  Ng_GetVertexElements (vnr+1, &elnrs[0]);
  for (int j = 0; j < nel; j++)
    elnrs[j]--;
}






  


#ifdef NETGEN_ELTRANS

void MeshAccess ::
GetElementTransformation (int elnr, ElementTransformation & eltrans,
			  LocalHeap & lh) const
{
  int elind = Ng_GetElementIndex (elnr+1)-1;
  eltrans.SetElement (0, elnr, elind);
  eltrans.SetElementType (GetElType(elnr));
  if (pts.Size())
    eltrans.SetGeometryData (&pts, &dxdxis, &first_of_element);

  if(higher_integration_order.Size() == GetNE() && higher_integration_order[elnr])
    eltrans.SetHigherIntegrationOrder();
  else
    eltrans.UnSetHigherIntegrationOrder();
}

#else

void MeshAccess :: 
GetElementTransformation (int elnr, ElementTransformation & eltrans,
			  LocalHeap & lh) const
{
  int elind = Ng_GetElementIndex (elnr+1)-1;
  int pnums[NG_ELEMENT_MAXPOINTS];
  int np;

  NG_ELEMENT_TYPE et = Ng_GetElement (elnr+1, pnums);

  switch (et)
    {
    case NG_TRIG:
        eltrans.SetElement (&trig1, elnr, elind); break;
    case NG_TRIG6:
      eltrans.SetElement (&trig2, elnr, elind); break;
    case NG_QUAD:
      eltrans.SetElement (&quad1, elnr, elind); break;
    case NG_QUAD6:
      eltrans.SetElement (&quad2, elnr, elind); break;

    case NG_TET:
      eltrans.SetElement (&tet1, elnr, elind); break;
    case NG_TET10:
      eltrans.SetElement (&tet2, elnr, elind); break;
    case NG_PRISM:
      eltrans.SetElement (&prism1, elnr, elind); break;
    case NG_PRISM12:
      eltrans.SetElement (&prism2, elnr, elind); break;
    case NG_PYRAMID:
      eltrans.SetElement (&pyramid1, elnr, elind); break;
    case NG_HEX:
      eltrans.SetElement (&hex1, elnr, elind); break;

    default:
      cerr << "element transformation for element " << int(et)
	   << " not defined" << endl;
    }

  np = eltrans.GetElement().GetNDof();

  double point[3];
  int i, j;

  int dim = GetDimension();
  
  FlatMatrix<> & pmat = const_cast<FlatMatrix<>&> (eltrans.PointMatrix());
  pmat.AssignMemory (dim, np, lh);

  for (i = 0; i < np; i++)
    {
      Ng_GetPoint (pnums[i], point);
      for (j = 0; j < dim; j++)
	pmat(j, i) = point[j];
    }


  if(higher_integration_order.Size() == GetNE() && higher_integration_order[elnr])
    eltrans.SetHigherIntegrationOrder();
  else
    eltrans.UnSetHigherIntegrationOrder();
}

#endif



#ifdef NETGEN_ELTRANS

void MeshAccess ::
GetSurfaceElementTransformation (int elnr, ElementTransformation & eltrans,
				 LocalHeap & lh) const
{
  int elind = Ng_GetSurfaceElementIndex (elnr+1)-1;
  eltrans.SetElement (1, elnr, elind);
  eltrans.SetElementType (GetSElType(elnr));
}    

#else

void MeshAccess :: 
GetSurfaceElementTransformation (int elnr, ElementTransformation & eltrans, 
				 LocalHeap & lh) const
{
  try
    {
      int elind = Ng_GetSurfaceElementIndex (elnr+1)-1;
      
      int np;
      int pnums[NG_SURFACE_ELEMENT_MAXPOINTS];
      
      NG_ELEMENT_TYPE et = Ng_GetSurfaceElement (elnr+1, pnums);
      switch (et)
	{
	case NG_SEGM:
	  eltrans.SetElement (&segm1, elnr, elind); break;
	case NG_SEGM3:
	  eltrans.SetElement (&segm2, elnr, elind); break;
	  
	case NG_TRIG:
	  eltrans.SetElement (&trig1, elnr, elind); break;
	case NG_TRIG6:
	  eltrans.SetElement (&trig2, elnr, elind); break;
	case NG_QUAD:
	  eltrans.SetElement (&quad1, elnr, elind); break;
	case NG_QUAD6:
	  eltrans.SetElement (&quad2, elnr, elind); break;
	  
	default:
	  cerr << "surface element transformation for element " 
	       << int(et)
	       << " not defined" << endl;
	}

      np = eltrans.GetElement().GetNDof();

	
      FlatMatrix<> & pmat = const_cast<FlatMatrix<>&> (eltrans.PointMatrix());
      FlatMatrix<> & nvmat = const_cast<FlatMatrix<>&> (eltrans.NVMatrix());
      
      int dim = GetDimension();
      
      pmat.AssignMemory (dim, np, lh);
      nvmat.AssignMemory (dim, np, lh);
      
      double point[3], nv[3];
      int i, j;
      
      for (i = 0; i < np; i++)
	{
	  Ng_GetPoint (pnums[i], point);
	  for (j = 0; j < dim; j++)
	    pmat(j, i) = point[j];
	  
	  Ng_GetNormalVector (elnr+1, i+1, nv);
	  for (j = 0; j < dim; j++)
	    nvmat(j, i) = nv[j];
	}
      
      if (dim == 2)
	{
	  double tx = pmat(0,1)-pmat(0,0);
	  double ty = pmat(1,1)-pmat(1,0);
	  double len = sqrt (tx*tx+ty*ty);
	  if (len != 0)
	    { tx /= len; ty /= len; }
	  for (i = 0; i < np; i++)
	    {
	      nvmat(0,i) = ty;
	      nvmat(1,i) = -tx;
	    }
	}
    }
  catch (Exception & e)
    {
      stringstream ost;
      ost << "in GetSurfaceElementTransformation" << endl;
      e.Append (ost.str());
      throw;
    }
  catch (exception & e)
    {
      throw (Exception (string(e.what()) +
			string("\n in GetSurfaceElementTransformation\n")));
    }

}

#endif



  
double MeshAccess :: ElementVolume (int elnr) const
{

  const FiniteElement * fe;
  switch (GetElType (elnr))
    {
    case ET_TRIG: fe = &trig0; break;
    case ET_QUAD: fe = &quad0; break;
    case ET_TET: fe = &tet0; break;
    case ET_PYRAMID: fe = &pyramid0; break;
    case ET_PRISM: fe = &prism0; break;
    case ET_HEX: fe = &hex0; break;
    default:
      {
	cerr << "ElementVolume not implemented for el " << GetElType(elnr) << endl;
      }
    }
  
  char d[10000];
  LocalHeap lh(d, 10000);

  ElementTransformation trans;
  GetElementTransformation (elnr, trans, lh);
  ConstantCoefficientFunction ccf(1);


  if (GetDimension() == 2)
    {
      SourceIntegrator<2> si( &ccf );
      FlatVector<> elvec;
      si.AssembleElementVector (*fe, trans, elvec, lh);
      return elvec(0);
    }
  else
    {
      SourceIntegrator<3> si( &ccf );
      FlatVector<> elvec;
      si.AssembleElementVector (*fe, trans, elvec, lh);
      return elvec(0);
    }
}


  
double MeshAccess :: SurfaceElementVolume (int selnr) const
{
  const FiniteElement * fe;
  switch (GetSElType (selnr))
    {
    case ET_TRIG: fe = &trig0; break;
    case ET_QUAD: fe = &quad0; break;
    default:
      {
	cerr << "SurfaceElementVolume not implemented for el " << GetElType(selnr) << endl;
	return 0;
      }
    }

  char d[10000];
  LocalHeap lh(d, 10000);

  ElementTransformation trans;
  GetSurfaceElementTransformation (selnr, trans, lh);
  ConstantCoefficientFunction ccf(1);

  if (GetDimension() == 2)
    {
      NeumannIntegrator<2> si( &ccf );
      FlatVector<> elvec;
      si.AssembleElementVector (*fe, trans, elvec, lh);
      return elvec(0);
    }
  else
    {
      NeumannIntegrator<3> si( &ccf );
      FlatVector<> elvec;
      si.AssembleElementVector (*fe, trans, elvec, lh);
      return elvec(0);
    }
}

/*
int MeshAccess :: GetNLevels() const
{
  return Ng_GetNLevels();
}
*/



void MeshAccess :: SetPointSearchStartElement(const int el) const
{
  Ng_SetPointSearchStartElement(el+1);
}

int MeshAccess :: FindElementOfPoint (FlatVector<double> point,
				      IntegrationPoint & ip, 
				      bool build_searchtree,
				      const int index) const
{
  ARRAY<int> dummy(1);
  dummy[0] = index;
  return FindElementOfPoint(point,ip,build_searchtree,&dummy);
}

int MeshAccess :: FindElementOfPoint (FlatVector<double> point,
				      IntegrationPoint & ip,
				      bool build_searchtree,
				      const ARRAY<int> * const indices) const
{
  int elnr;
  double lami[3];
  if(indices != NULL && indices->Size()>0)
    {
      int * dummy = new int[indices->Size()];
      for(int i=0; i<indices->Size(); i++) dummy[i] = (*indices)[i]+1;

      elnr = Ng_FindElementOfPoint (&point(0), lami, build_searchtree,dummy,indices->Size());

      delete [] dummy;
    }
  else
    {  
      elnr = Ng_FindElementOfPoint (&point(0), lami, build_searchtree);
    }

  if (elnr == 0) return -1;

  for (int k = 0; k < 3; k++)
    ip(k) = lami[k];

  return elnr-1;
}


int MeshAccess :: FindSurfaceElementOfPoint (FlatVector<double> point,
					     IntegrationPoint & ip, 
					     bool build_searchtree,
					     const int index) const
{
  ARRAY<int> dummy(1);
  dummy[0] = index;
  return FindSurfaceElementOfPoint(point,ip,build_searchtree,&dummy);
}

int MeshAccess :: FindSurfaceElementOfPoint (FlatVector<double> point,
					     IntegrationPoint & ip,
					     bool build_searchtree,
					     const ARRAY<int> * const indices) const
{
  int elnr;
  double lami[3];
  if(indices != NULL && indices->Size()>0)
    {
      int * dummy = new int[indices->Size()];
      for(int i=0; i<indices->Size(); i++) dummy[i] = (*indices)[i]+1;

      elnr = Ng_FindSurfaceElementOfPoint (&point(0), lami, build_searchtree,dummy,indices->Size());

      delete [] dummy;
    }
  else
    {  
      elnr = Ng_FindSurfaceElementOfPoint (&point(0), lami, build_searchtree);
    }

  if (elnr == 0) return -1;

  for (int k = 0; k < 3; k++)
    ip(k) = lami[k];


  return elnr-1;
}


int MeshAccess :: GetNPairsPeriodicVertices () const 
{
  return Ng_GetNPeriodicVertices(0);
}

int MeshAccess :: GetNPairsPeriodicVertices (int idnr) const 
{
  return Ng_GetNPeriodicVertices(idnr);
}
 
void MeshAccess :: GetPeriodicVertices ( ARRAY<INT<2> > & pairs) const
{
  int npairs;

  npairs = Ng_GetNPeriodicVertices (0);
  pairs.SetSize (npairs);

  Ng_GetPeriodicVertices (0,&pairs[0][0]);
  for (int i = 0; i < pairs.Size(); i++)
    {
      pairs[i][0]--;
      pairs[i][1]--;
    }
}
 
void MeshAccess :: GetPeriodicVertices (int idnr, ARRAY<INT<2> > & pairs) const
{
  int npairs;

  npairs = Ng_GetNPeriodicVertices (idnr);

  pairs.SetSize (npairs);

  Ng_GetPeriodicVertices (idnr,&pairs[0][0]);

  for (int i = 0; i < pairs.Size(); i++)
    {
      pairs[i][0]--;
      pairs[i][1]--;
    }
}




int MeshAccess :: GetNPairsPeriodicEdges () const 
{
  return Ng_GetNPeriodicEdges(0);
}
 
void MeshAccess :: GetPeriodicEdges ( ARRAY<INT<2> > & pairs) const
{
  int npairs;

  npairs = Ng_GetNPeriodicEdges (0);
  pairs.SetSize (npairs);

  Ng_GetPeriodicEdges (0,&pairs[0][0]);
  for (int i = 0; i < pairs.Size(); i++)
    {
      pairs[i][0]--;
      pairs[i][1]--;
    }
}

int MeshAccess :: GetNPairsPeriodicEdges (int idnr) const 
{
  return Ng_GetNPeriodicEdges(idnr);
}
 
void MeshAccess :: GetPeriodicEdges (int idnr, ARRAY<INT<2> > & pairs) const
{
  int npairs;

  npairs = Ng_GetNPeriodicEdges (idnr);
  pairs.SetSize (npairs);

  Ng_GetPeriodicEdges (idnr,&pairs[0][0]);
  for (int i = 0; i < pairs.Size(); i++)
    {
      pairs[i][0]--;
      pairs[i][1]--;
    }
}

/*
///// Added by Roman Stainko ....
void MeshAccess::GetVertexElements( int vnr, ARRAY<int>& elems) const
{
  int nel = Ng_GetVertex_NElements( vnr+1 );
  elems.SetSize( nel );
  
  Ng_GetVertex_Elements( vnr+1, &elems[0] );
  
  for( int i=0; i<nel; i++ )
    elems[i]--;
 
}
*/

///// Added by Roman Stainko ....
void MeshAccess::GetVertexSurfaceElements( int vnr, ARRAY<int>& elems) const
{
  int nel = Ng_GetVertex_NSurfaceElements( vnr+1 );
  elems.SetSize( nel );
  
  Ng_GetVertex_SurfaceElements( vnr+1, &elems[0] );

  for( int i=0; i<nel; i++ )
    elems[i]--;
}



void MeshAccess::SetHigherIntegrationOrder(int elnr)
{
  if(higher_integration_order.Size() != GetNE())
    {
      higher_integration_order.SetSize(GetNE());
      higher_integration_order = false;
    }
  higher_integration_order[elnr] = true;
}
void MeshAccess::UnSetHigherIntegrationOrder(int elnr)
{
  if(higher_integration_order.Size() != GetNE())
    {
      higher_integration_order.SetSize(GetNE());
      higher_integration_order = false;
    }
  higher_integration_order[elnr] = false;
}



void MeshAccess :: PrecomputeGeometryData(int intorder)
{
  return;

  cout << "Precompute Geometry Data, intorder = " << intorder << " ..." << flush;

  int npts = 0;
  ARRAY<int> vnums;
  for (int i = 0; i < GetNE(); i++)
    {
      const IntegrationRule & ir1d =
	GetIntegrationRules().SelectIntegrationRule (ET_SEGM, intorder);
      GetElVertices (i, vnums);
      
      /*
      TensorProductIntegrationRule ir (GetElType (i), &vnums[0], ir1d);
      npts += ir.GetNIP();
      */
    }

  cout << ", npts = " << npts << ", memory = " << npts * 8 * 12 << " ... " << flush;

  pts.SetSize(npts);
  dxdxis.SetSize (npts);
  first_of_element.SetSize (GetNE());

  npts = 0;
  for (int i = 0; i < GetNE(); i++)
    {
      const IntegrationRule & ir1d =
	GetIntegrationRules().SelectIntegrationRule (ET_SEGM, intorder);
      GetElVertices (i, vnums);

      /*
      TensorProductIntegrationRule ir (GetElType (i), &vnums[0], ir1d);

      first_of_element[i] = npts;

      for (int j = 0; j < ir.GetNIP(); j++, npts++)
	Ng_GetElementTransformation (i+1, &ir[j](0), &pts[npts](0), &dxdxis[npts](0,0));
      */
    }

  // (*testout) << "computed pts = " << pts << endl;
  cout << " done" << endl;
}



void MeshAccess :: InitPointCurve(double red, double green, double blue) const
{
  Ng_InitPointCurve(red, green, blue);
}

void MeshAccess :: AddPointCurvePoint(const Vec<3> & point) const
{
  Ng_AddPointCurvePoint(&(point(0)));
}

}





