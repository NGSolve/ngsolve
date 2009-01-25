#ifdef PARALLEL


#include <meshing.hpp>
#include "paralleltop.hpp"


namespace netgen
{



  void ParallelMeshTopology :: Reset ()
  {
    *testout << "ParallelMeshTopology::Reset" << endl;
    
    if ( ntasks == 1 ) return;
    PrintMessage ( 4, "RESET");
    
    int nvold = nv;
    int nedold = ned;
    int nfaold = nfa;
    
    ne = mesh.GetNE();
    nv = mesh.GetNV();
    nseg = mesh.GetNSeg();
    nsurfel = mesh.GetNSE();
    
    ned = mesh.GetTopology().GetNEdges();
    nfa = mesh.GetTopology().GetNFaces();
    

    loc2distedge.ChangeSize (ned);
    for (int i = 0; i < ned; i++)
      if (loc2distedge[i].Size() == 0)
        loc2distedge.Add (i, -1);  // will be the global nr

    loc2distface.ChangeSize (nfa);
    for (int i = 0; i < nfa; i++)
      if (loc2distface[i].Size() == 0)
        loc2distface.Add (i, -1);  // will be the global nr


  if ( !isexchangevert ) 
    {
      isexchangevert = new BitArray (nv * ( ntasks+1 ));
      isexchangevert->Clear();
    }
  if ( !isexchangeedge ) 
    {
      isexchangeedge = new BitArray (ned*(ntasks+1) );
      isexchangeedge->Clear();
    }
  if ( !isexchangeface ) 
    {
      isexchangeface = new BitArray (nfa*(ntasks+1) );
      isexchangeface->Clear();
    }
  if ( !isexchangeel ) 
    {
      isexchangeel = new BitArray (ne*(ntasks+1) );
      isexchangeel->Clear();
    }


  // if the number of vertices did not change, return
  if ( nvold == nv ) return;

  // faces and edges get new numbers -> delete 
  isexchangeface -> SetSize(nfa*(ntasks+1) );
  isexchangeedge -> SetSize(ned*(ntasks+1) );
  isexchangeface -> Clear();
  isexchangeedge -> Clear();


  SetNV(nv);
  SetNE(ne);


  if ( !isghostedge.Size() )
    {
      isghostedge.SetSize(ned);
      isghostedge.Clear();
    }
  if ( !isghostface.Size() )
    {
      isghostface.SetSize(nfa);
      isghostface.Clear();
    }

}


 ParallelMeshTopology :: ~ParallelMeshTopology ()
 {
   delete isexchangeface;
   delete isexchangevert;
   delete isexchangeedge;
   delete isexchangeel;
 }



 ParallelMeshTopology :: ParallelMeshTopology ( const netgen::Mesh & amesh )
  : mesh(amesh)
{
  ned = 0; //mesh.GetTopology().GetNEdges();
  nfa = 0; //mesh.GetTopology().GetNFaces();
  nv = 0;
  ne = 0;
  np = 0;
  nseg = 0;
  nsurfel = 0;
  neglob = 0;
  nvglob = 0;

  nparel = 0;

  isexchangeface = 0;
  isexchangevert = 0;
  isexchangeel = 0;
  isexchangeedge = 0;

  coarseupdate = 0;

  isghostedge.SetSize(0);
  isghostface.SetSize(0);

  overlap = 0;
}


  int ParallelMeshTopology :: Glob2Loc_Vert (int globnum )
  {
    for (int i = 1; i <= nv; i++)
      if ( globnum == loc2distvert[i][0] )
	return i;

    return -1;
  }

int ParallelMeshTopology :: Glob2Loc_VolEl (int globnum )
{
  int locnum = -1;
  for (int i = 0; i < ne; i++)
    {
      if ( globnum == loc2distel[i][0] )
	{
	  locnum = i+1;
	}
    }
  return locnum;
}

int ParallelMeshTopology :: Glob2Loc_SurfEl (int globnum )
  {
    int locnum = -1;
    for (int i = 0; i < nsurfel; i++)
      {
	if ( globnum == loc2distsurfel[i][0] )
	  {
	    locnum = i+1;
	  }
      }
    return locnum;
  }

int ParallelMeshTopology :: Glob2Loc_Segm (int globnum )
  {
    int locnum = -1;
    for (int i = 0; i < nseg; i++)
      {
	if ( globnum == loc2distsegm[i][0] )
	  {
	    locnum = i+1;
	  }
      }
    return locnum;
  }

 

  void ParallelMeshTopology ::  Print() const
  {

      (*testout) << endl <<  "TOPOLOGY FOR PARALLEL MESHES" << endl << endl;


      for ( int i = 1; i <= nv; i++ )
       	if ( IsExchangeVert (i) )
	  {
	    (*testout) << "exchange point  " << i << ":  global " << GetLoc2Glob_Vert(i) << endl;
 	    for ( int dest = 0; dest < ntasks; dest ++)
	      if ( dest != id )
 	      if ( GetDistantPNum( dest, i  ) > 0 )
 		(*testout) << "   p" << dest << ": " << GetDistantPNum ( dest, i ) << endl; 
	  }

      for ( int i = 1; i <= ned; i++ )
	if ( IsExchangeEdge ( i ) )
	  {
	    int v1, v2;
	    mesh . GetTopology().GetEdgeVertices(i, v1, v2);
	    (*testout) << "exchange edge  " << i << ":  global vertices "  << GetLoc2Glob_Vert(v1) << "  " 
		       << GetLoc2Glob_Vert(v2) << endl;
	    for ( int dest = 0; dest < ntasks; dest++)
	      if ( GetDistantEdgeNum ( dest, i ) > 0 )
		if ( dest != id )
		{
		  (*testout) << "   p" << dest << ": " << GetDistantEdgeNum ( dest, i ) << endl;
		}
	  }


      for ( int i = 1; i <= nfa; i++ )
	if ( IsExchangeFace(i) )
	  {
	    Array<int> facevert;
	    mesh . GetTopology().GetFaceVertices(i, facevert);
	    
	    (*testout) << "exchange face  " << i << ":  global vertices " ;
	    for ( int fi=0; fi < facevert.Size(); fi++)
	      (*testout) << GetLoc2Glob_Vert(facevert[fi]) << "  ";
	    (*testout) << endl; 
 	    for ( int dest = 0; dest < ntasks; dest++)
	      if ( dest != id )
	      {
		if ( GetDistantFaceNum ( dest, i ) >= 0 )
 		(*testout) << "   p" << dest << ": " << GetDistantFaceNum ( dest, i ) << endl;
	      }
	  }


      for ( int i = 1; i < mesh.GetNE(); i++)
	{
	  if ( !IsExchangeElement(i) ) continue;
	  Array<int> vert;
	  const Element & el = mesh.VolumeElement(i);

	  (*testout) << "parallel local element " << i << endl;
	  
	  (*testout) << "vertices " ;
	  for ( int j = 0; j < el.GetNV(); j++)
	    (*testout) << el.PNum(j+1)  << "  ";
	  (*testout) << "is ghost " << IsGhostEl(i) << endl;
	  (*testout) << endl;


	}
  }





int  ParallelMeshTopology :: GetDistantPNum ( int proc, int locpnum ) const
  {
    if ( proc == 0 )
      return loc2distvert[locpnum][0];

    for (int i = 1; i < loc2distvert[locpnum].Size(); i += 2)
      if ( loc2distvert[locpnum][i] == proc )
	return loc2distvert[locpnum][i+1];

    return -1;
  } 

int  ParallelMeshTopology :: GetDistantFaceNum ( int proc, int locfacenum ) const
  {
    if ( proc == 0 )
      return loc2distface[locfacenum-1][0];

    for ( int i = 1; i < loc2distface[locfacenum-1].Size(); i+=2 )
      if ( loc2distface[locfacenum-1][i] == proc )
	return loc2distface[locfacenum-1][i+1];

    return -1;
  } 

int  ParallelMeshTopology :: GetDistantEdgeNum ( int proc, int locedgenum ) const
  {
    if ( proc == 0 )
      return loc2distedge[locedgenum-1][0];

    for ( int i = 1; i < loc2distedge[locedgenum-1].Size(); i+=2 )
      if ( loc2distedge[locedgenum-1][i] == proc )
	return loc2distedge[locedgenum-1][i+1];

    return -1;
  } 

int  ParallelMeshTopology :: GetDistantElNum ( int proc, int locelnum ) const
  {
    if ( proc == 0 )
      return loc2distel[locelnum-1][0];

    for ( int i = 1; i < loc2distel[locelnum-1].Size(); i+=2 )
      if ( loc2distel[locelnum-1][i] == proc )
	return loc2distel[locelnum-1][i+1];

    return -1;
  } 



  /*
//
// gibt anzahl an distant pnums zurueck
int  ParallelMeshTopology :: GetNDistantPNums ( int locpnum ) const
{
  return loc2distvert[locpnum].Size() / 2 + 1;
} 

int  ParallelMeshTopology :: GetNDistantFaceNums ( int locfacenum ) const
  {
    int size = loc2distface[locfacenum-1].Size() / 2 + 1;
    return size;
  } 

int  ParallelMeshTopology :: GetNDistantEdgeNums ( int locedgenum ) const
  {
    int size = loc2distedge[locedgenum-1].Size() / 2 + 1;
    return size;
  } 

int  ParallelMeshTopology :: GetNDistantElNums ( int locelnum ) const
  {
    int size = loc2distel[locelnum-1].Size() / 2 + 1;
    return size;
  } 
*/


  // gibt anzahl an distant pnums zurueck
  // * pnums entspricht Array<int[2] >
int  ParallelMeshTopology :: GetDistantPNums ( int locpnum, int * distpnums ) const
  {
//     distpnums[0] = loc2distvert[locpnum][0];

//     for (int i = 1; i < loc2distvert[locpnum].Size(); i += 2)
//        distpnums[ loc2distvert[locpnum][i] ] = loc2distvert[locpnum][i+1];
    distpnums[0] = 0;
    distpnums[1] = loc2distvert[locpnum][0];
    for ( int i = 1; i < loc2distvert[locpnum].Size(); i++ )
      distpnums[i+1] = loc2distvert[locpnum][i];

     int size = loc2distvert[locpnum].Size() / 2 + 1;
     return size;
  } 

int  ParallelMeshTopology :: GetDistantFaceNums ( int locfacenum, int * distfacenums ) const
  {
//     distfacenums[0] = loc2distface[locfacenum-1][0];

//     for ( int i = 1; i < loc2distface[locfacenum-1].Size(); i+=2 )
//       distfacenums[loc2distface[locfacenum-1][i]] = loc2distface[locfacenum-1][i+1];

    distfacenums[0] = 0;
    distfacenums[1] = loc2distface[locfacenum-1][0];

    for ( int i = 1; i < loc2distface[locfacenum-1].Size(); i++ )
      distfacenums[i+1] = loc2distface[locfacenum-1][i];

    int size = loc2distface[locfacenum-1].Size() / 2 + 1;
    return size;
  } 

int  ParallelMeshTopology :: GetDistantEdgeNums ( int locedgenum, int * distedgenums ) const
  {
//     distedgenums[0] = loc2distedge[locedgenum-1][0];

//     for ( int i = 1; i < loc2distedge[locedgenum-1].Size(); i+=2 )
//       distedgenums[loc2distedge[locedgenum-1][i]] = loc2distedge[locedgenum-1][i+1];

    distedgenums[0] = 0;
    distedgenums[1] = loc2distedge[locedgenum-1][0];

    for ( int i = 1; i < loc2distedge[locedgenum-1].Size(); i++ )
      distedgenums[i+1] = loc2distedge[locedgenum-1][i];

    int size = loc2distedge[locedgenum-1].Size() / 2 + 1;
    return size;
  } 

int  ParallelMeshTopology :: GetDistantElNums ( int locelnum, int * distelnums ) const
  {
//     distelnums[0] = loc2distel[locelnum-1][0];

//     for ( int i = 1; i < loc2distel[locelnum-1].Size(); i+=2 )
//       distelnums[loc2distel[locelnum-1][i]] = loc2distel[locelnum-1][i+1];

    distelnums[0] = 0;
    distelnums[1] = loc2distel[locelnum-1][0];

    for ( int i = 1; i < loc2distel[locelnum-1].Size(); i++ )
      distelnums[i+1] = loc2distel[locelnum-1][i];

    int size = loc2distel[locelnum-1].Size() / 2 + 1;
    return size;
  } 





  void ParallelMeshTopology :: SetDistantFaceNum ( int dest, int locnum, int distnum )
  {
    if ( dest == 0 )
      {
	loc2distface[locnum-1][0] = distnum;
	return;
      }

    for ( int i = 1; i < loc2distface[locnum-1].Size(); i+=2 )
      if ( loc2distface[locnum-1][i] == dest )
	{
	  loc2distface[locnum-1][i+1] = distnum;
	  return;
	}

    loc2distface.Add(locnum-1, dest);
    loc2distface.Add(locnum-1, distnum);
  }

  void ParallelMeshTopology :: SetDistantPNum ( int dest, int locnum, int distnum )
  {
    if ( dest == 0 )
      {
	loc2distvert[locnum][0] = distnum;  // HERE
	return;
      }

    for ( int i = 1;  i < loc2distvert[locnum].Size(); i+=2 )
      if ( loc2distvert[locnum][i] == dest )
	{
	  loc2distvert[locnum][i+1] = distnum;
	  return;
	}

    loc2distvert.Add (locnum, dest);  
    loc2distvert.Add (locnum, distnum); 
  }


  void ParallelMeshTopology :: SetDistantEdgeNum ( int dest, int locnum, int distnum )
  {
    if ( dest == 0 )
      {
	loc2distedge[locnum-1][0] = distnum;
	return;
      }

    for ( int i = 1; i < loc2distedge[locnum-1].Size(); i+=2 )
      if ( loc2distedge[locnum-1][i] == dest )
	{
	  loc2distedge[locnum-1][i+1] = distnum;
	  return;
	}

    loc2distedge.Add (locnum-1, dest);
    loc2distedge.Add (locnum-1, distnum);
  }

  void ParallelMeshTopology :: SetDistantEl ( int dest, int locnum, int distnum )
  {
    if ( dest == 0 )
      {
	loc2distel[locnum-1][0] = distnum;
	return;
      }

    for ( int i = 1; i < loc2distel[locnum-1].Size(); i+=2 )
      if ( loc2distel[locnum-1][i] == dest )
	{
	  loc2distel[locnum-1][i+1] = distnum;
	  return;
	}


    loc2distel.Add (locnum-1, dest);
    loc2distel.Add (locnum-1, distnum);
  }

  void ParallelMeshTopology :: SetDistantSurfEl ( int dest, int locnum, int distnum )
  {
     if ( dest == 0 )
      {
	loc2distsurfel[locnum-1][0] = distnum;
	return;
      }

     for ( int i = 1;  i < loc2distsurfel[locnum-1].Size(); i+=2 )
      if ( loc2distsurfel[locnum-1][i] == dest )
	{
	  loc2distsurfel[locnum-1][i+1] = distnum;
	  return;
	}

    loc2distsurfel.Add (locnum-1, dest);
    loc2distsurfel.Add (locnum-1, distnum);
  }

  void ParallelMeshTopology :: SetDistantSegm ( int dest, int locnum, int distnum )
  {
    if ( dest == 0 )
      {
	loc2distsegm[locnum-1][0] = distnum;
	return;
      }

    for (int i = 1; i < loc2distsegm[locnum-1].Size(); i+=2 )
      if ( loc2distsegm[locnum-1][i] == dest )
	{
	  loc2distsegm[locnum-1][i+1] = distnum;
	  return;
	}

    loc2distsegm.Add (locnum-1, dest);
    loc2distsegm.Add (locnum-1, distnum);
  }

  void ParallelMeshTopology :: GetVertNeighbours ( int vnum, Array<int> & dests ) const
  {
    dests.SetSize(0);
    int i = 1;
    while ( i < loc2distvert[vnum].Size() )
      {
	dests.Append ( loc2distvert[vnum][i] );
	i+=2;
      }
  }


  void ParallelMeshTopology :: Update ()
  {
    ne = mesh.GetNE();
    nv = mesh.GetNV();
    nseg = mesh.GetNSeg();
    nsurfel = mesh.GetNSE();
    
    ned = mesh.GetTopology().GetNEdges();
    nfa = mesh.GetTopology().GetNFaces();
  }








void ParallelMeshTopology :: UpdateRefinement ()
{
  ;
}




void ParallelMeshTopology :: UpdateCoarseGridGlobal ()
{
  PrintMessage ( 1, "UPDATE GLOBAL COARSEGRID STARTS" );      // JS

  // MPI_Barrier (MPI_COMM_WORLD);
  //   PrintMessage ( 1, "all friends are here " );      // JS
  //  MPI_Barrier (MPI_COMM_WORLD);



  int timer = NgProfiler::CreateTimer ("UpdateCoarseGridGlobal");
  NgProfiler::RegionTimer reg(timer);




  *testout << "ParallelMeshTopology :: UpdateCoarseGridGlobal" << endl;
  const MeshTopology & topology = mesh.GetTopology();
  
  Array<int> sendarray, recvarray;

  nfa = topology . GetNFaces();
  ned = topology . GetNEdges();
  np = mesh . GetNP();
  nv = mesh . GetNV();
  ne = mesh . GetNE();
  nseg = mesh.GetNSeg();
  nsurfel = mesh.GetNSE();
  

  // low order processor - save mesh partition
  if ( id == 0 )
    {
      if ( !isexchangeel )
        {
          isexchangeel = new BitArray ( (ntasks+1) * ne );
          isexchangeel -> Clear();
        }


      for ( int eli = 1; eli <= ne; eli++ )
        {
          loc2distel[eli-1][0] = eli;
          SetExchangeElement ( eli );
          const Element & el = mesh . VolumeElement ( eli );
          int dest = el . GetPartition ( );
          SetExchangeElement ( dest, eli );
	  
          for ( int i = 0; i < el.GetNP(); i++ )
            {
              SetExchangeVert ( dest, el.PNum(i+1) );
              SetExchangeVert ( el.PNum(i+1) );
            }
          Array<int> edges;
          topology . GetElementEdges ( eli, edges );
          for ( int i = 0; i < edges.Size(); i++ )
            {
              SetExchangeEdge ( dest, edges[i] );
              SetExchangeEdge ( edges[i] );
            }
          topology . GetElementFaces ( eli, edges );
          for ( int i = 0; i < edges.Size(); i++ )
            {
              SetExchangeFace ( dest, edges[i] );
              SetExchangeFace ( edges[i] );
            }
        }
      
      // HERE
      for ( int i = 1; i <= mesh .GetNV(); i++)
        loc2distvert[i][0] = i;
      
      for ( int i = 0; i < mesh . GetNSeg(); i++)
        loc2distsegm[i][0] = i+1;
      
      for ( int i = 0; i < mesh . GetNSE(); i++)
        loc2distsurfel[i][0] = i+1;
      
      for ( int i = 0; i < topology .GetNEdges(); i++)
        loc2distedge[i][0] = i+1;

      for ( int i = 0; i < topology .GetNFaces(); i++)
        loc2distface[i][0] = i+1;
    }
  

    if ( id == 0 )
      sendarray.Append (nfa);

    BitArray recvface(nfa);
    recvface.Clear();
   
    /*
    Array<int> edges, pnums, faces;
    for ( int el = 1; el <= ne; el++ )
      {
	topology.GetElementFaces (el, faces);
	int globeli = GetLoc2Glob_VolEl(el);

	for ( int fai = 0; fai < faces.Size(); fai++)
	  {
	    int fa = faces[fai];

            topology.GetElementEdges ( el, edges );
	    topology.GetFaceVertices ( fa, pnums );
            
            // send : 
            // localfacenum
            // np
            // ned
            // globalpnums
            // localpnums
            
            // localedgenums mit globalv1, globalv2
                
            sendarray. Append ( fa );
            sendarray. Append ( globeli );
            sendarray. Append ( pnums.Size() );
            sendarray. Append ( edges.Size() );

            if (id == 0)
              for ( int i = 0; i < pnums.Size(); i++ )
                sendarray. Append( pnums[i] );
            else
              for ( int i = 0; i < pnums.Size(); i++ )
                sendarray. Append( GetLoc2Glob_Vert(pnums[i]) );

            for ( int i = 0; i < pnums.Size(); i++ )
              sendarray. Append(pnums[i] );
                
            for ( int i = 0; i < edges.Size(); i++ )
              {
                sendarray. Append(edges[i] );
                int v1, v2;
                topology . GetEdgeVertices ( edges[i], v1, v2 );
                int dv1 = GetLoc2Glob_Vert ( v1 );
                int dv2 = GetLoc2Glob_Vert ( v2 );
                
                if (id > 0) if ( dv1 > dv2 ) swap ( dv1, dv2 );
                sendarray . Append ( dv1 );
                sendarray . Append ( dv2 );
              }
	  }
      }
    */

      // new version
    Array<int> edges, pnums, faces, elpnums;
    sendarray.Append (ne);
    for ( int el = 1; el <= ne; el++ )
      {
	topology.GetElementFaces (el, faces);
	topology.GetElementEdges ( el, edges );
	const Element & volel = mesh.VolumeElement (el);

	int globeli = GetLoc2Glob_VolEl(el);

	sendarray. Append ( globeli );
	sendarray. Append ( faces.Size() );
	sendarray. Append ( edges.Size() );
	sendarray. Append ( volel.GetNP() );

	for ( int i = 0; i < faces.Size(); i++ )
	  sendarray. Append(faces[i] );
	for ( int i = 0; i < edges.Size(); i++ )
	  sendarray. Append(edges[i] );

	for ( int i = 0; i < volel.GetNP(); i++ )
	  if (id == 0)
	    sendarray. Append(volel[i] );
	  else
	    sendarray. Append(GetLoc2Glob_Vert (volel[i]));
      }
    // end new version



    BitArray edgeisinit(ned), vertisinit(np);
    edgeisinit.Clear();
    vertisinit.Clear();
    
    // Array for temporary use, to find local from global element fast
    Array<int,1> glob2loc_el;
    if ( id != 0 ) 
      {
	glob2loc_el.SetSize (neglob);  
	glob2loc_el = -1;
	for ( int locel = 1; locel <= mesh.GetNE(); locel++)
	  glob2loc_el[GetLoc2Glob_VolEl(locel)] = locel;
      }


    // MPI_Barrier (MPI_COMM_WORLD);

    MPI_Request sendrequest;

    if (id == 0)
      {
 	PrintMessage (4, "UpdateCoarseGridGlobal : bcast, size = ", int (sendarray.Size()*sizeof(int)) );
	MyMPI_Bcast ( sendarray );
      }
    else
      MyMPI_ISend ( sendarray, 0, sendrequest );


    int nloops = (id == 0) ? ntasks-1 : 1;
    for (int hi = 0; hi < nloops; hi++)
      {
	int sender;

	if (id == 0)
	  {
	    sender = MyMPI_Recv ( recvarray );
	    PrintMessage (4, "have received from ", sender);
	  }
	else
	  {
	    MyMPI_Bcast ( recvarray );
	    sender = 0;
	  }
	
	// compare received vertices with own ones
	int ii = 0;
	int cntel = 0;
	int volel = 1;

	if ( id != 0 )
	  nfaglob = recvarray[ii++];
	


	Array<int> faces, edges;
	Array<int> pnums, globalpnums;

	int recv_ne = recvarray[ii++];
	for (int hi = 0; hi < recv_ne; hi++)
	  {
	    int globvolel   = recvarray[ii++];
	    int distnfa = recvarray[ii++];
	    int distned = recvarray[ii++];
	    int distnp  = recvarray[ii++];

	    if ( id > 0 )      
	      volel = glob2loc_el[globvolel];
	    else
	      volel = globvolel;
	    
	    if (volel != -1)
	      {
		topology.GetElementFaces( volel, faces);
		topology.GetElementEdges ( volel, edges);
		const Element & volelement = mesh.VolumeElement (volel);

		for ( int i = 0; i  < faces.Size(); i++)
		  SetDistantFaceNum ( sender, faces[i], recvarray[ii++]);
		
		for ( int i = 0; i  < edges.Size(); i++)
		  SetDistantEdgeNum ( sender, edges[i], recvarray[ii++]);

		for ( int i = 0; i  < distnp; i++)
		  SetDistantPNum ( sender, volelement[i], recvarray[ii++]);
	      }
	    else
	      ii += distnfa + distned + distnp;
	  }
      }


    coarseupdate = 1;
    
    if (id != 0)
      {
	MPI_Status status;
	MPI_Wait (&sendrequest, &status);
      }

#ifdef SCALASCA
#pragma pomp inst end(updatecoarsegrid)
#endif
}






  void ParallelMeshTopology :: UpdateCoarseGrid ()
  {
    int timer = NgProfiler::CreateTimer ("UpdateCoarseGrid");
    NgProfiler::RegionTimer reg(timer);

#ifdef SCALASCA
#pragma pomp inst begin(updatecoarsegrid)
#endif
    (*testout) << "UPDATE COARSE GRID PARALLEL TOPOLOGY " << endl;
    PrintMessage ( 1, "UPDATE COARSE GRID PARALLEL TOPOLOGY " );

    // find exchange edges - first send exchangeedges locnum, v1, v2
    // receive distant distnum, v1, v2
    // find matching
    const MeshTopology & topology = mesh.GetTopology();

    UpdateCoarseGridGlobal();


    if ( id == 0 ) return;


    Array<int> * sendarray, *recvarray;
    sendarray = new Array<int> (0);
    recvarray = new Array<int>;

    nfa = topology . GetNFaces();
    ned = topology . GetNEdges();
    np = mesh . GetNP();
    nv = mesh . GetNV();
    ne = mesh . GetNE();
    nseg = mesh.GetNSeg();
    nsurfel = mesh.GetNSE();
    

    sendarray->SetSize (0);

    BitArray recvface(nfa);
    recvface.Clear();
    for ( int fa = 1; fa <= nfa; fa++ )
      {

	if ( !IsExchangeFace ( fa ) ) continue;

	Array<int> edges, pnums;
	int globfa = GetDistantFaceNum ( 0, fa );

	topology.GetFaceEdges ( fa, edges );
	topology.GetFaceVertices ( fa, pnums );
	
	
	// send : 
	// localfacenum
	// globalfacenum
	// np
	// ned
	// globalpnums
	// localpnums
	
	// localedgenums mit globalv1, globalv2
	// 
	
	sendarray -> Append ( fa );
	sendarray -> Append ( globfa );
	sendarray -> Append ( pnums.Size() );
	sendarray -> Append ( edges.Size() );
	for ( int i = 0; i < pnums.Size(); i++ )
	  {
	    sendarray -> Append( GetLoc2Glob_Vert(pnums[i]) );
	  }
	for ( int i = 0; i < pnums.Size(); i++ )
	  {
	    sendarray -> Append(pnums[i] );
	  }
	for ( int i = 0; i < edges.Size(); i++ )
	  {
	    sendarray -> Append(edges[i] );
	    int v1, v2;
	    topology . GetEdgeVertices ( edges[i], v1, v2 );
	    int dv1 = GetLoc2Glob_Vert ( v1 );
	    int dv2 = GetLoc2Glob_Vert ( v2 );
	    sendarray -> Append ( dv1 );
	    sendarray -> Append ( dv2 );
	  }
      }

    BitArray edgeisinit(ned), vertisinit(np);
    edgeisinit.Clear();
    vertisinit.Clear();
    
    // Array for temporary use, to find local from global element fast
    // only for not too big meshes 
    // seems ok, as low-order space is treated on one proc
    Array<int,1> * glob2locfa;
    glob2locfa = new Array<int,1> ( nfaglob );
    (*glob2locfa) = -1;

    for ( int locfa = 1; locfa <= nfa; locfa++)
      {
	if ( !IsExchangeFace ( locfa ) ) continue;
	(*glob2locfa)[GetDistantFaceNum(0, locfa) ] = locfa;
      }

    for ( int sender = 1; sender < ntasks; sender ++ )
      {
		
	if ( id == sender )
	  MyMPI_Bcast ( *sendarray, sender-1, MPI_HIGHORDER_COMM);

	if ( id != sender )
	  {
	    MyMPI_Bcast ( *recvarray, sender-1, MPI_HIGHORDER_COMM);
	    // compare received vertices with own ones
	    int ii = 0;
	    int cntel = 0;
	    int locfa = 1;

	    while ( ii< recvarray -> Size() )
	      {
		
		// receive list : 
		// distant facenum
		// global facenum
		// np
		// ned
		// globalpnums
		// distant pnums
		
		// distant edgenums mit globalv1, globalv2
		
		int distfa = (*recvarray)[ii++];
		int globfa = (*recvarray)[ii++];
		int distnp = (*recvarray)[ii++];
		int distned =(*recvarray)[ii++];
		
		int locfa = (*glob2locfa) [globfa];
		
		if ( locfa == -1 ) 
		  {
		    ii += 2*distnp + 3*distned;
		    locfa = 1;
		    continue;
		  }
		
		Array<int> edges;
		int fa = locfa;
		
		Array<int> pnums, globalpnums;
		topology.GetFaceEdges ( fa, edges );
		topology.GetFaceVertices ( fa, pnums );
		
		
		globalpnums.SetSize ( distnp );
		for ( int i = 0; i < distnp; i++)
		  globalpnums[i] = GetLoc2Glob_Vert ( pnums[i] );
		
		SetDistantFaceNum ( sender, fa, distfa );
		
		// find exchange points
		for ( int i = 0;  i < distnp; i++)
		  {
		    int distglobalpnum = (*recvarray)[ii+i];
		    for ( int j = 0; j < distnp; j++ )
		      if ( globalpnums[j] == distglobalpnum )
			{
			  // set sender -- distpnum  ---- locpnum
			  int distpnum = (*recvarray)[ii + i +distnp];
			  SetDistantPNum ( sender, pnums[j], distpnum );
			}
		    
		  }
			
		int * distedgenums  = new int [distned];
		// find exchange edges
		for ( int i = 0; i  < edges.Size(); i++)
		  {
		    int v1, v2;
		    topology . GetEdgeVertices ( edges[i], v1, v2 );
		    int dv1 = GetLoc2Glob_Vert ( v1 );
		    int dv2 = GetLoc2Glob_Vert ( v2 );
		    if ( dv1 > dv2 ) swap ( dv1, dv2 );
		    for ( int ed = 0; ed < distned; ed++)
		      {
			distedgenums[ed] = (*recvarray)[ii + 2*distnp + 3*ed];
			int ddv1 = (*recvarray)[ii + 2*distnp + 3*ed + 1];
			int ddv2 = (*recvarray)[ii + 2*distnp + 3*ed + 2];
			if ( ddv1 > ddv2 ) swap ( ddv1, ddv2 );
			if ( dv1 == ddv1 && dv2 == ddv2 )
			  {
			    // set sender -- distednum -- locednum
			    SetDistantEdgeNum ( sender, edges[i], distedgenums[ed] );
			  }
		      }
		    
		    
		  }
		delete [] distedgenums; 
		
		ii += 2*distnp + 3*distned;
		
	      }
	  }
      }
  

    // set which elements are where for the master processor

    delete sendarray; delete recvarray;
    if ( id > 0 )
      delete glob2locfa;
    coarseupdate = 1;
	
#ifdef SCALASCA
#pragma pomp inst end(updatecoarsegrid)
#endif
  }

  void ParallelMeshTopology :: UpdateCoarseGridOverlap ()
  {

    UpdateCoarseGridGlobal();

#ifdef SCALASCA
#pragma pomp inst begin(updatecoarsegrid)
#endif
    (*testout) << "UPDATE COARSE GRID PARALLEL TOPOLOGY, OVERLAP " << endl;
    PrintMessage ( 1, "UPDATE COARSE GRID PARALLEL TOPOLOGY, OVERLAP " );

    const MeshTopology & topology = mesh.GetTopology();

    nfa = topology . GetNFaces();
    ned = topology . GetNEdges();
    np = mesh . GetNP();
    nv = mesh . GetNV();
    ne = mesh . GetNE();
    nseg = mesh.GetNSeg();
    nsurfel = mesh.GetNSE();


    if ( id != 0 )
      {

    // find exchange edges - first send exchangeedges locnum, v1, v2
    // receive distant distnum, v1, v2
    // find matching

    Array<int> * sendarray, *recvarray;
    sendarray = new Array<int> (0);
    recvarray = new Array<int>;
    
    sendarray->SetSize (0);

    BitArray recvface(nfa);
    recvface.Clear();
   
    for ( int el = 1; el <= ne; el++ )
      {
	Array<int> edges, pnums, faces;
	topology.GetElementFaces (el, faces);
	int globeli = GetLoc2Glob_VolEl(el);
	for ( int fai = 0; fai < faces.Size(); fai++)
	  {
	    int fa = faces[fai];

	    topology.GetFaceEdges ( fa, edges );
	    topology.GetFaceVertices ( fa, pnums );
	    
	    if ( !IsExchangeElement ( el ) ) continue;

	    int globfa = GetDistantFaceNum(0, fa) ;

	    // send : 
	    // localfacenum
	    // globalfacenum
	    // globalelnum
	    // np
	    // ned
	    // globalpnums
	    // localpnums
	    
	    // localedgenums mit globalelnums mit globalv1, globalv2
	    // 

	    sendarray -> Append ( fa );
	    sendarray -> Append ( globfa );
	    sendarray -> Append ( globeli );
	    sendarray -> Append ( pnums.Size() );
	    sendarray -> Append ( edges.Size() );
	    for ( int i = 0; i < pnums.Size(); i++ )
	      {
		sendarray -> Append( GetLoc2Glob_Vert(pnums[i]) );
	      }
	    for ( int i = 0; i < pnums.Size(); i++ )
	      {
		sendarray -> Append(pnums[i] );
  	      }
	    for ( int i = 0; i < edges.Size(); i++ )
	      {
		int globedge = GetDistantEdgeNum(0, edges[i] );
		int v1, v2;
		topology . GetEdgeVertices ( edges[i], v1, v2 );
		int dv1 = GetLoc2Glob_Vert ( v1 );
		int dv2 = GetLoc2Glob_Vert ( v2 );

		sendarray -> Append(edges[i] );
		sendarray -> Append (globedge);
		sendarray -> Append ( dv1 );
		sendarray -> Append ( dv2 );
	      }
	  }
      }
    
    BitArray edgeisinit(ned), vertisinit(np);
    edgeisinit.Clear();
    vertisinit.Clear();
    
    // Array for temporary use, to find local from global element fast
    // only for not too big meshes 
    // seems ok, as low-order space is treated on one proc
    Array<int,1> * glob2loc_el;

    glob2loc_el = new Array<int,1> ( neglob );
    (*glob2loc_el) = -1;
    for ( int locel = 1; locel <= mesh.GetNE(); locel++)
      (*glob2loc_el)[GetLoc2Glob_VolEl(locel)] = locel;
      
    for ( int sender = 1; sender < ntasks; sender ++ )
      {
	if ( id == sender )
	  MyMPI_Bcast (*sendarray, sender-1, MPI_HIGHORDER_COMM);

// 	  {
// 	    for ( int dest = 1; dest < ntasks; dest ++ )
// 	      if ( dest != id)
// 		{
// 		  MyMPI_Send (*sendarray, dest);
// 		}
// 	  }
	    
	if ( id != sender )
	  {
// 	    MyMPI_Recv ( *recvarray, sender);
	    MyMPI_Bcast (*recvarray, sender-1, MPI_HIGHORDER_COMM);
	    // compare received vertices with own ones
	    int ii = 0;
	    int cntel = 0;
	    int volel = 1;

	    while ( ii< recvarray -> Size() )
	      {

		// receive list : 
		// distant facenum
		// np
		// ned
		// globalpnums
		// distant pnums

		// distant edgenums mit globalv1, globalv2

		int distfa = (*recvarray)[ii++];
		int globfa = (*recvarray)[ii++];
		int globvolel = (*recvarray)[ii++];
		int distnp = (*recvarray)[ii++];
		int distned =(*recvarray)[ii++];

		if ( id > 0 ) // GetLoc2Glob_VolEl ( volel ) != globvolel )
		  volel = (*glob2loc_el)[globvolel]; //Glob2Loc_VolEl ( globvolel );
		else
		  volel = globvolel;

		if ( volel == -1 ) 
		  {
		    ii += 2*distnp + 4*distned;
		    volel = 1;
		    continue;
		  }

		Array<int> faces, edges;
		topology.GetElementFaces( volel, faces);
		topology.GetElementEdges ( volel, edges);
		for ( int fai= 0; fai < faces.Size(); fai++ )
		  {
		    int fa = faces[fai];
		    if ( !IsExchangeFace ( fa ) && sender != 0 ) continue;
		    //		    if ( recvface.Test ( fa-1 ) ) continue;

		    Array<int> pnums, globalpnums;
		    //topology.GetFaceEdges ( fa, edges );
		    topology.GetFaceVertices ( fa, pnums );


		    // find exchange faces ...
		    // have to be of same type
		    if ( pnums.Size () != distnp ) continue;


		    globalpnums.SetSize ( distnp );
		    for ( int i = 0; i < distnp; i++)
		      globalpnums[i] = GetLoc2Glob_Vert ( pnums[i] );





		    // test if 3 vertices match
		    bool match = 1;
		    for ( int i = 0;  i < distnp; i++)
		      if ( !globalpnums.Contains ( (*recvarray)[ii+i] ) )
			match = 0;
		    
		    if ( !match ) continue;

		    //  recvface.Set(fa-1);

		    SetDistantFaceNum ( sender, fa, distfa );

		    SetDistantFaceNum ( 0, fa, globfa );

		    // find exchange points
		    for ( int i = 0;  i < distnp; i++)
		      {
			int distglobalpnum = (*recvarray)[ii+i];
			for ( int j = 0; j < distnp; j++ )
			  if ( globalpnums[j] == distglobalpnum )
			    {
			      // set sender -- distpnum  ---- locpnum
			      int distpnum = (*recvarray)[ii + i +distnp];
			      SetDistantPNum ( sender, pnums[j], distpnum );
			    }
			
		      }

		    int * distedgenums  = new int [distned];
		    // find exchange edges
      		    for ( int i = 0; i  < edges.Size(); i++)
		      {
			int v1, v2;
			topology . GetEdgeVertices ( edges[i], v1, v2 );
			int dv1 = GetLoc2Glob_Vert ( v1 );
			int dv2 = GetLoc2Glob_Vert ( v2 );
			if ( dv1 > dv2 ) swap ( dv1, dv2 );
			for ( int ed = 0; ed < distned; ed++)
			  {
			    distedgenums[ed] = (*recvarray)[ii + 2*distnp + 4*ed];
			    int globedgenum = (*recvarray)[ii + 2*distnp + 4*ed + 1];
			    int ddv1 = (*recvarray)[ii + 2*distnp + 4*ed + 2];
			    int ddv2 = (*recvarray)[ii + 2*distnp + 4*ed + 3];
			    if ( ddv1 > ddv2 ) swap ( ddv1, ddv2 );
			    if ( dv1 == ddv1 && dv2 == ddv2 )
			      {
				// set sender -- distednum -- locednum
				SetDistantEdgeNum ( sender, edges[i], distedgenums[ed] );
				SetDistantEdgeNum ( 0, edges[i], globedgenum );
			      }
			  }


		      }
		    delete [] distedgenums; 


		  }
		ii += 2*distnp + 4*distned;
	      }
		
	    
	    
	     
          }
	
      }

    // set which elements are where for the master processor

    delete sendarray; delete recvarray;
    if ( id > 0 )
      delete glob2loc_el;
    coarseupdate = 1;
	
      }


    // send global-local el/face/edge/vert-info to id 0


//     nfa = topology . GetNFaces();
//     ned = topology . GetNEdges();
//     np = mesh . GetNP();
//     nv = mesh . GetNV();
//     ne = mesh . GetNE();
//     nseg = mesh.GetNSeg();
//     nsurfel = mesh.GetNSE();
    if ( id != 0 )
      {
	Array<int> * sendarray;
	sendarray = new Array<int> (4);

	int sendnfa = 0, sendned = 0;

	(*sendarray)[0] = ne;
	(*sendarray)[1] = nfa;
	(*sendarray)[2] = ned;
	(*sendarray)[3] = np;

	int ii = 4;
	for ( int el = 1; el <= ne; el++ )
	  (*sendarray).Append ( GetLoc2Glob_VolEl (el ) );

	for ( int fa = 1; fa <= nfa; fa++ )
	  {
	    if ( !IsExchangeFace (fa) ) continue;
	    sendnfa++;
	    (*sendarray).Append ( fa );
	    (*sendarray).Append ( GetDistantFaceNum (0, fa) );
	  }

	for ( int ed = 1; ed <= ned; ed++ )
	  {
	    if ( !IsExchangeEdge (ed) ) continue;
	    sendned++;
	    sendarray->Append ( ed );
	    sendarray->Append ( GetDistantEdgeNum(0, ed) );
	  }

	for ( int vnum = 1; vnum <= np; vnum++ )
	  sendarray->Append ( GetLoc2Glob_Vert(vnum) );

	(*sendarray)[1] = sendnfa;
	(*sendarray)[2] = sendned;

	MyMPI_Send (*sendarray, 0);

	delete sendarray;
      }

    else
      {
	Array<int> * recvarray = new Array<int>;

	for ( int sender = 1; sender < ntasks; sender++ )
	  {
	    MyMPI_Recv ( *recvarray, sender);

	    int distnel = (*recvarray)[0];
	    int distnfa = (*recvarray)[1];
	    int distned = (*recvarray)[2];
	    int distnp = (*recvarray)[3];

	    int ii = 4;

	    for ( int el = 1; el <= distnel; el++ )
	      SetDistantEl ( sender, (*recvarray)[ii++], el );

	    for ( int fa = 1; fa <= distnfa; fa++ )
	      {
		int distfa = (*recvarray)[ii++];
		SetDistantFaceNum ( sender, (*recvarray)[ii++], distfa );
	      }
	    for ( int ed = 1; ed <= distned; ed++ )
	      {
		int disted = (*recvarray)[ii++];
		SetDistantEdgeNum ( sender, (*recvarray)[ii++], disted );
	      }
	    for ( int vnum = 1; vnum <= distnp; vnum++ )
	      SetDistantPNum ( sender, (*recvarray)[ii++], vnum );
	  }

	delete recvarray;
      }    
#ifdef SCALASCA
#pragma pomp inst end(updatecoarsegrid)
#endif
  }



  void ParallelMeshTopology :: UpdateTopology () 
  {
    // loop over parallel faces and edges, find new local face/edge number, 

    const MeshTopology & topology = mesh.GetTopology();
    int nfa = topology.GetNFaces();
    int ned = topology.GetNEdges();

    isghostedge.SetSize(ned);
    isghostface.SetSize(nfa);
    isghostedge.Clear();
    isghostface.Clear();

    for ( int ed = 1; ed <= ned; ed++)
      {
	int v1, v2;
	topology.GetEdgeVertices ( ed, v1, v2 );
	if ( IsGhostVert(v1) || IsGhostVert(v2) )
	  SetGhostEdge ( ed );
      }


    Array<int> pnums;
    for ( int fa = 1; fa <= nfa; fa++)
      {
	topology.GetFaceVertices ( fa, pnums );
	for ( int i = 0; i < pnums.Size(); i++)
	  if ( IsGhostVert( pnums[i] ) )
	    {
	      SetGhostFace ( fa );
	      break;
	    }
      }
  }


void ParallelMeshTopology :: UpdateExchangeElements()
{
  (*testout) << "UPDATE EXCHANGE ELEMENTS " << endl;
  const MeshTopology & topology = mesh.GetTopology();

  isexchangeedge->SetSize ( (ntasks+1) * topology.GetNEdges() );
  isexchangeface->SetSize ( (ntasks+1) * topology.GetNFaces() );

  isexchangeedge->Clear();
  isexchangeface->Clear();

  for ( int eli = 1; eli <= mesh.GetNE(); eli++)
    {
      if ( ! IsExchangeElement ( eli ) ) continue;
      const Element & el = mesh.VolumeElement(eli);
      Array<int> faces, edges;
      int np = el.NP();

      topology.GetElementEdges ( eli, edges );
      topology.GetElementFaces ( eli, faces );
      for ( int i = 0; i < edges.Size(); i++)
	{
	  SetExchangeEdge ( edges[i] );
	}
      for ( int i = 0; i < faces.Size(); i++)
	{
	  SetExchangeFace ( faces[i] );
	}
      for ( int i = 0; i < np; i++)
	{
	  SetExchangeVert ( el[i] );
	}
    }

  if ( id == 0 ) return;



  Array<int> ** elementonproc, ** recvelonproc;
  elementonproc = new Array<int>*[ntasks];
  recvelonproc = new Array<int>*[ntasks];

  for ( int i = 1; i < ntasks; i++ )
    {
      elementonproc[i] = new Array<int>(0);
      recvelonproc[i] = new Array<int>(0);
    }


  for ( int eli = 1; eli <= mesh.GetNE(); eli++ )
    {
      if ( !IsExchangeElement(eli) ) continue;
      for ( int i = 1; i < ntasks; i++ )
	if ( GetDistantElNum(i, eli) != -1 && i != id ) 
	  {
	    elementonproc[i] -> Append(eli);
	    elementonproc[i] -> Append(GetDistantElNum(i, eli));
	  }

    }

  for ( int sender = 1; sender < ntasks; sender ++ )
    {
      if ( id == sender )
	for ( int dest = 1; dest < ntasks; dest ++ )
	  if ( dest != id)
	    {
	      MyMPI_Send ( *(elementonproc[dest]), dest);
	      elementonproc[dest] -> SetSize(0);
	    }
	
      
      if ( id != sender )
	{
	  MyMPI_Recv (*( recvelonproc[sender]), sender);
	}
    }


  int ii = 0;
  for ( int sender = 1; sender < ntasks; sender++ )
    {
      if ( sender == id ) continue;

      ii = 0;
      while ( recvelonproc[sender]->Size() > ii )
	{
	  int distelnum = (*recvelonproc[sender])[ii++];
	  int locelnum =  (*recvelonproc[sender])[ii++];
	  SetDistantEl ( sender, locelnum, distelnum);
	}
      recvelonproc[sender]->SetSize(0);
    }


  BitArray procs(ntasks);
  procs.Clear();
  for ( int eli = 1; eli <= mesh.GetNE(); eli++)
    {
      if ( IsGhostEl(eli) ) continue;
      if ( !IsExchangeElement(eli) ) continue;

      procs.Clear();
      int sumprocs = 0;
      for ( int i = 1; i < ntasks; i++ )
	if ( GetDistantElNum(i, eli) != -1 && i != id ) 
	  {
	    procs.Set(i);
	    sumprocs++;
	  }

      for ( int dest = 1; dest < ntasks; dest++)
	{
	  if ( !procs.Test(dest) ) continue;
	  elementonproc[dest]->Append(GetDistantElNum(dest, eli));
	  elementonproc[dest]->Append(sumprocs);
	  for ( int i = 1; i < ntasks; i++ )
	    if ( procs.Test(i) )
	      { 
		elementonproc[dest]->Append(i);
		elementonproc[dest]->Append(GetDistantElNum(i, eli));
	      }
	}
    }

  for ( int sender = 1; sender < ntasks; sender ++ )
    {
      if ( id == sender )
	for ( int dest = 1; dest < ntasks; dest ++ )
	  if ( dest != id)
	    {
	      MyMPI_Send ( *(elementonproc[dest]), dest);
	      delete elementonproc[dest];
	    }
	
      
      if ( id != sender )
	{
	  MyMPI_Recv (*( recvelonproc[sender]), sender);
	}
    }

  for ( int sender = 1; sender < ntasks; sender++ )
    {
      if ( sender == id ) continue;
      ii = 0;
      while ( recvelonproc[sender]->Size() > ii )
	{
	  int locelnum = (*recvelonproc[sender])[ii++];
	  int nprocs =  (*recvelonproc[sender])[ii++];
	  for ( int iproc = 0; iproc < nprocs; iproc++)
	    {
	      int proc = (*recvelonproc[sender])[ii++];
 	      int distelnum = (*recvelonproc[sender])[ii++];
	      if ( id == proc ) continue;
	      SetExchangeElement (locelnum, proc);
 	      SetDistantEl( proc, locelnum, distelnum );
	    }
	}
      delete recvelonproc[sender];
    }

  delete [] elementonproc;
  delete [] recvelonproc;
}




void ParallelMeshTopology :: SetNV ( const int anv )
  {
    *testout << "called setnv"  << endl
             << "old size: " << loc2distvert.Size() << endl
             << "new size: " << anv << endl;

    loc2distvert.ChangeSize (anv);
    for (int i = 1; i <= anv; i++)
      if (loc2distvert.EntrySize(i) == 0)
        loc2distvert.Add (i, -1);  // will be the global nr

    BitArray * isexchangevert2 = new BitArray( (ntasks+1) * anv );
    isexchangevert2->Clear();
    if ( isexchangevert )
      {
	for ( int i = 0; i < min2( isexchangevert->Size(), isexchangevert2->Size() ); i++ )
	  if ( isexchangevert->Test(i) ) isexchangevert2->Set(i);
	delete isexchangevert;
      }

    isexchangevert = isexchangevert2;
    nv = anv;

  }

void ParallelMeshTopology :: SetNE ( const int ane )
  {
    loc2distel.ChangeSize (ane);
    for (int i = 0; i < ane; i++)
      {
	if (loc2distel[i].Size() == 0)
	  loc2distel.Add (i, -1);   // will be the global nr
      }
    BitArray * isexchangeel2 = new BitArray ( (ntasks+1) * ane );
    isexchangeel2->Clear();
    if ( isexchangeel )
      {
	for ( int i = 0; i < min2(isexchangeel->Size(), isexchangeel2->Size() ) ; i++ )
	  if ( isexchangeel->Test(i) ) isexchangeel2->Set(i);
	delete isexchangeel;
      }

    ne = ane;
    isexchangeel = isexchangeel2;

  }

void ParallelMeshTopology :: SetNSE ( int anse )
  {
    loc2distsurfel.ChangeSize (anse);
    for (int i = 0; i < anse; i++)
      if (loc2distsurfel[i].Size() == 0)
        loc2distsurfel.Add (i, -1);  // will be the global nr

   nsurfel = anse;
  }

void ParallelMeshTopology :: SetNSegm ( int anseg )
{
  loc2distsegm.ChangeSize (anseg);
  for (int i = 0; i < anseg; i++)
    if (loc2distsegm[i].Size() == 0)
        loc2distsegm.Add (i, -1);  // will be the global nr

   nseg = anseg;
  }


}




#endif
