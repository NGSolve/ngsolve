#ifdef PARALLEL


#include <meshing.hpp>
#include "paralleltop.hpp"


namespace netgen
{

  static MPI_Group MPI_HIGHORDER_WORLD;
  static MPI_Comm MPI_HIGHORDER_COMM;

  void MyMPI_ExchangeTable (TABLE<int> & send_verts, 
			    TABLE<int> & recv_verts, int tag,
			    MPI_Comm comm = MPI_COMM_WORLD)
  {
    int ntasks, rank;
    MPI_Comm_size(comm, &ntasks);
    MPI_Comm_rank(comm, &rank);

    Array<MPI_Request> requests;
    for (int dest = 0; dest < ntasks; dest++)
      if (dest != rank)
	requests.Append (MyMPI_ISend (send_verts[dest], dest, tag, comm));

    for (int i = 0; i < ntasks-1; i++)
      {
	MPI_Status status;
	MPI_Probe (MPI_ANY_SOURCE, tag, comm, &status);
	int size, src = status.MPI_SOURCE;
	MPI_Get_count (&status, MPI_INT, &size);
	recv_verts.SetEntrySize (src, size, sizeof(int));
	requests.Append (MyMPI_IRecv (recv_verts[src], src, tag, comm));
      }
    MPI_Waitall (requests.Size(), &requests[0], MPI_STATUS_IGNORE);
  }




  void ParallelMeshTopology :: Reset ()
  {
    *testout << "ParallelMeshTopology::Reset" << endl;
    
    if ( ntasks == 1 ) return;
    
    int nvold = nv;
    
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

    if ( nvold == nv ) return;

    SetNV(nv);
    SetNE(ne);
  }


  ParallelMeshTopology :: ~ParallelMeshTopology ()
  {
    ;
  }



  ParallelMeshTopology :: ParallelMeshTopology ( const netgen::Mesh & amesh )
    : mesh(amesh)
  {
    ned = 0; 
    nfa = 0; 
    nv = 0;
    ne = 0;
    np = 0;
    nseg = 0;
    nsurfel = 0;
    neglob = 0;
    nvglob = 0;

    nparel = 0;

    coarseupdate = 0;
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
		(*testout) << "   p" << dest << ": " << GetDistantEdgeNum ( dest, i ) << endl;
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


    /*
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
    */
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



  // gibt anzahl an distant pnums zurueck
  // * pnums entspricht Array<int[2] >
  int  ParallelMeshTopology :: GetDistantPNums ( int locpnum, int * distpnums ) const
  {
    distpnums[0] = 0;
    distpnums[1] = loc2distvert[locpnum][0];
    for ( int i = 1; i < loc2distvert[locpnum].Size(); i++ )
      distpnums[i+1] = loc2distvert[locpnum][i];

    int size = loc2distvert[locpnum].Size() / 2 + 1;
    return size;
  } 

  int  ParallelMeshTopology :: GetDistantFaceNums ( int locfacenum, int * distfacenums ) const
  {
    distfacenums[0] = 0;
    distfacenums[1] = loc2distface[locfacenum-1][0];

    for ( int i = 1; i < loc2distface[locfacenum-1].Size(); i++ )
      distfacenums[i+1] = loc2distface[locfacenum-1][i];

    int size = loc2distface[locfacenum-1].Size() / 2 + 1;
    return size;
  } 

  int  ParallelMeshTopology :: GetDistantEdgeNums ( int locedgenum, int * distedgenums ) const
  {
    distedgenums[0] = 0;
    distedgenums[1] = loc2distedge[locedgenum-1][0];

    for ( int i = 1; i < loc2distedge[locedgenum-1].Size(); i++ )
      distedgenums[i+1] = loc2distedge[locedgenum-1][i];

    int size = loc2distedge[locedgenum-1].Size() / 2 + 1;
    return size;
  } 

  int  ParallelMeshTopology :: GetDistantElNums ( int locelnum, int * distelnums ) const
  {
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
    if (id == 0)
      PrintMessage ( 3, "UPDATE GLOBAL COARSEGRID STARTS" );      // JS


    MPI_Group MPI_GROUP_WORLD;
    int process_ranks[] = { 0 };
    MPI_Comm_group (MPI_COMM_WORLD, &MPI_GROUP_WORLD);
    MPI_Group_excl (MPI_GROUP_WORLD, 1, process_ranks, &MPI_HIGHORDER_WORLD);
    MPI_Comm_create (MPI_COMM_WORLD, MPI_HIGHORDER_WORLD, &MPI_HIGHORDER_COMM);


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
	for ( int eli = 1; eli <= ne; eli++ )
	  loc2distel[eli-1][0] = eli;
      
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
      {
	sendarray.Append (nfa);
	sendarray.Append (ned);

	Array<int> edges, pnums, faces, elpnums;
	sendarray.Append (ne);
	for ( int el = 1; el <= ne; el++ )
	  {
	    topology.GetElementFaces (el, faces);
	    topology.GetElementEdges ( el, edges );
	    const Element & volel = mesh.VolumeElement (el);

	    int globeli = GetLoc2Glob_VolEl(el);
	    // cout << "el = " << el << ", globeli = " << globeli << endl;

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


	sendarray.Append (nsurfel);
	for (int el = 1; el <= nsurfel; el++)
	  {
	    topology.GetSurfaceElementEdges ( el, edges );
	    const Element2d & surfel = mesh.SurfaceElement (el);

	    sendarray. Append ( el );
	    sendarray. Append ( edges.Size() );
	    sendarray. Append ( surfel.GetNP() );

	    for ( int i = 0; i < edges.Size(); i++ )
	      sendarray. Append(edges[i] );

	    for ( int i = 0; i < surfel.GetNP(); i++ )
	      sendarray. Append(surfel[i] );
	  }

      }


    if (id == 0)
      MyMPI_Bcast ( sendarray );
    else
      MyMPI_Bcast ( recvarray );


    if (id != 0)
      {
	Array<int,1> glob2loc_el;

	glob2loc_el.SetSize (neglob);  
	glob2loc_el = -1;
	for ( int locel = 1; locel <= mesh.GetNE(); locel++)
	  glob2loc_el[GetLoc2Glob_VolEl(locel)] = locel;

	int ii = 0;
	nfaglob = recvarray[ii++];
	nedglob = recvarray[ii++];

	Array<int> faces, edges;
	Array<int> pnums, globalpnums;

	int recv_ne = recvarray[ii++];
	for (int hi = 0; hi < recv_ne; hi++)
	  {
	    int globvolel = recvarray[ii++];
	    int distnfa = recvarray[ii++];
	    int distned = recvarray[ii++];
	    int distnp  = recvarray[ii++];

	    int volel = glob2loc_el[globvolel];
	    if (volel != -1)
	      {
		topology.GetElementFaces( volel, faces);
		topology.GetElementEdges ( volel, edges);
		const Element & volelement = mesh.VolumeElement (volel);

		for ( int i = 0; i  < faces.Size(); i++)
		  SetDistantFaceNum ( 0, faces[i], recvarray[ii++]);
		
		for ( int i = 0; i  < edges.Size(); i++)
		  SetDistantEdgeNum ( 0, edges[i], recvarray[ii++]);

		for ( int i = 0; i  < volelement.GetNP(); i++)
		  SetDistantPNum ( 0, volelement[i], recvarray[ii++]);
	      }
	    else
	      ii += distnfa + distned + distnp;
	  }

	
	Array<int,1> glob2loc_sel;

	int recv_nse = recvarray[ii++];
	nseglob = recv_nse;

	glob2loc_sel.SetSize (nseglob);  
	glob2loc_sel = -1;
	for ( int locel = 1; locel <= mesh.GetNSE(); locel++)
	  glob2loc_sel[GetLoc2Glob_SurfEl(locel)] = locel;


	for (int hi = 0; hi < recv_nse; hi++)
	  {
	    int globvolel = recvarray[ii++];
	    int distned = recvarray[ii++];
	    int distnp  = recvarray[ii++];

	    int surfel = glob2loc_sel[globvolel];
	    if (surfel != -1)
	      {
		topology.GetSurfaceElementEdges ( surfel, edges);
		const Element2d & element = mesh.SurfaceElement (surfel);
		
		for ( int i = 0; i  < edges.Size(); i++)
		  SetDistantEdgeNum ( 0, edges[i], recvarray[ii++]);

		for ( int i = 0; i  < element.GetNP(); i++)
		  SetDistantPNum ( 0, element[i], recvarray[ii++]);
	      }
	    else
	      ii += distned + distnp;
	  }


      }
    
    if (id != 0)
      {
	*testout << "l2d - vert = " << loc2distvert << endl;
	*testout << "l2d - edge = " << loc2distedge << endl;
	*testout << "l2d - el = " << loc2distel << endl;
	*testout << "l2d - sel = " << loc2distsurfel << endl;
      }

    coarseupdate = 1;
  }






  void ParallelMeshTopology :: UpdateCoarseGrid ()
  {
    static int timer = NgProfiler::CreateTimer ("UpdateCoarseGrid");
    NgProfiler::RegionTimer reg(timer);


    (*testout) << "UPDATE COARSE GRID PARALLEL TOPOLOGY " << endl;
    if (id == 0)
      PrintMessage (1, "UPDATE COARSE GRID PARALLEL TOPOLOGY ");


    // find exchange edges - first send exchangeedges locnum, v1, v2
    // receive distant distnum, v1, v2
    // find matching
    const MeshTopology & topology = mesh.GetTopology();

    UpdateCoarseGridGlobal();
    
    MPI_Barrier (MPI_COMM_WORLD);

    if ( id == 0 ) 
      {
	return;
      }


    Array<int> sendarray, recvarray;

    nfa = topology . GetNFaces();
    ned = topology . GetNEdges();
    np = mesh . GetNP();
    nv = mesh . GetNV();
    ne = mesh . GetNE();
    nseg = mesh.GetNSeg();
    nsurfel = mesh.GetNSE();
    

    static int timerv = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex vertices");
    static int timere = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex edges");
    static int timerf = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex faces");


    Array<int,1> glob2loc;
    Array<int> cnt_send(ntasks-1);

    NgProfiler::StartTimer (timere);

    // exchange edges
    int maxedge = 0;
    for (int edge = 1; edge <= ned; edge++)
      maxedge = max (maxedge, GetDistantEdgeNum (0, edge));

    glob2loc.SetSize (maxedge);
    glob2loc = -1;
    
    for (int edge = 1; edge <= ned; edge++)
      glob2loc[GetDistantEdgeNum(0, edge)] = edge;

    cnt_send = 0;
    int v1, v2;
    for (int edge = 1; edge <= ned; edge++)
      {
	topology.GetEdgeVertices (edge, v1, v2);
	for (int dest = 1; dest < ntasks; dest++)
	  if (IsExchangeVert (dest, v1) && 
	      IsExchangeVert (dest, v2))
	    {
	      cnt_send[dest-1]+=2;
	    }
      }
    
    TABLE<int> send_edges(cnt_send);
    for (int edge = 1; edge <= ned; edge++)
      {
	topology.GetEdgeVertices (edge, v1, v2);
	for (int dest = 1; dest < ntasks; dest++)
	  {
	    if (IsExchangeVert (dest, v1) && 
		IsExchangeVert (dest, v2))
	      {
		send_edges.Add (dest-1, GetDistantEdgeNum(0, edge));
		send_edges.Add (dest-1, edge);
	      }
	  }
      }


    // *testout << "send exchange edges: " << send_edges << endl;

    TABLE<int> recv_edges(ntasks-1);
    MyMPI_ExchangeTable (send_edges, recv_edges, MPI_TAG_MESH+9, MPI_HIGHORDER_COMM);

    // *testout << "recv exchange edges: " << recv_edges << endl;

    for (int sender = 1; sender < ntasks; sender ++)
      if (id != sender)
	{
	  FlatArray<int> recvarray = recv_edges[sender-1];

	  for (int ii = 0; ii < recvarray.Size(); )
	    { 
	      int globe = recvarray[ii++];
	      int diste = recvarray[ii++];

	      if (globe <= maxedge)
		{
		  int loce = glob2loc[globe];
		  if (loce != -1)
		    SetDistantEdgeNum (sender, loce, diste);
		}
	    }
	} 


 
    NgProfiler::StopTimer (timere);


    MPI_Barrier (MPI_HIGHORDER_COMM);

    if (mesh.GetDimension() == 3)
      {

    NgProfiler::StartTimer (timerf);

    glob2loc.SetSize (nfaglob);
    glob2loc = -1;
    
    for (int loc = 1; loc <= nfa; loc++)
      glob2loc[GetDistantFaceNum(0, loc)] = loc;

    cnt_send = 0;
    Array<int> verts;
    for (int face = 1; face <= nfa; face++)
      {
	topology.GetFaceVertices (face, verts);
	for (int dest = 1; dest < ntasks; dest++)
	  if (IsExchangeVert (dest, verts[0]) && 
	      IsExchangeVert (dest, verts[1]) &&
	      IsExchangeVert (dest, verts[2]))
	    {
	      cnt_send[dest-1]+=2;
	    }
      }
    
    TABLE<int> send_faces(cnt_send);
    for (int face = 1; face <= nfa; face++)
      {
	topology.GetFaceVertices (face, verts);
	for (int dest = 1; dest < ntasks; dest++)
	  {
	    if (IsExchangeVert (dest, verts[0]) && 
		IsExchangeVert (dest, verts[1]) &&
		IsExchangeVert (dest, verts[2]))
	      {
		send_faces.Add (dest-1, GetDistantFaceNum(0, face));
		send_faces.Add (dest-1, face);
	      }
	  }
      }
    TABLE<int> recv_faces(ntasks-1);
    MyMPI_ExchangeTable (send_faces, recv_faces, MPI_TAG_MESH+8, MPI_HIGHORDER_COMM);

    // *testout << "send exchange faces: " << send_faces << endl;
    // *testout << "recv exchange faces: " << recv_faces << endl;

    for (int sender = 1; sender < ntasks; sender ++)
      if (id != sender)
	{
	  FlatArray<int> recvarray = recv_faces[sender-1];

	  for (int ii = 0; ii < recvarray.Size(); )
	    { 
	      int globf = recvarray[ii++];
	      int distf = recvarray[ii++];
	      int locf = glob2loc[globf];

	      // *testout << "set distant face, sender = " << sender << ", locf = " << locf << "; distf = " << distf << endl;
	      if (locf != -1)
		SetDistantFaceNum (sender, locf, distf);
	    }
	} 


    NgProfiler::StopTimer (timerf);
      }


    // set which elements are where for the master processor

    coarseupdate = 1;
  }





  void ParallelMeshTopology :: UpdateTopology () 
  {
    ;
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

    nv = anv;
  }

  void ParallelMeshTopology :: SetNE ( const int ane )
  {
    loc2distel.ChangeSize (ane);
    for (int i = 0; i < ane; i++)
      if (loc2distel[i].Size() == 0)
	loc2distel.Add (i, -1);   // will be the global nr
    ne = ane;
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
