#ifdef PARALLEL


#include <meshing.hpp>
#include "paralleltop.hpp"


namespace netgen
{

  ParallelMeshTopology :: ParallelMeshTopology (const Mesh & amesh)
    : mesh(amesh)
  {
    is_updated = false;
  }
 

  ParallelMeshTopology :: ~ParallelMeshTopology ()
  {
    ;
  }


  void ParallelMeshTopology :: Reset ()
  {
    *testout << "ParallelMeshTopology::Reset" << endl;
    
    if ( ntasks == 1 ) return;

    int ned = mesh.GetTopology().GetNEdges();
    int nfa = mesh.GetTopology().GetNFaces();

    if (glob_edge.Size() != ned)
      {
	glob_edge.SetSize(ned);
	glob_face.SetSize(nfa);
	glob_edge = -1;
	glob_face = -1;

	loc2distedge.ChangeSize (ned);
	loc2distface.ChangeSize (nfa);
      }

    if (glob_vert.Size() != mesh.GetNV())
      {
	SetNV(mesh.GetNV());
	SetNE(mesh.GetNE());
      }
  }


  void ParallelMeshTopology ::  Print() const
  {
    ;
  }


  void ParallelMeshTopology :: SetDistantFaceNum (int dest, int locnum)
  {
    for ( int i = 0; i < loc2distface[locnum-1].Size(); i+=1 )
      if ( loc2distface[locnum-1][i] == dest )
	return;
    loc2distface.Add(locnum-1, dest);
  }

  void ParallelMeshTopology :: SetDistantPNum (int dest, int locnum)
  {
    for ( int i = 0;  i < loc2distvert[locnum-1].Size(); i+=1 )
      if ( loc2distvert[locnum-1][i] == dest )
	return;
    loc2distvert.Add (locnum-1, dest);  
  }


  void ParallelMeshTopology :: SetDistantEdgeNum (int dest, int locnum)
  {
    for ( int i = 0; i < loc2distedge[locnum-1].Size(); i+=1 )
      if ( loc2distedge[locnum-1][i] == dest )
	return;
    loc2distedge.Add (locnum-1, dest);
  }

  void ParallelMeshTopology :: SetNV (int anv)
  {
    glob_vert.SetSize(anv);
    glob_vert = -1;
    loc2distvert.ChangeSize (anv);
  }

  void ParallelMeshTopology :: SetNE ( int ane )
  {
    glob_el.SetSize (ane);
    glob_el = -1;
  }

  void ParallelMeshTopology :: SetNSE ( int anse )
  {
    glob_surfel.SetSize(anse);
    glob_surfel = -1;
  }

  void ParallelMeshTopology :: SetNSegm ( int anseg )
  {
    glob_segm.SetSize (anseg);
    glob_segm = -1;
  }






  void ParallelMeshTopology :: UpdateCoarseGridGlobal ()
  {
    cout << "updatecoarsegridglobal called" << endl;
    if (id == 0)
      PrintMessage ( 3, "UPDATE GLOBAL COARSEGRID STARTS" );      

    int timer = NgProfiler::CreateTimer ("UpdateCoarseGridGlobal");
    NgProfiler::RegionTimer reg(timer);

    *testout << "ParallelMeshTopology :: UpdateCoarseGridGlobal" << endl;

    const MeshTopology & topology = mesh.GetTopology();

    if ( id == 0 )
      {
	Array<Array<int>*> sendarrays(ntasks);
	for (int dest = 1; dest < ntasks; dest++)
	  sendarrays[dest] = new Array<int>;

	Array<int> edges, faces;
	for (int el = 1; el <= mesh.GetNE(); el++)
	  {
	    topology.GetElementFaces (el, faces);
	    topology.GetElementEdges (el, edges);
	    const Element & volel = mesh.VolumeElement (el);

	    Array<int> & sendarray = *sendarrays[volel.GetPartition()];

	    for ( int i = 0; i < edges.Size(); i++ )
	      sendarray.Append (edges[i]);
	    for ( int i = 0; i < faces.Size(); i++ )
	      sendarray.Append (faces[i]);
	  }

	for (int el = 1; el <= mesh.GetNSE(); el++)
	  {
	    topology.GetSurfaceElementEdges (el, edges);
	    const Element2d & surfel = mesh.SurfaceElement (el);
	    Array<int> & sendarray = *sendarrays[surfel.GetPartition()];

	    for ( int i = 0; i < edges.Size(); i++ )
	      sendarray.Append (edges[i]);
	    sendarray.Append (topology.GetSurfaceElementFace (el));
	  }

	Array<MPI_Request> sendrequests;
	for (int dest = 1; dest < ntasks; dest++)
	  sendrequests.Append (MyMPI_ISend (*sendarrays[dest], dest, MPI_TAG_MESH+10));
	MPI_Waitall (sendrequests.Size(), &sendrequests[0], MPI_STATUS_IGNORE);

	for (int dest = 1; dest < ntasks; dest++)
	  delete sendarrays[dest];
      }

    else

      {
	Array<int> recvarray;
	MyMPI_Recv (recvarray, 0, MPI_TAG_MESH+10);

	int ii = 0;

	Array<int> faces, edges;

	for (int volel = 1; volel <= mesh.GetNE(); volel++)
	  {
	    topology.GetElementEdges ( volel, edges);
	    for ( int i = 0; i  < edges.Size(); i++)
	      SetLoc2Glob_Edge ( edges[i], recvarray[ii++]);

	    topology.GetElementFaces( volel, faces);
	    for ( int i = 0; i  < faces.Size(); i++)
	      SetLoc2Glob_Face ( faces[i], recvarray[ii++]);
	  }

	for (int surfel = 1; surfel <= mesh.GetNSE(); surfel++)
	  {
	    topology.GetSurfaceElementEdges (surfel, edges);
	    for (int i = 0; i  < edges.Size(); i++)
	      SetLoc2Glob_Edge (edges[i], recvarray[ii++]);
	    int face = topology.GetSurfaceElementFace (surfel);
	    SetLoc2Glob_Face ( face, recvarray[ii++]);
	  }
      }
    
    is_updated = true;
  }




  void ParallelMeshTopology :: UpdateCoarseGrid ()
  {
    cout << "UpdateCoarseGrid" << endl;
    // if (is_updated) return;

    Reset();
    static int timer = NgProfiler::CreateTimer ("UpdateCoarseGrid");
    NgProfiler::RegionTimer reg(timer);


    (*testout) << "UPDATE COARSE GRID PARALLEL TOPOLOGY " << endl;
    if (id == 0)
      PrintMessage (1, "update parallel topology");


    // UpdateCoarseGridGlobal();


    
    // MPI_Barrier (MPI_COMM_WORLD);

    MPI_Group MPI_GROUP_WORLD;
    MPI_Group MPI_LocalGroup;
    MPI_Comm MPI_LocalComm;

    int process_ranks[] = { 0 };
    MPI_Comm_group (MPI_COMM_WORLD, &MPI_GROUP_WORLD);
    MPI_Group_excl (MPI_GROUP_WORLD, 1, process_ranks, &MPI_LocalGroup);
    MPI_Comm_create (MPI_COMM_WORLD, MPI_LocalGroup, &MPI_LocalComm);

    if (id == 0) return;

    const MeshTopology & topology = mesh.GetTopology();

    Array<int> cnt_send(ntasks-1);


    // update new vertices after mesh-refinement
    if (mesh.mlbetweennodes.Size() > 0)
      {
	cout << "UpdateCoarseGrid - vertices" << endl;
        int newnv = mesh.mlbetweennodes.Size();
        loc2distvert.ChangeSize(mesh.mlbetweennodes.Size());
	/*
        for (PointIndex pi = PointIndex::BASE; pi < newnv+PointIndex::BASE; pi++)
          {
            PointIndex v1 = mesh.mlbetweennodes[pi][0];
            PointIndex v2 = mesh.mlbetweennodes[pi][1];
            if (mesh.mlbetweennodes[pi][0] != PointIndex::BASE-1)
              for (int dest = 1; dest < ntasks; dest++)
                if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
                  SetDistantPNum(dest, pi);
          }
	*/

	bool changed = true;
	while (changed)
	  {
	    changed = false;

	    // build exchange vertices
	    cnt_send = 0;
	    for (PointIndex pi : mesh.Points().Range())
	      for (int dist : GetDistantPNums(pi-PointIndex::BASE))
		cnt_send[dist-1]++;
	    TABLE<int> dest2vert(cnt_send);    
	    for (PointIndex pi : mesh.Points().Range())
	      for (int dist : GetDistantPNums(pi-PointIndex::BASE))
		dest2vert.Add (dist-1, pi);

	    
	    for (PointIndex pi = PointIndex::BASE; pi < newnv+PointIndex::BASE; pi++)
	      {
		PointIndex v1 = mesh.mlbetweennodes[pi][0];
		PointIndex v2 = mesh.mlbetweennodes[pi][1];
		if (mesh.mlbetweennodes[pi][0] != PointIndex::BASE-1)
		  // for (int dest = 1; dest < ntasks; dest++)
		  for (int dest : GetDistantPNums(v1-PointIndex::BASE))
		    if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
		      cnt_send[dest-1]++;
	      }
	    
	    TABLE<int> dest2pair(cnt_send);
	    // for (int dest = 1; dest < ntasks; dest++)
	      for (PointIndex pi = PointIndex::BASE; pi < newnv+PointIndex::BASE; pi++)
		{
		  PointIndex v1 = mesh.mlbetweennodes[pi][0];
		  PointIndex v2 = mesh.mlbetweennodes[pi][1];
		  if (mesh.mlbetweennodes[pi][0] != PointIndex::BASE-1)
		    for (int dest : GetDistantPNums(v1-PointIndex::BASE))
		      if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
			dest2pair.Add (dest-1, pi);
		}

	    cnt_send = 0;
	    int v1, v2;
	    for (PointIndex pi = PointIndex::BASE; pi < newnv+PointIndex::BASE; pi++)
	      {
		PointIndex v1 = mesh.mlbetweennodes[pi][0];
		PointIndex v2 = mesh.mlbetweennodes[pi][1];
		if (mesh.mlbetweennodes[pi][0] != PointIndex::BASE-1)
		  for (int dest : GetDistantPNums(v1-PointIndex::BASE))
		    if (IsExchangeVert(dest, v2))
		      cnt_send[dest-1]+=2;
	      }
	    
	    TABLE<int> send_verts(cnt_send);
	    
	    Array<int, PointIndex::BASE> loc2exchange(mesh.GetNV());
	    for (int dest = 1; dest < ntasks; dest++)
	      if (dest != id)
		{
		  loc2exchange = -1;
		  int cnt = 0;
		  /*
		  for (PointIndex pi : mesh.Points().Range())
		    if (IsExchangeVert(dest, pi))
		      loc2exchange[pi] = cnt++;
		  */
		  for (PointIndex pi : dest2vert[dest-1])
		    loc2exchange[pi] = cnt++;
		  
		  // for (PointIndex pi = PointIndex::BASE; pi < newnv+PointIndex::BASE; pi++)
		  for (PointIndex pi : dest2pair[dest-1])
		    {
		      PointIndex v1 = mesh.mlbetweennodes[pi][0];
		      PointIndex v2 = mesh.mlbetweennodes[pi][1];
		      if (mesh.mlbetweennodes[pi][0] != PointIndex::BASE-1)
			if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
			  {
			    send_verts.Add (dest-1, loc2exchange[v1]);
			    send_verts.Add (dest-1, loc2exchange[v2]);
			  }
		    }
		}
	    
	    TABLE<int> recv_verts(ntasks-1);
	    MyMPI_ExchangeTable (send_verts, recv_verts, MPI_TAG_MESH+9, MPI_LocalComm);
	    
	    for (int dest = 1; dest < ntasks; dest++)
	      if (dest != id)
		{
		  loc2exchange = -1;
		  int cnt = 0;
		  /*
		  for (PointIndex pi : mesh.Points().Range())
		    if (IsExchangeVert(dest, pi))
		      loc2exchange[pi] = cnt++;
		  */
		  for (PointIndex pi : dest2vert[dest-1])
		    loc2exchange[pi] = cnt++;
		  
		  FlatArray<int> recvarray = recv_verts[dest-1];
		  for (int ii = 0; ii < recvarray.Size(); ii+=2)
		    for (PointIndex pi : dest2pair[dest-1])
		      // for (PointIndex pi = PointIndex::BASE; pi < newnv+PointIndex::BASE; pi++)
		      {
			PointIndex v1 = mesh.mlbetweennodes[pi][0];
			PointIndex v2 = mesh.mlbetweennodes[pi][1];
			if (mesh.mlbetweennodes[pi][0] != PointIndex::BASE-1)
			  {
			    INDEX_2 re(recvarray[ii], recvarray[ii+1]);
			    INDEX_2 es(loc2exchange[v1], loc2exchange[v2]);
			    if (es == re && !IsExchangeVert(dest, pi))
			      {
				SetDistantPNum(dest, pi);
				changed = true;
			      }
			  }
		      }
		}
	  }
      }

    Array<int> sendarray, recvarray;
    cout << "UpdateCoarseGrid - edges" << endl;

    // static int timerv = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex vertices");
    static int timere = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex edges");
    static int timerf = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex faces");


    NgProfiler::StartTimer (timere);


    int nfa = topology . GetNFaces();
    int ned = topology . GetNEdges();
    
    // build exchange vertices
    cnt_send = 0;
    for (PointIndex pi : mesh.Points().Range())
      for (int dist : GetDistantPNums(pi-PointIndex::BASE))
	cnt_send[dist-1]++;
    TABLE<int> dest2vert(cnt_send);    
    for (PointIndex pi : mesh.Points().Range())
      for (int dist : GetDistantPNums(pi-PointIndex::BASE))
	dest2vert.Add (dist-1, pi);

    // exchange edges
    cnt_send = 0;
    int v1, v2;
    for (int edge = 1; edge <= ned; edge++)
      {
	topology.GetEdgeVertices (edge, v1, v2);
	for (int dest = 1; dest < ntasks; dest++)
	  if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
	    cnt_send[dest-1]+=1;
      }
    
    TABLE<int> dest2edge(cnt_send);
    for (int & v : cnt_send) v *= 2;
    TABLE<int> send_edges(cnt_send);

    for (int edge = 1; edge <= ned; edge++)
      {
	topology.GetEdgeVertices (edge, v1, v2);
	for (int dest = 1; dest < ntasks; dest++)
	  if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
	    dest2edge.Add (dest-1, edge);
      }


    Array<int, PointIndex::BASE> loc2exchange(mesh.GetNV());
    for (int dest = 1; dest < ntasks; dest++)
      {
        loc2exchange = -1;
        int cnt = 0;
        for (PointIndex pi : dest2vert[dest-1])
	  loc2exchange[pi] = cnt++;

	for (int edge : dest2edge[dest-1])
          {
            topology.GetEdgeVertices (edge, v1, v2);
            if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
              {
                send_edges.Add (dest-1, loc2exchange[v1]);
                send_edges.Add (dest-1, loc2exchange[v2]);
              }
          }
      }

    cout << "UpdateCoarseGrid - edges mpi-exchange" << endl;
    TABLE<int> recv_edges(ntasks-1);
    MyMPI_ExchangeTable (send_edges, recv_edges, MPI_TAG_MESH+9, MPI_LocalComm);
    cout << "UpdateCoarseGrid - edges mpi-exchange done" << endl;

    /*
    for (int dest = 1; dest < ntasks; dest++)
      {
	auto ex2loc = dest2vert[dest-1];
	FlatArray<int> recvarray = recv_edges[dest-1];
        for (int ii = 0; ii < recvarray.Size(); ii+=2)
	  for (int edge : dest2edge[dest-1])
	    {
	      topology.GetEdgeVertices (edge, v1, v2);
	      INDEX_2 re(ex2loc[recvarray[ii]], 
			 ex2loc[recvarray[ii+1]]);
	      INDEX_2 es(v1, v2);
	      if (es == re)
		SetDistantEdgeNum(dest, edge);
	    }
      }
    */

    for (int dest = 1; dest < ntasks; dest++)
      {
	auto ex2loc = dest2vert[dest-1];
	if (ex2loc.Size() == 0) continue;

	INDEX_2_CLOSED_HASHTABLE<int> vert2edge(2*dest2edge[dest-1].Size()+10); 
	for (int edge : dest2edge[dest-1])
	  {
	    topology.GetEdgeVertices (edge, v1, v2);
	    vert2edge.Set(INDEX_2(v1,v2), edge);
	  }

	FlatArray<int> recvarray = recv_edges[dest-1];
        for (int ii = 0; ii < recvarray.Size(); ii+=2)
	  {
	    INDEX_2 re(ex2loc[recvarray[ii]], 
		       ex2loc[recvarray[ii+1]]);
	    if (vert2edge.Used(re))
	      SetDistantEdgeNum(dest, vert2edge.Get(re));
	  }
      }



    NgProfiler::StopTimer (timere);

    // MPI_Barrier (MPI_LocalComm);

    cout << "UpdateCoarseGrid - faces" << endl;
    if (mesh.GetDimension() == 3)
      {
	NgProfiler::StartTimer (timerf);
	Array<int> verts;

	// exchange faces
	cnt_send = 0;
	for (int face = 1; face <= nfa; face++)
	  {
	    topology.GetFaceVertices (face, verts);
	    for (int dest = 1; dest < ntasks; dest++)
	      if (dest != id)
		if (IsExchangeVert (dest, verts[0]) && 
		    IsExchangeVert (dest, verts[1]) &&
		    IsExchangeVert (dest, verts[2]))
		  cnt_send[dest-1]++;
	  }
	
	TABLE<int> dest2face(cnt_send);
	for (int face = 1; face <= nfa; face++)
	  {
	    topology.GetFaceVertices (face, verts);
	    for (int dest = 1; dest < ntasks; dest++)
	      if (dest != id)
		if (IsExchangeVert (dest, verts[0]) && 
		    IsExchangeVert (dest, verts[1]) &&
		    IsExchangeVert (dest, verts[2]))
		  dest2face.Add(dest-1, face);
	  }

	for (int & c : cnt_send) c*=3;
	TABLE<int> send_faces(cnt_send);
	Array<int, PointIndex::BASE> loc2exchange(mesh.GetNV());
	for (int dest = 1; dest < ntasks; dest++)
	  if (dest != id)
	    {
	      /*
	      loc2exchange = -1;
	      int cnt = 0;
	      for (PointIndex pi : mesh.Points().Range())
		if (IsExchangeVert(dest, pi))
		  loc2exchange[pi] = cnt++;
	      */
	      if (dest2vert[dest-1].Size() == 0) continue;

	      loc2exchange = -1;
	      int cnt = 0;
	      for (PointIndex pi : dest2vert[dest-1])
		loc2exchange[pi] = cnt++;
	      
	      for (int face : dest2face[dest-1])
		{
		  topology.GetFaceVertices (face, verts);
		  if (IsExchangeVert (dest, verts[0]) && 
		      IsExchangeVert (dest, verts[1]) &&
		      IsExchangeVert (dest, verts[2]))
		    {
		      send_faces.Add (dest-1, loc2exchange[verts[0]]);
		      send_faces.Add (dest-1, loc2exchange[verts[1]]);
		      send_faces.Add (dest-1, loc2exchange[verts[2]]);
		    }
		}
	    }
	
	cout << "UpdateCoarseGrid - faces mpi-exchange" << endl;
	TABLE<int> recv_faces(ntasks-1);
	MyMPI_ExchangeTable (send_faces, recv_faces, MPI_TAG_MESH+9, MPI_LocalComm);
	cout << "UpdateCoarseGrid - faces mpi-exchange done" << endl;

	/*
	for (int dest = 1; dest < ntasks; dest++)
	  if (dest != id)
	    {
	      loc2exchange = -1;
	      int cnt = 0;
	      for (PointIndex pi : dest2vert[dest-1])
		loc2exchange[pi] = cnt++;
	      
	      FlatArray<int> recvarray = recv_faces[dest-1];
	      for (int ii = 0; ii < recvarray.Size(); ii+=3)
		for (int face : dest2face[dest-1])
		  {
		    topology.GetFaceVertices (face, verts);
		    INDEX_3 re(recvarray[ii], recvarray[ii+1], recvarray[ii+2]);
		    INDEX_3 es(loc2exchange[verts[0]], loc2exchange[verts[1]], loc2exchange[verts[2]]);
		    if (es == re)
		      SetDistantFaceNum(dest, face);
		  }
	    }
	*/

	
	for (int dest = 1; dest < ntasks; dest++)
	  {
	    auto ex2loc = dest2vert[dest-1];
	    if (ex2loc.Size() == 0) continue;
	    
	    INDEX_3_CLOSED_HASHTABLE<int> vert2face(2*dest2face[dest-1].Size()+10); 
	    for (int face : dest2face[dest-1])
	      {
		topology.GetFaceVertices (face, verts);
		vert2face.Set(INDEX_3(verts[0], verts[1], verts[2]), face);
	      }
	    
	    FlatArray<int> recvarray = recv_faces[dest-1];
	    for (int ii = 0; ii < recvarray.Size(); ii+=3)
	      {
		INDEX_3 re(ex2loc[recvarray[ii]], 
			   ex2loc[recvarray[ii+1]],
			   ex2loc[recvarray[ii+2]]);
		if (vert2face.Used(re))
		  SetDistantFaceNum(dest, vert2face.Get(re));
	      }
	  }
	



	

	/*
	  Array<int,1> glob2loc;

	int maxface = 0;
	for (int face = 1; face <= nfa; face++)
	  maxface = max (maxface, GetGlobalFaceNum (face));
	
	// glob2loc.SetSize (nfaglob);
	glob2loc.SetSize (maxface);
	glob2loc = -1;
	
	for (int loc = 1; loc <= nfa; loc++)
	  glob2loc[GetGlobalFaceNum(loc)] = loc;
	
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
		    send_faces.Add (dest-1, GetGlobalFaceNum(face));
		    send_faces.Add (dest-1, face);
		  }
	      }
	  }
	TABLE<int> recv_faces(ntasks-1);
	MyMPI_ExchangeTable (send_faces, recv_faces, MPI_TAG_MESH+8, MPI_LocalComm);
	
	for (int sender = 1; sender < ntasks; sender ++)
	  if (id != sender)
	    {
	      FlatArray<int> recvarray = recv_faces[sender-1];
	      
	      for (int ii = 0; ii < recvarray.Size(); )
		{ 
		  int globf = recvarray[ii++];
		  int distf = recvarray[ii++];
		  
		  if (globf <= maxface)
		    {
		      int locf = glob2loc[globf];
		      if (locf != -1)
			SetDistantFaceNum (sender, locf);
		    }
		}
	    } 
	*/	
	
	NgProfiler::StopTimer (timerf);
      }
    cout << "UpdateCoarseGrid - done" << endl;
    
    is_updated = true;
  }
}


#endif
