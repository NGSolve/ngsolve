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
    if (is_updated) return;

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



    Array<int> sendarray, recvarray;


    static int timerv = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex vertices");
    static int timere = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex edges");
    static int timerf = NgProfiler::CreateTimer ("UpdateCoarseGrid - ex faces");


    Array<int> cnt_send(ntasks-1);

    NgProfiler::StartTimer (timere);


    const MeshTopology & topology = mesh.GetTopology();
    int nfa = topology . GetNFaces();
    int ned = topology . GetNEdges();
    

    // exchange edges
    cnt_send = 0;
    int v1, v2;
    for (int edge = 1; edge <= ned; edge++)
      {
	topology.GetEdgeVertices (edge, v1, v2);
	for (int dest = 1; dest < ntasks; dest++)
	  if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
	    cnt_send[dest-1]+=2;
      }
    
    TABLE<int> send_edges(cnt_send);
    INDEX_2_HASHTABLE<int> gv2e(2*ned);

    for (int edge = 1; edge <= ned; edge++)
      {
	topology.GetEdgeVertices (edge, v1, v2);
	INDEX_2 es(GetGlobalPNum(v1), GetGlobalPNum(v2));
	es.Sort();

	gv2e.Set (es, edge);

	for (int dest = 1; dest < ntasks; dest++)
	  if (IsExchangeVert (dest, v1) && IsExchangeVert (dest, v2))
	    {
	      send_edges.Add (dest-1, es[0]);
	      send_edges.Add (dest-1, es[1]);
	    }
      }

    TABLE<int> recv_edges(ntasks-1);
    MyMPI_ExchangeTable (send_edges, recv_edges, MPI_TAG_MESH+9, MPI_LocalComm);

    for (int sender = 1; sender < ntasks; sender ++)
      if (id != sender)
	{
	  FlatArray<int> recvarray = recv_edges[sender-1];
	  for (int ii = 0; ii < recvarray.Size(); ii+=2)
	    { 
	      INDEX_2 gv12 (recvarray[ii],recvarray[ii+1]);
	      if (gv2e.Used (gv12))
		SetDistantEdgeNum (sender, gv2e.Get(gv12));
	    }
	}

 
    NgProfiler::StopTimer (timere);

    // MPI_Barrier (MPI_LocalComm);


    if (mesh.GetDimension() == 3)
      {
	NgProfiler::StartTimer (timerf);

	// exchange faces
	cnt_send = 0;
	Array<int> verts;
	for (int face = 1; face <= nfa; face++)
	  {
	    topology.GetFaceVertices (face, verts);
	    for (int dest = 1; dest < ntasks; dest++)
	      if (IsExchangeVert (dest, verts[0]) && 
		  IsExchangeVert (dest, verts[1]) &&
		  IsExchangeVert (dest, verts[2]))
		cnt_send[dest-1]+=3;
	  }
    
	TABLE<int> send_faces(cnt_send);
	INDEX_3_HASHTABLE<int> gv2f(2*nfa);

	for (int face = 1; face <= nfa; face++)
	  {
	    topology.GetFaceVertices (face, verts);
	    INDEX_3 fs (GetGlobalPNum(verts[0]),
			GetGlobalPNum(verts[1]),
			GetGlobalPNum(verts[2]));
	    fs.Sort();

	    gv2f.Set (fs, face);

	    for (int dest = 1; dest < ntasks; dest++)
	      if (IsExchangeVert (dest, verts[0]) && 
		  IsExchangeVert (dest, verts[1]) &&
		  IsExchangeVert (dest, verts[2]))
		{
		  send_faces.Add (dest-1, fs[0]);
		  send_faces.Add (dest-1, fs[1]);
		  send_faces.Add (dest-1, fs[2]);
		}
	  }
	
	TABLE<int> recv_faces(ntasks-1);
	MyMPI_ExchangeTable (send_faces, recv_faces, MPI_TAG_MESH+9, MPI_LocalComm);
	
	for (int sender = 1; sender < ntasks; sender ++)
	  if (id != sender)
	    {
	      FlatArray<int> recvarray = recv_faces[sender-1];
	      for (int ii = 0; ii < recvarray.Size(); ii+=3)
		{ 
		  INDEX_3 gv123 (recvarray[ii],recvarray[ii+1],recvarray[ii+2]);
		  if (gv2f.Used (gv123))
		    SetDistantFaceNum (sender, gv2f.Get(gv123));
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
    
    is_updated = true;
  }
}


#endif
