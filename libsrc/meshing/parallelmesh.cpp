#ifdef PARALLEL

#include <meshing.hpp>
#include "paralleltop.hpp"

// #define METIS4


#ifdef METIS
namespace metis {
  extern "C" {

#include <metis.h>

#if METIS_VER_MAJOR >= 5
#define METIS5
    typedef idx_t idxtype;   
#else
#define METIS4
    typedef idxtype idx_t;  
#endif
  } 
}

using namespace metis;
#endif

namespace netgen
{

  template <>
  inline MPI_Datatype MyGetMPIType<PointIndex> ( )
  { return MPI_INT; }


  void Mesh :: SendRecvMesh ()
  {
    if (id == 0)
      PrintMessage (1, "Send/Receive mesh");

    {
      // distribute global information
      int nelglob, nseglob, nvglob;
      if (id == 0)
	{
	  paralleltop -> SetNV (GetNV());
	  paralleltop -> SetNE (GetNE());
	  paralleltop -> SetNSegm (GetNSeg());
	  paralleltop -> SetNSE (GetNSE());

	  nelglob = GetNE();
	  nseglob = GetNSE();
	  nvglob = GetNV();
	}
      
      MyMPI_Bcast (dimension);
      MyMPI_Bcast (nelglob);
      MyMPI_Bcast (nseglob);
      MyMPI_Bcast (nvglob);
    }



    {
      // distribute number of local elements
      if (id == 0)
	{
	  Array<int> num_els_on_proc(ntasks);
	  num_els_on_proc = 0;
	  for (ElementIndex ei = 0; ei < GetNE(); ei++)
	    num_els_on_proc[(*this)[ei].GetPartition()]++;
	  
	  MPI_Scatter (&num_els_on_proc[0], 1, MPI_INT,
		       MPI_IN_PLACE, -1, MPI_INT, 0, MPI_COMM_WORLD);
	}
      else
	{
	  int nelloc;
	  MPI_Scatter (NULL, 0, MPI_INT,
		       &nelloc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	  
	  paralleltop -> SetNE (nelloc);
	}
    }

    if (id == 0)
      SendMesh ();
    else
      ReceiveParallelMesh();

    paralleltop -> UpdateCoarseGrid();
  }







  void Mesh :: SendMesh () const   
  {
    Array<MPI_Request> sendrequests;

    PrintMessage ( 3, "Sending vertices");


    Array<int> num_els_on_proc(ntasks);
    num_els_on_proc = 0;
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      num_els_on_proc[(*this)[ei].GetPartition()]++;

    TABLE<ElementIndex> els_of_proc (num_els_on_proc);
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      els_of_proc.Add ( (*this)[ei].GetPartition(), ei);


    Array<int> num_sels_on_proc(ntasks);
    num_sels_on_proc = 0;
    for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
      num_sels_on_proc[(*this)[ei].GetPartition()]++;

    TABLE<SurfaceElementIndex> sels_of_proc (num_sels_on_proc);
    for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
      sels_of_proc.Add ( (*this)[ei].GetPartition(), ei);

    Array<int> num_segs_on_proc(ntasks);
    num_segs_on_proc = 0;
    for (SegmentIndex ei = 0; ei < GetNSeg(); ei++)
      num_segs_on_proc[(*this)[ei].GetPartition()]++;

    TABLE<SegmentIndex> segs_of_proc (num_segs_on_proc);
    for (SegmentIndex ei = 0; ei < GetNSeg(); ei++)
      segs_of_proc.Add ( (*this)[ei].GetPartition(), ei);




    Array<int, PointIndex::BASE> vert_flag (GetNV());
    Array<int, PointIndex::BASE> num_procs_on_vert (GetNV());
    Array<int> num_verts_on_proc (ntasks);

    num_verts_on_proc = 0;
    num_procs_on_vert = 0;
    vert_flag = -1;


    for (int dest = 1; dest < ntasks; dest++)
      {
	FlatArray<ElementIndex> els = els_of_proc[dest];
	for (int hi = 0; hi < els.Size(); hi++)
	  {
	    const Element & el = (*this) [ els[hi] ];
	    for (int i = 0; i < el.GetNP(); i++)
	      {
		PointIndex epi = el[i]; 
		if (vert_flag[epi] < dest)
		  {
		    vert_flag[epi] = dest;

		    num_verts_on_proc[dest]++;
		    num_procs_on_vert[epi]++;

		    paralleltop -> SetDistantPNum (dest, epi); 
		  }
	      }
	  }


	FlatArray<SurfaceElementIndex> sels = sels_of_proc[dest];
	for (int hi = 0; hi < sels.Size(); hi++)
	  {
	    const Element2d & el = (*this) [ sels[hi] ];
	    for (int i = 0; i < el.GetNP(); i++)
	      {
		PointIndex epi = el[i]; 
		if (vert_flag[epi] < dest)
		  {
		    vert_flag[epi] = dest;

		    num_verts_on_proc[dest]++;
		    num_procs_on_vert[epi]++;

		    paralleltop -> SetDistantPNum (dest, epi);
		  }
	      }
	  }

	FlatArray<SegmentIndex> segs = segs_of_proc[dest];
	for (int hi = 0; hi < segs.Size(); hi++)
	  {
	    const Segment & el = (*this) [segs[hi]];
	    for (int i = 0; i < 2; i++)
	      {
		PointIndex epi = el[i]; 
		if (vert_flag[epi] < dest)
		  {
		    vert_flag[epi] = dest;

		    num_verts_on_proc[dest]++;
		    num_procs_on_vert[epi]++;

		    paralleltop -> SetDistantPNum (dest, epi);
		  }
	      }
	  }
      }

    TABLE<PointIndex> verts_of_proc (num_verts_on_proc);
    TABLE<int, PointIndex::BASE> procs_of_vert (num_procs_on_vert);
    TABLE<int, PointIndex::BASE> loc_num_of_vert (num_procs_on_vert);

    vert_flag = -1;
    for (int dest = 1; dest < ntasks; dest++)
      {
	FlatArray<ElementIndex> els = els_of_proc[dest];
	for (int hi = 0; hi < els.Size(); hi++)
	  {
	    const Element & el = (*this) [ els[hi] ];
	    for (int i = 0; i < el.GetNP(); i++)
	      {
		PointIndex epi = el[i];
		if (vert_flag[epi] < dest)
		  {
		    vert_flag[epi] = dest;
		    procs_of_vert.Add (epi, dest);
		  }
	      }
	  }

	FlatArray<SurfaceElementIndex> sels = sels_of_proc[dest];
	for (int hi = 0; hi < sels.Size(); hi++)
	  {
	    const Element2d & el = (*this) [ sels[hi] ];
	    for (int i = 0; i < el.GetNP(); i++)
	      {
		PointIndex epi = el[i];
		if (vert_flag[epi] < dest)
		  {
		    vert_flag[epi] = dest;
		    procs_of_vert.Add (epi, dest);
		  }
	      }
	  }

	FlatArray<SegmentIndex> segs = segs_of_proc[dest];
	for (int hi = 0; hi < segs.Size(); hi++)
	  {
	    const Segment & el = (*this) [ segs[hi] ];
	    for (int i = 0; i < 2; i++)
	      {
		PointIndex epi = el[i];
		if (vert_flag[epi] < dest)
		  {
		    vert_flag[epi] = dest;
		    procs_of_vert.Add (epi, dest);
		  }
	      }
	  }

      }

    for (int vert = 1; vert <= GetNP(); vert++ )
      {
	FlatArray<int> procs = procs_of_vert[vert];
	for (int j = 0; j < procs.Size(); j++)
	  {
	    int dest = procs[j];
	    verts_of_proc.Add (dest, vert);
	    loc_num_of_vert.Add (vert, verts_of_proc[dest].Size());
	  }
      }

    for (int dest = 1; dest < ntasks; dest++)
      {
	FlatArray<PointIndex> verts = verts_of_proc[dest];
	sendrequests.Append (MyMPI_ISend (verts, dest, MPI_TAG_MESH+1));

	MPI_Datatype mptype = MeshPoint::MyGetMPIType();

	int numv = verts.Size();

	MPI_Datatype newtype;
	Array<int> blocklen (numv);  
	blocklen = 1;
	
	MPI_Type_indexed (numv, &blocklen[0], 
			  reinterpret_cast<int*> (&verts[0]), 
			  mptype, &newtype);
	MPI_Type_commit (&newtype);

	MPI_Request request;
	MPI_Isend( &points[0], 1, newtype, dest, MPI_TAG_MESH+1, MPI_COMM_WORLD, &request);
	sendrequests.Append (request);
      }

    Array<int> num_distpnums(ntasks);
    num_distpnums = 0;

    for (int vert = 1; vert <= GetNP(); vert++)
      {
	FlatArray<int> procs = procs_of_vert[vert];
	for (int j = 0; j < procs.Size(); j++)
	  num_distpnums[procs[j]] += 3 * (procs.Size()-1);
      }

    TABLE<int> distpnums (num_distpnums);

    for (int vert = 1; vert <= GetNP(); vert++)
      {
	FlatArray<int> procs = procs_of_vert[vert];
	for (int j = 0; j < procs.Size(); j++)
	  for (int k = 0; k < procs.Size(); k++)
	    if (j != k)
	      {
		distpnums.Add (procs[j], loc_num_of_vert[vert][j]);
		distpnums.Add (procs[j], procs_of_vert[vert][k]);
		distpnums.Add (procs[j], loc_num_of_vert[vert][k]);
	      }
      }
    
    for ( int dest = 1; dest < ntasks; dest ++ )
      sendrequests.Append (MyMPI_ISend (distpnums[dest], dest, MPI_TAG_MESH+1));



    PrintMessage ( 3, "Sending elements" );

    Array<int> elarraysize (ntasks);
    elarraysize = 0;
    for ( int ei = 1; ei <= GetNE(); ei++)
      {
	const Element & el = VolumeElement (ei);
	int dest = el.GetPartition();
	elarraysize[dest] += 3 + el.GetNP();
      }

    TABLE<int> elementarrays(elarraysize);

    for (int ei = 1; ei <= GetNE(); ei++)
      {
	const Element & el = VolumeElement (ei);
	int dest = el.GetPartition();

	elementarrays.Add (dest, ei);
	elementarrays.Add (dest, el.GetIndex());
	elementarrays.Add (dest, el.GetNP());
	for (int i = 0; i < el.GetNP(); i++)
	  elementarrays.Add (dest, el[i]);
      }

    for (int dest = 1; dest < ntasks; dest ++ )
      sendrequests.Append (MyMPI_ISend (elementarrays[dest], dest, MPI_TAG_MESH+2));


    PrintMessage ( 3, "Sending Face Descriptors" );

    Array<double> fddata (6 * GetNFD());
    for (int fdi = 1; fdi <= GetNFD(); fdi++)
      {
	fddata[6*fdi-6] = GetFaceDescriptor(fdi).SurfNr();
	fddata[6*fdi-5] = GetFaceDescriptor(fdi).DomainIn();	
	fddata[6*fdi-4] = GetFaceDescriptor(fdi).DomainOut();
	fddata[6*fdi-3] = GetFaceDescriptor(fdi).BCProperty();
	fddata[6*fdi-2] = GetFaceDescriptor(fdi).domin_singular;
	fddata[6*fdi-1] = GetFaceDescriptor(fdi).domout_singular;
	
      }
    for (int dest = 1; dest < ntasks; dest++)
      sendrequests.Append (MyMPI_ISend (fddata, dest, MPI_TAG_MESH+3));

    PrintMessage ( 3, "Sending Surface elements" );
  
    Array <int> nlocsel(ntasks), bufsize(ntasks); 
    nlocsel = 0;
    bufsize = 1;

    for (int sei = 1; sei <= GetNSE(); sei++ )
      {
	const Element2d & sel = SurfaceElement (sei);
	int dest = sel.GetPartition();
	nlocsel[dest] ++;
	bufsize[dest] += 4 + 2*sel.GetNP();
      }
    
    TABLE<int> selbuf(bufsize);

    for (int dest = 1; dest < ntasks; dest++ )
      selbuf.Add (dest, nlocsel[dest]);

    for (int sei = 1; sei <= GetNSE(); sei ++ )
      {
	const Element2d & sel = SurfaceElement (sei);
	int dest = sel.GetPartition();

	selbuf.Add (dest, sei);
	selbuf.Add (dest, sel.GetIndex());
	// selbuf.Add (dest, 0);
	selbuf.Add (dest, sel.GetNP());
	
	for ( int ii = 1; ii <= sel.GetNP(); ii++)
	  {
	    selbuf.Add (dest, sel.PNum(ii));
	    selbuf.Add (dest, sel.GeomInfoPi(ii).trignum);
	  }
      }

    for (int dest = 1; dest < ntasks; dest++)
      sendrequests.Append (MyMPI_ISend(selbuf[dest], dest, MPI_TAG_MESH+4));


    PrintMessage ( 3, "Sending Edge Segments");

    
    Array <int> nloc(ntasks); 
    nloc = 0;
    bufsize = 1;
      
    for (int i = 1; i <= GetNSeg(); i++ )
      {
	const Segment & seg = LineSegment (i);
	int dest = seg.GetPartition();
	nloc[dest] ++;
	bufsize[dest] += 14;
      }
    
    TABLE<double> segm_buf(bufsize);
    
    /*
    for (int dest = 1; dest < ntasks; dest++ )
      segm_buf.Add (dest, nloc[dest]);
    */
    for (int i = 1; i <= GetNSeg(); i++ )
      {
	const Segment & seg = LineSegment (i);
	int dest = seg.GetPartition();
	
	segm_buf.Add (dest, i);
	segm_buf.Add (dest, seg.si);
	segm_buf.Add (dest, seg.pnums[0]);
	segm_buf.Add (dest, seg.pnums[1]);
	segm_buf.Add (dest, seg.geominfo[0].trignum);
	segm_buf.Add (dest, seg.geominfo[1].trignum);
	segm_buf.Add (dest, seg.surfnr1);
	segm_buf.Add (dest, seg.surfnr2);
	segm_buf.Add (dest, seg.edgenr);
	segm_buf.Add (dest, seg.epgeominfo[0].dist);
	segm_buf.Add (dest, seg.epgeominfo[1].edgenr);
	segm_buf.Add (dest, seg.epgeominfo[1].dist);
	segm_buf.Add (dest, seg.singedge_right);
	segm_buf.Add (dest, seg.singedge_left);
      }
    
    for (int dest = 1; dest < ntasks; dest++)
      sendrequests.Append (MyMPI_ISend(segm_buf[dest], dest, MPI_TAG_MESH+5));

    /*
    Array <int> nlocseg(ntasks), segi(ntasks);
    for ( int i = 0; i < ntasks; i++)
      {
	nlocseg[i] = 0;
	bufsize[i] = 0;
	segi[i] = 0;
      }

    for (int segi = 1; segi <= GetNSeg(); segi ++ )
      {
	Array<int> volels;
	const MeshTopology & topol = GetTopology();
	topol . GetSegmentSurfaceElements ( segi, volels );
        for (int j = 0; j < volels.Size(); j++)
          {
            int ei = volels[j];
            if ( ei > 0 && ei <= GetNSE() )
              {
                const Element2d & el = SurfaceElement (ei);
                int dest = el.GetPartition();
		nlocseg[dest] ++;
		bufsize[dest] += 14;
              }
	  }
      }

    TABLE<double> segmbuf(bufsize);

    for ( int ls=1; ls <= GetNSeg(); ls++)
      {
	Array<int> volels;
	GetTopology().GetSegmentSurfaceElements ( ls, volels );
	const Segment & seg = LineSegment (ls);

        for (int j = 0; j < volels.Size(); j++)
          {
            int ei = volels[j];
            if ( ei > 0 && ei <= GetNSE() )
              {
                const Element2d & el = SurfaceElement (ei);
                int dest = el.GetPartition();

		if ( dest > 0 )
		  {
		    segmbuf.Add (dest, ls);
		    segmbuf.Add (dest, seg.si);
		    segmbuf.Add (dest, seg.pnums[0]);
		    segmbuf.Add (dest, seg.pnums[1]);
		    segmbuf.Add (dest, seg.geominfo[0].trignum);
		    segmbuf.Add (dest, seg.geominfo[1].trignum);
		    segmbuf.Add (dest, seg.surfnr1);
		    segmbuf.Add (dest, seg.surfnr2);
		    segmbuf.Add (dest, seg.edgenr);
		    segmbuf.Add (dest, seg.epgeominfo[0].dist);
		    segmbuf.Add (dest, seg.epgeominfo[1].edgenr);
		    segmbuf.Add (dest, seg.epgeominfo[1].dist);
		    segmbuf.Add (dest, seg.singedge_right);
		    segmbuf.Add (dest, seg.singedge_left);
		    segi[dest] += 14;
		  }
		// paralleltop -> SetDistantSegm ( dest, ls, int ( segi[dest] / 14 ) );
	      }
	  }
      }

    for ( int dest = 1; dest < ntasks; dest++)
      sendrequests.Append (MyMPI_ISend(segmbuf[dest], dest, MPI_TAG_MESH+5));
    */

    PrintMessage ( 3, "now wait ...");

    MPI_Waitall (sendrequests.Size(), &sendrequests[0], MPI_STATUS_IGNORE);

    MPI_Barrier(MPI_COMM_WORLD);
  }








  // slaves receive the mesh from the master
  void Mesh :: ReceiveParallelMesh ( )
  {
    int timer = NgProfiler::CreateTimer ("ReceiveParallelMesh");
    int timer_pts = NgProfiler::CreateTimer ("Receive points");
    int timer_els = NgProfiler::CreateTimer ("Receive elements");
    int timer_sels = NgProfiler::CreateTimer ("Receive surface elements");
    NgProfiler::RegionTimer reg(timer);


    // string st;

    // receive vertices
    NgProfiler::StartTimer (timer_pts);

    Array<int> verts;
    MyMPI_Recv (verts, 0, MPI_TAG_MESH+1);

    int numvert = verts.Size();
    paralleltop -> SetNV (numvert);
    
    // INDEX_CLOSED_HASHTABLE<int> glob2loc_vert_ht (3*numvert+1);
    INDEX_HASHTABLE<int> glob2loc_vert_ht (3*numvert+1);

    for (int vert = 0; vert < numvert; vert++)
      {
	int globvert = verts[vert];
	paralleltop->SetLoc2Glob_Vert ( vert+1, globvert  );
	glob2loc_vert_ht.Set (globvert, vert+1);
      }
    
    for (int i = 0; i < numvert; i++)
      AddPoint (netgen::Point<3> (0,0,0));
    
    MPI_Datatype mptype = MeshPoint::MyGetMPIType();
    MPI_Status status;
    MPI_Recv( &points[1], numvert, mptype, 0, MPI_TAG_MESH+1, MPI_COMM_WORLD, &status);
    
    Array<int> dist_pnums;
    MyMPI_Recv (dist_pnums, 0, MPI_TAG_MESH+1);
    
    for (int hi = 0; hi < dist_pnums.Size(); hi += 3)
      paralleltop ->
	SetDistantPNum (dist_pnums[hi+1], dist_pnums[hi]); // , dist_pnums[hi+2]);
    
    NgProfiler::StopTimer (timer_pts);
    *testout << "got " << numvert << " vertices" << endl;

    {
      Element el;
      
      Array<int> elarray;
      MyMPI_Recv (elarray, 0, MPI_TAG_MESH+2);
      
      NgProfiler::RegionTimer reg(timer_els);

      for (int ind = 0, elnum = 1; ind < elarray.Size(); elnum++)
	{
	  paralleltop->SetLoc2Glob_VolEl ( elnum,  elarray[ind++]);
	  
	  el.SetIndex(elarray[ind++]);
	  el.SetNP(elarray[ind++]);
	  
	  for ( int j = 0; j < el.GetNP(); j++)
	    el[j] = glob2loc_vert_ht.Get (elarray[ind++]); 
	  
	  AddVolumeElement (el);
	}
    }

    {
      Array<double> fddata;
      MyMPI_Recv (fddata, 0, MPI_TAG_MESH+3);
      for (int i = 0; i < fddata.Size(); i += 6)
	{
	  int faceind = AddFaceDescriptor 
	    (FaceDescriptor(int(fddata[i]), int(fddata[i+1]), int(fddata[i+2]), 0));
	  GetFaceDescriptor(faceind).SetBCProperty (int(fddata[i+3]));
	  GetFaceDescriptor(faceind).domin_singular = fddata[i+4];
	  GetFaceDescriptor(faceind).domout_singular = fddata[i+5];
	}
    }

    {
      NgProfiler::RegionTimer reg(timer_sels);
      Array<int> selbuf;

      MyMPI_Recv ( selbuf, 0, MPI_TAG_MESH+4);
      
      int ii = 0;
      int sel = 0;

      int nlocsel = selbuf[ii++];
      paralleltop -> SetNSE ( nlocsel );
      
      while (ii < selbuf.Size()-1)
	{
	  int globsel = selbuf[ii++];
	  int faceind = selbuf[ii++];
	  //bool isghost = selbuf[ii++];
	  int nep = selbuf[ii++];
	  Element2d tri(nep);
	  tri.SetIndex(faceind);
	  for(int j = 1; j <= nep; j++)
	    {
	      tri.PNum(j) = glob2loc_vert_ht.Get (selbuf[ii++]);
	      tri.GeomInfoPi(j).trignum = selbuf[ii++];
	    }
	  
	  paralleltop->SetLoc2Glob_SurfEl ( sel+1, globsel );
	  AddSurfaceElement (tri);
	  sel ++;
	}
    }
    


    {
      Array<double> segmbuf;
      MyMPI_Recv ( segmbuf, 0, MPI_TAG_MESH+5);

      Segment seg;
      int globsegi;
      int ii = 0;
      int segi = 1;
      int nsegloc = int ( segmbuf.Size() / 14 ) ;
      paralleltop -> SetNSegm ( nsegloc );

      while ( ii < segmbuf.Size() )
	{
	  globsegi = int (segmbuf[ii++]);
	  seg.si = int (segmbuf[ii++]);
	  
	  seg.pnums[0] = glob2loc_vert_ht.Get (int(segmbuf[ii++]));
	  seg.pnums[1] = glob2loc_vert_ht.Get (int(segmbuf[ii++]));
	  seg.geominfo[0].trignum = int( segmbuf[ii++] );
	  seg.geominfo[1].trignum = int ( segmbuf[ii++]);
	  seg.surfnr1 = int ( segmbuf[ii++]);
	  seg.surfnr2 = int ( segmbuf[ii++]);
	  seg.edgenr = int ( segmbuf[ii++]);
	  seg.epgeominfo[0].dist = segmbuf[ii++];
	  seg.epgeominfo[1].edgenr = int (segmbuf[ii++]);
	  seg.epgeominfo[1].dist = segmbuf[ii++];
	  
	  seg.singedge_left = segmbuf[ii++];
	  seg.singedge_right = segmbuf[ii++];
	  
	  seg.epgeominfo[0].edgenr = seg.epgeominfo[1].edgenr;
	  
	  seg.domin = seg.surfnr1;
	  seg.domout = seg.surfnr2;
	  if ( seg.pnums[0] >0 && seg.pnums[1] > 0 )
	    {
	      paralleltop-> SetLoc2Glob_Segm ( segi,  globsegi );
	      
	      AddSegment (seg);
	      segi++;
	    }
	  
	}
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    int timerloc = NgProfiler::CreateTimer ("Update local mesh");
    int timerloc2 = NgProfiler::CreateTimer ("CalcSurfacesOfNode");

    NgProfiler::RegionTimer regloc(timerloc);
    stringstream str;
    str << "p" << id << ": got " << GetNE() << " elements and " 
	 << GetNSE() << " surface elements";
    cout << str.str() << endl;
    // PrintMessage (2, "Got ", GetNE(), " elements and ", GetNSE(), " surface elements");
    // PrintMessage (2, "Got ", GetNSE(), " surface elements");

    NgProfiler::StartTimer (timerloc2);

    CalcSurfacesOfNode ();

    NgProfiler::StopTimer (timerloc2);

    topology -> Update();
    clusters -> Update();

    // paralleltop -> UpdateCoarseGrid();
    
    SetNextMajorTimeStamp();
  }
  


  
  
  // distribute the mesh to the slave processors
  // call it only for the master !
  void Mesh :: Distribute ()
  {
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if (id != 0 || ntasks == 1 ) return;

#ifdef METIS
    ParallelMetis ();
#else
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      (*this)[ei].SetPartition(ntasks * ei/GetNE() + 1);
#endif

    /*
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      *testout << "el(" << ei << ") is in part " << (*this)[ei].GetPartition() << endl;
    for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
      *testout << "sel(" << int(ei) << ") is in part " << (*this)[ei].GetPartition() << endl;
      */
    
    // MyMPI_SendCmd ("mesh");
    SendRecvMesh (); 
  }
  

#ifdef METIS5
  void Mesh :: ParallelMetis ( )  
  {
    PrintMessage (3, "call metis 5 ...");

    int timer = NgProfiler::CreateTimer ("Mesh::Partition");
    NgProfiler::RegionTimer reg(timer);

    idx_t ne = GetNE() + GetNSE() + GetNSeg();
    idx_t nn = GetNP();

    Array<idx_t> eptr, eind;
    for (int i = 0; i < GetNE(); i++)
      {
	eptr.Append (eind.Size());
	const Element & el = VolumeElement(i+1);
	for (int j = 0; j < el.GetNP(); j++)
	  eind.Append (el[j]-1);
      }
    for (int i = 0; i < GetNSE(); i++)
      {
	eptr.Append (eind.Size());
	const Element2d & el = SurfaceElement(i+1);
	for (int j = 0; j < el.GetNP(); j++)
	  eind.Append (el[j]-1);
      }
    for (int i = 0; i < GetNSeg(); i++)
      {
	eptr.Append (eind.Size());
	const Segment & el = LineSegment(i+1);
	eind.Append (el[0]-1);
	eind.Append (el[1]-1);
      }
    eptr.Append (eind.Size());
    Array<idx_t> epart(ne), npart(nn);

    idxtype nparts = ntasks-1;
    idxtype edgecut;

    idxtype ncommon = 3;
    METIS_PartMeshDual (&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &ncommon, &nparts,
			NULL, NULL,
			&edgecut, &epart[0], &npart[0]);

    /*
    METIS_PartMeshNodal (&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &nparts,
			 NULL, NULL,
			 &edgecut, &epart[0], &npart[0]);
    */
    PrintMessage (3, "metis complete");
    // cout << "done" << endl;

    for (int i = 0; i < GetNE(); i++)
      VolumeElement(i+1).SetPartition(epart[i] + 1);
    for (int i = 0; i < GetNSE(); i++)
      SurfaceElement(i+1).SetPartition(epart[i+GetNE()] + 1);
    for (int i = 0; i < GetNSeg(); i++)
      LineSegment(i+1).SetPartition(epart[i+GetNE()+GetNSE()] + 1);



    if (GetDimension() == 3)
      {
	
	// surface elements attached to volume elements
	Array<bool, PointIndex::BASE> boundarypoints (GetNP());
	boundarypoints = false;
	
	for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  {
	    const Element2d & el = (*this)[sei];
	    for (int j = 0; j < el.GetNP(); j++)
	      boundarypoints[el[j]] = true;
	  }
	// Build Pnt2Element table, boundary points only
	Array<int, PointIndex::BASE> cnt(GetNP());
	cnt = 0;
	for (ElementIndex ei = 0; ei < GetNE(); ei++)
	  {
	    const Element & el = (*this)[ei];
	    for (int j = 0; j < el.GetNP(); j++)
	      if (boundarypoints[el[j]])
		cnt[el[j]]++;
	  }
	TABLE<ElementIndex, PointIndex::BASE> pnt2el(cnt);
	cnt = 0;
	for (ElementIndex ei = 0; ei < GetNE(); ei++)
	  {
	    const Element & el = (*this)[ei];
	    for (int j = 0; j < el.GetNP(); j++)
	      if (boundarypoints[el[j]])
		pnt2el.Add (el[j], ei);
	  }
	
	for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
	  {
	    Element2d & sel = (*this)[sei];
	    PointIndex pi1 = sel[0];
	    FlatArray<ElementIndex> els = pnt2el[pi1];
	    
	    sel.SetPartition (-1);
	    
	    for (int j = 0; j < els.Size(); j++)
	      {
		const Element & el = (*this)[els[j]];
		
		bool hasall = true;
		
		for (int k = 0; k < sel.GetNP(); k++)
		  {
		    bool haspi = false;
		    for (int l = 0; l < el.GetNP(); l++)
		      if (sel[k] == el[l])
			haspi = true;

		    if (!haspi) hasall = false;
		  }
		
		if (hasall)
		  {
		    sel.SetPartition (el.GetPartition());
		    break;
		  }
	      }
	    if (sel.GetPartition() == -1)
	      cerr << "no volume element found" << endl;
	  }


	for (SegmentIndex si = 0; si < GetNSeg(); si++)
	  {
	    Segment & sel = (*this)[si];
	    PointIndex pi1 = sel[0];
	    FlatArray<ElementIndex> els = pnt2el[pi1];
	    
	    sel.SetPartition (-1);
	    
	    for (int j = 0; j < els.Size(); j++)
	      {
		const Element & el = (*this)[els[j]];
		
		bool haspi[9] = { false };  // max surfnp
		
		for (int k = 0; k < 2; k++)
		  for (int l = 0; l < el.GetNP(); l++)
		    if (sel[k] == el[l])
		      haspi[k] = true;
		
		bool hasall = true;
		for (int k = 0; k < sel.GetNP(); k++)
		  if (!haspi[k]) hasall = false;
		
		if (hasall)
		  {
		    sel.SetPartition (el.GetPartition());
		    break;
		  }
	      }
	    if (sel.GetPartition() == -1)
	      cerr << "no volume element found" << endl;
	  }
	
	
      }
  }

#endif





//========================== weights =================================================================



  // distribute the mesh to the slave processors
  // call it only for the master !
  void Mesh :: Distribute (Array<int> & volume_weights , Array<int>  & surface_weights, Array<int>  & segment_weights)
  {
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if (id != 0 || ntasks == 1 ) return;

#ifdef METIS
    ParallelMetis (volume_weights, surface_weights, segment_weights);
#else
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      (*this)[ei].SetPartition(ntasks * ei/GetNE() + 1);
#endif

    /*
    for (ElementIndex ei = 0; ei < GetNE(); ei++)
      *testout << "el(" << ei << ") is in part " << (*this)[ei].GetPartition() << endl;
    for (SurfaceElementIndex ei = 0; ei < GetNSE(); ei++)
      *testout << "sel(" << int(ei) << ") is in part " << (*this)[ei].GetPartition() << endl;
      */
    
    // MyMPI_SendCmd ("mesh");
    SendRecvMesh (); 
  }
  

#ifdef METIS5
  void Mesh :: ParallelMetis (Array<int> & volume_weights , Array<int> & surface_weights, Array<int> & segment_weights)  
  {
    PrintMessage (3, "call metis 5 with weights ...");
    
    // cout << "segment_weights " << segment_weights << endl;
    // cout << "surface_weights " << surface_weights << endl;
    // cout << "volume_weights " << volume_weights << endl;

    int timer = NgProfiler::CreateTimer ("Mesh::Partition");
    NgProfiler::RegionTimer reg(timer);

    idx_t ne = GetNE() + GetNSE() + GetNSeg();
    idx_t nn = GetNP();
    
    Array<idx_t> eptr, eind , nwgt;
    for (int i = 0; i < GetNE(); i++)
      {
	eptr.Append (eind.Size());
	
	const Element & el = VolumeElement(i+1);
	
	int ind = el.GetIndex();	
	if (volume_weights.Size()<ind)
	    nwgt.Append(0);
	else
	    nwgt.Append (volume_weights[ind -1]);
	
	for (int j = 0; j < el.GetNP(); j++)
	  eind.Append (el[j]-1);
      }
    for (int i = 0; i < GetNSE(); i++)
      {
	eptr.Append (eind.Size());
	const Element2d & el = SurfaceElement(i+1);
	
	
	int ind = el.GetIndex(); 
	ind = GetFaceDescriptor(ind).BCProperty();
	if (surface_weights.Size()<ind)
	    nwgt.Append(0);
	else
	    nwgt.Append (surface_weights[ind -1]);

	
	for (int j = 0; j < el.GetNP(); j++)
	  eind.Append (el[j]-1);
      }
    for (int i = 0; i < GetNSeg(); i++)
      {
	eptr.Append (eind.Size());
	
	const Segment & el = LineSegment(i+1);	
	
	int ind = el.si;
	if (segment_weights.Size()<ind)
	    nwgt.Append(0);
	else
	    nwgt.Append (segment_weights[ind -1]);
	
	eind.Append (el[0]);
	eind.Append (el[1]);
      }
      
    eptr.Append (eind.Size());
    Array<idx_t> epart(ne), npart(nn);

    idxtype nparts = ntasks-1;
    idxtype edgecut;


    idxtype ncommon = 3;
    METIS_PartMeshDual (&ne, &nn, &eptr[0], &eind[0], &nwgt[0], NULL, &ncommon, &nparts,
			NULL, NULL,
			&edgecut, &epart[0], &npart[0]);
    /*
    METIS_PartMeshNodal (&ne, &nn, &eptr[0], &eind[0], NULL, NULL, &nparts,
			 NULL, NULL,
			 &edgecut, &epart[0], &npart[0]);
    */
    PrintMessage (3, "metis complete");
    // cout << "done" << endl;

    for (int i = 0; i < GetNE(); i++)
      VolumeElement(i+1).SetPartition(epart[i] + 1);
    for (int i = 0; i < GetNSE(); i++)
      SurfaceElement(i+1).SetPartition(epart[i+GetNE()] + 1);
    for (int i = 0; i < GetNSeg(); i++)
      LineSegment(i+1).SetPartition(epart[i+GetNE()+GetNSE()] + 1);
  }
#endif 



//===========================================================================================









#ifdef METIS4
  void Mesh :: ParallelMetis ( )  
  {
    int timer = NgProfiler::CreateTimer ("Mesh::Partition");
    NgProfiler::RegionTimer reg(timer);

    PrintMessage (3, "Metis called");
      
    if (GetDimension() == 2) 
      {
	PartDualHybridMesh2D ( ); // neloc );
	return;
      }


    idx_t ne = GetNE();
    idx_t nn = GetNP();

    if (ntasks <= 2 || ne <= 1)
      {
        if (ntasks == 1) return;
        
        for (int i=1; i<=ne; i++)
          VolumeElement(i).SetPartition(1);

        for (int i=1; i<=GetNSE(); i++)
          SurfaceElement(i).SetPartition(1);

        return;
      }


    bool uniform_els = true;

    ELEMENT_TYPE elementtype = TET; 
    for (int el = 1; el <= GetNE(); el++)
      if (VolumeElement(el).GetType() != elementtype)
	{
	  uniform_els = false;
	  break;
	}


    if (!uniform_els)
      {
	PartHybridMesh ();  
      }
    else
      {
	
	// uniform (TET) mesh,  JS
	int npe = VolumeElement(1).GetNP();
	Array<idxtype> elmnts(ne*npe);
	
	int etype;
	if (elementtype == TET)
	  etype = 2;
	else if (elementtype == HEX)
	  etype = 3;
	
    
	for (int i=1; i<=ne; i++)
	  for (int j=1; j<=npe; j++)
	    elmnts[(i-1)*npe+(j-1)] = VolumeElement(i).PNum(j)-1;
	
	int numflag = 0;
	int nparts = ntasks-1;
	int ncommon = 3;
	int edgecut;
	Array<idxtype> epart(ne), npart(nn);
	
	//     if ( ntasks == 1 ) 
	//       {
	// 	(*this) = *mastermesh;
	// 	nparts = 4;	   
	// 	metis :: METIS_PartMeshDual (&ne, &nn, elmnts, &etype, &numflag, &nparts,
	// 				     &edgecut, epart, npart);
	// 	cout << "done" << endl;
	
	// 	cout << "edge-cut: " << edgecut << ", balance: " << metis :: ComputeElementBalance(ne, nparts, epart) << endl;
	
	// 	for (int i=1; i<=ne; i++)
	// 	  {
	// 	    mastermesh->VolumeElement(i).SetPartition(epart[i-1]);
	// 	  }
	
	// 	return;
	//       }
	
	
	int timermetis = NgProfiler::CreateTimer ("Metis itself");
	NgProfiler::StartTimer (timermetis);
	
#ifdef METIS4
	cout << "call metis(4)_PartMeshDual ... " << flush;
	METIS_PartMeshDual (&ne, &nn, &elmnts[0], &etype, &numflag, &nparts,
			    &edgecut, &epart[0], &npart[0]);
#else
	cout << "call metis(5)_PartMeshDual ... " << endl;
	// idx_t options[METIS_NOPTIONS];
	
	Array<idx_t> eptr(ne+1);
	for (int j = 0; j < ne+1; j++)
	  eptr[j] = 4*j;
	
	METIS_PartMeshDual (&ne, &nn, &eptr[0], &elmnts[0], NULL, NULL, &ncommon, &nparts,
			    NULL, NULL,
			    &edgecut, &epart[0], &npart[0]);
#endif
	
	NgProfiler::StopTimer (timermetis);
	
	cout << "complete" << endl;
#ifdef METIS4
	cout << "edge-cut: " << edgecut << ", balance: " 
	     << ComputeElementBalance(ne, nparts, &epart[0]) << endl;
#endif
	
	// partition numbering by metis : 0 ...  ntasks - 1
	// we want:                       1 ...  ntasks
	for (int i=1; i<=ne; i++)
	  VolumeElement(i).SetPartition(epart[i-1] + 1);
      }
    

    for (int sei = 1; sei <= GetNSE(); sei++ )
      {
	int ei1, ei2;
	GetTopology().GetSurface2VolumeElement (sei, ei1, ei2);
	Element2d & sel = SurfaceElement (sei);

        for (int j = 0; j < 2; j++)
          {
            int ei = (j == 0) ? ei1 : ei2;
            if ( ei > 0 && ei <= GetNE() )
              {
		sel.SetPartition (VolumeElement(ei).GetPartition());
		break;
	      }
	  }	
      }
    
  }
#endif


  void Mesh :: PartHybridMesh () 
  {
#ifdef METIS
    int ne = GetNE();
    
    int nn = GetNP();
    int nedges = topology->GetNEdges();

    idxtype  *xadj, * adjacency, *v_weights = NULL, *e_weights = NULL;

    int weightflag = 0;
    int numflag = 0;
    int nparts = ntasks - 1;

    int options[5];
    options[0] = 0;
    int edgecut;
    idxtype * part;

    xadj = new idxtype[nn+1];
    part = new idxtype[nn];

    Array<int> cnt(nn+1);
    cnt = 0;

    for ( int edge = 1; edge <= nedges; edge++ )
      {
	int v1, v2;
	topology->GetEdgeVertices ( edge, v1, v2);
	cnt[v1-1] ++;
	cnt[v2-1] ++;
      }

    xadj[0] = 0;
    for ( int n = 1; n <= nn; n++ )
      {
	xadj[n] = idxtype(xadj[n-1] + cnt[n-1]); 
      }

    adjacency = new idxtype[xadj[nn]];
    cnt = 0;

    for ( int edge = 1; edge <= nedges; edge++ )
      {
	int v1, v2;
	topology->GetEdgeVertices ( edge, v1, v2);
	adjacency[ xadj[v1-1] + cnt[v1-1] ] = v2-1;
	adjacency[ xadj[v2-1] + cnt[v2-1] ] = v1-1;
	cnt[v1-1]++;
	cnt[v2-1]++;
      }

    for ( int vert = 0; vert < nn; vert++ )
      {
	FlatArray<idxtype> array ( cnt[vert], &adjacency[ xadj[vert] ] );
	BubbleSort(array);
      }

#ifdef METIS4
    METIS_PartGraphKway ( &nn, xadj, adjacency, v_weights, e_weights, &weightflag, 
			  &numflag, &nparts, options, &edgecut, part );
#else
    cout << "currently not supported (metis5), A" << endl;
#endif

    Array<int> nodesinpart(ntasks);

    for ( int el = 1; el <= ne; el++ )
      {
	Element & volel = VolumeElement(el);
	nodesinpart = 0;

	
	int el_np = volel.GetNP();
	int partition = 0; 
	for ( int i = 0; i < el_np; i++ )
	  nodesinpart[ part[volel[i]-1]+1 ] ++;

	for ( int i = 1; i < ntasks; i++ )
	  if ( nodesinpart[i] > nodesinpart[partition] ) 
	    partition = i;

	volel.SetPartition(partition);
      }

    delete [] xadj;
    delete [] part;
    delete [] adjacency;
#else
    cout << "parthybridmesh not available" << endl;
#endif
  }


  void Mesh :: PartDualHybridMesh ( ) // Array<int> & neloc ) 
  {
#ifdef METIS
    int ne = GetNE();
    
    // int nn = GetNP();
    // int nedges = topology->GetNEdges();
    int nfaces = topology->GetNFaces();

    idxtype  *xadj, * adjacency, *v_weights = NULL, *e_weights = NULL;

    int weightflag = 0;
    // int numflag = 0;
    int nparts = ntasks - 1;

    int options[5];
    options[0] = 0;
    int edgecut;
    idxtype * part;

    Array<int, 0> facevolels1(nfaces), facevolels2(nfaces);
    facevolels1 = -1;
    facevolels2 = -1;

    Array<int, 0> elfaces;
    xadj = new idxtype[ne+1];
    part = new idxtype[ne];

    Array<int, 0> cnt(ne+1);
    cnt = 0;

    for ( int el=1; el <= ne; el++ )
      {
	Element volel = VolumeElement(el);
	topology->GetElementFaces(el, elfaces);
	for ( int i = 0; i < elfaces.Size(); i++ )
	  {
	    if ( facevolels1[elfaces[i]-1] == -1 )
	      facevolels1[elfaces[i]-1] = el;
	    else
	      {
		facevolels2[elfaces[i]-1] = el;
		cnt[facevolels1[elfaces[i]-1]-1]++;
		cnt[facevolels2[elfaces[i]-1]-1]++;
	      }
	  }
      }

    xadj[0] = 0;
    for ( int n = 1; n <= ne; n++ )
      {
	xadj[n] = idxtype(xadj[n-1] + cnt[n-1]); 
      }

    adjacency = new idxtype[xadj[ne]];
    cnt = 0;

    for ( int face = 1; face <= nfaces; face++ )
      {
	int e1, e2;
	e1 = facevolels1[face-1];
	e2 = facevolels2[face-1];
	if ( e2 == -1 ) continue;
	adjacency[ xadj[e1-1] + cnt[e1-1] ] = e2-1;
	adjacency[ xadj[e2-1] + cnt[e2-1] ] = e1-1;
	cnt[e1-1]++;
	cnt[e2-1]++;
      }

    for ( int el = 0; el < ne; el++ )
      {
	FlatArray<idxtype> array ( cnt[el], &adjacency[ xadj[el] ] );
	BubbleSort(array);
      }

    int timermetis = NgProfiler::CreateTimer ("Metis itself");
    NgProfiler::StartTimer (timermetis);

#ifdef METIS4
    METIS_PartGraphKway ( &ne, xadj, adjacency, v_weights, e_weights, &weightflag, 
			  &numflag, &nparts, options, &edgecut, part );
#else
    cout << "currently not supported (metis5), B" << endl;
#endif


    NgProfiler::StopTimer (timermetis);

    Array<int> nodesinpart(ntasks);

    for ( int el = 1; el <= ne; el++ )
      {
	// Element & volel = VolumeElement(el);
	nodesinpart = 0;

	VolumeElement(el).SetPartition(part[el-1 ] + 1);
	
      }

    /*    
    for ( int i=1; i<=ne; i++)
      {
	neloc[ VolumeElement(i).GetPartition() ] ++;
      }
    */

    delete [] xadj;
    delete [] part;
    delete [] adjacency;
#else
    cout << "partdualmesh not available" << endl;
#endif

  }





  void Mesh :: PartDualHybridMesh2D ( ) 
  {
#ifdef METIS
    idxtype ne = GetNSE();
    int nv = GetNV();

    Array<idxtype> xadj(ne+1);
    Array<idxtype> adjacency(ne*4);

    // first, build the vertex 2 element table:
    Array<int, PointIndex::BASE> cnt(nv);
    cnt = 0;
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      for (int j = 0; j < (*this)[sei].GetNP(); j++)
	cnt[ (*this)[sei][j] ] ++;
    
    TABLE<SurfaceElementIndex, PointIndex::BASE> vert2els(cnt);
    for (SurfaceElementIndex sei = 0; sei < GetNSE(); sei++)
      for (int j = 0; j < (*this)[sei].GetNP(); j++)
	vert2els.Add ((*this)[sei][j], sei);
    

    // find all neighbour elements
    int cntnb = 0;
    Array<int> marks(ne);   // to visit each neighbour just once
    marks = -1;
    for (SurfaceElementIndex sei = 0; sei < ne; sei++)
      {
	xadj[sei] = cntnb;
	for (int j = 0; j < (*this)[sei].GetNP(); j++)
	  {
	    PointIndex vnr = (*this)[sei][j];

	    // all elements with at least one common vertex
	    for (int k = 0; k < vert2els[vnr].Size(); k++)   
	      {
		SurfaceElementIndex sei2 = vert2els[vnr][k];
		if (sei == sei2) continue;
		if (marks[sei2] == sei) continue;
		
		// neighbour, if two common vertices
		int common = 0;
		for (int m1 = 0; m1 < (*this)[sei].GetNP(); m1++)
		  for (int m2 = 0; m2 < (*this)[sei2].GetNP(); m2++)
		    if ( (*this)[sei][m1] == (*this)[sei2][m2])
		      common++;
		
		if (common >= 2)
		  {
		    marks[sei2] = sei;     // mark as visited
		    adjacency[cntnb++] = sei2;
		  }
	      }
	  }
      }
    xadj[ne] = cntnb;

    idxtype *v_weights = NULL, *e_weights = NULL;

    idxtype weightflag = 0;
    // int numflag = 0;
    idxtype nparts = ntasks - 1;

    idxtype edgecut;
    Array<idxtype> part(ne);

    for ( int el = 0; el < ne; el++ )
      BubbleSort (adjacency.Range (xadj[el], xadj[el+1]));

#ifdef METIS4	
    int options[5];
    options[0] = 0;
    METIS_PartGraphKway ( &ne, &xadj[0], &adjacency[0], v_weights, e_weights, &weightflag, 
			  &numflag, &nparts, options, &edgecut, &part[0] );
#else
    idx_t ncon = 1;
    METIS_PartGraphKway ( &ne, &ncon, &xadj[0], &adjacency[0], 
			  v_weights, NULL, e_weights, 
			  &nparts, 
			  NULL, NULL, NULL,
			  &edgecut, &part[0] );
#endif


    for (SurfaceElementIndex sei = 0; sei < ne; sei++)
      (*this) [sei].SetPartition (part[sei]+1);
#else
    cout << "partdualmesh not available" << endl;
#endif

  }



}



#endif
