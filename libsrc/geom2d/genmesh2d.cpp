#include <meshing.hpp>
#include <geometry2d.hpp>

namespace netgen
{
  extern DLL_HEADER MeshingParameters mparam;

  extern void Optimize2d (Mesh & mesh, MeshingParameters & mp);



  void CalcPartition (double l, double h, double h1, double h2,
		      double hcurve, double elto0, Array<double> & points);

  // partitionizes spline curve
  void Partition (const SplineSegExt & spline,
		  double h, double elto0,
		  Mesh & mesh, Point3dTree & searchtree, int segnr) 
  {
    enum { D = 2 };
    int i, j;
    double l; // , r1, r2, ra;
    double lold, dt, frac;
    int n = 100;
    Point<D> p, pold, mark, oldmark;
    Array<double> curvepoints;
    double edgelength, edgelengthold;
    l = spline.Length();

    double h1 = min (spline.StartPI().hmax, h/spline.StartPI().refatpoint);
    double h2 = min (spline.EndPI().hmax, h/spline.EndPI().refatpoint);
    double hcurve = min (spline.hmax, h/spline.reffak);

    CalcPartition (l, h, h1, h2, hcurve, elto0, curvepoints);
    //  cout << "curvepoints = " << curvepoints << endl;

    dt = 1.0 / n;

    l = 0;
    j = 1;

    pold = spline.GetPoint (0);
    lold = 0;
    oldmark = pold;
    edgelengthold = 0;
    Array<int> locsearch;
    
    for (i = 1; i <= n; i++)
      {
	p = spline.GetPoint (i*dt);
	l = lold + Dist (p, pold);
	while (j < curvepoints.Size() && (l >= curvepoints[j] || i == n))
	  {
	    frac = (curvepoints[j]-lold) / (l-lold);
	    edgelength = i*dt + (frac-1)*dt;
	    // mark = pold + frac * (p-pold);
	    mark = spline.GetPoint (edgelength);
	  
	    // cout << "mark = " << mark << " =?= " << GetPoint (edgelength) << endl;

	    {
	      PointIndex pi1 = -1, pi2 = -1;
	  
	      Point3d mark3(mark(0), mark(1), 0);
	      Point3d oldmark3(oldmark(0), oldmark(1), 0);

	      Vec<3> v (1e-4*h, 1e-4*h, 1e-4*h);
	      searchtree.GetIntersecting (oldmark3 - v, oldmark3 + v, locsearch);

	      for (int k = 0; k < locsearch.Size(); k++)
		if ( mesh[PointIndex(locsearch[k])].GetLayer() == spline.layer)
		  pi1 = locsearch[k];
	      // if (locsearch.Size()) pi1 = locsearch[0];
	      
	      searchtree.GetIntersecting (mark3 - v, mark3 + v, locsearch);
	      for (int k = 0; k < locsearch.Size(); k++)
		if ( mesh[PointIndex(locsearch[k])].GetLayer() == spline.layer)
		  pi2 = locsearch[k];
	      // if (locsearch.Size()) pi2 = locsearch[0];

	      /*	    
		for (PointIndex pk = PointIndex::BASE; 
		pk < mesh.GetNP()+PointIndex::BASE; pk++)
		{
		if (Dist (mesh[pk], oldmark3) < 1e-4 * h) pi1 = pk;
		if (Dist (mesh[pk], mark3) < 1e-4 * h) pi2 = pk;
		}
	      */
	    
	      
	      //	    cout << "pi1 = " << pi1 << endl;
	      //	    cout << "pi2 = " << pi2 << endl;
	    
	      if (pi1 == -1)
		{
		  pi1 = mesh.AddPoint(oldmark3, spline.layer);
		  searchtree.Insert (oldmark3, pi1);
		}
	      if (pi2 == -1)
		{
		  pi2 = mesh.AddPoint(mark3, spline.layer);
		  searchtree.Insert (mark3, pi2);
		}

	      Segment seg;
	      seg.edgenr = segnr;
	      seg.si = spline.bc; // segnr;
	      seg[0] = pi1;
	      seg[1] = pi2;
	      seg.domin = spline.leftdom;
	      seg.domout = spline.rightdom;
	      seg.epgeominfo[0].edgenr = segnr;
	      seg.epgeominfo[0].dist = edgelengthold;
	      seg.epgeominfo[1].edgenr = segnr;
	      seg.epgeominfo[1].dist = edgelength;
	      seg.singedge_left = spline.hpref_left;
	      seg.singedge_right = spline.hpref_right;
	      mesh.AddSegment (seg);
	    }
	
	    oldmark = mark;
	    edgelengthold = edgelength;
	    j++;
	  }
    
	pold = p;
	lold = l;
      }
  }



  void SplineGeometry2d :: PartitionBoundary (double h, Mesh & mesh2d)
  {
    enum { D = 2 };
    Box<D> bbox;
    GetBoundingBox (bbox);
    double dist = Dist (bbox.PMin(), bbox.PMax());
    Point<3> pmin;
    Point<3> pmax;
  
    pmin(2) = -dist; pmax(2) = dist;
    for(int j=0;j<D;j++)
      {
	pmin(j) = bbox.PMin()(j);
	pmax(j) = bbox.PMax()(j);
      }

    Point3dTree searchtree (pmin, pmax);

    for (int i = 0; i < splines.Size(); i++)
      for (int side = 0; side <= 1; side++)
	{
	  int dom = (side == 0) ? GetSpline(i).leftdom : GetSpline(i).rightdom;
	  if (dom != 0) GetSpline(i).layer = GetDomainLayer (dom);
	}

    for (int i = 0; i < splines.Size(); i++)
      if (GetSpline(i).copyfrom == -1)
	{
	  // astrid - set boundary meshsize to  domain meshsize h
	  // if no domain mesh size is given, the max h value from the bounding box is used
	  double minimum = min2 ( GetDomainMaxh ( GetSpline(i).leftdom ), GetDomainMaxh ( GetSpline(i).rightdom ) );
	  double maximum = max2 ( GetDomainMaxh ( GetSpline(i).leftdom ), GetDomainMaxh ( GetSpline(i).rightdom ) );
	  minimum = min2 ( minimum, h );
	  maximum = min2 ( maximum, h);
	  if ( minimum > 0 )
	    // GetSpline(i).Partition(minimum, elto0, mesh2d, searchtree, i+1);
	    Partition(GetSpline(i), minimum, elto0, mesh2d, searchtree, i+1);	    
	  else if ( maximum > 0 )
	    // GetSpline(i).Partition(maximum, elto0, mesh2d, searchtree, i+1);
	    Partition(GetSpline(i), maximum, elto0, mesh2d, searchtree, i+1);
	  else
	    // GetSpline(i).Partition(h, elto0, mesh2d, searchtree, i+1);
	    Partition(GetSpline(i), h, elto0, mesh2d, searchtree, i+1);
	}
      else
	{
	  CopyEdgeMesh (GetSpline(i).copyfrom, i+1, mesh2d, searchtree);
	}
  }






  void SplineGeometry2d :: CopyEdgeMesh (int from, int to, Mesh & mesh, Point3dTree & searchtree)
  {
    const int D = 2;
    int i;

    Array<int, PointIndex::BASE> mappoints (mesh.GetNP());
    Array<double, PointIndex::BASE> param (mesh.GetNP());
    mappoints = -1;
    param = 0;

    Point3d pmin, pmax;
    mesh.GetBox (pmin, pmax);
    double diam2 = Dist2(pmin, pmax);

    if (printmessage_importance>0)
      cout << "copy edge, from = " << from << " to " << to << endl;
  
    for (i = 1; i <= mesh.GetNSeg(); i++)
      {
	const Segment & seg = mesh.LineSegment(i);
	if (seg.edgenr == from)
	  {
	    mappoints.Elem(seg[0]) = 1;
	    param.Elem(seg[0]) = seg.epgeominfo[0].dist;

	    mappoints.Elem(seg[1]) = 1;
	    param.Elem(seg[1]) = seg.epgeominfo[1].dist;
	  }
      }

    bool mapped = false;
    for (i = 1; i <= mappoints.Size(); i++)
      {
	if (mappoints.Get(i) != -1)
	  {
	    Point<D> newp = splines.Get(to)->GetPoint (param.Get(i));
	    Point<3> newp3;
	    for(int j=0; j<min2(D,3); j++)
	      newp3(j) = newp(j);
	    for(int j=min2(D,3); j<3; j++)
	      newp3(j) = 0;
	  
	    int npi = -1;
	  
	    for (PointIndex pi = PointIndex::BASE; 
		 pi < mesh.GetNP()+PointIndex::BASE; pi++)
	      if (Dist2 (mesh.Point(pi), newp3) < 1e-12 * diam2)
		npi = pi;
	  
	    if (npi == -1)
	      {
		npi = mesh.AddPoint (newp3);
		searchtree.Insert (newp3, npi);
	      }

	    mappoints.Elem(i) = npi;

	    mesh.GetIdentifications().Add (i, npi, to);
	    mapped = true;
	  }
      }
    if(mapped)
      mesh.GetIdentifications().SetType(to,Identifications::PERIODIC);

    // copy segments
    int oldnseg = mesh.GetNSeg();
    for (i = 1; i <= oldnseg; i++)
      {
	const Segment & seg = mesh.LineSegment(i);
	if (seg.edgenr == from)
	  {
	    Segment nseg;
	    nseg.edgenr = to;
	    nseg.si = GetSpline(to-1).bc;      // splines.Get(to)->bc;
	    nseg[0] = mappoints.Get(seg[0]);
	    nseg[1] = mappoints.Get(seg[1]);
	    nseg.domin = GetSpline(to-1).leftdom;
	    nseg.domout = GetSpline(to-1).rightdom;
	  
	    nseg.epgeominfo[0].edgenr = to;
	    nseg.epgeominfo[0].dist = param.Get(seg[0]);
	    nseg.epgeominfo[1].edgenr = to;
	    nseg.epgeominfo[1].dist = param.Get(seg[1]);
	    mesh.AddSegment (nseg);
	  }
      }
  }





  void MeshFromSpline2D (SplineGeometry2d & geometry,
			 Mesh *& mesh, 
			 MeshingParameters & mp)
  {
    PrintMessage (1, "Generate Mesh from spline geometry");

    double h = mp.maxh;

    Box<2> bbox = geometry.GetBoundingBox ();

    if (bbox.Diam() < h) 
      {
	h = bbox.Diam();
	mp.maxh = h;
      }

    mesh = new Mesh;
    mesh->SetDimension (2);

    geometry.PartitionBoundary (h, *mesh);


    // marks mesh points for hp-refinement
    for (int i = 0; i < geometry.GetNP(); i++)
      if (geometry.GetPoint(i).hpref)
	{
	  double mindist = 1e99;
	  PointIndex mpi(0);
	  Point<2> gp = geometry.GetPoint(i);
	  Point<3> gp3(gp(0), gp(1), 0);
	  for (PointIndex pi = PointIndex::BASE; 
	       pi < mesh->GetNP()+PointIndex::BASE; pi++)
	    if (Dist2(gp3, (*mesh)[pi]) < mindist)
	      {
		mpi = pi;
		mindist = Dist2(gp3, (*mesh)[pi]);
	      }
	  (*mesh)[mpi].Singularity(1.);
	}


    int maxdomnr = 0;
    for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
      {
	if ( (*mesh)[si].domin > maxdomnr) maxdomnr = (*mesh)[si].domin;
	if ( (*mesh)[si].domout > maxdomnr) maxdomnr = (*mesh)[si].domout;
      }

    mesh->ClearFaceDescriptors();
    for (int i = 1; i <= maxdomnr; i++)
      mesh->AddFaceDescriptor (FaceDescriptor (i, 0, 0, i));

    // set Array<string*> bcnames... 
    // number of bcnames
    int maxsegmentindex = 0;
    for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
      {
	if ( (*mesh)[si].si > maxsegmentindex) maxsegmentindex = (*mesh)[si].si;
      }

    mesh->SetNBCNames(maxsegmentindex);

    for ( int sindex = 0; sindex < maxsegmentindex; sindex++ )
      mesh->SetBCName ( sindex, geometry.GetBCName( sindex+1 ) );

    for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
      (*mesh)[si].SetBCName ( (*mesh).GetBCNamePtr( (*mesh)[si].si-1 ) );

    Point3d pmin(bbox.PMin()(0), bbox.PMin()(1), -bbox.Diam());
    Point3d pmax(bbox.PMax()(0), bbox.PMax()(1), bbox.Diam());

    mesh->SetLocalH (pmin, pmax, mparam.grading);
    mesh->SetGlobalH (h);
  
    mesh->CalcLocalH(mparam.grading);

    int bnp = mesh->GetNP(); // boundary points

    int hquad = mparam.quad;


    for (int domnr = 1; domnr <= maxdomnr; domnr++)
      if (geometry.GetDomainTensorMeshing (domnr))
        { // tensor product mesh
          
          Array<PointIndex, PointIndex::BASE> nextpi(bnp);
          Array<int, PointIndex::BASE> si1(bnp), si2(bnp);
          PointIndex firstpi;
          
          nextpi = -1;
          si1 = -1;
          si2 = -1;
          for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
            {
              int p1 = -1, p2 = -2;

              if ( (*mesh)[si].domin == domnr)
                { p1 = (*mesh)[si][0]; p2 = (*mesh)[si][1]; }
              if ( (*mesh)[si].domout == domnr)
                { p1 = (*mesh)[si][1]; p2 = (*mesh)[si][0]; }
              
              if (p1 == -1) continue;

              nextpi[p1] = p2;       // counter-clockwise
              
              int index = (*mesh)[si].si;
              if (si1[p1] != index && si2[p1] != index)
                { si2[p1] = si1[p1]; si1[p1] = index; }
              if (si1[p2] != index && si2[p2] != index)
                { si2[p2] = si1[p2]; si1[p2] = index; }
            }

          PointIndex c1(0), c2, c3, c4;  // 4 corner points
          int nex = 1, ney = 1;

          for (PointIndex pi = 1; pi <= si2.Size(); pi++)
            if (si2[pi] != -1)
              { c1 = pi; break; }      

          for (c2 = nextpi[c1]; si2[c2] == -1; c2 = nextpi[c2], nex++);
          for (c3 = nextpi[c2]; si2[c3] == -1; c3 = nextpi[c3], ney++);
          for (c4 = nextpi[c3]; si2[c4] == -1; c4 = nextpi[c4]);



          Array<PointIndex> pts ( (nex+1) * (ney+1) );   // x ... inner loop
          pts = -1;

          for (PointIndex pi = c1, i = 0; pi != c2; pi = nextpi[pi], i++)
            pts[i] = pi;
          for (PointIndex pi = c2, i = 0; pi != c3; pi = nextpi[pi], i++)
            pts[(nex+1)*i+nex] = pi;
          for (PointIndex pi = c3, i = 0; pi != c4; pi = nextpi[pi], i++)
            pts[(nex+1)*(ney+1)-i-1] = pi;
          for (PointIndex pi = c4, i = 0; pi != c1; pi = nextpi[pi], i++)
            pts[(nex+1)*(ney-i)] = pi;


          for (PointIndex pix = nextpi[c1], ix = 0; pix != c2; pix = nextpi[pix], ix++)
            for (PointIndex piy = nextpi[c2], iy = 0; piy != c3; piy = nextpi[piy], iy++)
              {
                Point<3> p = (*mesh)[pix] + ( (*mesh)[piy] - (*mesh)[c2] );
                pts[(nex+1)*(iy+1) + ix+1] = mesh -> AddPoint (p , 1, FIXEDPOINT);
              }

          for (int i = 0; i < ney; i++)
            for (int j = 0; j < nex; j++)
              {
                Element2d el(QUAD);
                el[0] = pts[i*(nex+1)+j];
                el[1] = pts[i*(nex+1)+j+1];
                el[2] = pts[(i+1)*(nex+1)+j+1];
                el[3] = pts[(i+1)*(nex+1)+j];
                el.SetIndex (domnr);

                mesh -> AddSurfaceElement (el);
              }
        }




    for (int domnr = 1; domnr <= maxdomnr; domnr++)
      {
        if (geometry.GetDomainTensorMeshing (domnr)) continue;

	if ( geometry.GetDomainMaxh ( domnr ) > 0 )
	  h = geometry.GetDomainMaxh(domnr);


	PrintMessage (3, "Meshing domain ", domnr, " / ", maxdomnr);

	int oldnf = mesh->GetNSE();

        mparam.quad = hquad || geometry.GetDomainQuadMeshing (domnr);

	Meshing2 meshing (mparam, Box<3> (pmin, pmax));

	Array<int, PointIndex::BASE> compress(bnp);
	compress = -1;
	int cnt = 0;
	for (PointIndex pi = PointIndex::BASE; pi < bnp+PointIndex::BASE; pi++)
	  if ( (*mesh)[pi].GetLayer() == geometry.GetDomainLayer(domnr))
	    {
	      meshing.AddPoint ( (*mesh)[pi], pi);
	      cnt++;
	      compress[pi] = cnt;
	    }

	PointGeomInfo gi;
	gi.trignum = 1;
	for (SegmentIndex si = 0; si < mesh->GetNSeg(); si++)
	  {
	    if ( (*mesh)[si].domin == domnr)
	      {
		meshing.AddBoundaryElement ( compress[(*mesh)[si][0]], 
					     compress[(*mesh)[si][1]], gi, gi);
	      }
	    if ( (*mesh)[si].domout == domnr)
	      {
		meshing.AddBoundaryElement ( compress[(*mesh)[si][1]],
					     compress[(*mesh)[si][0]], gi, gi);
		
	      }

	  }


	mparam.checkoverlap = 0;
	
	/*
	if (!mparam.quad)
	  meshing.Delaunay (*mesh, domnr, mparam);
	else
	*/
	meshing.GenerateMesh (*mesh, mparam, h, domnr);

	for (SurfaceElementIndex sei = oldnf; sei < mesh->GetNSE(); sei++)
	  (*mesh)[sei].SetIndex (domnr);


	// astrid
	char * material;
	geometry.GetMaterial( domnr, material );
	if ( material )
	  {
	    (*mesh).SetMaterial ( domnr,  material );
	  }

      }

    mparam.quad = hquad;



    int hsteps = mp.optsteps2d;

    mp.optimize2d = "smcm"; 
    mp.optsteps2d = hsteps/2;
    Optimize2d (*mesh, mp);

    mp.optimize2d = "Smcm"; 
    mp.optsteps2d = (hsteps+1)/2;
    Optimize2d (*mesh, mp);

    mp.optsteps2d = hsteps;

    mesh->Compress();
    mesh -> SetNextMajorTimeStamp();


    extern DLL_HEADER void Render();
    Render();
  }








}
