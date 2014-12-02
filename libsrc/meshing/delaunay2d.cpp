#include <mystdlib.h>
#include "meshing.hpp"

// not yet working ....

namespace netgen
{



  class DelaunayTrig
  {
    PointIndex pnums[3];
    Point<3> c;
    double r;
    double rad2;
  public:    
    DelaunayTrig () { ; }
    DelaunayTrig (int p1, int p2, int p3)
    { pnums[0] = p1; pnums[1] = p2; pnums[2] = p3; }

    PointIndex & operator[] (int j) { return pnums[j]; }
    const PointIndex & operator[] (int j) const { return pnums[j]; }

    void CalcCenter (Mesh & mesh)
    {
      Point<3> p1 = mesh[pnums[0]];
      Point<3> p2 = mesh[pnums[1]];
      Point<3> p3 = mesh[pnums[2]];
      Vec<3> v1 = p2-p1;
      Vec<3> v2 = p3-p1;
      Mat<2,2> mat, inv;
      mat(0,0) = v1*v1;
      mat(0,1) = v1*v2;
      mat(1,0) = v2*v1;
      mat(1,1) = v2*v2;
      Vec<2> rhs, sol;
      rhs(0) = 0.5 * v1*v1;
      rhs(1) = 0.5 * v2*v2;
      CalcInverse (mat, inv);
      sol = inv * rhs;
      
      c = p1 + sol(0) * v1 + sol(1) * v2;
      rad2 = Dist2(c, p1);
      r = sqrt(rad2); 
    }

    Point<3> Center() const { return c; }
    double Radius2() const { return rad2; }
    Box<3> BoundingBox() const { return Box<3> (c-Vec<3>(r,r,0.1), c+Vec<3>(r,r,0.1)); }
  };

  ostream & operator<< (ostream & ost, DelaunayTrig trig)
  {
    ost << trig[0] << "-" << trig[1] << "-" << trig[2] << endl;
    return ost;
  }

  
  void Meshing2 :: BlockFillLocalH (Mesh & mesh, const MeshingParameters & mp)
  {
    static int timer = NgProfiler::CreateTimer ("Meshing2::BlockFill");
    static int timer1 = NgProfiler::CreateTimer ("Meshing2::BlockFill 1");
    static int timer2 = NgProfiler::CreateTimer ("Meshing2::BlockFill 2");
    static int timer3 = NgProfiler::CreateTimer ("Meshing2::BlockFill 3");
    static int timer4 = NgProfiler::CreateTimer ("Meshing2::BlockFill 4");
    NgProfiler::RegionTimer reg (timer);

    NgProfiler::StartTimer (timer1);

    double filldist = mp.filldist;
    
    cout << "blockfill local h" << endl;
    cout << "rel filldist = " << filldist << endl;
    PrintMessage (3, "blockfill local h");

    Array<Point<3> > npoints;
    
    // adfront -> CreateTrees();

    Box<3> bbox ( Box<3>::EMPTY_BOX );
    double maxh = 0;
    
    for (int i = 0; i < adfront->GetNFL(); i++)
      {
	const FrontLine & line = adfront->GetLine (i);

	const Point<3> & p1 = adfront->GetPoint(line.L().I1());
	const Point<3> & p2 = adfront->GetPoint(line.L().I2());
	
        maxh = max (maxh, Dist (p1, p2));
	
	bbox.Add (p1);
	bbox.Add (p2);
      }

    
    cout << "bbox = " << bbox << endl;


    // Point<3> mpc = bbox.Center();
    bbox.Increase (bbox.Diam()/2);
    Box<3> meshbox = bbox;

    NgProfiler::StopTimer (timer1);    
    NgProfiler::StartTimer (timer2);    


    LocalH loch2 (bbox, 1, 2);
    
    if (mp.maxh < maxh) maxh = mp.maxh;
    
    bool changed;
    do 
      {
	mesh.LocalHFunction().ClearFlags();
	
	for (int i = 0; i < adfront->GetNFL(); i++)
	  {
	    const FrontLine & line = adfront->GetLine(i);
	    
	    Box<3> bbox (adfront->GetPoint (line.L().I1()));
	    bbox.Add (adfront->GetPoint (line.L().I2()));

	    
	    double filld = filldist * bbox.Diam();
	    bbox.Increase (filld);
	    
	    mesh.LocalHFunction().CutBoundary (bbox); 
	  }
	

	mesh.LocalHFunction().FindInnerBoxes (adfront, NULL);
	
	npoints.SetSize(0);
	mesh.LocalHFunction().GetInnerPoints (npoints);

	changed = false;
	for (int i = 0; i < npoints.Size(); i++)
	  {
	    if (mesh.LocalHFunction().GetH(npoints[i]) > 1.2 * maxh)
	      {
		mesh.LocalHFunction().SetH (npoints[i], maxh);
		changed = true;
	      }
	  }
      }
    while (changed);

    NgProfiler::StopTimer (timer2);    
    NgProfiler::StartTimer (timer3);    


    if (debugparam.slowchecks)
      (*testout) << "Blockfill with points: " << endl;
    *testout << "loch = " << mesh.LocalHFunction() << endl;
    
    *testout << "npoints = " << endl << npoints << endl;

    for (int i = 1; i <= npoints.Size(); i++)
      {
	if (meshbox.IsIn (npoints.Get(i)))
	  {
	    PointIndex gpnum = mesh.AddPoint (npoints.Get(i));
	    adfront->AddPoint (npoints.Get(i), gpnum);
	    
	    if (debugparam.slowchecks)
	      {
		(*testout) << npoints.Get(i) << endl;

		Point<2> p2d (npoints.Get(i)(0), npoints.Get(i)(1));
		if (!adfront->Inside(p2d))
		  {
		    cout << "add outside point" << endl;
		    (*testout) << "outside" << endl;
		  }
	      }
	    
	  }
      }
    
    NgProfiler::StopTimer (timer3);    
    NgProfiler::StartTimer (timer4);    
  

  // find outer points
  
    loch2.ClearFlags();

    for (int i = 0; i < adfront->GetNFL(); i++)
      {
	const FrontLine & line = adfront->GetLine(i);
	
	Box<3> bbox (adfront->GetPoint (line.L().I1()));
	bbox.Add (adfront->GetPoint (line.L().I2()));
	
	loch2.SetH (bbox.Center(), bbox.Diam());
      }


    for (int i = 0; i < adfront->GetNFL(); i++)
      {
	const FrontLine & line = adfront->GetLine(i);
	
	Box<3> bbox (adfront->GetPoint (line.L().I1()));
	bbox.Add (adfront->GetPoint (line.L().I2()));

	bbox.Increase (filldist * bbox.Diam());
	loch2.CutBoundary (bbox);
      }
    
    loch2.FindInnerBoxes (adfront, NULL);

      // outer points : smooth mesh-grading
    npoints.SetSize(0);
    loch2.GetOuterPoints (npoints);
    
    for (int i = 1; i <= npoints.Size(); i++)
      {
	if (meshbox.IsIn (npoints.Get(i)))
	  {
	    PointIndex gpnum = mesh.AddPoint (npoints.Get(i));
	    adfront->AddPoint (npoints.Get(i), gpnum);
	  }
      }  

    NgProfiler::StopTimer (timer4);    
  }



  
  void Meshing2 :: Delaunay (Mesh & mesh, int domainnr, const MeshingParameters & mp)
  {
    static int timer = NgProfiler::CreateTimer ("Meshing2::Delaunay - total");
    static int timerstart = NgProfiler::CreateTimer ("Meshing2::Delaunay - start");
    static int timerfinish = NgProfiler::CreateTimer ("Meshing2::Delaunay - finish");
    static int timer1 = NgProfiler::CreateTimer ("Meshing2::Delaunay - incremental");
    static int timer1a = NgProfiler::CreateTimer ("Meshing2::Delaunay - incremental a");
    static int timer1b = NgProfiler::CreateTimer ("Meshing2::Delaunay - incremental b");
    static int timer1c = NgProfiler::CreateTimer ("Meshing2::Delaunay - incremental c");
    static int timer1d = NgProfiler::CreateTimer ("Meshing2::Delaunay - incremental d");
    NgProfiler::RegionTimer reg (timer);



    cout << "2D Delaunay meshing (in progress)" << endl;


    BlockFillLocalH (mesh, mp);

    NgProfiler::StartTimer (timerstart);

    // do the delaunay


    // face bounding box:
    Box<3> bbox (Box<3>::EMPTY_BOX);

    for (int i = 0; i < adfront->GetNFL(); i++)
      {
	const FrontLine & line = adfront->GetLine(i);
        bbox.Add (Point<3> (adfront->GetPoint (line.L()[0])));
        bbox.Add (Point<3> (adfront->GetPoint (line.L()[1])));
      }

    for (int i = 0; i < mesh.LockedPoints().Size(); i++)
      bbox.Add (mesh.Point (mesh.LockedPoints()[i]));

    cout << "bbox = " << bbox << endl;

    // external point
    Vec<3> vdiag = bbox.PMax()-bbox.PMin();
    
    auto old_points = mesh.Points().Range();
    DelaunayTrig startel;
    startel[0] = mesh.AddPoint (bbox.PMin() + Vec<3> (-8*vdiag(0), -8*vdiag(1), 0));
    startel[1] = mesh.AddPoint (bbox.PMin() + Vec<3> (+8*vdiag(0), -8*vdiag(1), 0));
    startel[2] = mesh.AddPoint (bbox.PMin() + Vec<3> (0, 8*vdiag(1), 0));

    Box<3> hbox;
    hbox.Set (mesh[startel[0]]);
    hbox.Add (mesh[startel[1]]);
    hbox.Add (mesh[startel[2]]);
    Point<3> hp = mesh[startel[0]];
    hp(2) = 1;  hbox.Add (hp);
    hp(2) = -1; hbox.Add (hp);
    Box3dTree searchtree(hbox);

    Array<DelaunayTrig> tempels;
    startel.CalcCenter (mesh);

    tempels.Append (startel);
    searchtree.Insert(startel.BoundingBox(), 0);

    Array<int> closeels;
    Array<int> intersecting;
    Array<INDEX_2> edges;




    // reorder points
    Array<PointIndex, PointIndex::BASE, PointIndex> mixed(old_points.Size());
    int prims[] = { 11, 13, 17, 19, 23, 29, 31, 37 };
    int prim;
  
    {
      int i = 0;
      while (old_points.Size() % prims[i] == 0) i++;
      prim = prims[i];
    }

    for (PointIndex pi : old_points)
      mixed[pi] = PointIndex ( (prim * pi) % old_points.Size() + PointIndex::BASE );
    
    NgProfiler::StopTimer (timerstart);
    NgProfiler::StartTimer (timer1);


    for (PointIndex i1 : old_points)
      {
        PointIndex i = mixed[i1];

        NgProfiler::StartTimer (timer1a);
        Point<3> newp = mesh[i];
        intersecting.SetSize(0);
        edges.SetSize(0);

        searchtree.GetIntersecting (newp, newp, closeels);
        // for (int jj = 0; jj < closeels.Size(); jj++)
        // for (int j = 0; j < tempels.Size(); j++)
        for (int j : closeels)
          {
            if (tempels[j][0] < 0) continue;
            Point<3> c = tempels[j].Center();
            double r2 = tempels[j].Radius2();

            bool inside = Dist2(mesh[i], c) < r2;
            if (inside) intersecting.Append (j);
          }

        NgProfiler::StopTimer (timer1a);
        NgProfiler::StartTimer (timer1b);

        // find outer edges
        for (auto j : intersecting)
          {
            const DelaunayTrig & trig = tempels[j];
            for (int k = 0; k < 3; k++)
              {
                int p1 = trig[k];
                int p2 = trig[(k+1)%3];
                INDEX_2 edge(p1,p2);
                edge.Sort();
                bool found = false;
                for (int l = 0; l < edges.Size(); l++)
                  if (edges[l] == edge)
                    {
                      edges.Delete(l); 
                      found = true;
                      break;
                    }
                if (!found) edges.Append (edge);
              }
          }

        NgProfiler::StopTimer (timer1b);
        NgProfiler::StartTimer (timer1c);

        /*
        for (int j = intersecting.Size()-1; j >= 0; j--)
          tempels.Delete (intersecting[j]);
        */
        for (int j : intersecting)
          {
            searchtree.DeleteElement (j);
            tempels[j][0] = -1;
            tempels[j][1] = -1;
            tempels[j][2] = -1;
          }

        NgProfiler::StopTimer (timer1c);
        NgProfiler::StartTimer (timer1d);

        for (auto edge : edges)
          {
            DelaunayTrig trig (edge[0], edge[1], i);
            trig.CalcCenter (mesh);
            tempels.Append (trig);
            searchtree.Insert(trig.BoundingBox(), tempels.Size()-1);
          }

        NgProfiler::StopTimer (timer1d);
      }

    NgProfiler::StopTimer (timer1);
    NgProfiler::StartTimer (timerfinish);

    for (DelaunayTrig & trig : tempels)
      {
        if (trig[0] < 0) continue;

        Point<3> c = Center (mesh[trig[0]], mesh[trig[1]], mesh[trig[2]]);
        if (!adfront->Inside (Point<2> (c(0),c(1)))) continue;

        Vec<3> n = Cross (mesh[trig[1]]-mesh[trig[0]], 
                          mesh[trig[2]]-mesh[trig[0]]);
        if (n(2) < 0) Swap (trig[1], trig[2]);

        Element2d el(trig[0], trig[1], trig[2]);
        el.SetIndex (domainnr);
        mesh.AddSurfaceElement (el);
      }

    for (PointIndex pi : mesh.Points().Range())
      *testout << pi << ": " << mesh[pi].Type() << endl;

    NgProfiler::StopTimer (timerfinish);
  }

}
