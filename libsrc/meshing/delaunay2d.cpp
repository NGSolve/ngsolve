#include <mystdlib.h>
#include "meshing.hpp"

// not yet working ....

namespace netgen
{

  void Meshing2 :: BlockFillLocalH (Mesh & mesh, const MeshingParameters & mp)
  {
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
	
	double hi = Dist (p1, p2);
	if (hi > maxh) maxh = hi;
	
	bbox.Add (p1);
	bbox.Add (p2);
      }

    
    cout << "bbox = " << bbox << endl;


    // Point<3> mpc = bbox.Center();
    bbox.Increase (bbox.Diam()/2);
    Box<3> meshbox = bbox;
    
    LocalH loch2 (bbox, 1);
    
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
	    if (mesh.LocalHFunction().GetH(npoints[i]) > 1.5 * maxh)
	      {
		mesh.LocalHFunction().SetH (npoints[i], maxh);
		changed = true;
	      }
	  }
      }
    while (changed);

    if (debugparam.slowchecks)
      (*testout) << "Blockfill with points: " << endl;
    *testout << "loch = " << mesh.LocalHFunction() << endl;
    
    *testout << "npoints = " << endl << npoints << endl;

    for (int i = 1; i <= npoints.Size(); i++)
      {
	if (meshbox.IsIn (npoints.Get(i)))
	  {
	    int gpnum = mesh.AddPoint (npoints.Get(i));
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
    
    npoints.SetSize(0);
    loch2.GetOuterPoints (npoints);
    
    for (int i = 1; i <= npoints.Size(); i++)
      {
	if (meshbox.IsIn (npoints.Get(i)))
	  {
	    int gpnum = mesh.AddPoint (npoints.Get(i));
	    adfront->AddPoint (npoints.Get(i), gpnum);
	  }
      }  

  }
  
  void Meshing2 :: Delaunay (Mesh & mesh, int domainnr, const MeshingParameters & mp)
  {
    cout << "2D Delaunay meshing (in progress)" << endl;

    // int oldnp = mesh.GetNP();

    cout << "np, old = " << mesh.GetNP() << endl;

    BlockFillLocalH (mesh, mp);


    cout << "np, now = " << mesh.GetNP() << endl;
    
  }

}
