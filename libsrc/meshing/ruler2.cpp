#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{


  static double CalcElementBadness (const Array<Point2d> & points,
				    const Element2d & elem)
  {
    // badness = sqrt(3) /36 * circumference^2 / area - 1 +
    //           h / li + li / h - 2

    Vec2d v12, v13, v23;
    double l12, l13, l23, cir, area;
    static const double c = sqrt(3.0) / 36;

    v12 = points.Get(elem.PNum(2)) - points.Get(elem.PNum(1));
    v13 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(1));
    v23 = points.Get(elem.PNum(3)) - points.Get(elem.PNum(2));

    l12 = v12.Length();
    l13 = v13.Length();
    l23 = v23.Length();

    cir = l12 + l13 + l23;
    area = 0.5 * (v12.X() * v13.Y() - v12.Y() * v13.X());
    if (area < 1e-6)
      {
	return 1e8;
      }

    if (testmode)
      {
	(*testout) << "l = " << l12 << " + " << l13 << " + " << l23 << " = " 
		   << cir << ", area = " << area << endl;
	(*testout) << "shapeerr = " << 10 * (c * cir * cir / area - 1) << endl
		   << "sizeerr = " << 1/l12 + l12 + 1/l13 + l13 + 1/l23 + l23 - 6
		   << endl;
      }

    return 10 * (c * cir * cir / area - 1)
      + 1/l12 + l12 + 1/l13 + l13 + 1/l23 + l23 - 6;
  }



  int Meshing2 ::ApplyRules (Array<Point2d> & lpoints, 
			     Array<int> & legalpoints,
			     int maxlegalpoint,
			     Array<INDEX_2> & llines1,
			     int maxlegalline,
			     Array<Element2d> & elements,
			     Array<INDEX> & dellines, int tolerance,
			     const MeshingParameters & mp)
  {
    static int timer = NgProfiler::CreateTimer ("meshing2::ApplyRules");
    NgProfiler::RegionTimer reg (timer);



    double maxerr = 0.5 + 0.3 * tolerance;
    double minelerr = 2 + 0.5 * tolerance * tolerance;

    int noldlp = lpoints.Size();
    int noldll = llines1.Size();


    ArrayMem<int,100> pused(maxlegalpoint), lused(maxlegalline);
    ArrayMem<int,100> pnearness(noldlp), lnearness(llines1.Size());

    ArrayMem<int, 20> pmap, pfixed, lmap;
  
    ArrayMem<Point2d,100> tempnewpoints;
    ArrayMem<INDEX_2,100> tempnewlines;
    ArrayMem<int,100> tempdellines;
    ArrayMem<Element2d,100> tempelements;


    elements.SetSize (0);
    dellines.SetSize (0);

    testmode = debugparam.debugoutput;

#ifdef LOCDEBUG
    int loctestmode = testmode;

    if (loctestmode)
      {
	(*testout) << endl << endl << "Check new environment" << endl;
	(*testout) << "tolerance = " << tolerance << endl;
	for (int i = 1; i <= lpoints.Size(); i++)
	  (*testout) << "P" << i << " = " << lpoints.Get(i) << endl;
	(*testout) << endl;
	for (int i = 1; i <= llines1.Size(); i++)
	  (*testout) << "(" << llines1.Get(i).I1() << "-" << llines1.Get(i).I2() << ")" << endl;
      }
#endif

    // check every rule

    int found = 0;   // rule number

    pnearness = 1000;
  
    for (int j = 0; j < 2; j++)
      pnearness.Set(llines1[0][j], 0);



    enum { MAX_NEARNESS = 3 };

    for (int cnt = 0; cnt < MAX_NEARNESS; cnt++)
      {
	bool ok = true;
	for (int i = 0; i < maxlegalline; i++)
	  {
	    const INDEX_2 & hline = llines1[i];

	    int minn = min2 (pnearness.Get(hline[0]),  pnearness.Get(hline[1]));

	    for (int j = 0; j < 2; j++)
	      if (pnearness.Get(hline[j]) > minn+1)
		{
		  ok = false;
		  pnearness.Set(hline[j], minn+1);
		}
	  }
	if (!ok) break;
      }


    for (int i = 0; i < maxlegalline; i++)
      lnearness[i] = pnearness.Get(llines1[i][0]) + pnearness.Get(llines1[i][1]);


    // resort lines after lnearness
    Array<INDEX_2> llines(llines1.Size());
    Array<int> sortlines(llines1.Size());
    int lnearness_class[MAX_NEARNESS];

    for (int j = 0; j < MAX_NEARNESS; j++)
      lnearness_class[j] = 0;
    for (int i = 0; i < maxlegalline; i++)
      if (lnearness[i] < MAX_NEARNESS)
	lnearness_class[lnearness[i]]++;
    
    int cumm = 0;
    for (int j = 0; j < MAX_NEARNESS; j++)
      {
	int hcnt = lnearness_class[j];
	lnearness_class[j] = cumm;
	cumm += hcnt;
      }

    for (int i = 0; i < maxlegalline; i++)
      if (lnearness[i] < MAX_NEARNESS)
	{
	  llines[lnearness_class[lnearness[i]]] = llines1[i];
	  sortlines[lnearness_class[lnearness[i]]] = i+1;
	  lnearness_class[lnearness[i]]++;
	}
      else
	{
	  llines[cumm] = llines1[i];
	  sortlines[cumm] = i+1;
	  cumm++;
	}

    for (int i = maxlegalline; i < llines1.Size(); i++)
      {
	llines[cumm] = llines1[i];
	sortlines[cumm] = i+1;
	cumm++;
      }

    for (int i = 0; i < maxlegalline; i++)
      lnearness[i] = pnearness.Get(llines[i][0]) + pnearness.Get(llines[i][1]);




    static bool firsttime = true;
    static int timers[100];
    static int timers2[100];
    static int timers3[100];
    if (firsttime)
      {
	for (int ri = 0; ri < rules.Size(); ri++)
	  timers[ri] = NgProfiler::CreateTimer (string("netrule ")+rules[ri]->Name());
	for (int ri = 0; ri < rules.Size(); ri++)
	  timers2[ri] = NgProfiler::CreateTimer (string("netrule,mapped ")+rules[ri]->Name());
	for (int ri = 0; ri < rules.Size(); ri++)
	  timers3[ri] = NgProfiler::CreateTimer (string("netrule,lines mapped ")+rules[ri]->Name());
	firsttime = false;
      }

    lused = 0;
    pused = 0;


    static int timer1 = NgProfiler::CreateTimer ("meshing2::ApplyRules 1");
    NgProfiler::RegionTimer reg1 (timer1);


    for (int ri = 1; ri <= rules.Size(); ri++)
      {
	NgProfiler::RegionTimer reg(timers[ri-1]);
	netrule * rule = rules.Get(ri);

#ifdef LOCDEBUG
	if (loctestmode)
	  (*testout) << "Rule " << rule->Name() << endl;
#endif

	if (rule->GetQuality() > tolerance) continue;

	pmap.SetSize (rule->GetNP());
	lmap.SetSize (rule->GetNL());
      
	pmap = 0;
	lmap = 0;

	lused[0] = 1; 
	lmap[0] = 1;  

	for (int j = 0; j < 2; j++)
	  {
	    pmap.Elem(rule->GetLine(1)[j]) = llines[0][j];
	    pused.Elem(llines[0][j])++;
	  }



	int nlok = 2;


	bool ok = false;

	while (nlok >= 2)
	  {

	    if (nlok <= rule->GetNOldL())

	      {
		ok = 0;
		
		int maxline = (rule->GetLNearness(nlok) < MAX_NEARNESS) ? lnearness_class[rule->GetLNearness(nlok)] : maxlegalline;
		// int maxline = maxlegalline;

		while (!ok && lmap.Get(nlok) < maxline)
		  {
		    lmap.Elem(nlok)++;
		    int locli = lmap.Get(nlok);

		    if (lnearness.Get(locli) > rule->GetLNearness (nlok) ) continue;
		    if (lused.Get(locli)) continue;


		    ok = 1;

		    INDEX_2 loclin = llines.Get(locli);
		    Vec2d linevec = lpoints.Get(loclin.I2()) - lpoints.Get(loclin.I1());

		    if (rule->CalcLineError (nlok, linevec) > maxerr)
		      {
			ok = 0;
#ifdef LOCDEBUG
			if(loctestmode)
			  (*testout) << "not ok pos1" << endl;
#endif
			continue;
		      }

		    for (int j = 0; j < 2; j++)
		      {
			int refpi = rule->GetLine(nlok)[j];

			if (pmap.Get(refpi) != 0)
			  {
			    if (pmap.Get(refpi) != loclin[j])
			      {
				ok = 0;
#ifdef LOCDEBUG
				if(loctestmode)
				  (*testout) << "not ok pos2" << endl;
#endif
				break;
			      }
			  }
			else
			  {
			    if (rule->CalcPointDist (refpi, lpoints.Get(loclin[j])) > maxerr
				|| !legalpoints.Get(loclin[j])
				|| pused.Get(loclin[j]))
			      {
				ok = 0;
#ifdef LOCDEBUG
				if(loctestmode)
				  {
				    (*testout) << "nok pos3" << endl;
				    //if(rule->CalcPointDist (refpi, lpoints.Get(loclin[j])) > maxerr)
				    //(*testout) << "r1" << endl;
				    //if(!legalpoints.Get(loclin[j]))
				    //(*testout) << "r2 legalpoints " << legalpoints << " loclin " << loclin << " j " << j << endl;
				    //if(pused.Get(loclin[j]))
				    //(*testout) << "r3" << endl;
				  }
#endif
				break;
			      }
			  }
		      }
		  }

		if (ok)
		  {
		    int locli = lmap.Get(nlok);
		    INDEX_2 loclin = llines.Get(locli);

		    lused.Elem (locli) = 1;
		    for (int j = 0; j < 2; j++)
		      {
			pmap.Set(rule->GetLine (nlok)[j], loclin[j]);
			pused.Elem(loclin[j])++;
		      }

		    nlok++;
		  }
		else
		  {
		    lmap.Elem(nlok) = 0;
		    nlok--;

		    lused.Elem (lmap.Get(nlok)) = 0;
		    for (int j = 0; j < 2; j++)
		      {
			pused.Elem(llines.Get(lmap.Get(nlok))[j]) --;
			if (! pused.Get (llines.Get (lmap.Get (nlok))[j]))
			  pmap.Set (rule->GetLine (nlok)[j], 0);
		      }
		  }
	      }

	    else

	      {
		NgProfiler::RegionTimer reg(timers3[ri-1]);

		// all lines are mapped !!

		// map also all points:

		int npok = 1;
		int incnpok = 1;

		pfixed.SetSize (pmap.Size());
		for (int i = 0; i < pmap.Size(); i++)
		  pfixed[i] = (pmap[i] >= 1);
 
		while (npok >= 1)
		  {

		    if (npok <= rule->GetNOldP())

		      {
			if (pfixed.Get(npok))

			  {
			    if (incnpok)
			      npok++;
			    else
			      npok--;
			  }

			else

			  {
			    ok = 0;

			    if (pmap.Get(npok))
			      pused.Elem(pmap.Get(npok))--;

			    while (!ok && pmap.Get(npok) < maxlegalpoint)
			      {
				ok = 1;

				pmap.Elem(npok)++;

				if (pused.Get(pmap.Get(npok)))
				  {
				    ok = 0;
				  }
				else
				  {
				    if (rule->CalcPointDist (npok, lpoints.Get(pmap.Get(npok))) > maxerr 
					|| !legalpoints.Get(pmap.Get(npok))) 
                                    
				      ok = 0;
				  }
			      }

			    if (ok)
			      {
				pused.Elem(pmap.Get(npok))++;
				npok++;
				incnpok = 1;
			      }

			    else

			      {
				pmap.Elem(npok) = 0;
				npok--;
				incnpok = 0;
			      }
			  }
		      }

		    else

		      {
			NgProfiler::RegionTimer reg(timers2[ri-1]);

			npok = rule->GetNOldP();
			incnpok = 0;

			if (ok)
			  foundmap.Elem(ri)++; 

#ifdef LOCDEBUG
			if (loctestmode)
			  (*testout) << "lines and points mapped" << endl;
#endif

			ok = 1;

			// check orientations

			for (int i = 1; i <= rule->GetNOrientations(); i++)
			  {
			    if (CW (lpoints.Get(pmap.Get(rule->GetOrientation(i).i1)),
				    lpoints.Get(pmap.Get(rule->GetOrientation(i).i2)),
				    lpoints.Get(pmap.Get(rule->GetOrientation(i).i3))) )
			      {
				ok = 0;
#ifdef LOCDEBUG
				if (loctestmode)
				  (*testout) << "Orientation " << i << " not ok" << endl;
#endif
				break;
			      }
			  }


			if (!ok) continue;

			Vector oldu (2 * rule->GetNOldP());
		      
			for (int i = 1; i <= rule->GetNOldP(); i++)
			  {
			    Vec2d ui(rule->GetPoint(i), lpoints.Get(pmap.Get(i)));
			    oldu (2*i-2) = ui.X();
			    oldu (2*i-1) = ui.Y();
			  }
		      
			rule -> SetFreeZoneTransformation (oldu, tolerance);

		      
			if (!ok) continue;
			if (!rule->ConvexFreeZone())
			  {
			    ok = 0;
#ifdef LOCDEBUG
			    if (loctestmode) 
			      (*testout) << "freezone not convex" << endl;
#endif
			    /*
			      static int cnt = 0;
			      cnt++;
			      if (cnt % 100 == 0)
			      {
			      cout << "freezone not convex, cnt = " << cnt << "; rule = " << rule->Name() << endl;
			      (*testout) << "freezone not convex, cnt = " << cnt << "; rule = " << rule->Name() << endl;
			      (*testout) << "tol = " << tolerance << endl;
			      (*testout) << "maxerr = " << maxerr << "; minerr = " << minelerr << endl;
			      (*testout) << "freezone = " << rule->GetTransFreeZone() << endl;
			      }
			    */
			  }

			// check freezone:
			if (!ok) continue;
			for (int i = 1; i <= maxlegalpoint && ok; i++)
			  {
			    if ( !pused.Get(i) &&
				 rule->IsInFreeZone (lpoints.Get(i)) )
			      {
				ok = 0;
#ifdef LOCDEBUG
				if (loctestmode)
				  (*testout) << "Point " << i << " in freezone" << endl;
#endif
				break;
			      }
			  }

			if (!ok) continue;
			for (int i = maxlegalpoint+1; i <= lpoints.Size(); i++)
			  {
			    if ( rule->IsInFreeZone (lpoints.Get(i)) )
			      {
				ok = 0;
#ifdef LOCDEBUG
				if (loctestmode)
				  (*testout) << "Point " << i << " in freezone" << endl;
#endif
				break;
			      }
			  }


			if (!ok) continue;
			for (int i = 1; i <= maxlegalline; i++)
			  {
			    if (!lused.Get(i) && 
				rule->IsLineInFreeZone (lpoints.Get(llines.Get(i).I1()),
							lpoints.Get(llines.Get(i).I2())))
			      {
				ok = 0;
#ifdef LOCDEBUG
				if (loctestmode)
				  (*testout) << "line " << llines.Get(i).I1() << "-"
					     << llines.Get(i).I2() << " in freezone" << endl;
#endif
				break;
			      }
			  }

			if (!ok) continue;

			for (int i = maxlegalline+1; i <= llines.Size(); i++)
			  {
			    if (rule->IsLineInFreeZone (lpoints.Get(llines.Get(i).I1()),
							lpoints.Get(llines.Get(i).I2())))
			      {
				ok = 0;
#ifdef LOCDEBUG
				if (loctestmode)
				  (*testout) << "line " << llines.Get(i).I1() << "-"
					     << llines.Get(i).I2() << " in freezone" << endl;
#endif
				break;
			      }
			  }


			/*
			// check orientations

			for (i = 1; i <= rule->GetNOrientations() && ok; i++)
			{
			if (CW (lpoints.Get(pmap.Get(rule->GetOrientation(i).i1)),
			lpoints.Get(pmap.Get(rule->GetOrientation(i).i2)),
			lpoints.Get(pmap.Get(rule->GetOrientation(i).i3))) )
			{
			ok = 0;
			if (loctestmode)
			(*testout) << "Orientation " << i << " not ok" << endl;
			}
			}
			*/


			if (!ok) continue;

#ifdef LOCDEBUG
			if (loctestmode)
			  (*testout) << "rule ok" << endl;
#endif

			// Setze neue Punkte:
			if (rule->GetNOldP() < rule->GetNP())
			  {
			    Vector newu(rule->GetOldUToNewU().Height());
			    rule->GetOldUToNewU().Mult (oldu, newu);
			    
			    int oldnp = rule->GetNOldP();
			    for (int i = oldnp + 1; i <= rule->GetNP(); i++)
			      {
				Point2d np = rule->GetPoint(i);
				np.X() += newu (2 * (i-oldnp) - 2);
				np.Y() += newu (2 * (i-oldnp) - 1);
				
				pmap.Elem(i) = lpoints.Append (np);
			      }
			  }

			// Setze neue Linien:

			for (int i = rule->GetNOldL() + 1; i <= rule->GetNL(); i++)
			  {
			    llines.Append (INDEX_2 (pmap.Get(rule->GetLine (i)[0]),
						    pmap.Get(rule->GetLine (i)[1])));
			  }


			// delete old lines:
			for (int i = 1; i <= rule->GetNDelL(); i++)
			  dellines.Append (sortlines.Elem (lmap.Get(rule->GetDelLine(i))));
			// dellines.Append (lmap.Get(rule->GetDelLine(i))));

			// dellines.Append (lmap.Elem(rule->GetDelLines()));
			// lmap[rule->GetDelLines()];


			// insert new elements:

			for (int i = 1; i <= rule->GetNE(); i++)
			  {
			    elements.Append (rule->GetElement(i));
			    for (int j = 1; j <= elements.Get(i).GetNP(); j++)
			      elements.Elem(i).PNum(j) = pmap.Get(elements.Get(i).PNum(j));
			  }


			double elerr = 0;
			for (int i = 1; i <= elements.Size(); i++)
			  {
			    double hf;
			    if (!mp.quad)
			      hf = CalcElementBadness (lpoints, elements.Get(i));
			    else
			      hf = elements.Get(i).CalcJacobianBadness (lpoints) * 5;
#ifdef LOCDEBUG
			    if (loctestmode)
			      (*testout) << "r " << rule->Name() << "bad = " << hf << endl;
#endif
			    if (hf > elerr) elerr = hf;
			  }

#ifdef LOCDEBUG
			if (loctestmode)
			  (*testout) << "error = " << elerr;
#endif

			canuse.Elem(ri) ++;

			if (elerr < 0.99*minelerr)
			  {
#ifdef LOCDEBUG
			    if (loctestmode)
			      {
				(*testout) << "rule = " << rule->Name() << endl;
				(*testout) << "class = " << tolerance << endl;
				(*testout) << "lpoints: " << endl;
				for (int i = 1; i <= lpoints.Size(); i++)
				  (*testout) << lpoints.Get(i) << endl;
				(*testout) << "llines: " << endl;
				for (int i = 1; i <= llines.Size(); i++)
				  (*testout) << llines.Get(i).I1() << " " << llines.Get(i).I2() << endl;

				(*testout) << "Freezone: ";
				for (int i = 1; i <= rule -> GetTransFreeZone().Size(); i++)
				  (*testout) << rule->GetTransFreeZone().Get(i) << endl;
			      }
#endif

			    minelerr = elerr;
			    found = ri;

			    tempnewpoints = lpoints.Range (noldlp, lpoints.Size());
			    tempnewlines = llines.Range (noldll, llines.Size());
			    tempdellines = dellines;
			    tempelements = elements;
			  }

			lpoints.SetSize (noldlp);
			llines.SetSize (noldll);
			dellines.SetSize (0);
			elements.SetSize (0);
			ok = 0;
		      }
		  }

		nlok = rule->GetNOldL();

		lused.Set (lmap.Get(nlok), 0);

		for (int j = 1; j <= 2; j++)
		  {
		    int refpi = rule->GetPointNr (nlok, j);
		    pused.Elem(pmap.Get(refpi))--;

		    if (pused.Get(pmap.Get(refpi)) == 0)
		      pmap.Set(refpi, 0);
		  }
	      }
	  }
      }


    if (found)
      {
	lpoints.Append (tempnewpoints);
	llines1.Append (tempnewlines);
	dellines.Append (tempdellines);
	elements.Append (tempelements);
      }


    return found;
  }





}
