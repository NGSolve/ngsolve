/*

2d Spline curve for Mesh generator

*/

#include <meshing.hpp>
#include <geometry2d.hpp>

namespace netgen
{



  SplineGeometry2d :: ~SplineGeometry2d()
  {
    for ( int i = 0; i < bcnames.Size(); i++ )
      delete bcnames[i];
    for (int i=0; i<materials.Size(); i++)
      delete [] materials[i];
  }


  void SplineGeometry2d :: Load (const char * filename)
  {

    ifstream infile;
    Point<2> x;
    char buf[50];


    infile.open (filename);
  
    if ( ! infile.good() )
      throw NgException(string ("Input file '") + 
			string (filename) +
			string ("' not available!"));

    TestComment ( infile );
  
    infile >> buf;   // file recognition

    tensormeshing.SetSize(0);
    quadmeshing.SetSize(0);

    TestComment ( infile );
    if ( strcmp (buf, "splinecurves2dnew") == 0 )
      {
	LoadDataNew ( infile );
      }
    else if ( strcmp (buf, "splinecurves2dv2") == 0 )
      {
	LoadDataV2 ( infile );
      }
    else
      {
	LoadData(infile );
      }
    infile.close();
  }



  // herbert: fixed TestComment
  void SplineGeometry2d :: TestComment ( ifstream & infile )
  {
    bool comment = true;
    char ch;
    while ( comment == true && !infile.eof() ) {
      infile.get(ch);
      if ( ch == '#' ) { // skip comments
	while (  ch != '\n' && !infile.eof() ) {
	  infile.get(ch);
	}
      }
      else if ( ch == '\n' )  { // skip empty lines
	;
      }
      else if ( isspace(ch) ) { // skip whitespaces
	; 
      }
      else { // end of comment
	infile.putback(ch);
	comment = false;
      }
    }
    return;
  }




  void SplineGeometry2d :: LoadData ( ifstream & infile )
  {      
    enum { D = 2 };

    int nump, numseg, leftdom, rightdom;
    Point<D> x;
    int hi1, hi2, hi3;
    double hd;
    char buf[50], ch;

    materials.SetSize(0);
    maxh.SetSize(0);
    infile >> elto0;

    TestComment ( infile );

    infile >> nump;
    for (int i = 0; i < nump; i++)
      {
	TestComment ( infile );
	for(int j=0; j<D; j++)
	  infile >> x(j);
	infile >> hd;

	Flags flags;

	ch = 'a';
	// infile >> ch;
	do {
	  infile.get (ch);
	} while (isspace(ch) && ch != '\n');
	while (ch == '-')
	  {
	    char flag[100];
	    flag[0]='-';
	    infile >> (flag+1);
	    flags.SetCommandLineFlag (flag);
	    ch = 'a';
	    do {
	      infile.get (ch);
	    } while (isspace(ch) && ch != '\n');
	  }
    
	if (infile.good())
	  infile.putback (ch);

	geompoints.Append (GeomPoint<D>(x, hd));
	geompoints.Last().hpref = flags.GetDefineFlag ("hpref");
	geompoints.Last().hmax = 1e99;
      }

    PrintMessage (3, nump, " points loaded");
    TestComment ( infile );

    infile >> numseg;
    bcnames.SetSize(numseg);
    for ( int i = 0; i < numseg; i++ )
      bcnames[i] = 0; // "default";

    SplineSeg<D> * spline = 0;

    PrintMessage (3, numseg, " segments loaded");
    for (int i = 0; i < numseg; i++)
      {
	TestComment ( infile );
      
	infile >> leftdom >> rightdom;

	// cout << "add spline " << i << ", left = " << leftdom << ", right = " << rightdom << endl;
      
	infile >> buf;
	// type of spline segement
	if (strcmp (buf, "2") == 0)
	  { // a line
	    infile >> hi1 >> hi2;
	    spline = new LineSeg<D>(geompoints[hi1-1],
				    geompoints[hi2-1]);
	  }
	else if (strcmp (buf, "3") == 0)
	  { // a rational spline
	    infile >> hi1 >> hi2 >> hi3;
	    spline = new SplineSeg3<D> (geompoints[hi1-1],
					geompoints[hi2-1],
					geompoints[hi3-1]);
	  }
	else if (strcmp (buf, "4") == 0)
	  { // an arc
	    infile >> hi1 >> hi2 >> hi3;
	    spline = new CircleSeg<D> (geompoints[hi1-1],
				       geompoints[hi2-1],
				       geompoints[hi3-1]);
	    // 	  break;
	  }
	else if (strcmp (buf, "discretepoints") == 0)
	  {
	    int npts;
	    infile >> npts;
	    Array< Point<D> > pts(npts);
	    for (int j = 0; j < npts; j++)
	      for(int k=0; k<D; k++)
		infile >> pts[j](k);

	    spline = new DiscretePointsSeg<D> (pts);
	  }
    

	SplineSegExt * spex = new SplineSegExt (*spline);
	
	infile >> spex->reffak;
	spex -> leftdom = leftdom;
	spex -> rightdom = rightdom;
	spex -> hmax = 1e99;
	splines.Append (spex);


	Flags flags;
	ch = 'a';
	infile >> ch;
	while (ch == '-')
	  {
	    char flag[100];
	    flag[0]='-';
	    infile >> (flag+1);
	    flags.SetCommandLineFlag (flag);
	    ch = 'a';
	    infile >> ch;
	  }
    
	if (infile.good())
	  infile.putback (ch);
    
	spex->bc = int (flags.GetNumFlag ("bc", i+1));
	spex->hpref_left = int (flags.GetDefineFlag ("hpref")) || 
	  int (flags.GetDefineFlag ("hprefleft"));
	spex->hpref_right = int (flags.GetDefineFlag ("hpref")) || 
	  int (flags.GetDefineFlag ("hprefright"));
	spex->copyfrom = int (flags.GetNumFlag ("copy", -1));
	if ( flags.StringFlagDefined("bcname") )
	  {
	    int mybc = spex->bc-1;
	    delete bcnames[mybc];
	    bcnames[mybc] = new string (flags.GetStringFlag("bcname","") );
	  }
      }
  }




  void SplineGeometry2d :: LoadDataNew ( ifstream & infile )
  {
    enum { D = 2 };
    int nump, numseg, leftdom, rightdom;
    Point<D> x;
    int hi1, hi2, hi3;
    double hd;
    char buf[50], ch;
    int pointnr;


    TestComment ( infile );
    infile >> elto0;
    TestComment ( infile );
      
    infile >> nump;
    geompoints.SetSize(nump);
      
    for (int i = 0; i < nump; i++)
      {
	TestComment ( infile );
	infile >> pointnr;
	if ( pointnr > nump )
	  {
	    throw NgException(string ("Point number greater than total number of points") );
	  }
	for(int j=0; j<D; j++)
	  infile >> x(j);


	// hd is now optional, default 1
	//  infile >> hd;
	hd = 1;

	Flags flags;


	// get flags, 
	ch = 'a';
	// infile >> ch;
	do 
	  {

	    infile.get (ch);
	    // if another int-value, set refinement flag to this value
	    // (corresponding to old files)
	    if ( int (ch) >= 48 && int(ch) <= 57 )
	      {
		infile.putback(ch);
		infile >> hd;
		infile.get(ch);
	      }
	  } 
	while (isspace(ch) && ch != '\n');
	while (ch == '-')
	  {
	    char flag[100];
	    flag[0]='-';
	    infile >> (flag+1);
	    flags.SetCommandLineFlag (flag);
	    ch = 'a';
	    do {
	      infile.get (ch);
	    } while (isspace(ch) && ch != '\n');
	  }
    
	if (infile.good())
	  infile.putback (ch);

	if ( hd == 1 )
	  hd = flags.GetNumFlag ( "ref", 1.0);
	//       geompoints.Append (GeomPoint<D>(x, hd));
	geompoints[pointnr-1] = GeomPoint<D>(x, hd);
	geompoints[pointnr-1].hpref = flags.GetDefineFlag ("hpref");
      }

    TestComment ( infile );

    infile >> numseg;
    bcnames.SetSize(numseg);
    for ( int i = 0; i < numseg; i++ )
      bcnames[i] = 0;//new"default";

    SplineSeg<D> * spline = 0;
    for (int i = 0; i < numseg; i++)
      {
	TestComment ( infile );
      
	infile >> leftdom >> rightdom;

	// cout << "add spline " << i << ", left = " << leftdom << endl;

	infile >> buf;
	// type of spline segement
	if (strcmp (buf, "2") == 0)
	  { // a line
	    infile >> hi1 >> hi2;
	    spline = new LineSeg<D> (geompoints[hi1-1],
				     geompoints[hi2-1]);
	  }
	else if (strcmp (buf, "3") == 0)
	  { // a rational spline
	    infile >> hi1 >> hi2 >> hi3;
	    spline = new SplineSeg3<D> (geompoints[hi1-1],
					geompoints[hi2-1],
					geompoints[hi3-1]);
	  }
	else if (strcmp (buf, "4") == 0)
	  { // an arc
	    infile >> hi1 >> hi2 >> hi3;
	    spline = new CircleSeg<D> (geompoints[hi1-1],
				       geompoints[hi2-1],
				       geompoints[hi3-1]);
	    // 	  break;
	  }
	else if (strcmp (buf, "discretepoints") == 0)
	  {
	    int npts;
	    infile >> npts;
	    Array< Point<D> > pts(npts);
	    for (int j = 0; j < npts; j++)
	      for(int k=0; k<D; k++)
		infile >> pts[j](k);

	    spline = new DiscretePointsSeg<D> (pts);
	  }
    
	//      infile >> spline->reffak;

	SplineSegExt * spex = new SplineSegExt (*spline);

	spex -> leftdom = leftdom;
	spex -> rightdom = rightdom;
	splines.Append (spex);

	// hd is now optional, default 1
	//  infile >> hd;
	hd = 1;
	infile >> ch;
      
	// get refinement parameter, if it is there
	// infile.get (ch);
	// if another int-value, set refinement flag to this value
	// (corresponding to old files)
	if ( int (ch) >= 48 && int(ch) <= 57 )
	  {
	    infile.putback(ch);
	    infile >> hd;
	    infile >> ch ;
	  }
      
	Flags flags;
	while (ch == '-')
	  {
	    char flag[100];
	    flag[0]='-';
	    infile >> (flag+1);
	    flags.SetCommandLineFlag (flag);
	    ch = 'a';
	    infile >> ch;
	  }
    
	if (infile.good())
	  infile.putback (ch);
    
	spex->bc = int (flags.GetNumFlag ("bc", i+1));
	spex->hpref_left = int (flags.GetDefineFlag ("hpref")) || 
	  int (flags.GetDefineFlag ("hprefleft"));
	spex->hpref_right = int (flags.GetDefineFlag ("hpref")) || 
	  int (flags.GetDefineFlag ("hprefright"));
	spex->copyfrom = int (flags.GetNumFlag ("copy", -1));
	spex->reffak = flags.GetNumFlag ("ref", 1 );
	spex->hmax = flags.GetNumFlag ("maxh", 1e99 );

	if ( flags.StringFlagDefined("bcname") )
	  {
	    int mybc = spex->bc-1;
	    if ( bcnames[mybc] ) delete bcnames[mybc];
	    bcnames[mybc] = new string (flags.GetStringFlag("bcname","") );
	  }

	if ( hd != 1 )
	  spex->reffak = hd;
      }
    if ( !infile.good() )
      return;
    TestComment ( infile );
    int numdomains;
    int domainnr;
    char material[100];

    if ( !infile.good() ) 
      return;

    infile >> numdomains;
    materials.SetSize(numdomains) ;
    maxh.SetSize ( numdomains ) ;
    maxh = 1e99;

    TestComment ( infile );

    for ( int i=0; i<numdomains; i++)
      materials [ i ] = new char (100);

    for ( int i=0; i<numdomains && infile.good(); i++)
      {
	TestComment ( infile );
	infile >> domainnr;
	infile >> material;
	strcpy(materials[domainnr-1], material);

	Flags flags;
	ch = 'a';
	infile >> ch;
	while (ch == '-')
	  {
	    char flag[100];
	    flag[0]='-';
	    infile >> (flag+1);
	    flags.SetCommandLineFlag (flag);
	    ch = 'a';
	    infile >> ch;
	  }
    
	if (infile.good())
	  infile.putback (ch);
	 
	maxh[domainnr-1] = flags.GetNumFlag ( "maxh", 1e99);
      }
    return;
  }




  void SplineGeometry2d :: LoadDataV2 ( ifstream & infile )
  { 
    enum { D = 2 };
    // new parser by Astrid Sinwel
    
    PrintMessage (1, "Load 2D Geometry V2");
    int nump, leftdom, rightdom;
    Point<D> x;
    int hi1, hi2, hi3;
    double hd;
    char buf[50], ch;
    int pointnr;

    string keyword;

    Array < GeomPoint<D> > infilepoints (0);
    Array <int> pointnrs (0);
    nump = 0;
    int numdomains = 0;


    TestComment ( infile );
    // refinement factor
    infile >> elto0;
    TestComment ( infile );
      

    // test if next ch is a letter, i.e. new keyword starts
    bool ischar = false;

    while ( infile.good() )
      {
	infile >> keyword;

	ischar = false;

	if ( keyword == "points" )
	  {
	    PrintMessage (3, "load points");
	    infile.get(ch);
	    infile.putback(ch);

	    // test if ch is a letter
	    if ( int(ch) >= 65 && int(ch) <=90 )
	      ischar = true;
	    if ( int(ch) >= 97 && int(ch) <= 122 )
	      ischar = true;

	    while ( ! ischar )
	      {
		TestComment ( infile );
		infile >> pointnr;
		// pointnrs 1-based
		if ( pointnr > nump ) nump = pointnr; 
		pointnrs.Append(pointnr);
	      
		for(int j=0; j<D; j++)
		  infile >> x(j);
		// hd is now optional, default 1
		//  infile >> hd;
		hd = 1;
	      
		Flags flags;
	      
	      
		// get flags, 
		ch = 'a';
		// infile >> ch;
		do 
		  {
		    infile.get (ch);
		    // if another int-value, set refinement flag to this value
		    // (corresponding to old files)
		    if ( int (ch) >= 48 && int(ch) <= 57 )
		      {
			infile.putback(ch);
			infile >> hd;
			infile.get(ch);
		      }
		  } 
		while (isspace(ch) && ch != '\n');
		while (ch == '-')
		  {
		    char flag[100];
		    flag[0]='-';
		    infile >> (flag+1);
		    flags.SetCommandLineFlag (flag);
		    ch = 'a';
		    do {
		      infile.get (ch);
		    } while (isspace(ch) && ch != '\n');
		  }
		if (infile.good())
		  infile.putback (ch);
	      
		if ( hd == 1 )
		  hd = flags.GetNumFlag ( "ref", 1.0);
		//       geompoints.Append (GeomPoint<D>(x, hd));

		infilepoints.Append ( GeomPoint<D>(x, hd) );
		infilepoints.Last().hpref = flags.GetDefineFlag ("hpref");
		infilepoints.Last().hmax = flags.GetNumFlag ("maxh", 1e99);

		TestComment(infile);
		infile.get(ch);
		infile.putback(ch);

		// test if letter
		if ( int(ch) >= 65 && int(ch) <=90 )
		  ischar = true;
		if ( int(ch) >= 97 && int(ch) <= 122 )
		  ischar = true;
	      }

	    //	  infile.putback (ch);

	    geompoints.SetSize(nump);
	    for ( int i = 0; i < nump; i++ )
	      {
		geompoints[pointnrs[i] - 1] = infilepoints[i];
		geompoints[pointnrs[i] - 1].hpref = infilepoints[i].hpref; 
	      }
	    TestComment(infile);
	  }

	else if ( keyword == "segments" )
	  {
	    PrintMessage (3, "load segments");

	    bcnames.SetSize(0);
	    infile.get(ch);
	    infile.putback(ch);
	    int i = 0;

	    // test if ch is a letter
	    if ( int(ch) >= 65 && int(ch) <=90 )
	      ischar = true;
	    if ( int(ch) >= 97 && int(ch) <= 122 )
	      ischar = true;

	    while ( !ischar ) //ch != 'p' && ch != 'm' )
	      {
		i++;
		TestComment ( infile );

		SplineSeg<D> * spline = 0;
		TestComment ( infile );
		  
		infile >> leftdom >> rightdom;
	      
		if ( leftdom > numdomains ) numdomains = leftdom;
		if ( rightdom > numdomains ) numdomains = rightdom;

	      
		infile >> buf;
		// type of spline segement
		if (strcmp (buf, "2") == 0)
		  { // a line
		    infile >> hi1 >> hi2;
		    spline = new LineSeg<D>(geompoints[hi1-1],
					    geompoints[hi2-1]);
		  }
		else if (strcmp (buf, "3") == 0)
		  { // a rational spline
		    infile >> hi1 >> hi2 >> hi3;
		    spline = new SplineSeg3<D> (geompoints[hi1-1],
						geompoints[hi2-1],
						geompoints[hi3-1]);
		  }
		else if (strcmp (buf, "4") == 0)
		  { // an arc
		    infile >> hi1 >> hi2 >> hi3;
		    spline = new CircleSeg<D> (geompoints[hi1-1],
					       geompoints[hi2-1],
					       geompoints[hi3-1]);
		  }
		else if (strcmp (buf, "discretepoints") == 0)
		  {
		    int npts;
		    infile >> npts;
		    Array< Point<D> > pts(npts);
		    for (int j = 0; j < npts; j++)
		      for(int k=0; k<D; k++)
			infile >> pts[j](k);
		  
		    spline = new DiscretePointsSeg<D> (pts);
		  }
		else if (strcmp (buf, "bsplinepoints") == 0)
		  {
		    int npts,order;
		    infile >> npts;    
		    infile >> order;
		    Array< Point<D> > pts(npts);
		    for (int j = 0; j < npts; j++)
		      for(int k=0; k<D; k++)
			infile >> pts[j](k);	    		    
		    if(order<2)		      
			cerr<<"Minimum order of 2 is required!!"<<endl;
		    else if(order==2)
		      spline = new BSplineSeg<D,2> (pts);
		      else if(order==3)
			spline = new BSplineSeg<D,3> (pts);
		      else if(order==4)
			spline = new BSplineSeg<D,4> (pts);
		      else if(order>4)		      
			cerr<<"Maximum allowed order is 4!!"<<endl;
		  }
	      
		//      infile >> spline->reffak;
		SplineSegExt * spex = new SplineSegExt (*spline);

		spex -> leftdom = leftdom;
		spex -> rightdom = rightdom;
		splines.Append (spex);
	      
	      
		// hd is now optional, default 1
		//  infile >> hd;
		hd = 1;

		infile >> ch;
	      
		// get refinement parameter, if it is there
		//infile.get (ch);
		// if another int-value, set refinement flag to this value
		// (corresponding to old files)
                /*
		if ( int (ch) >= 48 && int(ch) <= 57 )
		  {
		    infile.putback(ch);
		    infile >> hd;
		    infile >> ch ;
		  }
                */
		// get flags, 
		Flags flags;
		while (ch == '-')
		  {
		    char flag[100];
		    flag[0]='-';
		    infile >> (flag+1);
		    flags.SetCommandLineFlag (flag);
		    ch = 'a';
		    infile >> ch;
		  }
	      
		if (infile.good())
		  infile.putback (ch);
	      
		spex->bc = int (flags.GetNumFlag ("bc", i+1));
		spex->hpref_left = int (flags.GetDefineFlag ("hpref")) || 
		  int (flags.GetDefineFlag ("hprefleft"));
		spex->hpref_right = int (flags.GetDefineFlag ("hpref")) || 
		  int (flags.GetDefineFlag ("hprefright"));
		spex->copyfrom = int (flags.GetNumFlag ("copy", -1));
		spex->reffak = flags.GetNumFlag ("ref", 1 );
		spex->hmax = flags.GetNumFlag ("maxh", 1e99 );
		if ( hd != 1 ) spex->reffak = hd;

		if ( flags.StringFlagDefined("bcname") )
		  {
		    int mybc = spex->bc-1;
		    for ( int ii = bcnames.Size(); ii <= mybc; ii++ )
		      bcnames.Append ( new string ("default"));
		    if ( bcnames[mybc] ) delete bcnames[mybc];
		    bcnames[mybc] = new string (flags.GetStringFlag("bcname","") );
		  }

		TestComment(infile);
		infile.get(ch);
		infile.putback(ch);

		// test if ch is a letter
		if ( int(ch) >= 65 && int(ch) <=90 )
		  ischar = true;
		if ( int(ch) >= 97 && int(ch) <= 122 )
		  ischar = true;

	      }
	  
	    infile.get(ch);
	    infile.putback(ch);
	

	  }
	else if ( keyword == "materials" )
	  {
	    TestComment ( infile );
	    int domainnr;
	    char material[100];
	  
	    if ( !infile.good() ) 
	      return;
	  
	    materials.SetSize(numdomains) ;
	    maxh.SetSize ( numdomains ) ;
	    for ( int i = 0; i < numdomains; i++)
	      maxh[i] = 1000;
	    quadmeshing.SetSize ( numdomains );
	    quadmeshing = false;
	    tensormeshing.SetSize ( numdomains );
	    tensormeshing = false;
	    layer.SetSize ( numdomains );
	    layer = 1;

	  
	    TestComment ( infile );
	  
	    for ( int i=0; i<numdomains; i++)
	      materials [ i ] = new char[100];
	  
	    for ( int i=0; i<numdomains && infile.good(); i++)
	      {
		TestComment ( infile );
		infile >> domainnr;
		infile >> material;

		strcpy (materials[domainnr-1], material);
	      
		Flags flags;
		ch = 'a';
		infile >> ch;
		while (ch == '-')
		  {
		    char flag[100];
		    flag[0]='-';
		    infile >> (flag+1);
		    flags.SetCommandLineFlag (flag);
		    ch = 'a';
		    infile >> ch;
		  }
	      
		if (infile.good())
		  infile.putback (ch);
	      
		maxh[domainnr-1] = flags.GetNumFlag ( "maxh", 1000);
		if (flags.GetDefineFlag("quad")) quadmeshing[domainnr-1] = true;
		if (flags.GetDefineFlag("tensor")) tensormeshing[domainnr-1] = true;
		layer[domainnr-1] = int(flags.GetNumFlag ("layer", 1));
	      }
	  }
      }
    return;
  }









  /*
  void CalcPartition (const SplineSegExt & spline,
		      double l, double h, double h1, double h2,
		      double hcurve, double elto0, Array<double> & points)
  {
    double fperel, oldf, f;

    int n = 1000;

    points.SetSize (0);

    double dt = l / n;

    double sum = 0;
    for (int i = 1; i <= n; i++)
      {
	double t = (i-0.5)*dt;
	double fun = min3 (hcurve, t/elto0 + h1, (l-t)/elto0 + h2);
	double curv = spline.CalcCurvature (t/l);
	cout << "curv = " << curv << endl;
	if (curv < 1e-10) curv = 1e-10;
	fun = min2 (fun, 0.1/curv);
	sum += dt / fun;
      }

    int nel = int (sum+1);
    fperel = sum / nel;

    points.Append (0);

    int i = 1;
    oldf = 0;
    // t = 0.5 * dt;
    for (int j = 1; j <= n && i < nel; j++)
      {
	double t = (j-0.5)*dt;
	double fun = min3 (hcurve, t/elto0 + h1, (l-t)/elto0 + h2);
	double curv = spline.CalcCurvature (t/l);
	if (curv < 1e-10) curv = 1e-10;
	fun = min2 (fun, 0.1/curv);

	f = oldf + dt / fun;

	while (i * fperel < f && i < nel)
	  {
	    points.Append ( dt * (j-1) +  (i * fperel - oldf) * fun);
	    i++;
	  }
	oldf = f;
	t += dt;
      }
    points.Append (l);
  }
  */







  string SplineGeometry2d :: GetBCName( const int  bcnr ) const
  {
    if ( bcnames.Size() >= bcnr)
      if (bcnames[bcnr-1] )
	return *bcnames[bcnr-1];
    return "default";
  }

  string * SplineGeometry2d :: BCNamePtr( const int bcnr ) 
  {
    if ( bcnr > bcnames.Size() )
      return 0;
    else
      return bcnames[bcnr-1];
  }





  void SplineGeometry2d :: GetMaterial (int  domnr, char* & material )
  {
    if ( materials.Size() >= domnr)
      material =  materials[domnr-1];
    else 
      material = 0;
  }

  void SplineGeometry2d :: SetMaterial (int  domnr, const string & material)
  {
    int oldsize = materials.Size();
    if (domnr > materials.Size()) materials.SetSize (domnr);
    for (int i = oldsize; i < domnr; i++)
      materials[i] = nullptr;
    
    if (domnr >= 1) //  && domnr <= materials.Size())
      {
        delete materials[domnr-1];
        materials[domnr-1] = new char[material.size()+1];
        strcpy(materials[domnr-1], material.c_str());
      }
    else
      throw NgException ("material index out of range");
  }

  

  double SplineGeometry2d :: GetDomainMaxh (const int domnr )
  {
    if ( maxh.Size() >= domnr  && domnr > 0)
      return maxh[domnr-1];
    else
      return -1;
  }

  void SplineGeometry2d :: SetDomainMaxh (int domnr, double h)
  {
    int oldsize = maxh.Size();
    if (domnr > maxh.Size()) maxh.SetSize (domnr);
    for (int i = oldsize; i < domnr; i++)
      maxh[i] = 1e99;
    
    if (domnr >= 1) 
      maxh[domnr] = h;
    else
      throw NgException ("material index out of range");
  }

  

  extern void MeshFromSpline2D (SplineGeometry2d & geometry,
				shared_ptr<Mesh> & mesh, 
				MeshingParameters & mp);


  int SplineGeometry2d :: GenerateMesh (shared_ptr<Mesh> & mesh, MeshingParameters & mparam,
					int perfstepsstart, int perfstepsend)
  {
    MeshFromSpline2D (*this, mesh, mparam);
    return 0;
  }


  Refinement & SplineGeometry2d :: GetRefinement () const
  {
    return * new Refinement2d (*this);
  }



  class SplineGeometryRegister : public GeometryRegister
  {
  public:
    virtual NetgenGeometry * Load (string filename) const;
  };

  NetgenGeometry *  SplineGeometryRegister :: Load (string filename) const
  {
    const char * cfilename = filename.c_str();
    if (strcmp (&cfilename[strlen(cfilename)-4], "in2d") == 0)
      {
	PrintMessage (1, "Load 2D-Spline geometry file ", cfilename);
	

	ifstream infile(cfilename);

	SplineGeometry2d * hgeom = new SplineGeometry2d();
	hgeom -> Load (cfilename);
	return hgeom;
      }
    
    return NULL;
  }

  class SplineGeoInit
  {
  public:
    SplineGeoInit()
    {
      geometryregister.Append (new SplineGeometryRegister);
    }
  };

  SplineGeoInit sginit;

}
