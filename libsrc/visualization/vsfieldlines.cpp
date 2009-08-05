#ifndef NOTCL

#include <mystdlib.h>
#include "incvis.hpp"


#include <myadt.hpp>
#include <meshing.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>

#include <visual.hpp>


namespace netgen
{

  extern AutoPtr<Mesh> mesh;


  
  RKStepper :: ~RKStepper() 
  {
    delete a;
  }
    
  RKStepper :: RKStepper(int type) : a(NULL), tolerance(1e100)
  {
    notrestarted = 0;

    if (type == 0) // explicit Euler
      {
	c.SetSize(1); c[0] = 0;
	b.SetSize(1); b[0] = 1;
	steps = order = 1;
      }
    else if (type == 1) // Euler-Cauchy
      {
	c.SetSize(2); c[0] = 0; c[1] = 0.5;
	b.SetSize(2); b[0] = 0; b[1] = 1;
	Array<int> size(2);
	size[0] = 0; size[1] = 1;
	a = new TABLE<double>(size);
	a->Set(2,1,0.5);  // Set, Get: 1-based!
	steps = order = 2;
      }
    else if (type == 2) // Simpson
      {
	c.SetSize(3); c[0] = 0; c[1] = 1; c[2] = 0.5;
	b.SetSize(3); b[0] = b[1] = 1./6.; b[2] = 2./3.;
	Array<int> size(3);
	size[0] = 0; size[1] = 1; size[2] = 2;
	a = new TABLE<double>(size);
	a->Set(2,1,1);
	a->Set(3,1,0.25); a->Set(3,2,0.25); 
	steps = order = 3;
      }
    else if (type == 3) // classical Runge-Kutta
      {
	c.SetSize(4); c[0] = 0; c[1] = c[2] = 0.5; c[3] = 1;
	b.SetSize(4); b[0] = b[3] = 1./6.; b[1] = b[2] = 1./3.;
	Array<int> size(4);
	size[0] = 0; size[1] = 1; size[2] = 2; size[3] = 3;
	a = new TABLE<double>(size);
	a->Set(2,1,0.5);
	a->Set(3,1,0); a->Set(3,2,0.5); 
	a->Set(4,1,0); a->Set(4,2,0); a->Set(4,3,1); 
	steps = order = 4;
      }
    
    K.SetSize(steps);
  }

  void RKStepper :: StartNextValCalc(const Point3d & astartval, const double astartt, const double ah, const bool aadaptive)
  {
    //cout << "Starting RK-Step with h=" << ah << endl;

    stepcount = 0;
    h = ah;
    startt = astartt;
    startval = astartval;
    adaptive = aadaptive;
    adrun = 0;
  }

  bool RKStepper :: GetNextData(Point3d & val, double & t, double & ah)
  {
    bool finished(false);
    
    
    //cout << "stepcount " << stepcount << endl;
    
    if(stepcount <= steps)
      {
	t = startt + c[stepcount-1]*h;
	val = startval;
	for(int i=0; i<stepcount-1; i++)
	  val += h * a->Get(stepcount,i+1) * K[i];
      }
    
    
    if(stepcount == steps)
      {
	val = startval;
	for(int i=0; i<steps; i++)
	  val += h * b[i] * K[i];
	
	if(adaptive)
	  {
	    if(adrun == 0)
	      {
		stepcount = 0;
		h *= 0.5;
		adrun = 1;
		valh = val;
	      }
	    else if (adrun == 1)
	      {
		stepcount = 0;
		startval_bak = startval;
		startval = val;
		startt_bak = startt;
		startt += h;//0.5*h;
		adrun = 2;
	      }
	    else if (adrun == 2)
	      {
		Point3d valh2 = val;
		val = valh2 + 1./(pow(2.,order)-1.) * (valh2 - valh);
		Vec3d errvec = val - valh;
		
		double err = errvec.Length();
		
		double fac = 0.7 * pow(tolerance/err,1./(order+1.));
		if(fac > 1.3) fac = 1.3;
		
		if(fac < 1 || notrestarted >= 2)
		  ah = 2.*h * fac;
		
		if(err < tolerance) 
		  {
		    finished = true;
		    notrestarted++;
		    //(*testout) << "finished RK-Step, new h=" << ah << " tolerance " << tolerance << " err " << err << endl;
		  }
		else
		  {
		    //ah *= 0.9;
		    notrestarted = 0;
		    //(*testout) << "restarting h " << 2.*h << " ah " << ah << " tolerance " << tolerance << " err " << err << endl;
		    StartNextValCalc(startval_bak,startt_bak, ah, adaptive);
		  }
	      }
	  }
	else 
	  {
	    t = startt + h;
	    finished = true;
	  }
	
      }
    
    if(stepcount == 0)
      {
	t = startt + c[stepcount]*h;
	val = startval;
	for(int i=0; i<stepcount; i++)
	  val += h * a->Get(stepcount,i) * K[i];
      }
    
    return finished;
  }


  bool RKStepper :: FeedNextF(const Vec3d & f)
  {
    K[stepcount] = f;
    stepcount++;
    return true;
  }
  


  void FieldLineCalc :: GenerateFieldLines(Array<Point3d> & potential_startpoints, const int numlines, const int gllist,
					   const double minval, const double maxval, const int logscale, double phaser, double phasei)
  {

    
    Array<Point3d> points;
    Array<double> values;
    Array<bool> drawelems;
    Array<int> dirstart;


    if(vsol -> iscomplex)
      SetPhase(phaser,phasei);

    double crit = 1.0;

    if(randomized)
      {
	double sum = 0;
	double lami[3];
	double values[6];
	Vec3d v;
	
	for(int i=0; i<potential_startpoints.Size(); i++)
	  {
	    int elnr = mesh.GetElementOfPoint(potential_startpoints[i],lami,true) - 1;
	    if(elnr == -1)
	      continue;

	    mesh.SetPointSearchStartElement(elnr);
	    
	    if (mesh.GetDimension()==3)
              vss.GetValues ( vsol, elnr, lami[0], lami[1], lami[2], values);
            else
              vss.GetSurfValues ( vsol, elnr, lami[0], lami[1], values);
              
	    
	    VisualSceneSolution::RealVec3d ( values, v, vsol->iscomplex, phaser, phasei);
	    
	    sum += v.Length();
	  }

	crit = sum/double(numlines);
      }


    int calculated = 0;

    cout << endl;




    for(int i=0; i<potential_startpoints.Size(); i++)
      {
	cout << "\rFieldline Calculation " << int(100.*i/potential_startpoints.Size()) << "%"; cout.flush();

	if(randomized)
	  SetCriticalValue((double(rand())/RAND_MAX)*crit);

	if(calculated >= numlines) break;

	Calc(potential_startpoints[i],points,values,drawelems,dirstart);

	bool usable = false;

	for(int j=1; j<dirstart.Size(); j++)
	  for(int k=dirstart[j-1]; k<dirstart[j]-1; k++)
	    {
	      if(!drawelems[k] || !drawelems[k+1]) continue;
	     
	      usable = true;
 
	      // vss.SetOpenGlColor  (0.5*(values[k]+values[k+1]), minval, maxval, logscale);
              
              /*
              if (vss.usetexture == 1)
                glTexCoord1f ( 0.5*(values[k]+values[k+1]) );
              else
              */
              vss.SetOpenGlColor  (0.5*(values[k]+values[k+1]) );
	      vss.DrawCylinder (points[k], points[k+1], thickness);
	    }

	if(usable) calculated++;
      }
    cout << "\rFieldline Calculation " << 100 << "%" << endl;
    
  }



  FieldLineCalc :: FieldLineCalc(const Mesh & amesh, VisualSceneSolution & avss, const VisualSceneSolution::SolData * solution, 
				 const double rel_length, const int amaxpoints, 
				 const double rel_thickness, const double rel_tolerance, const int rk_type, const int adirection) :
    mesh(amesh), vss(avss), vsol(solution), stepper(rk_type)
  {
    mesh.GetBox (pmin, pmax);
    rad = 0.5 * Dist (pmin, pmax);
    

    maxlength = (rel_length > 0) ? rel_length : 0.5;
    maxlength *= 2.*rad;

    thickness = (rel_thickness > 0) ? rel_thickness : 0.0015;
    thickness *= 2.*rad;
    
    double auxtolerance = (rel_tolerance > 0) ? rel_tolerance : 1.5e-3;
    auxtolerance *= 2.*rad;

    stepper.SetTolerance(auxtolerance);
    
    direction = adirection;
    
    
    maxpoints = amaxpoints;

    if(direction == 0)
      {
	maxlength *= 0.5;
	maxpoints /= 2;
      }
    

    phaser = 1;
    phasei = 0;
    
    critical_value = -1;

    randomized = false;
    
  }
  


  
  void FieldLineCalc :: Calc(const Point3d & startpoint, Array<Point3d> & points, Array<double> & vals, Array<bool> & drawelems, Array<int> & dirstart)
  {
    double lami[3], startlami[3];
    double values[6];
    double dummyt(0);
    Vec3d v;
    Vec3d startv;
    Point3d newp;
    double h;
    
    double startval;
    bool startdraw;
    bool drawelem = false;
    int elnr;

    for (int i=0; i<6; i++) values[i]=0.0;
    for (int i=0; i<3; i++) lami[i]=0.0;
    for (int i=0; i<3; i++) startlami[i]=0.0;
    
    points.SetSize(0);
    vals.SetSize(0);
    drawelems.SetSize(0);

    dirstart.SetSize(0);
    dirstart.Append(0);


    int startelnr = mesh.GetElementOfPoint(startpoint,startlami,true) - 1;
    (*testout) << "p = " << startpoint << "; elnr = " << startelnr << endl;
    if (startelnr == -1)
      return;
      
    mesh.SetPointSearchStartElement(startelnr);

    if (mesh.GetDimension()==3)
      startdraw = vss.GetValues ( vsol, startelnr, startlami[0], startlami[1], startlami[2], values);
    else
      startdraw = vss.GetSurfValues ( vsol, startelnr, startlami[0], startlami[1], values);

    VisualSceneSolution::RealVec3d ( values, startv, vsol->iscomplex, phaser, phasei);

    startval = startv.Length();

    if(critical_value > 0 && fabs(startval) < critical_value)
      return;

    //cout << "p = " << startpoint << "; elnr = " << startelnr << endl;


      
    for(int dir = 1; dir >= -1; dir -= 2)
      {
	if(dir*direction < 0) continue;
	  
	points.Append(startpoint);
	vals.Append(startval);
	drawelems.Append(startdraw);
	  
	h = 0.001*rad/startval; // otherwise no nice lines; should be made accessible from outside
	
	v = startv;
	if(dir == -1) v *= -1.;

	elnr = startelnr;
	lami[0] = startlami[0]; lami[1] = startlami[1]; lami[2] = startlami[2]; 
	  

	for(double length = 0; length < maxlength; length += h*vals.Last())
	  {
	    if(v.Length() < 1e-12*rad)
	      {
		(*testout) << "Current fieldlinecalculation came to a stillstand at " << points.Last() << endl;
		break;
	      }

	    stepper.StartNextValCalc(points.Last(),dummyt,h,true);
	    stepper.FeedNextF(v);

	    while(!stepper.GetNextData(newp,dummyt,h) && elnr != -1)
	      {
		elnr = mesh.GetElementOfPoint(newp,lami,true) - 1;
		if(elnr != -1)
		  {
		    mesh.SetPointSearchStartElement(elnr);
                    if (mesh.GetDimension()==3)
                        drawelem = vss.GetValues (vsol, elnr, lami[0], lami[1], lami[2], values);
                    else
                      drawelem = vss.GetSurfValues (vsol, elnr, lami[0], lami[1], values);

		    VisualSceneSolution::RealVec3d (values, v, vsol->iscomplex, phaser, phasei);
		    if(dir == -1) v *= -1.;
		    stepper.FeedNextF(v);
		  }
	      }

	    if (elnr == -1)
	      {
		//cout << "direction " <<dir << " reached the wall." << endl;
		break;
	      }

	    points.Append(newp);
	    vals.Append(v.Length());
	    drawelems.Append(drawelem);

	    if(points.Size() % 40 == 0 && points.Size() > 1)
	      (*testout) << "Points in current fieldline: " << points.Size() << ", current position: " << newp << endl;

	    if(maxpoints > 0 && points.Size() >= maxpoints)
	      {
		break;
	      }

	    //cout << "length " << length << " h " << h << " vals.Last() " << vals.Last()  << " maxlength " << maxlength << endl;
	  }
	dirstart.Append(points.Size());
      }
  }




  
  void VisualSceneSolution :: BuildFieldLinesFromBox(Array<Point3d> & startpoints)
  {
    if(fieldlines_startarea_parameter[0] > fieldlines_startarea_parameter[3] ||
       fieldlines_startarea_parameter[1] > fieldlines_startarea_parameter[4] ||
       fieldlines_startarea_parameter[2] > fieldlines_startarea_parameter[5])
      {
	Point3d pmin, pmax;
	mesh->GetBox (pmin, pmax);
	
	fieldlines_startarea_parameter[0] = pmin.X();	
	fieldlines_startarea_parameter[1] = pmin.Y();
	fieldlines_startarea_parameter[2] = pmin.Z();
	fieldlines_startarea_parameter[3] = pmax.X();	
	fieldlines_startarea_parameter[4] = pmax.Y();
	fieldlines_startarea_parameter[5] = pmax.Z();
      }
    
    for (int i = 1; i <= startpoints.Size(); i++)
      {
	Point3d p (fieldlines_startarea_parameter[0] + double (rand()) / RAND_MAX * (fieldlines_startarea_parameter[3]-fieldlines_startarea_parameter[0]),
		   fieldlines_startarea_parameter[1] + double (rand()) / RAND_MAX * (fieldlines_startarea_parameter[4]-fieldlines_startarea_parameter[1]),
		   fieldlines_startarea_parameter[2] + double (rand()) / RAND_MAX * (fieldlines_startarea_parameter[5]-fieldlines_startarea_parameter[2]));
	
	startpoints[i-1] = p;
      }
  }

  void VisualSceneSolution :: BuildFieldLinesFromLine(Array<Point3d> & startpoints)
  {
    for (int i = 1; i <= startpoints.Size(); i++)
      {
	double s = double (rand()) / RAND_MAX;

	Point3d p (fieldlines_startarea_parameter[0] + s * (fieldlines_startarea_parameter[3]-fieldlines_startarea_parameter[0]),
		   fieldlines_startarea_parameter[1] + s * (fieldlines_startarea_parameter[4]-fieldlines_startarea_parameter[1]),
		   fieldlines_startarea_parameter[2] + s * (fieldlines_startarea_parameter[5]-fieldlines_startarea_parameter[2]));
	
	startpoints[i-1] = p;
      }
  }


  void VisualSceneSolution :: BuildFieldLinesFromFile(Array<Point3d> & startpoints)
  {
    ifstream * infile;

    infile = new ifstream(fieldlines_filename.c_str());

    //cout << "reading from file " << fieldlines_filename << endl;

    int numpoints = 0;

    string keyword;

    
    double dparam;
    int iparam;

    while(infile->good())
      {
	(*infile) >> keyword;

	if(keyword == "point") numpoints++;
	else if(keyword == "line" || keyword == "box")
	  {
	    for(int i=0; i<6; i++) (*infile) >> dparam;
	    (*infile) >> iparam;
	    numpoints += iparam;
	  }
      }

    delete infile;


    //cout << numpoints << " startpoints" << endl;

    startpoints.SetSize(numpoints);
    
    infile = new ifstream(fieldlines_filename.c_str());

    numpoints = 0;

    while(infile->good())
      {
	(*infile) >> keyword;

	if (keyword == "point")
	  {
	    (*infile) >> startpoints[numpoints].X(); (*infile) >> startpoints[numpoints].Y(); (*infile) >> startpoints[numpoints].Z();
	    numpoints++;
	  }
	else if (keyword == "line" || keyword == "box")
	  {
	    for(int i=0; i<6; i++) (*infile) >> fieldlines_startarea_parameter[i];
	    (*infile) >> iparam;

	    Array<Point3d> auxpoints(iparam);
	    
	    if (keyword == "box")
	      BuildFieldLinesFromBox(auxpoints);
	    else if (keyword == "line")
	      BuildFieldLinesFromLine(auxpoints);
	    
	    for(int i=0; i<iparam; i++)
	      {
		startpoints[numpoints] = auxpoints[i];
		numpoints++;
	      }
	  }

	//cout << "startpoints " << startpoints << endl;
      }

    delete infile;
    
    

    
  }

  
  void VisualSceneSolution :: BuildFieldLinesFromFace(Array<Point3d> & startpoints)
  {
    Array<SurfaceElementIndex> elements_2d;
    
    //cout << "fieldlines_startface " << fieldlines_startface << endl;
    mesh->GetSurfaceElementsOfFace(fieldlines_startface,elements_2d);
    if(elements_2d.Size() == 0)
      {
	cerr << "No Elements on selected face (?)" << endl;
	return;
      }
    Vec3d v1,v2,cross;
    
    double area = 0;

	int i;
    for(i=0; i<elements_2d.Size(); i++)
      {
	const Element2d & elem = mesh->SurfaceElement(elements_2d[i]);
	
	v1 = mesh->Point(elem[1]) - mesh->Point(elem[0]);
	v2 = mesh->Point(elem[2]) - mesh->Point(elem[0]);
	cross = Cross(v1,v2);
	area += cross.Length();
	
	if(elem.GetNV() == 4)
	  {
	    v1 = mesh->Point(elem[2]) - mesh->Point(elem[0]);
	    v2 = mesh->Point(elem[3]) - mesh->Point(elem[0]);
	    cross = Cross(v1,v2);
	    area += cross.Length();
	  }
      }
    
    int startpointsp = 0;
    i = 0;
    
    while(startpointsp < startpoints.Size())
      {
	const Element2d & elem = mesh->SurfaceElement(elements_2d[i]);
	
	int numtri = (elem.GetNV() == 3) ? 1 : 2;
	
	for(int tri = 0; startpointsp < startpoints.Size() && tri<numtri; tri++)
	  {
	    
	    if(tri == 0)
	      {
		v1 = mesh->Point(elem[1]) - mesh->Point(elem[0]);
		v2 = mesh->Point(elem[2]) - mesh->Point(elem[0]);
		cross = Cross(v1,v2);
	      }
	    else if(tri == 1)
	      {
		v1 = mesh->Point(elem[2]) - mesh->Point(elem[0]);
		v2 = mesh->Point(elem[3]) - mesh->Point(elem[0]);
		cross = Cross(v1,v2);
	      }
	    
	    double thisarea = cross.Length();
	    
	    int numloc = int(startpoints.Size()*thisarea/area);
	    if(double (rand()) / RAND_MAX < startpoints.Size()*thisarea/area - numloc)
	      numloc++;
	    
	    for(int j=0; startpointsp < startpoints.Size() && j<numloc; j++)
	      {
		double s = double (rand()) / RAND_MAX;
		double t = double (rand()) / RAND_MAX;
		if(s+t > 1)
		  {
		    s = 1.-s; t = 1.-t;
		  }
		startpoints[startpointsp] = mesh->Point(elem[0]) + s*v1 +t*v2;
		startpointsp++;
	      }
	  }
	i++;
	if(i == elements_2d.Size()) i = 0;
      } 
    
  }


  void VisualSceneSolution :: BuildFieldLinesPlot ()
  {
    if (fieldlinestimestamp >= solutiontimestamp) 
      return;
    fieldlinestimestamp = solutiontimestamp;
    

    if (fieldlineslist)
      glDeleteLists (fieldlineslist, num_fieldlineslists);

    if (vecfunction == -1)
      return;

    const SolData * vsol = soldata[fieldlines_vecfunction];

    num_fieldlineslists = (vsol -> iscomplex && !fieldlines_fixedphase) ? 100 : 1;
   

    FieldLineCalc linecalc(*mesh,*this,vsol,
			   fieldlines_rellength,fieldlines_maxpoints,fieldlines_relthickness,fieldlines_reltolerance,fieldlines_rktype);

    if(fieldlines_randomstart) 
      linecalc.Randomized();

    fieldlineslist = glGenLists (num_fieldlineslists);

    int num_startpoints = num_fieldlines / num_fieldlineslists;
    if (num_fieldlines % num_fieldlineslists != 0) num_startpoints++;

    if(fieldlines_randomstart)
      num_startpoints *= 10;

    
    Array<Point3d> startpoints(num_startpoints);
    

    for (int ln = 0; ln < num_fieldlineslists; ln++)
      {
	if(fieldlines_startarea == 0)
	  BuildFieldLinesFromBox(startpoints);
	else if(fieldlines_startarea == 1)
	  BuildFieldLinesFromFile(startpoints);
	else if(fieldlines_startarea == 2)
	  BuildFieldLinesFromFace(startpoints);


	    
	double phi;
	
	if(vsol -> iscomplex)
	  {
	    if(fieldlines_fixedphase)
	      phi = fieldlines_phase;
	    else
	      phi = 2*M_PI*ln / num_fieldlineslists;
	  }
	else
	  phi = 0;

	cout << "phi = " << phi << endl;

	double phaser = cos(phi), phasei = sin(phi);
	

        glNewList(fieldlineslist+ln, GL_COMPILE);
        
        SetTextureMode (usetexture);
	linecalc.GenerateFieldLines(startpoints,num_fieldlines / num_fieldlineslists+1,
				    fieldlineslist+ln,minval,maxval,logscale,phaser,phasei);

        glEndList ();

      }
  }



  
}


#endif // NOTCL
