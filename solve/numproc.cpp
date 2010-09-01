#include <solve.hpp>

#include <tcl.h>
#if TCL_MAJOR_VERSION==8 && TCL_MINOR_VERSION>=4
#define tcl_const const
#else
#define tcl_const
#endif


#include <parallelngs.hpp>


namespace ngsolve
{
  using namespace ngsolve;
  using namespace ngparallel;



  NumProc :: NumProc (PDE & apde, const int acallposition)
    : NGS_Object (apde.GetMeshAccess(), "numproc"), pde(apde)
  {
    callposition = acallposition;
  }
  
  NumProc :: ~NumProc()
  {
    ;
  }

  /*
  void NumProc :: Do()
  {
    ;
  }

  void NumProc :: Do (LocalHeap & lh)
  {
    Do();
  }
  */

  void NumProc :: PrintReport (ostream & ost)
  {
    ost << "Base-class NumProc" << endl;
  }

  void NumProc :: PrintDoc (ostream & ost)
  {
    ost << "No documentation available" << endl;
  }








  /* ***************************** Numproc CalcFlux ************************** */



  ///
  class NumProcCalcFlux : public NumProc
  {
  protected:
    ///
    BilinearForm * bfa;
    ///
    GridFunction * gfu;
    ///
    GridFunction * gfflux;
    /// compute flux, not gradient
    bool applyd;
    ///
    // bool useall;
    ///
    int domain;
  public:
    ///
    NumProcCalcFlux (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcCalcFlux();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcCalcFlux (pde, flags);
    }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual string GetClassName () const
    {
      return "Calc Flux";
    }


    virtual void PrintReport (ostream & ost)
    {
      ost << GetClassName() << endl
	  << "Bilinear-form    = " << bfa->GetName() << endl
	  << "Differential-Op  = " << bfa->GetIntegrator(0)->Name() << endl
	  << "Gridfunction-In  = " << gfu->GetName() << endl
	  << "Gridfunction-Out = " << gfflux->GetName() << endl
	  << "apply coeffs     = " << applyd << endl;
      // << "use all int'rs   = " << useall << endl;
    }
  };


  NumProcCalcFlux :: NumProcCalcFlux (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL));
    if (bfa->NumIntegrators()==0)
      throw Exception ("bilinearform used for CalcFlux needs at least one integrator");

    gfu = pde.GetGridFunction (flags.GetStringFlag ("solution", NULL));
    gfflux = pde.GetGridFunction (flags.GetStringFlag ("flux", NULL));
    applyd = flags.GetDefineFlag ("applyd");
    // useall = flags.GetDefineFlag ("useall");
    domain = static_cast<int>(flags.GetNumFlag("domain",0))-1;
  }


  NumProcCalcFlux :: ~NumProcCalcFlux()
  {
    ;
  }

  void NumProcCalcFlux :: PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc CalcFlux:\n" \
      "-----------------\n" \
      "Computes the natural flux of the bvp:\n\n"\
      "- Heat flux for thermic problems\n"\
      "- Stresses for mechanical problems\n"\
      "- Induction for magnetostatic problems\n\n"\
      "Required flags:\n" 
      "-bilinearform=<bfname>\n" 
      "    the first integrator for the bf computes the flux\n" \
      "-solution=<gfname>\n" \
      "    grid-function providing the primal solution field\n" \
      "-flux=<gfname>\n" \
      "    grid-function used for storing the flux (e.g., vector-valued L2)\n\n" \
      "\nOptional flags:\n" \
      "-applyd   apply coefficient matrix (compute either strains or stresses, B-field or H-field,..\n"
	<< endl;

    //      "-useall   use all integrators for computing the flux, and add up result\n" 
  }


  void NumProcCalcFlux :: Do(LocalHeap & lh)
  {
    CalcFluxProject (pde.GetMeshAccess(), *gfu, *gfflux,
		     *bfa->GetIntegrator(0),
		     applyd, domain, lh);
    
    /*
    // useall - version currently not supported
    gfflux->GetVector() = 0;
    for (int k = 0; k < bfa->NumIntegrators(); k++)
    {
    if (!bfa->GetFESpace().IsComplex())
    CalcFlux (pde.GetMeshAccess(),
    dynamic_cast<const S_GridFunction<double>&> (*gfu), 
    dynamic_cast<S_GridFunction<double>&> (*gfflux), 
    *bfa->GetIntegrator(k),
    applyd, 1, domain);
    else
    CalcFlux (pde.GetMeshAccess(),
    dynamic_cast<const S_GridFunction<Complex>&> (*gfu), 
    dynamic_cast<S_GridFunction<Complex>&> (*gfflux), 
    *bfa->GetIntegrator(k),
    applyd, 1, domain);
    }
    */
  }





  /* ***************************** Numproc CalcFlux ************************** */



  ///
  class NumProcSetValues : public NumProc
  {
  protected:
    ///
    GridFunction * gfu;
    ///
    CoefficientFunction * coef;
    ///
    bool boundary;
    ///
    int component;
  public:
    ///
    NumProcSetValues (PDE & apde, const Flags & flags)
    : NumProc (apde)
    {
      gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
      coef = pde.GetCoefficientFunction (flags.GetStringFlag ("coefficient", ""));
      boundary = flags.GetDefineFlag ("boundary");
      component = int (flags.GetNumFlag ("component", 0))-1;
    }

    ///
    virtual ~NumProcSetValues() { ; }

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcSetValues (pde, flags);
    }

    static void PrintDoc (ostream & ost)
    {
      ost << 
	"\n\nNumproc setvalues:\n"		\
	"-----------------\n"				\
	"Set a gridfunction to given values\n\n" \
	"Required flags:\n"
 	"" \
	"-gridfunction=<gfname>\n"						\
	"    grid-function to be set\n"	\
	"-coefficient=<coefname>\n"						\
	"    coefficient providing values\n\n" \
	"\nOptional flags:\n"						\
	"-boundary\n only boundary values are set\n" \
	"-component=<comp>\n set only this component (of CompoundFESpace)\n" \
	  << endl;
    }
    
    ///
    virtual void Do(LocalHeap & lh)
    {
      GridFunction * hgfu = gfu;
      if (component != -1)
	hgfu = gfu->GetComponent(component);

      SetValues (pde.GetMeshAccess(), *coef, 
		 *hgfu, boundary, 0, lh);

      if (component != -1)
	delete hgfu;
    }

    ///
    virtual string GetClassName () const
    {
      return "SetValues";
    }

    virtual void PrintReport (ostream & ost)
    {
      ost << GetClassName() << endl
	  << "Gridfunction-Out = " << gfu->GetName() << endl;
    }
  };



  /* **************************** Numproc Draw Flux ********************************* */



  ///
  class NumProcDrawFlux : public NumProc
  {
  protected:
    netgen::SolutionData * vis;
    ///
    BilinearForm * bfa;
    ///
    GridFunction * gfu;
    /// compute flux, not gradient
    bool applyd;
    ///
    bool useall;

    //BilinearFormIntegrator * bfi2d;
    //BilinearFormIntegrator * bfi3d;
    string label;

  public:
    ///
    NumProcDrawFlux (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcDrawFlux();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcDrawFlux (pde, flags);
    }
    ///
    virtual void Do (LocalHeap & lh);

    static void PrintDoc (ostream & ost);
    ///
    virtual string GetClassName () const
    {
      return "Draw Flux";
    }


    virtual void PrintReport (ostream & ost)
    {
      ost << GetClassName() << endl;
      if (bfa) ost << "Bilinear-form    = " << bfa->GetName() << endl;
      if (bfa) ost << "Differential-Op  = " << bfa->GetIntegrator(0)->Name() << endl;
      if (gfu) ost << "Gridfunction-In  = " << gfu->GetName() << endl;
      ost << "apply coeffs     = " << applyd << endl;
    }
  };


  NumProcDrawFlux :: NumProcDrawFlux (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = NULL;
    gfu = NULL;
    
    if (flags.GetDefineFlag ("order"))
      {
	Ng_SolutionData soldata;
	Ng_InitSolutionData (&soldata);
	soldata.name = "order";
	soldata.soltype = NG_SOLUTION_ELEMENT_ORDER;
	Ng_SetSolutionData (&soldata);
	return;
      }
    if (flags.GetDefineFlag ("marked"))
      {
	Ng_SolutionData soldata;
	Ng_InitSolutionData (&soldata);
	soldata.name = "marked";
	soldata.soltype = NG_SOLUTION_MARKED_ELEMENTS;
	Ng_SetSolutionData (&soldata);
	return;
      }

    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));

    gfu = pde.GetGridFunction (flags.GetStringFlag ("solution", ""));
    if(!gfu)
      gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""));

    applyd = flags.GetDefineFlag ("applyd");
    label = flags.GetStringFlag ("label", "");
    useall = flags.GetDefineFlag ("useall");

    Array<BilinearFormIntegrator *> bfi2d, bfi3d;

    for (int i = 0; i < bfa->NumIntegrators(); i++)
      {
	if ((!bfi3d.Size() || useall) && bfa->GetIntegrator(i)->DimElement() == 3)
	  bfi3d.Append(bfa->GetIntegrator(i));
	if ((!bfi2d.Size() || useall) && bfa->GetIntegrator(i)->DimElement() == 2)
	  bfi2d.Append(bfa->GetIntegrator(i));
      }

    /*
      if (bfi2d) cout << "bfi2d = " 
      << bfi2d->Name()
      << ", dim = " << bfi2d->DimFlux()
      << endl;
      if (bfi3d) cout << "bfi3d = " << bfi3d->Name()
      << ", dim = " << bfi3d->DimFlux()
      << endl;
    */

    if (!bfa->GetFESpace().IsComplex())
      vis = new VisualizeGridFunction<double> (ma, gfu, bfi2d, bfi3d, applyd);
    else
      vis = new VisualizeGridFunction<Complex> (ma, gfu, bfi2d, bfi3d, applyd);

    Ng_SolutionData soldata;
    Ng_InitSolutionData (&soldata);
  
    // soldata.name = const_cast<char*> (gfu->GetName().c_str());
    soldata.name = (char*)label.c_str();
    soldata.data = 0;
    soldata.components = vis->GetComponents();
    soldata.iscomplex = vis->IsComplex();
    soldata.draw_surface = bfi2d.Size() != 0;
    soldata.draw_volume  = bfi3d.Size() != 0;
    soldata.dist = 1;
    soldata.soltype = NG_SOLUTION_VIRTUAL_FUNCTION;
    soldata.solclass = vis;
    Ng_SetSolutionData (&soldata);
  }

  NumProcDrawFlux :: ~NumProcDrawFlux()
  {
    ;
  }

  void NumProcDrawFlux :: Do(LocalHeap & lh)
  {
    cout << "Num-proc draw flux" << endl;
  }










  void NumProcDrawFlux :: PrintDoc (ostream & ost)

  {
    ost <<
      "\n\nNumproc DrawFlux:\n" \
      "-----------------\n" \
      "Adds the natural flux to the visualization dialogbox:\n"\
      "It takes the first integrator of the bilinear-form\n" \
      "- Heat flux for thermic problems\n"\
      "- Stresses for mechanical problems\n"\
      "- Induction for magnetostatic problems\n\n"\
      "Required flags:\n"
      "-bilinearform=<bfname>\n" \
      "    the first integrator for the bf computes the flux\n" \
      "-solution=<gfname>\n" \
      "    grid-function providing the primal solution field\n" \
      "\nOptional flags:\n" \
      "-applyd\n" \
      "    apply coefficient matrix (compute either strains or stresses, B-field or H-field,..\n"\
      "-label=<name>\n" \
      "    label printed in the visualization dialogbox\n\n";

  }



  /* **************************** Numproc Evaluate ********************************* */




  ///
  class NumProcEvaluate : public NumProc
  {
  protected:
    ///
    BilinearForm * bfa;
    ///
    LinearForm * lff;
    ///
    GridFunction * gfu;
    ///
    GridFunction * gfv;
    ///
    Vector<double> point;
    ///
    Array<int> domains;
    ///
    Vector<double> point2;
    ///
    Vector<double> point3;
    ///
    Vector<double> point4;
    ///
    bool integrateonplanes;
    ///
    bool usepoint3and4;
    ///
    int variabledirection;
    ///
    int n[3];
    ///
    string filename, text;
    ///
    string variablename;
    ///
    bool applyd;
    ///
    bool hermitsch;
    ///
    int component;
    ///
    int outputprecision;
  public:
    ///
    NumProcEvaluate (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcEvaluate();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcEvaluate (pde, flags);
    }
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "Evaluate";
    }

  };





  NumProcEvaluate :: NumProcEvaluate (PDE & apde, const Flags & flags)
    : NumProc (apde), point(1), point2(1), point3(1), point4(1)
  {
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""), 1); 
    lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", ""), 1);
    gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""), 0); 
    gfv = pde.GetGridFunction (flags.GetStringFlag ("gridfunction2", ""), 1); 

    variablename = flags.GetStringFlag ("resultvariable", "");

    if (flags.NumListFlagDefined ("point"))
      {
	const Array<double> & p = flags.GetNumListFlag ("point");
	point.SetSize(p.Size());
	for (int i = 0; i < p.Size(); i++)
	  point(i) = p[i];
	//      cout << "point = " << point << endl;
      }

    if (flags.NumListFlagDefined ("domains"))
      {
	const Array<double> & ds = flags.GetNumListFlag("domains");
	domains.SetSize(ds.Size());
	for (int i = 0; i < ds.Size(); i++)
	  domains[i] = int(ds[i])-1;
      }

    if (flags.NumListFlagDefined ("point2"))
      {
	const Array<double> & p = flags.GetNumListFlag ("point2");
	point2.SetSize(p.Size());
	for (int i = 0; i < p.Size(); i++)
	  point2(i) = p[i];
	//      cout << "point2 = " << point2 << endl;
      }

    usepoint3and4 = (flags.NumListFlagDefined ("point3") &&
		     flags.NumListFlagDefined ("point4"));
    if (usepoint3and4)
      {
	const Array<double> & p3 = flags.GetNumListFlag ("point3");
	point3.SetSize(p3.Size());
	for (int i = 0; i < p3.Size(); i++)
	  point3(i) = p3[i];
	const Array<double> & p4 = flags.GetNumListFlag ("point4");
	point4.SetSize(p4.Size());
	for (int i = 0; i < p4.Size(); i++)
	  point4(i) = p4[i];
      }
    
    integrateonplanes = flags.GetDefineFlag("integrateonplanes");

    variabledirection = static_cast<int>(flags.GetNumFlag("variabledirection",0))-1;
    
    n[0] = static_cast<int>(flags.GetNumFlag("n1",0));
    n[1] = static_cast<int>(flags.GetNumFlag("n2",0));
    n[2] = static_cast<int>(flags.GetNumFlag("n3",0));

    text = flags.GetStringFlag ("text","value");

    if(flags.StringFlagDefined("filename"))
      filename = pde.GetDirectory() + dirslash + flags.GetStringFlag("filename","");
    else
      filename = "err.out";

    

    applyd = flags.GetDefineFlag ("applyd");
    hermitsch = flags.GetDefineFlag ("hermitsch");

    outputprecision = (pde.ConstantUsed("outputprecision")) ? int(pde.GetConstant("outputprecision")) : -1;
    if(flags.NumFlagDefined("outputprecision"))
      outputprecision = int(flags.GetNumFlag("outputprecision",-1));

    component = static_cast<int>(flags.GetNumFlag("cachecomp",1))-1;
  }

  NumProcEvaluate :: ~NumProcEvaluate()
  {
    ;
  }



  void NumProcEvaluate :: PrintDoc (ostream & ost)

  {
    ost <<
      "\n\nNumproc Evaluate:\n" \
      "-----------------\n" \
      "Evaluates linearform or bilinearform, or pointvalues:\n"\
      "Required flags:\n" \
      "-gridfunction=<gfname>\n" \
      "    gridfunction to evaluate\n" \
      "\nOptional flags:\n" \
      "-bilinearform=<bfname>\n" \
      "    evaluates bilinear-form bfname(gfname,gfname2)\n"\
      "    needs second gridfunction <gfname2>\n" \
      "-gridfunction2=<gfname2>\n" \
      "    second gridfunciton for bilinear-form evaluation\n" \
      "-cachecomp=<n>"\
      "    for gridfunctions with cachesize > 1 use the given component"
      "-linearform=<lfname>\n" \
      "    evaluates linearform lfname(gfname)\n" \
      "-point=[x,y,z]\n" \
      "    evaluates diffop applied to gridfunction in point p\n" \
      "    diffop taken from first term in bilinear-form\n" \
      "-domains=[d1,...,dn]\n" \
      "    the point has to lie in one of the specified domains\n" \
      "-point2=[x2,y2,z2]\n" \
      "    evaluates diffop applied to gridfunction in line p1-p2, and writes to file\n" \
      "-applyd\n" \
      "    evaluates flux instead of derivatives\n" \
      "-integrateonplanes\n" \
      "    diffop applied to gridfunction is integrated on a sequence of rectangles bounded by a brick\n" \
      "       -point=... and -point2=... now set the minimal and maximal point of the brick\n" \
      "       -point3=..., -point4=...\n"\
      "            if these additional points are set, the bounds are not given by a brick but by\n" \
      "            a parallelepiped\n"\
      "       -variabledirection=1|2|3\n" \
      "            sets the direction of the normal vector of the rectangles (x|y|z)\n"\
      "       -n1=<number of points>\n" \
      "       -n2=<number of points>\n" \
      "       -n3=<number of points>\n" \
      "            set the number of integration points resp. the number of rectangles in the 3 directions\n" \
      "-text=blabla \n" \
      "    prints text \n" \
      " -resultvariabe=<varname> \n" \
      "    stores scalar, real results in variable <varname>\n" \
      "-filename=<fn> \n" \
      "    writes values in file \n" \
      "-outputprecision=<n> \n" \
      "    sets the number of digits for the output of the results\n";
  }




  void NumProcEvaluate :: Do(LocalHeap & lh)
  {
    int i, j, k;
    double result(0);
    ofstream ofile (filename.c_str());

    int old_cout_precision = cout.precision();
    int old_ofile_precision = ofile.precision();

    if(outputprecision > 0)
      {
	cout.precision(outputprecision);
	ofile.precision(outputprecision);
      }

    if (lff)
      {
	if (!gfu)
	  throw Exception ("evaluate linear-form needs an argument -gridfunction=u");

	cout << "<" << lff->GetName() << ", " << gfu->GetName() << "> = " << flush;
	if (!lff->GetFESpace().IsComplex())
	  {
	    result = S_InnerProduct<double>(lff->GetVector(), gfu->GetVector());
	    cout << result << endl;
	  }
	else
	  cout << S_InnerProduct<Complex>(lff->GetVector(), gfu->GetVector()) << endl;
      }
    else if (point.Size() >= 2)
      {
	const BilinearFormIntegrator & bfi = (bfa) ? *bfa->GetIntegrator(0) : *gfu->GetFESpace().GetEvaluator();

	if (point2.Size() >= 2)
	  {
	    pde.GetEvaluateFiles() += " ";
	    pde.GetEvaluateFiles() += filename;
	    pde.GetEvaluateFiles() += " ";
	    pde.GetEvaluateFiles() += text;
	    
	    if(!integrateonplanes)
	      {
		// plot values along line p1-p2

		int numpoints = 0;
		for(i = 0; i<3; i++) 
		  if(n[i] > numpoints) numpoints = n[i];
		if(numpoints == 0)
		  numpoints = 1000;

		
		for (i = 0; i <= numpoints; i++)
		  {
		    lh.CleanUp();
		    FlatVector<double> p(point.Size(), lh);
		    p = point + (double(i)/numpoints) * (point2-point);
		    
		    if (!gfu->GetFESpace().IsComplex())
		      {
			FlatVector<double> pflux(bfi.DimFlux(), lh);
			bool ok =
			  CalcPointFlux<double> (ma, *gfu, p,
						 pflux, bfi, applyd, lh, component);
			
			if (i==0)
			  {
			    ofile << "#nglineplotinfo ";
			    for (j = 0; j < p.Size(); j++) ofile << "coord" << j+1 << " ";
			    for (j = 0; j < pflux.Size(); j++) ofile << "flux" << j+1 << " ";
			    ofile << endl;
			  }

			for (j = 0; j < p.Size(); j++)
			  ofile << p(j) << " ";
			for (j = 0; ok && j < pflux.Size(); j++)
			  ofile << pflux(j) << " ";
			ofile << endl;
		      }
		    else
		      {
			FlatVector<Complex> pflux(bfi.DimFlux(), lh);
			bool ok =
			  CalcPointFlux (ma, *gfu, p,
					 pflux, bfi, applyd, lh, component);
			
			if (i==0)
			  {
			    ofile << "#nglineplotinfo ";
			    for (j = 0; j < p.Size(); j++) ofile << "coord" << j+1 << " ";
			    for (j = 0; j < pflux.Size(); j++) ofile << "real(flux" << j+1 << ") imag(flux" << j+1 << ") ";
			    ofile << endl;
			  }

			
			for (j = 0; j < p.Size(); j++)
			  ofile << p(j) << " ";
			for (j = 0; ok && j < pflux.Size(); j++)
			  ofile << pflux(j).real() << " " << pflux(j).imag() << " ";
			ofile << endl;
		      }
		  }
	      }
	    else
	      {
		// plot values in parallelepiped (p,p2,p3,p4) in direction [variabledirection] where the values on each plane are integrated

		if(!usepoint3and4)
		  {
		    point3.SetSize(point.Size());
		    point4.SetSize(point.Size());
		    
		    point3(0) = point(0); point3(1) = point2(1); 
		    point4(0) = point(0); point4(1) = point(1); 
		    point2(1) = point(1);
		    if(point.Size() > 2)
		      {
			point3(2) = point(2);	
			point4(2) = point2(2);
			point2(2) = point(2);
		      }
		  }

		Vector<double> dirvecv(point.Size()),dirvec1(point.Size()),dirvec2(point.Size());

		int dir1,dir2;

		if(variabledirection == 0) 
		  { 
		    dir1 = 1; dir2 = 2;
		    dirvecv = point2 - point;
		    dirvec1 = point3 - point;
		    dirvec2 = point4 - point;
		  }
		else if(variabledirection == 1) 
		  { 
		    dir1 = 0; dir2 = 2; 
		    dirvecv = point3 - point;
		    dirvec1 = point2 - point;
		    dirvec2 = point4 - point;
		  }
		else if(variabledirection == 2) 
		  {
		    dir1 = 0; dir2 = 1; 
		    dirvecv = point4 - point;
		    dirvec1 = point2 - point;
		    dirvec2 = point3 - point;
		  }
		else
		  throw Exception ("numproc evaluate: flag \"variabledirection\" not correctly set");
		
		if(n[dir1] < 1) n[dir1] = 1;
		if(n[dir2] < 1) n[dir2] = 1;
		if(n[variabledirection] < 1) n[variabledirection] = 1000;

		double pos1,pos2;

		Array < Array < double > * > values(n[variabledirection]+1);
		values = NULL;

		FlatVector<double> p(point.Size(), lh);
		
		void * heapp = lh.GetPointer();

		bool complexflux (false);

		bool firstone(true);
		
		for(i=0, pos1 = 0.5/n[dir1]; i<n[dir1]; i++, pos1+=1./n[dir1])
		  for(j=0, pos2 = 0.5/n[dir2]; j<n[dir2]; j++, pos2+=1./n[dir2])
		    for(k=0; k <= n[variabledirection]; k++)
		      {
			p = point + pos1*dirvec1 + pos2*dirvec2 + (double(k)/n[variabledirection]) * dirvecv;
			
			if(firstone)
			  {
			    IntegrationPoint dummyip;
			    ma.FindElementOfPoint(p,dummyip,true);
			    firstone = false;
			  }
			
			if (!gfu->GetFESpace().IsComplex())
			  {
			    FlatVector<double> pflux(bfi.DimFlux(), lh);
                            CalcPointFlux (ma, *gfu, p,
                                           pflux, bfi, applyd, lh, component);
			    
			    if(values[k] == NULL)
			      {
				values[k] = new Array<double>(pflux.Size());
				(*values[k]) = 0;
			      }
			    
			    for(int ii=0; ii<pflux.Size(); ii++) (*values[k])[ii] += pflux[ii];
			  }
			else
			  {
			    FlatVector<Complex> pflux(bfi.DimFlux(), lh);
			    complexflux = true;
                            CalcPointFlux (ma, *gfu, p,
                                           pflux, bfi, applyd, lh, component);
			    
			    if(values[k] == NULL)
			      {
				values[k] = new Array<double>(2*pflux.Size());
				(*values[k]) = 0;
			      }
			    
			    for(int ii=0; ii<pflux.Size(); ii++)
			      {
				(*values[k])[2*ii] += pflux[ii].real();
				(*values[k])[2*ii+1] += pflux[ii].imag();
			      }
			  }
			lh.CleanUp(heapp);
		      }
	      
		lh.CleanUp();

		Vector<double> cross(3); cross = 0;
		if(point.Size() > 2)
		  {
		    cross(0) = dirvec1(1)*dirvec2(2) - dirvec1(2)*dirvec2(1);
		    cross(1) = dirvec1(2)*dirvec2(0) - dirvec1(0)*dirvec2(2);
		  }
		cross(2) = dirvec1(0)*dirvec2(1) - dirvec1(1)*dirvec2(0);

		double factor = L2Norm(cross)/(n[dir1]*n[dir2]);

		ofile << "#nglineplotinfo coord" << variabledirection+1 << " ";
		
		if(!complexflux)
		  for(j = 0; j<(*values[0]).Size(); j++) ofile << "flux" << j+1 << " ";
		else
		  for(j = 0; j<(*values[0]).Size()/2; j++) ofile << "real(flux" << j+1 << ") imag(flux" << j+1 << ") ";

		ofile << endl;

		double lv = L2Norm(dirvecv);
		double pv = InnerProduct(dirvecv,point)/lv;
		for(i=0; i<values.Size(); i++)
		  {
		    ofile << pv + (double(i)/n[variabledirection]) * lv << " ";
		    for(j = 0; j<(*values[i]).Size(); j++)
		      ofile << factor*(*values[i])[j] << " ";
		    ofile << endl;

		    delete values[i];
		  }
	      }
	  }
	else
	  {
	    cout << text << " (" << point(0);
	    for (int i = 1; i < point.Size(); i++)
	      cout << "," << point(i);
	    cout << ") = " << flush;

	    if (!gfu->GetFESpace().IsComplex())
	      {
		FlatVector<double> pflux(bfi.DimFlux(), lh);
		CalcPointFlux (ma, *gfu, point, domains,
			       pflux, bfi, applyd, lh, component);
		for (int i = 0; i < pflux.Size(); i++)
		  cout << pflux(i) << "   ";
		cout << endl;
	      }
	    else
	      {
		FlatVector<Complex> pflux(bfi.DimFlux(), lh);
		CalcPointFlux (ma, *gfu, point, domains,
			       pflux, bfi, applyd, lh, component);
		for (int i = 0; i < pflux.Size(); i++)
		  cout << pflux(i) << "   ";
		cout << endl;
	      }
	  }
      }
    else if (bfa)
      {

	cout << " bfa is ok " << endl; 
	if (!gfu)
	  throw Exception ("evaluate bilinear-form needs an argument -gridfunction=u");
	if (!gfv)
	  throw Exception ("evaluate bilinear-form needs an argument -gridfunction2=v");

	cout << bfa->GetName() << "(" << gfu->GetName() << ", " << gfv->GetName() << ") = " << flush;
	BaseVector & vecu = gfu->GetVector();
	BaseVector & vecv = gfv->GetVector();
	BaseVector & hv = *vecu.CreateVector();
	bfa->GetMatrix().Mult (vecu, hv);
	if (!bfa->GetFESpace().IsComplex())
	  {
	    result = S_InnerProduct<double>(vecv, hv);
	    //cout << setprecision(16) << result << endl;

	    cout << "bf(gf,gf2) = " << result << endl;
	    ofile << "err = " << result << endl;
	    //cout << "bf(gf,gf2) = " << setprecision(16) << result << endl;
	    //ofile << "err = " << setprecision(16) << result << endl;

	  }
	else
	  {
	    if (!hermitsch)
	      cout << S_InnerProduct<Complex>(vecv, hv) << endl;
	    else
	      {
		FlatVector<Complex> fvecv = vecv.FVComplex();
		FlatVector<Complex> fhv = hv.FVComplex();
		Complex sum = 0.0;
		for (int i = 0; i < fvecv.Size(); i++)
		  sum += fvecv(i) * conj (fhv(i));
		cout << sum << endl;
	      }
	  }
	delete &hv;
      }
    
    pde.GetVariable(variablename,true) = result;
    
    cout.precision(old_cout_precision);
    ofile.precision(old_ofile_precision);
    ofile.close();
  }

  void NumProcEvaluate :: PrintReport (ostream & ost)
  {
    ost << "NumProcEvaluate:" << endl;
  }






  /////////////////////////////////////////////////////////////////////////////
  ///
  class NumProcAnalyze : public NumProc
  {
  protected:
    GridFunction * gfu;
    ///
    string variablename;
    ///
    bool nodistinction;
    ///
    bool volanalyze;
    ///
    bool surfanalyze;
    ///
    int component;
    ///
    Array<int> surfdomains;
    ///
    Array<int> voldomains;
  public:
    ///
    NumProcAnalyze (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcAnalyze();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcAnalyze (pde, flags);
    }
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "Analyze";
    }

  };





  NumProcAnalyze :: NumProcAnalyze (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""), 0); 

    variablename = flags.GetStringFlag ("resultvariable", "");


    volanalyze = flags.GetDefineFlag("volume");
    surfanalyze = flags.GetDefineFlag("surface");
    nodistinction = flags.GetDefineFlag("nodistinction");

    if(!volanalyze && !surfanalyze) volanalyze = true;

    component = static_cast<int>(flags.GetNumFlag("comp",0)) - 1;

    if(flags.NumListFlagDefined("voldomains"))
      {
	voldomains.SetSize(flags.GetNumListFlag("voldomains").Size());
	for(int i=0; i<voldomains.Size(); i++)
	  {
	    voldomains[i] = static_cast<int>(flags.GetNumListFlag("voldomains")[i]);
	  }
      }
    
    if(flags.NumListFlagDefined("surfdomains"))
      {
	surfdomains.SetSize(flags.GetNumListFlag("surfdomains").Size());
	for(int i=0; i<surfdomains.Size(); i++)
	  {
	    surfdomains[i] = static_cast<int>(flags.GetNumListFlag("surfdomains")[i]);
	  }
      }
	    
  }

  NumProcAnalyze :: ~NumProcAnalyze()
  {
    ;
  }



  void NumProcAnalyze :: PrintDoc (ostream & ost)

  {
    ost <<
      "\n\nNumproc Analyze:\n" \
      "--------------------------\n" \
      "\nAnalyzes a gridfunction, i.e. calculates maximum, minimum, and average value\n"\
      "Required flags:\n" \
      "-gridfunction=<gfname>\n" \
      "    gridfunction to analyze\n" \
      "\nOptional flags:\n" \
      "-volume\n"\
      "    analyze inside volumes\n" \
      "-surface\n" \
      "    analyze on surfaces\n" \
      "-comp=<component>\n" \
      "    analyze only one given component\n" \
      "-voldomains=<domainslist>\n" \
      "    analyze only inside given domains\n" \
      "-surfdomains=<domainslist>\n" \
      "    analyze only on given surfaces\n" \
      "-resultvariable=<basename>\n" \
      "    the results are saved to the variables basename.<vol|surf><number>.comp<number>.<min|max|av>\n" \
      "-nodistinction\n" \
      "    if set then there is no distinction of the different components and domains. Results are saved to basename.<min|max|av>\n";
  }




  void NumProcAnalyze :: Do(LocalHeap & lh)
  {
    int dom;

    bool writevar = (strcmp(variablename.c_str(), "") != 0);
    string actvarname;

    const BilinearFormIntegrator * Evaluator_ptr;
    const BilinearFormIntegrator * BoundaryEvaluator_ptr;
    

    const FESpace & fes = gfu->GetFESpace();

    const int components = fes.GetEvaluator()->DimFlux();
    int ndomains;
    
    string typestring;

    for(int count=0; count < 2; count++)
      {
	if(count == 0)
	  {
	    if(!volanalyze) continue;

	    typestring = ".vol";

	    Evaluator_ptr = fes.GetEvaluator();
	    BoundaryEvaluator_ptr = NULL;
	    ndomains = pde.GetMeshAccess().GetNDomains();

	  }
	if(count == 1)
	  {
	    if(volanalyze) continue;

	    typestring = ".surf";

	    Evaluator_ptr = NULL;
	    BoundaryEvaluator_ptr = fes.GetBoundaryEvaluator();
	    ndomains = pde.GetMeshAccess().GetNBoundaries();
	  }

	Array<int> & domains = ((count == 0) ? voldomains : surfdomains);


	if(domains.Size() == 0)
	  {
	    domains.SetSize(ndomains);
	    for(int i=0; i<domains.Size(); i++)
	      {
		domains[i] = i+1;
	      }
	  }

	// const FESpace & fes = gfu->GetFESpace();

	VisualizeGridFunction<double> vgfu(pde.GetMeshAccess(),gfu,
					   BoundaryEvaluator_ptr,
					   Evaluator_ptr,false);


	Array<double> mini, maxi, average, vol;

	if(component == -1)
	  {
	    mini.SetSize(components*ndomains);
	    maxi.SetSize(components*ndomains);
	    average.SetSize(components*ndomains);
	  }
	else
	  {
	    mini.SetSize(ndomains);
	    maxi.SetSize(ndomains);
	    average.SetSize(ndomains);
	  }
	  

	if(nodistinction)
	  {
	    vol.SetSize(ndomains);
	    vgfu.Analyze(mini,maxi,average,vol,component);
	  }
	else
	  {
	    vgfu.Analyze(mini,maxi,average,component);
	  }	

	(*testout) << variablename << "min " << mini << endl
		   << variablename << "max " << maxi << endl;


	cout << endl << gfu->GetName();

	if(nodistinction)
	  {
	    if (count == 0) cout << " (volume)" << endl;
	    else cout << " (surface)" << endl;

	    double gmin(1e100), gmax(-1e100), gav(0), gvol(0);
	    
	    for(int i=0; i<domains.Size(); i++)
	      {
		dom = domains[i]-1;
		gvol += vol[dom];
		
		if(component == -1)
		  {
		    for(int j=0; j<components; j++)
		      {
			if(mini[dom*components+j] < gmin) gmin = mini[dom*components+j];
			if(maxi[dom*components+j] > gmax) gmax = maxi[dom*components+j];
			gav += average[(dom)*components+j];
		      }
		  }
		else
		  {
		    if(mini[dom] < gmin) gmin = mini[dom];
		    if(maxi[dom] > gmax) gmax = maxi[dom];
		    gav += average[dom];
		  }
	      }
	    gav /= gvol;

	    cout << "min: " << gmin << endl
		 << "max: " << gmax << endl
		 << "av:  " << gav << endl;
	    if(writevar)
	      {
		ostringstream avnmin,avnmax,avnav;

		avnmin << variablename << ".min";
		actvarname = avnmin.str();
		if(!pde.VariableUsed(actvarname)) pde.AddVariable(actvarname,gmin);
		else pde.GetVariable(actvarname) = gmin;

		avnmax << variablename << ".max";
		actvarname = avnmax.str();
		if(!pde.VariableUsed(actvarname)) pde.AddVariable(actvarname,gmax);
		else pde.GetVariable(actvarname) = gmax;

		avnav << variablename << ".av";
		actvarname = avnav.str();
		if(!pde.VariableUsed(actvarname)) pde.AddVariable(actvarname,gav);
		else pde.GetVariable(actvarname) = gav;
	      }
	  }
	else
	  {
	    if(component != -1) cout << ", component " << component+1;
	    cout << endl;
	    for(int i=0; i<domains.Size(); i++)
	      {
		dom = domains[i];
		
		if(count == 0) cout << "on domain " << dom<<": " << endl;
		else cout << "on surfacedomain " << dom<<": " << endl;
		
		dom--;

		if(component == -1)
		  {
		    cout << "min:" << endl;
		    for(int j=0; j<components; j++)
		      {
			cout << mini[dom*components+j] << " ";
			if(writevar)
			  {
			    ostringstream avn;
			    avn << variablename << typestring << dom+1 << ".comp" << j+1 << ".min";
			    actvarname = avn.str();
			    if(!pde.VariableUsed(actvarname)) pde.AddVariable(actvarname,mini[dom*components+j]);
			    else pde.GetVariable(actvarname) = mini[dom*components+j];
			  }
		      }
		    cout << endl << "max:" << endl;
		    for(int j=0; j<components; j++)
		      {
			cout << maxi[dom*components+j] << " ";
			if(writevar)
			  {
			    ostringstream avn;
			    avn << variablename << typestring << dom+1 << ".comp" << j+1 << ".max";
			    actvarname = avn.str();
			    if(!pde.VariableUsed(actvarname)) pde.AddVariable(actvarname,maxi[dom*components+j]);
			    else pde.GetVariable(actvarname) = maxi[dom*components+j];
			  }
		      }
		    cout << endl << "av:" << endl;
		    for(int j=0; j<components; j++)
		      {
			cout << average[dom*components+j] << " ";
			if(writevar)
			  {
			    ostringstream avn;
			    avn << variablename << typestring << dom+1 << ".comp" << j+1 << ".av";
			    actvarname = avn.str();
			    if(!pde.VariableUsed(actvarname)) pde.AddVariable(actvarname,average[dom*components+j]);
			    else pde.GetVariable(actvarname) = average[dom*components+j];
			  }
		      }
		    cout << endl;
		  }
		else
		  {
		    cout << "min: " << mini[dom] << " max: " << maxi[dom] << " av: " << average[dom] << endl;
		    if(writevar)
		      {
			ostringstream avn;
			avn << variablename << typestring << dom+1 << ".comp" << component+1;
			actvarname = avn.str()+".min";
			if(!pde.VariableUsed(actvarname)) pde.AddVariable(actvarname,mini[dom]);
			else pde.GetVariable(actvarname) = mini[dom];
			actvarname = avn.str()+".max";
			if(!pde.VariableUsed(actvarname)) pde.AddVariable(actvarname,maxi[dom]);
			else pde.GetVariable(actvarname) = maxi[dom];
			actvarname = avn.str()+".av";
			if(!pde.VariableUsed(actvarname)) pde.AddVariable(actvarname,average[dom]);
			else pde.GetVariable(actvarname) = average[dom];
		      }
		  }
	      }
	  }
      }
  }

  void NumProcAnalyze :: PrintReport (ostream & ost)
  {
    ost << "NumProcAnalyze:" << endl;
  }









  //////////////////////////////////////////////////////////////////////////////


  class NumProcWarn : public NumProc
  {
  protected:
    ///
    string variablename1,variablename2;
    ///
    double val1,val2;
    ///
    bool less,lessorequal,greater,greaterorequal;
    ///
    string text;
    
  public:
    ///
    NumProcWarn (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcWarn();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcWarn (pde, flags);
    }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "Warn";
    }
  };


  
  NumProcWarn :: NumProcWarn (PDE & apde, const Flags & flags)
    : NumProc(apde)
  {
    text = flags.GetStringFlag("text","");

    variablename1 = flags.GetStringFlag("var1","");
    variablename2 = flags.GetStringFlag("var2","");

    val1 = flags.GetNumFlag("val1",0);
    val2 = flags.GetNumFlag("val2",0);

    less = flags.GetDefineFlag("less");
    lessorequal = flags.GetDefineFlag("lessorequal");
    greater = flags.GetDefineFlag("greater");
    greaterorequal = flags.GetDefineFlag("greaterorequal");

    //if(!less && !lessorequal && !greater && !greaterorequal) less = true;
  }

  NumProcWarn :: ~NumProcWarn(){;}

  void NumProcWarn :: PrintDoc (ostream & ost)
  {
    ost <<
      "\n\nNumproc Warn:\n" \
      "--------------------------\n" \
      "Checks if variables fulfill a given inequality\n"\
      "\nFlags:\n"\
      " -var1=<variablename> or -val1=<value>\n"\
      " -var2=<variablename> or -val2=<value>\n"\
      " -less or -lessorequal or -greater or -greaterorequal\n"\
      "\nOptional flags:\n"\
      " -text=<text>\n"\
      "     prints the text \n"\
      "\nIf var1 resp. val1 is less, lessorequal, greater, greaterorequal than var2 resp. val2, "\
      "then a warning is given.\n";
  }

  void NumProcWarn :: Do(LocalHeap & lh)
  {
    double value1,value2;
    ostringstream warnleft, warnright;
    string warnmiddle;

    if(strcmp(variablename1.c_str(),"") != 0)
      { 
	value1 = pde.GetVariable(variablename1);
	warnleft << variablename1 << " ("<< value1 <<")";
      }
    else 
      {
	value1 = val1;
	warnleft << value1;
      }
    if(strcmp(variablename2.c_str(),"") != 0)
      { 
	value2 = pde.GetVariable(variablename2);
	warnright << variablename2 << " ("<< value2 <<")";
      }
    else 
      {
	value2 = val2;
	warnright << value2;
      }

    bool warn;

    if(less)
      {
	warn = (value1 < value2);
	warnmiddle = " < ";
      }
    else if(lessorequal)
      {
	warn = (value1 <= value2);
	warnmiddle = " <= ";
      }
    else if(greater)
      {
	warn = (value1 > value2);
	warnmiddle = " > ";
      }
    else if(greaterorequal)
      {
	warn = (value1 >= value2);
	warnmiddle = " >= ";
      }
    else
      {
	throw Exception ("no sensible warning-condition");
      }

    if(warn)
      {
	cout << "Warning: " << text << endl << warnleft.str() << warnmiddle << warnright.str() << endl;
	ostringstream tclstring;
	

	tclstring << "printwarning \"" << text << "\\n" << warnleft.str() << warnmiddle << warnright.str() << "\"" << endl;

	char *dummy; dummy = new char[tclstring.str().size()+1];
	strcpy(dummy,tclstring.str().c_str());

	Tcl_Eval(pde.tcl_interpreter,dummy);

	delete [] dummy;
		 
      }
  }


  void NumProcWarn :: PrintReport (ostream & ost)
  {
    ost << "NumProcWarn:" << endl;
  }



  //////////////////////////////////////////////////////////////////////////////


  class NumProcTclTable : public NumProc
  {
  protected:
    
    int rows, columns;

    Array < string > tableentries;

    string title;

    bool noprint;
        
  public:
    ///
    NumProcTclTable (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcTclTable();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcTclTable (pde, flags);
    }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "TclTable";
    }
  };


  
  NumProcTclTable :: NumProcTclTable (PDE & apde, const Flags & flags)
    : NumProc(apde)
  {
    noprint = flags.GetDefineFlag("noprint");

    rows = static_cast<int>(flags.GetNumFlag("rows",0));
    columns = static_cast<int>(flags.GetNumFlag("columns",0));

    tableentries.SetSize(rows*columns);
    tableentries = "empty";

    title = flags.GetStringFlag("title","");

    const Array<char*> & textarray = flags.GetStringListFlag("entries");

    for(int i=0; i<tableentries.Size() && i<textarray.Size(); i++)
      tableentries[i] = textarray[i];

  }

  NumProcTclTable :: ~NumProcTclTable(){;}

  void NumProcTclTable :: PrintDoc (ostream & ost)
  {
    ost <<
      "\n\nNumproc TclTable:\n" \
      "--------------------------\n" \
      "opens a window with text and variable- and constant-values in a simple table format\n" \
      "  and saves the table at the last position of the tcl variable \"tablesforoutput\"\n\n" \
      "Flags:\n"
      " -rows=<num>\n" \
      "    sets the number of rows\n" \
      " -columns=<num>\n" \
      "    sets the number of columns\n" \
      " -entries=[<e11>,<e12>,...<erc>]\n" \
      "    <eij> can be:\n" \
      "       the keyword \"empty\": then this cell stays empty\n" \
      "       the name of a variable or a constant: then the cell gets its value\n" \
      "       some other text\n\n" \
      "Optional Flags:\n" \
      " -title=<name>\n" \
      "    the window-title\n" \
      " -noprint\n" \
      "    do not print the table now\n\n";
  }

  void NumProcTclTable :: Do(LocalHeap & lh)
  {
    ostringstream tclstring;

    //tclstring << "printtable {\" " << title << "\" " << rows << " " << columns << endl;
    tclstring << "lappend tablesforoutput {\" " << title << "\" " << rows << " " << columns << endl;
    for(int i=0; i<tableentries.Size(); i++)
      {
     if (tableentries[i] == "empty")
#ifdef SCALASCA
		   // scalasca does not compile strings containing " \"\" "
	  tclstring << "\" \" ";
#else  	
	  tclstring << "\"\" ";
#endif
	else if (pde.VariableUsed(tableentries[i]))
	  tclstring << pde.GetVariable(tableentries[i]) << " ";
	else if (pde.ConstantUsed(tableentries[i]))
	  tclstring << pde.GetConstant(tableentries[i]) << " ";
	else
	  tclstring << "\"" << tableentries[i] << "\" ";
      }
    tclstring << "}" << endl;
    
    if(!noprint) tclstring << "printtable [ lindex $tablesforoutput end]" << endl;


    
    char *dummy; dummy = new char[tclstring.str().size()+1];
    strcpy(dummy,tclstring.str().c_str());
    
    //cout << "tclstring: " << endl << dummy << endl;

    Tcl_Eval(pde.tcl_interpreter,dummy);

    delete [] dummy;
		 
  }


  void NumProcTclTable :: PrintReport (ostream & ost)
  {
    ost << "NumProcTclTable:" << endl;
  }




  //////////////////////////////////////////////////////////////////////////////
  
  class NumProcSaveSolution : public NumProc
  {
  protected:
    string filename;
    bool ascii;

  public:
    NumProcSaveSolution (PDE & apde, const Flags & flags);

    
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcSaveSolution (pde, flags);
    }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "SaveSolution";
    }

  };
  
  NumProcSaveSolution :: NumProcSaveSolution (PDE & apde, const Flags & flags)
    : NumProc(apde)
  {
    filename = pde.GetDirectory()+dirslash+flags.GetStringFlag("filename","");
    ascii = flags.GetDefineFlag("ascii");
  }

  
  void NumProcSaveSolution :: Do(LocalHeap & lh)
  {
    if(filename != "")
      pde.SaveSolution(filename,ascii);
  }

  
  void NumProcSaveSolution :: PrintDoc (ostream & ost)
  {
    ost <<
      "\n\nNumproc SaveSolution:\n" \
      "--------------------------\n" \
      "Saves the current solution to a file\n"\
      "\nFlags:\n"\
      " -filename=<name>\n"\
      "      file where to save the solution\n"\
      " -ascii\n"\
      "      the file is not binary\n\n";

  }

  void NumProcSaveSolution :: PrintReport (ostream & ost)
  {
    ost << "NumProcSaveSolution:" << endl;
  }
    

  //////////////////////////////////////////////////////////////////////////////
  
#ifdef ASTRID
  class NumProcSaveZipSolution : public NumProc
  {
  protected:
    string dirname;
    bool ascii;

  public:
    NumProcSaveZipSolution (PDE & apde, const Flags & flags);

    
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcSaveZipSolution (pde, flags);
    }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "SaveSolution";
    }

  };
  
  NumProcSaveZipSolution :: NumProcSaveZipSolution (PDE & apde, const Flags & flags)
    : NumProc(apde)
  {
    dirname = flags.GetStringFlag("filename","solution");
    ascii = flags.GetDefineFlag("ascii");
  }

  
  void NumProcSaveZipSolution :: Do(LocalHeap & lh)
  {
    if ( dirname == "" ) return;
    string suffix = "";
    if ( dirname.size() > 7 )
      suffix = dirname.substr (dirname.length() - 7, 7);
    string filename;
    if ( suffix == ".tar.gz" )
      filename = pde.GetDirectory()+dirslash+dirname;
    else
      filename = pde.GetDirectory()+dirslash+dirname + ".tar.gz";
    pde.SaveZipSolution(filename,ascii);

  }

  
  void NumProcSaveZipSolution :: PrintDoc (ostream & ost)
  {
    ost <<
      "\n\nNumproc SaveZipSolution:\n" \
      "--------------------------\n" \
      "Saves the current solution, mesh, geometry and \n"\
      "pdefile to a .tar.gz file\n"\
      "\nFlags:\n"\
      " -filename=<name>\n"\
      "      file where to save the solution\n"\
      " -ascii\n"\
      "      the solution file is not binary\n\n";

  }

  void NumProcSaveZipSolution :: PrintReport (ostream & ost)
  {
    ost << "NumProcSaveZipSolution:" << endl;
  }
    
#endif


  //////////////////////////////////////////////////////////////////////////////


  class NumProcLoadSolution : public NumProc
  {
  protected:
    string filename;
    bool ascii;

  public:
    NumProcLoadSolution (PDE & apde, const Flags & flags);

    
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcLoadSolution (pde, flags);
    }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "LoadSolution";
    }
  };


  NumProcLoadSolution :: NumProcLoadSolution (PDE & apde, const Flags & flags)
    : NumProc(apde)
  {
    filename = pde.GetDirectory()+dirslash+flags.GetStringFlag("filename","");
    ascii = flags.GetDefineFlag("ascii");
  }

  void NumProcLoadSolution :: Do(LocalHeap & lh)
  {
    if(filename != "")
      pde.LoadSolution(filename,ascii);
  }

  
  void NumProcLoadSolution :: PrintDoc (ostream & ost)
  {
    ost <<
      "\n\nNumproc LoadSolution:\n" \
      "--------------------------\n" \
      "Loads a previously computed solution from a file\n"\
      "\nFlags:\n"\
      " -filename=<name>\n"\
      "      file from where to load the solution\n"\
      " -ascii\n"\
      "      the file is not binary\n\n";

  }

  void NumProcLoadSolution :: PrintReport (ostream & ost)
  {
    ost << "NumProcLoadSolution:" << endl;
  }
    

  class NumProcLoadZipSolution : public NumProc
  {
  protected:
    string dirname;
    bool ascii;

  public:
    NumProcLoadZipSolution (PDE & apde, const Flags & flags);

    
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcLoadZipSolution (pde, flags);
    }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "LoadZipSolution";
    }
  };

#ifdef ASTRID
  NumProcLoadZipSolution :: NumProcLoadZipSolution (PDE & apde, const Flags & flags)
    : NumProc(apde)
  {
    dirname = flags.GetStringFlag("filename","");
    ascii = flags.GetDefineFlag("ascii");
  }

  void NumProcLoadZipSolution :: Do(LocalHeap & lh)
  {
    if ( dirname == "" ) return;
    string suffix = "";
    if ( dirname.size() > 7 )
      suffix = dirname.substr (dirname.length() - 7, 7);
    string filename;
    if ( suffix == ".tar.gz" )
      filename = pde.GetDirectory()+dirslash+dirname;
    else
      filename = pde.GetDirectory()+dirslash+dirname + ".tar.gz";
    pde.LoadZipSolution(filename,ascii);

  }

  
  void NumProcLoadZipSolution :: PrintDoc (ostream & ost)
  {
    ost <<
      "\n\nNumproc LoadZipSolution:\n" \
      "--------------------------\n" \
      "Loads a previously computed pdefile, mesh, geometry\n"\
      "and solution from a .tar.gz file\n"\
      "\nFlags:\n"\
      " -filename=<name>\n"\
      "      .tar.gz file where the data is stored\n"\
      " -ascii\n"\
      "      the solution file is not binary\n\n";

  }

  void NumProcLoadZipSolution :: PrintReport (ostream & ost)
  {
    ost << "NumProcLoadZipSolution:" << endl;
  }
    
#endif

  //////////////////////////////////////////////////////////////////////////////


  class NumProcTclMenu : public NumProc
  {
  protected:


  public:
    ///
    NumProcTclMenu (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcTclMenu();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcTclMenu (pde, flags);
    }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "TclMenu";
    }
  };


  
  NumProcTclMenu :: NumProcTclMenu (PDE & apde, const Flags & flags)
    : NumProc(apde)
  {

    bool newmenu = flags.GetDefineFlag("newmenu");

    string menuname (flags.GetStringFlag("menuname",""));

    string text (flags.GetStringFlag("text",""));


    Array<double> centerpoint;
    bool center = flags.NumListFlagDefined("centerpoint");
    if(center) centerpoint = flags.GetNumListFlag("centerpoint");


    Array<double> rotation_pars; // (alpha, v1, v2, v3)
    bool rotation = flags.NumListFlagDefined("rotation");
    if(rotation) rotation_pars = flags.GetNumListFlag("rotation");
    

    Array<double> clipvec;
    bool clip = flags.NumListFlagDefined("clipvec");
    if(clip) clipvec = flags.GetNumListFlag("clipvec");

    bool noclipsol = flags.GetDefineFlag("noclipsol");


    bool vectorfunction = flags.GetDefineFlag("vectorfunction");
    
    string fieldname (flags.GetStringFlag("fieldname",""));
    int component = static_cast<int>(flags.GetNumFlag("comp",1));
    string evaluate (flags.GetStringFlag("evaluate",""));
    if(evaluate != "") component = 0;


    double deformationscale = flags.GetNumFlag("deformationscale",0);
    bool deformationoff = ( flags.NumFlagDefined("deformationscale") && fabs(deformationscale) < 1.e-6 );
    bool deformationon = ( !deformationoff && flags.NumFlagDefined("deformationscale") );

    double light = flags.GetNumFlag("light",-1);
    if(light > 1) light = 1;

    double minval = 0, maxval = 1;
    bool autoscale = flags.GetDefineFlag("autoscale");
    bool noautoscale = ( flags.NumFlagDefined("minval") && flags.NumFlagDefined("maxval") );
    if(noautoscale)
      {
	minval = flags.GetNumFlag("minval",0);
	maxval = flags.GetNumFlag("maxval",0);
      }

    bool stopsolutiondrawing = flags.GetDefineFlag("stopsolutiondrawing");
    bool solutiondrawing = flags.GetDefineFlag("solutiondrawing");


    int printtable = static_cast<int>(flags.GetNumFlag("printtcltable",0));
    bool printlasttable = flags.GetDefineFlag("printlasttcltable");


    string systemcommand (flags.GetStringFlag("systemcommand",""));
    string systemcommandflag1 (flags.GetStringFlag("systemcommandflag1",""));
    string systemcommandflag2 (flags.GetStringFlag("systemcommandflag2",""));
    string systemcommandflag3 (flags.GetStringFlag("systemcommandflag3",""));

    

    // Build menus

    
    ostringstream tclstring;

    bool ng_setvispar(false), ng_vissetpar(false);

    if(newmenu)
      {
	tclstring << ".ngmenu add cascade -label \"" << text 
		  <<"\" -menu .ngmenu." << menuname <<" -underline 0\n"
		  << "menu .ngmenu." << menuname << endl;
      }
    else
      {
	tclstring << ".ngmenu." << menuname << " add command -label \"" << text << "\" \\" << endl
		  << "-command {" << endl;

	if(stopsolutiondrawing)
	  {
#ifdef SCALASCA
	    tclstring << "set selectvisual \" \"" << endl;
#else
	    tclstring << "set selectvisual \"\"" << endl;
#endif
	    ng_setvispar = true;
	  }

	if(solutiondrawing)
	  {
	    tclstring << "set selectvisual \"solution\"" << endl;
	    ng_setvispar = true;
	  }
	    

	if(center)
	  {
	    for(int i = centerpoint.Size()-1; i<3; i++) centerpoint.Append(0);
	    tclstring << "set viewoptions.usecentercoords 1" << endl
		      << "set viewoptions.centerx " << centerpoint[0] << endl
		      << "set viewoptions.centery " << centerpoint[1] << endl
		      << "set viewoptions.centerz " << centerpoint[2] << endl
		      << "set dummy $selectvisual" << endl
		      << "set selectvisual \"mesh\"" << endl
		      << "Ng_SetVisParameters; Ng_Center" << endl;
	    if(!stopsolutiondrawing)
	      tclstring << "if {$dummy != \"mesh\"} {set selectvisual $dummy; Ng_SetVisParameters}" << endl;
	  }
	if(clip)
	  {
	    for(int i = centerpoint.Size()-1; i<3; i++) clipvec.Append(0);
	    tclstring << "set viewoptions.clipping.enable 1" << endl
		      << "set viewoptions.clipping.nx " << clipvec[0] << endl
		      << "set viewoptions.clipping.ny " << clipvec[1] << endl
		      << "set viewoptions.clipping.nz " << clipvec[2] << endl
		      << "set viewoptions.clipping.dist 0" << endl;
	    if(noclipsol)
	      tclstring << "set visoptions.clipsolution none" << endl;
	    ng_setvispar = true;
	  }
	if(rotation)
	  {
	    for(int i = rotation_pars.Size(); i<4; i++) rotation_pars.Append(0);
	    tclstring << "Ng_ArbitraryRotation";
	    for(int i=0; i<rotation_pars.Size(); i++) tclstring << " " << rotation_pars[i];
	    tclstring << ";" << endl;
	  }


	if(fieldname != "")
	  {
	    if(deformationon)
	      {
		tclstring << "set visoptions.deformation 1" << endl
		      << "set visoptions.scaledeform1 " << deformationscale << endl
		      << "set visoptions.scaledeform2 1" << endl
		      << "set visoptions.vecfunction " << fieldname << endl;
		ng_vissetpar = true;
	      }
	    else
	      {
		if(vectorfunction)
		  {
		    tclstring << "set visoptions.vecfunction " << fieldname << endl;
		    if(clip && !noclipsol)
		      tclstring << "set visoptions.clipsolution vec" << endl;
		    ng_vissetpar = true;
		  }
		else
		  {
		    if(evaluate != "") 
		      tclstring << "set visoptions.evaluate " << evaluate<< endl;
		    tclstring << "set visoptions.scalfunction " << fieldname<< "." << component << endl;
		    if(clip && !noclipsol)
		      tclstring << "set visoptions.clipsolution scal" << endl;
		    ng_vissetpar = true;
		  }
	      }
	  }


	if(deformationoff)
	  {
	    tclstring << "set visoptions.deformation 0" << endl;
	    ng_vissetpar = true;
	  }

	if(light >= 0)
	  {
	    tclstring << "set viewoptions.light.amb " << light << endl;
	    ng_setvispar = true;
	  }


	if(autoscale)
	  {
	    tclstring << "set visoptions.autoscale 1" << endl;
	    ng_vissetpar = true;
	  }

	if(noautoscale)
	  {
	    tclstring << "set visoptions.autoscale 0" << endl
		      << "set visoptions.mminval " << minval << endl
		      << "set visoptions.mmaxval " << maxval << endl;
	    ng_vissetpar = true;
	  }


	if(printtable > 0)
	  tclstring << "printtable [lindex $tablesforoutput " << printtable-1 << "]" << endl;

	if(printlasttable)
	  tclstring << "printtable [lindex $tablesforoutput end]" << endl;

	if(ng_setvispar)
	  tclstring << "Ng_SetVisParameters" << endl;
	if(ng_vissetpar)
	  tclstring << "Ng_Vis_Set parameters" << endl;

	if(systemcommand != "")
	  {
	    tclstring << "exec " << systemcommand;
	    if(systemcommandflag1 != "") tclstring << " " << systemcommandflag1;
	    if(systemcommandflag2 != "") tclstring << " " << systemcommandflag2;
	    if(systemcommandflag3 != "") tclstring << " " << systemcommandflag3;
	    tclstring << " &" <<endl;
	  }
	

	tclstring << "redraw" << endl
		  << "}" << endl;
	//cout << tclstring.str() << endl;
	
      }


    // the following is a lengthy workaround for a visual-c++ problem
    char *dummy; dummy = new char[tclstring.str().size()+1];
    strcpy(dummy,tclstring.str().c_str());

    Tcl_Eval(pde.tcl_interpreter,dummy);
    
    delete [] dummy;
  }

  NumProcTclMenu :: ~NumProcTclMenu(){;}

  void NumProcTclMenu :: PrintDoc (ostream & ost)
  {
    ost <<
      "\n\nNumproc TclMenu:\n" \
      "--------------------------\n" \
      "Adds menu entries for specific visualization\n"\
      "\nFlags:\n"\
      " -menuname=<name>\n"\
      "      menu identifier\n"\
      "\nOptional Flags:\n"\
      " -newmenu\n"\
      "      adds a new menu\n"\
      " -text=<text>\n"\
      "      name of the menu or the menu-entry\n"\
      " -centerpoint=[x,y,z]\n"\
      "      sets the center point\n"\
      " -clipvec=[x,y,z]\n"\
      "      turns on a clipping plane with the given normal-vector and opens the clipping dialog-box\n"\
      " -noclipsol\n"\
      "      no solution is visualized on the clipping plane\n"\
      " -rotation=[alpha1,vx1,vy1,vz1,alpha2,vx2,vy2,vz2,...]\n"\
      "      rotates the standard view by the angle alpha1 round the vector (vx1,vy1,vz1), then by alpha2 round (vx2,vy2,vz2) etc.\n"\
      " -fieldname=<name>\n"\
      "      name of the field to be displayed\n"\
      " -comp=<component>\n"\
      "      component\n"\
      " -vectorfunction\n"\
      "      the field is a vectorfield, and is visualized as such\n"\
      " -evaluate = < abs | abstens | mises | main >\n"
      " -deformationscale = <value>\n"\
      "      display deformation scaled by given value\n"\
      " -light = <value>\n"\
      "      sets the ambient light to a value between 0 and 1\n"\
      " -autoscale\n"\
      "      turns on autoscale\n"\
      " -minval=<value> -maxval=<value>\n"\
      "      turns off autoscale to use the given minimal and maximal values\n" \
      " -printtcltable=<n>\n"\
      "      print table generated by nth instance of numproc tcltable\n" \
      " -printlasttcltable\n"\
      "      print table generated by latest instance of numproc tcltable\n" \
      " -systemcommand=<name>\n"\
      "      a system command is executed, parameters can be set using -systemcommandflagN=<par>, with N=1,2,3\n\n" ;
  }

  void NumProcTclMenu :: Do(LocalHeap & lh)
  {
    
    
  }


  void NumProcTclMenu :: PrintReport (ostream & ost)
  {
    ost << "NumProcTclMenu:" << endl;
  }






  
  //////////////////////////////////////////////////////////////////////////////


  class NumProcGenerateOne : public NumProc
  {
  protected:
    GridFunction * gfone;
    

  public:
     NumProcGenerateOne (PDE & apde, const Flags & flags);

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcGenerateOne (pde, flags);
    }
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "GenerateOne";
    }
  };

  
  NumProcGenerateOne ::  NumProcGenerateOne (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    gfone = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", ""), 0); 
  }

  void NumProcGenerateOne ::  PrintDoc (ostream & ost)

  {
    ost <<
      "\n\nNumproc GenerateOne:\n";
  }


  void NumProcGenerateOne :: Do(LocalHeap & lh)
  {
    const HCurlHighOrderFESpace * hcurlhofespace = dynamic_cast<const HCurlHighOrderFESpace *>(&gfone->GetFESpace());
    const H1HighOrderFESpace * h1hofespace = dynamic_cast<const H1HighOrderFESpace *>(&gfone->GetFESpace());
    const NodalFESpace * h1fespace = dynamic_cast<const NodalFESpace *>(&gfone->GetFESpace());
    const L2HighOrderFESpace * l2hofespace = dynamic_cast<const L2HighOrderFESpace *>(&gfone->GetFESpace());
    
    S_GridFunction<double> * gfoned = dynamic_cast<S_GridFunction<double> *>(gfone);
    S_GridFunction<Complex> * gfonec = dynamic_cast<S_GridFunction<Complex> *>(gfone);

    if(hcurlhofespace)
      {

      }
    else if(h1hofespace || h1fespace)
      {
	Array<int> dnums(1);
	Vector<Complex> vecc(1);
	Vector<double> vecd(1);
	
	vecd(0) = 1.;
	vecc(0) = 1.;

	for(dnums[0]=0; dnums[0] < ma.GetNV(); dnums[0]++)
	  {
	    if(gfoned)
	      gfoned->SetElementVector(dnums,vecd);
	    else if(gfonec)
	      gfonec->SetElementVector(dnums,vecc);
	  }
      }
    else if(l2hofespace)
      {

      }
    else
      {
	throw Exception("NumProcGenerateOne: gridfunction has unknown fespace.");
      }


  }

  void NumProcGenerateOne  :: PrintReport (ostream & ost)
  {
    ost << "NumProcGenerateOne:" << endl;
  }

  //////////////////////////////////////////////////////////////////////////////






  /////////////////////////////////////
  //   num proc visualization
  /////////////////////////////////////


  class NumProcVisualization : public NumProc
  {
  protected:


  public:
    ///
    NumProcVisualization (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcVisualization();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcVisualization (pde, flags);
    }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "Visualization";
    }
  };


  
  NumProcVisualization :: NumProcVisualization (PDE & apde, const Flags & flags)
    : NumProc(apde)
  {

    if ( id > 0 ) return;

    Array<double> centerpoint;
    bool center = flags.NumListFlagDefined("centerpoint");
    if(center) centerpoint = flags.GetNumListFlag("centerpoint");


    Array<double> rotation_pars; // (alpha, v1, v2, v3)
    bool rotation = flags.NumListFlagDefined("rotation");
    if(rotation) rotation_pars = flags.GetNumListFlag("rotation");
    

    Array<double> clipvec;
    bool clip = flags.NumListFlagDefined("clipvec");
    if(clip) clipvec = flags.GetNumListFlag("clipvec");

    string scalarfun (flags.GetStringFlag("scalarfunction",""));
    int scalarcomp = int(flags.GetNumFlag("comp",1));
    string vecfun  (flags.GetStringFlag("vectorfunction",""));
    string clipfun (flags.GetStringFlag("clipsolution",""));
    double clipdist = flags.GetNumFlag("clipdist", 0);

    string evaluate (flags.GetStringFlag("evaluate",""));
    if(evaluate != "") scalarcomp = 0;

    double deformationscale = flags.GetNumFlag("deformationscale",0);
    bool deformationoff = ( flags.NumFlagDefined("deformationscale") && fabs(deformationscale) < 1.e-6 );
    bool deformationon = ( !deformationoff && flags.NumFlagDefined("deformationscale") );

    double light = flags.GetNumFlag("light",-1);
    if(light > 1) light = 1;

    double minval = 0, maxval = 1;
    bool autoscale = flags.GetDefineFlag("autoscale");
    bool noautoscale = ( flags.NumFlagDefined("minval") && flags.NumFlagDefined("maxval") );
    if(noautoscale)
      {
	minval = flags.GetNumFlag("minval",0);
	maxval = flags.GetNumFlag("maxval",0);
      }

    bool stopsolutiondrawing = flags.GetDefineFlag("stopsolutiondrawing");
    bool solutiondrawing = flags.GetDefineFlag("solutiondrawing");


    int printtable = static_cast<int>(flags.GetNumFlag("printtcltable",0));
    bool printlasttable = flags.GetDefineFlag("printlasttcltable");


    string systemcommand (flags.GetStringFlag("systemcommand",""));
    string systemcommandflag1 (flags.GetStringFlag("systemcommandflag1",""));
    string systemcommandflag2 (flags.GetStringFlag("systemcommandflag2",""));
    string systemcommandflag3 (flags.GetStringFlag("systemcommandflag3",""));

    // HERBERT:
    int subdivision = (int) flags.GetNumFlag("subdivision",1);
    bool notexture = (bool) flags.GetDefineFlag("notexture");
    bool nooutline = (bool) flags.GetDefineFlag("nooutline");

    bool nolineartexture = (bool) flags.GetDefineFlag("nolineartexture");
    // Build menus

    
    
    ostringstream tclstring;

    bool ng_setvispar(false), ng_vissetpar(false);


    if(stopsolutiondrawing)
      {
#ifdef SCALASCA
	 	    tclstring << "set ::selectvisual \" \"" << endl;
#else
		    tclstring << "set ::selectvisual \"\"" << endl;
#endif
	ng_setvispar = true;
      }

    if(solutiondrawing)
      {
        tclstring << "set ::selectvisual \"solution\"" << endl;
	ng_setvispar = true;
      }
	    

    if(center)
      {
	for(int i = centerpoint.Size()-1; i<3; i++) centerpoint.Append(0);
        tclstring << "set ::viewoptions.usecentercoords 1" << endl
            << "set ::viewoptions.centerx " << centerpoint[0] << endl
            << "set ::viewoptions.centery " << centerpoint[1] << endl
            << "set ::viewoptions.centerz " << centerpoint[2] << endl
            << "set ::dummy $selectvisual" << endl
		  << "set selectvisual \"mesh\"" << endl
		  << "Ng_SetVisParameters; Ng_Center" << endl;
	if(!stopsolutiondrawing)
          tclstring << "if {$dummy != \"mesh\"} {set ::selectvisual $dummy; Ng_SetVisParameters}" << endl;
      }
    if(clip)
      {
	for(int i = centerpoint.Size()-1; i<3; i++) clipvec.Append(0);
        tclstring << "set ::viewoptions.clipping.enable 1" << endl
            << "set ::viewoptions.clipping.nx " << clipvec[0] << endl
            << "set ::viewoptions.clipping.ny " << clipvec[1] << endl
            << "set ::viewoptions.clipping.nz " << clipvec[2] << endl
		  << "set ::viewoptions.clipping.dist " << clipdist << endl;
	ng_setvispar = true;
      }
    if(rotation)
      {
	for(int i = rotation_pars.Size(); i<4; i++) rotation_pars.Append(0);
	tclstring << "Ng_ArbitraryRotation";
	for(int i=0; i<rotation_pars.Size(); i++) tclstring << " " << rotation_pars[i];
	tclstring << ";" << endl;
      }
    

      // herbert: allow deformation for all sorts of functions!
      if(deformationon)
        {
          tclstring << "set ::visoptions.deformation 1" << endl
              << "set ::visoptions.scaledeform1 " << deformationscale << endl
              << "set ::visoptions.scaledeform2 1" << endl;
        }
      else
        tclstring << "set ::visoptions.deformation 0" << endl;
    
      if(vecfun != "")
        {
          if(deformationon)
            {
              tclstring << "set ::visoptions.deformation 1" << endl
                        << "set ::visoptions.scaledeform1 " << deformationscale << endl
                        << "set ::visoptions.scaledeform2 1" << endl;
            }
          else 
            tclstring << "set ::visoptions.showsurfacesolution 1" << endl;
          
          tclstring  << "set ::visoptions.vecfunction " << vecfun << endl;
          ng_vissetpar = true;
	  
      }
	
    if(scalarfun != "")
      {
        tclstring << "set ::visoptions.scalfunction " << scalarfun << "." << scalarcomp << endl;
	ng_vissetpar = true;
	  
      }
	
    if(evaluate != "") 
      tclstring << "set ::visoptions.evaluate " << evaluate<< endl;
    if(clipfun == "scalar")
      tclstring << "set ::visoptions.clipsolution scal" << endl;
    else if (clipfun == "vector")
      tclstring << "set ::visoptions.clipsolution vec" << endl;
	ng_vissetpar = true;
      
 


    if(deformationoff)
      {
        tclstring << "set ::visoptions.deformation 0" << endl;
	ng_vissetpar = true;
      }
    
    if(light >= 0)
      {
        tclstring << "set ::viewoptions.light.amb " << light << endl;
	ng_setvispar = true;
      }
    
    
    if(autoscale)
      {
        tclstring << "set ::visoptions.autoscale 1" << endl;
	ng_vissetpar = true;
      }
    
    if(noautoscale)
      {
        tclstring << "set ::visoptions.autoscale 0" << endl
            << "set ::visoptions.mminval " << minval << endl
            << "set ::visoptions.mmaxval " << maxval << endl;
	ng_vissetpar = true;
      }
    
    
    if(printtable > 0)
      tclstring << "printtable [lindex $::tablesforoutput " << printtable-1 << "]" << endl;
    
    if(printlasttable)
      tclstring << "printtable [lindex $::tablesforoutput end]" << endl;
    
    // HERBERT:
    tclstring << "set ::visoptions.subdivisions " << subdivision <<  endl;
    tclstring << "set ::visoptions.usetexture " << (notexture ? 0 : 1) <<  endl;
    tclstring << "set ::viewoptions.drawoutline " << (nooutline ? 0 : 1) <<  endl;

    tclstring << "set ::visoptions.lineartexture " << ( nolineartexture ? 0 : 1 ) << endl;
 
    if(ng_setvispar)
      tclstring << "Ng_SetVisParameters" << endl;
    if(ng_vissetpar)
      tclstring << "Ng_Vis_Set parameters" << endl;
    
    if(systemcommand != "")
      {
	tclstring << "exec " << systemcommand;
	if(systemcommandflag1 != "") tclstring << " " << systemcommandflag1;
	if(systemcommandflag2 != "") tclstring << " " << systemcommandflag2;
	if(systemcommandflag3 != "") tclstring << " " << systemcommandflag3;
	tclstring << " &" <<endl;
      }
      
     
    
    
    tclstring << "redraw" << endl
	      << "}" << endl;
    //cout << tclstring.str() << endl;
    
  
    // the following is a lengthy workaround for a visual-c++ problem
    char *dummy; dummy = new char[tclstring.str().size()+1];
    strcpy(dummy,tclstring.str().c_str());
    
    Tcl_Eval(pde.tcl_interpreter,dummy);
    
    delete [] dummy;
  }
  
  NumProcVisualization :: ~NumProcVisualization(){;}
  
  void NumProcVisualization :: PrintDoc (ostream & ost)
  {
    ost <<
      "\n\nNumproc Visualization:\n" \
      "--------------------------\n" \
      "Adds menu entries for specific visualization\n"\
      "\nFlags:\n"\
      "\nOptional Flags:\n"\
      " -centerpoint=[x,y,z]\n"\
      "      sets the center point\n"\
      " -clipvec=[x,y,z]\n"\
      "      turns on a clipping plane with the given normal-vector and opens the clipping dialog-box\n"\
      " -rotation=[alpha1,vx1,vy1,vz1,alpha2,vx2,vy2,vz2,...]\n"\
      "      rotates the standard view by the angle alpha1 round the vector (vx1,vy1,vz1), then by alpha2 round (vx2,vy2,vz2) etc.\n"\
      " -scalarfunction=<name>\n" \
      "      name of the scalar function to be displayed\n"\
      " -scalarcomp=<value>\n"				   \
      "      component of the scalar function, which is displayed (default 1)\n"\
      " -vectorfunction=<name>\n"\
      "      name of the vector function to be displayes\n" \
      " -clipsolution=< scalar | vector > \n"\
      "      specifies clipping plane solution\n" \
      " -evaluate = < abs | abstens | mises | main >\n"\
      "      evaluation for scalar function, if specified, it overrides the scalar component\n"\
      " -deformationscale = <value>\n"\
      "      display deformation scaled by given value\n"\
      " -light = <value>\n"\
      "      sets the ambient light to a value between 0 and 1\n"\
      " -autoscale\n"\
      "      turns on autoscale\n"\
      " -minval=<value> -maxval=<value>\n"\
      "      turns off autoscale to use the given minimal and maximal values\n" \
      " -printtcltable=<n>\n"\
      "      print table generated by nth instance of numproc tcltable\n" \
      " -printlasttcltable\n"\
      "      print table generated by latest instance of numproc tcltable\n" \
      " -systemcommand=<name>\n"\
      "      a system command is executed, parameters can be set using -systemcommandflagN=<par>, with N=1,2,3\n"\
      " -nolineartexture\n"\
      "      turns off linear texture\n"\
      " -notextures\n"\
      "      turns usetextures off\n"\
      " -subdivision = <value>\n"\
      "      sets subdivision to <value>\n\n" ;
  }

  void NumProcVisualization :: Do(LocalHeap & lh)
  {
    
    
  }


  void NumProcVisualization :: PrintReport (ostream & ost)
  {
    ost << "NumProcVisualization:" << endl;
  }





  class NumProcClearGridFunctions : public NumProc
  {
  protected:
    Array<GridFunction *> gf;

  public:
    NumProcClearGridFunctions (PDE & apde, const Flags & flags) : NumProc(apde)
    {
      if(flags.StringFlagDefined("gridfunction"))
	gf.Append(pde.GetGridFunction(flags.GetStringFlag("gridfunction","")));

      if(flags.StringListFlagDefined("gridfunctions"))
	for(int i=0; i<flags.GetStringListFlag("gridfunctions").Size(); i++)
	  gf.Append(pde.GetGridFunction(flags.GetStringListFlag("gridfunctions")[i]));
    }
    
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcClearGridFunctions (pde, flags);
    }
    
    
    virtual string GetClassName () const
    {
      return " Clear GridFunction";
    }

    static void PrintDoc (ostream & ost)
    {
      ost <<
	"\n\nNumproc cleargridfunctions:\n" \
	"-------------------------------\n" \
	"sets gridfunction-values to zero\n\n"\
	"flags:\n"\
	"-gridfunction=<name>\n"\
	"   gridfunction to clear\n"\
	"-gridfunctions=[<name1>,...,<namen>]\n"\
	"   gridfunctions to clear\n";
    }
    ///
    virtual void Do(LocalHeap & lh)
    {
      for(int i=0; i<gf.Size(); i++)
	gf[i]->GetVector() = 0.;
    }
  };


  //////////////////////////////////////////////////////////////////////////////


  class NumProcQuit : public NumProc
  {
  public:
    NumProcQuit (PDE & apde, const Flags & flags) : NumProc(apde)
    {
      if(flags.GetDefineFlag("immedeately"))
	 exit(0);
    }
    
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcQuit (pde, flags);
    }

    virtual string GetClassName () const
    {
      return " Quit";
    }

    static void PrintDoc (ostream & ost)
    {
      ost <<
	"\n\nNumproc quit:\n" \
	"--------------------------\n" \
	"Quits Netgen/NGSolve\n\n";
    }
    ///
    virtual void Do(LocalHeap & lh)
    {
      char exstr[] = "Ng_Exit\n";
      Tcl_Eval(pde.tcl_interpreter,exstr);
      exit(0);
    }
  };


  //////////////////////////////////////////////////////////////////////////////


//   class NumProcDirectSolverRegion : public NumProc
//   {
//   public:
//     NumProcDirectSolverRegion (PDE & apde, const Flags & flags);

    

//     static NumProc * Create (PDE & pde, const Flags & flags)
//     {
//       return new NumProcDirectSolverRegion (pde, flags);
//     }

    
//     virtual string GetClassName () const
//     {
//       return " DirectSolverRegion";
//     }
//   };


//   NumProcDirectSolverRegion :: NumProcDirectSolverRegion (PDE & apde, const Flags & flags) : 
//     NumProc(apde)
//   {
    
//     FESpace * space = pde.GetFESpace(flags.GetStringFlag ("fespace", NULL));

//     CoefficientFunction * cf = pde.GetCoefficientFunction(flags.GetStringFlag ("coefficient", NULL));

//     space -> SetDirectSolverRegionCoefficientFunction(cf);
//   }




  //////////////////////////////////////////////////////////////////////////////










  NumProcs::NumProcInfo::
  NumProcInfo (const string & aname,
	       NumProc* (*acreator)(PDE & pde, const Flags & flags),
	       void (*aprintdoc) (ostream & ost) )
    : name(aname), creator(acreator), printdoc(aprintdoc)
  {
    ;
  }
  
  NumProcs :: NumProcs ()
  {
    ;
  }

  NumProcs :: ~NumProcs()
  {
    for (int i = 0; i < npa.Size(); i++)
      delete npa[i];
  }
  
  void NumProcs :: 
  AddNumProc (const string & aname,
	      NumProc* (*acreator)(PDE & pde, const Flags & flags),
	      void (*printdoc) (ostream & ost) )
  {
    npa.Append (new NumProcInfo(aname, acreator, printdoc));
  }

  const NumProcs::NumProcInfo * 
  NumProcs::GetNumProc(const string & name)
  {
    for (int i = 0; i < npa.Size(); i++)
      {
	if (name == npa[i]->name)
	  return npa[i];
      }
    return 0;
  }

  void NumProcs :: Print (ostream & ost) const
  {
    ost << endl << "NumProcs:" << endl;
    ost <<         "---------" << endl;
    ost << setw(20) << "Name" << endl;
    for (int i = 0; i < npa.Size(); i++)
      ost << setw(20) << npa[i]->name << endl;
  }



 
  NumProcs & GetNumProcs ()
  {
    static NumProcs nps;
    return nps;
  }


  // standard numprocs:




  namespace numproc_cpp
 {


    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetNumProcs().AddNumProc ("calcflux", NumProcCalcFlux::Create, NumProcCalcFlux::PrintDoc);
      GetNumProcs().AddNumProc ("setvalues", NumProcSetValues::Create, NumProcSetValues::PrintDoc);
      GetNumProcs().AddNumProc ("drawflux", NumProcDrawFlux::Create, NumProcDrawFlux::PrintDoc);
      GetNumProcs().AddNumProc ("evaluate", NumProcEvaluate::Create, NumProcEvaluate::PrintDoc);
      GetNumProcs().AddNumProc ("analyze", NumProcAnalyze::Create, NumProcAnalyze::PrintDoc);
      GetNumProcs().AddNumProc ("warn", NumProcWarn::Create, NumProcWarn::PrintDoc);
      GetNumProcs().AddNumProc ("tcltable", NumProcTclTable::Create, NumProcTclTable::PrintDoc);
      GetNumProcs().AddNumProc ("tclmenu", NumProcTclMenu::Create, NumProcTclMenu::PrintDoc);
      GetNumProcs().AddNumProc ("visualization", NumProcVisualization::Create, NumProcVisualization::PrintDoc);
      GetNumProcs().AddNumProc ("loadsolution", NumProcLoadSolution::Create, NumProcLoadSolution::PrintDoc);
#ifdef ASTRID
      GetNumProcs().AddNumProc ("loadzipsolution", NumProcLoadZipSolution::Create, NumProcLoadZipSolution::PrintDoc);
#endif
      GetNumProcs().AddNumProc ("savesolution", NumProcSaveSolution::Create, NumProcSaveSolution::PrintDoc);
#ifdef ASTRID
      GetNumProcs().AddNumProc ("savezipsolution", NumProcSaveZipSolution::Create, NumProcSaveZipSolution::PrintDoc);
#endif
      GetNumProcs().AddNumProc ("quit", NumProcQuit::Create, NumProcQuit::PrintDoc);
      GetNumProcs().AddNumProc ("generateone", NumProcGenerateOne::Create, NumProcGenerateOne::PrintDoc);
      GetNumProcs().AddNumProc ("cleargridfunctions", NumProcClearGridFunctions::Create, NumProcClearGridFunctions::PrintDoc);
      
      //GetNumProcs().AddNumProc ("directsolverregion", NumProcDirectSolverRegion::Create, NumProcDirectSolverRegion::PrintDoc);
      

      //      GetNumProcs().AddNumProc ("setvisual", NumProcSetVisual::Create);



      /*
      cout << "test lapack" << endl;
      Matrix<double> a(1000), b(1000), c(1000);
      a = 1.0;
      b = 1.0;
      for (int i = 0; i < 1000; i++)
	a(i,i) = b(i,i) = 10;
      c = 0.0;

      for (int i = 0; i < 100; i++)
	{
	  cout << "i = " << i << endl;
	  LapackMultAddABt (a, b, 1.0, c);
	}
      */
    }

    
    Init init;
  }
}
  
