#include <solve.hpp>
#include <ctime>
#include <parallelngs.hpp>

#ifdef VTUNE
#include <libittnotify.h>
#endif

namespace ngsolve
{

  /*
  NumProc :: NumProc (PDE & apde) // , const int acallposition)
    : NGS_Object (apde.GetMeshAccess(), "numproc"), pde(apde)
                     // , callposition(acallposition) 
  { ; }
  */


  /*
  NumProc :: NumProc (PDE & apde, const Flags & flags) // , const int acallposition)
    : NGS_Object (apde.GetMeshAccess(int(flags.GetNumFlag("mesh",1))-1), "numproc"), 
      pde(apde) // , callposition(acallposition) 
  {
    if (flags.StringFlagDefined ("name"))
      SetName (flags.GetStringFlag ("name",""));
  }
  
  NumProc :: ~NumProc()
  {
    ;
  }

  void NumProc :: PrintReport (ostream & ost)
  {
    ost << "Base-class NumProc" << endl;
  }

  void NumProc :: PrintDoc (ostream & ost)
  {
    ost << "No documentation available" << endl;
  }

  */






  /* ***************************** Numproc CalcFlux ************************** */



  ///
  class NumProcCalcFlux : public NumProc
  {
  protected:
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<GridFunction> gfu;
    ///
    shared_ptr<GridFunction> gfflux;
    /// compute flux, not gradient
    bool applyd;
    ///
    // bool useall;
    ///
    int domain;
  public:
    ///
    NumProcCalcFlux (shared_ptr<PDE> apde, const Flags & flags);
    ///
    NumProcCalcFlux (shared_ptr<PDE> apde, shared_ptr<BilinearForm> abfa,
                     shared_ptr<GridFunction> agfu, shared_ptr<GridFunction> agfflux,
                     bool aapplyd);
    ///
    virtual ~NumProcCalcFlux() { ; }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual string GetClassName () const
    {
      return "Calc Flux";
    }


    virtual void PrintReport (ostream & ost) const
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


  NumProcCalcFlux :: NumProcCalcFlux (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL));
    if (bfa->NumIntegrators()==0)
      throw Exception ("bilinearform used for CalcFlux needs at least one integrator");

    gfu = apde->GetGridFunction (flags.GetStringFlag ("solution", NULL));
    gfflux = apde->GetGridFunction (flags.GetStringFlag ("flux", NULL));
    applyd = flags.GetDefineFlag ("applyd");
    // useall = flags.GetDefineFlag ("useall");
    domain = static_cast<int>(flags.GetNumFlag("domain",0))-1;
  }

  NumProcCalcFlux :: NumProcCalcFlux (shared_ptr<PDE> apde,
                     shared_ptr<BilinearForm> abfa,
                     shared_ptr<GridFunction> agfu,
                     shared_ptr<GridFunction> agfflux,
                     bool aapplyd)
      : NumProc(apde), bfa(abfa), gfu(agfu), gfflux(agfflux), applyd(aapplyd), domain(-1)
  {
    if (bfa->NumIntegrators()==0)
      throw Exception ("bilinearform used for CalcFlux needs at least one integrator");
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
    CalcFluxProject (*gfu, *gfflux,
		     bfa->GetIntegrator(0),
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





  /* ***************************** Numproc SetValues ************************** */



  ///
  class NumProcSetValues : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfu;
    shared_ptr<CoefficientFunction> coef;
    VorB vb;
    bool coarsegridonly;
    int component;
    bool print;
  public:
    ///
    NumProcSetValues (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
    {
      gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
      coef = apde->GetCoefficientFunction (flags.GetStringFlag ("coefficient", ""));
      if(flags.GetDefineFlag ("boundary"))
	vb = BND;
      else
	vb = VOL;
      coarsegridonly = flags.GetDefineFlag ("coarsegridonly");
      component = int (flags.GetNumFlag ("component", 0))-1;
      print = flags.GetDefineFlag ("print");

      if (flags.NumFlagDefined ("component"))
	{
	  cerr << "!!!!     numproc setvalues   ... -component   is depreciated and will be removed soon" << endl
	       << "!!!!     please use  -gridfuncion=" << gfu->GetName() << "." << component << " instead" << endl;
	}
    }

    ///
    virtual ~NumProcSetValues() { ; }


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
        "-coarsegridonly\n" \
        "    set values only on coarsest grid \n" 
	  << endl;
    }
    
    ///
    virtual void Do(LocalHeap & lh)
    {
      if (coarsegridonly && ma->GetNLevels() > 1) return;
      shared_ptr<GridFunction> hgfu = gfu;
      if (component != -1)
	hgfu = gfu->GetComponent(component);

      SetValues (coef, *hgfu, vb, 0, lh);
      if (print) 
        *testout << "setvalues result:" << endl << hgfu->GetVector() << endl;
    }

    ///
    virtual string GetClassName () const
    {
      return "SetValues";
    }

    virtual void PrintReport (ostream & ost) const
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
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<GridFunction> gfu;
    /// compute flux, not gradient
    bool applyd;
    ///
    bool useall;

    //BilinearFormIntegrator * bfi2d;
    //BilinearFormIntegrator * bfi3d;
    string label;

  public:
    ///
    NumProcDrawFlux (shared_ptr<PDE> apde, const Flags & flags);

    ///
    NumProcDrawFlux (shared_ptr<BilinearForm> abfa, 
                     shared_ptr<GridFunction> agfu,
                     string alabel,
                     bool aapplyd,
                     bool auseall)
      : bfa(abfa), gfu(agfu), applyd(aapplyd), useall(auseall), label(alabel)
    { 
      Array<shared_ptr<BilinearFormIntegrator>> bfi2d, bfi3d;

      for (int i = 0; i < bfa->NumIntegrators(); i++)
        {
          if ((!bfi3d.Size() || useall) && bfa->GetIntegrator(i)->DimElement() == 3)
            bfi3d.Append(bfa->GetIntegrator(i));
          if ((!bfi2d.Size() || useall) && bfa->GetIntegrator(i)->DimElement() == 2)
            bfi2d.Append(bfa->GetIntegrator(i));
        }

      if (!bfa->GetFESpace()->IsComplex())
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
    
    ///
    virtual ~NumProcDrawFlux() { ; }

    ///
    virtual void Do (LocalHeap & lh);

    static void PrintDoc (ostream & ost);
    ///
    virtual string GetClassName () const
    {
      return "Draw Flux";
    }


    virtual void PrintReport (ostream & ost) const
    {
      ost << GetClassName() << endl;
      if (bfa) ost << "Bilinear-form    = " << bfa->GetName() << endl;
      if (bfa) ost << "Differential-Op  = " << bfa->GetIntegrator(0)->Name() << endl;
      if (gfu) ost << "Gridfunction-In  = " << gfu->GetName() << endl;
      ost << "apply coeffs     = " << applyd << endl;
    }
  };


  NumProcDrawFlux :: NumProcDrawFlux (shared_ptr<PDE> apde, const Flags & flags)
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

    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));

    gfu = apde->GetGridFunction (flags.GetStringFlag ("solution", ""));
    if(!gfu)
      gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));

    applyd = flags.GetDefineFlag ("applyd");
    label = flags.GetStringFlag ("label", "");
    useall = flags.GetDefineFlag ("useall");

    Array<shared_ptr<BilinearFormIntegrator>> bfi2d, bfi3d;

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

    if (!bfa->GetFESpace()->IsComplex())
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


  void NumProcDrawFlux :: Do(LocalHeap & lh)
  {
    // cout << "Num-proc draw flux" << endl;
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




  /* **************************** Numproc DrawCoefficient ********************************* */


  class NumProcDrawCoefficient : public NumProc
  {
  protected:
    netgen::SolutionData * vis;
    shared_ptr<CoefficientFunction> cf;
    string label;

  public:
    ///
    NumProcDrawCoefficient (shared_ptr<PDE> apde, const Flags & flags);
    virtual ~NumProcDrawCoefficient() { ; }

    virtual void Do (LocalHeap & lh) { ; }

    static void PrintDoc (ostream & ost);
    virtual string GetClassName () const {return "Draw Flux";}
    virtual void PrintReport (ostream & ost) const { ; }
  };


  NumProcDrawCoefficient :: NumProcDrawCoefficient (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    cf = apde->GetCoefficientFunction (flags.GetStringFlag ("coefficient", ""));
    label = flags.GetStringFlag ("label", "");

    vis = new VisualizeCoefficientFunction (ma, cf);

    Ng_SolutionData soldata;
    Ng_InitSolutionData (&soldata);
  
    soldata.name = (char*)label.c_str();
    soldata.data = 0;
    soldata.components = cf -> Dimension();
    if (cf->IsComplex()) soldata.components *= 2;
    soldata.iscomplex = cf -> IsComplex();
    soldata.draw_surface = true;
    soldata.draw_volume  = true; 
    if (flags.GetDefineFlag("volume"))
      soldata.draw_surface = false;
    if (flags.GetDefineFlag("boundary"))
      soldata.draw_volume = false;
    soldata.dist = 1;
    soldata.soltype = NG_SOLUTION_VIRTUAL_FUNCTION;
    soldata.solclass = vis;
    Ng_SetSolutionData (&soldata);
  }



  void NumProcDrawCoefficient :: PrintDoc (ostream & ost)

  {
    ost <<
      "\n\nNumproc draw:\n" \
      "-----------------\n" \
      "Draws a coefficient function:\n"\
      "Required flags:\n"
      "-coefficient=<cfname>\n" \
      "    name of the coefficient\n" \
      "-label=<label>\n" \
      "    name displayed in the visualization dialog\n";
  }





  /* **************************** Numproc Evaluate ********************************* */




  ///
  class NumProcEvaluate : public NumProc
  {
  protected:
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<LinearForm> lff;
    ///
    shared_ptr<GridFunction> gfu;
    ///
    shared_ptr<GridFunction> gfv;
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
    NumProcEvaluate (shared_ptr<PDE> apde, const Flags & flags);
    ///
    virtual ~NumProcEvaluate() { ; }

    /*
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcEvaluate (pde, flags);
    }
    */
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost) const; 

    virtual string GetClassName () const
    {
      return "Evaluate";
    }

  };





  NumProcEvaluate :: NumProcEvaluate (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde), point(1), point2(1), point3(1), point4(1)
  {
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", ""), 1); 
    lff = apde->GetLinearForm (flags.GetStringFlag ("linearform", ""), 1);
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""), 0); 
    gfv = apde->GetGridFunction (flags.GetStringFlag ("gridfunction2", ""), 1); 

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
      filename = apde->GetDirectory() + dirslash + flags.GetStringFlag("filename","");
    else
      filename = "err.out";

    

    applyd = flags.GetDefineFlag ("applyd");
    hermitsch = flags.GetDefineFlag ("hermitsch");

    outputprecision = (apde->ConstantUsed("outputprecision")) ? int(apde->GetConstant("outputprecision")) : -1;
    if(flags.NumFlagDefined("outputprecision"))
      outputprecision = int(flags.GetNumFlag("outputprecision",-1));

    component = static_cast<int>(flags.GetNumFlag("cachecomp",1))-1;
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
      "       -n3=<number of points>\n"                                 \
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

	cout << IM(1) << "<" << lff->GetName() << ", " << gfu->GetName() << "> = " << flush;
	if (!lff->GetFESpace()->IsComplex())
	  {
	    result = S_InnerProduct<double>(lff->GetVector(), gfu->GetVector());
	    cout << IM(1) << result << endl;
	  }
	else
	  cout << IM(1) << S_InnerProduct<Complex>(lff->GetVector(), gfu->GetVector()) << endl;
      }
    else if (point.Size() >= 2)
      {
	auto bfi = (bfa) ? bfa->GetIntegrator(0) : gfu->GetFESpace()->GetIntegrator();

	if (point2.Size() >= 2)
	  {
	    shared_ptr<PDE>(pde)->GetEvaluateFiles() += " ";
	    shared_ptr<PDE>(pde)->GetEvaluateFiles() += filename;
	    shared_ptr<PDE>(pde)->GetEvaluateFiles() += " ";
	    shared_ptr<PDE>(pde)->GetEvaluateFiles() += text;
	    
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
		    
		    if (!gfu->GetFESpace()->IsComplex())
		      {
			FlatVector<double> pflux(bfi->DimFlux(), lh);
			bool ok =
			  CalcPointFlux<double> (*gfu, p,
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
			FlatVector<Complex> pflux(bfi->DimFlux(), lh);
			bool ok =
			  CalcPointFlux (*gfu, p,
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
			    ma->FindElementOfPoint(p,dummyip,true);
			    firstone = false;
			  }
			
			if (!gfu->GetFESpace()->IsComplex())
			  {
			    FlatVector<double> pflux(bfi->DimFlux(), lh);
                            CalcPointFlux (*gfu, p,
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
			    FlatVector<Complex> pflux(bfi->DimFlux(), lh);
			    complexflux = true;
                            CalcPointFlux (*gfu, p,
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

	    if (!gfu->GetFESpace()->IsComplex())
	      {
		FlatVector<double> pflux(bfi->DimFlux(), lh);
		CalcPointFlux (*gfu, point, domains,
			       pflux, bfi, applyd, lh, component);
		for (int i = 0; i < pflux.Size(); i++)
		  cout << pflux(i) << "   ";
		cout << endl;
	      }
	    else
	      {
		FlatVector<Complex> pflux(bfi->DimFlux(), lh);
		CalcPointFlux (*gfu, point, domains,
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
	auto hv = vecu.CreateVector();
	bfa->GetMatrix().Mult (vecu, *hv);
	if (!bfa->GetFESpace()->IsComplex())
	  {
	    result = S_InnerProduct<double>(vecv, *hv);
	    //cout << setprecision(16) << result << endl;

	    cout << "bf(gf,gf2) = " << result << endl;
	    ofile << "err = " << result << endl;
	    //cout << "bf(gf,gf2) = " << setprecision(16) << result << endl;
	    //ofile << "err = " << setprecision(16) << result << endl;

	  }
	else
	  {
	    if (!hermitsch)
	      cout << S_InnerProduct<Complex>(vecv, *hv) << endl;
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
      }
    
    shared_ptr<PDE>(pde)->GetVariable(variablename,true) = result;
    
    cout.precision(old_cout_precision);
    ofile.precision(old_ofile_precision);
    ofile.close();
  }

  void NumProcEvaluate :: PrintReport (ostream & ost) const
  {
    ost << "NumProcEvaluate:" << endl;
  }






  /////////////////////////////////////////////////////////////////////////////
  ///
  class NumProcAnalyze : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfu;
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
    NumProcAnalyze (shared_ptr<PDE> apde, const Flags & flags);
    ///
    virtual ~NumProcAnalyze();

    /*
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcAnalyze (pde, flags);
    }
    */
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost) const;

    virtual string GetClassName () const
    {
      return "Analyze";
    }

  };





  NumProcAnalyze :: NumProcAnalyze (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""), 0); 

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

    shared_ptr<BilinearFormIntegrator> Integrator_ptr;
    shared_ptr<BilinearFormIntegrator> BoundaryIntegrator_ptr;
    

    const FESpace & fes = *gfu->GetFESpace();

    const int components = fes.GetIntegrator()->DimFlux();
    int ndomains;
    
    string typestring;

    for(int count=0; count < 2; count++)
      {
	if(count == 0)
	  {
	    if(!volanalyze) continue;

	    typestring = ".vol";

	    Integrator_ptr = fes.GetIntegrator();
	    BoundaryIntegrator_ptr = NULL;
	    ndomains = shared_ptr<PDE>(pde)->GetMeshAccess()->GetNDomains();

	  }
	if(count == 1)
	  {
	    if(volanalyze) continue;

	    typestring = ".surf";

	    Integrator_ptr = NULL;
	    BoundaryIntegrator_ptr = fes.GetIntegrator(BND);
	    ndomains = shared_ptr<PDE>(pde)->GetMeshAccess()->GetNBoundaries();
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

	VisualizeGridFunction<double> vgfu(shared_ptr<PDE>(pde)->GetMeshAccess(),gfu,
					   BoundaryIntegrator_ptr,
					   Integrator_ptr,false);


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
		if(!shared_ptr<PDE>(pde)->VariableUsed(actvarname)) shared_ptr<PDE>(pde)->AddVariable(actvarname,gmin);
		else shared_ptr<PDE>(pde)->GetVariable(actvarname) = gmin;

		avnmax << variablename << ".max";
		actvarname = avnmax.str();
		if(!shared_ptr<PDE>(pde)->VariableUsed(actvarname)) shared_ptr<PDE>(pde)->AddVariable(actvarname,gmax);
		else shared_ptr<PDE>(pde)->GetVariable(actvarname) = gmax;

		avnav << variablename << ".av";
		actvarname = avnav.str();
		if(!shared_ptr<PDE>(pde)->VariableUsed(actvarname)) shared_ptr<PDE>(pde)->AddVariable(actvarname,gav);
		else shared_ptr<PDE>(pde)->GetVariable(actvarname) = gav;
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
			    if(!shared_ptr<PDE>(pde)->VariableUsed(actvarname)) shared_ptr<PDE>(pde)->AddVariable(actvarname,mini[dom*components+j]);
			    else shared_ptr<PDE>(pde)->GetVariable(actvarname) = mini[dom*components+j];
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
			    if(!shared_ptr<PDE>(pde)->VariableUsed(actvarname)) shared_ptr<PDE>(pde)->AddVariable(actvarname,maxi[dom*components+j]);
			    else shared_ptr<PDE>(pde)->GetVariable(actvarname) = maxi[dom*components+j];
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
			    if(!shared_ptr<PDE>(pde)->VariableUsed(actvarname)) shared_ptr<PDE>(pde)->AddVariable(actvarname,average[dom*components+j]);
			    else shared_ptr<PDE>(pde)->GetVariable(actvarname) = average[dom*components+j];
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
			if(!shared_ptr<PDE>(pde)->VariableUsed(actvarname)) shared_ptr<PDE>(pde)->AddVariable(actvarname,mini[dom]);
			else shared_ptr<PDE>(pde)->GetVariable(actvarname) = mini[dom];
			actvarname = avn.str()+".max";
			if(!shared_ptr<PDE>(pde)->VariableUsed(actvarname)) shared_ptr<PDE>(pde)->AddVariable(actvarname,maxi[dom]);
			else shared_ptr<PDE>(pde)->GetVariable(actvarname) = maxi[dom];
			actvarname = avn.str()+".av";
			if(!shared_ptr<PDE>(pde)->VariableUsed(actvarname)) shared_ptr<PDE>(pde)->AddVariable(actvarname,average[dom]);
			else shared_ptr<PDE>(pde)->GetVariable(actvarname) = average[dom];
		      }
		  }
	      }
	  }
      }
  }

  void NumProcAnalyze :: PrintReport (ostream & ost) const
  {
    ost << "NumProcAnalyze:" << endl;
  }



  class NumProcIntegrate : public NumProc
  {
    shared_ptr<CoefficientFunction> coef;
    int order;
  public:
    NumProcIntegrate (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde, flags)
    {
      order = int (flags.GetNumFlag ("order", 2));
      coef = apde->GetCoefficientFunction (flags.GetStringFlag ("coefficient", "") );

      if (!coef->IsComplex())
	apde->AddVariable (string("integrate.")+GetName()+".value", 0.0, 6);
      else
	{
	  apde->AddVariable (string("integrate.")+GetName()+".value.real", 0.0, 6);
	  apde->AddVariable (string("integrate.")+GetName()+".value.imag", 0.0, 6);
	}
    }

    virtual string GetClassName () const
    {
      return "Integrate";
    }

    static void PrintDoc (ostream & ost)
    {
      ost <<
	"\n\nNumproc integrate:\n"					\
	"--------------------------\n"					\
	"\nIntegrates a coefficient function\n"				\
	"\nRequired flags:\n"						\
	"-coefficient=<cfname>\n"					\
	"    coefficientfunction to integrate\n"					\
	"\nOptional flags:\n"						\
	"-order\n"							\
	"    integration order\n"; 
    }


    template <typename SCAL>
    SCAL DoScal (LocalHeap & lh)
    {
      mutex npintegrate_mutex;
      SCAL sum = 0;
      cout << "np integrate,ne = " << ma->GetNE() << endl;
      ParallelForRange( IntRange(ma->GetNE()), [&] ( IntRange r )
      {
	LocalHeap slh = lh.Split(), &lh = slh;
	SCAL lsum = 0;
	for (int i : r)
	  {
	    HeapReset hr(lh);
            ElementId ei(VOL, i);
	    ElementTransformation & eltrans = ma->GetTrafo (ei, lh);
	    IntegrationRule ir (eltrans.GetElementType(), order);
	    const BaseMappedIntegrationRule & mir = eltrans(ir, lh);
	      
	    FlatMatrix<SCAL> result(mir.Size(), 1, lh);
	    coef -> Evaluate (mir, result);
	    SCAL hsum = 0;
	    for (int j = 0; j < mir.Size(); j++)
	      hsum += mir[j].GetWeight() * result(j);

	    lsum += hsum;
	  }
	{
          lock_guard<mutex> guard(npintegrate_mutex);
	  sum += lsum;
	}
      });

      sum = MyMPI_AllReduce (sum);
      return sum;
    }

    virtual void Do (LocalHeap & lh)
    {
      if (!coef -> IsComplex())
	{
	  double sum = DoScal<double> (lh);
	  cout << IM(1) << "Integral = " << sum << endl;
	  shared_ptr<PDE>(pde)->AddVariable (string("integrate.")+GetName()+".value", sum, 6);
	}
      else
	{
	  Complex sum = DoScal<Complex> (lh);
	  cout << IM(1) << "Integral = " << sum << endl;
	  shared_ptr<PDE>(pde)->AddVariable (string("integrate.")+GetName()+".value.real", sum.real(), 6);
	  shared_ptr<PDE>(pde)->AddVariable (string("integrate.")+GetName()+".value.imag", sum.imag(), 6);
	}
    }
  };



  class NumProcWriteFile : public NumProc
  {
    ofstream * outfile;
    int outputprecision;
    Array<string> output_vars;
  public:
    NumProcWriteFile (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    {
      outfile = NULL;

      string filename = flags.GetStringFlag ("filename","");

      outputprecision = (apde->ConstantUsed("outputprecision")) ? int(apde->GetConstant("outputprecision")) : -1;
      if(flags.NumFlagDefined("outputprecision"))
        outputprecision = int(flags.GetNumFlag("outputprecision",-1));
      
      if (filename.length() && (MyMPI_GetId() == 0) )
	{
	  filename = apde->GetDirectory() + dirslash + filename;
	  cout << "NP WriteFile: outputfile is " << filename << endl;
          if (!flags.GetDefineFlag ("append"))
            outfile = new ofstream (filename.c_str());
          else
            outfile = new ofstream (filename.c_str(), ios_base::app);

          if(outputprecision > 0)
            outfile->precision(outputprecision);
	}
      else
	outfile = 0;

      output_vars = flags.GetStringListFlag ("variables");
      // cout << "variables = " << endl << output_vars;

      if (outfile && !flags.GetDefineFlag("append"))
	{
	  *outfile << "# ";
	  for (int i = 0; i < output_vars.Size(); i++)
	    *outfile << output_vars[i] << " ";
	  *outfile << endl;
	}
    }

    ~NumProcWriteFile () 
    {
      delete outfile;
    }

    virtual string GetClassName () const
    {
      return "WriteFile";
    }

    static void PrintDoc (ostream & ost)
    {
      ost <<
	"\n\nNumproc writefile:\n"					\
	"--------------------------\n"					\
	"\nDumps variables/constants to a file\n"			\
	"\nRequired flags:\n"						\
	"-filename=<name>\n"						\
	"    specify the filename\n"					\
	"-append\n"                                                     \
	"    append to file\n"                                          \
	"-variables=[var1,var2...]\n"					\
	"    variables for output\n";					
    }

    virtual void Do (LocalHeap & lh)
    {
      for (int i = 0; i < output_vars.Size(); i++)
	{
	  if (shared_ptr<PDE>(pde)->StringConstantUsed (output_vars[i]))
            {
              string sval = shared_ptr<PDE>(pde)->GetStringConstant(output_vars[i]);
              cout << IM(3) << output_vars[i] << " = " << sval << endl;
              if (outfile) *outfile << sval << " ";
            }
          else
            {
              double val = -1e99;
              if (shared_ptr<PDE>(pde)->ConstantUsed (output_vars[i])) val = shared_ptr<PDE>(pde)->GetConstant(output_vars[i]);
              if (shared_ptr<PDE>(pde)->VariableUsed (output_vars[i])) val = shared_ptr<PDE>(pde)->GetVariable(output_vars[i]);
              cout << IM(3) << output_vars[i] << " = " << val << endl;
              if (outfile) *outfile << val << " ";
            }
	}
      if (outfile) *outfile << endl;
    }
  };






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
    NumProcWarn (shared_ptr<PDE> apde, const Flags & flags);
    ///
    virtual ~NumProcWarn();

    /*
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcWarn (pde, flags);
    }
    */
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost) const;

    virtual string GetClassName () const
    {
      return "Warn";
    }
  };


  
  NumProcWarn :: NumProcWarn (shared_ptr<PDE> apde, const Flags & flags)
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
	value1 = shared_ptr<PDE>(pde)->GetVariable(variablename1);
	warnleft << variablename1 << " ("<< value1 <<")";
      }
    else 
      {
	value1 = val1;
	warnleft << value1;
      }
    if(strcmp(variablename2.c_str(),"") != 0)
      { 
	value2 = shared_ptr<PDE>(pde)->GetVariable(variablename2);
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

	shared_ptr<PDE>(pde)->Tcl_Eval(tclstring.str());

	delete [] dummy;
		 
      }
  }


  void NumProcWarn :: PrintReport (ostream & ost) const
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
    NumProcTclTable (shared_ptr<PDE> apde, const Flags & flags);
    ///
    virtual ~NumProcTclTable();

    /*
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcTclTable (pde, flags);
    }
    */
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost) const;

    virtual string GetClassName () const
    {
      return "TclTable";
    }
  };


  
  NumProcTclTable :: NumProcTclTable (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc(apde)
  {
    noprint = flags.GetDefineFlag("noprint");

    rows = static_cast<int>(flags.GetNumFlag("rows",0));
    columns = static_cast<int>(flags.GetNumFlag("columns",0));

    tableentries.SetSize(rows*columns);
    tableentries = "empty";

    title = flags.GetStringFlag("title","");

    const Array<string> & textarray = flags.GetStringListFlag("entries");

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
	else if (shared_ptr<PDE>(pde)->VariableUsed(tableentries[i]))
	  tclstring << shared_ptr<PDE>(pde)->GetVariable(tableentries[i]) << " ";
	else if (shared_ptr<PDE>(pde)->ConstantUsed(tableentries[i]))
	  tclstring << shared_ptr<PDE>(pde)->GetConstant(tableentries[i]) << " ";
	else
	  tclstring << "\"" << tableentries[i] << "\" ";
      }
    tclstring << "}" << endl;
    
    if(!noprint) tclstring << "printtable [ lindex $tablesforoutput end]" << endl;


    
    char *dummy; dummy = new char[tclstring.str().size()+1];
    strcpy(dummy,tclstring.str().c_str());
    
    //cout << "tclstring: " << endl << dummy << endl;

    shared_ptr<PDE>(pde)->Tcl_Eval(tclstring.str()); // shared_ptr<PDE>(pde)->GetTclInterpreter(), dummy);

    delete [] dummy;
		 
  }


  void NumProcTclTable :: PrintReport (ostream & ost) const
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
    NumProcSaveSolution (shared_ptr<PDE> apde, const Flags & flags);

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost) const; 

    virtual string GetClassName () const
    {
      return "SaveSolution";
    }

  };
  
  NumProcSaveSolution :: NumProcSaveSolution (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc(apde)
  {
    filename = apde->GetDirectory()+dirslash+flags.GetStringFlag("filename","");
    ascii = flags.GetDefineFlag("ascii");
  }

  
  void NumProcSaveSolution :: Do(LocalHeap & lh)
  {
    if(filename != "")
      shared_ptr<PDE>(pde)->SaveSolution(filename,ascii);
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

  void NumProcSaveSolution :: PrintReport (ostream & ost) const
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
    virtual void PrintReport (ostream & ost) const;

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

  void NumProcSaveZipSolution :: PrintReport (ostream & ost) const
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
    NumProcLoadSolution (shared_ptr<PDE> apde, const Flags & flags);

    /*
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcLoadSolution (pde, flags);
    }
    */
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost) const; 

    virtual string GetClassName () const
    {
      return "LoadSolution";
    }
  };


  NumProcLoadSolution :: NumProcLoadSolution (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc(apde)
  {
    filename = apde->GetDirectory()+dirslash+flags.GetStringFlag("filename","");
    ascii = flags.GetDefineFlag("ascii");
  }

  void NumProcLoadSolution :: Do(LocalHeap & lh)
  {
    if(filename != "")
      shared_ptr<PDE>(pde)->LoadSolution(filename,ascii);
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

  void NumProcLoadSolution :: PrintReport (ostream & ost) const
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
    virtual void PrintReport (ostream & ost) const;

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

  void NumProcLoadZipSolution :: PrintReport (ostream & ost) const
  {
    ost << "NumProcLoadZipSolution:" << endl;
  }
    
#endif

  //////////////////////////////////////////////////////////////////////////////






  
  
  
  class NumProcLoadSolution2 : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfu;
    string filename;

  public:
    NumProcLoadSolution2 (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    {
      gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));	
      filename = flags.GetStringFlag("filename", "solution.out");	
    }
    
    virtual ~NumProcLoadSolution2() {;}
    
    virtual void Do(LocalHeap & lh)
    {
      ifstream infile(filename.c_str(), ios::binary);
      gfu -> Load(infile);
    }
    
    virtual string GetClassName () const
    {
      return "NumProcLoadSolution2";
    }
  };


  class NumProcSaveSolution2 : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfu;
    string filename;

  public:
    NumProcSaveSolution2 (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    {
      gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
      filename = flags.GetStringFlag("filename", "solution.out");
    }
    
    virtual ~NumProcSaveSolution2() {};		
    
    virtual void Do(LocalHeap & lh)
    {
      ofstream out(filename.c_str(), ios::binary);
      gfu -> Save(out);
    }
    
    virtual string GetClassName () const
    {
      return "NumProcSaveSolution2";
    }
  };


  class NumProcAssembleLinearization : public NumProc
  {
  protected:
    shared_ptr<BilinearForm> bf;
    shared_ptr<GridFunction> gfu;

  public:
    NumProcAssembleLinearization (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    {
      bf = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", NULL));
      gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""));
    }
    
    virtual ~NumProcAssembleLinearization() {};		
    
    virtual void Do(LocalHeap & lh)
    {
      BaseVector & vecu = gfu -> GetVector();
      cout << " assemble linearization:" << endl;
      bf -> AssembleLinearization(vecu, lh);
    }
    
    virtual string GetClassName () const
    {
      return "NumProcAssembleLinearization";
    }



    static void PrintDoc (ostream & ost)
    {
      ost <<
        "\n\nNumproc AssembleLinearization:\n" \
        "--------------------------\n" \
        "Assembles Bilinearform depending on a gridfunction \n"\
        "\nFlags:\n"\
        " -bilinearform=<name>\n"\
        "      bilinearform that depends on gridfunction\n"\
        " -gridfunction\n"\
        "      gridfunction (the argument of the linearization)\n";
    }
  
  };





  class NumProcTclMenu : public NumProc
  {
  protected:


  public:
    ///
    NumProcTclMenu (shared_ptr<PDE> apde, const Flags & flags);
    ///
    virtual ~NumProcTclMenu();
    
    /*
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcTclMenu (pde, flags);
    }
    */
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost) const;

    virtual string GetClassName () const
    {
      return "TclMenu";
    }
  };


  
  NumProcTclMenu :: NumProcTclMenu (shared_ptr<PDE> apde, const Flags & flags)
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
		    tclstring << "set visoptions.scalfunction " << fieldname<< ":" << component << endl;
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

    apde->Tcl_Eval(tclstring.str());
    
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


  void NumProcTclMenu :: PrintReport (ostream & ost) const
  {
    ost << "NumProcTclMenu:" << endl;
  }






  
  //////////////////////////////////////////////////////////////////////////////


  class NumProcGenerateOne : public NumProc
  {
  protected:
    shared_ptr<GridFunction> gfone;
    

  public:
    NumProcGenerateOne (shared_ptr<PDE> apde, const Flags & flags);
    /*
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcGenerateOne (pde, flags);
    }
    */
    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost) const;

    virtual string GetClassName () const
    {
      return "GenerateOne";
    }
  };

  
  NumProcGenerateOne ::  NumProcGenerateOne (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    gfone = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", ""), 0); 
  }

  void NumProcGenerateOne ::  PrintDoc (ostream & ost)

  {
    ost <<
      "\n\nNumproc GenerateOne:\n";
  }


  void NumProcGenerateOne :: Do(LocalHeap & lh)
  {
    const HCurlHighOrderFESpace * hcurlhofespace = dynamic_cast<const HCurlHighOrderFESpace *>(gfone->GetFESpace().get());
    const H1HighOrderFESpace * h1hofespace = dynamic_cast<const H1HighOrderFESpace *>(gfone->GetFESpace().get());
    const NodalFESpace * h1fespace = dynamic_cast<const NodalFESpace *>(gfone->GetFESpace().get());
    const L2HighOrderFESpace * l2hofespace = dynamic_cast<const L2HighOrderFESpace *>(gfone->GetFESpace().get());
    
    S_GridFunction<double> * gfoned = dynamic_cast<S_GridFunction<double> *>(gfone.get());
    S_GridFunction<Complex> * gfonec = dynamic_cast<S_GridFunction<Complex> *>(gfone.get());

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

	for(dnums[0]=0; dnums[0] < ma->GetNV(); dnums[0]++)
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

  void NumProcGenerateOne  :: PrintReport (ostream & ost) const
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
    NumProcVisualization (shared_ptr<PDE> apde, const Flags & flags);
    ///
    virtual ~NumProcVisualization() { ; }

    static void PrintDoc (ostream & ost);

    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual void PrintReport (ostream & ost) const;

    virtual string GetClassName () const
    {
      return "Visualization";
    }
  };


  
  NumProcVisualization :: NumProcVisualization (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc(apde)
  {

    if ( MyMPI_GetId() != 0 ) return;

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
        tclstring << "set ::visoptions.scalfunction " << scalarfun << ":" << scalarcomp << endl;
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

    apde->Tcl_Eval(tclstring.str());
    
    delete [] dummy;
  }
  
  
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


  void NumProcVisualization :: PrintReport (ostream & ost) const
  {
    ost << "NumProcVisualization:" << endl;
  }





  class NumProcClearGridFunctions : public NumProc
  {
  protected:
    Array<shared_ptr<GridFunction>> gf;

  public:
    NumProcClearGridFunctions (shared_ptr<PDE> apde, const Flags & flags) : NumProc(apde)
    {
      if(flags.StringFlagDefined("gridfunction"))
	gf.Append(apde->GetGridFunction(flags.GetStringFlag("gridfunction","")));

      if(flags.StringListFlagDefined("gridfunctions"))
	for(int i=0; i<flags.GetStringListFlag("gridfunctions").Size(); i++)
	  gf.Append(apde->GetGridFunction(flags.GetStringListFlag("gridfunctions")[i]));
    }
    
    /*
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcClearGridFunctions (pde, flags);
    }
    */
    
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
    NumProcQuit (shared_ptr<PDE> apde, const Flags & flags) : NumProc(apde)
    {
      if(flags.GetDefineFlag("immedeately"))
	 exit(0);
    }
    
    /*
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcQuit (pde, flags);
    }
    */

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
      shared_ptr<PDE>(pde)->Tcl_Eval(exstr);
      exit(0);
    }
  };



class NumProcPause : public NumProc
{
  double seconds;
public:
  NumProcPause (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    seconds = flags.GetNumFlag ("seconds", 10);
  }

  virtual void Do(LocalHeap & lh)
  {
#ifndef WIN32
    sleep (seconds);
#else
	  Sleep(1000*seconds);
#endif
  }


  virtual string GetClassName () const
  {
    return "NumProcPause";
  }

  virtual void PrintReport (ostream & ost) const
  {
    ost << GetClassName() << endl
        << "pause for " << seconds << " seconds" << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc Pause:\n"                
	<< endl;
  }
};


class NumProcVtuneProfiling : public NumProc
{
  bool pause;
  bool resume;
public:
  NumProcVtuneProfiling  (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    pause = flags.GetDefineFlag ("pause");
    resume = flags.GetDefineFlag ("resume");
    if(pause==resume) cout << "NumProcVtuneProfiling: only exactly one flag (pause, resume) allowed!!!" << endl;
  }

  virtual void Do(LocalHeap & lh)
  {
#ifdef VTUNE
      if(pause)
          __itt_pause();
      if(resume)
          __itt_resume();
#else
      if(pause==resume) cout << "NumProcVtuneProfiling: VTUNE not supported, compile with -DUSE_VTUNE=ON" << endl;
#endif // VTUNE
  }


  virtual string GetClassName () const
  {
    return "NumProcVtuneProfiling";
  }

  virtual void PrintReport (ostream & ost) const
  {
    ost << GetClassName() << endl
        << "pause/resume vtune profiling" << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc VtuneProfiling:\n"                
	<< endl;
  }
};



//////////////////////////////////////////////////////////////////////////////


class NumProcTestVariable : public NumProc
{
  string varname;
  Array<double> refvalues;
  double tolerance;
  bool abstol;
  bool cdash;
  int calls = 0;
public:
  NumProcTestVariable (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    varname = flags.GetStringFlag("variable");

    if (flags.NumFlagDefined("refvalue"))
    {
      const double refvalue = flags.GetNumFlag("refvalue",0.0);
      refvalues.Append(refvalue);
    }
    else if (flags.NumListFlagDefined("refvalues"))
    {
      refvalues = flags.GetNumListFlag("refvalues");
    }
    else
    {
      cout << "WARNING: no reference values given, will not compare anything" << endl;
    }
    
    tolerance = flags.GetNumFlag("tolerance",0.0);
    abstol = flags.GetDefineFlag("abstol");
    cdash = flags.GetDefineFlag("cdash");
  }

  virtual void Do(LocalHeap & lh)
  {
    if (calls >= refvalues.Size())
      return;
    const double refvalue = refvalues[calls];
    double variable = pde.lock()->GetVariable(varname, false);

    if (cdash)
    {
      string dart_varname = varname;
      for (unsigned int i=0; i< dart_varname.length(); i++) {
        char & c = dart_varname[i];
        if(c==' ' || c==':' || c=='-' || c=='.') {
          dart_varname.erase(i,1);
          i--;
        }
      }
      cout << "<DartMeasurement name=" << '"' << dart_varname << '"' << endl;
      cout << "type=\"numeric/double\">" << variable << "</DartMeasurement>" << endl;
    }


    if(abstol)
    {
      if (abs(variable - refvalue) > tolerance)
      {
        ostringstream exctext;
        exctext << "NumProcTestVariable(" << GetName();
        exctext << "NumProcTestVariable(" << GetName();
        exctext << ": Violated absolute tolerance: ";
        exctext << "value = " << variable;
        exctext << ", refvalue = " << refvalue;
        exctext << ", tolerance = " << tolerance;
        throw Exception (exctext.str()); 
      }
    }
    else 
    {
      if (abs(variable - refvalue)/abs(refvalue) > tolerance)
      {
        ostringstream exctext;
        exctext << "NumProcTestVariable(" << GetName();
        exctext << "NumProcTestVariable(" << GetName();
        exctext << ": Violated relative tolerance: ";
        exctext << "value = " << variable;
        exctext << ", refvalue = " << refvalue;
        exctext << ", tolerance = " << tolerance;
        throw Exception (exctext.str()); 
      }
    }
    std::cout << " variable " << varname << " withtin tolerance: " << std::endl;
    std::cout << " value = " << variable << ", refvalue = " << refvalue << std::endl;
    std::cout << " abs. error. = " << abs(variable - refvalue)  << std::endl;
    std::cout << " rel. error. = " << abs(variable - refvalue) / refvalue  << std::endl;
    calls++;
  }


  virtual string GetClassName () const
  {
    return "NumProcTestVariable";
  }

  virtual void PrintReport (ostream & ost) const
  {
    ost << GetClassName() << endl
        << "Compare variable" << varname << " with reference values " 
        << refvalues << "and (";
    if (abstol) ost << "absolute)"; else ost << "relative)";
    ost << " tolerance of " << tolerance << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost <<
      "\n\nNumproc NumProcTestVariable:\n" \
      "-------------------------------\n" \
      "Compare variables with reference\n" \
      "value and throw exc. on failure\n\n"\
      "flags:\n"\
      "-variable=<name>\n"\
      "   variables (for diff. levels) to compare\n"\
      "-refvalues=[<val>,<val>,..]\n"\
      "   variable to compare (only on first level)\n"\
      "-refvalue=<val>\n"\
      "   reference value\n"\
      "-tolerance=<val>\n"\
      "   tolerance\n"\
      "-abstol\n"\
      "   use absolute tolerance \n"\
      "   (default is relative tol.)\n";
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










  // standard numprocs:

  static RegisterNumProc<NumProcSetValues> npinitsetvlues("setvalues");
  static RegisterNumProc<NumProcCalcFlux> npinitcalcflux("calcflux");
  static RegisterNumProc<NumProcVisualization> npinitvisual("visualization");
  static RegisterNumProc<NumProcIntegrate> npinitintegrate("integrate");
  static RegisterNumProc<NumProcWriteFile> npinitwf ("writefile");
  static RegisterNumProc<NumProcDrawFlux> npinitdf ("drawflux");
  static RegisterNumProc<NumProcDrawCoefficient> npinitdc ("draw");
  static RegisterNumProc<NumProcPause> npinitpause ("pause");
  static RegisterNumProc<NumProcTestVariable> nptestvar ("testvariable");

  static RegisterNumProc<NumProcLoadSolution2> npload ("loadgridfunction2");
  static RegisterNumProc<NumProcSaveSolution2> npsave ("savegridfunction2");
  static RegisterNumProc<NumProcAssembleLinearization> npassnonl ("assemblelinearization");



  static RegisterNumProc<NumProcEvaluate> npeval("evaluate");
  static RegisterNumProc<NumProcAnalyze> npanalyze ("analyze");
  static RegisterNumProc<NumProcWarn> npwarn("warn");
  static RegisterNumProc<NumProcTclTable> nptcltable("tcltable");
  static RegisterNumProc<NumProcTclMenu> nptclmenu("tclmenu");
  static RegisterNumProc<NumProcLoadSolution> nploadsol("loadsolution");
  
  static RegisterNumProc<NumProcSaveSolution> npsavesol ("savesolution");
  static RegisterNumProc<NumProcQuit> npquad ("quit");
  static RegisterNumProc<NumProcGenerateOne> npgen1 ("generateone");
  static RegisterNumProc<NumProcClearGridFunctions> npcleargf ("cleargridfunctions");
  
#ifdef ASTRID
      GetNumProcs().AddNumProc ("loadzipsolution", NumProcLoadZipSolution::Create, NumProcLoadZipSolution::PrintDoc);
      GetNumProcs().AddNumProc ("savezipsolution", NumProcSaveZipSolution::Create, NumProcSaveZipSolution::PrintDoc);
#endif
  
#ifdef VTUNE
  static RegisterNumProc<NumProcVtuneProfiling> npvtune ("vtune");
#endif
}
  




#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"

using namespace ngsolve;
void ExportDrawFlux(py::module &m)
{
  // cout << "exporting CalcFlux and DrawFlux numproc" << endl;

  m.def ("CalcFlux", FunctionPointer
           ([](shared_ptr<PDE> pde,
               shared_ptr<BilinearForm> bfa,
               shared_ptr<GridFunction> gfu,
               shared_ptr<GridFunction> gfflux,
               bool applyd) -> shared_ptr<NumProc>

            {
              return make_shared<NumProcCalcFlux> (pde, bfa, gfu, gfflux, applyd);
            }),
            py::arg("pde"), py::arg("bf"), py::arg("gf"),
            py::arg("flux"), py::arg("applyd")=false
	   );
  m.def ("DrawFlux", FunctionPointer
           ([](shared_ptr<BilinearForm> bfa,
               shared_ptr<GridFunction> gfu,
               const string & alabel,
               bool applyd,
               bool useall) -> shared_ptr<NumProc>
            
            {
              return make_shared<NumProcDrawFlux> (bfa, gfu, alabel, applyd, useall);
            }),
            py::arg("bf"), py::arg("gf"), 
            py::arg("label")="flux", py::arg("applyd")=false, py::arg("useall")=false
	   );


}
#endif
