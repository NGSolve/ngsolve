#include "../include/solve.hpp"

namespace ngsolve
{
  using namespace ngsolve;
  using namespace ngcomp; 
  // using namespace ngbla;
  // using namespace ngla;

  ///
  class NumProcZZErrorEstimator : public NumProc
  {
  private:
    ///
    BilinearForm * bfa;
    ///
    GridFunction * gfu;
    ///
    GridFunction * gferr;
    string filename; 

  public:
    NumProcZZErrorEstimator (PDE & apde, const Flags & flags);
    virtual ~NumProcZZErrorEstimator();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcZZErrorEstimator (pde, flags);
    }

    static void PrintDoc (ostream & ost);
    virtual void Do(LocalHeap & lh);
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "ZZ Error Estimator";
    }
  };



  NumProcZZErrorEstimator :: NumProcZZErrorEstimator (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    gfu = pde.GetGridFunction (flags.GetStringFlag ("solution", ""));
    gferr = pde.GetGridFunction (flags.GetStringFlag ("error", ""));
    filename = flags.GetStringFlag ("filename","");
  }

  NumProcZZErrorEstimator :: ~NumProcZZErrorEstimator()
  {
    ;
  }
  
  void NumProcZZErrorEstimator :: PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc ZZ-error estimator:\n" \
      "---------------------------\n" \
      "Computes the Zienkiewicz-Zhu error estimator\n\n" \
      "Required flags:\n" \
      "-bilinearform=<bfname>\n"
      "    takes first integrator of bilinearform to compute the flux\n"
      "-solution=<solname>\n"
      "    gridfunction storing the finite element solution\n"
      "-error=<errname>\n"
      "    piece-wise constant gridfuntion to store the computed element-wise error\n";
  }
  
  

  void NumProcZZErrorEstimator :: Do(LocalHeap & lh)
  {
    cout << "ZZ error-estimator" << endl;

    if (bfa->NumIntegrators() == 0)
      throw Exception ("ZZErrorEstimator: Bilinearform needs an integrator");

    BilinearFormIntegrator * bfi = bfa->GetIntegrator(0);

    Flags fesflags;
  
     
    if(bfa->GetFESpace().VarOrder())
      {
	cout << " Set Flux Space Relorder "  << bfa->GetFESpace().GetRelOrder() << endl; 
	fesflags.SetFlag("relorder",bfa->GetFESpace().GetRelOrder()); 
      }
    else 
      {
	int order = bfa->GetFESpace().GetOrder();
	if(order == 0)
	  order = 1;
	cout << "Set Flux Space order " << order << endl;
	fesflags.SetFlag("order", order);
      } 

    
    fesflags.SetFlag ("dim", bfi->DimFlux());
    if (bfa->GetFESpace().IsComplex())
      fesflags.SetFlag ("complex");
    

    *testout << " ************ ZZ ErrorEstimator fesflux " << endl; 
    FESpace * fesflux = 
      new H1HighOrderFESpace (ma, fesflags);
    // new NodalFESpace (ma, fesflags);   // for low order fe
 
    
    
    fesflux -> Update(lh);

    Flags flags;
    //    flags.SetFlag ("novisual");
    GridFunction * flux = CreateGridFunction (fesflux, "fluxzz", flags);
    flux->Update();

    FlatVector<double> err = 
      dynamic_cast<T_BaseVector<double>&> (gferr->GetVector()).FV();

    err = 0;
  
    int ndom = ma.GetNDomains();
    //    for (int k = 4; k <= 4; k++)
    for (int k = 0; k < ndom; k++)
      {
	if (!bfa->GetFESpace().IsComplex())
	  {
	    CalcFluxProject (ma, 
			     dynamic_cast<const S_GridFunction<double>&> (*gfu), 
			     dynamic_cast<S_GridFunction<double>&> (*flux), 
			     *bfi,
			     1, k, lh);
	  
	    CalcError (ma, 
		       dynamic_cast<const S_GridFunction<double>&> (*gfu), 
		       dynamic_cast<const S_GridFunction<double>&> (*flux), 
		       *bfi,
		       err, k, lh);
	  }
	else
	  {
	    CalcFluxProject (ma, 
			     dynamic_cast<const S_GridFunction<Complex>&> (*gfu), 
			     dynamic_cast<S_GridFunction<Complex>&> (*flux), 
			     *bfi,
			     1, k, lh);
	  
	    CalcError (ma,
		       dynamic_cast<const S_GridFunction<Complex>&> (*gfu), 
		       dynamic_cast<const S_GridFunction<Complex>&> (*flux), 
		       *bfi,
		       err, k, lh);
	  }
      
      }
    // delete flux;
    double sum = 0;
    for (int i = 0; i < err.Size(); i++)
      sum += err(i);
    cout << " estimated error = " << sqrt (sum) << endl;
    static ofstream errout ("error.out");
    errout << ma.GetNLevels() 
	   << "  " << bfa->GetFESpace().GetNDof() 
	   << "  " << sqrt(double (bfa->GetFESpace().GetNDof())) 
	   << " " << sqrt(sum) << endl;

    static ofstream ofile (filename.c_str());
    ofile << ma.GetNLevels() 
	   << "  " << bfa->GetFESpace().GetNDof() 
	   << "  " << sqrt(double (bfa->GetFESpace().GetNDof())) 
	   << " " << sqrt(sum) << endl;
       
  }




  void NumProcZZErrorEstimator :: PrintReport (ostream & ost)
  {
    ost << "NumProcZZErrorEstimator:" << endl;
    ost << "Bilinear-form = " << endl;
  }














  class NumProcDifference : public NumProc
  {
  private:
    /// use flux from bfa1
    BilinearForm * bfa1;
    /// first gridfunction
    GridFunction * gfu1;
    /// use flux from bfa2
    BilinearForm * bfa2;
    /// second gridfunction
    GridFunction * gfu2;
    /// use coefficient function ( coef_real, i * coef_imag )
    CoefficientFunction * coef_real;
    /// imaginary part of function
    CoefficientFunction * coef_imag;

    /// difference function
    GridFunction * gfdiff;
    /// output to file
    string filename; 
    ofstream * file;

  public:
    NumProcDifference (PDE & apde, const Flags & flags);
    virtual ~NumProcDifference();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcDifference (pde, flags);
    }

    static void PrintDoc (ostream & ost);
    virtual void Do(LocalHeap & lh);
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "Calc Difference";
    }
  };



  NumProcDifference :: NumProcDifference (PDE & apde, const Flags & flags)
    : NumProc (apde), bfa2(0), gfu2(0), coef_real(0), coef_imag(0)
  {
  
    bfa1 = pde.GetBilinearForm 
      (flags.GetStringFlag ("bilinearform1", flags.GetStringFlag ("bilinearform", "")));
    gfu1 = pde.GetGridFunction 
      (flags.GetStringFlag ("solution1", flags.GetStringFlag("solution","")));

    if ( flags.StringFlagDefined("bilinearform2") )
      {
	bfa2 = pde.GetBilinearForm 
	  (flags.GetStringFlag ("bilinearform2", flags.GetStringFlag ("bilinearform", "")));
	gfu2 = pde.GetGridFunction (flags.GetStringFlag ("solution2", ""));
      }
    else
      {
	coef_real = pde.GetCoefficientFunction(flags.GetStringFlag("function",""));
	if ( flags.StringFlagDefined("function_imag"))
	  coef_imag = pde.GetCoefficientFunction(flags.GetStringFlag("function_imag",""));
      }

    gfdiff = pde.GetGridFunction (flags.GetStringFlag ("diff", ""), 1);

    filename = flags.GetStringFlag ("filename","");
    if (filename.length())
      file = new ofstream (filename.c_str());
    else
      file = 0;
  }

  NumProcDifference :: ~NumProcDifference()
  {
    delete file;
  }
  
  void NumProcDifference :: PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc Difference:\n" \
      "---------------------------\n" \
      "Computes the difference of fluxes of 2 gridfunctions\n"
			"or of flux of gridfunction and coefficient function\n\n" \
      "Required flags:\n" \
      "-bilinearform1=<bfname>\n"
      "    takes first integrator of bilinearform1 to compute the flux\n"
      "-solution1=<solname>\n"
      "    gridfunction storing the finite element solution 1\n"
      "-bilinearform2=<bfname>\n"
      "    takes first integrator of bilinearform2 to compute the flux for solution2\n"
      "-solution2=<solname>\n"
      "    gridfunction storing the finite element solution 2\n"
      "-function=<coeff_name>\n"
      "    coefficient storing the real part of the function\n"
      "-function_imag=<coeff_name>\n"
      "    coefficient storing the imaginary part of the function\n"
      "-diff=<errname>\n"
      "    piece-wise constant gridfuntion to store the computed element-wise error\n"
      "    between flux1 and flux2 or flux1 and function\n";
  }
  
  

  void NumProcDifference :: Do(LocalHeap & lh)
  {
    cout << "Compute difference" << endl;

    

    if (bfa1->NumIntegrators() == 0)
      throw Exception ("Difference: Bilinearform1 needs an integrator");

    BilinearFormIntegrator * bfi1 = bfa1->GetIntegrator(0);
    FlatVector<double> diff =
      dynamic_cast<T_BaseVector<double>&> (gfdiff->GetVector()).FV();
    diff = 0;

    int ndom = ma.GetNDomains();

    if ( bfa2 )
      {
	BilinearFormIntegrator * bfi2 = bfa2->GetIntegrator(0);
    
	for (int k = 0; k < ndom; k++)
	  {
	    if (!bfa1->GetFESpace().IsComplex())
	      {
		CalcDifference (ma, 
				dynamic_cast<const S_GridFunction<double>&> (*gfu1), 
				dynamic_cast<const S_GridFunction<double>&> (*gfu2), 
				*bfi1, *bfi2,
				diff, k, lh);
	      }
	    else
	      {
		CalcDifference (ma,
				dynamic_cast<const S_GridFunction<Complex>&> (*gfu1), 
				dynamic_cast<const S_GridFunction<Complex>&> (*gfu2), 
				*bfi1, *bfi2,
				diff, k, lh);
	      }
	    
	  }
      }
    else
      {
	for (int k = 0; k < ndom; k++)
	  {
	    if (!bfa1->GetFESpace().IsComplex())
	      {
		CalcDifference (ma, 
				dynamic_cast<const S_GridFunction<double>&> (*gfu1), 
				*bfi1,
				coef_real, coef_imag,
				diff, k, lh);
	      }
	    else
	      {
		CalcDifference (ma,
				dynamic_cast<const S_GridFunction<Complex>&> (*gfu1), 
				*bfi1, 
				coef_real, coef_imag,
				diff, k, lh);
	      }
	  }

      }
   
    double sum = 0;
    for (int i = 0; i < diff.Size(); i++)
      sum += diff(i);

    cout << " total difference = " << sqrt (sum) << endl;

    if (file)
      {
	(*file) << ma.GetNLevels() 
		<< "  " << bfa1->GetFESpace().GetNDof() 
		<< "  " << sqrt(double (bfa1->GetFESpace().GetNDof())) 
		<< " " << sqrt(sum) << endl;
      }
  }



  void NumProcDifference :: PrintReport (ostream & ost)
  {
    ost << "NumProcDifference:" << endl;
    ost << "Bilinear-form = " << endl;
  }















  ///
  class NumProcRTZZErrorEstimator : public NumProc
  {
  private:
    ///
    BilinearForm * bfa;
    ///
    GridFunction * gfu;
    ///
    GridFunction * gferr;
  public:
    NumProcRTZZErrorEstimator (PDE & apde, const Flags & flags);
    virtual ~NumProcRTZZErrorEstimator();
  
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcRTZZErrorEstimator (pde, flags);
    }

    virtual void Do(LocalHeap & lh);
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "RTZZ Error Estimator";
    }
  };



  NumProcRTZZErrorEstimator :: NumProcRTZZErrorEstimator (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    gfu = pde.GetGridFunction (flags.GetStringFlag ("solution", ""));
    gferr = pde.GetGridFunction (flags.GetStringFlag ("error", ""));
  }

  NumProcRTZZErrorEstimator :: ~NumProcRTZZErrorEstimator()
  {
    ;
  }

  void NumProcRTZZErrorEstimator :: Do(LocalHeap & lh)
  {
    cout << "RTZZ error-estimator" << endl;

    if (bfa->NumIntegrators() == 0)
      throw Exception ("RTZZErrorEstimator: Bilinearform needs an integrator");

    BilinearFormIntegrator * bfi = bfa->GetIntegrator(0);

    Flags fesflags;
    fesflags.SetFlag ("order", bfa->GetFESpace().GetOrder()+8);
    // fesflags.SetFlag ("dim", bfi->DimFlux());
    if (bfa->GetFESpace().IsComplex())
      fesflags.SetFlag ("complex");

    HDivHighOrderFESpace & fesflux = 
      *new HDivHighOrderFESpace (ma, fesflags);

    fesflux.Update(lh);

    Flags flags;
    GridFunction * flux = CreateGridFunction (&fesflux, "fluxzz", flags);
    flux->Update();

    FlatVector<double> err = 
      dynamic_cast<T_BaseVector<double>&> (gferr->GetVector()).FV();

    err = 0;
  
    int i, j, k;

    if (!bfa->GetFESpace().IsComplex())
      {
	CalcFluxProject (ma, 
			 dynamic_cast<const S_GridFunction<double>&> (*gfu), 
			 dynamic_cast<S_GridFunction<double>&> (*flux), 
			 *bfi,
			 1, -1, lh);
	  
	CalcError (ma, 
		   dynamic_cast<const S_GridFunction<double>&> (*gfu), 
		   dynamic_cast<const S_GridFunction<double>&> (*flux), 
		   *bfi,
		   err, -1, lh);
      }
    else
      {
	CalcFluxProject (ma, 
			 dynamic_cast<const S_GridFunction<Complex>&> (*gfu), 
			 dynamic_cast<S_GridFunction<Complex>&> (*flux), 
			 *bfi,
			 1, -1, lh);
	
	CalcError (ma,
		   dynamic_cast<const S_GridFunction<Complex>&> (*gfu), 
		   dynamic_cast<const S_GridFunction<Complex>&> (*flux), 
		   *bfi,
		   err, -1, lh);
      }

    // delete flux;
    double sum = 0;
    for (i = 0; i < err.Size(); i++)
      sum += err(i);
    cout << "estimated error = " << sqrt (sum) << endl;
    static ofstream errout ("error.out");
    errout << ma.GetNLevels() 
	   << "  " << bfa->GetFESpace().GetNDof() 
	   << "  " << sqrt(double (bfa->GetFESpace().GetNDof())) 
	   << " " << sqrt(sum) << endl;
  }




  void NumProcRTZZErrorEstimator :: PrintReport (ostream & ost)
  {
    ost << "NumProcRTZZErrorEstimator:" << endl;
    ost << "Bilinear-form = " << endl;
  }

  ///
  class NumProcHierarchicalErrorEstimator : public NumProc
  {
  private:
    ///
    BilinearForm * bfa;
    ///
    BilinearForm * bfa2;
    ///
    LinearForm * lff;
    ///
    GridFunction * gfu;
    ///
    GridFunction * gferr;
    ///
    FESpace * vtest;
  public:
    NumProcHierarchicalErrorEstimator (PDE & apde, const Flags & flags)
      : NumProc (apde)
    {
      bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
      bfa2 = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform2", ""), 1);
      if (!bfa2) bfa2 = bfa;
      lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", ""));
      gfu = pde.GetGridFunction (flags.GetStringFlag ("solution", ""));
      vtest = pde.GetFESpace (flags.GetStringFlag ("testfespace", ""));
      gferr = pde.GetGridFunction (flags.GetStringFlag ("error", ""));
    }

    virtual ~NumProcHierarchicalErrorEstimator()
    {
      ;
    }
  
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcHierarchicalErrorEstimator (pde, flags);
    }

    virtual void Do(LocalHeap & lh)
    {
      cout << "Hierarchical error-estimator" << endl;
      
      FlatVector<double> err = gferr->GetVector().FVDouble();
      if (!bfa->GetFESpace().IsComplex())
	{
	  CalcErrorHierarchical (ma, 
				 dynamic_cast<const S_BilinearForm<double>&> (*bfa), 
				 dynamic_cast<const S_BilinearForm<double>&> (*bfa2), 
				 dynamic_cast<const S_LinearForm<double>&> (*lff), 
				 dynamic_cast<S_GridFunction<double>&> (*gfu), 
				 *vtest, err, lh);
	}

      // delete flux;
      double sum = 0;
      for (int i = 0; i < err.Size(); i++)
	sum += err(i);
      cout << "estimated error = " << sqrt (sum) << endl;
    }


    virtual void PrintReport (ostream & ost)
    {
      ost << "NumProcHierarchicalErrorEstimator:" << endl;
      ost << "Bilinear-form = " << endl;
    }
    
    virtual string GetClassName () const
    {
      return "Hierarchical Error Estimator";
    }
  };












  /**
     Mark elements for refinement
  */
  class NumProcMarkElements : public NumProc
  {
  protected:
    ///
    GridFunction * gferr;
    ///
    GridFunction * gferr2;
    ///
    int minlevel;
    ///
    double fac;
    ///
    double factor;

  public:
    ///
    NumProcMarkElements (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcMarkElements();

    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcMarkElements (pde, flags);
    }
    ///
    virtual void Do();
    ///
    virtual string GetClassName () const
    {
      return "Element Marker";
    }
    virtual void PrintReport (ostream & ost);
  };







  NumProcMarkElements :: NumProcMarkElements (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    gferr = pde.GetGridFunction (flags.GetStringFlag ("error", ""));
    gferr2 = pde.GetGridFunction (flags.GetStringFlag ("error2", ""), 1);
    minlevel = int(flags.GetNumFlag ("minlevel", 0));
    fac = flags.GetNumFlag ("fac", -1);
    if (fac != -1)
      throw Exception ("numproc markelements:\n Flag 'fac' not supported anymore\nNew one is -factor=xxx, default 0.5, lower: less refinement");
    factor = flags.GetNumFlag ("factor", 0.5);
  }

  NumProcMarkElements :: ~NumProcMarkElements()
  {
    ;
  }

  void NumProcMarkElements :: Do()
  {
    cout << "Element marker, " << flush;

    if (ma.GetNLevels() < minlevel) 
      {
	cout << endl;
	return;
      }


    int i;
    FlatVector<double> err =
      dynamic_cast<T_BaseVector<double>&> (gferr->GetVector()).FV();

    double maxerr = 0;
    double toterr = 0;

    if (gferr2)
      {
	const FlatVector<double> & err2 =
	  dynamic_cast<T_BaseVector<double>&> (gferr2->GetVector()).FV();
      
	for (i = 0; i < err.Size(); i++)
	  {
	    err(i) = sqrt (err(i) * err2(i));

	    if (err(i) > maxerr) maxerr = err(i);
	    toterr += err(i);
	  }

	cout << "goal driven error estimator, est.err. = " << toterr << endl;
      }
    else
      {
	for (i = 0; i < err.Size(); i++)
	  {
	    toterr += err(i);
	    if (err(i) > maxerr) maxerr = err(i);
	  }
      }

    //cout << "maxerr = " << maxerr << "; toterr= " << toterr << endl;

    int nref;
    fac = 1;
  
    if (factor > 0.9999) factor = 0.9999;

    while (1)
      {
	fac *= 0.8;
	double markerr = 0;
	nref = 0;

	for (i = 0; i < err.Size(); i++)
	  {
	    if (err(i) > fac * maxerr)
	      {
		nref++;
		Ng_SetRefinementFlag (i+1, 1);
				  
		markerr += err(i);
	      }
	    else
	      {
		Ng_SetRefinementFlag (i+1, 0);
	      }
	  }


	//cout << "fac = " << fac << ", nmark = " << nref 
	//  << ", markerr = " << markerr << ", toterr = " << toterr << endl;


	if (markerr >= factor * toterr) break;
	//if(nref >= factor*err.Size()) break;
      }

    /*
    for(i = 0; i< ma.GetNSE(); i++)
      {
	Array<int> elts;
	ma.GetFaceElements(ma.GetSElFace(i),elts);
	for(int j=0; j<elts.Size(); j++)
	  if(refine[elts[j]])
	    Ng_SetSurfaceRefinementFlag(i+1,1);
      }
    */

	
    

    

    /*
    (*testout) << "toterr = " << toterr << " marked: " << endl;
    for (i = 0; i < err.Size(); i++)
      {
	if (err(i) > fac * maxerr)
	  (*testout) << "cell " << i << " err " << err(i) << endl;
      }
    */    

    cout << nref << "/" << err.Size() << " elements marked." << endl;

    if (ma.GetDimension() == 3)
      {
	int nse = ma.GetNSE();
	for (int i = 0; i < nse; i++)
	  Ng_SetSurfaceRefinementFlag (i+1, 0);
      }
	  
  }




  void NumProcMarkElements :: PrintReport (ostream & ost)
  {
    ost << "NumProcMarkElements:" << endl;
  }








  class NumProcSetVisual : public NumProc
  {
    ///
Flags visflags;
  public:
    ///
    NumProcSetVisual (PDE & apde, const Flags & flags);
    ///
    virtual ~NumProcSetVisual ();

    ///
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcSetVisual (pde, flags);
    }

    ///
    virtual void Do();
  };






  NumProcSetVisual ::   
  NumProcSetVisual (PDE & apde, const Flags & flags)
    : NumProc (apde), visflags (flags)
  {
    cout << "SetVisual has flags" << endl;
    visflags.PrintFlags(cout);
  }

  NumProcSetVisual ::  ~NumProcSetVisual ()
  {
    ;
  }


  void NumProcSetVisual :: Do()
  {
    int i;

    /*
      cout << "Set Visualization Flag:" << endl;

      for (i = 0; i < visflags.GetNStringFlags(); i++)
      {
      const char * name;
      const char * str;
      str = visflags.GetStringFlag (i, name);
      Ng_SetVisualizationParameter (name, str);
      }

      for (i = 0; i < visflags.GetNNumFlags(); i++)
      {
      const char * name;
      double val;
      char str[100];
      val = visflags.GetNumFlag (i, name);
      sprintf (str, "%f", val);
      cout << "set flag " << name << " to " << str << endl;
      Ng_SetVisualizationParameter (name, str);
      }
    */
  }







  ///
  class NumProcPrimalDualErrorEstimator : public NumProc
  {
  private:
    ///
    BilinearForm * bfa;
    ///
    GridFunction * gfu;
    ///
    GridFunction * gfflux;
    ///
    GridFunction * gferr;
  public:
    NumProcPrimalDualErrorEstimator (PDE & apde, const Flags & flags);
    virtual ~NumProcPrimalDualErrorEstimator();
  
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcPrimalDualErrorEstimator (pde, flags);
    }

    virtual void Do(LocalHeap & lh);
    virtual void PrintReport (ostream & ost);

    virtual string GetClassName () const
    {
      return "PrimalDual Error Estimator";
    }
  };



  NumProcPrimalDualErrorEstimator :: NumProcPrimalDualErrorEstimator (PDE & apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    gfu = pde.GetGridFunction (flags.GetStringFlag ("solution", ""));
    gfflux = pde.GetGridFunction (flags.GetStringFlag ("flux", ""));
    gferr = pde.GetGridFunction (flags.GetStringFlag ("error", ""));
  }

  NumProcPrimalDualErrorEstimator :: ~NumProcPrimalDualErrorEstimator()
  {
    ;
  }

  void NumProcPrimalDualErrorEstimator :: Do(LocalHeap & lh)
  {
    cout << "PrimalDual error-estimator" << endl;

    if (bfa->NumIntegrators() == 0)
      throw Exception ("PrimalDualErrorEstimator: Bilinearform needs an integrator");

    BilinearFormIntegrator * bfi = bfa->GetIntegrator(0);

    FlatVector<double> err = 
      dynamic_cast<T_BaseVector<double>&> (gferr->GetVector()).FV();

    err = 0;
  
    int i, j, k;
    if (!bfa->GetFESpace().IsComplex())
      {
	CalcError (ma, 
		   dynamic_cast<const S_GridFunction<double>&> (*gfu), 
		   dynamic_cast<const S_GridFunction<double>&> (*gfflux), 
		   *bfi,
		   err, -1, lh);
      }
    else
      {
	CalcError (ma,
		   dynamic_cast<const S_GridFunction<Complex>&> (*gfu), 
		   dynamic_cast<const S_GridFunction<Complex>&> (*gfflux), 
		   *bfi,
		   err, -1, lh);
      }
      
      
    double sum = 0;
    for (i = 0; i < err.Size(); i++)
      sum += err(i);
    cout << "estimated error = " << sqrt (sum) << endl;
  }




  void NumProcPrimalDualErrorEstimator :: PrintReport (ostream & ost)
  {
    ost << "NumProcPrimalDualErrorEstimator:" << endl;
    ost << "Bilinear-form = " << endl;
  }









  


  namespace
#ifdef MACOS
  numprocee_cpp
#endif
  {
    class Init
    { 
    public: 
      Init ();
    };
    
    Init::Init()
    {
      GetNumProcs().AddNumProc ("zzerrorestimator", NumProcZZErrorEstimator::Create, NumProcZZErrorEstimator::PrintDoc);
      GetNumProcs().AddNumProc ("difference", NumProcDifference::Create, NumProcDifference::PrintDoc);
//       GetNumProcs().AddNumProc ("difference_exact", NumProcDifferenceExact::Create, NumProcDifferenceExact::PrintDoc);
      GetNumProcs().AddNumProc ("rtzzerrorestimator", NumProcRTZZErrorEstimator::Create);
      GetNumProcs().AddNumProc ("hierarchicalerrorestimator", 
				NumProcHierarchicalErrorEstimator::Create);
      GetNumProcs().AddNumProc ("primaldualerrorestimator", 
				NumProcPrimalDualErrorEstimator::Create);
      GetNumProcs().AddNumProc ("markelements", NumProcMarkElements::Create);    }
    
    
    Init init;
    
  }
  


}
