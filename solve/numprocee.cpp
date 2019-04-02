#include "../include/solve.hpp"

namespace ngsolve
{

  ///
  class NumProcZZErrorEstimator : public NumProc
  {
  private:
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<GridFunction> gfu;
    ///
    shared_ptr<GridFunction> gferr;
    ///
    string filename; 
    ///
    ofstream outfile;

  public:
    NumProcZZErrorEstimator (shared_ptr<PDE> apde, const Flags & flags);
    virtual ~NumProcZZErrorEstimator() { ; }

    static void PrintDoc (ostream & ost);
    virtual void Do(LocalHeap & lh);
    virtual void PrintReport (ostream & ost) const;

    virtual string GetClassName () const
    {
      return "ZZ Error Estimator";
    }
  };


  NumProcZZErrorEstimator :: NumProcZZErrorEstimator (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde, flags)
  {
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("solution", ""));
    gferr = apde->GetGridFunction (flags.GetStringFlag ("error", ""));
    filename = flags.GetStringFlag ("filename","error.out");
    outfile.open (filename.c_str());
    apde->AddVariable (string("ZZerrest.")+GetName()+".err", 1e99);
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
      "    piece-wise constant gridfuntion to store the computed element-wise error\n"
      "-filename=<name>\n"
      "    file to store level, unknowns, estimated error\n";
  }
  
  

  void NumProcZZErrorEstimator :: Do(LocalHeap & lh)
  {
    cout << "ZZ error-estimator" << endl;

    if (bfa->NumIntegrators() == 0)
      throw Exception ("ZZErrorEstimator: Bilinearform needs an integrator");

    auto bfi = bfa->GetIntegrator(0);

    Flags fesflags;
     
    if(bfa->GetFESpace()->VarOrder())
      {
	cout << " Set Flux Space Relorder "  << bfa->GetFESpace()->GetRelOrder() << endl; 
	fesflags.SetFlag("relorder",bfa->GetFESpace()->GetRelOrder()); 
      }
    else 
      {
	int order = bfa->GetFESpace()->GetOrder();
	if(order == 0)
	  order = 1;
	cout << "Set Flux Space order " << order << endl;
	fesflags.SetFlag("order", order);
      } 
    
    fesflags.SetFlag ("dim", bfi->DimFlux());
    if (bfa->GetFESpace()->IsComplex())
      fesflags.SetFlag ("complex");
    

    *testout << " ************ ZZ ErrorEstimator fesflux " << endl; 

    shared_ptr<FESpace> fesflux = make_shared<H1HighOrderFESpace> (ma, fesflags);
    fesflux -> Update(lh);

    shared_ptr<GridFunction> flux = CreateGridFunction (fesflux, "fluxzz", 
                                                        Flags().SetFlag("novisual"));
    flux->Update();

    FlatVector<double> err = gferr->GetVector().FV<double>();
    err = 0;
  
    int ndom = ma->GetNDomains();
    for (int k = 0; k < ndom; k++)
      {
	CalcFluxProject (*gfu, *flux, bfi, 1, k, lh);
	CalcError (*gfu, *flux, bfi, err, k, lh);
      }

    // delete flux;
    // delete fesflux;

    double sum = 0;
    for (int i = 0; i < err.Size(); i++) 
      sum += err(i);

    cout << " estimated error = " << sqrt (sum) << endl;
    shared_ptr<PDE>(pde)->AddVariable (string("ZZerrest.")+GetName()+".err", sqrt(sum));
    

    outfile << ma->GetNLevels() 
	    << "  "<< bfa->GetFESpace()->GetNDof() 
	    << " " << sqrt(sum) << endl;
  }

  void NumProcZZErrorEstimator :: PrintReport (ostream & ost) const
  {
    ost << "NumProcZZErrorEstimator:" << endl;
    ost << "Bilinear-form = " << endl;
  }














  class NumProcDifference : public NumProc
  {
  private:
    /// use flux from bfa1
    shared_ptr<BilinearForm> bfa1;
    /// first gridfunction
    shared_ptr<GridFunction> gfu1;
    /// use flux from bfa2
    shared_ptr<BilinearForm> bfa2;
    /// second gridfunction
    shared_ptr<GridFunction> gfu2;
    /// use coefficient function ( coef_real, i * coef_imag )
    shared_ptr<CoefficientFunction> coef_real;
    /// imaginary part of function
    shared_ptr<CoefficientFunction> coef_imag;

    /// difference function
    shared_ptr<GridFunction> gfdiff;
    /// output to file
    string filename; 
    ofstream * file;

  public:
    NumProcDifference (shared_ptr<PDE> apde, const Flags & flags);
    virtual ~NumProcDifference();


    static void PrintDoc (ostream & ost);
    virtual void Do(LocalHeap & lh);
    virtual void PrintReport (ostream & ost) const;

    virtual string GetClassName () const
    {
      return "Calc Difference";
    }
  };



  NumProcDifference :: NumProcDifference (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde), bfa2(0), gfu2(0), coef_real(0), coef_imag(0)
  {
  
    bfa1 = apde->GetBilinearForm 
      (flags.GetStringFlag ("bilinearform1", flags.GetStringFlag ("bilinearform", "")));
    gfu1 = apde->GetGridFunction 
      (flags.GetStringFlag ("solution1", flags.GetStringFlag("solution","")));

    if ( flags.StringFlagDefined("bilinearform2") )
      {
	bfa2 = apde->GetBilinearForm 
	  (flags.GetStringFlag ("bilinearform2", flags.GetStringFlag ("bilinearform", "")));
	gfu2 = apde->GetGridFunction (flags.GetStringFlag ("solution2", ""));
      }
    else
      {
	coef_real = apde->GetCoefficientFunction(flags.GetStringFlag("function",""));
	if ( flags.StringFlagDefined("function_imag"))
	  coef_imag = apde->GetCoefficientFunction(flags.GetStringFlag("function_imag",""));
      }

    gfdiff = apde->GetGridFunction (flags.GetStringFlag ("diff", ""), 1);

    filename = flags.GetStringFlag ("filename","");
    if (filename.length() && ma->GetCommunicator().Rank() == 0)
      {
        if (!flags.GetDefineFlag ("append"))
          file = new ofstream (filename.c_str());
        else
          file = new ofstream (filename.c_str(), ios_base::app);
      }
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
    cout << IM(3) << "Compute difference ... " << flush;

    double sum = 0;
    // if (working_proc)
      {
	if (bfa1->NumIntegrators() == 0)
	  throw Exception ("Difference: Bilinearform1 needs an integrator");

	auto bfi1 = bfa1->GetIntegrator(0);
	FlatVector<double> diff = gfdiff->GetVector().FV<double> ();
	diff = 0;

	int ndom = ma->GetNDomains();

	if ( bfa2 )
	  {
	    auto bfi2 = bfa2->GetIntegrator(0);
    
	    for (int k = 0; k < ndom; k++)
	      {
		if (!bfa1->GetFESpace()->IsComplex())
		  {
		    CalcDifference (dynamic_cast<const S_GridFunction<double>&> (*gfu1), 
                                    dynamic_cast<const S_GridFunction<double>&> (*gfu2), 
                                    bfi1, bfi2, diff, k, lh);
		  }
		else
		  {
		    CalcDifference (dynamic_cast<const S_GridFunction<Complex>&> (*gfu1), 
                                    dynamic_cast<const S_GridFunction<Complex>&> (*gfu2), 
                                    bfi1, bfi2, diff, k, lh);
		  }	    
	      }
	  }
	else
	  {
	    for (int k = 0; k < ndom; k++)
	      CalcDifference (*gfu1, bfi1, coef_real, diff, k, lh);
	  }

	for (int i = 0; i < diff.Size(); i++)
	  sum += diff(i);
      }
    
      sum = ma->GetCommunicator().AllReduce (sum, MPI_SUM);

    cout << IM(1) << " total difference = " << sqrt (sum) << endl;
    shared_ptr<PDE>(pde)->AddVariable (string("calcdiff.")+GetName()+".diff", sqrt(sum), 6);
    
    int ndof = bfa1 -> GetFESpace()->GetNDofGlobal();

    if (file)
      {
	(*file) << ma->GetNLevels() 
		<< "  " << ndof
		<< "  " << sqrt(double (ndof))
		<< " " << sqrt(sum) << endl;
      }
  }



  void NumProcDifference :: PrintReport (ostream & ost) const
  {
    ost << "NumProcDifference:" << endl;
    ost << "Bilinear-form = " << endl;
  }















  ///
  class NumProcRTZZErrorEstimator : public NumProc
  {
  private:
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<GridFunction> gfu;
    ///
    shared_ptr<GridFunction> gferr;
  public:
    NumProcRTZZErrorEstimator (shared_ptr<PDE> apde, const Flags & flags);
    virtual ~NumProcRTZZErrorEstimator();
  
    virtual void Do(LocalHeap & lh);
    virtual void PrintReport (ostream & ost) const;

    virtual string GetClassName () const
    {
      return "RTZZ Error Estimator";
    }
  };



  NumProcRTZZErrorEstimator :: NumProcRTZZErrorEstimator (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("solution", ""));
    gferr = apde->GetGridFunction (flags.GetStringFlag ("error", ""));
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

    auto bfi = bfa->GetIntegrator(0);

    Flags fesflags;
    fesflags.SetFlag ("order", bfa->GetFESpace()->GetOrder());
    if (bfa->GetFESpace()->IsComplex())
      fesflags.SetFlag ("complex");

    auto fesflux = make_shared<HDivHighOrderFESpace> (ma, fesflags);
    fesflux -> Update(lh);

    // Flags flags;
    // flags.SetFlag ("novisual");
    auto flux = CreateGridFunction (fesflux, "fluxzz",
                                    Flags().SetFlag("novisual"));
    flux->Update();

    FlatVector<double> err = gferr->GetVector().FV<double>();

    err = 0;

    CalcFluxProject (*gfu, *flux, bfi, 1, -1, lh);
    CalcError (*gfu, *flux, bfi, err, -1, lh);

    double sum = 0;
    for (int i = 0; i < err.Size(); i++)
      sum += err(i);
    cout << "estimated error = " << sqrt (sum) << endl;

    shared_ptr<PDE>(pde)->AddVariable (string("RTZZerrest.")+GetName()+".err", sqrt(sum));

    static ofstream errout ("error.out");
    errout << ma->GetNLevels() 
	   << "  " << bfa->GetFESpace()->GetNDof() 
	   << " " << sqrt(sum) << endl;
  }

  void NumProcRTZZErrorEstimator :: PrintReport (ostream & ost) const
  {
    ost << "NumProcRTZZErrorEstimator:" << endl;
    ost << "Bilinear-form = " << endl;
  }






  
  ///
  class NumProcHierarchicalErrorEstimator : public NumProc
  {
  private:
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<BilinearForm> bfa2;
    ///
    shared_ptr<LinearForm> lff;
    ///
    shared_ptr<GridFunction> gfu;
    ///
    shared_ptr<GridFunction> gferr;
    ///
    shared_ptr<FESpace> vtest;
  public:
    NumProcHierarchicalErrorEstimator (shared_ptr<PDE> apde, const Flags & flags)
      : NumProc (apde)
    {
      bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
      bfa2 = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform2", ""), 1);
      if (!bfa2) bfa2 = bfa;
      lff = apde->GetLinearForm (flags.GetStringFlag ("linearform", ""));
      gfu = apde->GetGridFunction (flags.GetStringFlag ("solution", ""));
      vtest = apde->GetFESpace (flags.GetStringFlag ("testfespace", ""));
      gferr = apde->GetGridFunction (flags.GetStringFlag ("error", ""));
    }

    virtual ~NumProcHierarchicalErrorEstimator()
    {
      ;
    }
  
    /*
    static NumProc * Create (PDE & pde, const Flags & flags)
    {
      return new NumProcHierarchicalErrorEstimator (pde, flags);
    }
    */

    virtual void Do(LocalHeap & lh)
    {
      cout << "Hierarchical error-estimator" << endl;
      
      FlatVector<double> err = gferr->GetVector().FVDouble();
      if (!bfa->GetFESpace()->IsComplex())
	{
	  CalcErrorHierarchical (dynamic_cast<const S_BilinearForm<double>&> (*bfa), 
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


    virtual void PrintReport (ostream & ost) const
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
    shared_ptr<GridFunction> gferr;
    ///
    shared_ptr<GridFunction> gferr2;
    ///
    int minlevel;
    ///
    double fac;
    ///
    double factor;

  public:
    ///
    NumProcMarkElements (shared_ptr<PDE> apde, const Flags & flags);
    ///
    virtual ~NumProcMarkElements() { ; }
    ///
    virtual void Do(LocalHeap & lh);
    ///
    virtual string GetClassName () const
    {
      return "Element Marker";
    }
    virtual void PrintReport (ostream & ost) const;
  };







  NumProcMarkElements :: NumProcMarkElements (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    gferr = apde->GetGridFunction (flags.GetStringFlag ("error", ""));
    gferr2 = apde->GetGridFunction (flags.GetStringFlag ("error2", ""), 1);
    minlevel = int(flags.GetNumFlag ("minlevel", 0));
    fac = flags.GetNumFlag ("fac", -1);
    if (fac != -1)
      throw Exception ("numproc markelements:\n Flag 'fac' not supported anymore\nNew one is -factor=xxx, default 0.5, lower: less refinement");
    factor = flags.GetNumFlag ("factor", 0.5);
  }

  void NumProcMarkElements :: Do(LocalHeap & lh)
  {
    cout << "Element marker, " << flush;

    if (ma->GetNLevels() < minlevel) 
      {
	cout << endl;
	return;
      }

    FlatVector<double> err = gferr->GetVector().FV<double> ();

    double maxerr = 0;
    double toterr = 0;

    if (gferr2)
      {
	FlatVector<double> err2 = gferr2->GetVector().FV<double>();
      
	for (int i = 0; i < err.Size(); i++)
	  {
	    err(i) = sqrt (err(i) * err2(i));

	    if (err(i) > maxerr) maxerr = err(i);
	    toterr += err(i);
	  }

	cout << "goal driven error estimator, est.err. = " << toterr << endl;
      }
    else
      {
	for (int i = 0; i < err.Size(); i++)
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

	for (int i = 0; i < err.Size(); i++)
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

	if (markerr >= factor * toterr) break;
      }

    /*
    for(i = 0; i< ma->GetNSE(); i++)
      {
	Array<int> elts;
	ma->GetFaceElements(ma->GetSElFace(i),elts);
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

    if (ma->GetDimension() == 3)
      {
	int nse = ma->GetNSE();
	for (int i = 0; i < nse; i++)
	  Ng_SetSurfaceRefinementFlag (i+1, 0);
      }
	  
  }

  void NumProcMarkElements :: PrintReport (ostream & ost) const
  {
    ost << "NumProcMarkElements:" << endl;
  }








  class NumProcSetVisual : public NumProc
  {
    ///
    Flags visflags;
  public:
    ///
    NumProcSetVisual (shared_ptr<PDE> apde, const Flags & flags);
    ///
    virtual ~NumProcSetVisual ();
    
    virtual void Do(LocalHeap & lh);
  };






  NumProcSetVisual ::   
  NumProcSetVisual (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde), visflags (flags)
  {
    cout << "SetVisual has flags" << endl;
    visflags.PrintFlags(cout);
  }

  NumProcSetVisual ::  ~NumProcSetVisual ()
  {
    ;
  }


  void NumProcSetVisual :: Do(LocalHeap & lh)
  {
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
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<GridFunction> gfu;
    ///
    shared_ptr<GridFunction> gfflux;
    ///
    shared_ptr<GridFunction> gferr;
  public:
    NumProcPrimalDualErrorEstimator (shared_ptr<PDE> apde, const Flags & flags);
    virtual ~NumProcPrimalDualErrorEstimator() { ; }

    virtual void Do(LocalHeap & lh);
    virtual void PrintReport (ostream & ost) const;
    static void PrintDoc (ostream & ost);

    virtual string GetClassName () const
    {
      return "PrimalDual Error Estimator";
    }
  };



  NumProcPrimalDualErrorEstimator :: NumProcPrimalDualErrorEstimator (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearform", ""));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("solution", ""));
    gfflux = apde->GetGridFunction (flags.GetStringFlag ("flux", ""));
    gferr = apde->GetGridFunction (flags.GetStringFlag ("error", ""));
  }


  void NumProcPrimalDualErrorEstimator :: Do(LocalHeap & lh)
  {
    cout << "PrimalDual error-estimator" << endl;

    if (bfa->NumIntegrators() == 0)
      throw Exception ("PrimalDualErrorEstimator: Bilinearform needs an integrator");

    auto bfi = bfa->GetIntegrator(0);

    FlatVector<double> err = gferr->GetVector().FV<double>();
      // dynamic_cast<T_BaseVector<double>&> (gferr->GetVector()).FV();
    
    err = 0;
    CalcError (*gfu, *gfflux, bfi, err, -1, lh);
  
    double sum = 0;
    for (int i = 0; i < err.Size(); i++)
      sum += err(i);
    cout << "estimated error = " << sqrt (sum) << endl;
  }


  void NumProcPrimalDualErrorEstimator :: PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc PrimalDual - error estimator:\n" \
      "---------------------------\n" \
      "Computes the error comparing flux of a primal unknown and a dual unknown\n\n" \
      "Required flags:\n" \
      "-bilinearform=<bfname>\n"
      "    takes first integrator of bilinearform to compute the flux\n"
      "-solution=<solname>\n"
      "    gridfunction storing the finite element solution\n"
      "-flux=<fluxname>\n"
      "    gridfunction storing the flux solution\n"
      "-error=<errname>\n"
      "    piece-wise constant gridfuntion to store the computed element-wise error\n";
  }
  
  

  void NumProcPrimalDualErrorEstimator :: PrintReport (ostream & ost) const
  {
    ost << "NumProcPrimalDualErrorEstimator:" << endl;
    ost << "Bilinear-form = " << endl;
  }


  
  static RegisterNumProc<NumProcZZErrorEstimator> npinitzz("zzerrorestimator");
  static RegisterNumProc<NumProcRTZZErrorEstimator> npinitrtzz("rtzzerrorestimator");

  static RegisterNumProc<NumProcHierarchicalErrorEstimator> npinithieree("hierarchicalerrorestimator");
  static RegisterNumProc<NumProcPrimalDualErrorEstimator> npinitpdee("primaldualerrorestimator");  

  static RegisterNumProc<NumProcDifference> npinitdiff("difference");
  static RegisterNumProc<NumProcMarkElements> npinitmark("markelements");


  namespace numprocee_cpp
  {
    int link_it;
  }
}
