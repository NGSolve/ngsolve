#include <solve.hpp>

#include <parallelngs.hpp>

namespace ngsolve
{
  using namespace ngsolve;

  using namespace ngparallel;

  EvalVariable :: EvalVariable(const MeshAccess & ama, const string & aname)
    : NGS_Object(ama,aname)
  {
    variable = NULL;
  }

  void EvalVariable :: SetVariable(double & avariable)
  {
    variable = &avariable;
  }
  
  double EvalVariable :: Evaluate(void)
  {
    if(variable)
      {
	*variable = evaluator.Eval((double*)(0));
	return *variable;
      }
    else
      return evaluator.Eval((double*)(0));
  }




  PDE :: PDE (MeshAccess & ama)
    : ma (ama)
  {
    levelsolved = -1;
    SetGood (true);
  }
  
  PDE :: ~PDE()
  {
    for (int i = 0; i < coefficients.Size(); i++)
      delete coefficients[i];
    coefficients.DeleteAll();
    for (int i = 0; i < spaces.Size(); i++)
      delete spaces[i];
    spaces.DeleteAll();
    for (int i = 0; i < gridfunctions.Size(); i++)
      delete gridfunctions[i];
    gridfunctions.DeleteAll();
    for (int i = 0; i < bilinearforms.Size(); i++)
      delete bilinearforms[i];
    bilinearforms.DeleteAll();
    for (int i = 0; i < linearforms.Size(); i++)
      delete linearforms[i];
    linearforms.DeleteAll();
    for (int i = 0; i < preconditioners.Size(); i++)
      delete preconditioners[i];
    preconditioners.DeleteAll();
    for (int i = 0; i < numprocs.Size(); i++)
      delete numprocs[i];
    numprocs.DeleteAll();

    for (int i = 0; i < string_constants.Size(); i++)
      delete string_constants[i];
    string_constants.DeleteAll();

    for(int i = 0; i < evaluators.Size(); i++)
      delete evaluators[i];
    evaluators.DeleteAll();

    for(int i = 0; i < CurvePointIntegratorFilenames.Size(); i++)
      delete CurvePointIntegratorFilenames[i];
    CurvePointIntegratorFilenames.DeleteAll();

    CurvePointIntegrators.DeleteAll();

    Ng_ClearSolutionData ();


  }

  void PDE :: SavePDE (const string & filename)
  {
    ;
  } 
  
  
  void PDE :: SetFilename(const string str)
  { 
    filename = str;
#ifdef WIN32
    for(int i=0; filename[i]!=0 && i<filename.size(); i++)
      if(filename[i] == '/')
	filename[i] = '\\';
#endif
  }

  void PDE :: SaveSolution (const string & filename, const bool ascii)
  {
    // ofstream outfile(filename.c_str());
    ofstream outfile;
    if(ascii)
      outfile.open(filename.c_str());
    else
      outfile.open(filename.c_str(), ios_base::binary);

    for (int i = 0; i < gridfunctions.Size(); i++)
      {
        if (printmessage_importance>0)
          cout << "Writing gridfunction " << gridfunctions.GetName(i) << endl;

	if(gridfunctions[i]->IsUpdated())
	  {
	    if(ascii)
	      gridfunctions[i]->GetVector().SaveText (outfile);
	    else
	      gridfunctions[i]->GetVector().Save (outfile);
	  }
	else
	  {
	    cerr << "ERROR: Gridfunction \"" << gridfunctions.GetName(i) << "\" is not initialized" << endl
		 << "       => not in the output file" << endl;
	  }
      }
    outfile.close();
  }

  ///
  void PDE :: LoadSolution (const string & filename, const bool ascii)
  {
#ifdef NETGEN_ELTRANS
    int geometryorder = 1;
    if (constants.Used ("geometryorder"))
      geometryorder = int (constants["geometryorder"]);

    bool rational = false;
    if (constants.Used ("rationalgeometry"))
      rational = bool (constants["rationalgeometry"]);

    Ng_HighOrder (geometryorder, rational);
#endif    

    ifstream infile;
    if(ascii)
      infile.open(filename.c_str());
    else
      infile.open(filename.c_str(), ios_base::binary);

    LocalHeap lh(100009);
    for (int i = 0; i < spaces.Size(); i++)
      spaces[i]->Update(lh);
    for (int i = 0; i < gridfunctions.Size(); i++)
      {
	gridfunctions[i]->Update();
        if (printmessage_importance>0)
          cout << "Loading gridfunction " << gridfunctions.GetName(i) << endl;
	if(ascii)
	  gridfunctions[i]->GetVector().LoadText (infile);
	else
	  gridfunctions[i]->GetVector().Load (infile);
      }
    infile.close();
  }


  void PDE :: PrintReport (ostream & ost)
  { 
    ost << endl << "PDE Description:" << endl;


    for (int i = 0; i < constants.Size(); i++)
      ost << "constant " << constants.GetName(i) << " = " << constants[i] << endl;
    for (int i = 0; i < string_constants.Size(); i++)
      ost << "string constant " << string_constants.GetName(i) << " = " << string_constants[i] << endl;
    for (int i = 0; i < variables.Size(); i++)
      ost << "variable " << variables.GetName(i) << " = " << variables[i] << endl;

    ost << endl;
  
    ost << "Coefficients:" << endl
	<< "-------------" << endl;
    for (int i = 0; i < coefficients.Size(); i++)
      {
	ost << "coefficient " << coefficients.GetName(i) << ":" << endl;
	coefficients[i]->PrintReport (ost);
      }
    
    ost << endl
	<< "Spaces:" << endl
	<< "-------" << endl;
    for (int i = 0; i < spaces.Size(); i++)
      {
	ost << "space " << spaces.GetName(i) << ":" << endl;
	spaces[i]->PrintReport (ost);
      }

    ost << endl
	<< "Bilinear-forms:" << endl
	<< "---------------" << endl;
    for (int i = 0; i < bilinearforms.Size(); i++)
      {
	ost << "bilinear-form " << bilinearforms.GetName(i) << ":" << endl;
	bilinearforms[i]->PrintReport (ost);
      }

    ost << endl 
	<< "Linear-forms:" << endl
	<< "-------------" << endl;
    for (int i = 0; i < linearforms.Size(); i++)
      {
	ost << "linear-form " << linearforms.GetName(i) << ":" << endl;
	linearforms[i]->PrintReport (ost);
      }

    ost << endl 
	<< "Grid-functions:" << endl
	<< "---------------" << endl;
    for (int i = 0; i < gridfunctions.Size(); i++)
      {
	ost << "grid-function " << gridfunctions.GetName(i) << ":" << endl;
	gridfunctions[i]->PrintReport (ost);
      }

    ost << endl
	<< "Preconditioners:" << endl
	<< "----------------" << endl;
    for (int i = 0; i < preconditioners.Size(); i++)
      {
	ost << "preconditioner " << preconditioners.GetName(i) << ":" << endl;
	preconditioners[i]->PrintReport (ost);
      }

    ost << endl 
	<< "Numprocs:" << endl
	<< "---------" << endl;
    for (int i = 0; i < numprocs.Size(); i++)
      {
	ost << "numproc " << numprocs.GetName(i) << ":" << endl;
	numprocs[i]->PrintReport (ost);
      }
  }
  

  void PDE :: PrintMemoryUsage (ostream & ost)
  {
    Array<MemoryUsageStruct*> memuse;
    for (int i = 0; i < spaces.Size(); i++)
      spaces[i]->MemoryUsage (memuse);
    for (int i = 0; i < bilinearforms.Size(); i++)
      bilinearforms[i]->MemoryUsage (memuse);
    for (int i = 0; i < linearforms.Size(); i++)
      linearforms[i]->MemoryUsage (memuse);
    for (int i = 0; i < gridfunctions.Size(); i++)
      gridfunctions[i]->MemoryUsage (memuse);
    for (int i = 0; i < preconditioners.Size(); i++)
      preconditioners[i]->MemoryUsage (memuse);

    int sumbytes = 0, sumblocks = 0;
    for (int i = 0; i < memuse.Size(); i++)
      {
	ost << memuse[i]->Name() << ": " << memuse[i]->NBytes()
	    << " bytes in " << memuse[i]->NBlocks() << " blocks." << endl;
	sumbytes += memuse[i]->NBytes();
	sumblocks += memuse[i]->NBlocks();
      }
    if (printmessage_importance>0)
      cout << "total bytes " << sumbytes << " in " << sumblocks << " blocks." << endl;
  }













  bool PDE ::
  ConstantUsed (const string & aname) const
  {
    return constants.Used(aname);
  }
  
  double PDE ::
  GetConstant (const string & name, bool opt) const
  { 
    if (constants.Used(name))
      return constants[name]; 
    if (opt) return 0;

    stringstream str;
    str << "Constant '" << name << "' not defined\n";
    throw Exception (str.str());
  }

  bool PDE ::
  StringConstantUsed (const string & aname) const
  {
    return string_constants.Used(aname);
  }
  
  string PDE ::
  GetStringConstant (const string & name, bool opt) const
  { 
    if (string_constants.Used(name))
      return *string_constants[name]; 
    if (opt) return string("");

    stringstream str;
    str << "String constant '" << name << "' not defined\n";
    throw Exception (str.str());
  }

  bool PDE ::
  VariableUsed (const string & aname) const
  {
    return variables.Used(aname);
  }
  
  double & PDE :: 
  GetVariable (const string & name, bool opt)
  { 
    if (variables.Used(name))
      return variables[name]; 

    static double dummy;
    if (opt) return dummy;

    stringstream str;
    str << "Variable '" << name << "' not defined\n";
    throw Exception (str.str());
  }

  CoefficientFunction * PDE :: 
  GetCoefficientFunction (const string & name, bool opt)
  { 
    if (coefficients.Used(name))
      return coefficients[name]; 

    if (opt) return 0;
    throw Exception (string ("CoefficientFunction '") + name + "' not defined\n");
  }

  FESpace * PDE :: 
  GetFESpace (const string & name, bool opt)
  { 
    if (spaces.Used(name))
      return spaces[name]; 

    if (opt) return 0;
    throw Exception (string("FESpace '") + name + "' not defined\n");
  }

  GridFunction * PDE :: 
  GetGridFunction (const string & name, bool opt)
  { 
    if (gridfunctions.Used(name))
      return gridfunctions[name]; 

    if (opt) return 0;
    throw Exception (string("GridFunction '") + name + "' not defined\n");
  }

  BilinearForm * PDE :: 
  GetBilinearForm (const string & name, bool opt)
  { 
    if (bilinearforms.Used(name))
      return bilinearforms[name]; 

    if (opt) return 0;
    throw Exception (string("Bilinear-form '") + name + "' not defined\n");
  }

  LinearForm * PDE :: 
  GetLinearForm (const string & name, bool opt)
  {
    if (linearforms.Used(name))
      return linearforms[name]; 

    if (opt) return 0;
    throw Exception (string("Linear-form '") + name + "' not defined\n");
  }

  Preconditioner * PDE :: 
  GetPreconditioner (const string & name, bool opt)
  { 
    if (preconditioners.Used(name))
      return preconditioners[name]; 

    if (opt) return 0;
    throw Exception (string("Preconditioner '") + name + "' not defined\n");
  }

  NumProc * PDE :: 
  GetNumProc (const string & name, bool opt)
  { 
    if (numprocs.Used(name))
      return numprocs[name]; 

    if (opt) return 0;
    throw Exception (string("Numproc '") + name + "' not defined\n");
  }

  const CoefficientFunction * PDE :: 
  GetCoefficientFunction (const string & name, bool opt) const
  { 
    if (coefficients.Used(name))
      return coefficients[name]; 

    if (opt) return 0;
    stringstream str;
    str << "CoefficientFunction '" << name << "' not defined\n";
    throw Exception (str.str());
  }

  const FESpace * PDE :: 
  GetFESpace (const string & name, bool opt) const
  { 
    if (spaces.Used(name))
      return spaces[name]; 

    if (opt) return 0;
    stringstream str;
    str << "FESpace '" << name << "' not defined\n";
    throw Exception (str.str());
  }

  const GridFunction * PDE :: 
  GetGridFunction (const string & name, bool opt) const
  { 
    if (gridfunctions.Used(name))
      return gridfunctions[name]; 

    if (opt) return 0;
    stringstream str;
    str << "Grid-function '" << name << "' not defined\n";
    throw Exception (str.str());
  }

  const BilinearForm * PDE :: 
  GetBilinearForm (const string & name, bool opt) const
  { 
    if (bilinearforms.Used(name))
      return bilinearforms[name]; 

    if (opt) return 0;
    cout << "name = (" << name << ")" << endl;
    stringstream str;
    str << "Bilinear-form '" << name << "' not defined\n";
    throw Exception (str.str());
  }

  const LinearForm * PDE :: 
  GetLinearForm (const string & name, bool opt) const
  { 
    if (linearforms.Used(name))
      return linearforms[name]; 

    if (opt) return 0;
    stringstream str;
    str << "Linear-form '" << name << "' not defined\n";
    throw Exception (str.str());
  }

  const Preconditioner * PDE :: 
  GetPreconditioner (const string & name, bool opt) const
  { 
    if (preconditioners.Used(name))
      return preconditioners[name]; 

    if (opt) return 0;
    stringstream str;
    str << "Preconditioner '" << name << "' not defined\n";
    throw Exception (str.str());
  }

  const NumProc * PDE :: 
  GetNumProc (const string & name, bool opt) const
  { 
    if (numprocs.Used(name))
      return numprocs[name]; 

    if (opt) return 0;
    stringstream str;
    str << "Numproc '" << name << "' not defined\n";
    throw Exception (str.str());
  }






  void PDE :: SolveBVP ()
  {
    static int timer = NgProfiler::CreateTimer ("Solver - Total");
    NgProfiler::RegionTimer reg (timer);

#ifdef PARALLEL
    MPI_Comm_size ( MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank ( MPI_COMM_WORLD, &id );
    hoprocs.SetSize(ntasks-1);
    for ( int i = 0; i<ntasks-1; i++ )
      hoprocs[i] = i+1;

    if ( id == 0 )
      for ( int dest = 1; dest < ntasks; dest++)
        MyMPI_Send("ngs_solvepde" , dest);
#endif

    size_t heapsize = 1000000;
    if (constants.Used ("heapsize"))
      heapsize = size_t(constants["heapsize"]);

#ifdef _OPENMP
    heapsize *= omp_get_max_threads();
#endif
    

    LocalHeap lh(heapsize);

    clock_t starttime, endtime;
    starttime = clock();

#ifdef PARALLEL
    double startwtime, endwtime;

    MPI_Barrier( MPI_COMM_WORLD );
    if ( id == 0 ) 
      startwtime = MPI_Wtime();
    MPI_Barrier( MPI_COMM_WORLD );
#endif

    bool solvebvpstd = true;


    if(solvebvpstd)
      {
	static int meshtimer = NgProfiler::CreateTimer ("Mesh adaption");
	NgProfiler::StartTimer (meshtimer);
    
	if (levelsolved >= 0)
	  {
	    if (constants.Used ("refinehp"))
	      Ng_Refine(NG_REFINE_HP);
	    else if (constants.Used ("refinep"))
	      Ng_Refine(NG_REFINE_P);
	    else
	      Ng_Refine(NG_REFINE_H);
	  }

	if (constants.Used ("secondorder"))
	  throw Exception ("secondorder is obsolete \n Please use  'define constant geometryorder = 2' instead");
	
	if (constants.Used ("hpref") && levelsolved == -1)
	  {
	    if( constants.Used("hpref_setnoorders") )
	      {
		if( constants.Used("hpref_geom_factor")) 
		  Ng_HPRefinement (int (constants["hpref"]),constants["hpref_geom_factor"],false);
		else
		  Ng_HPRefinement (int (constants["hpref"]),0.125 ,false);
	      }
	    else
	      {
		if( constants.Used("hpref_geom_factor")) 
		  Ng_HPRefinement (int (constants["hpref"]),constants["hpref_geom_factor"]);
		else
		  Ng_HPRefinement (int (constants["hpref"]));
	      }
	  }


	int geometryorder = 1;
	if (constants.Used ("geometryorder"))
	  geometryorder = int (constants["geometryorder"]);

        bool rational = false;
        if (constants.Used ("rationalgeometry"))
          rational = bool (constants["rationalgeometry"]);


	if (constants.Used("common_integration_order"))
	  {
            if (printmessage_importance>0)
              cout << " !!! comminintegrationorder = " << int (constants["common_integration_order"]) << endl;
	    Integrator::SetCommonIntegrationOrder (int (constants["common_integration_order"]));
	  }

	if (geometryorder)
          Ng_HighOrder (geometryorder, rational);


	NgProfiler::StopTimer (meshtimer);

        ma.UpdateBuffers();   // update global mesh infos

        if (printmessage_importance>0)
	  {
	    cout << "Solve at level " << ma.GetNLevels()-1
		 << ", NE = " << ma.GetNE() 
		 << ", NP = " << ma.GetNP() << endl;
	  }
    
	// line-integrator curve points can only be built if
	// element-curving has been done
	for(int i=0; i<CurvePointIntegrators.Size(); i++)
	  BuildLineIntegratorCurvePoints(*CurvePointIntegratorFilenames[i],
					 ma,
					 *CurvePointIntegrators[i]);
		    


	for(int i=0; i<todo.Size(); i++)
	  {
	    EvalVariable * ev = dynamic_cast<EvalVariable *>(todo[i]);
	    FESpace * fes = dynamic_cast<FESpace *>(todo[i]);
	    GridFunction * gf = dynamic_cast<GridFunction *>(todo[i]);
	    BilinearForm * bf = dynamic_cast<BilinearForm *>(todo[i]);
	    LinearForm * lf = dynamic_cast<LinearForm *>(todo[i]);
	    Preconditioner * pre = dynamic_cast<Preconditioner *>(todo[i]);
	    NumProc * np = dynamic_cast<NumProc *>(todo[i]);

	    if (ev)
	      {
                if (printmessage_importance>0)
                  cout << "evaluate variable " << ev->GetName() << " = " << ev->Evaluate() << endl;
	      }

	    else if (fes)
	      {
		try
		  {
		    NgProfiler::RegionTimer timer(fes->GetTimer());
                    if (printmessage_importance>0)
		      {
			cout << "Update "
			     << fes -> GetClassName()
			     << " " << fes -> GetName () << flush;
		      }

		    fes -> Update(lh);
		    lh.CleanUp();

		    if (fes->GetDimension() == 1)
                      {
                        if (printmessage_importance>0)
                          cout << ", ndof = " << fes -> GetNDof() << endl;
                      }
                    else
                      if (printmessage_importance>0)
                        {
                          cout << ", ndof = " 
                               << fes -> GetDimension() << " x " 
                               << fes -> GetNDof() << endl;
                        }
		  }
		catch (exception & e)
		  {
		    throw Exception (e.what() + 
				     string ("\nthrown by update space ") +
				     string(fes->GetName()));
		  }
#ifdef _MSC_VER
# ifndef MSVC_EXPRESS
		catch (CException * e)
		  {
		    TCHAR msg[255];
		    e->GetErrorMessage(msg, 255);
		    throw Exception (msg + 
				     string ("\nthrown by update space ") +
				     string(fes->GetName()));
		  }
# endif // MSVC_EXPRESS
#endif
		catch (Exception & e)
		  {
		    e.Append (string ("\nthrown by update space ") +
			      string (fes->GetName()));
		    throw;
		  }
	      }

	    else if (gf)
	      {
		try
		  {
                    if (printmessage_importance>0)
                      cout << "Update gridfunction " << gf->GetName() << endl;
		    NgProfiler::RegionTimer timer(gf->GetTimer());
		    gf->Update();

		  }
		catch (exception & e)
		  {
		    throw Exception (e.what() + 
				     string ("\nthrown by update grid-function ") +
				     string (gf->GetName()));
		  }
#ifdef _MSC_VER
# ifndef MSVC_EXPRESS
		catch (CException * e)
		  {
		    TCHAR msg[255];
		    e->GetErrorMessage(msg, 255);
		    throw Exception (msg + 
				     string ("\nthrown by update grid-function ") +
				     string(gf->GetName()));
		  }
# endif // MSVC_EXPRESS
#endif
		catch (Exception & e)
		  {
		    e.Append (string ("\nthrown by update grid-function ") +
			      string (gf->GetName()));
		    throw;
		  }
	      }

	    else if (bf)
	      {
		try
		  {
                    if (printmessage_importance>0)
		      {
			cout << "update bilinear-form " << bf->GetName() << endl;
			(*testout) << "update bilinear-form " << bf->GetName() << endl;
		      }
		    NgProfiler::RegionTimer timer(bf->GetTimer());
		    bf->Assemble(lh);
		    lh.CleanUp();
		  }
      

		catch (exception & e)
		  {
		    throw Exception (e.what() + 
				     string ("\nthrown by update bilinear-form ") +
				     string (bf->GetName()));
		  }
#ifdef _MSC_VER
# ifndef MSVC_EXPRESS
		catch (CException * e)
		  {
		    TCHAR msg[255];
		    e->GetErrorMessage(msg, 255);
		    throw Exception (msg + 
				     string ("\nthrown by update bilinear-form ") +
				     string(bf->GetName()));
		  }
# endif // MSVC_EXPRESS
#endif
		catch (Exception & e)
		  {
		    e.Append (string ("\nthrown by update bilinear-form ") +
			      string (bf->GetName()));
		    throw;
		  }
	      }

	    else if (lf)
	      {
		if( lf->InitialAssembling() )
		  {
		    try
		      {
                        if (printmessage_importance>0)
                          cout << "Update linear-form " << lf->GetName() << endl;
			NgProfiler::RegionTimer timer(lf->GetTimer());
			
			lf->Assemble(lh);
			lh.CleanUp();
			
		      }
		    catch (exception & e)
		      {
			throw Exception (e.what() + 
					 string ("\nthrown by update linear-form ") +
					 string (lf->GetName()));
		      }
#ifdef _MSC_VER
# ifndef MSVC_EXPRESS
		    catch (CException * e)
		      {
			TCHAR msg[255];
			e->GetErrorMessage(msg, 255);
			throw Exception (msg + 
					 string ("\nthrown by update linear-form ") +
					 string(lf->GetName()));
		      }
# endif // MSVC_EXPRESS
#endif
		    catch (Exception & e)
		      {
			e.Append (string ("\nthrown by update linear-form ") +
				  string (lf->GetName()));
			throw;
		      }
		  }
	      }

	    else if (pre)
	      {
		try
		  {

		    if ( pre->LaterUpdate() )
		      {  
			if (printmessage_importance>0)
			  {  
			    cout << endl << "WARNING: Update of " << pre->ClassName() 
				 << "  " << pre->GetName() << " postponed!" << endl;
			  }
		      }
		    else
		      {	    
                        if (printmessage_importance>0)
			  {
			    cout << "Update " << pre->ClassName() 
				 << "  " << pre->GetName() << endl;
			  }
                        NgProfiler::RegionTimer timer(pre->GetTimer());
			pre->Update();
			//	  preconditioners[i]->Test();
		      }
		  }

		catch (exception & e)
		  {
		    throw Exception (e.what() + 
				     string ("\nthrown by update preconditioner ") +
				     string (pre->GetName()));
		  }
#ifdef _MSC_VER
# ifndef MSVC_EXPRESS
		catch (CException * e)
		  {
		    TCHAR msg[255];
		    e->GetErrorMessage(msg, 255);
		    throw Exception (msg + 
				     string ("\nthrown by update preconditioner ") +
				     string(pre->GetName()));
		  }
# endif // MSVC_EXPRESS
#endif
		catch (Exception & e)
		  {
		    e.Append (string ("\nthrown by update preconditioner ") +
			      string (pre->GetName()));
		    throw;
		  }
	      }

	    else if (np)
	      {
		try
		  {
                    if (printmessage_importance>0)
		      {
			cout << "Call numproc " << np->GetClassName() 
			     << "  " << np->GetName() << endl;
		      }
		    NgProfiler::RegionTimer timer(np->GetTimer());
		    np->Do(lh);
		    lh.CleanUp();
		  }
		catch (exception & e)
		  {
		    throw Exception (e.what() + 
				     string ("\nthrown by update numproc ") +
				     string (np->GetName()));
		  }
#ifdef _MSC_VER
# ifndef MSVC_EXPRESS
		catch (CException * e)
		  {
		    TCHAR msg[255];
		    e->GetErrorMessage(msg, 255);
		    throw Exception (msg + 
				     string ("\nthrown by update numproc ") +
				     string(np->GetName()));
		  }
# endif // MSVC_EXPRESS
#endif
		catch (Exception & e)
		  {
		    e.Append (string ("\nthrown by update numproc ") +
			      string (np->GetName()));
		    throw;
		  }
	      }
	    else
	      cerr << "???????????????" << endl;


	  }




	  for (int i = 0; i < preconditioners.Size(); i++)
	    if(!preconditioners[i]->SkipCleanUp())
              preconditioners[i]->CleanUpLevel();
          for (int i = 0; i < bilinearforms.Size(); i++)
	    if(!bilinearforms[i]->SkipCleanUp())
	      bilinearforms[i]->CleanUpLevel();
	  for (int i = 0; i < linearforms.Size(); i++)
	    if(!linearforms[i]->SkipCleanUp())
	      linearforms[i]->CleanUpLevel();
	  

	  // set solution data
	  for (int i = 0; i < gridfunctions.Size(); i++)
	    gridfunctions[i]->Visualize(gridfunctions.GetName(i));

	  Ng_Redraw();
	  levelsolved++;
      }
    
    if (printmessage_importance>0)
      cout << "Equation Solved" << endl;
    endtime = clock();
    if (printmessage_importance>0)
      cout << "Total Time = " << double(endtime - starttime)/CLOCKS_PER_SEC << endl << endl;
    
#ifdef PARALLEL
    MPI_Barrier( MPI_COMM_WORLD );
    if ( id == 0 ) {
      endwtime = MPI_Wtime();
      printf( "WALL CLOCK TIME = %f\n", double(endwtime-startwtime) );
    }
#endif
  }










  void PDE :: AddConstant (const string & name, double val)
  {
    if (printmessage_importance>0)
      cout << "add constant " << name << " = " << val << endl;
    constants.Set (name.c_str(), val);
  }

  void PDE :: AddStringConstant (const string & name, const string & val)
  {
    if (printmessage_importance>0)
      cout << "add string constant " << name << " = " << val << endl;
    if(string_constants.Used(name))
      delete string_constants[name];

    string_constants.Set (name.c_str(), new string(val));
  }

  void PDE :: AddVariable (const string & name, double val)
  {
    if (printmessage_importance>0)
      cout << "add variable " << name << " = " << val << endl;
    variables.Set (name.c_str(), val);
  }

  
  void PDE :: AddVariable (const string & name, EvalVariable * eval)
  {
    evaluators.Append(eval);
    todo.Append(eval);
    variables.Set (name, 0);
    eval->SetVariable(variables[name]);
    if (printmessage_importance>0)
      cout << "add variable " << name << " = " << eval->Evaluate() << endl;
  }

  

  void PDE :: AddCoefficientFunction (const string & name, CoefficientFunction* fun)
  {
    if (printmessage_importance>0)
      cout << "add coefficient-function, name = " << name << endl;
    coefficients.Set (name.c_str(), fun);
  }




  FESpace * PDE :: AddFESpace (const string & name, Flags & flags)
  {
    if (printmessage_importance>0)
      cout << "add fespace " << name << endl;

    FESpace * space = 0;
    if (flags.GetDefineFlag ("vec")) 
      flags.SetFlag ("dim", ma.GetDimension());
    if (flags.GetDefineFlag ("tensor")) 
      flags.SetFlag ("dim", sqr (ma.GetDimension()));
    if (flags.GetDefineFlag ("symtensor")) 
      flags.SetFlag ("dim", ma.GetDimension()*(ma.GetDimension()+1) / 2);

    int order = int(flags.GetNumFlag ("order", -1));
    string type = flags.GetStringFlag("type", "");
    
    if (type != "") // this should become standard
      {
	for (int i = 0; i < GetFESpaceClasses().GetFESpaces().Size(); i++)
	  {
	    if (type == GetFESpaceClasses().GetFESpaces()[i]->name)
	      {
		if ( id == 0 && ntasks > 1)
		  {
		    FESpace * hospace = GetFESpaceClasses().GetFESpaces()[i]-> creator ( ma, flags) ;
		    // low order space if existent
		    space = & (hospace -> LowOrderFESpace()) ;
		    // else space, but with  order 0
		    if ( space == 0 )
		      {
			flags.SetFlag("order",0.0);
			if ( hospace->IsComplex() )
			  flags.SetFlag("complex");
			
			space = GetFESpaceClasses().GetFESpaces()[i]->creator (ma, flags);
		      }
		  }
		else
		  space = GetFESpaceClasses().GetFESpaces()[i]->creator (ma, flags);
	      }
	  }
	
	if (type == "compound")
	  {
	    const Array<char*> & spacenames = flags.GetStringListFlag ("spaces");
	    Array<const FESpace*> cspaces (spacenames.Size());
	    for (int i = 0; i < cspaces.Size(); i++)
	      cspaces[i] = GetFESpace (spacenames[i]);
	    space = new CompoundFESpace (GetMeshAccess(), cspaces, flags);
	  }
	if (!space) 
	  cout << "warning: unknown fespace type " << type << endl;
      }
    else // old method
      {
	for (int i = 0; i < GetFESpaceClasses().GetFESpaces().Size(); i++)
	  {
	    if (flags.GetDefineFlag (GetFESpaceClasses().GetFESpaces()[i]->name))
#ifdef PARALLEL
	      {
		if ( id == 0 && ntasks > 1)
		  {
		    FESpace * hospace = GetFESpaceClasses().GetFESpaces()[i]-> creator ( ma, flags) ;
		    // low order space if existent
		    space = & (hospace -> LowOrderFESpace()) ;
		    // else space, but with  order 0
		    if ( space == 0 )
		      {
			flags.SetFlag("order",0.0);
			flags.PrintFlags(*testout);
			space = GetFESpaceClasses().GetFESpaces()[i]->creator (ma, flags);
		      }
		    //delete hospace;
		  }
		else
		  space = GetFESpaceClasses().GetFESpaces()[i]->creator (ma, flags);
	      }
#else
	    space = GetFESpaceClasses().GetFESpaces()[i]->creator (ma, flags);
#endif
	    
	  }
      }
    
    if (space) 
      ;
    else if (flags.GetDefineFlag ("compound"))
      {
	const Array<char*> & spacenames = flags.GetStringListFlag ("spaces");
        if (printmessage_importance>0)
          cout << "   spaces=" << spacenames << endl;
	Array<const FESpace*> cspaces (spacenames.Size());
	for (int i = 0; i < cspaces.Size(); i++)
	  cspaces[i] = GetFESpace (spacenames[i]);
	space = new CompoundFESpace (GetMeshAccess(), cspaces, flags);
      }
    else
      {
	if (order <= 1)
	  space = new NodalFESpace (GetMeshAccess(), flags);
	else
	  space = new H1HighOrderFESpace (GetMeshAccess(), flags);
      }
    
    // if (flags.GetDefineFlag ("bem")) space->SetBEM (1);

    if (flags.NumListFlagDefined ("dirichletboundaries"))
      {
	BitArray dirbnds(ma.GetNBoundaries());
	dirbnds.Clear();
	const Array<double> & array = flags.GetNumListFlag ("dirichletboundaries");
	for (int i = 0; i < array.Size(); i++)
	  dirbnds.Set (int(array[i])-1);
	space->SetDirichletBoundaries (dirbnds);
      }
  
    if (flags.NumListFlagDefined ("domains"))
      {
	BitArray definedon(ma.GetNDomains());
	definedon.Clear();
	const Array<double> & domains = flags.GetNumListFlag ("domains");
	for (int i = 0; i < domains.Size(); i++)
	  definedon.Set (int(domains[i])-1);
	space->SetDefinedOn (definedon);
      }

    if (flags.NumListFlagDefined ("boundaries"))
      {
	BitArray definedon(ma.GetNBoundaries());
	definedon.Clear();
	const Array<double> & boundaries = flags.GetNumListFlag ("boundaries");
	for (int i = 0; i < boundaries.Size(); i++)
	  definedon.Set (int(boundaries[i])-1);
	space->SetDefinedOnBoundary (definedon);
      }

    space->SetName (name);
    spaces.Set (name, space);
    todo.Append(space);

    return space;
  }

  void PDE :: AddFESpace (const string & name, FESpace * space)
  {
    space->SetName (name);
    spaces.Set (name, space);
    todo.Append(space);
  }




  GridFunction * PDE :: AddGridFunction (const string & name, Flags & flags)
  {
    if (printmessage_importance>0)
      cout << "add grid-function " << name << endl;

    string spacename = flags.GetStringFlag ("fespace", "");

    if (!spaces.Used (spacename))
      {
	throw Exception (string ("Gridfuncton '") + name +
			 "' uses undefined space '" + spacename + "'");
      }

    const FESpace * space = GetFESpace(spacename);
 
    GridFunction * gf = CreateGridFunction (space, name, flags);
    gridfunctions.Set (name, gf);

    todo.Append(gf);

    return gf;
  }

  void PDE :: AddGridFunction (const string & name, GridFunction * gf)
  {
    gf -> SetName (name);
    gridfunctions.Set (name, gf);
    todo.Append(gf);
  }


  BilinearForm * PDE :: AddBilinearForm (const string & name, Flags & flags)
  {
    if (printmessage_importance>0)
      cout << "add bilinear-form " << name << endl;
    string spacename = flags.GetStringFlag ("fespace", "");

    if (!spaces.Used (spacename))
      {
	cerr << "space " << spacename << " not defined " << endl;
	return 0;
      }
    const FESpace * space = spaces[spacename];

    const FESpace * space2 = NULL;
    if (flags.StringFlagDefined ("fespace2"))
      space2 = spaces[flags.GetStringFlag ("fespace2", "")];
  
    if (!space2)
      {
	bilinearforms.Set (name, CreateBilinearForm (space, name, flags));
      }

    //  else
    //    bilinearforms.Set (name, new T_BilinearForm<double > (*space, *space2, name));
    
    if (flags.StringFlagDefined ("linearform"))
      bilinearforms[name] -> SetLinearForm (GetLinearForm (flags.GetStringFlag ("linearform", 0)));

    todo.Append(bilinearforms[name]);

    return bilinearforms[name];
  }


 
  LinearForm * PDE :: AddLinearForm (const string & name, Flags & flags)
  {
    if (printmessage_importance>0)
      cout << "add linear-form " << name << endl;

    string spacename = flags.GetStringFlag ("fespace", "");

    if (!spaces.Used (spacename))
      {
	throw Exception (string ("Linear-form '") + name +
			 "' uses undefined space '" + spacename + "'");
      }

    const FESpace * space = spaces[spacename];

    linearforms.Set (name, CreateLinearForm (space, name, flags));
    todo.Append(linearforms[name]);

    return linearforms[name];
  }



  Preconditioner * PDE :: AddPreconditioner (const string & name, Flags & flags)
  {
    if (printmessage_importance>0)
      cout << "add preconditioner " << name << flush;

    //  flags.PrintFlags (cout);
    Preconditioner * pre = NULL;
    const char * type = flags.GetStringFlag ("type", NULL);

    if ( ntasks == 1 )
      {

	if (strcmp (type, "multigrid") == 0)
	  pre = new MGPreconditioner (this, flags, name);
	
	else if (strcmp (type, "direct") == 0)
	  pre = new DirectPreconditioner (this, flags, name);
	
	else if (strcmp (type, "local") == 0)
	  pre = new LocalPreconditioner (this, flags, name);
	
	//  else if (strcmp (type, "vefc") == 0)
	// pre = new VEFC_Preconditioner (this, flags, name);
	
	else if (strcmp (type, "twolevel") == 0)
	  pre = new TwoLevelPreconditioner (this, flags, name);
	
	else if (strcmp (type, "complex") == 0)
	  pre = new ComplexPreconditioner (this, flags, name);
	
	else if (strcmp (type, "chebychev") == 0)
	  pre = new ChebychevPreconditioner (this, flags, name);
	
	//  else if (strcmp (type, "constrained") == 0)
	//    pre = new ConstrainedPreconditioner (this, flags);
	
	else if (strcmp (type, "amg") == 0)
	  pre = new CommutingAMGPreconditioner (this, flags, name);
	
	else if (strcmp (type, "dndd") == 0)
	  pre = new DNDDPreconditioner (this, flags, name);
	
	// changed 08/19/2003, Bachinger
	else if (strcmp (type, "nonsymmetric") == 0)
	  pre = new NonsymmetricPreconditioner (this, flags, name);
        /*
	else if (strcmp (type, "wirebasket" ) == 0 )
	  {
	    const BilinearForm * bfa = 
	      GetBilinearForm(flags.GetStringFlag("bilinearform", "") );
	    if ( bfa -> GetFESpace() . IsComplex() )
	      pre = new WireBasketPreconditioner<Complex>(this, flags, name);
	    else
	      pre = new WireBasketPreconditioner<double>(this, flags, name);
	  }
        */
	else
	  for (int i = 0; i < GetPreconditionerClasses().GetPreconditioners().Size(); i++)
	    {
	      if (flags.GetDefineFlag (GetPreconditionerClasses().GetPreconditioners()[i]->name))
		pre = GetPreconditionerClasses().GetPreconditioners()[i]->creator (*this, flags);
	      
	      if (string(type) == GetPreconditionerClasses().GetPreconditioners()[i]->name)
		pre = GetPreconditionerClasses().GetPreconditioners()[i]->creator (*this, flags);
	    }
      }
    else  // parallel computations
      {
	if (strcmp (type, "multigrid") == 0)
	  {
	      pre = new MGPreconditioner (this, flags, name);
	  }
	else if (strcmp (type, "local") == 0)
	  {
	    if ( id >= 0 )
	      pre = new LocalPreconditioner (this, flags, name);
	  }
	
	else if (strcmp (type, "direct") == 0)
	  pre = new DirectPreconditioner (this, flags, name);
	
	else if (strcmp (type, "twolevel") == 0)
	  pre = new TwoLevelPreconditioner (this, flags, name);

	else if (strcmp (type, "complex") == 0)
	  pre = new ComplexPreconditioner (this, flags, name);
	
	else if (strcmp (type, "chebychev") == 0)
	  pre = new ChebychevPreconditioner (this, flags, name);
	
	//  else if (strcmp (type, "constrained") == 0)
	//    pre = new ConstrainedPreconditioner (this, flags);
	
	else if (strcmp (type, "amg") == 0)
	  pre = new CommutingAMGPreconditioner (this, flags, name);
	
	else if (strcmp (type, "dndd") == 0)
	  pre = new DNDDPreconditioner (this, flags, name);
	
	// changed 08/19/2003, Bachinger
	else if (strcmp (type, "nonsymmetric") == 0)
	  pre = new NonsymmetricPreconditioner (this, flags, name);
	else
	  for (int i = 0; i < GetPreconditionerClasses().GetPreconditioners().Size(); i++)
	    {
	      if (flags.GetDefineFlag (GetPreconditionerClasses().GetPreconditioners()[i]->name))
		pre = GetPreconditionerClasses().GetPreconditioners()[i]->creator (*this, flags);

	      if (string(type) == GetPreconditionerClasses().GetPreconditioners()[i]->name)
		pre = GetPreconditionerClasses().GetPreconditioners()[i]->creator (*this, flags);
	      
	      
	    }
      }

	
    if (!pre)
      throw Exception ("Unknown preconditioner type " + string(type));
    
    preconditioners.Set (name, pre);
    if (printmessage_importance>0)
      cout << ", type = " << pre->ClassName() << endl;
    
    if (!flags.GetDefineFlag ("ondemand"))
      todo.Append(pre);
    
    return pre;
  }



  void PDE :: AddNumProc (const string & name, NumProc * np)
  {
    if (printmessage_importance>0)
      cout << "add numproc " << name << ", type = " << np->GetClassName() << endl;
    np->SetName (name);
    numprocs.Set (name, np);

    todo.Append(np);
  }



  
  void PDE :: AddBilinearFormIntegrator (const string & name, BilinearFormIntegrator * part,
					 const bool deletable)
  {
    BilinearForm * form = GetBilinearForm (name);
    if (form && part)
      {
	form->AddIntegrator (part,deletable);
        if (printmessage_importance>0)
          cout << "integrator " << part->Name() << endl;
      }
    else
      {
        if (printmessage_importance>0)
          cerr << "Bilinearform = " << form << ", part = " << part << endl;
      }
  }
  
  void PDE :: AddIndependentBilinearFormIntegrator (const string & name, BilinearFormIntegrator * part,
						    const int master, const int slave,
						    const bool deletable)
  {
    BilinearForm * form = GetBilinearForm (name);
    if (form && part)
      {
	form->AddIndependentIntegrator (part,master,slave,deletable);
        if (printmessage_importance>0)
          cout << "integrator " << part->Name() << endl;
      }
    else
      {
        if (printmessage_importance>0)
          cerr << "Bilinearform = " << form << ", part = " << part << endl;
      }
  }

 
  void PDE :: AddLinearFormIntegrator (const string & name, LinearFormIntegrator * part)
  {
    LinearForm * form = GetLinearForm (name);
    if (form && part)
      {
	form->AddIntegrator (part);
        if (printmessage_importance>0)
          cout << "integrator " << part->Name() << endl;
      }
    else
      {
        if (printmessage_importance>0)
          cerr << "Linearform = " << form << ", part = " << part << endl;
      }
  }
  
  void PDE :: SetLineIntegratorCurvePointInfo(const string & filename,
					      Integrator * integrator)
  {
    CurvePointIntegrators.Append(integrator);
    CurvePointIntegratorFilenames.Append(new string(filename));
  }

  void PDE :: WritePDEFile ( string abspdefile, string geofile, 
		      string meshfile, string matfile, string oldpdefile )
  {
    ofstream pdeout ( abspdefile.c_str() );
    ifstream pdein ( oldpdefile.c_str() );

    pdeout << "geometry = " << geofile << endl;
    pdeout << "mesh = " << meshfile << endl;
    if ( matfile != "" )
      pdeout << "matfile = " << matfile << endl;

    string token;
    char ch;

    bool init = true;
    while ( init )
      {
	pdein.get(ch);
	if ( ch == '\n' )
	  continue;
	else if ( ch == '#' )
	  {
	    while ( ch != '\n' )
	      pdein.get(ch);
	    continue;
	  }
	pdein.putback(ch);
	pdein >> token;
	if ( token == "mesh" || token == "geometry" || token == "matfile" )
	  {
	    while ( ch != '\n' )
	      pdein.get(ch);
	    continue;
	  }
	pdeout << token;
	init = false;
      }

    // copy rest of file
    
    while ( pdein.good())
      {
	pdein.get(ch);
	pdeout.put(ch);
      }
  }

#ifdef ASTRID
  void PDE :: SaveZipSolution (const string & filename, const bool ascii )
  {
    string::size_type pos1 = filename.rfind('\\');
    string::size_type pos2 = filename.rfind('/');
    
    if (pos1 == filename.npos) pos1 = 0;
    if (pos2 == filename.npos) pos2 = 0;

    string oldmatfile = this -> GetMatfile();
    string oldpdefile = this -> GetFilename();
    string oldgeofile    = this -> GetGeoFileName();

    string pde_directory = filename.substr (0, max2(pos1, pos2));
    string dirname = filename.substr (max2(pos1, pos2) +1, filename.length() - max2(pos1,pos2) - 8);
    string absdirname = pde_directory +"/" + dirname;

    string meshfile = dirname + ".vol";
    string absmeshfile = absdirname + "/" + meshfile;
    string pdefile = dirname + ".pde";
    string geofile = dirname + ".geo";
    string solfile = dirname + ".sol";
    string abspdefile = absdirname + "/" + pdefile;
    string absgeofile = absdirname + "/" + geofile;
    string abssolfile = absdirname + "/" + solfile;
    string absmatfile = "";
    string matfile = "";


    string mkdir = "mkdir " + absdirname;
    system ( mkdir.c_str() );


    string cpgeo = "cp " + oldgeofile + " " + absgeofile;

    system ( cpgeo.c_str() );


    if ( oldmatfile != "" )
      {
	matfile = dirname + ".mat";
	absmatfile = absdirname + "/" + matfile;
	string cpmat = "cp " + oldmatfile + " " + absmatfile;
	system (cpmat.c_str() );
      }

    // save current mesh
    Ng_SaveMesh(absmeshfile.c_str());
    // write pdefile with changed geo, mesh and matfile
    (*this).WritePDEFile ( abspdefile, geofile, meshfile, matfile, oldpdefile );
    // save solution to .sol file
    (*this).SaveSolution ( abssolfile, ascii ); 

    string tarfiles = "cd " + pde_directory + "\ntar -cf " + dirname + ".tar " + dirname;
    //pdefile + " " + geofile + " " + matfile + " " + solfile + " " + meshfile; 
    string gzipfiles = "cd " + pde_directory + "\ngzip -f " + dirname + ".tar";
    //\nrm "  + pdefile + " " + geofile + " " + matfile + " " + solfile + " " + meshfile; 
    string rmdir = "cd " + absdirname + "\nrm " + pdefile + " " +  geofile + " " + matfile + " " + solfile + " " + meshfile + "\ncd ..\nrmdir " + dirname; 
    system (tarfiles.c_str() );
    system (gzipfiles.c_str() );
    system ( rmdir.c_str() );

    if (printmessage_importance>0)
    {
      cout << "saved geometry, mesh, pde and solution to file " << endl << pde_directory << "/" 
	 << dirname << ".tar.gz" << endl;
    }

  }
///
void PDE :: LoadZipSolution (const string & filename, const bool ascii)
{ 
  string::size_type pos1 = filename.rfind('\\');
  string::size_type pos2 = filename.rfind('/');
  
  if (pos1 == filename.npos) pos1 = 0;
  if (pos2 == filename.npos) pos2 = 0;

  string pde_directory = filename.substr (0, max2(pos1, pos2));
  string dirname = filename.substr (max2(pos1, pos2) +1, filename.length() - max2(pos1,pos2) - 8);
  string absdirname = pde_directory +"/" + dirname;
  
  string pdefile = dirname + "/" + dirname + ".pde";
  string solfile = dirname + "/" + dirname + ".sol";
  string abspdefile = pde_directory + "/" + pdefile;
  string abssolfile = pde_directory + "/" + solfile;
  
  string unzipfiles =  "cd " + pde_directory + "\ngzip -d " + dirname + ".tar.gz";
  string untarfiles =  "cd " + pde_directory + "\ntar -xf " + dirname + ".tar";
  string rmtar = "rm " + absdirname + ".tar";

  system (unzipfiles.c_str() );
  system (untarfiles.c_str() );
  system (rmtar.c_str());

  (*this).LoadPDE(abspdefile, 0, 0);
  (*this).LoadSolution(abssolfile, ascii );
}
#endif
  ///

}
