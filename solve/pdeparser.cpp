#include <solve.hpp>
#include <parallelngs.hpp>


#ifdef HAVE_DLFCN_H 
#include <dlfcn.h>
#endif


namespace ngsolve
{
  using namespace ngsolve;
  using namespace ngparallel;

  // parser for pde file


  enum TOKEN_TYPE 
    { 
      UNDEF = 0,
      NUMBER = 1, 
      STRING = 2,
      END = 3,  
  
      PLUS = '+', 
      MINUS = '-', 
      MUL = '*', DIV = '/',
      LP = '(', RP = ')', EQUAL = '=', COMMA = ',',
      LSB = '[', RSB = ']',
      COMMENT = '#',
      KW_DEFINE = 100, KW_GEOMETRY, KW_MESH, KW_SHARED,
      KW_CONSTANT, KW_VARIABLE, KW_COEFFICIENT, KW_FESPACE, KW_GRIDFUNCTION, 
      KW_BILINEARFORM, KW_LINEARFORM, KW_PRECONDITIONER, KW_BEMELEMENT,
      KW_INTEGRATOR = 200, KW_NUMPROC_ID,
      KW_NUMPROC = 300, 
      KW_MATFILE,
      KW_OVERLAP
    };

  
  struct kwstruct
  {
    TOKEN_TYPE kw; 
    const char * name;
  };

  static kwstruct defkw[] =
    {
      { KW_DEFINE,      "define" },
      { KW_GEOMETRY,    "geometry" },
      { KW_MESH,        "mesh" },
      { KW_SHARED,      "shared" },
      { KW_CONSTANT,    "constant" },
      { KW_VARIABLE,    "variable" },
      { KW_COEFFICIENT, "coefficient" },
      { KW_FESPACE,     "fespace" },
      { KW_GRIDFUNCTION, "gridfunction" },
      { KW_BILINEARFORM, "bilinearform" },
      { KW_BEMELEMENT,  "bemelement" },
      { KW_LINEARFORM,  "linearform" },
      { KW_PRECONDITIONER, "preconditioner" },
      { KW_NUMPROC,     "numproc" },
      { KW_MATFILE,     "matfile"},
      { KW_OVERLAP,     "overlap"},
      { TOKEN_TYPE(0) }
    };



  class PDEScanner
  {
    TOKEN_TYPE token;
    double num_value;
    string string_value;
    SymbolTable<TOKEN_TYPE> keywords;
    SymbolTable<int> integrators;
    SymbolTable<int> numprocs;

    int linenum;

    string copy_of_stream;

  private:
    void HandleStringConstants(void);

  public:
    istream * scanin;

    PDEScanner (istream * ascanin);

    ~PDEScanner();

    TOKEN_TYPE GetToken() const
    { return token; }

    double GetNumValue() const
    { return num_value; }

    const char * GetStringValueC() const
    { return string_value.c_str(); }
    const string & GetStringValue() const
    { return string_value; }

    void ReadNext();

    void Error (const string & err);
  };


  PDEScanner :: PDEScanner (istream * ascanin)
  {

#ifdef PARALLEL
      MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
      MPI_Comm_rank(MPI_COMM_WORLD, &id);
      hoprocs.SetSize(ntasks-1);
      for ( int i = 0; i<ntasks-1; i++ )
	hoprocs[i] = i+1;
#endif

    scanin = ascanin;
    token = END;
    num_value = 0;
    linenum = 1;

    // initialize keywords:

    kwstruct * kwp;
    kwp = defkw;
    while (kwp->kw)
      {
	keywords.Set (kwp->name, kwp->kw);
	kwp++;
      }

    // initialize integrator and numproc keywords:

    Integrators & itgs = GetIntegrators();
    for (int i = 0; i < itgs.GetBFIs().Size(); i++)
      integrators.Set (itgs.GetBFIs()[i]->name, 1);

    for (int i = 0; i < itgs.GetLFIs().Size(); i++)
      integrators.Set (itgs.GetLFIs()[i]->name, 1);

    NumProcs & nps = GetNumProcs();
    for (int i = 0; i < nps.GetNumProcs().Size(); i++)
      numprocs.Set (nps.GetNumProcs()[i]->name, 1);


    HandleStringConstants();
    scanin = new stringstream(copy_of_stream);
  }

  PDEScanner :: ~PDEScanner()
  {
    delete scanin;
  }

  void PDEScanner :: HandleStringConstants(void)
  {
    copy_of_stream = "";
    char ch;

    while(!scanin->eof())
      {
	scanin->get(ch);
	if(!scanin->eof())
	  copy_of_stream += ch;
      }

    size_t pos = 0;
    int ok = 0;
    string constname,constvalue;
    while(pos < copy_of_stream.size())
      {
	if(copy_of_stream[pos] == '#')
	  {
	    while(pos < copy_of_stream.size() && copy_of_stream[pos] != '\n')
	      pos++;
	  }
	else if(copy_of_stream[pos] == ' ' || copy_of_stream[pos] == '\t' || copy_of_stream[pos] == '\r' || copy_of_stream[pos] == '\n')
	  pos++;
	else if(ok == 2)
	  {
	    constname = "";
	    while(pos < copy_of_stream.size() && 
		  !(copy_of_stream[pos] == ' ' || copy_of_stream[pos] == '\t' || copy_of_stream[pos] == '\r' || copy_of_stream[pos] == '\n'))
	      {	      
		constname += copy_of_stream[pos];
		pos++;
	      }
	    ok = 3;
	  }
	else if(copy_of_stream.substr(pos,6) == "define")
	  {
	    pos += 6;
	    ok = 1;
	  }
	else if(ok == 1 && copy_of_stream.substr(pos,8) == "constant")
	  {
	    pos += 8;
	    ok = 2;
	  }
	else if(ok == 3 && copy_of_stream[pos] == '=')
	  {
	    pos++;
	    ok = 4;
	  }
	else if(ok == 4 && copy_of_stream[pos] == '\"')
	  {
	    pos++;
	    constvalue = "";
	    while(pos < copy_of_stream.size() && copy_of_stream[pos] != '\"')
	      {
		if(copy_of_stream[pos] == '\n' || copy_of_stream[pos] == '\t')
		  constvalue += ' ';
		else if(copy_of_stream[pos] == '\r')
		  ;
		else
		  constvalue += copy_of_stream[pos];
		pos++;
	      }

	    string fullname = string("$(")+constname+string(")");
	    size_t spos = pos;
	    while(spos != string::npos)
	      {
		spos = copy_of_stream.find(fullname,spos);

		if(spos != string::npos)
		  copy_of_stream.replace(spos,
					 strlen(fullname.c_str()),
					 constvalue);
	      }
	    ok = 0;
	  }
	else
	  {
	    ok = 0;
	    pos++;
	  }
      }
  }




  void PDEScanner :: ReadNext ()
  {
    char ch(0);
  

    // skip whitespaces
    do
      { 
	scanin->get(ch);

	if (ch == '\n') 
	  linenum++;

	// end of file reached
	if (scanin->eof())
	  {
	    token = END;
	    return;
	  }

	// skip comment line
	if (ch == '#')
	  {
	    while (ch != '\n')
	      {
		scanin->get(ch);
		if (scanin->eof())
		  {
		    token = END;
		    return;
		  }
	      }
	    linenum++;
	  }	
      }
    while (isspace(ch));
  
  
    //cout << "ch = \"" << ch << "\"" << endl;

    switch (ch)
      {
	// case '*': case '/': 
	// case '+': 
      case '-': 
      case '(': case ')':
      case '[': case ']':
      case '=': case ',':
	{
	  token = TOKEN_TYPE (ch);
	  break;
	}
      
      default:
	{
	  if (isdigit (ch) || ch == '.')
 	    {
	      scanin->putback (ch);

	      int pos = scanin->tellg();

	      bool isdig;
	      if(isdigit(ch))
		isdig = true;
	      else
		{
		  scanin->get(ch); // '.'
		  scanin->get(ch);
		  isdig = isdigit(ch);
		}

	      scanin->seekg(pos);

	      if(isdig)
		{
		  (*scanin) >> num_value;
		  token = NUMBER;
		  return;
		}
	      else
		{
		  scanin->get(ch);
		}
	    }


	  if (ch == '\"')
	    {
	      scanin->get(ch);
	      string_value = "";
	      while(ch != '\"')
		{
		  if(ch != '\r')
		    string_value += ch;
		  scanin->get(ch);
		  if(ch == '\n' || ch == '\t')
		    ch = ' ';
		}
	    }
	  else
	    {
	      (*scanin).putback (ch);
	      (*scanin) >> string_value;
	    }

          
          
          Integrators & itgs = GetIntegrators();
          for (int i = 0; i < itgs.GetBFIs().Size(); i++)
            if (itgs.GetBFIs()[i]->name == string_value)
              {
                token = KW_INTEGRATOR;
                return;
              }
          for (int i = 0; i < itgs.GetLFIs().Size(); i++)
            if (itgs.GetLFIs()[i]->name == string_value)
              {
                token = KW_INTEGRATOR;
                return;
              }

          // for (int i = 0; i < itgs.GetLFIs().Size(); i++)
          // integrators.Set (itgs.GetLFIs()[i]->name, 1);

          /*
	  if (integrators.Used (string_value))
	    {
	      token = KW_INTEGRATOR;
	      return;
	    }
          */

	  if (numprocs.Used (string_value))
	    {
	      token = KW_NUMPROC_ID;
	      return;
	    }
	
	  if (keywords.Used (string_value))
	    {
	      token = keywords[string_value];
	      return;
	    }
	
	  token = STRING;
	}
      }
  }

  void PDEScanner :: Error (const string & err)
  {
    stringstream errstr;
    errstr << "Parsing error in line " << linenum << ": " << endl << err << endl;
    errstr << "input continues with <<<";
    for (int i = 0; i < 50; i++)
      {
	char ch;
	scanin->get(ch);
	errstr << ch;
	if (scanin->eof())
	  {
	    errstr << "(end of file)";
	    break;
	  }
      }
    errstr << endl << ">>> stop parsing" << endl;
    throw Exception (errstr.str());
  }


  


  void CommandList (bool nomeshload = false, const bool nogeometryload = false);
  void DefineCommand ();
  void NumProcCommand ();
  void CheckFlags (Flags & flags);

  
  static PDEScanner * scan;
  static PDE * pde;
  

  void CommandList (bool nomeshload, const bool nogeometryload)
  {
    while (scan->GetToken() != END)
      {
	switch (scan->GetToken())
	  {

	  case KW_GEOMETRY:
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != '=')
		{
		  //cout << "Got Token: " << scan->GetToken() << endl;
		  scan->Error ("Expected =");
		}
	      scan->ReadNext();
	      	      
	      string geofile = pde->GetDirectory()+dirslash+scan->GetStringValue();

	      if(!nogeometryload)
		{
		  if (ifstream (geofile.c_str()))
		    {
		      Ng_LoadGeometry (const_cast<char*> (geofile.c_str()));
		      pde->SetGeoFileName(geofile);
		    }
		  else
		    {
		      Ng_LoadGeometry (const_cast<char*> (scan->GetStringValueC()));
		      pde->SetGeoFileName(scan->GetStringValue());
                    }
		}
	      scan->ReadNext();
	      break;	    }

	  case KW_MESH:
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != '=')
		scan->Error ("Expected '='");
	      scan->ReadNext();

	      // cout << "Load Mesh from File " << scan->GetStringValue() << endl;
	      // Ng_LoadMesh ((char*)scan->GetStringValueC());

	      string meshfile = pde->GetDirectory()+dirslash+scan->GetStringValue();

	      pde->SetMeshFileName(meshfile);
	      if(!nomeshload)
		{
		  if (ifstream (meshfile.c_str()))
		    {
                      if (printmessage_importance>0)
                        cout << "Load mesh from file " << meshfile << endl;
		      Ng_LoadMesh (const_cast<char*> (meshfile.c_str()));
		    }
		  else
		    {
                      if (printmessage_importance>0)
                        cout << "Load mesh from file " << scan->GetStringValue() << endl;
		      Ng_LoadMesh (const_cast<char*> (scan->GetStringValueC()));
		    }
		}

              pde->GetMeshAccess().UpdateBuffers();

	      if (!pde->GetMeshAccess().GetNP())
		throw Exception ("No mesh or empty mesh file\n");
	      scan->ReadNext();
	      break;
	    }

	  case KW_SHARED:
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != '=')
		scan->Error ("Expected '='");
	      scan->ReadNext();

	      // cout << "Load Mesh from File " << scan->GetStringValue() << endl;
	      // Ng_LoadMesh ((char*)scan->GetStringValueC());

	      // string shared = pde->GetDirectory()+dirslash+scan->GetStringValue();
              string shared = scan->GetStringValue() + ".so";
	      scan->ReadNext();

              cout << "load shared library '" << shared << "'" << endl;

#ifdef HAVE_DLFCN_H 
              void * handle = dlopen (shared.c_str(), RTLD_LAZY);
              if (!handle)
                {
                  stringstream err;
                  err << "Cannot load shared library '" << shared << "' \nerrmsg: "  << dlerror();
                  throw Exception (err.str());
                }
#else
              throw Exception ("cannot handle shared libraries");
#endif

              /*
              cout << "handle = " << handle << endl;
              if (!handle)
                cout << "dlerr = " << dlerror() << endl;
              */

              /*
              void (*symbolhandle)();
              symbolhandle =  ( void (*)() ) dlsym (handle, "My_DL_Init");
              cout << "symbolhandle = " << (void*)symbolhandle << endl;
              if (symbolhandle) (*symbolhandle)();
              */

              break;
            }

	  case KW_OVERLAP:
	    {
	      cout << "juhuuuuuuuuuu ... " << scan->GetStringValue() << endl;
	      scan->ReadNext();
	      if (scan->GetToken() != '=')
		scan->Error ("Expected '='");
	      scan->ReadNext();
#ifdef PARALLEL
	      int overlap = int(scan -> GetNumValue());
	      while (pde->GetMeshAccess().Overlap() < overlap && id == 0 && ntasks != 1)
		{
		  for ( int dest = 1; dest < ntasks; dest++)
		    {
		      MyMPI_Send ( "overlap++", dest );
		    }
		  Ng_UpdateOverlap();
		}
#endif
	      scan->ReadNext();
	      break;
	    }
	    
	  case KW_MATFILE:
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != '=')
		scan->Error ("Expected '='");
	      scan->ReadNext();
	      
	      string matfile  = scan->GetStringValue();
	      if (ifstream (matfile.c_str()).good())
		{
                  if (printmessage_importance>0)
                    cout << "Materialdata from file " << matfile << endl;
		  pde->SetMatfile(matfile); 
		}
	      else
		{
                  if (printmessage_importance>0)
                    cout << "Materialdata from file " << matfile << endl;
		  throw Exception("**** Materialdata-File not found!! \n ");  
		}
	      
	      scan->ReadNext();
	      
	      break;
	    } 
	    

	  case KW_DEFINE:
	    {
	      DefineCommand ();
	      break;
	    }

	  case KW_NUMPROC:
	    {
	      scan->ReadNext();
	      string npid = scan->GetStringValue();

	      scan->ReadNext();
	      string name = scan->GetStringValue();

	      scan->ReadNext();
	      Flags flags;
	      CheckFlags (flags);
	      
	      if (GetNumProcs().GetNumProc(npid))
		{
		  pde -> AddNumProc (name, GetNumProcs().GetNumProc(npid)->creator(*pde, flags));
// #ifdef SOCKETS
// 		  if(pde -> ConstantUsed ("clientserver") && pde -> GetConstant("clientserver") > 0.5)
// 		    pde -> GetClientSocketAccess().CheckNumProc(name,flags,pde->GetNumProc(name)->GetCallPosition());
// #endif
		}
	      else
		throw Exception (string("Undefined numproc ") + npid);
	      break;
	    }

	  default:
	    {
	      stringstream errstr;
	      errstr << "Unknown command";
	      if (scan->GetToken() == STRING)
		errstr << ": " << scan->GetStringValue() << endl;
	      else
		errstr << ", token = " << scan->GetToken() << endl;
	      scan -> Error (errstr.str());

	      scan->ReadNext();
	    }
	  }
      }
      if (printmessage_importance>0)
        cout << "End of file reached" << endl << endl << endl;
  }


		      
  



  void DefineCommand ()
  {
    scan->ReadNext();
    switch (scan->GetToken())
      {
      case KW_CONSTANT:
	{
	  scan->ReadNext();
	  string name = scan->GetStringValue ();

	  scan->ReadNext();
	  if (scan->GetToken() != '=')
	    scan->Error ("Expected =");

	  scan->ReadNext();

	  double val = 0;
	  string sval;

	  if (scan->GetToken() == LP)
	    {
	      EvalFunction * fun = new EvalFunction ();
	      for (int i = 0; i < pde->GetConstantTable().Size(); i++)
		fun->DefineConstant (pde->GetConstantTable().GetName(i),
				     pde->GetConstantTable()[i]);
	      for (int i = 0; i < pde->GetVariableTable().Size(); i++)
		fun->DefineGlobalVariable (pde->GetVariableTable().GetName(i),
					   &pde->GetVariableTable()[i]);
	      
	      fun->Parse (*scan->scanin);
	      if (fun->IsConstant())
		val = fun->Eval ((double*)(0));
	      else
		scan->Error  ("Expression not constant");
	    }
	  else if (scan->GetToken() == NUMBER)
	    val = scan->GetNumValue();
	  else if (scan->GetToken() == MINUS)
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != NUMBER)
		scan -> Error ("expected a number");
	      val = -scan->GetNumValue();
	    }
	  else if (scan->GetToken() == STRING)
	    {
	      sval = scan->GetStringValue();
	    }
	  else
	    scan->Error ("Expected a number");

	  if(sval == "")
	    pde->AddConstant (name, val);
	  else
	    pde->AddStringConstant (name, sval);
	    
	  scan->ReadNext();
	  break;
	}

      case KW_VARIABLE:
	{
	  scan->ReadNext();

	  string name = scan->GetStringValue ();

	  scan->ReadNext();
	  if (scan->GetToken() != '=')
	    scan->Error ("Expected =");

	  scan->ReadNext();

	  double val = 0;
	  
	  EvalVariable * eval = NULL;

	  if (scan->GetToken() == LP)
	    {
	      eval = new EvalVariable(pde->GetMeshAccess(),name);
	      for (int i = 0; i < pde->GetConstantTable().Size(); i++)
		eval->GetEvaluator().DefineConstant (pde->GetConstantTable().GetName(i),
						     pde->GetConstantTable()[i]);
	      for (int i = 0; i < pde->GetVariableTable().Size(); i++)
		eval->GetEvaluator().DefineGlobalVariable (pde->GetVariableTable().GetName(i),
							   &pde->GetVariableTable()[i]);
	      
	      
	      eval->GetEvaluator().Parse (*scan->scanin);
	    }
	  else if (scan->GetToken() == NUMBER)
	    val = scan->GetNumValue();
	  else if (scan->GetToken() == MINUS)
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != NUMBER)
		scan -> Error ("expected a number");
	      val = -scan->GetNumValue();
	    }
	  else
	    scan->Error ("Expected a number");

	  if(eval)
	    pde->AddVariable(name,eval);
	  else
	    pde->AddVariable (name, val);
	  scan->ReadNext();
	  break;
	}

      case KW_COEFFICIENT:
	{
	  Array<double> dcoeffs;
	  Array<EvalFunction*> coeffs;
	  Array < Array < Array<double>* >* > polycoeffs;
	  Array < Array<double>* > polybounds;
	  scan->ReadNext();

	  string name = scan->GetStringValue ();
	  scan->ReadNext();

	  if (strcmp (scan->GetStringValueC(), "material") == 0)
	    {
	      // define values by material

	      scan->ReadNext();
	      //SymbolTable<double> values;
	      SymbolTable<EvalFunction * > funs;
	      while (scan->GetToken() == STRING)
		{
		  string mat = scan->GetStringValue();
		  scan->ReadNext();

		  //double val = 1;

		  if (scan->GetToken() == LP)
		    {
		      EvalFunction * fun = new EvalFunction ();
		      for (int i = 0; i < pde->GetConstantTable().Size(); i++)
			fun->DefineConstant (pde->GetConstantTable().GetName(i),
					     pde->GetConstantTable()[i]);
		      for (int i = 0; i < pde->GetVariableTable().Size(); i++)
			fun->DefineGlobalVariable (pde->GetVariableTable().GetName(i),
						   &pde->GetVariableTable()[i]);
		      
		      
		      fun->Parse (*scan->scanin);
		      funs.Set (mat,fun);
// 		      if (fun->IsConstant())
// 			val = fun->Eval (NULL);
// 		      else
// 			scan->Error  ("Expression not constant");
		    }
		  else 
		    {
		      double val = 1;
		      if (scan->GetToken() == MINUS)
			{
			  val = -1;
			  scan->ReadNext();
			}
		      if (scan->GetToken() != NUMBER)
			scan->Error ("number expected");
		      val *= scan->GetNumValue();

		      EvalFunction * fun = new EvalFunction ();
		      fun->AddConstant (val);
		      funs.Set (mat,fun);
		      
		    }
		  scan->ReadNext();

		  if (scan->GetToken() == COMMA) 
		    scan->ReadNext();
		  //values.Set (mat, val);
		}

	      //cout << "values = " << endl << values << endl;

	      int maxdom = -1;
	      int ne = pde->GetMeshAccess().GetNE();
	      for (int i = 0; i < ne; i++)
		maxdom = max2 (maxdom, pde->GetMeshAccess().GetElIndex(i));
	      maxdom++;

	      dcoeffs.SetSize(maxdom);
	      dcoeffs = 0;
	      coeffs.SetSize(maxdom);
	      coeffs = NULL;
	      bool only_constant = true;
	      for (int i = 0; i < ne; i++)
		{
		  EvalFunction * fun = NULL;
		  string mat = pde->GetMeshAccess().GetElMaterial(i);
		  int index = pde->GetMeshAccess().GetElIndex(i);
		  // cout << "mat = " << mat << ", ind = " << index << endl;


		  bool used = false;
		  for(int j=0; !used && j<funs.Size(); j++)
		    {
		      used = StringFitsPattern(mat,funs.GetName(j));
		      if(used)
			{
			  //(*testout) << "\"" << mat << "\" fits \"" << funs.GetName(j) << "\"" << endl;
			  fun = funs[j];
			}
		    }
		  if(!used)
		    {
		      if (funs.Used ("default"))
			fun = funs["default"];
		      else
			throw Exception (string ("No value defined for material ")+mat);
		    }
		  
			 
		  if(coeffs[index] == NULL) coeffs[index] = new EvalFunction(*fun);
		  if(fun->IsConstant())
		    dcoeffs[index] = fun->Eval( (double*)(0));
		  else
		    only_constant = false;
		}

	      for(int i=0; i<funs.Size(); i++)
		delete funs[i];

	      if(only_constant)
		{
		  (*testout) << "material coefficients = " << endl << dcoeffs << endl;
		  pde->AddCoefficientFunction
		    (name, new DomainConstantCoefficientFunction(dcoeffs));
		}
	      else
		{
		  (*testout) << "material coefficients variable " << endl;
		  if (pde->GetMeshAccess().GetDimension() == 2)
		    pde->AddCoefficientFunction
		      (name, new DomainVariableCoefficientFunction<2>(coeffs));
		  else
		    pde->AddCoefficientFunction
		      (name, new DomainVariableCoefficientFunction<3>(coeffs));
		}
	      
	      for (int hi = 0; hi < coeffs.Size(); hi++)
		delete coeffs[hi];
	      
	    }

	  else if (strcmp (scan->GetStringValueC(), "files") == 0)
	    {
	      string ipfilename,infofilename,valuesfilename;

// 	      scan -> ReadNext();
// 	      ipfilename = scan->GetStringValue();
// 	      scan -> ReadNext();
// 	      infofilename = scan->GetStringValue();
// 	      scan -> ReadNext();
// 	      valuesfilename = scan->GetStringValue();
	      
 	      scan -> ReadNext();
 	      Flags flags;
 	      CheckFlags (flags);

	      ipfilename = pde->GetDirectory()+dirslash+flags.GetStringFlag("ipfile","ipfile");
	      infofilename = pde->GetDirectory()+dirslash+flags.GetStringFlag("infofile","infofile");
	      valuesfilename = pde->GetDirectory()+dirslash+flags.GetStringFlag("valuesfile","valuesfile");

	      
	      const bool loadvalues = flags.GetDefineFlag("loadvalues");
	      


	      pde->AddCoefficientFunction
		(name, new FileCoefficientFunction(ipfilename,
						   infofilename,
						   valuesfilename,
						   loadvalues));
	      
	      //scan -> ReadNext();
	    }

	  else if (strcmp (scan->GetStringValueC(), "bcnames") == 0)

	    {
	      // define values by names of boundary conditions

	      scan->ReadNext();
	      //SymbolTable<double> values;
	      SymbolTable<EvalFunction * > funs;
	      while (scan->GetToken() == STRING)
		{
		  string bcname = scan->GetStringValue();
		  scan->ReadNext();

		  //double val = 1;

		  if (scan->GetToken() == LP)
		    {
		      EvalFunction * fun = new EvalFunction ();
		      for (int i = 0; i < pde->GetConstantTable().Size(); i++)
			fun->DefineConstant (pde->GetConstantTable().GetName(i),
					     pde->GetConstantTable()[i]);
		      for (int i = 0; i < pde->GetVariableTable().Size(); i++)
			fun->DefineGlobalVariable (pde->GetVariableTable().GetName(i),
						   &pde->GetVariableTable()[i]);
		      
		      
		      fun->Parse (*scan->scanin);
		      funs.Set (bcname,fun);
// 		      if (fun->IsConstant())
// 			val = fun->Eval (NULL);
// 		      else
// 			scan->Error  ("Expression not constant");
		    }
		  else 
		    {
		      double val = 1;
		      if (scan->GetToken() == MINUS)
			{
			  val = -1;
			  scan->ReadNext();
			}
		      if (scan->GetToken() != NUMBER)
			scan->Error ("number expected");
		      val *= scan->GetNumValue();

		      EvalFunction * fun = new EvalFunction ();
		      fun->AddConstant (val);
		      funs.Set (bcname,fun);
		      
		    }
		  scan->ReadNext();

		  if (scan->GetToken() == COMMA) 
		    scan->ReadNext();
		  //values.Set (bcname, val);
		}

	      //cout << "values = " << endl << values << endl;

	      int maxbc = -1;
	      int nse = pde->GetMeshAccess().GetNSE();
	      for (int i = 0; i < nse; i++)
		maxbc = max2 (maxbc, pde->GetMeshAccess().GetSElIndex(i));
	      maxbc++;

	      dcoeffs.SetSize(maxbc);
	      dcoeffs = 0;
	      coeffs.SetSize(maxbc);
	      coeffs = NULL;
	      bool only_constant = true;
	      for (int i = 0; i < nse; i++)
		{
		  EvalFunction * fun = NULL;
		  string bcname = pde->GetMeshAccess().GetSElBCName(i);
		  int index = pde->GetMeshAccess().GetSElIndex(i);
		  // cout << "bcname = " << bcname << ", ind = " << index << endl;

		  bool used = false;
		  for(int j=0; !used && j<funs.Size(); j++)
		    {
		      used = StringFitsPattern(bcname,funs.GetName(j));
		      if(used)
			fun = funs[j];
		    }
		  if(!used)
		    {
		      if (funs.Used ("default"))
			fun = funs["default"];
		      else
			throw Exception (string ("No value defined for boundary condition ")+bcname);
		    }
			 
		  if(coeffs[index] == NULL) coeffs[index] = new EvalFunction(*fun);
		  if(fun->IsConstant())
		    dcoeffs[index] = fun->Eval( (double*)(0) );
		  else
		    only_constant = false;
		}

	      for(int i=0; i<funs.Size(); i++)
		delete funs[i];

	      if(only_constant)
		{
		  (*testout) << "material coefficients = " << endl << dcoeffs << endl;
		  pde->AddCoefficientFunction
		    (name, new DomainConstantCoefficientFunction(dcoeffs));
		}
	      else
		{
		  (*testout) << "material coefficients variable " << endl;
		  if (pde->GetMeshAccess().GetDimension() == 2)
		    pde->AddCoefficientFunction
		      (name, new DomainVariableCoefficientFunction<2>(coeffs));
		  else
		    pde->AddCoefficientFunction
		      (name, new DomainVariableCoefficientFunction<3>(coeffs));
		}
	      
	      for (int hi = 0; hi < coeffs.Size(); hi++)
		delete coeffs[hi];
	      
	    }

	  else
	    
	    {
	      while (scan->GetToken() == NUMBER || 
		     scan->GetToken() == MINUS ||
		     scan->GetToken() == LP ||
		     scan->GetToken() == LSB)
		{
		  if (scan->GetToken() == LP)
		    {
		      cout << "a" << endl;
		      EvalFunction * fun = new EvalFunction ();
		      for (int i = 0; i < pde->GetConstantTable().Size(); i++)
			fun->DefineConstant (pde->GetConstantTable().GetName(i),
					     pde->GetConstantTable()[i]);
		      for (int i = 0; i < pde->GetVariableTable().Size(); i++)
			fun->DefineGlobalVariable (pde->GetVariableTable().GetName(i),
						   &pde->GetVariableTable()[i]);
		      
		      cout << "a1" << endl;
		      fun->Parse (*scan->scanin);
		      cout << "a2" << endl;
		      coeffs.Append (fun);
		      cout << "a3" << endl;
		      fun -> Print (cout);
		      if (fun->IsConstant())
			dcoeffs.Append (fun->Eval ( (double*)(0) ));
		      scan->ReadNext();
		      cout << "b" << endl;
		      // fun -> Print(cout);
		    }
		  else if (scan->GetToken() == LSB) // polynomial [MW]
		    {
                      if (printmessage_importance>0)
                        cout << "polynomial: ";
		      scan->ReadNext();
		      double val;
		      Array< Array<double>* > *polyco = new Array< Array<double>* >;
		      Array<double> *polyb = new Array<double>;
		      while(scan->GetToken() != RSB)
			{
			  Array<double> *polyc = new Array<double>;
			  while(scan->GetToken() != RSB && scan->GetToken() != LP)
			    {
			      while (scan->GetToken() == COMMA) 
				scan->ReadNext();
			      val = scan->GetNumValue();
			      //cout << "val " << val << endl;
			      if (scan->GetToken() == MINUS)
				{
				  scan->ReadNext();
				  if (scan->GetToken() != NUMBER)
				    scan -> Error ("expected a number");
				  
				  val = - scan->GetNumValue();
				}
			      polyc->Append(val);
			      scan->ReadNext();
			    }
                          if (printmessage_importance>0)
                            cout << (*polyc)[0];
                          
			  for(int i = 1; i<polyc->Size(); i++)
			    {
                              if (printmessage_importance>0)
                                cout << " + " << (*polyc)[i] << "*t^"<<i;
			    }
			  

			  polyco->Append(polyc);
			  if(scan->GetToken() == LP)
			    {
			      scan->ReadNext();
			      val = scan->GetNumValue();
			      if (scan->GetToken() == MINUS)
				{
				  scan->ReadNext();
				  if (scan->GetToken() != NUMBER)
				    scan -> Error ("expected a number");
				  
				  val = - scan->GetNumValue();
				}
			      polyb->Append(val);
                              if (printmessage_importance>0)
                                cout << " until " << val <<", then ";
			      scan->ReadNext(); // RP
			      scan->ReadNext();
			    }
			}
		      polybounds.Append(polyb);
		      polycoeffs.Append(polyco);

                      if (printmessage_importance>0)
                        cout << endl;

		      scan->ReadNext();
		    }
		  else
		    {
		      double val = scan->GetNumValue();
		      if (scan->GetToken() == MINUS)
			{
			  scan->ReadNext();
			  if (scan->GetToken() != NUMBER)
			    scan -> Error ("expected a number");
			  
			  val = - scan->GetNumValue();
			}
		      
		      EvalFunction * fun = new EvalFunction();
	 	      fun->AddConstant (val);
		      coeffs.Append (fun);
		      dcoeffs.Append (val);
		      scan->ReadNext();
		    }
		  
		  if (scan->GetToken() == COMMA) 
		    scan->ReadNext();
		}

	      bool allconst = 1;
	      for (int i = 0; i < coeffs.Size(); i++)
		if (!coeffs[i]->IsConstant())
		  allconst = 0;
	      
	      if (polycoeffs.Size() > 0)
		{
		  pde->AddCoefficientFunction
		    (name, new PolynomialCoefficientFunction(polycoeffs,polybounds));
		}
	      else if (allconst)
		{
		  pde->AddCoefficientFunction
		    (name, new DomainConstantCoefficientFunction(dcoeffs));
		  for (int hi = 0; hi < coeffs.Size(); hi++)
		    delete coeffs[hi];
		}
	      else
		{
		  if (pde->GetMeshAccess().GetDimension() == 2)
		    {
		      pde->AddCoefficientFunction
			(name, new DomainVariableCoefficientFunction<2>(coeffs));
		    }
		  else
		    pde->AddCoefficientFunction
		      (name, new DomainVariableCoefficientFunction<3>(coeffs));

		  for (int hi = 0; hi < coeffs.Size(); hi++)
		    delete coeffs[hi];
		}
	    }

	  break;
	}

      case KW_FESPACE:
	{
	  scan->ReadNext();

	  string name = scan->GetStringValue ();
	  scan->ReadNext();
	  Flags flags;
	  CheckFlags (flags);
	  pde->AddFESpace (name, flags);
	  break;
	}

      case KW_GRIDFUNCTION:
	{
	  scan->ReadNext();
	  string name = scan->GetStringValue ();
	  scan->ReadNext();

	  Flags flags;
	  CheckFlags (flags);
	  pde->AddGridFunction (name, flags);
	  break;
	}

      case KW_BILINEARFORM:
	{
	  scan->ReadNext();

	  string name = scan->GetStringValue ();
	  scan->ReadNext();

	  Flags flags;
	  CheckFlags (flags);
	  pde->AddBilinearForm (name, flags);


	  //
	  if(flags.StringFlagDefined("useintegratorsof"))
	    {
	      BilinearForm * source = pde->GetBilinearForm(flags.GetStringFlag("useintegratorsof",""));
	      
	      for(int i = 0; i < source->NumIntegrators(); i++)
		pde->AddBilinearFormIntegrator(name,source->GetIntegrator(i),false);
	    
	      for(int i = 0; i < source->NumIndependentIntegrators(); i++)
		pde->AddBilinearFormIntegrator(name,source->GetIndependentIntegrator(i),false);
	    }

	  //
	
	  // read bilinear-form components
	  int ncoeffs;
	  while (1)
	    {
	      ncoeffs = -1;

	      TOKEN_TYPE integrator_token = scan->GetToken();

	      if (integrator_token == KW_INTEGRATOR)
		{
		  string integrator_name = scan->GetStringValue();
		  const ngfem::Integrators::IntegratorInfo * info =
		    ngfem::GetIntegrators() 
		    . GetBFI(integrator_name, 
			     pde->GetMeshAccess().GetDimension());

		  if (info)
		    {
		      Array<CoefficientFunction*> coeffs(info->numcoeffs);
		      scan->ReadNext();
		      for (int i = 0; i < info->numcoeffs; i++)
			{
			  coeffs[i] = pde->GetCoefficientFunction (scan->GetStringValue(), 1);
			  if (!coeffs[i])
			    {
			      T_GridFunction<double> * gf = 
				dynamic_cast<T_GridFunction<double>*> (pde->GetGridFunction (scan->GetStringValue(), 1));
			      if (gf)
				coeffs[i] = new GridFunctionCoefficientFunction (*gf);
			      else
				{
				  throw Exception (string("undefined coefficient ") + scan->GetStringValue());
				}
			    }

			  scan->ReadNext();
			}
		    
		      Flags partflags;
		      CheckFlags (partflags);

		      ngfem::BilinearFormIntegrator * integrator = 
			dynamic_cast<ngfem::BilinearFormIntegrator*> (info->creator(coeffs));
		      integrator -> SetName (integrator_name);

		      if (partflags.NumFlagDefined ("order"))
			{
			  integrator -> 
			    SetIntegrationOrder
			    (int (partflags.GetNumFlag ("order",0)));
			}

		      if (partflags.NumFlagDefined ("higherorder"))
			{
			  integrator -> 
			    SetHigherIntegrationOrder
			    (int (partflags.GetNumFlag ("higherorder",0)));
			}

		      if (partflags.StringFlagDefined ("filename"))
			{
#ifdef WIN32
			  char slash = '\\';
#else
			  char slash = '/';
#endif
			  integrator -> SetFileName(pde->GetDirectory() + slash + partflags.GetStringFlag("filename",""));
			}


		      if (partflags.NumFlagDefined ("comp"))
			{
			  if (dynamic_cast<const CompoundFESpace*> 
			      (&pde->GetBilinearForm (name)->GetFESpace()))
			    {
			      integrator = new CompoundBilinearFormIntegrator
				(*integrator, 
				 int(partflags.GetNumFlag ("comp", 1))-1);
			    }
			  else
			    {
			      integrator = new BlockBilinearFormIntegrator
				(*integrator, 
				 pde->GetBilinearForm (name)->GetFESpace().GetDimension(),
				 int(partflags.GetNumFlag ("comp", 1))-1);
			    }
			}

		      if (partflags.GetDefineFlag ("normal"))
			{
			  throw Exception ("-normal flag currently not available");
			  /*
			  integrator = new NormalBilinearFormIntegrator
			    (*integrator);
			  */
			}

		      if (partflags.GetDefineFlag ("imag"))
			{
			  integrator = new ComplexBilinearFormIntegrator
			    (*integrator, Complex(0,1));
			}


		      if (partflags.NumFlagDefined ("definedon"))
			{
			  int domain = int(partflags.GetNumFlag("definedon", 0));
			  int size = max2 (pde->GetMeshAccess().GetNDomains(), 
					   pde->GetMeshAccess().GetNBoundaries());
			  BitArray definedon(size);
			  definedon.Clear();
			  definedon.Set(domain-1);
			  (*testout) << "definedon = " << definedon << endl;
			  integrator->SetDefinedOn (definedon);
			}

		      // integrator -> SetFastIntegration (partflags.GetDefineFlag("fast"),
		      // partflags.GetDefineFlag("checkfast"));

		      integrator -> SetConstantCoefficient (partflags.GetDefineFlag("const"));


// 		      bool indep(false);
// 		      if (partflags.NumFlagDefined ("master") && partflags.NumFlagDefined ("slave"))
// 			{
// 			  indep = true;
// 			  pde->AddIndependentBilinearFormIntegrator(name, integrator, 
// 								    static_cast<int>(partflags.GetNumFlag("master",0))-1,
// 								    static_cast<int>(partflags.GetNumFlag("slave",0))-1);
// 			}

// 		      if(!indep)
//                      {

		      pde->AddBilinearFormIntegrator (name, integrator);

		      continue;
		    }
		  else
		    {
		      throw Exception ("Should have integrator");
		    }
		

		}
	      break;
	 
	    }
		   
	  break;
	}

      case KW_LINEARFORM:
	{
	  scan->ReadNext();
	  string name = scan->GetStringValue ();

	  scan->ReadNext();
	  Flags flags;
	  CheckFlags (flags);

	  pde->AddLinearForm (name, flags);
	
	  // read linear-form components
	  int ncoeffs;
	  while (1)
	    {
	      ncoeffs = -1;

	      TOKEN_TYPE integrator_token = scan->GetToken();

	      if (integrator_token == KW_INTEGRATOR)
		{

		  const ngfem::Integrators::IntegratorInfo * info =
		    ngfem::GetIntegrators() 
		    . GetLFI(scan->GetStringValue(), 
			     pde->GetMeshAccess().GetDimension());
		
		  if (info)
		    {
		      Array<CoefficientFunction*> coeffs(info->numcoeffs);
		      scan->ReadNext();
		      for (int i = 0; i < info->numcoeffs; i++)
			{
			  coeffs[i] = pde->GetCoefficientFunction (scan->GetStringValue(),true);
			  if (!coeffs[i])
			    {
			      T_GridFunction<double> * gf = 
				dynamic_cast<T_GridFunction<double>*> (pde->GetGridFunction (scan->GetStringValue(), 1));
			      if (gf)
				coeffs[i] = new GridFunctionCoefficientFunction (*gf);
			      else
				throw Exception (string("undefined coefficient ") + scan->GetStringValue());
			    }
			  scan->ReadNext();

			  // for vector valued coefficient functions
			  if (i == 0 &&
			      coeffs[i]->Dimension() == info->numcoeffs)
			    {
			      coeffs.SetSize(1);
			      break;
			    }
			}
		    
		      Flags partflags;
		      CheckFlags (partflags);
		    
		      ngfem::LinearFormIntegrator * integrator = 
			dynamic_cast<ngfem::LinearFormIntegrator*> (info->creator(coeffs));

		      (*testout) << "partflags = " << endl << partflags << endl;
		    
		      if (partflags.NumFlagDefined ("order"))
			{
			  integrator -> 
			    SetIntegrationOrder
			    (int (partflags.GetNumFlag ("order",0)));
			}



		      if (partflags.NumFlagDefined ("comp"))
			{
			  if (dynamic_cast<const CompoundFESpace*> 
			      (&pde->GetLinearForm (name)->GetFESpace()))
			    {
			      integrator = new CompoundLinearFormIntegrator
				(*integrator, 
				 int(partflags.GetNumFlag ("comp", 1))-1);
			    }
			  else
			    {
			      integrator = new BlockLinearFormIntegrator
				(*integrator, 
				 pde->GetLinearForm (name)->GetFESpace().GetDimension(),
				 int(partflags.GetNumFlag ("comp", 1))-1);
			    }
			}

		      if (partflags.NumFlagDefined ("normal"))
			{
			  throw Exception ("-normal flag currently not available");
			}
		    
		      if (partflags.GetDefineFlag ("imag"))
			{
			  
			  integrator = new ComplexLinearFormIntegrator
			    (*integrator, Complex(0,1));
			} 


		      if (partflags.NumFlagDefined ("definedon") || partflags.NumListFlagDefined("definedon"))
			{
			  int size = max2 (pde->GetMeshAccess().GetNDomains(), 
					   pde->GetMeshAccess().GetNBoundaries());
			  BitArray definedon(size);
			  definedon.Clear();

			  if(partflags.NumFlagDefined ("definedon"))
			    definedon.Set(int(partflags.GetNumFlag("definedon", 0))-1);
			  
			  if(partflags.NumListFlagDefined("definedon"))
			    for(int i=0; i<partflags.GetNumListFlag("definedon").Size(); i++)
			      definedon.Set(int(partflags.GetNumListFlag("definedon")[i])-1);
			      

			  integrator->SetDefinedOn (definedon);
			}

		      if (partflags.NumFlagDefined ("cachecomp"))
			{
			  integrator->SetCacheComp(int(partflags.GetNumFlag("cachecomp",0)));
			}
		    
		      if (partflags.StringFlagDefined ("curvefile"))
			{		
			  pde->SetLineIntegratorCurvePointInfo(pde->GetDirectory() + dirslash + partflags.GetStringFlag("curvefile",""),
			  				       integrator);

			  //BuildLineIntegratorCurvePoints ( pde->GetDirectory() + dirslash + partflags.GetStringFlag("curvefile",""),
			  //				   pde->GetMeshAccess(),
			  //				   *integrator);
			}

		      // integrator -> SetFastIntegration (partflags.GetDefineFlag("fast"),
		      // partflags.GetDefineFlag("checkfast"));

		      integrator -> SetConstantCoefficient (partflags.GetDefineFlag("const"));

		      pde->AddLinearFormIntegrator (name, integrator);
		      continue;
		    }
		  else
		    {
		      throw Exception ("Should have integrator");
		    }

		}

	      break;

	    }
		   
	  break;
	}


      case KW_PRECONDITIONER:
	{
	  scan->ReadNext();
	  string name = scan->GetStringValue ();
	
	  scan->ReadNext();
	  Flags flags;

	  CheckFlags (flags);

	  pde->AddPreconditioner (name, flags);

	  break;	
	}


	/*
      case KW_BEMELEMENT:
	{
	  scan->ReadNext();
	  string name = scan->GetStringValue ();

	  scan->ReadNext();
	  Flags flags;

	  CheckFlags (flags);
	  pde->AddBEMElement (name, flags);

	  break;
	}
	*/

      case STRING:
	{
	  scan -> Error ("Illegal define command " + string(scan->GetStringValue()));
	  break;
	}
      default:
	{
	  stringstream errstr;
	  errstr << "Illegal define token " << scan->GetStringValue();
	  scan -> Error (errstr.str());
	}
      }
  }





  void CheckFlags (Flags & flags)
  {
    while (scan->GetToken() == '-')
      {
	scan->ReadNext();
	string flag = string("-") + scan->GetStringValue();
	flags.SetCommandLineFlag (flag.c_str());
	scan->ReadNext();
      }
  }


  
  void BuildLineIntegratorCurvePoints ( const string filename,
					const MeshAccess & ma,
					Integrator & integrator,
					bool draw)
  {
    ifstream infile(filename.c_str());
    

    if(!infile)
      {
	string errstring = string("Error for integration along curve: could not open \"") +
	  filename + string("\"\n");
	cerr << errstring;
	throw Exception(errstring);
      }
    
    BuildLineIntegratorCurvePoints(infile,ma,integrator,draw);

    infile.close();
  }

  void BuildLineIntegratorCurvePoints ( istream & infile,
					const MeshAccess & ma,
					Integrator & integrator,
					bool draw)
  {
    Array<int> domains;

    if(integrator.DefinedOnSubdomainsOnly())
      {
	for(int i=0; i<ma.GetNDomains(); i++)
	  if(integrator.DefinedOn(i))
	    domains.Append(i);
      }

    integrator.UnSetIntegrationAlongCurve();

    if(draw)
      ma.InitPointCurve();

    Vector<> current(3);

    bool latestpoint = true;

    string inp;
    infile >> inp;
    while(infile)
      {
	if(inp[0] == '#')
	  {
	    infile.ignore(1000,'\n');
	  }
	else if(inp == "p")
	  {
	    latestpoint = true;
	    for(int i=0; i<3; i++)
	      infile >> current[i];

	    Vector<> tangent(3);
	    for(int i=0; i<3; i++)
	      infile >> tangent[i];

	    integrator.AppendCurvePoint(current,tangent);
	    if(draw)
	      ma.AddPointCurvePoint(current);
	  }
	else if (inp == "l")
	  {
	    if(!latestpoint)
	      {
		if(draw)
		  ma.InitPointCurve();
		integrator.SetCurveClearance();
	      }
	    latestpoint = false;
	    double secfac;
	    infile >> secfac;
	    Vector<> start(3), end(3), v2(3);
	    
	    for(int i=0; i<3; i++)
	      infile >> start[i];
	    for(int i=0; i<3; i++)
	      infile >> end[i];

	    v2 = end - start;
	    double l = L2Norm(v2);

	    v2 *= 1./l;

	    IntegrationPoint dummyip;

	    current = start;
	    
	    Array<int> verts;

	    int numpoints = 0;

	    int oldelement=-1;
	    double h = 0;
	    if(draw)
	      ma.AddPointCurvePoint(start);

	    while(L2Norm(current - start) < l)
	      {
		integrator.AppendCurvePoint(current,v2);
		numpoints++;

		int element;
		if(domains.Size() > 0)
		  element = ma.FindElementOfPoint(current,dummyip,true,&domains);
		else
		  element = ma.FindElementOfPoint(current,dummyip,true);


		if(element == -1)
		  {
		    cerr << endl << "WARNING: Integration line may lie outside domain! (at point "<<current <<")" << endl << endl;
		  }

		if(element != oldelement && element != -1)
		  {		
		    ma.GetElVertices(element,verts);
		    Vec<3> center = 0;
		    for(int i = 0; i<verts.Size(); i++)
		      {
			Vec<3> vertex;
			ma.GetPoint(verts[i],vertex);
			center += vertex;
		      }
		    center *= 1./double(verts.Size());
		    h=1e10;
		    for(int i = 0; i<verts.Size(); i++)
		      {
			Vec<3> vertex;
			ma.GetPoint(verts[i],vertex);
			double auxh = L2Norm(vertex-center);
			if(auxh < h) h = auxh;
		      }
		    
		    h *= secfac;
		  }
		oldelement = element;

		
		current += h*v2;
	      }
	    integrator.AppendCurvePoint(end,v2);
	    if(draw)
	      ma.AddPointCurvePoint(end);
	    numpoints++;

	    (*testout) << "curve-line with " << numpoints << " integration points" << endl;
	  }
	else if (inp == "s")
	  {
	    if(!latestpoint)
	      {
		if(draw)
		  ma.InitPointCurve();
		integrator.SetCurveClearance();
	      }
	    latestpoint = false;
	    double secfac;
	    infile >> secfac;
	    Vector<> p1(3), p2(3), p3(3);
	    for(int i=0; i<3; i++)
	      infile >> p1[i];
	    for(int i=0; i<3; i++)
	      infile >> p2[i];
	    for(int i=0; i<3; i++)
	      infile >> p3[i];

	    
	    IntegrationPoint dummyip;
	    double t = 0;
	    double tstep = 0.1;
	    int oldelement=-1;
	    double h = 0;
            
	    Array<int> verts;
	    

	    Vector<> oldp(3);
	    oldp = p1;

	    Vector<> tangent(3);
	    tangent = p2-p1;
	    integrator.AppendCurvePoint(p1,tangent);
	    if(draw)
	      ma.AddPointCurvePoint(p1);
	    int numpoints = 1;

	    while (t < 1)
	      {		
		int element;
		if(domains.Size() > 0)
		  element = ma.FindElementOfPoint(oldp,dummyip,true,&domains);
		else
		  element = ma.FindElementOfPoint(oldp,dummyip,true);

		if(element != oldelement)
		  {		
		    ma.GetElVertices(element,verts);
		    Vec<3> center = 0;
		    for(int i = 0; i<verts.Size(); i++)
		      {
			Vec<3> vertex;
			ma.GetPoint(verts[i],vertex);
			center += vertex;
		      }
		    center *= 1./double(verts.Size());
		    h=1e10;
		    for(int i = 0; i<verts.Size(); i++)
		      {
			Vec<3> vertex;
			ma.GetPoint(verts[i],vertex);
			double auxh = L2Norm(vertex-center);
			if(auxh < h) h = auxh;
		      }
		    h *= secfac;
		  }
		oldelement = element;

		tstep *= 2;
		do
		  {
		    tstep *= 0.75;
		    double testt = t+tstep;
		    double b1 = (1.-testt)*(1.-testt);
		    double b2 = sqrt(2.)*testt*(1.-testt);
		    double b3 = testt*testt;
		    double w = b1+b2+b3;
		    current = (b1/w)*p1 + (b2/w)*p2 + (b3/w)*p3;
		  } 
		while(L2Norm(current-oldp) > h);
		
		t+=tstep;
		if(t < 1)
		  {
		    tangent = (1.-t)*((1.-sqrt(2.))*t-1.)*p1 +
		      (1.-2*t)*p2 +
		      (t*(sqrt(2.)+(1.-sqrt(2.))*t))*p3;
		    integrator.AppendCurvePoint(current,tangent);
		    if(draw)
		      ma.AddPointCurvePoint(current);
		    numpoints++;
		  }

		oldp = current;
	      }
	    tangent = p3-p2;
	    integrator.AppendCurvePoint(p3,tangent);
	    if(draw)
	      ma.AddPointCurvePoint(p3);
	    
	    (*testout) << "curve-spline with " << numpoints << " integration points" << endl;

	  }

	else if (inp == "end")
	  {
	    break;
	  }
	else
	  {
	    throw Exception("Error for integration along curve: something wrong in file\n");
	  }
	infile >> inp;
      }
  }


  void PDE :: LoadPDE (istream & input, const bool nomeshload, const bool nogeometryload)
  {
    pde = this;
    
    // Reset geometries 
    Ng_LoadGeometry("");
    
#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    hoprocs.SetSize(ntasks-1);
    for ( int i = 0; i<ntasks-1; i++ )
      hoprocs[i] = i+1;
#endif
    
    scan = new PDEScanner(&input);
    scan->ReadNext();
    CommandList(nomeshload,nogeometryload);
    delete scan;
  }


  void PDE :: LoadPDE (const string & filename, const bool nomeshload, const bool nogeometryload)
  {
    if (printmessage_importance>0)
      cout << "Load PDE from file " << filename << endl;
 

    string::size_type pos1 = filename.rfind('\\');
    string::size_type pos2 = filename.rfind('/');
    
    if (pos1 == filename.npos) pos1 = 0;
    if (pos2 == filename.npos) pos2 = 0;
    
    string pde_directory = filename.substr (0, max2(pos1, pos2));
    (*testout) << "pdefile " ;//<< pde->GetFilename() << endl;

    if(pde_directory == "")
      pde_directory = ".";

    if (printmessage_importance>0)
      cout << "dir = " << pde_directory << endl;

    pde = this;

#ifdef WIN32
	for(int i=0; pde_directory[i]!=0 && i<pde_directory.size(); i++)
		if(pde_directory[i] == '/')
			pde_directory[i] = '\\';
#endif
	

    pde->SetDirectory(pde_directory);
    pde->SetFilename(filename);
  
    ifstream infile (filename.c_str());
    if (!infile.good())
      {
	throw Exception (string ("PDE file " + filename + " not found"));
      }

      
    LoadPDE(infile,nomeshload,nogeometryload);

#ifdef PARALLEL
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    hoprocs.SetSize(ntasks-1);
    for ( int i = 0; i<ntasks-1; i++ )
      hoprocs[i] = i+1;
    
    if ( id == 0 )
       for ( int dest = 1; dest < ntasks; dest ++)
	  {
	    MyMPI_Send ("ngs_pdefile", dest );
	    MyMPI_Send(filename, dest);
	  }
#endif
  }


} // namespace

