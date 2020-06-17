#include <ngstd.hpp>
#include <nginterface.h>

#include <solve.hpp>
#include <parallelngs.hpp>


#ifdef HAVE_DLFCN_H 
#include <dlfcn.h>
#else
#include <windows.h>
#endif


#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
extern PythonEnvironment pyenv;
#endif



namespace ngcomp
{

  // parser for pde file

  enum TOKEN_TYPE 
    { 
      UNDEF = 0, NUMBER = 1, STRING = 2, END = 3,  
  
      PLUS = '+', MINUS = '-', 
      MUL = '*', DIV = '/',
      LP = '(', RP = ')', EQUAL = '=', COMMA = ',',
      LSB = '[', RSB = ']',
      COMMENT = '#',
      KW_DEFINE = 100, KW_GEOMETRY, KW_MESH, KW_SHARED, KW_PYMODULE,
      KW_CONSTANT, KW_VARIABLE, KW_FLAGS,
      KW_COEFFICIENT, KW_FESPACE, KW_GRIDFUNCTION, 
      KW_BILINEARFORM, KW_LINEARFORM, KW_PRECONDITIONER, 
      KW_INTEGRATOR = 200, KW_NUMPROC_ID,
      KW_NUMPROC = 300, KW_PYNUMPROC,
      KW_MATFILE,
      KW_LABEL = 1000, KW_GOTO, KW_IF, KW_STOP
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
      { KW_PYMODULE,    "pymodule" },
      { KW_CONSTANT,    "constant" },
      { KW_VARIABLE,    "variable" },
      { KW_FLAGS,       "flags" },
      { KW_COEFFICIENT, "coefficient" },
      { KW_FESPACE,     "fespace" },
      { KW_GRIDFUNCTION, "gridfunction" },
      { KW_BILINEARFORM, "bilinearform" },
      { KW_LINEARFORM,  "linearform" },
      { KW_PRECONDITIONER, "preconditioner" },
      { KW_NUMPROC,     "numproc" },
      { KW_PYNUMPROC,   "pynumproc" },
      { KW_MATFILE,     "matfile"},
      { KW_LABEL,       "label"},
      { KW_GOTO,        "goto"},
      { KW_STOP,        "stop"},
      { KW_IF,          "if"},
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
    streampos lastpos;

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
    void WriteBack()
    {
      scanin->seekg(lastpos);
    }

    void Error (const string & err);
  };


  PDEScanner :: PDEScanner (istream * ascanin)
  {
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
    for (auto bfi_register : itgs.GetBFIs())
      integrators.Set (bfi_register->name, 1);
    for (auto lfi_register : itgs.GetLFIs())
      integrators.Set (lfi_register->name, 1);

    NumProcs & nps = GetNumProcs();
    for (auto np_register : nps.GetNumProcs())
      numprocs.Set (np_register->name, 1);

    HandleStringConstants();
    scanin = new stringstream(copy_of_stream);
  }

  PDEScanner :: ~PDEScanner()
  {
    delete scanin;
  }

  // replace $(variable)
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

    int argc;
    char ** argv;
    Ng_GetArgs (argc, argv);
    
    for (int i = 0; i < argc; i++)
      {
	string flag = argv[i];
	if (flag.substr(0,2) == "-D")
	  {
	    int pos = flag.find ("=");
	    string varname = flag.substr (2, pos-2);
	    string value = flag.substr (pos+1, flag.length());
	    string fullname = string("$(")+varname+string(")");

	    size_t spos = 0;
	    while(spos != string::npos)
	      {
		spos = copy_of_stream.find(fullname,spos);

		if(spos != string::npos)
		  copy_of_stream.replace(spos,
					 fullname.length(),
					 value);
	      }
	  }
      }

    // cout << "new stream = " << copy_of_stream << endl;



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
        lastpos = scanin->tellg();    
  
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
	      string_value = "";
	      while ((*scanin) && ((isalnum (ch) 
                                    || ch == '_' || ch == '.'
                                    || ch == '?' || ch == '*'    // regex
                                    ) ))
		{
		  string_value += ch;
                  scanin->get(ch);
                  // if (scanin->eof()) break;
		}
	      scanin->putback (ch);

	      // (*scanin) >> string_value;
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


  class PDEEvalFunction : public EvalFunction
  {
    PDE & pde;
  public:
    PDEEvalFunction (PDE & apde)
      : pde(apde)
    {
      for (int i = 0; i < pde.GetConstantTable().Size(); i++)
        DefineConstant (pde.GetConstantTable().GetName(i),
                        pde.GetConstantTable()[i]);
      for (int i = 0; i < pde.GetVariableTable().Size(); i++)
        DefineGlobalVariable (pde.GetVariableTable().GetName(i),
                              pde.GetVariableTable()[i].get());
      for (int i = 0; i < pde.GenericVariables().Size(); i++)
        DefineGlobalVariable (pde.GenericVariables().GetName(i),
                              &pde.GenericVariables()[i]);
    }

    bool Parse (istream & aist, Array<shared_ptr<CoefficientFunction>> & depends)
    {
      for (int i = 0; i < pde.GetCoefficientTable().Size(); i++)
        DefineArgument (pde.GetCoefficientTable().GetName(i),-1,
                        pde.GetCoefficientTable()[i]->Dimension(),
                        pde.GetCoefficientTable()[i]->IsComplex());
      

      bool ok = EvalFunction::Parse(aist);
      
      int tot = 3;
      for (int i = 0; i < pde.GetCoefficientTable().Size(); i++)
        tot += pde.GetCoefficientTable()[i]->Dimension();
      Array<const char*> names(tot);
      names = NULL;
      for (int i = 0; i < arguments.Size(); i++)
        {
          int num = arguments[i].argnum;
          if (num >= 3)
            names[num] = arguments.GetName(i).c_str();
        }
      for (int i = 0; i < names.Size(); i++)
        if (names[i])
          {
            cout << "depend on " << names[i] << endl;
            depends.Append (pde.GetCoefficientFunction(names[i]));
          }
      return ok;
    }

  };


  void CommandList (bool nomeshload = false, const bool nogeometryload = false);
  void DefineCommand ();
  void NumProcCommand ();
  void CheckFlags (Flags & flags);
  Flags ParseFlags ();

  
  static PDEScanner * scan;
  NGS_DLL_HEADER shared_ptr<PDE> pde;  // valid douring pde-parsing
  
  
  void CommandList (bool nomeshload, const bool nogeometryload)
  {
    while (scan->GetToken() != END)
      {
	switch (scan->GetToken())
	  {

	  case KW_GEOMETRY:
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != '=') scan->Error ("Expected =");
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
	      if (scan->GetToken() != '=') scan->Error ("Expected '='");
	      scan->ReadNext();

	      string meshfile = pde->GetDirectory()+dirslash+scan->GetStringValue();
	      pde->SetMeshFileName(meshfile);
	      if(!nomeshload)
		{
                  if (ifstream (meshfile.c_str()))
		    {
		      cout << IM(1) << "Load mesh from file " << meshfile << endl;
                      auto ma = make_shared<MeshAccess>(meshfile);
                      // ma -> LoadMesh (meshfile);
                      pde -> AddMeshAccess(ma);
		    }
		  else
		    {
		      cout << IM(1) << "Load mesh from file " << scan->GetStringValue() << endl;
                      auto ma = make_shared<MeshAccess>(scan->GetStringValue());
                      // ma -> LoadMesh (scan->GetStringValue());
                      pde -> AddMeshAccess(ma);
		    }
		}

              // pde->GetMeshAccess().UpdateBuffers();
	      pde->AddVariable ("mesh.levels", pde->GetMeshAccess()->GetNLevels(), 6);
	      pde->AddVariable ("mesh.ne", pde->GetMeshAccess()->GetNE(), 6);
	      pde->AddVariable ("mesh.nv", pde->GetMeshAccess()->GetNV(), 6);
	      
	      // if (!pde->GetMeshAccess().GetNP())
	      if (pde->GetMeshAccess()->GetDimension() == -1)
		throw Exception ("No mesh or empty mesh file\n");
	      scan->ReadNext();
	      break;
	    }

	  case KW_SHARED:
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != '=') scan->Error ("Expected '='");
	      scan->ReadNext();
		  
              string shared = scan->GetStringValue();
	      scan->ReadNext();

              
#ifdef HAVE_DLFCN_H
#ifdef __APPLE__
	      shared += ".dylib";
#else
	      shared += ".so";
#endif
              cout << IM(1) << "load shared library '" << shared << "'" << endl;

              void * handle = dlopen (shared.c_str(), RTLD_LAZY | RTLD_GLOBAL);
              if (!handle)
                {
                  stringstream err;
                  err << "Cannot load shared library '" << shared << "' \nerrmsg: "  << dlerror();
                  throw Exception (err.str());
                }
#else
	      shared += ".dll";
              cout << IM(1) << "load shared library '" << shared << "'" << endl;
	      
	      HINSTANCE handle = LoadLibrary (shared.c_str());
	      if (!handle)
		{
		  stringstream err;
                  err << "Cannot load shared library '" << shared << "' \nerrmsg: "; //   << dlerror();
                  throw Exception (err.str());
		}
#endif
              break;
            }

	  case KW_PYMODULE:
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != '=') scan->Error ("Expected '='");
	      scan->ReadNext();
		  
              string module_name = scan->GetStringValue();
	      scan->ReadNext();
#ifdef NGS_PYTHON
              {
                cout << "*********** load py-module " << module_name << endl;
                AcquireGIL gil_lock;
                string command = string("from ") + module_name + " import*";
                pyenv.exec (command.c_str());
                // pyenv.exec("from module1 import *");
              }

#else
              cout << "want to load py-module, but python not enabled" << endl;
#endif
              break;
            }
	    
	  case KW_MATFILE:
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != '=') scan->Error ("Expected '='");
	      scan->ReadNext();
	      
	      string matfile  = scan->GetStringValue();
	      if (ifstream (matfile.c_str()).good())
		{
		  cout << IM(1) << "Materialdata from file " << matfile << endl;
		  pde->SetMatfile(matfile); 
		}
	      else
		{
		  cout << IM(1) << "Materialdata from file " << matfile << endl;
		  throw Exception("**** Materialdata-File not found!! \n ");  
		}
	      
	      scan->ReadNext();
	      
	      break;
	    } 
	    

	  case KW_DEFINE:
	    {
              scan->ReadNext();
	      DefineCommand ();
	      break;
	    }

          case KW_CONSTANT:
          case KW_VARIABLE:
          case KW_FLAGS:
          case KW_COEFFICIENT:
          case KW_FESPACE:
          case KW_GRIDFUNCTION:
          case KW_BILINEARFORM:
          case KW_LINEARFORM:
          case KW_PRECONDITIONER:
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
	      flags.SetFlag ("name", name.c_str());

	      int dim = pde->GetMeshAccess()->GetDimension();
	      if (GetNumProcs().GetNumProc(npid, dim))
		{
		  pde -> AddNumProc (name, GetNumProcs().GetNumProc(npid, dim)->creator(pde, flags));
// #ifdef SOCKETS
// 		  if(pde -> ConstantUsed ("clientserver") && pde -> GetConstant("clientserver") > 0.5)
// 		    pde -> GetClientSocketAccess().CheckNumProc(name,flags,pde->GetNumProc(name)->GetCallPosition());
// #endif
		}
	      else
		throw Exception (string("Undefined numproc ") + npid);
	      break;
	    }

	  case KW_PYNUMPROC:
            {
	      scan->ReadNext();
	      string npid = scan->GetStringValue();

	      scan->ReadNext();
	      string name = scan->GetStringValue();

	      scan->ReadNext();
	      Flags flags;
	      CheckFlags (flags);
	      flags.SetFlag ("name", name.c_str());
              
#ifdef NGS_PYTHON
              {
                cout << "*********** create py-numproc " << npid << ", name = " << name << endl;
                AcquireGIL gil_lock;

                static int pynp_cnt = 0;
                string pyname = "_pynumproc"+ToString(pynp_cnt++);

		try{
                    // get the python class
                    py::object npclass = pyenv[npid.c_str()];
                    // call constructor and release the python instance to avoid reference counting on python side
                    py::handle np = npclass(py::cast(pde),py::cast(flags)).release();
                    pde->AddNumProc (name, np.cast<shared_ptr<NumProc>>());
                }
		catch (py::error_already_set const &e) {
                    cerr << e.what() << endl;
                    PyErr_Print();
		}
              }
#else
              cerr << "sorry, python not enabled" << endl;
#endif
              break;
            }

          case KW_LABEL:
            {
              scan -> ReadNext();
              string label = scan->GetStringValue();
              scan -> ReadNext();
              pde -> AddControlStatement (make_shared<LabelStatement> (pde->GetMeshAccess(), "dummy", label));
              break;
            }

          case KW_GOTO:
            {
              scan -> ReadNext();
              string label = scan->GetStringValue();
              scan -> ReadNext();
              pde -> AddControlStatement (make_shared<GotoStatement> (pde->GetMeshAccess(), "dummy", label));
              break;
            }

          case KW_IF:
            {
              EvalFunction * fun = new PDEEvalFunction(*pde);
              fun->Parse (*scan->scanin);
              
              scan -> ReadNext();
              if (scan->GetToken() != KW_GOTO) cerr << "goto expected" << endl;

              scan -> ReadNext();
              string label = scan->GetStringValue();
              scan -> ReadNext();

              pde -> AddControlStatement (make_shared<ConditionalGotoStatement> (pde->GetMeshAccess(), "dummy", fun, label));              
              break;
            }

          case KW_STOP:
            {
              scan -> ReadNext();
              pde -> AddControlStatement (make_shared<StopStatement> (pde->GetMeshAccess(), "dummy"));
              break;
            }

	  default:
	    {
              if (scan -> GetToken() == STRING &&
                  pde -> GetVariableTable().Used (scan->GetStringValue()))
                {
                  cout << "add EvalFunction" << endl;
                  double & var = pde->GetVariable(scan->GetStringValue());
                  scan -> ReadNext(); // '='

                  EvalFunction * fun = new PDEEvalFunction (*pde);
                  fun -> Parse (*scan->scanin);
                  /*
                  for (int i = 0; i < pde->GetConstantTable().Size(); i++)
                    ev->GetEvaluator().DefineConstant (pde->GetConstantTable().GetName(i),
                                                        pde->GetConstantTable()[i]);
                  for (int i = 0; i < pde->GetVariableTable().Size(); i++)
                    ev->GetEvaluator().DefineGlobalVariable (pde->GetVariableTable().GetName(i),
                                                              pde->GetVariableTable()[i]);
                  ev -> GetEvaluator().Parse (*scan->scanin);
                  */
                  
                  auto ev = make_shared<EvalVariable> (pde->GetMeshAccess(), scan->GetStringValue(), fun);
                  ev -> SetVariable (var);
                  pde -> AddVariableEvaluation (ev);
                  scan->ReadNext();
                }
              
              else if (scan -> GetToken() == STRING &&
                       pde -> GetGridFunctionTable().Used (scan->GetStringValue()))
                {
                  cout << "add GridFunction assignment" << endl;
                  string gfname = scan->GetStringValue();
                  shared_ptr<GridFunction> gf = pde->GetGridFunction(gfname);
                  scan -> ReadNext(); // '='

                  PDEEvalFunction * fun = new PDEEvalFunction(*pde);
                  Array<shared_ptr<CoefficientFunction>> depends;
                  fun->Parse (*scan->scanin, depends);
                  scan -> ReadNext();
                  cout << "set gridfunction " << gf -> GetName() << " to evalfunction ";

                  static int cnt = 0;
                  cnt++;
                  string coefname = "assigncoef" + ToString(cnt);
                  pde->AddCoefficientFunction
                    (coefname, make_shared<DomainVariableCoefficientFunction>(*fun, depends));
                  Flags flags = ParseFlags();
                  flags.SetFlag ("gridfunction", gfname.c_str());
                  flags.SetFlag ("coefficient", coefname.c_str());

                  string npname = "assignnp" + ToString(cnt);
                  pde -> AddNumProc (npname, GetNumProcs().GetNumProc("setvalues", pde->GetMeshAccess()->GetDimension())
                                     -> creator (pde, flags));
                }

              
              else
                {
                  stringstream errstr;
                  errstr << "Unknown command";
                  if (scan->GetToken() == STRING)
                    errstr << ": " 
                           << "\"" << scan->GetStringValue() << "\"" << endl;
                  else
                    errstr << ", token = " << scan->GetToken() 
                           << ", char = " << char(scan->GetToken()) << endl;
                  scan -> Error (errstr.str());
                  
                  scan->ReadNext();
                }
	    }
	  }
      }


    cout << IM(1) << "End of file reached" << endl << endl << endl;
  }


		      
  



  void DefineCommand ()
  {
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

	  if (scan->GetToken() == LP ||
              scan->GetToken() == NUMBER ||
              scan->GetToken() == MINUS)
	    {
	      EvalFunction * fun = new PDEEvalFunction (*pde);
              scan->WriteBack();
	      fun->Parse (*scan->scanin);

	      if (fun->IsConstant())
		val = fun->EvalConstant();
	      else
		scan->Error  ("Expression not constant");
              delete fun;
	    }
          /*
	  else if (scan->GetToken() == NUMBER)
	    val = scan->GetNumValue();
	  else if (scan->GetToken() == MINUS)
	    {
	      scan->ReadNext();
	      if (scan->GetToken() != NUMBER)
		scan -> Error ("expected a number");
	      val = -scan->GetNumValue();
	    }
          */
	  else if (scan->GetToken() == STRING)
	    {
	      sval = scan->GetStringValue();
	    }
	  else
	    scan->Error ("Expected a number");

	  if(sval == "")
	    {
	      pde->AddConstant (name, val);
	      pde->AddCoefficientFunction
		(name, shared_ptr<CoefficientFunction> (new ConstantCoefficientFunction(val)));
	    }
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

	  // double val = 0;
	  
          EvalFunction * fun = new PDEEvalFunction(*pde);
          scan->WriteBack();
          fun->Parse (*scan->scanin);

          // cout << "function parsing complete, fun = ";
          // fun->Print(cout);
          // cout << endl;

          if (! fun -> IsResultComplex() && fun->Dimension() == 1)
            {
              if (fun->IsConstant())
                pde->AddVariable (name, fun->EvalConstant(), 1);
              else
                pde->AddVariable (name, make_shared<EvalVariable> (pde->GetMeshAccess(), name, fun));
            }
          else
            {
              GenericVariable var (fun->IsResultComplex(), fun->Dimension());
              if (var.IsComplex())
                fun -> Eval ( (Complex*)NULL, &var.ValueDouble(), var.Dimension());
              else
                fun -> Eval ( (double*)NULL, &var.ValueDouble(), var.Dimension());
              
              cout << "generic variable " << name << " = " << var << endl;
              pde->GenericVariables().Set(name, var);
            }
	  scan->ReadNext();

          /*
	  EvalVariable * eval = NULL;

	  if (scan->GetToken() == LP)
	    {
	      eval = new EvalVariable(pde->GetMeshAccess(),name);
	      for (int i = 0; i < pde->GetConstantTable().Size(); i++)
		eval->GetEvaluator().DefineConstant (pde->GetConstantTable().GetName(i),
						     pde->GetConstantTable()[i]);
	      for (int i = 0; i < pde->GetVariableTable().Size(); i++)
		eval->GetEvaluator().DefineGlobalVariable (pde->GetVariableTable().GetName(i),
							   pde->GetVariableTable()[i]);
              scan->WriteBack();
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
	    pde->AddVariable (name, val, 1);
	  scan->ReadNext();
          */
	  break;
	}


      case KW_FLAGS:
	{
	  scan->ReadNext();

	  string name = scan->GetStringValue ();

	  scan->ReadNext();
	  if (scan->GetToken() != '=')
	    scan->Error ("Expected =");

	  scan->ReadNext();

          Flags flags;
          CheckFlags (flags);

          pde->AddFlags (name, flags);
          break;
        }

      case KW_COEFFICIENT:
	{
	  Array<double> dcoeffs;
	  Array<shared_ptr<EvalFunction>> coeffs;
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
		      EvalFunction * fun = new PDEEvalFunction (*pde);
                      /*
		      for (int i = 0; i < pde->GetConstantTable().Size(); i++)
                       fun->DefineConstant (pde->GetConstantTable().GetName(i),
					     pde->GetConstantTable()[i]);
		      for (int i = 0; i < pde->GetVariableTable().Size(); i++)
			fun->DefineGlobalVariable (pde->GetVariableTable().GetName(i),
						   pde->GetVariableTable()[i]);
                      */
		      
                      scan->WriteBack();
		      fun->Parse (*scan->scanin);
		      funs.Set (mat,fun);
		      *testout << "set material " << mat << " to evalfunc" << endl;
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
		      *testout << "set material " << mat << " to constant evalfunc" << endl;
		    }
		  scan->ReadNext();

		  if (scan->GetToken() == COMMA) 
		    scan->ReadNext();
		  //values.Set (mat, val);
		}

	      //cout << "values = " << endl << values << endl;

	      int maxdom = -1;
	      int ne = pde->GetMeshAccess()->GetNE();
	      for (int i = 0; i < ne; i++)
		maxdom = max2 (maxdom, pde->GetMeshAccess()->GetElIndex(ElementId(VOL,i)));
	      maxdom++;

	      dcoeffs.SetSize(maxdom);
	      dcoeffs = 0;
	      coeffs.SetSize(maxdom);
	      coeffs = nullptr;
	      bool only_constant = true;
	      for (int i = 0; i < ne; i++)
		{
		  int index = pde->GetMeshAccess()->GetElIndex(ElementId(VOL,i));
		  if (coeffs[index]) continue;

		  string mat = pde->GetMeshAccess()->GetMaterial(ElementId(VOL,i));
		  // cout << "mat = " << mat << ", ind = " << index << endl;

		  EvalFunction * fun = NULL;
		  bool used = false;

		  for(int j = 0; j < funs.Size(); j++)
		    {
		      // *testout << "check pattern, mat = '" << mat << "', name = '" << funs.GetName(j) << "'" << endl;
		      used = StringFitsPattern(mat,funs.GetName(j));
		      if(used)
			{
			  //(*testout) << "\"" << mat << "\" fits \"" << funs.GetName(j) << "\"" << endl;
			  fun = funs[j];
			  // *testout << "pattern match, mat = " << mat << " fun.name = " << funs.GetName(j) << endl;
			  break;
			}
		    }
		  if(!used)
		    {
		      if (funs.Used ("default"))
			{
			  fun = funs["default"];
			  // *testout << "use default value for mat " << mat << endl;
			}
		      else
			throw Exception (string ("No value defined for material ")+mat);
		    }
		  
			 
		  coeffs[index] = make_shared<EvalFunction> (*fun);
		  if(fun->IsConstant())
		    dcoeffs[index] = fun->Eval( (double*)(0));
		  else
		    only_constant = false;
		}

	      // iff not all indices are in use ...
	      shared_ptr<EvalFunction> some_fun;
	      for (int i = 0; i < coeffs.Size(); i++)
		if (coeffs[i]) some_fun = coeffs[i];
	      for (int i = 0; i < coeffs.Size(); i++)
		if (!coeffs[i])
		  coeffs[i] = some_fun; // new EvalFunction (*some_fun);

	      for(int i=0; i<funs.Size(); i++)
		delete funs[i];

	      if(only_constant)
		{
		  (*testout) << "material coefficients = " << endl << dcoeffs << endl;
		  pde->AddCoefficientFunction
		    (name, shared_ptr<CoefficientFunction> (new DomainConstantCoefficientFunction(dcoeffs)));
		}
	      else
		{
		  (*testout) << "material coefficients variable " << endl;
                  pde->AddCoefficientFunction
                    (name, make_shared<DomainVariableCoefficientFunction>(coeffs));
		}

              /*
	      for (int hi = 0; hi < coeffs.Size(); hi++)
		delete coeffs[hi];
              */
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
		(name, shared_ptr<CoefficientFunction> 
                 (new FileCoefficientFunction(ipfilename, infofilename,
                                              valuesfilename, loadvalues)));
	      
	      //scan -> ReadNext();
	    }

	  else if (strcmp (scan->GetStringValueC(), "bcnames") == 0)

	    {
	      // define values by names of boundary conditions

	      scan->ReadNext();
	      //SymbolTable<double> values;
	      SymbolTable<shared_ptr<EvalFunction>> funs;
	      while (scan->GetToken() == STRING)
		{
		  string bcname = scan->GetStringValue();
		  scan->ReadNext();

		  //double val = 1;

		  if (scan->GetToken() == LP)
		    {
		      auto fun = make_shared<EvalFunction> ();
		      for (int i = 0; i < pde->GetConstantTable().Size(); i++)
			fun->DefineConstant (pde->GetConstantTable().GetName(i),
					     pde->GetConstantTable()[i]);
		      for (int i = 0; i < pde->GetVariableTable().Size(); i++)
			fun->DefineGlobalVariable (pde->GetVariableTable().GetName(i),
						   pde->GetVariableTable()[i].get());
		      
		      
                      scan->WriteBack();
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

		      auto fun = make_shared<EvalFunction> ();
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
	      int nse = pde->GetMeshAccess()->GetNSE();
	      for (int i = 0; i < nse; i++)
		maxbc = max2 (maxbc, pde->GetMeshAccess()->GetElIndex(ElementId(BND,i)));
	      maxbc++;

	      dcoeffs.SetSize(maxbc);
	      dcoeffs = 0;
	      coeffs.SetSize(maxbc);
	      coeffs = nullptr;
	      bool only_constant = true;
	      for (int i = 0; i < nse; i++)
		{
		  int index = pde->GetMeshAccess()->GetElIndex(ElementId(BND,i));
		  if (coeffs[index]) continue;

		  shared_ptr<EvalFunction> fun = NULL;
		  string bcname = pde->GetMeshAccess()->GetMaterial(ElementId(BND, i));
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
			 
		  if(coeffs[index] == NULL) coeffs[index] = fun; // new EvalFunction(*fun);
		  if(fun->IsConstant())
		    dcoeffs[index] = fun->Eval( (double*)(0) );
		  else
		    only_constant = false;
		}

	      // for(int i=0; i<funs.Size(); i++)
              // delete funs[i];

	      if(only_constant)
		{
		  (*testout) << "material coefficients = " << endl << dcoeffs << endl;
		  pde->AddCoefficientFunction
		    (name, shared_ptr<CoefficientFunction> (new DomainConstantCoefficientFunction(dcoeffs)));
		}
	      else
		{
		  (*testout) << "material coefficients variable " << endl;
                  pde->AddCoefficientFunction
                    (name, make_shared<DomainVariableCoefficientFunction>(coeffs));
		}
	      
	      // for (int hi = 0; hi < coeffs.Size(); hi++)
              // delete coeffs[hi];
	    }

	  else
	    
	    {
	      while (scan->GetToken() == NUMBER || 
		     scan->GetToken() == MINUS ||
		     scan->GetToken() == LP ||
		     scan->GetToken() == LSB ||
                     (scan->GetToken() == STRING &&
                      (
                       pde->GetConstantTable().Used(scan->GetStringValue()) ||
                       pde->GetVariableTable().Used(scan->GetStringValue()) ||
                       pde->GetCoefficientTable().Used(scan->GetStringValue())
                       )
                      )
                     )
		{
		  // if (scan->GetToken() == LP)
                  if (scan->GetToken() != LSB)
		    {
		      auto fun = make_shared<EvalFunction> ();
		      for (int i = 0; i < pde->GetConstantTable().Size(); i++)
			fun->DefineConstant (pde->GetConstantTable().GetName(i),
					     pde->GetConstantTable()[i]);
		      for (int i = 0; i < pde->GetVariableTable().Size(); i++)
			fun->DefineGlobalVariable (pde->GetVariableTable().GetName(i),
						   pde->GetVariableTable()[i].get());
		      for (int i = 0; i < pde->GetCoefficientTable().Size(); i++)
			fun->DefineArgument (pde->GetCoefficientTable().GetName(i),-1,
					     pde->GetCoefficientTable()[i]->Dimension(),
					     pde->GetCoefficientTable()[i]->IsComplex());
                      
                      scan->WriteBack();
		      fun->Parse (*scan->scanin);
		      coeffs.Append (fun);
		      if (fun->IsConstant())
			dcoeffs.Append (fun->Eval ( (double*)(0) ));
		      scan->ReadNext();
		    }
		  else //  Token() == LSB   polynomial [MW]
		    {
		      cout << IM(1) << "polynomial: ";
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
			  cout << IM(1) << (*polyc)[0];
                          
			  for(int i = 1; i<polyc->Size(); i++)
			    cout << IM(1) << " + " << (*polyc)[i] << "*t^"<<i;
			  

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
			      cout << IM(1) << " until " << val <<", then ";
			      scan->ReadNext(); // RP
			      scan->ReadNext();
			    }
			}
		      polybounds.Append(polyb);
		      polycoeffs.Append(polyco);

		      cout << IM(1) << endl;

		      scan->ReadNext();
		    }
                  /*
		  else
		    {  // now combined with parser
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
                  */
		  if (scan->GetToken() == COMMA) 
		    scan->ReadNext();
		}

	      bool allconst = true;
	      for (int i = 0; i < coeffs.Size(); i++)
		if (!coeffs[i]->IsConstant())
		  allconst = false;
	      
	      if (polycoeffs.Size() > 0)
		{
		  pde->AddCoefficientFunction
		    (name, shared_ptr<CoefficientFunction> (new PolynomialCoefficientFunction(polycoeffs,polybounds)));
		}
	      else if (allconst)
		{
		  if (dcoeffs.Size() == 1)
		    pde->AddCoefficientFunction
		      (name, shared_ptr<CoefficientFunction> (new ConstantCoefficientFunction(dcoeffs[0])));
		  else
		    pde->AddCoefficientFunction
		      (name, shared_ptr<CoefficientFunction> (new DomainConstantCoefficientFunction(dcoeffs)));
		  // for (int hi = 0; hi < coeffs.Size(); hi++)
                  // delete coeffs[hi];
		}
	      else
		{
		  Array<shared_ptr<CoefficientFunction>> depends; // pde->GetCoefficientTable().Size()+3);
		  if (coeffs[0])
		    {
		      int tot = 3;
		      for (int i = 0; i < pde->GetCoefficientTable().Size(); i++)
			tot += pde->GetCoefficientTable()[i]->Dimension();
		      Array<const char*> names(tot);
		      names = NULL;
		      for (int i = 0; i < coeffs[0]->arguments.Size(); i++)
			{
			  int num = coeffs[0]->arguments[i].argnum;
			  if (num >= 3)
			    names[num] = coeffs[0]->arguments.GetName(i).c_str();
			}
		      for (int i = 0; i < names.Size(); i++)
			if (names[i])
			  {
			    cout << "depend on " << names[i] << endl;
			    depends.Append (pde->GetCoefficientFunction(names[i]));
			  }
		    }

                  pde->AddCoefficientFunction
                    (name, make_shared<DomainVariableCoefficientFunction>(coeffs, depends));
                  
		  // for (int hi = 0; hi < coeffs.Size(); hi++)
                  // delete coeffs[hi];
		}
	    }

	  break;
	}

      case KW_FESPACE:
	{
	  scan->ReadNext();

	  string name = scan->GetStringValue ();
	  scan->ReadNext();

	  pde->AddFESpace (name, ParseFlags());
	  break;
	}

      case KW_GRIDFUNCTION:
	{
	  scan->ReadNext();
	  string name = scan->GetStringValue ();
	  scan->ReadNext();

	  pde->AddGridFunction (name, ParseFlags());
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
	      auto source = pde->GetBilinearForm(flags.GetStringFlag("useintegratorsof",""));
	      for(int i = 0; i < source->NumIntegrators(); i++)
		pde->AddBilinearFormIntegrator(name,source->GetIntegrator(i),false);
	    
	      // for(int i = 0; i < source->NumIndependentIntegrators(); i++)
	      // pde->AddBilinearFormIntegrator(name,source->GetIndependentIntegrator(i),false);
	    }

	  //
	
	  // read bilinear-form components
	  while (1)
	    {
	      TOKEN_TYPE integrator_token = scan->GetToken();

	      if (integrator_token == KW_INTEGRATOR)
		{
		  string integrator_name = scan->GetStringValue();
		  //const ngfem::Integrators::IntegratorInfo
                  auto info = ngfem::GetIntegrators() 
		    . GetBFI(integrator_name, pde->GetMeshAccess()->GetDimension());

		  if (info)
		    {
		      Array<shared_ptr<CoefficientFunction>> coeffs(info->numcoeffs);
		      scan->ReadNext();
		      for (int i = 0; i < info->numcoeffs; i++)
			{
			  coeffs[i] = pde->GetCoefficientFunction (scan->GetStringValue(), 1);
			  if (!coeffs[i])
			    {
                              shared_ptr<GridFunction> gf = pde->GetGridFunction (scan->GetStringValue(), 1);
			      if (gf) coeffs[i] = make_shared<GridFunctionCoefficientFunction> (gf);
			    }

                          if (!coeffs[i])
                            {
                              EvalFunction * fun = new PDEEvalFunction(*pde);
                              scan-> WriteBack();
                              if (fun->Parse(*scan->scanin))
                                {
                                  if (fun->IsConstant())
                                    coeffs[i] = shared_ptr<CoefficientFunction> (new ConstantCoefficientFunction (fun->EvalConstant()));
                                  else
                                    coeffs[i] = make_shared<DomainVariableCoefficientFunction> (*fun);
                                }
                              delete fun;
                            }

                          if (!coeffs[i])
                            throw Exception (string("undefined coefficient ") + scan->GetStringValue());

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
                      
                      /*
		      ngfem::BilinearFormIntegrator * integrator = 
			dynamic_cast<ngfem::BilinearFormIntegrator*> (info->creator(coeffs));
                      */
                      auto integrator = info->creator(coeffs);
		      integrator -> SetName (integrator_name);
		      integrator -> SetFlags (partflags);



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
			      (pde->GetBilinearForm (name)->GetFESpace().get()))
			    {
			      integrator = make_shared<CompoundBilinearFormIntegrator>
                                (integrator, int(partflags.GetNumFlag ("comp", 1))-1);
                        }
			  else
			    {
			      integrator = make_shared<BlockBilinearFormIntegrator>
                                (integrator, 
                                 pde->GetBilinearForm (name)->GetFESpace()->GetDimension(),
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
			  integrator = make_shared<ComplexBilinearFormIntegrator> (integrator, Complex(0,1));
			}

		      if (partflags.GetDefineFlag ("transpose"))
			{
			  integrator = make_shared<TransposeBilinearFormIntegrator> (integrator);
			}

		      if (partflags.GetDefineFlag ("real"))
			{
			  integrator = make_shared<ComplexBilinearFormIntegrator> (integrator, Complex(1,0));
			}


		      if (partflags.NumFlagDefined ("definedon") || partflags.NumListFlagDefined("definedon"))
			{
			  int size = max2 (pde->GetMeshAccess()->GetNDomains(), 
					   pde->GetMeshAccess()->GetNBoundaries());
			  BitArray definedon(size);
			  definedon.Clear();

			  if(partflags.NumFlagDefined ("definedon"))
			    {
			      int defon = int(partflags.GetNumFlag("definedon", 0))-1;
			      if (defon < size) definedon.SetBit(defon);
			    }
			  if(partflags.NumListFlagDefined("definedon"))
			    for(int i=0; i<partflags.GetNumListFlag("definedon").Size(); i++)
			      {
				int defon = int(partflags.GetNumListFlag("definedon")[i])-1;
				if (defon < size) definedon.SetBit(defon);
			      }
			  integrator->SetDefinedOn (definedon);
			}


		      integrator -> SetConstantCoefficient (partflags.GetDefineFlag("const"));

		      int numregions = integrator -> BoundaryForm() ? 
			pde->GetMeshAccess()->GetNBoundaries() : pde->GetMeshAccess()->GetNDomains();
		      while (numregions > 0 && !integrator->DefinedOn(numregions-1))
			numregions--;

		      for (int i = 0; i < coeffs.Size(); i++)
			{
			  if (coeffs[i] -> NumRegions () < numregions)
			    {
			      stringstream str;
			      str << "coefficient " << i 
				  << " of integrator " << integrator_name << " has " 
				  << coeffs[i]->NumRegions() << " regions, but should be " << numregions;
			      
			      throw Exception (str.str());
			    }
			}




// 		      bool indep(false);
// 		      if (partflags.NumFlagDefined ("master") && partflags.NumFlagDefined ("other"))
// 			{
// 			  indep = true;
// 			  pde->AddIndependentBilinearFormIntegrator(name, integrator, 
// 								    static_cast<int>(partflags.GetNumFlag("master",0))-1,
// 								    static_cast<int>(partflags.GetNumFlag("other",0))-1);
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
	  while (1)
	    {
	      TOKEN_TYPE integrator_token = scan->GetToken();

	      if (integrator_token == KW_INTEGRATOR)
		{
		  string integrator_name = scan->GetStringValue();

		  auto info =
		    ngfem::GetIntegrators().GetLFI(integrator_name, pde->GetMeshAccess()->GetDimension());
		
		  if (info)
		    {
		      Array<shared_ptr<CoefficientFunction>> coeffs(info->numcoeffs);
		      scan->ReadNext();
		      for (int i = 0; i < info->numcoeffs; i++)
			{
			  coeffs[i] = pde->GetCoefficientFunction (scan->GetStringValue(),true);
			  if (!coeffs[i])
			    {
			      shared_ptr<GridFunction> gf = pde->GetGridFunction (scan->GetStringValue(), 1);
			      if (gf)
				coeffs[i] = make_shared<GridFunctionCoefficientFunction> (gf);
			    }

                          if (!coeffs[i])
                            {
                              EvalFunction * fun = new PDEEvalFunction(*pde);
                              scan-> WriteBack();
                              if (fun->Parse(*scan->scanin))
                                {
                                  if (fun->IsConstant())
                                    coeffs[i] = make_shared<ConstantCoefficientFunction> (fun->EvalConstant());
                                  else
                                    coeffs[i] = make_shared<DomainVariableCoefficientFunction> (*fun);
                                }
                              delete fun;
                            }

                          if (!coeffs[i])
                            throw Exception (string("undefined coefficient ") + scan->GetStringValue());

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

                      /*
		      ngfem:: integrator = 
			dynamic_cast<ngfem::LinearFormIntegrator*> (info->creator(coeffs));
                      */
                      auto integrator = info->creator(coeffs);
		      int numregions = integrator -> BoundaryForm() ? 
			pde->GetMeshAccess()->GetNBoundaries() : pde->GetMeshAccess()->GetNDomains();
		      for (int i = 0; i < coeffs.Size(); i++)
			{
			  if (coeffs[i] -> NumRegions () < numregions)
			    {
			      stringstream str;
			      str << "coefficient " << i 
				  << " of integrator " << integrator_name << " has " 
				  << coeffs[i]->NumRegions() << " regions, but should be " << numregions;
			      
			      throw Exception (str.str());
			    }
			}


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
			      (pde->GetLinearForm (name)->GetFESpace().get()))
			    {
			      integrator = make_shared<CompoundLinearFormIntegrator>
                                (integrator, int(partflags.GetNumFlag ("comp", 1))-1);
			    }
			  else
			    {
			      integrator = make_shared<BlockLinearFormIntegrator>
                                (integrator, pde->GetLinearForm (name)->GetFESpace()->GetDimension(),
                                 int(partflags.GetNumFlag ("comp", 1))-1);
                        }
			}

		      if (partflags.NumFlagDefined ("normal"))
			{
			  throw Exception ("-normal flag currently not available");
			}
		    
		      if (partflags.GetDefineFlag ("imag"))
			{
			  integrator = shared_ptr<LinearFormIntegrator> (new ComplexLinearFormIntegrator
                                                                         (integrator, Complex(0,1)));
			} 

		      if (partflags.GetDefineFlag ("real"))
			{
			  integrator = shared_ptr<LinearFormIntegrator> (new ComplexLinearFormIntegrator (integrator, Complex(1,0)));
			} 

		      if (partflags.NumFlagDefined ("definedon") || partflags.NumListFlagDefined("definedon"))
			{
			  int size = max2 (pde->GetMeshAccess()->GetNDomains(), 
					   pde->GetMeshAccess()->GetNBoundaries());
			  BitArray definedon(size);
			  definedon.Clear();

			  if(partflags.NumFlagDefined ("definedon"))
			    {
			      int defon = int(partflags.GetNumFlag("definedon", 0))-1;
			      if (defon < size) definedon.SetBit(defon);
			    }
			  if(partflags.NumListFlagDefined("definedon"))
			    for(int i=0; i<partflags.GetNumListFlag("definedon").Size(); i++)
			      {
				int defon = int(partflags.GetNumListFlag("definedon")[i])-1;
				if (defon < size) definedon.SetBit(defon);
			      }

			  integrator->SetDefinedOn (definedon);
			}

		      if (partflags.NumFlagDefined ("cachecomp"))
			{
			  integrator->SetCacheComp(int(partflags.GetNumFlag("cachecomp",0)));
			}
		    
		      if (partflags.StringFlagDefined ("curvefile"))
			{		
			  pde->SetLineIntegratorCurvePointInfo(pde->GetDirectory() + dirslash + partflags.GetStringFlag("curvefile",""),
			  				       integrator.get());

			  //BuildLineIntegratorCurvePoints ( pde->GetDirectory() + dirslash + partflags.GetStringFlag("curvefile",""),
			  //				   pde->GetMeshAccess(),
			  //				   *integrator);
			}

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
        /*
	scan->ReadNext();
	string flag = string("-") + scan->GetStringValue();
	flags.SetCommandLineFlag (flag.c_str());
	scan->ReadNext();
        */
        scan -> WriteBack();
        string str;
        *(scan->scanin) >> str;
	flags.SetCommandLineFlag (str.c_str(), &pde->GetFlagsTable());
	scan->ReadNext();
      }
  }


  Flags ParseFlags()
  {
    Flags flags;
    CheckFlags (flags);
    return flags;
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
		    verts = ma.GetElVertices(ElementId(VOL,element));
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
		    verts = ma.GetElVertices(ElementId(VOL,element));
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


  void LoadPDE (shared_ptr<PDE> apde, istream & input, const bool nomeshload, const bool nogeometryload)
  {
    pde = apde;
    
    // Reset geometries 
    Ng_LoadGeometry("");
    
    scan = new PDEScanner(&input);
    scan->ReadNext();
    CommandList(nomeshload,nogeometryload);
    delete scan;
    pde = nullptr;
  }
  

  shared_ptr<PDE> LoadPDE (istream & input, const bool nomeshload, const bool nogeometryload)
  {
    shared_ptr<PDE> apde = make_shared<PDE>();
    LoadPDE (apde, input, nomeshload, nogeometryload);
    return apde;
  }


  void LoadPDE (shared_ptr<PDE> apde, const string & filename, 
                const bool nomeshload, const bool nogeometryload)
  {
    static Timer timer("LoadPDE");
    RegionTimer reg (timer);

    cout << IM(1) << "Load PDE from file " << filename << endl;
    string data;
    pde = apde;
    NgMPI_Comm ngs_comm; //  = MPI_COMM_WORLD;
    if (ngs_comm.Rank() == 0)
      {
	string::size_type pos1 = filename.rfind('\\');
	string::size_type pos2 = filename.rfind('/');
	
	if (pos1 == filename.npos) pos1 = 0;
	if (pos2 == filename.npos) pos2 = 0;
	
	string pde_directory = filename.substr (0, max2(pos1, pos2));
	(*testout) << "pdefile " ;//<< pde->GetFilename() << endl;
	
	if(pde_directory == "")
	  pde_directory = ".";
	
	cout << IM(1) << "dir = " << pde_directory << endl;
#ifdef WIN32
	for(int i=0; pde_directory[i]!=0 && i<pde_directory.size(); i++)
	  if(pde_directory[i] == '/')
	    pde_directory[i] = '\\';
#endif
	pde->SetDirectory(pde_directory);
	pde->SetFilename(filename);

	ifstream infile (filename.c_str());
	if (!infile.good())
	  throw Exception (string ("PDE file " + filename + " not found"));
	
	while (!infile.eof())
	  {
	    char ch;
	    infile.get(ch);
	    data += ch;
	  }

	string hfilename = filename;
	ngs_comm.Bcast (hfilename);
	ngs_comm.Bcast (pde_directory);
      }

    else

      {
	string filename, pde_directory;

	ngs_comm.Bcast (filename);
	ngs_comm.Bcast (pde_directory);
	pde->SetDirectory(pde_directory);
	pde->SetFilename(filename);
      }
    
    // MyMPI_Bcast (data, ngs_comm);
    ngs_comm.Bcast (data);

    stringstream strdata(data);
    LoadPDE(pde, strdata, nomeshload, nogeometryload);

    pde = nullptr;
  }


  shared_ptr<PDE> LoadPDE (const string & filename, 
                           const bool nomeshload, const bool nogeometryload)
  {
    shared_ptr<PDE> apde = make_shared<PDE>();
    LoadPDE (apde, filename, nomeshload, nogeometryload);
    return apde;
  }


} // namespace

