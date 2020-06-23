#ifndef FILE_PDE
#define FILE_PDE

/*********************************************************************/
/* File:   pde.hpp                                                   */
/* Author: Joachim Schoeberl                                         */
/* Date:   08. Jul. 2000                                             */
/*********************************************************************/


// #include <../ngstd/evalfunc.hpp>

struct Tcl_Interp;
namespace ngcomp
{
  

  class EvalVariable : public NGS_Object
  {
  private:
    double * variable;
    EvalFunction * evaluator;
  public:
    EvalVariable(shared_ptr<MeshAccess> ama, const string & aname,
                 EvalFunction * aevaluator)
      : NGS_Object(ama,aname), variable(NULL), evaluator(aevaluator)
    { ; }
      
    void SetVariable(double & avariable)
    { 
      variable = &avariable; 
    }
    /*
    EvalFunction & GetEvaluator(void)
    {
      return evaluator;
    }
    */
    double Evaluate()
    {
      if(variable)
	return *variable = evaluator->Eval((double*)(0));
      else
	return evaluator->Eval((double*)(0));
    }

    virtual string GetClassName () const
    {
      return "EvalVariable";
    }
 
  };



  class LabelStatement : public NGS_Object
  {
    string label;
  public:
    LabelStatement(shared_ptr<MeshAccess> ama, const string & aname,
                  string alabel)
      : NGS_Object(ama,aname)
    { 
      label = alabel;
    }

    string GetLabel() const { return label; }
  };

  class GotoStatement : public NGS_Object
  {
    string target;
  public:
    GotoStatement(shared_ptr<MeshAccess> ama, const string & aname,
                  string atarget)
      : NGS_Object(ama,aname)
    { 
      target = atarget;
    }

    string GetTarget() const
    {
      return target;
    }
  };

  class ConditionalGotoStatement : public NGS_Object
  {
    EvalFunction * fun;
    string target;
  public:
    ConditionalGotoStatement(shared_ptr<MeshAccess> ama, const string & aname,
                  EvalFunction * afun, string atarget)
      : NGS_Object(ama,aname)
    { 
      fun = afun;
      target = atarget;
    }
    
    bool EvaluateCondition () const
    {
      return fun -> Eval ((double*)NULL);
    }

    string GetTarget() const
    {
      return target;
    }
  };

  class StopStatement : public NGS_Object
  {
  public:
    StopStatement(shared_ptr<MeshAccess> ama, const string & aname)
      : NGS_Object(ama,aname) 
    { ; }
  };


  NGS_DLL_HEADER void BuildLineIntegratorCurvePoints ( const string filename,
					const MeshAccess & ma,
					Integrator & integrator,
					bool draw = true);
  NGS_DLL_HEADER void BuildLineIntegratorCurvePoints(istream & infile,
					const MeshAccess & ma,
					Integrator & integrator,
					bool draw = true);



  /** 
      Description of partial differential equation
  */
  class NGS_DLL_HEADER PDE
  {
    ///
    Array<shared_ptr<MeshAccess>> mas;

    ///
    string geometryfilename;
    ///
    string meshfilename;
    ///
    SymbolTable<double> constants;
    ///
    SymbolTable<string*> string_constants;
    ///
    SymbolTable<shared_ptr<double>> variables;
    ///
    SymbolTable<GenericVariable> generic_variables;
    ///
    SymbolTable<Flags> flags;
    ///
    Array<shared_ptr<EvalVariable>> evaluators;
    ///
    SymbolTable<shared_ptr<CoefficientFunction>> coefficients;
    ///
    SymbolTable<shared_ptr<FESpace>> spaces;
    ///
    SymbolTable<shared_ptr<GridFunction>> gridfunctions;
    ///
    SymbolTable<shared_ptr<BilinearForm>> bilinearforms;
    ///
    SymbolTable<shared_ptr<LinearForm>> linearforms;
    ///
    SymbolTable<shared_ptr<Preconditioner>> preconditioners;
    ///
    SymbolTable<shared_ptr<NumProc>> numprocs;
  
    /// 
    bool isgood;
    ///
    Array<Integrator*> CurvePointIntegrators;
    Array<string*> CurvePointIntegratorFilenames;

    ///
    Array<shared_ptr<NGS_Object>> todo;

    ///
    int levelsolved;
    //
    string matfile;

    string filename;
    string workingdirectory;

    string evaluatefiles;

    /// program counter
    int pc;  

    /// a hack 
    Tcl_Interp * tcl_interpreter;

  public:
    ///
    PDE();
    ///
    /*
    PDE (const string & filename)
      : PDE() 
    {
      LoadPDE (filename);
    }
    */
    ///
    ~PDE();

    ///
    void SaveSolution (const string & filename, const bool ascii = false);
    ///
    void LoadSolution (const string & filename, const bool ascii = false);
    ///  
    void Solve ();
    // void SolveBVP () { Solve(); }
    
    void DoArchive (Archive & archive);

    ///
    void PrintReport (ostream & ost) const;

    ///
    void PrintMemoryUsage (ostream & ost);

    ///
    shared_ptr<MeshAccess> GetMeshAccess (int nr = 0) const  
    {
      if (nr >= mas.Size()) throw Exception ("PDE::GetMeshAccess, no mesh");
      return mas[nr]; 
    }
    ///
    void AddMeshAccess (shared_ptr<MeshAccess> ma) { mas.Append (ma); }
    ///
    bool ConstantUsed (const string & aname) const;
    ///
    double GetConstant (const string & aname, bool opt = 0) const;
    ///
    bool StringConstantUsed (const string & aname) const;
    ///
    string GetStringConstant (const string & aname, bool opt = 0) const;
    ///
    bool VariableUsed (const string & aname) const;

    ///
    double & GetVariable (const string & aname, bool opt = 0);

    bool FlagsUsed (const string & name) const { return flags.Used(name); }
    
    const Flags & GetFlags (const string & name) const { return flags[name]; }

    ///
    // CoefficientFunction * GetCoefficientFunction (const string & name, bool opt = 0);
    shared_ptr<CoefficientFunction> GetCoefficientFunction (const string & name, bool opt = 0) const;
    ///
    // FESpace * GetFESpace (const string & name, bool opt = 0);
    shared_ptr<FESpace> GetFESpace (const string & name, bool opt = 0) const;
    ///
    /*
    GridFunction * GetGridFunction (const string & name, bool opt = 0)
    {
      return GridFunctionPtr(name, opt).get();
    }
    */
    ///
    shared_ptr<GridFunction> GetGridFunction (const string & name, bool opt = 0) const;
    ///
    shared_ptr<BilinearForm> GetBilinearForm (const string & name, bool opt = 0) const;
    ///
    shared_ptr<LinearForm> GetLinearForm (const string & name, bool opt = 0) const;
    ///
    shared_ptr<Preconditioner> GetPreconditioner (const string & name, bool opt = 0) const;
    ///
    shared_ptr<NumProc> GetNumProc (const string & name, bool opt = 0) const;
    ///
    // shared_ptr<CoefficientFunction> GetCoefficientFunction (const string & name, bool opt = 0) const;
    ///
    // const FESpace * GetFESpace (const string & name, bool opt = 0) const;
    ///
    // const GridFunction * GetGridFunction (const string & name, bool opt = 0) const;
    ///
    // const BilinearForm * GetBilinearForm (const string & name, bool opt = 0) const;
    ///
    // const LinearForm * GetLinearForm (const string & name, bool opt = 0) const;
    ///
    // const Preconditioner * GetPreconditioner (const string & name, bool opt = 0) const;
    ///
    // const NumProc * GetNumProc (const string & name, bool opt = 0) const;





    ///
    void AddConstant (const string & name, double val);
    ///
    void AddStringConstant (const string & name, const string & val);
    ///
    void AddVariable (const string & name, double val, int im = 5);
    ///
    void AddVariable (const string & name, shared_ptr<EvalVariable> eval);
    ///
    void AddFlags (const string & name, const Flags & aflags);
    ///
    void AddVariableEvaluation (shared_ptr<EvalVariable> eval);
    ///
    void AddCoefficientFunction (const string & name, shared_ptr<CoefficientFunction> fun);
    ///
    shared_ptr<FESpace> AddFESpace (const string & name, const Flags & flags);
    ///
    void AddFESpace (const string & name, shared_ptr<FESpace> space);
    ///
    shared_ptr<GridFunction> AddGridFunction (const string & name, const Flags & flags);
    ///
    void AddGridFunction (const string & name, shared_ptr<GridFunction> gf, bool addcf = false);
    ///
    shared_ptr<BilinearForm> AddBilinearForm (const string & name, const Flags & flags);
    ///
    void AddBilinearForm (const string & name, shared_ptr<BilinearForm> bf);
    ///
    shared_ptr<LinearForm> AddLinearForm (const string & name, const Flags & flags);
    void AddLinearForm (const string & name, shared_ptr<LinearForm> lf);
    ///
    shared_ptr<Preconditioner> AddPreconditioner (const string & name, const Flags & flags);
    void AddPreconditioner (const string & name, shared_ptr<Preconditioner> pre);
    ///
    void AddNumProc (const string & name, shared_ptr<NumProc> np);


    ///
    void AddBilinearFormIntegrator (const string & name, shared_ptr<BilinearFormIntegrator> part,
				    const bool deletable = true);

    /*
    void AddIndependentBilinearFormIntegrator (const string & name, BilinearFormIntegrator * part,
					       const int master, const int other,
					       const bool deletable = true);
    */

    ///
    void AddLinearFormIntegrator (const string & name, shared_ptr<LinearFormIntegrator> part);

    ///
    void SetLineIntegratorCurvePointInfo(const string & filename,
					 Integrator * integrator);

    ///
    SymbolTable<double> & GetConstantTable ()
    { return constants; }
    ///
    SymbolTable<string*> & GetStringConstantTable ()
    { return string_constants; }
    ///
    SymbolTable<shared_ptr<double>> & GetVariableTable ()
    { return variables; }

    SymbolTable<GenericVariable> & GenericVariables() { return generic_variables; }
    const SymbolTable<GenericVariable> & GenericVariables() const { return generic_variables; }

    ///
    SymbolTable<shared_ptr<CoefficientFunction>> & GetCoefficientTable ()
    { return coefficients; }
    ///
    SymbolTable<shared_ptr<FESpace>> & GetSpaceTable ()
    { return spaces; }
    ///
    SymbolTable<Flags> & GetFlagsTable ()
    { return flags; }
    ///
    SymbolTable<shared_ptr<GridFunction>> & GetGridFunctionTable()
    { return gridfunctions; }
    ///
    SymbolTable<shared_ptr<BilinearForm>> & GetBilinearFormTable()
    { return bilinearforms; }
    ///
    SymbolTable<shared_ptr<LinearForm>> & GetLinearFormTable()
    { return linearforms; }
    ///
    SymbolTable<shared_ptr<Preconditioner>> & GetPreconditionerTable ()
    { return preconditioners; }
    ///
    SymbolTable<shared_ptr<NumProc>> & GetNumProcTable()
    { return numprocs; }  


    int GetNStatements () { return todo.Size(); }
    const shared_ptr<NGS_Object> & GetStatement (int nr) { return todo[nr]; }
    void AddControlStatement (shared_ptr<NGS_Object> obj)
    {
      todo.Append (obj);
    }

    ///
    bool IsGood () { return isgood; }
    ///
    void SetGood (bool agood) { isgood = agood; }
  
    int GetPC() const { return pc; }
    
  
    string GetMatfile() const
    { return matfile;} 
  
    void SetMatfile(const string str)
    { matfile = str; } 

    string GetDirectory() const
    { return workingdirectory; }
  
    void SetDirectory(const string str)
    { workingdirectory = str; }

    string GetFilename() const
    { return filename; }

    void SetFilename(const string str);


    string & GetEvaluateFiles(void) { return evaluatefiles; }
    // #ifdef SOCKETS
    //   ClientSocketAccess & GetClientSocketAccess(void);
    //   void SolveBVPClientServer (LocalHeap & lh);
    // private:
    //   void CallNumProcsClientServer(const int position, LocalHeap & lh);
    // #endif
    void SetMeshFileName(const string ameshfilename) { meshfilename = ameshfilename; }
    void SetGeoFileName ( const string ageofilename) { geometryfilename = ageofilename; }

    string GetMeshFileName() const { return meshfilename; }
    string GetGeoFileName () const { return geometryfilename; }

    void WritePDEFile ( string abspdefile, string geofile, 
			string meshfile, string matfile, string oldpdefile );

    Tcl_Interp * GetTclInterpreter() const { return tcl_interpreter; }
    void SetTclInterpreter(Tcl_Interp * inter) { tcl_interpreter = inter; }


    void Tcl_Eval (string str);
    /*
    void Tcl_CreateCommand(Tcl_Interp *interp,
			   CONST char *cmdName, Tcl_CmdProc *proc,
			   ClientData clientData,
			   Tcl_CmdDeleteProc *deleteProc);
    */

#ifdef ASTRID
    void SaveZipSolution (const string & filename, const bool ascii = false);
    ///
    void LoadZipSolution (const string & filename, const bool ascii = false);
#endif
    ///
  };


  ///
  NGS_DLL_HEADER extern shared_ptr<PDE> LoadPDE (const string & filename, const bool nomeshload = false, const bool nogeometryload = false);
  ///
  NGS_DLL_HEADER extern shared_ptr<PDE> LoadPDE (istream & input, const bool nomeshload = false, const bool nogeometryload = false);

  
  NGS_DLL_HEADER extern void LoadPDE (shared_ptr<PDE> pde, const string & filename, const bool nomeshload = false, const bool nogeometryload = false);
  ///
  NGS_DLL_HEADER extern void LoadPDE (shared_ptr<PDE> pde, istream & input, const bool nomeshload = false, const bool nogeometryload = false);



  inline ostream & operator<< (ostream & ost, const PDE & pde)
  {
    pde.PrintReport (ost);
    return ost;
  }

}

#endif
