#ifndef FILE_PDE
#define FILE_PDE

/*********************************************************************/
/* File:   pde.hpp                                                   */
/* Author: Joachim Schoeberl                                         */
/* Date:   08. Jul. 2000                                             */
/*********************************************************************/


// #include <../ngstd/evalfunc.hpp>

namespace ngsolve
{


  class EvalVariable : public NGS_Object
  {
  private:
    double * variable;
    EvalFunction * evaluator;
  public:
    EvalVariable(const MeshAccess & ama, const string & aname,
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
    LabelStatement(const MeshAccess & ama, const string & aname,
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
    GotoStatement(const MeshAccess & ama, const string & aname,
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
    ConditionalGotoStatement(const MeshAccess & ama, const string & aname,
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
    StopStatement(const MeshAccess & ama, const string & aname)
      : NGS_Object(ama,aname) 
    { ; }
  };


  void BuildLineIntegratorCurvePoints ( const string filename,
					const MeshAccess & ma,
					Integrator & integrator,
					bool draw = true);
  void BuildLineIntegratorCurvePoints ( istream & infile,
					const MeshAccess & ma,
					Integrator & integrator,
					bool draw = true);



  /** 
      Description of partial differential equation
  */
  class NGS_DLL_HEADER PDE
  {
    ///
    Array<MeshAccess*> mas;

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
    Array<EvalVariable*> evaluators;
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
    Array<NGS_Object*> todo;

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
    PDE (const string & filename)
      : PDE() 
    {
      LoadPDE (filename);
    }
    ///
    ~PDE();

    ///
    void LoadPDE (const string & filename, const bool nomeshload = false, const bool nogeometryload = false);
    ///
    void LoadPDE (istream & input, const bool nomeshload = false, const bool nogeometryload = false);
    ///
    void SavePDE (const string & filename);

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
    const MeshAccess & GetMeshAccess (int nr = 0) const  
    { 
      return *mas[nr]; 
    }
    ///
    MeshAccess & GetMeshAccess (int nr = 0) 
    {
      return *mas[nr]; 
    }
    ///
    void AddMeshAccess (MeshAccess * ma) { mas.Append (ma); }
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
    ///
    shared_ptr<CoefficientFunction> GetCoefficientFunction (const string & name, bool opt = 0);
    ///
    FESpace * GetFESpace (const string & name, bool opt = 0);
    ///
    GridFunction * GetGridFunction (const string & name, bool opt = 0);
    ///
    BilinearForm * GetBilinearForm (const string & name, bool opt = 0);
    ///
    LinearForm * GetLinearForm (const string & name, bool opt = 0);
    ///
    Preconditioner * GetPreconditioner (const string & name, bool opt = 0);
    ///
    NumProc * GetNumProc (const string & name, bool opt = 0);


    ///
    const CoefficientFunction * GetCoefficientFunction (const string & name, bool opt = 0) const;
    ///
    const FESpace * GetFESpace (const string & name, bool opt = 0) const;
    ///
    const GridFunction * GetGridFunction (const string & name, bool opt = 0) const;
    ///
    const BilinearForm * GetBilinearForm (const string & name, bool opt = 0) const;
    ///
    const LinearForm * GetLinearForm (const string & name, bool opt = 0) const;
    ///
    const Preconditioner * GetPreconditioner (const string & name, bool opt = 0) const;
    ///
    const NumProc * GetNumProc (const string & name, bool opt = 0) const;





    ///
    void AddConstant (const string & name, double val);
    ///
    void AddStringConstant (const string & name, const string & val);
    ///
    void AddVariable (const string & name, double val, int im = 5);
    ///
    void AddVariable (const string & name, EvalVariable * eval);
    ///
    void AddVariableEvaluation (EvalVariable * eval);
    ///
    void AddCoefficientFunction (const string & name, shared_ptr<CoefficientFunction> fun);
    ///
    FESpace * AddFESpace (const string & name, const Flags & flags);
    ///
    void AddFESpace (const string & name, FESpace * space);
    ///
    GridFunction * AddGridFunction (const string & name, const Flags & flags);
    ///
    void AddGridFunction (const string & name, shared_ptr<GridFunction> gf, bool addcf = false);
    ///
    BilinearForm * AddBilinearForm (const string & name, const Flags & flags);
    ///
    LinearForm * AddLinearForm (const string & name, const Flags & flags);
    ///
    Preconditioner * AddPreconditioner (const string & name, const Flags & flags);
    ///
    void AddNumProc (const string & name, shared_ptr<NumProc> np);


    ///
    void AddBilinearFormIntegrator (const string & name, shared_ptr<BilinearFormIntegrator> part,
				    const bool deletable = true);

    /*
    void AddIndependentBilinearFormIntegrator (const string & name, BilinearFormIntegrator * part,
					       const int master, const int slave,
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
    const NGS_Object * GetStatement (int nr) { return todo[nr]; }
    void AddControlStatement (NGS_Object * obj)
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

}

#endif
