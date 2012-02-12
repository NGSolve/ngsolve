#ifndef FILE_PDE
#define FILE_PDE

/*********************************************************************/
/* File:   pde.hpp                                                   */
/* Author: Joachim Schoeberl                                         */
/* Date:   08. Jul. 2000                                             */
/*********************************************************************/

namespace ngsolve
{

  class EvalVariable : public NGS_Object
  {
  private:
    double * variable;

    EvalFunction evaluator;


  public:
    EvalVariable(const MeshAccess & ama, const string & aname);

    void SetVariable(double & avariable);

    EvalFunction & GetEvaluator(void)
    {
      return evaluator;
    }

    double Evaluate(void);

    virtual string GetClassName () const
    {
      return "EvalVariable";
    }
 
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
    MeshAccess & ma;

    ///
    string geometryfilename;
    ///
    string meshfilename;
    ///
    SymbolTable<double> constants;
    ///
    SymbolTable<string*> string_constants;
    ///
    SymbolTable<double> variables;
    ///
    Array<EvalVariable*> evaluators;
    ///
    SymbolTable<CoefficientFunction*> coefficients;
    ///
    SymbolTable<FESpace*> spaces;
    ///
    SymbolTable<GridFunction*> gridfunctions;
    ///
    SymbolTable<BilinearForm*> bilinearforms;
    ///
    SymbolTable<LinearForm*> linearforms;
    ///
    SymbolTable<Preconditioner*> preconditioners;
    ///
    SymbolTable<NumProc*> numprocs;
  
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

    // #ifdef SOCKETS
    //   ClientSocketAccess sa;
    // #endif

  public:
    ///
    PDE(MeshAccess & ama);
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
    void SolveBVP ();

    ///
    void PrintReport (ostream & ost);

    ///
    void PrintMemoryUsage (ostream & ost);

    ///
    const MeshAccess & GetMeshAccess() const { return ma; }

    ///
    MeshAccess & GetMeshAccess() { return ma; }
  


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
    CoefficientFunction * GetCoefficientFunction (const string & name, bool opt = 0);
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
    void AddVariable (const string & name, double val);
    ///
    void AddVariable (const string & name, EvalVariable * eval);
    ///
    void AddCoefficientFunction (const string & name, CoefficientFunction* fun);
    ///
    FESpace * AddFESpace (const string & name, Flags & flags);
    ///
    void AddFESpace (const string & name, FESpace * space);
    ///
    GridFunction * AddGridFunction (const string & name, Flags & flags);
    ///
    void AddGridFunction (const string & name, GridFunction * gf, bool addcf = false);
    ///
    BilinearForm * AddBilinearForm (const string & name, Flags & flags);
    ///
    LinearForm * AddLinearForm (const string & name, Flags & flags);
    ///
    Preconditioner * AddPreconditioner (const string & name, Flags & flags);
    ///
    void AddNumProc (const string & name, NumProc * np);


    ///
    void AddBilinearFormIntegrator (const string & name, BilinearFormIntegrator * part,
				    const bool deletable = true);
    ///
    void AddIndependentBilinearFormIntegrator (const string & name, BilinearFormIntegrator * part,
					       const int master, const int slave,
					       const bool deletable = true);
    ///
    void AddLinearFormIntegrator (const string & name, LinearFormIntegrator * part);

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
    SymbolTable<double> & GetVariableTable ()
    { return variables; }
    ///
    SymbolTable<CoefficientFunction*> & GetCoefficientTable ()
    { return coefficients; }
    ///
    SymbolTable<FESpace*> & GetSpaceTable ()
    { return spaces; }
    ///
    SymbolTable<GridFunction*> & GetGridFunctionTable()
    { return gridfunctions; }
    ///
    SymbolTable<BilinearForm*> & GetBilinearFormTable()
    { return bilinearforms; }
    ///
    SymbolTable<LinearForm*> & GetLinearFormTable()
    { return linearforms; }
    ///
    SymbolTable<Preconditioner*> & GetPreconditionerTable ()
    { return preconditioners; }
    ///
    SymbolTable<NumProc*> & GetNumProcTable()
    { return numprocs; }  

    ///
    bool IsGood () { return isgood; }
    ///
    void SetGood (bool agood) { isgood = agood; }
  


    /// a hack 
    Tcl_Interp * tcl_interpreter;
  
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

#ifdef ASTRID
    void SaveZipSolution (const string & filename, const bool ascii = false);
    ///
    void LoadZipSolution (const string & filename, const bool ascii = false);
#endif
    ///
  };

}

#endif
