#ifndef FILE_PRECONDITIONER
#define FILE_PRECONDITIONER

/*********************************************************************/
/* File:   preconditioner.hh                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Jul. 2000                                             */
/*********************************************************************/



namespace ngsolve
{
  class PDE;
}

namespace ngcomp
{
  using ngsolve::PDE;

  /**
     Base class for preconditioners.
  */
  class NGS_DLL_HEADER Preconditioner : public NGS_Object, public BaseMatrix
  {
  protected:
    bool test;
    bool timing;
    bool print;

    /// if true, the update in SolveBVP() is ignored, Update() has to be called explicitely.
    bool laterupdate;

    double * testresult_ok;
    double * testresult_min;
    double * testresult_max;
  
    Flags flags; 

    // for calculation of eigenvalues
    bool uselapack;

    int on_proc;

  public:
    ///
    //Preconditioner ();
    ///
    //Preconditioner (const Flags & flags);
    ///
    Preconditioner (const PDE * const apde, const Flags & aflags,
		    const string aname = "precond");
    ///
    virtual ~Preconditioner ();
  
    ///
    virtual bool LaterUpdate (void) { return laterupdate; }
    ///
    virtual void Update () = 0;
    ///
    virtual void CleanUpLevel () { ; }
    ///
    virtual const BaseMatrix & GetMatrix() const
    {
      return *this; 
    }
    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const
    {
      GetMatrix().Mult(x, y);
    }

    virtual void InitLevel () { ; }
    virtual void FinalizeLevel () { ; }
    virtual void AddElementMatrix (FlatArray<int> dnums,
				   const FlatMatrix<double> & elmat,
				   bool inner_element, int elnr,
				   LocalHeap & lh) { ; }

    virtual void AddElementMatrix (FlatArray<int> dnums,
				   const FlatMatrix<Complex> & elmat,
				   bool inner_element, int elnr,
				   LocalHeap & lh) { ; }



    virtual const BaseMatrix & GetAMatrix() const
    { throw Exception ("Preconditioner, A-Matrix not available"); }

    ///
    virtual const char * ClassName() const
    { return "base-class Preconditioner"; }


    virtual void PrintReport (ostream & ost)
    {
      ost << "type = " << ClassName() << endl;
    }

    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const
    {
      cout << "MemoryUsage not implemented for preconditioner " << ClassName() << endl;
    }

    void Test () const;
    void Timing () const;
  };


  /**
     Multigrid preconditioner.
     High level objects, contains a \Ref{MultigridPreconditioner} 
  */
  class NGS_DLL_HEADER MGPreconditioner : public Preconditioner
  {
    ///
    ngmg::MultigridPreconditioner * mgp;
    ///
    ngmg::TwoLevelMatrix * tlp;
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    MGPreconditioner * low_order_preconditioner;
    ///
    shared_ptr<Preconditioner> coarse_pre;
    ///
    int finesmoothingsteps;
    ///
    string smoothertype;
    /// 
    bool mgtest; 
    string mgfile; 
    int mgnumber; 
  
    string inversetype;

  public:
    ///
    MGPreconditioner (const PDE & pde, const Flags & aflags,
		      const string aname = "mgprecond");
    ///
    virtual ~MGPreconditioner();

    void FreeSmootherMem(void);

    ///
    virtual void Update ();
    ///
    virtual void CleanUpLevel ();
    ///
    virtual const BaseMatrix & GetMatrix() const;
    ///
    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa->GetMatrix(); 
    }
    ///
    virtual const char * ClassName() const
    { return "Multigrid Preconditioner"; }

    virtual void PrintReport (ostream & ost);

    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;

    void MgTest () const;
  };



  ///


  /**
     Local (Block-Jacobi or Block-Gauss-Seidel) preconditioner
  */
  class LocalPreconditioner : public Preconditioner
  {
  protected:
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    BaseMatrix * jacobi;
    ///
    bool block;
    bool locprectest; 
    string locprecfile; 

    string ct;
    shared_ptr<Preconditioner> coarse_pre;
  public:
    ///
    LocalPreconditioner (PDE * pde, const Flags & aflags,
			 const string aname = "localprecond");
    ///
    virtual ~LocalPreconditioner();
    ///
    virtual void Update ();
    ///
    virtual const BaseMatrix & GetMatrix() const;
    ///
    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa->GetMatrix(); 
    }
    ///
    virtual const char * ClassName() const
    { return "Local Preconditioner"; }
    void LocPrecTest () const;
  };




  ///
  class TwoLevelPreconditioner : public Preconditioner
  {
    ///
    PDE * pde;
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<Preconditioner> cpre;
    ///
    ngmg::TwoLevelMatrix * premat;
    ///
    int smoothingsteps;
  public:
    ///
    TwoLevelPreconditioner (PDE * apde, const Flags & aflags,
			    const string aname = "twolevelprecond");
    ///
    virtual ~TwoLevelPreconditioner();

    ///
    virtual void Update ();
    ///
    virtual const BaseMatrix & GetMatrix() const
    { return *new SparseMatrix<double> (1,1); } // *premat; }
    ///
    virtual const char * ClassName() const
    { return "TwoLevel Preconditioner"; }
  };








  ///
  class NGS_DLL_HEADER ComplexPreconditioner : public Preconditioner
  {
  protected:
    ///
    shared_ptr<Preconditioner> creal;
    ///
    // Real2ComplexMatrix<double,Complex> cm;
    int dim;
    BaseMatrix * cm;
  public:
    ///
    ComplexPreconditioner (PDE * apde, const Flags & aflags,
			   const string aname = "complexprecond");
    ///
    virtual ~ComplexPreconditioner();
    ///
    virtual void Update ();
    ///
    virtual const BaseMatrix & GetMatrix() const
    { 
      return *cm; 
    }
    ///
    virtual const char * ClassName() const
    { return "Complex Preconditioner"; }
  };




  ///
  class ChebychevPreconditioner : public Preconditioner
  {
  protected:
    ///
    shared_ptr<Preconditioner> csimple;
    /// 
    ChebyshevIteration * cm;
    /// 
    shared_ptr<BilinearForm> bfa;
    ///
    int steps; 
  public:
    ///
    ChebychevPreconditioner (PDE * apde, const Flags & aflags,
			     const string aname = "chebychevprecond");
    ///
    virtual ~ChebychevPreconditioner();
    ///
    virtual void Update ();
    ///
    virtual const BaseMatrix & GetMatrix() const
    { 
      return *cm; 
    }
    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa->GetMatrix(); 
    }

    ///
    virtual const char * ClassName() const
    { return "Chebychev Preconditioner"; }
  };




  class CommutingAMGPreconditioner : public Preconditioner
  {
  protected:
    PDE * pde;
    shared_ptr<BilinearForm> bfa;
    // CommutingAMG * amg;
    BaseMatrix * amg;
    shared_ptr<CoefficientFunction> coefe, coeff, coefse;
    bool hcurl;
    bool coarsegrid;
    int levels;
  public:
    CommutingAMGPreconditioner (PDE * apde, const Flags & aflags,
				const string aname = "commutingamgprecond");

    virtual ~CommutingAMGPreconditioner ();

    virtual void Update ();
    ///

    virtual const BaseMatrix & GetAMatrix() const
    {
      return bfa->GetMatrix(); 
    }

    virtual const BaseMatrix & GetMatrix() const
    { 
      return *amg; 
    }
    ///
    virtual void CleanUpLevel ();

    ///
    virtual const char * ClassName() const
    { return "CommutingAMG Preconditioner"; }
  };






  // added 08/19/2003:

  ////////////////////////////////////////////////////////////////////////////////
  //
  // special preconditioner for system
  //   (  A   M  )
  //   ( -M   A  )
  //
  // 
  // C = (  1  1  ) (  A+M       )
  //     ( -1  1  ) (       A+M  )
  //
  ////////////////////////////////////////////////////////////////////////////////
  class NonsymmetricPreconditioner : public Preconditioner
  {
  protected:
    ///
    shared_ptr<Preconditioner> cbase;
    ///
    int dim;
    BaseMatrix * cm;
  public:
    ///
    NonsymmetricPreconditioner (PDE * apde, const Flags & aflags,
				const string aname = "nonsymmetricprecond");
    ///
    virtual ~NonsymmetricPreconditioner();
    ///
    virtual void Update ();
    ///
    virtual const BaseMatrix & GetMatrix() const
    { 
      return *cm; 
    }
    ///
    virtual const char * ClassName() const
    { return "Nonsymmetric Preconditioner"; }
  };







  /// Registered Preconditioner classes
  class NGS_DLL_HEADER PreconditionerClasses
  {
  public:
    struct PreconditionerInfo
    {
      string name;
      Preconditioner* (*creator)(const PDE & pde, const Flags & aflags, const string & name);
      PreconditionerInfo (const string & aname,
			  Preconditioner* (*acreator)(const PDE & pde, 
						      const Flags & aflags,
						      const string & name));
    };
  
    Array<PreconditionerInfo*> prea;
  public:
    PreconditionerClasses();
    ~PreconditionerClasses();  
    void AddPreconditioner (const string & aname, 
			    Preconditioner* (*acreator)(const PDE & pde, 
							const Flags & aflags, 
							const string & name));
  
    const Array<PreconditionerInfo*> & GetPreconditioners() { return prea; }
    const PreconditionerInfo * GetPreconditioner(const string & name);

    void Print (ostream & ost) const;
  };
 
  extern NGS_DLL_HEADER PreconditionerClasses & GetPreconditionerClasses ();

  template <typename PRECOND>
  class RegisterPreconditioner
  {
  public:
    RegisterPreconditioner (string label, bool isparallel = true)
    {
      GetPreconditionerClasses().AddPreconditioner (label, Create);
      // cout << "register preconditioner '" << label << "'" << endl;
    }
    
    static Preconditioner * Create (const PDE & pde, const Flags & flags, const string & name)
    {
      return new PRECOND (pde, flags, name);
    }
  };

}

#endif

