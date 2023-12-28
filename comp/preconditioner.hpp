#ifndef FILE_PRECONDITIONER
#define FILE_PRECONDITIONER

/*********************************************************************/
/* File:   preconditioner.hh                                         */
/* Author: Joachim Schoeberl                                         */
/* Date:   10. Jul. 2000                                             */
/*********************************************************************/

#include <la.hpp>
#include <mgpre.hpp>

namespace ngcomp
{
  using namespace ngla;
  
  /**
     Base class for preconditioners.
  */
  class NGS_DLL_HEADER Preconditioner : public BaseMatrix, public NGS_Object
  {
  private:
    weak_ptr<BilinearForm> bf;
    bool is_registered = false;
  protected:
    bool test;
    bool timing;
    bool print;

    /// if true, the update in SolveBVP() is ignored, Update() has to be called explicitly.
    bool laterupdate;

    double * testresult_ok;
    double * testresult_min;
    double * testresult_max;
  
    // for calculation of eigenvalues
    bool uselapack;

    int on_proc;

  public:
    Preconditioner (shared_ptr<BilinearForm> bfa, const Flags & aflags,
		    const string aname = "precond");
    ///
    virtual ~Preconditioner ();
  
    ///
    virtual bool LaterUpdate (void) { return laterupdate; }
    ///
    virtual void Update ()  override = 0;
    ///
    virtual void CleanUpLevel () { ; }
    ///
    virtual const BaseMatrix & GetMatrix() const = 0;
    /*
    {
      return *this; 
    }
    */
    virtual shared_ptr<BaseMatrix> GetMatrixPtr()
    {
      return BaseMatrix::SharedFromThis<BaseMatrix>();
    }

    shared_ptr<BilinearForm> GetBilinearForm() const
    {
      return bf.lock();
    }
    
    virtual bool IsComplex() const override { return GetMatrix().IsComplex(); }
        
    ///
    virtual void Mult (const BaseVector & x, BaseVector & y) const override
    {
      GetMatrix().Mult(x, y);
    }

    virtual void MultAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      GetMatrix().MultAdd(s, x, y);
    }

    virtual void MultTrans (const BaseVector & x, BaseVector & y) const override
    {
      GetMatrix().MultTrans(x, y);
    }

    virtual void MultTransAdd (double s, const BaseVector & x, BaseVector & y) const override
    {
      GetMatrix().MultTransAdd(s, x, y);
    }



    
    virtual void InitLevel (shared_ptr<BitArray> freedofs = NULL) { ; }
    virtual void FinalizeLevel (const ngla::BaseMatrix * mat = NULL) { ; }
    virtual void AddElementMatrix (FlatArray<int> dnums,
				   const FlatMatrix<double> & elmat,
				   ElementId ei, 
				   LocalHeap & lh) { ; }

    virtual void AddElementMatrix (FlatArray<int> dnums,
				   const FlatMatrix<Complex> & elmat,
				   ElementId ei, 
				   LocalHeap & lh) { ; }



    virtual const BaseMatrix & GetAMatrix() const
    { throw Exception ("Preconditioner, A-Matrix not available"); }

    ///
    virtual const char * ClassName() const
    { return "base-class Preconditioner"; }

    virtual AutoVector CreateRowVector () const override
    {
      return GetAMatrix().CreateColVector();
    }

    virtual AutoVector CreateColVector () const override
    {
      return GetAMatrix().CreateRowVector();
    }

    virtual void PrintReport (ostream & ost) const override
    {
      ost << "type = " << ClassName() << endl;
    }

    virtual Array<MemoryUsage> GetMemoryUsage () const override
    {
      throw Exception(string("MemoryUsage not implemented for preconditioner ")+ClassName());
    }

    virtual int VHeight() const override { return GetMatrix().VHeight();}
    virtual int VWidth() const override { return GetMatrix().VWidth();}

    void Test () const;
    void Timing () const;
    void ThrowPreconditionerNotReady() const;
    const Flags & GetFlags() const { return flags; }

    using BaseMatrix::shared_from_this;
    using NGS_Object::GetMemoryTracer;
  };



  ///




  ///
  class TwoLevelPreconditioner : public Preconditioner
  {
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    shared_ptr<Preconditioner> cpre;
    ///
    ngmg::TwoLevelMatrix * premat;
    ///
    // int smoothingsteps;
  public:
    ///
    virtual ~TwoLevelPreconditioner();

    ///
    virtual void Update ();
    ///
    virtual const BaseMatrix & GetMatrix() const;
    // { return *new SparseMatrix<double> (1,1); } // *premat; }
    // { return *premat; }
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



  /**
     Multigrid preconditioner.
     High level objects, contains a \Ref{MultigridPreconditioner}
  */
  class NGS_DLL_HEADER MGPreconditioner : public Preconditioner
  {
    ///
    shared_ptr<ngmg::MultigridPreconditioner> mgp;
    ///
    shared_ptr<ngmg::TwoLevelMatrix> tlp;
    ///
    shared_ptr<BilinearForm> bfa;
    ///
    // MGPreconditioner * low_order_preconditioner;
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
    MGPreconditioner (shared_ptr<BilinearForm> bfa, const Flags & aflags,
		      const string aname = "mgprecond");
    ///
    virtual ~MGPreconditioner() { ; }

    void FreeSmootherMem(void);

    virtual void FinalizeLevel (const BaseMatrix * mat) override
    {
      Update();
    }

    ///
    virtual void Update () override;
    ///
    virtual void CleanUpLevel () override;
    ///
    virtual const BaseMatrix & GetMatrix() const override;
    ///
    virtual const BaseMatrix & GetAMatrix() const override
    {
      return bfa->GetMatrix();
    }
    ///
    virtual const char * ClassName() const override
    { return "Multigrid Preconditioner"; }

    virtual void PrintReport (ostream & ost) const override;

    virtual Array<MemoryUsage> GetMemoryUsage () const override;

    void MgTest () const;

    void SetDirectSolverCluster(shared_ptr<Array<int>> cluster);
    void SetCoarsePreconditioner(shared_ptr<Preconditioner> prec);
  };

  class CommutingAMGPreconditioner : public Preconditioner
  {
  protected:
    shared_ptr<BilinearForm> bfa;
    // CommutingAMG * amg;
    BaseMatrix * amg;
    shared_ptr<CoefficientFunction> coefe, coeff, coefse;
    bool hcurl;
    bool coarsegrid;
    int levels;
  public:
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
    virtual ~NonsymmetricPreconditioner();
    ///
    virtual bool IsComplex() const { return cm->IsComplex(); }

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
      function<shared_ptr<Preconditioner>(shared_ptr<BilinearForm>,const Flags &,const string)> creatorbf;
      DocInfo docinfo;
      
      PreconditionerInfo (const string & aname,
                          function<shared_ptr<Preconditioner>(shared_ptr<BilinearForm>, const Flags &, const string)> acreateorbf,
                          DocInfo adocinfo);
    };
  
    Array<unique_ptr<PreconditionerInfo>> prea;
  public:
    PreconditionerClasses() = default;
    PreconditionerClasses(const PreconditionerClasses &) = default;
    ~PreconditionerClasses() = default;

    void AddPreconditioner (const string & aname, 
			    function<shared_ptr<Preconditioner>(shared_ptr<BilinearForm>, const Flags&, const string)> createorbf,
                            DocInfo docinfo = DocInfo());
      
    const Array<unique_ptr<PreconditionerInfo>> & GetPreconditioners() { return prea; }
    const PreconditionerInfo * GetPreconditioner(const string & name);

    void Print (ostream & ost) const;
    void Cleanup();
  };
 
  extern NGS_DLL_HEADER PreconditionerClasses & GetPreconditionerClasses ();

  template <typename PRECOND>
  class RegisterPreconditioner
  {
  public:
    RegisterPreconditioner (string label, bool isparallel = true)
    {
      GetPreconditionerClasses().AddPreconditioner (label, CreateBF);
    }
    
    static shared_ptr<Preconditioner> CreateBF (shared_ptr<BilinearForm> bfa, const Flags & flags, const string & name)
    {
      return make_shared<PRECOND> (bfa, flags, name);
    }
  };

}

#endif

