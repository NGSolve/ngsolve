#ifndef FILE_BILINEARFORM
#define FILE_BILINEARFORM

/*********************************************************************/
/* File:   bilinearform.hpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngcomp
{

  class LinearForm;
  class Preconditioner;

  /** 
      A bilinear-form.
      A bilinear-form provides the system matrix. 
      It knows about its definition in terms of integrators.
      In most cases, it is defined on two copies of the same space V,
      but it can also live on V x W.
  */
  class NGS_DLL_HEADER BilinearForm : public NGS_Object
  {
  protected:
    /// Finite element space
    const FESpace & fespace;
    /// Test-space if different from trial-space, otherwise NULL (non operational)
    const FESpace * fespace2;

    /// don't assemble matrix
    bool nonassemble;
    /// store only diagonal of matrix
    bool diagonal;
    /// store matrices on mesh hierarchy
    bool multilevel;
    /// galerkin projection of coarse grid matrices
    bool galerkin;
    /// complex forms are hermitean (non operational)
    bool hermitean;
    /// bilinear form is symmetric
    bool symmetric;
    /// bilinear form is symmetric and positive definite (experimental)
    bool spd;
    /// add epsilon for regularization
    double eps_regularization; 
    /// diagonal value for unused dofs
    double unuseddiag;

    /// low order bilinear-form, 0 if not used
    shared_ptr<BilinearForm> low_order_bilinear_form;

    /// modify linear form due to static condensation
    LinearForm * linearform;

    /// some preconditioners need element matrices 
    Array<Preconditioner*> preconditioners; 

    /// matrices (sparse, application, diagonal, ...)
    Array<shared_ptr<BaseMatrix>> mats;

    /// bilinearform-integrators
    Array<shared_ptr<BilinearFormIntegrator>> parts;

    /*
    Array<BilinearFormIntegrator*> independent_parts;
    Array<bool> independent_parts_deletable;
    Array< Vec<2,int> > independent_meshindex;
    */

    /// does some timing in assemble
    bool timing;
    /// prints assembled matrix to testout
    bool print;
    /// prints element matrices to testout
    bool printelmat;
    /// calculates eigenvalues of element matrices
    bool elmat_ev;
    /// does static condensation of internal dofs
    bool eliminate_internal;
    /// keeps matrices for reconstruction of internal dofs
    bool keep_internal;
    /// should A_ii itself be stored?!
    bool store_inner; 
    
    /// precomputes some data for each element
    bool precompute;
    /// precomputed element-wise data
    Array<void*> precomputed_data;
    /// output of norm of matrix entries
    bool checksum;

  public:
    /// generate a bilinear-form
    BilinearForm (const FESpace & afespace,
		  const string & aname, 
		  const Flags & flags);

    /// generate a bilinear-form
    BilinearForm (const FESpace & afespace, 
		  const FESpace & afespace2, 
		  const string & aname,
		  const Flags & flags);

    virtual ~BilinearForm ();
  

    ///
    void AddIntegrator (shared_ptr<BilinearFormIntegrator> bfi);

    /*
    void AddIndependentIntegrator (BilinearFormIntegrator * bfi,
				   const int master_surface,
				   const int slave)
    */
  
    ///
    int NumIntegrators () const 
    {
      return parts.Size(); 
    }

    /*
    int NumIndependentIntegrators(void) const {return independent_parts.Size();}
    */

    /// the i-th integrator
    shared_ptr<BilinearFormIntegrator> GetIntegrator (int i) const { return parts[i]; }

    // const BilinearFormIntegrator * GetIntegrator (int i) const { return parts[i]; }


    /*
    BilinearFormIntegrator * GetIndependentIntegrator (const int i) 
    { return independent_parts[i]; }

    const BilinearFormIntegrator * GetIndependentIntegrator (const int i) const 
    { return independent_parts[i]; }

    int GetIndependentMasterSurfaceIndex(const int i) const
    {
      return independent_meshindex[i](0);
    }
    int GetIndependentSlaveIndex(const int i) const
    {
      return independent_meshindex[i](1);
    }
    */


    /// for static condensation of internal bubbles
    void SetLinearForm (LinearForm * alf) { linearform = alf; }
    
    /// preconditioner gets element-matrix
    void SetPreconditioner (Preconditioner * pre) 
    { preconditioners.Append (pre); }

    /// generates matrix graph
    virtual MatrixGraph * GetGraph (int level, bool symmetric);

    /// assembles the matrix
    void Assemble (LocalHeap & lh);

    /// re-assembles the matrix.
    /// if reallocate is false, the existing matrix is reused
    void ReAssemble (LocalHeap & lh, bool reallocate = 0);

    /// assembles matrix at linearization point given by lin
    /// needed for Newton's method
    virtual void AssembleLinearization (const BaseVector & lin,
					LocalHeap & lh, 
					bool reallocate = 0) = 0;

    /// applies the matrix without assembling
    void ApplyMatrix (const BaseVector & x,
		      BaseVector & y) const
    {
      y = 0;
      AddMatrix (1, x, y);
    }

    /// y += val * Mat * x
    virtual void AddMatrix (double val, const BaseVector & x,
			    BaseVector & y) const = 0;
  
    /// y += val * Mat * x
    virtual void AddMatrix (Complex val, const BaseVector & x,
			    BaseVector & y) const = 0;
  
    /// y += val * lin.mat * x
    virtual void ApplyLinearizedMatrixAdd (double val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y) const = 0;
    /// y += val * lin.mat * x
    virtual void ApplyLinearizedMatrixAdd (Complex val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y) const = 0;

    /// evaulates internal energy (usually  1/2 x^T A x)
    virtual double Energy (const BaseVector & x) const = 0;

    /// returns the assembled matrix
    BaseMatrix & GetMatrix () const { return *mats.Last(); }
    /// returns the assembled matrix
    shared_ptr<BaseMatrix> MatrixPtr () const { return mats.Last(); }

    // operator const BaseMatrix& () const { return GetMatrix(); }

    /// returns the assembled matrix on a given level
    BaseMatrix & GetMatrix (int level) const { return *mats[level]; }


    // BaseMatrix & GetMatrix ()  { return *mats.Last(); }
    // BaseMatrix & GetMatrix (int level)  { return *mats[level]; }

    /// reconstruct internal dofs from static condensation 
    /// -A_ii^{-1} A_ib
    virtual BaseMatrix & GetHarmonicExtension () const = 0;

    /// modify rhs doue to static condensation.
    /// -A_bi A_ii^{-1} 
    /// stored only for non-symmetric biforms
    virtual BaseMatrix & GetHarmonicExtensionTrans () const = 0;

    /// returns inverse of A_ii
    virtual BaseMatrix & GetInnerSolve () const = 0;

    /// returns A_ii
    virtual BaseMatrix & GetInnerMatrix () const = 0;

    /// is there a low-order biform ?
    bool HasLowOrderBilinearForm () const { return low_order_bilinear_form != NULL; }

    /// returns the low-order biform
    BilinearForm & GetLowOrderBilinearForm() const
    {
      return *low_order_bilinear_form;
    }


    /// use static condensation ?
    bool UsesEliminateInternal () const { return eliminate_internal; }

    /// stores the matrices for reconstructing internal dofs ?
    bool UsesKeepInternal () const { return keep_internal; }

    /// does it store Aii ?
    bool UsesStoreInner () const { return store_inner; }


    /// the finite element space
    const FESpace & GetFESpace() const { return fespace; }


    /// uses mixed spaces (non operational)
    bool MixedSpaces () const { return fespace2 != NULL; }

    /// returns the second space (form mixed spaces)
    const FESpace & GetFESpace2() const { return *fespace2; }

    ///
    int GetNLevels() const { return mats.Size(); }

    /// is the form symmetric ?
    bool IsSymmetric(void) const { return symmetric; }

    /// don't assemble the matrix
    void SetNonAssemble (bool na = true) { nonassemble = na; }

    ///
    void SetGalerkin (bool agalerkin = true) { galerkin = agalerkin; }

    ///
    void SetDiagonal (bool adiagonal = true) { diagonal = adiagonal; }

    ///
    void SetSymmetric (bool asymmetric = true) { symmetric = asymmetric; }

    ///
    void SetHermitean (bool ahermitean = true) { hermitean = ahermitean; }

    ///
    void SetMultiLevel (bool amultilevel = 1) { multilevel = amultilevel; }

    ///
    void SetTiming (bool at) { timing = at; }

    void SetEliminateInternal (bool eliminate) 
    { eliminate_internal = eliminate; }

    void SetKeepInternal (bool keep)
    { keep_internal = keep; }

    void SetStoreInner (bool storei) 
    { store_inner = storei; }

    void SetPrint (bool ap);
    void SetPrintElmat (bool ap);
    void SetElmatEigenValues (bool ee);

    /// computes low-order matrices from fines matrix
    void GalerkinProjection ();

    /// reconstruct internal dofs
    virtual void ComputeInternal (BaseVector & u, const BaseVector & f, LocalHeap & lh) const = 0;

    /// modify rhs due to static condensation
    virtual void ModifyRHS (BaseVector & f) const = 0;
  

  
    /// add eps I to the assembled matrix
    void SetEpsRegularization(double val) { eps_regularization = val; }

    /// set matrix diagonal for unused dofs to this value
    void SetUnusedDiag (double val) { unuseddiag = val; }

    /// does it use Galerkin projection ?
    bool UseGalerkin () const { return galerkin; }

    /// biform object
    virtual string GetClassName () const
    {
      return "BilinearForm";
    }

    /// prints report to file
    virtual void PrintReport (ostream & ost) const;

    ///
    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;

    /// creates a compatible vector
    virtual shared_ptr<BaseVector> CreateVector() const = 0;

    /// frees matrix 
    virtual void CleanUpLevel() { ; }

  private:
    /// assemble matrix
    virtual void DoAssemble (LocalHeap & lh) = 0;

    /// allocates (sparse) matrix data-structure
    virtual void AllocateMatrix () = 0;
  };






  /**
     We specify the scalar (double or Complex) of the biform.
   */
  template <class SCAL>
  class NGS_DLL_HEADER S_BilinearForm : public BilinearForm
  {
  protected:

    ElementByElementMatrix<SCAL> * harmonicext;
    BaseMatrix * harmonicexttrans;
    ElementByElementMatrix<SCAL> * innersolve;
    ElementByElementMatrix<SCAL> * innermatrix;

        
  public:
    /// 
    S_BilinearForm (const FESpace & afespace, const string & aname,
		    const Flags & flags)
      : BilinearForm (afespace, aname, flags) 
    { 
      harmonicext = NULL;
      harmonicexttrans = NULL;
      innersolve = NULL;
      innermatrix = NULL;
    }

    ///
    S_BilinearForm (const FESpace & afespace, 
		    const FESpace & afespace2,
		    const string & aname, const Flags & flags)
      : BilinearForm (afespace, afespace2, aname, flags) 
    {
      harmonicext = NULL;
      harmonicexttrans = NULL;
      innersolve = NULL;
      innermatrix = NULL;
    }


    ~S_BilinearForm();

    ///
    void AddMatrix1 (SCAL val, const BaseVector & x,
		     BaseVector & y) const;

    virtual void AddMatrix (double val, const BaseVector & x,
			    BaseVector & y) const
    {
      AddMatrix1 (val, x, y);
    }


    virtual void AddMatrix (Complex val, const BaseVector & x,
			    BaseVector & y) const
    {
      AddMatrix1 (ConvertTo<SCAL> (val), x, y);
    }


    virtual void LapackEigenSystem(FlatMatrix<SCAL> & elmat, LocalHeap & lh) const 
    {
      ;
    }
  
    void ApplyLinearizedMatrixAdd1 (SCAL val,
				    const BaseVector & lin,
				    const BaseVector & x,
				    BaseVector & y) const;
  
    virtual void ApplyLinearizedMatrixAdd (double val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y) const
    {
      ApplyLinearizedMatrixAdd1 (val, lin, x, y);
    }
  
    virtual void ApplyLinearizedMatrixAdd (Complex val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y) const
    {
      ApplyLinearizedMatrixAdd1 (ConvertTo<SCAL> (val), lin, x, y);
    }
  

    virtual double Energy (const BaseVector & x) const;

    virtual void ComputeInternal (BaseVector & u, const BaseVector & f, LocalHeap & lh) const;

    virtual void ModifyRHS (BaseVector & fd) const;

    ///
    virtual void DoAssemble (LocalHeap & lh);
    ///
    // virtual void DoAssembleIndependent (BitArray & useddof, LocalHeap & lh);
    ///
    virtual void AssembleLinearization (const BaseVector & lin,
					LocalHeap & lh, 
					bool reallocate = 0);
    ///
    virtual void AddElementMatrix (FlatArray<int> dnums1,
                                   FlatArray<int> dnums2,
                                   FlatMatrix<SCAL> elmat,
				   ElementId id, 
				   LocalHeap & lh) = 0;

    virtual void ApplyElementMatrix(const BaseVector & x,
				    BaseVector & y,
				    const SCAL & val,
				    const Array<int> & dnums,
				    const ElementTransformation & eltrans,
				    const int elnum,
				    const int type,
				    int & cnt,
				    LocalHeap & lh,
				    const FiniteElement * fel,
				    const SpecialElement * sel = NULL) const
    { cerr << "ApplyElementMatrix called for baseclass" << endl;}

    virtual void AddDiagElementMatrix (const Array<int> & dnums1,
				       const FlatVector<SCAL> & diag,
				       bool inner_element, int elnr,
				       LocalHeap & lh);


    BaseMatrix & GetHarmonicExtension () const 
    { 
      return *harmonicext; 
    }
    ///  
    BaseMatrix & GetHarmonicExtensionTrans () const
    { 
      return *harmonicexttrans; 
    }
    ///  
    BaseMatrix & GetInnerSolve () const
    { 
      return *innersolve; 
    }
    ///  
    BaseMatrix & GetInnerMatrix () const
    { 
      return *innermatrix; 
    }


  };



  template <class TM, class TV = typename mat_traits<TM>::TV_COL>
  class NGS_DLL_HEADER T_BilinearForm : public S_BilinearForm<typename mat_traits<TM>::TSCAL>
  {
  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV TV_COL;
    typedef SparseMatrix<TM,TV,TV> TMATRIX;
    
  protected:

  public:
    ///
    T_BilinearForm (const FESpace & afespace, const string & aname, const Flags & flags);
    ///
    T_BilinearForm (const FESpace & afespace, 
		    const FESpace & afespace2,
		    const string & aname,
		    const Flags & flags);
    ///
    virtual ~T_BilinearForm ();

    ///
    virtual void AllocateMatrix ();

    ///
    virtual shared_ptr<BaseVector> CreateVector() const;

    ///
    virtual void CleanUpLevel();

    ///
    virtual void AddElementMatrix (FlatArray<int> dnums1,
				   FlatArray<int> dnums2,
				   FlatMatrix<TSCAL> elmat,
				   ElementId id, 
				   LocalHeap & lh);

    virtual void ApplyElementMatrix(const BaseVector & x,
				    BaseVector & y,
				    const TSCAL & val,
				    const Array<int> & dnums,
				    const ElementTransformation & eltrans,
				    const int elnum,
				    const int type,
				    int & cnt,
				    LocalHeap & lh,
				    const FiniteElement * fel,
				    const SpecialElement * sel = NULL) const;

    virtual void LapackEigenSystem(FlatMatrix<TSCAL> & elmat, LocalHeap & lh) const;
  };






  template <class TM, class TV = typename mat_traits<TM>::TV_COL>
  class NGS_DLL_HEADER T_BilinearFormSymmetric : public S_BilinearForm<typename mat_traits<TM>::TSCAL>
  {

  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV TV_COL;
    typedef SparseMatrixSymmetric<TM,TV> TMATRIX;
    
  protected:
    

  public:
    T_BilinearFormSymmetric (const FESpace & afespace, const string & aname,
			     const Flags & flags);
    virtual ~T_BilinearFormSymmetric ();

    virtual void AllocateMatrix ();
    virtual void CleanUpLevel();

    virtual shared_ptr<BaseVector> CreateVector() const;

    virtual void AddElementMatrix (FlatArray<int> dnums1,
				   FlatArray<int> dnums2,
                                   FlatMatrix<TSCAL> elmat,
				   ElementId id, 
				   LocalHeap & lh);
    virtual void ApplyElementMatrix(const BaseVector & x,
				    BaseVector & y,
				    const TSCAL & val,
				    const Array<int> & dnums,
				    const ElementTransformation & eltrans,
				    const int elnum,
				    const int type,
				    int & cnt,
				    LocalHeap & lh,
				    const FiniteElement * fel,
				    const SpecialElement * sel = NULL) const;

    virtual void LapackEigenSystem(FlatMatrix<TSCAL> & elmat, LocalHeap & lh) const;
  };







  template <class TM>
  class NGS_DLL_HEADER T_BilinearFormDiagonal : public S_BilinearForm<typename mat_traits<TM>::TSCAL>
  {

  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef typename mat_traits<TM>::TV_COL TV_COL;
    typedef SparseMatrixSymmetric<TM> TMATRIX;

  protected:

  public:
    T_BilinearFormDiagonal (const FESpace & afespace, const string & aname,
			    const Flags & flags);
    virtual ~T_BilinearFormDiagonal ();

    virtual void AllocateMatrix ();
    virtual shared_ptr<BaseVector> CreateVector() const;

    virtual void AddElementMatrix (FlatArray<int> dnums1,
				   FlatArray<int> dnums2,
				   FlatMatrix<TSCAL> elmat,
				   ElementId id, 
				   LocalHeap & lh);

    virtual void AddDiagElementMatrix (const Array<int> & dnums1,
				       const FlatVector<TSCAL> & diag,
				       bool inner_element, int elnr,
				       LocalHeap & lh);
    virtual void ApplyElementMatrix(const BaseVector & x,
				    BaseVector & y,
				    const TSCAL & val,
				    const Array<int> & dnums,
				    const ElementTransformation & eltrans,
				    const int elnum,
				    const int type,
				    int & cnt,
				    LocalHeap & lh,
				    const FiniteElement * fel,
				    const SpecialElement * sel = NULL) const;
  };




  /**
     Allocates a bilinearform.
     Some flags are:
     -symmetric   ... assembles a symmetric matrix
   */
  extern NGS_DLL_HEADER shared_ptr<BilinearForm> CreateBilinearForm (shared_ptr<FESpace> space,
                                                                     const string & name,
                                                                     const Flags & flags);



  /**
     This objects provide the bilinear-form application as matrix vector product.
     If the bilinearform is indeed non-linear in the first argumen, the operator*
     will perform the non-linear operator application.
   */
  class NGS_DLL_HEADER BilinearFormApplication : public BaseMatrix
  {
  protected:
    ///
    const BilinearForm * bf;

  public:
    ///
    BilinearFormApplication (const BilinearForm * abf);
    ///
    virtual void Mult (const BaseVector & v, BaseVector & prod) const;
    ///
    virtual void MultAdd (double val, const BaseVector & v, BaseVector & prod) const;
    ///
    virtual void MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const;
    ///
    virtual shared_ptr<BaseVector> CreateVector () const;
    ///
    virtual int VHeight() const
    {
      return bf->GetFESpace().GetNDof(); 
    }
    ///
    virtual int VWidth() const
    {
      return bf->GetFESpace().GetNDof(); 
    }
  };

  /**
     Applies the matrix-vector product of linearized matrix.
     Linearization point is given in the constructor
   */
  class NGS_DLL_HEADER LinearizedBilinearFormApplication : public BilinearFormApplication 
  {
  protected:
    const BaseVector * veclin;

  public:
    ///
    LinearizedBilinearFormApplication (const BilinearForm * abf,
				       const BaseVector * aveclin);

    ///
    virtual void Mult (const BaseVector & v, BaseVector & prod) const;
    ///
    virtual void MultAdd (double val, const BaseVector & v, BaseVector & prod) const;
    ///
    virtual void MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const;
  };


  /**
     This bilinearform stores the element-matrices
   */
  template<class SCAL>
  class ElementByElement_BilinearForm : public S_BilinearForm<SCAL>
  {
  public:
    ElementByElement_BilinearForm (const FESpace & afespace, const string & aname,
		    const Flags & flags);
    virtual ~ElementByElement_BilinearForm ();
    
    virtual void AllocateMatrix ();
    virtual shared_ptr<BaseVector> CreateVector() const;
    
    virtual void AddElementMatrix (FlatArray<int> dnums1,
                                   FlatArray<int> dnums2,
                                   FlatMatrix<SCAL> elmat,
                                   ElementId id, 
                                   LocalHeap & lh);    
  };
  
  
  
}


#endif
