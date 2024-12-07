#ifndef FILE_BILINEARFORM
#define FILE_BILINEARFORM

/*********************************************************************/
/* File:   bilinearform.hpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#include "fespace.hpp"
#include <specialelement.hpp>
#include <sparsematrix.hpp>
#include <elementbyelement.hpp>


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
    shared_ptr<FESpace> fespace;
    /// Test-space if different from trial-space, otherwise NULL
    shared_ptr<FESpace> fespace2;

    /// don't assemble matrix
    bool nonassemble;
    /// store only diagonal of matrix
    bool diagonal = false;
    /// element-matrix for ref-elements
    bool geom_free = false;
    /// stores geom-free B factors, and D factors in integration points
    bool matrix_free_bdb = false;
    /// stores geom-free B factors, and D factors in integration points, and compiled CF
    bool nonlinear_matrix_free_bdb = false;
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
    /// check if all dofs declared used are used in assemble
    bool check_unused = true;
    /// low order bilinear-form, 0 if not used
    shared_ptr<BilinearForm> low_order_bilinear_form;

    /// modify linear form due to static condensation
    LinearForm * linearform;

    /// some preconditioners need element matrices 
    Array<Preconditioner*> preconditioners; 

    /// matrices (sparse, application, diagonal, ...)
    Array<shared_ptr<BaseMatrix>> mats;
    size_t graph_timestamp = 0;
    
    /// bilinearform-integrators
    Array<shared_ptr<BilinearFormIntegrator>> parts;
    Array<shared_ptr<BilinearFormIntegrator>> VB_parts[4];

    // loop over facets, VB=0 .. inner facets, VB=1 .. boundary facets
    Array<shared_ptr<FacetBilinearFormIntegrator>> facetwise_skeleton_parts[2];

    // geometry-free parts (only for apply)
    Array<shared_ptr<BilinearFormIntegrator>> geom_free_parts;
    
    // loop over elements
    Array<shared_ptr<FacetBilinearFormIntegrator>> elementwise_skeleton_parts;

    // #ifdef PARALLEL
    Array<shared_ptr<FacetBilinearFormIntegrator> > mpi_facet_parts;
    // #endif

    /// special elements for hacks (used for contact, periodic-boundary-penalty-constraints, ...
    Array<unique_ptr<SpecialElement>> specialelements;
    mutable unique_ptr<Table<int>> special_element_coloring;
    
    size_t specialelements_timestamp = 0;

    
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
    /// does static condensation of hidden dofs
    bool eliminate_hidden;
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
    ///
    optional<double> delete_zero_elements;
    
    mutable std::map<size_t, Matrix<>> precomputed;
  public:
    /// generate a bilinear-form
    BilinearForm (shared_ptr<FESpace> afespace,
		  const string & aname, 
		  const Flags & flags);

    /// generate a bilinear-form
    BilinearForm (shared_ptr<FESpace> afespace, 
		  shared_ptr<FESpace> afespace2, 
		  const string & aname,
		  const Flags & flags);
    BilinearForm (const BilinearForm&) = delete;
    BilinearForm& operator= (const BilinearForm&) = delete;
    virtual ~BilinearForm ();
  

    ///
    virtual BilinearForm & AddIntegrator (shared_ptr<BilinearFormIntegrator> bfi);

    BilinearForm & operator+= (shared_ptr<BilinearFormIntegrator> bfi)
    {
      return AddIntegrator(bfi);
    }

    /*
    void AddIndependentIntegrator (BilinearFormIntegrator * bfi,
				   const int master_surface,
				   const int other)
    */
  
    const Array<shared_ptr<BilinearFormIntegrator>> & Integrators() const
    {
      return parts;
    }

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
    int GetIndependentOtherIndex(const int i) const
    {
      return independent_meshindex[i](1);
    }
    */

    void AddSpecialElement (unique_ptr<SpecialElement> spel);
    auto & GetSpecialElements() const { return specialelements; }
    void DeleteSpecialElement(size_t index);
    void DeleteSpecialElements();
    Table<int> & SpecialElementColoring() const;
    
    /// for static condensation of internal bubbles
    void SetLinearForm (LinearForm * alf) { linearform = alf; }
    
    /// preconditioner gets element-matrix
    void SetPreconditioner (Preconditioner * pre);

    /// unregister preconditioner
    void UnsetPreconditioner(Preconditioner* pre);

    /// generates matrix graph
    virtual MatrixGraph GetGraph (int level, bool symmetric);

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
		      BaseVector & y, LocalHeap & lh) const
    {
      x.Cumulate();
      y = 0;
      AddMatrix (1, x, y, lh);
      y.SetParallelStatus(DISTRIBUTED);
    }

    /// y += val * Mat * x
    virtual void AddMatrix (double val, const BaseVector & x,
			    BaseVector & y, LocalHeap & lh) const = 0;
  
    /// y += val * Mat * x
    virtual void AddMatrix (Complex val, const BaseVector & x,
			    BaseVector & y, LocalHeap & lh) const = 0;

    /// y += val * Mat^T * x
    virtual void AddMatrixTrans (double val, const BaseVector & x,
                                 BaseVector & y, LocalHeap & lh) const = 0;
    
    /// y += val * lin.mat * x
    virtual void ApplyLinearizedMatrixAdd (double val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y, LocalHeap & lh) const = 0;
    /// y += val * lin.mat * x
    virtual void ApplyLinearizedMatrixAdd (Complex val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y, LocalHeap & lh) const = 0;

    /// evaulates internal energy (usually  1/2 x^T A x)
    virtual double Energy (const BaseVector & x, LocalHeap & lh) const = 0;

    /// returns the assembled matrix
    const BaseMatrix & GetMatrix () const { return *mats.Last(); }
    const BaseMatrix & GetMatrix (int level) const { return *mats[level]; }
    void DeleteMatrix()
    {
      if(mats.Size())
        mats.DeleteLast();
    }
    /// returns the assembled matrix
    shared_ptr<BaseMatrix> GetMatrixPtr () const;

    // operator const BaseMatrix& () const { return GetMatrix(); }

    /// returns the assembled matrix on a given level
    shared_ptr<BaseMatrix> GetMatrixPtr (int level) const { return mats[level]; }

    BaseMatrix & GetMatrix ()  { return *mats.Last(); }
    BaseMatrix & GetMatrix (int level)  { return *mats[level]; }

    /// reconstruct internal dofs from static condensation 
    /// -A_ii^{-1} A_ib
    virtual shared_ptr<BaseMatrix> GetHarmonicExtension () const = 0;

    /// modify rhs doue to static condensation.
    /// -A_bi A_ii^{-1} 
    /// stored only for non-symmetric biforms
    virtual shared_ptr<BaseMatrix> GetHarmonicExtensionTrans () const = 0;

    /// returns inverse of A_ii
    virtual shared_ptr<BaseMatrix> GetInnerSolve () const = 0;

    /// returns A_ii
    virtual shared_ptr<BaseMatrix> GetInnerMatrix () const = 0;

    /// is there a low-order biform ?
    bool HasLowOrderBilinearForm () const { return low_order_bilinear_form != NULL; }

    /// returns the low-order biform, if we can provide it ...
    virtual shared_ptr<BilinearForm> GetLowOrderBilinearForm() 
    {
      return low_order_bilinear_form;
    }


    /// use static condensation ?
    bool UsesEliminateInternal () const { return eliminate_internal; }

    /// use static condensation for hidden?
    bool UsesEliminateHidden () const { return eliminate_hidden; }

    /// stores the matrices for reconstructing internal dofs ?
    bool UsesKeepInternal () const { return keep_internal; }

    /// does it store Aii ?
    bool UsesStoreInner () const { return store_inner; }


    /// the finite element space
    // const FESpace & GetFESpace() const { return *fespace; }


    /// uses mixed spaces (non operational)
    bool MixedSpaces () const { return fespace2 != NULL; }

    /// returns the second space (form mixed spaces)
    // const FESpace & GetFESpace2() const { return *fespace2; }

    /// the finite element space
    shared_ptr<FESpace> GetFESpace() const { return fespace; }
    /// the finite element test space
    shared_ptr<FESpace> GetFESpace2() const { return fespace2; }
    
    shared_ptr<FESpace> GetTrialSpace() const { return fespace; }
    shared_ptr<FESpace> GetTestSpace() const { return fespace2 ? fespace2 : fespace; }
    ///
    int GetNLevels() const { return mats.Size(); }

    /// is the form symmetric ?
    bool IsSymmetric() const { return symmetric; }

    /// is the form symmetric and positive definite ?
    bool IsSPD() const { return spd; }

    ///
    virtual bool SymmetricStorage() const { return false; }

    /// don't assemble the matrix
    void SetNonAssemble (bool na = true) { nonassemble = na; }
    bool NonAssemble() const { return nonassemble; }

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

    void SetEliminateHidden (bool eliminate) 
    { eliminate_hidden = eliminate; }

    void SetKeepInternal (bool keep)
    { keep_internal = keep; }

    void SetStoreInner (bool storei) 
    { store_inner = storei; }

    void SetPrint (bool ap);
    void SetPrintElmat (bool ap);
    void SetElmatEigenValues (bool ee);
    void SetCheckUnused (bool b);
    
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
    virtual Array<MemoryUsage> GetMemoryUsage () const;

    /// creates a compatible vector
    virtual AutoVector CreateRowVector() const;
    virtual AutoVector CreateColVector() const;

    /// frees matrix 
    virtual void CleanUpLevel() { ; }

  protected:
    /// assemble matrix
    virtual void DoAssemble (LocalHeap & lh) = 0;
    void AssembleGF (LocalHeap & lh);
    void AssembleBDB (LocalHeap & lh, bool linear);

    /// allocates (sparse) matrix data-structure
    virtual void AllocateMatrix () = 0;
    virtual void AllocateInternalMatrices () = 0;
  };



  


  /**
     We specify the scalar (double or Complex) of the biform.
   */
  template <class SCAL>
  class NGS_DLL_HEADER S_BilinearForm : public BilinearForm
  {
  protected:

    // have parallel wrapper for distributed meshes:
    shared_ptr<BaseMatrix> harmonicext;
    shared_ptr<BaseMatrix> harmonicexttrans; 
    shared_ptr<BaseMatrix> innersolve;
    shared_ptr<BaseMatrix> innermatrix;

    // local operators:
    ElementByElementMatrix<SCAL> *harmonicext_ptr, *harmonicexttrans_ptr, *innersolve_ptr, *innermatrix_ptr;

    
    //data for mpi-facets; only has data if there are relevant integrators in the BLF!
    mutable bool have_mpi_facet_data = false;
    mutable Array<int> os_per;
    mutable Table<SCAL> send_table;
    mutable Table<SCAL> recv_table;
    
        
  public:
    ///
    using BilinearForm::BilinearForm;
    
    virtual ~S_BilinearForm();

    ///
    void AddMatrix1 (SCAL val, const BaseVector & x,
                    BaseVector & y, LocalHeap & lh) const;

    void AddMatrixGF (SCAL val, const BaseVector & x,
                      BaseVector & y, bool transpose, LocalHeap & lh) const;

    virtual void AddMatrix (double val, const BaseVector & x,
                           BaseVector & y, LocalHeap & lh) const override
    {
      x.Cumulate();
      y.Distribute();

      AddMatrix1 (val, x, y, lh);
    }

    virtual void AddMatrixTP (SCAL val, const BaseVector & x,
                             BaseVector & y, LocalHeap & lh) const;

    virtual void AddMatrix (Complex val, const BaseVector & x,
                           BaseVector & y, LocalHeap & lh) const override
    {
      x.Cumulate();
      y.Distribute();

      // AddMatrix1 (ConvertTo<SCAL> (val), x, y, lh);
      if constexpr (std::is_constructible<SCAL,Complex>())
        AddMatrix1 (SCAL(val), x, y, lh);
      else
        throw Exception("BilinearForm::AddMatrix(complex) called for real BilinearForm");
    }

    virtual void AddMatrixTrans (double val, const BaseVector & x,
                                 BaseVector & y, LocalHeap & lh) const override;

    virtual void LapackEigenSystem(FlatMatrix<SCAL> & elmat, LocalHeap & lh) const;
    // { ; }
  
    void ApplyLinearizedMatrixAdd1 (SCAL val,
				    const BaseVector & lin,
				    const BaseVector & x,
				    BaseVector & y, LocalHeap & lh) const;
  
    virtual void ApplyLinearizedMatrixAdd (double val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y, LocalHeap & lh) const override
    {
      lin.Cumulate();
      x.Cumulate();
      y.Distribute();

      ApplyLinearizedMatrixAdd1 (val, lin, x, y, lh);
    }
  
    virtual void ApplyLinearizedMatrixAdd (Complex val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y, LocalHeap & lh) const override
    {
      lin.Cumulate();
      x.Cumulate();
      y.Distribute();

      if constexpr (std::is_constructible<SCAL,Complex>())
        ApplyLinearizedMatrixAdd1 (SCAL(val), lin, x, y, lh);
      else
        throw Exception("BilinearForm::ApplyLinearizedMatrix(complex) called for real BilinearForm");
    }
  

    virtual double Energy (const BaseVector & x, LocalHeap & lh) const override;

    virtual void ComputeInternal (BaseVector & u, const BaseVector & f, LocalHeap & lh) const override;

    virtual void ModifyRHS (BaseVector & fd) const override;

    ///
    virtual void DoAssemble (LocalHeap & lh) override;
    virtual void Assemble_facetwise_skeleton_parts_VOL (Array<bool>& useddof, size_t & gcnt, LocalHeap & lh, const BaseVector * lin = nullptr);
    ///
    // virtual void DoAssembleIndependent (BitArray & useddof, LocalHeap & lh);
    ///
    virtual void AssembleLinearization (const BaseVector & lin,
					LocalHeap & lh, 
					bool reallocate = 0) override;
    ///
    virtual void AddElementMatrix (FlatArray<int> dnums1,
                                   FlatArray<int> dnums2,
                                   BareSliceMatrix<SCAL> elmat,
				   ElementId id, bool addatomic,
				   LocalHeap & lh) = 0;


    virtual void AddDiagElementMatrix (FlatArray<int> dnums1,
                                       FlatVector<SCAL> diag,
				       bool inner_element, int elnr,
				       LocalHeap & lh);


    shared_ptr<BaseMatrix> GetHarmonicExtension () const override
    { 
      return harmonicext; 
    }
    ///  
    shared_ptr<BaseMatrix> GetHarmonicExtensionTrans () const override
    { 
      return harmonicexttrans; 
    }
    ///  
    shared_ptr<BaseMatrix> GetInnerSolve () const override
    { 
      return innersolve; 
    }
    ///  
    shared_ptr<BaseMatrix> GetInnerMatrix () const override
    { 
      return innermatrix; 
    }

    virtual void AllocateInternalMatrices () override;
  };



  template <class TM, class TV = typename mat_traits<TM>::TV_COL>
  class NGS_DLL_HEADER T_BilinearForm : public S_BilinearForm<typename mat_traits<TM>::TSCAL>
  {
  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV TV_COL;
    typedef SparseMatrix<TM,TV,TV> TMATRIX;
    shared_ptr<TMATRIX> mymatrix;
  protected:

  public:
    using S_BilinearForm<TSCAL> :: S_BilinearForm;
    ///
    virtual ~T_BilinearForm () { };

    ///
    virtual shared_ptr<BilinearForm> GetLowOrderBilinearForm() override;

    virtual void AllocateMatrix () override;

    ///
    virtual void CleanUpLevel() override;

    ///
    virtual void AddElementMatrix (FlatArray<int> dnums1,
				   FlatArray<int> dnums2,
				   BareSliceMatrix<TSCAL> elmat,
				   ElementId id, bool addatomic, 
				   LocalHeap & lh) override;
  };






  template <class TM, class TV = typename mat_traits<TM>::TV_COL>
  class NGS_DLL_HEADER T_BilinearFormSymmetric : public S_BilinearForm<typename mat_traits<TM>::TSCAL>
  {

  public:
    typedef typename mat_traits<TM>::TSCAL TSCAL;
    typedef TV TV_COL;
    typedef SparseMatrixSymmetric<TM,TV> TMATRIX;
    shared_ptr<TMATRIX> mymatrix;    
  protected:
    

  public:
    T_BilinearFormSymmetric (shared_ptr<FESpace> afespace, const string & aname,
			     const Flags & flags);
    virtual ~T_BilinearFormSymmetric () override;

    virtual void AllocateMatrix () override;
    virtual void CleanUpLevel() override;
    virtual shared_ptr<BilinearForm> GetLowOrderBilinearForm() override;
    
    virtual void AddElementMatrix (FlatArray<int> dnums1,
				   FlatArray<int> dnums2,
                                   BareSliceMatrix<TSCAL> elmat,
				   ElementId id, bool addatomic, 
				   LocalHeap & lh) override;

    virtual bool SymmetricStorage() const override { return true; }
    virtual void LapackEigenSystem(FlatMatrix<TSCAL> & elmat, LocalHeap & lh) const override;
  };



  template <class TSCAL>
  class NGS_DLL_HEADER T_BilinearFormDynBlocks : public S_BilinearForm<TSCAL>
  {
  public:
    typedef SparseBlockMatrix<TSCAL> TMATRIX;
    shared_ptr<TMATRIX> mymatrix;
    size_t blockheight, blockwidth;
  protected:

  public:
    T_BilinearFormDynBlocks (shared_ptr<FESpace> afespace, 
                             const string & aname, const Flags & flags)
      : S_BilinearForm<TSCAL> (afespace, aname, flags),
      blockheight(afespace->GetDimension()),
      blockwidth(afespace->GetDimension()) { } 
    
    T_BilinearFormDynBlocks (shared_ptr<FESpace> afespace, 
                             shared_ptr<FESpace> afespace2,
                             const string & aname, const Flags & flags)
      : S_BilinearForm<TSCAL> (afespace, afespace2, aname, flags),
      blockheight(afespace2->GetDimension()),
      blockwidth(afespace->GetDimension()) { } 
    
    virtual ~T_BilinearFormDynBlocks () { };

    virtual shared_ptr<BilinearForm> GetLowOrderBilinearForm() override;

    virtual void AllocateMatrix () override;

    virtual void CleanUpLevel() override;

    virtual void AddElementMatrix (FlatArray<int> dnums1,
				   FlatArray<int> dnums2,
				   BareSliceMatrix<TSCAL> elmat,
				   ElementId id, bool addatomic, 
				   LocalHeap & lh) override;
  };






  

  template <class TSCAL>
  class NGS_DLL_HEADER S_BilinearFormNonAssemble : public S_BilinearForm<TSCAL>
  {
  public:
    S_BilinearFormNonAssemble (shared_ptr<FESpace> afespace, const string & aname,
                               const Flags & flags);
    S_BilinearFormNonAssemble (shared_ptr<FESpace> afespace, shared_ptr<FESpace> afespace2,
                               const string & aname, const Flags & flags);
    // virtual ~T_BilinearFormSymmetric ();

    virtual void AllocateMatrix () { cout << "S_BilinearFormNonAssemble :: Allocate: nothing to do" << endl; }
    virtual void CleanUpLevel() { ; } 

    virtual void AddElementMatrix (FlatArray<int> dnums1,
				   FlatArray<int> dnums2,
                                   BareSliceMatrix<TSCAL> elmat,
				   ElementId id, bool addatomic,
				   LocalHeap & lh)
    {
      throw Exception ("AddElementMatrix for non-assemble biform called");
    }

    virtual bool SymmetricStorage() const { return true; }

    virtual void LapackEigenSystem(FlatMatrix<TSCAL> & elmat, LocalHeap & lh) const
    { cout << "no eigensystem available" << endl; }
  };







  class ComponentBilinearForm : public BilinearForm
  {
    shared_ptr<BilinearForm> base_blf;
    int comp; // , ncomp;
  public:
    ComponentBilinearForm (shared_ptr<BilinearForm> abase_blf, int acomp, int ancomp);
    virtual BilinearForm & AddIntegrator (shared_ptr<BilinearFormIntegrator> bfi);

    virtual void Assemble (LocalHeap & lh) { cerr << "comp - assemble is illegal" << endl; }

    virtual void AssembleLinearization (const BaseVector & lin,
					LocalHeap & lh, 
					bool reallocate = 0) 
    { throw Exception ("comp-bf - AssembleLinearization is illegal"); }

    virtual void AddMatrix (double val, const BaseVector & x,
			    BaseVector & y, LocalHeap & lh) const
      { throw Exception ("comp-bf - AddMatrix is illegal"); }

    virtual void AddMatrix (Complex val, const BaseVector & x,
			    BaseVector & y, LocalHeap & lh) const
    { throw Exception ("comp-bf - AddMatrix is illegal"); }

    virtual void AddMatrixTrans (double val, const BaseVector & x,
                                 BaseVector & y, LocalHeap & lh) const
    { throw Exception ("comp-bf - AddMatrixTrans is illegal"); }
    
    virtual void ApplyLinearizedMatrixAdd (double val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y, LocalHeap & lh) const 
    { throw Exception ("comp-bf - AddMatrix is illegal"); }

    virtual void ApplyLinearizedMatrixAdd (Complex val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y, LocalHeap & lh) const 
    { throw Exception ("comp-bf - AddMatrix is illegal"); }

    virtual shared_ptr<BaseMatrix> GetHarmonicExtension () const 
    { throw Exception ("comp-bf - GetHarmonicExt is illegal"); }

    virtual shared_ptr<BaseMatrix> GetHarmonicExtensionTrans () const
    { throw Exception ("comp-bf - GetHarmonicExtTrans is illegal"); } 
    virtual shared_ptr<BaseMatrix> GetInnerSolve () const 
    { throw Exception ("comp-bf - GetInnerSolve is illegal"); } 
    virtual shared_ptr<BaseMatrix> GetInnerMatrix () const
    { throw Exception ("comp-bf - GetInnerMatrix is illegal"); } 
    virtual void ComputeInternal (BaseVector & u, const BaseVector & f, LocalHeap & lh) const
    { throw Exception ("comp-bf - ComputeInternal is illegal"); } 
    virtual void ModifyRHS (BaseVector & f) const 
    { throw Exception ("comp-bf - ModifyRHS is illegal"); } 
    virtual AutoVector CreateRowVector() const 
    { throw Exception ("comp-bf - CreateRowVector is illegal"); } 
    virtual AutoVector CreateColVector() const 
    { throw Exception ("comp-bf - CreateColVector is illegal"); } 
    virtual void DoAssemble (LocalHeap & lh) 
    { throw Exception ("comp-bf - DoAssemble is illegal"); } 
    virtual void AllocateMatrix ()
    { throw Exception ("comp-bf - AllocateMatrix is illegal"); } 
    virtual void AllocateInternalMatrices ()
    { throw Exception ("comp-bf - AllocateInternalMatrices is illegal"); } 
    virtual double Energy (const BaseVector & x, LocalHeap & lh) const 
    { throw Exception ("comp-bf - Energy is illegal"); } 

    /*
    virtual shared_ptr<BaseVector> GetVectorPtr() const
    { throw Exception ("comp - GetVectorPtr is illegal"); }
    virtual BaseVector & GetVector () const 
    { throw Exception ("comp - GetVector is illegal"); }
    */
  };






  /**
     Allocates a bilinearform.
     Some flags are:
     -symmetric   ... assembles a symmetric matrix
   */
  extern NGS_DLL_HEADER shared_ptr<BilinearForm> CreateBilinearForm (shared_ptr<FESpace> space,
                                                                     const string & name,
                                                                     const Flags & flags);

  extern NGS_DLL_HEADER shared_ptr<BilinearForm> CreateBilinearForm (shared_ptr<FESpace> space,
                                                                     shared_ptr<FESpace> space2,
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
    shared_ptr<BilinearForm> bf;
    LocalHeap & lh;
  public:
    ///
    BilinearFormApplication (shared_ptr<BilinearForm> abf, LocalHeap & alh);

    virtual bool IsComplex() const override
    {
      return bf->GetFESpace()->IsComplex();
    }
    
    ///
    virtual void Mult (const BaseVector & v, BaseVector & prod) const override;
    ///
    virtual void MultAdd (double val, const BaseVector & v, BaseVector & prod) const override;
    ///
    virtual void MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const override;
    ///
    virtual void MultTransAdd (double val, const BaseVector & v, BaseVector & prod) const override;
    
    virtual AutoVector CreateVector () const override;
    virtual AutoVector CreateRowVector () const override;
    virtual AutoVector CreateColVector () const override;
    
    ///
    virtual int VHeight() const override
    {
      if (bf->GetFESpace2())
        return bf->GetFESpace2()->GetNDof();         
      return bf->GetFESpace()->GetNDof(); 
    }
    ///
    virtual int VWidth() const override
    {
      return bf->GetFESpace()->GetNDof(); 
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
    LinearizedBilinearFormApplication (shared_ptr<BilinearForm> abf,
				       const BaseVector * aveclin,
                                       LocalHeap & alh);

    ///
    virtual void Mult (const BaseVector & v, BaseVector & prod) const override;
    ///
    virtual void MultAdd (double val, const BaseVector & v, BaseVector & prod) const override;
    ///
    virtual void MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const override;
  };



  class ApplyIntegrationPoints : public BaseMatrix
  {
    Array<shared_ptr<CoefficientFunction>> coefs;
    Array<ProxyFunction*> trialproxies;
    
    typedef void (*lib_function)(size_t nip, double * input, size_t dist_input,
                                 double * output, size_t dist_output,
                                 size_t dist, double * points, double * normals);

    unique_ptr<SharedLibrary> library;
    lib_function compiled_function = nullptr;
    
    size_t dimx, dimy;
    size_t nip;
    Matrix<double> points;
    Matrix<double> normals;
    
  public:
    ApplyIntegrationPoints (Array<shared_ptr<CoefficientFunction>> acoefs,
                            const Array<ProxyFunction*> & atrialproxies,
                            Matrix<double> apoints, Matrix<double> anormals,
                            size_t adimx, size_t adimy, size_t anip);
    
    AutoVector CreateColVector() const override;
    AutoVector CreateRowVector() const override;
    
    virtual int VHeight() const override { return nip*dimy; }
    virtual int VWidth() const override { return nip*dimx; }
    
    virtual void Mult (const BaseVector & x, BaseVector & y) const override;

    const Array<shared_ptr<CoefficientFunction>> & GetCFs() const { return coefs; } 
    const Array<ProxyFunction*> & GetTrialProxies() const { return trialproxies; }
    size_t GetDimX() const { return dimx; }
    size_t GetDimY() const { return dimy; }
    size_t GetNIP() const { return nip; }
  };  
  

  
  
  
}

#endif
