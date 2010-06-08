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
    /// Test-space if different from trial-space, otherwise 0
    const FESpace * fespace2;

    ///
    bool nonassemble;
    /// 
    bool diagonal;
    ///
    bool multilevel;
    /// galerkin projection of coarse grid matrices
    bool galerkin;
    /// complex forms are hermitean
    bool hermitean;
    /// bilinear form is symmetric
    bool symmetric;
    /// add epsilon for regularization
    double eps_regularization; 
    /// diagonal value for unused dofs
    double unuseddiag;

    /// low order bilinear-form, 0 if not used
    BilinearForm * low_order_bilinear_form;

    /// modify linear form due to static condensation
    LinearForm * linearform;

    /// matrices (sparse, application, diagonal, ...)
    Array<BaseMatrix*> mats;
    /// bilinearform-integrators
    Array<BilinearFormIntegrator*> parts;
    /// is biform responsible for the deallocation ?
    Array<bool> parts_deletable;
    ///
    Array<BilinearFormIntegrator*> independent_parts;
    /// is biform responsible for the deallocation ?
    Array<bool> independent_parts_deletable;
    ///
    Array< Vec<2,int> > independent_meshindex;

    ///
    bool timing;
    bool print;
    bool printelmat;
    bool elmat_ev;
    bool eliminate_internal;
  

    bool precompute;
    Array<void*> precomputed_data;


  public:
    BilinearForm (const FESpace & afespace,
		  const string & aname, 
		  const Flags & flags);

    ///
    BilinearForm (const FESpace & afespace, 
		  const FESpace & afespace2, 
		  const string & aname,
		  const Flags & flags);
    ///
    virtual ~BilinearForm ();
  

    ///
    void AddIntegrator (BilinearFormIntegrator * bfi, const bool deletable = true);

    ///
    void AddIndependentIntegrator (BilinearFormIntegrator * bfi,
				   const int master_surface,
				   const int slave,
				   const bool deletable = true);
  
    ///
    int NumIntegrators () const 
    {
      return parts.Size(); 
    }

    int NumIndependentIntegrators(void) const
    {
      return independent_parts.Size();
    }

    ///
    BilinearFormIntegrator * GetIntegrator (const int i) 
    { return parts[i]; }

    ///
    const BilinearFormIntegrator * GetIntegrator (const int i) const 
    { return parts[i]; }

    ///
    BilinearFormIntegrator * GetIndependentIntegrator (const int i) 
    { return independent_parts[i]; }

    ///
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


    /// for static condensation of internal bubbles
    void SetLinearForm (LinearForm * alf)
    { linearform = alf; }

    ///
    void Assemble (LocalHeap & lh);
    ///
    void ReAssemble (LocalHeap & lh, bool reallocate = 0);
    ///
    virtual void AssembleLinearization (const BaseVector & lin,
					LocalHeap & lh, 
					bool reallocate = 0) = 0;

    ///
    void ApplyMatrix (const BaseVector & x,
		      BaseVector & y) const
    {
      y = 0;
      AddMatrix (1, x, y);
    }

    ///
    virtual void AddMatrix (double val, const BaseVector & x,
			    BaseVector & y) const = 0;
  
    virtual void AddMatrix (Complex val, const BaseVector & x,
			    BaseVector & y) const = 0;
  
    ///
    virtual void ApplyLinearizedMatrixAdd (double val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y) const = 0;
    ///
    virtual void ApplyLinearizedMatrixAdd (Complex val,
					   const BaseVector & lin,
					   const BaseVector & x,
					   BaseVector & y) const = 0;


    virtual double Energy (const BaseVector & x) const = 0;

    ///
    const BaseMatrix & GetMatrix () const
    { 
      return *mats.Last(); 
    }

    ///
    const BaseMatrix & GetMatrix (int level) const
    { 
      return *mats[level];
    }

    ///  
    BaseMatrix & GetMatrix () 
    { 
      return *mats.Last(); 
    }

    ///
    BaseMatrix & GetMatrix (int level) 
    { 
      return *mats[level];
    }


    bool HasLowOrderBilinearForm(void) const {return low_order_bilinear_form != NULL;}
    bool UsesEliminateInternal(void) const {return eliminate_internal;}

    const BilinearForm & GetLowOrderBilinearForm() const
    {
      return *low_order_bilinear_form;
    }

  
    BilinearForm & GetLowOrderBilinearForm() 
    {
      return *low_order_bilinear_form;
    }

  

    ///
    const FESpace & GetFESpace() const
    { return fespace; }
    ///
    int MixedSpaces () const
    { return fespace2 != NULL; }
    ///
    const FESpace & GetFESpace2() const
    { return *fespace2; }

    ///
    int GetNLevels() const
    { return mats.Size(); }

    bool IsSymmetric(void) const {return symmetric;}


    void SetNonAssemble (bool na = true)
    { nonassemble = na; }

    ///
    void SetGalerkin (bool agalerkin = true)
    { galerkin = agalerkin; }
    ///
    void SetDiagonal (bool adiagonal = true)
    { diagonal = adiagonal; }
    ///
    void SetSymmetric (bool asymmetric = true)
    { symmetric = asymmetric; }
    ///
    void SetHermitean (bool ahermitean = true)
    { hermitean = ahermitean; }
    ///
    void SetMultiLevel (bool amultilevel = 1)
    { multilevel = amultilevel; }

    void SetTiming (bool at) 
    { timing = at; }

    void SetEliminateInternal (bool eliminate) 
    { eliminate_internal = eliminate; }

    void SetPrint (bool ap);
    void SetPrintElmat (bool ap);
    void SetElmatEigenValues (bool ee);

    ///
    void GalerkinProjection ();


    virtual void ComputeInternal (BaseVector & u, LocalHeap & lh) const = 0;
  
    ///
    void SetEpsRegularization(double val)
    { eps_regularization = val; }
    ///
    void SetUnusedDiag (double val)
    { unuseddiag = val; }

    ///
    int UseGalerkin () const
    { return galerkin; }

    ///
    virtual string GetClassName () const
    {
      return "BilinearForm";
    }

    ///
    virtual void PrintReport (ostream & ost);

    ///
    virtual void MemoryUsage (Array<MemoryUsageStruct*> & mu) const;

    ///
    virtual BaseVector * CreateVector() const = 0;

    /// frees matrix 
    virtual void CleanUpLevel() { ; }
  private:
    ///
    virtual void DoAssemble (LocalHeap & lh) = 0;

    ///
    virtual void AllocateMatrix () = 0;


#ifdef PARALLEL

    virtual void AllocateConsistentMatrix () = 0;

#endif
  };








  template <class SCAL>
  class NGS_DLL_HEADER S_BilinearForm : public BilinearForm
  {
  protected:
  public:
    S_BilinearForm (const FESpace & afespace, const string & aname,
		    const Flags & flags)
      : BilinearForm (afespace, aname, flags) { ; }

    ///
    S_BilinearForm (const FESpace & afespace, 
		    const FESpace & afespace2,
		    const string & aname, const Flags & flags)
      : BilinearForm (afespace, afespace2, aname, flags) { ; }


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

    virtual void ComputeInternal (BaseVector & u, LocalHeap & lh) const;

    ///
    virtual void DoAssemble (LocalHeap & lh);
    ///
    virtual void DoAssembleIndependent (BitArray & useddof, LocalHeap & lh);
    ///
    virtual void AssembleLinearization (const BaseVector & lin,
					LocalHeap & lh, 
					bool reallocate = 0);
    ///
    virtual void AddElementMatrix (const Array<int> & dnums1,
				   const Array<int> & dnums2,
				   const FlatMatrix<SCAL> & elmat,
				   bool inner_element, int elnr,
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
    virtual BaseVector * CreateVector() const;

    ///
    virtual void CleanUpLevel();

    ///
    virtual void AddElementMatrix (const Array<int> & dnums1,
				   const Array<int> & dnums2,
				   const FlatMatrix<TSCAL> & elmat,
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

    virtual void LapackEigenSystem(FlatMatrix<TSCAL> & elmat, LocalHeap & lh) const;

#ifdef PARALLEL
    virtual void AllocateConsistentMatrix ();
#endif
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

    virtual BaseVector * CreateVector() const;

    virtual void AddElementMatrix (const Array<int> & dnums1,
				   const Array<int> & dnums2,
				   const FlatMatrix<TSCAL> & elmat,
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

    virtual void LapackEigenSystem(FlatMatrix<TSCAL> & elmat, LocalHeap & lh) const;
#ifdef PARALLEL
    //   ///
    virtual void AllocateConsistentMatrix ();

#endif
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
    virtual BaseVector * CreateVector() const;

    virtual void AddElementMatrix (const Array<int> & dnums1,
				   const Array<int> & dnums2,
				   const FlatMatrix<TSCAL> & elmat,
				   bool inner_element, int elnr,
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
#ifdef PARALLEL
    //   ///
    virtual void AllocateConsistentMatrix ();

#endif
  };





  extern NGS_DLL_HEADER BilinearForm * CreateBilinearForm (const FESpace * space,
					    const string & name,
					    const Flags & flags);



  ///
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
    // virtual void MultTransAdd (double val, const BaseVector & v, BaseVector & prod) const;
    ///
    virtual BaseVector * CreateVector () const;

    virtual int VHeight() const
    {
      return bf->GetFESpace().GetNDof(); 
    }
    virtual int VWidth() const
    {
      return bf->GetFESpace().GetNDof(); 
    }
  };


  class NGS_DLL_HEADER LinearizedBilinearFormApplication : public BilinearFormApplication 
  {
  protected:
    const BaseVector * veclin;

  public:
    LinearizedBilinearFormApplication (const BilinearForm * abf,
				       const BaseVector * aveclin);

    ///
    virtual void Mult (const BaseVector & v, BaseVector & prod) const;
    ///
    virtual void MultAdd (double val, const BaseVector & v, BaseVector & prod) const;
    ///
    virtual void MultAdd (Complex val, const BaseVector & v, BaseVector & prod) const;
    ///
    // virtual void MultTransAdd (double val, const BaseVector & v, BaseVector & prod) const;

  };

}


#endif
