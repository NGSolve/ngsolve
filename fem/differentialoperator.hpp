#ifndef FILE_DIFFERENTIALOPERATOR
#define FILE_DIFFERENTIALOPERATOR

/*********************************************************************/
/* File:   differentialoperator.hpp                                  */
/* Author: Joachim Schoeberl                                         */
/* Date:   24. Nov. 2009                                             */
/*********************************************************************/

#include "finiteelement.hpp"
#include "intrule.hpp"


namespace ngfem
{

  /*
    How one factor of the value transforms from the reference element (NGSolve technical note 2026-08).
    k-forms counted from the bottom (0, 1) or from the top (n-1, n).
  */
  enum FORM_TYPE : uint8_t
  {
    FORM_VERTEX = 0,    // 0-form:      1
    FORM_EDGE = 1,      // 1-form:      F^{-T}
    FORM_FACET = 2,     // (n-1)-form:  F/J
    FORM_ELEMENT = 3,   // n-form:      1/J
    FORM_COMP = 4       // Cartesian component index, not transformed
  };

  inline bool FromTop (FORM_TYPE f) { return f == FORM_FACET || f == FORM_ELEMENT; }
  // distance of the degree from 0 (bottom) or n (top)
  inline int FormOffset (FORM_TYPE f) { return (f == FORM_EDGE || f == FORM_FACET) ? 1 : 0; }
  inline int FormDegree (FORM_TYPE f, int dim) { return FromTop(f) ? dim-FormOffset(f) : FormOffset(f); }
  inline bool HasIndex (FORM_TYPE f) { return f == FORM_COMP || FormOffset(f) == 1; }

  struct FormSignature
  {
    Array<FORM_TYPE> factors;
    bool exact = true;      // false: holds only up to terms with grad F on curved elements
    int detpower = 0;       // value additionally multiplied by J^detpower
  };

  inline optional<FormSignature> CombineForms (const optional<FormSignature> & base,
                                               const Array<FORM_TYPE> & before,
                                               const Array<FORM_TYPE> & after = Array<FORM_TYPE>())
  {
    if (!base) return nullopt;
    return FormSignature { Array<FORM_TYPE>(before + base->factors + after), base->exact, base->detpower };
  }

  /**
     Differential Operator.
     Base-class for run-time polymorphismus.
     Provides application and transpose-application
  */
  class DifferentialOperator
  {
  private:
    int dim;
    int blockdim;
    Array<int> dimensions;
    
    int vsdim;    // symmetric 3x3 matrix has dim=9, but vector-space dim=6
    optional<Matrix<>> vsembedding;
    Matrix<> vsembedding_pinv;   // left inverse: vs = pinv * full
  protected:
    VorB vb;
    int difforder;
  public:
    /*
    [[deprecated("Use DifferentialOperator(int,int,VorB,int) instead")]]
    NGS_DLL_HEADER DifferentialOperator(int adim, int ablockdim, bool boundary, int adifforder)
      : dim(adim), blockdim(ablockdim), vb(boundary ? BND : VOL),  difforder(adifforder)
     {
       if (blockdim == 1)
         dimensions = Array<int> ( { dim } );
       else
         dimensions = Array<int> ( { dim/blockdim, blockdim });
     }
    */
    NGS_DLL_HEADER DifferentialOperator(int adim, int ablockdim, VorB avb, int adifforder)
      : dim(adim), blockdim(ablockdim), vb(avb), difforder(adifforder)
    { 
      if (blockdim == 1)
        dimensions = Array<int> ( { dim } );
      else if (dim == 1)
        dimensions = Array<int> ( { blockdim } );
      else
        dimensions = Array<int> ( { dim/blockdim, blockdim });

      vsdim = dim;
    }
    /// destructor
    NGS_DLL_HEADER virtual ~DifferentialOperator () = default;

    virtual void DoArchive(Archive & ar) { ; }

    void SetDimensions (const Array<int> & adims) { dimensions = adims; }
    NGS_DLL_HEADER void SetVectorSpaceEmbedding (Matrix <> emb);
    optional<FlatMatrix<>> GetVSEmbedding() const { return vsembedding; }
    
    ///
    NGS_DLL_HEADER virtual string Name() const; // { return typeid(*this).name(); }
    /// dimension of range
    int Dim() const { return dim; }
    int VSDim() const { return vsdim; }
    const Array<int> & Dimensions() const { return dimensions; }
    /// number of copies of finite element by BlockDifferentialOperator
    int BlockDim() const { return blockdim; }
    /// does it live on the boundary ?
    bool Boundary() const { return vb == BND; }
    VorB VB() const { return vb; }

    virtual bool SupportsVB (VorB checkvb) const { return checkvb == vb; }
    virtual shared_ptr<DifferentialOperator> GetTrace() const
    {
      return nullptr;
      // throw Exception("GetTrace not overloaded for DifferentialOperator"+string(typeid(*this).name()));
    }
    /// form types of the factors of the value, nullopt if not declared
    virtual optional<FormSignature> Forms() const { return nullopt; }
    /// total polynomial degree is reduced by this order (i.e. minimal difforder)
    int DiffOrder() const { return difforder; } 

    /// polynomial order of the result, used to choose integration rules
    virtual int IntOrder (const FiniteElement & fel) const
    {
      auto et = fel.ElementType();
      return (et == ET_TRIG || et == ET_TET) ? fel.Order()-difforder : fel.Order();
    }

    virtual IntRange UsedDofs(const FiniteElement & fel) const { return IntRange(0, fel.GetNDof()); }
    // dofs used on facet facetnr of the element
    virtual IntRange UsedDofs(const FiniteElement & fel, int facetnr) const { return UsedDofs(fel); }

    virtual bool operator== (const DifferentialOperator & diffop2) const { return false; }

    /// calculates the matrix
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		BareSliceMatrix<double,ColMajor> mat,   
		LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & bfel,
		const BaseMappedIntegrationPoint & bmip,
		BareSliceMatrix<Complex,ColMajor> mat, 
		LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		BareSliceMatrix<double,ColMajor> mat,   
		LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		BareSliceMatrix<Complex,ColMajor> mat,   
		LocalHeap & lh) const;
    
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const SIMD_BaseMappedIntegrationRule & mir,
		BareSliceMatrix<SIMD<double>> mat) const;

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const SIMD_BaseMappedIntegrationRule & mir,
		BareSliceMatrix<SIMD<Complex>> mat) const;

    /// Bmat = vs-embedding * BmatVS    (if vs-embedding is set)
    NGS_DLL_HEADER virtual void
    CalcMatrixVS (const FiniteElement & fel,
                  const BaseMappedIntegrationPoint & mip,
                  SliceMatrix<double,ColMajor> mat,   
                  LocalHeap & lh) const;

    /// VSDim rows per dof (if vs-embedding is set), default extracts from CalcMatrix
    NGS_DLL_HEADER virtual void
    CalcMatrixVS (const FiniteElement & fel,
                  const SIMD_BaseMappedIntegrationRule & mir,
                  BareSliceMatrix<SIMD<double>> mat) const;

    NGS_DLL_HEADER virtual void
    CalcLinearizedMatrix (const FiniteElement & fel,
                          const BaseMappedIntegrationRule & mir,
                          BareSliceVector<double> x,
                          SliceMatrix<double,ColMajor> mat,   
                          LocalHeap & lh) const;

    NGS_DLL_HEADER virtual bool IsNonlinear() const { return false; }

    // second derivative of \sum_ipt wprime * B(u) 
    NGS_DLL_HEADER virtual void
    CalcHessianAdd (const FiniteElement & fel,
                    const BaseMappedIntegrationRule & mir,
                    SliceMatrix<> wprime,
                    BareSliceVector<> elvecu,
                    SliceMatrix<> hessian,   
                    LocalHeap & lh) const { ; } 

    // simd versions, nonlinear diffops throw NOSIMD if not overloaded
    NGS_DLL_HEADER virtual void
    CalcLinearizedMatrix (const FiniteElement & fel,
                          const SIMD_BaseMappedIntegrationRule & mir,
                          BareSliceVector<double> x,
                          BareSliceMatrix<SIMD<double>> mat) const;

    NGS_DLL_HEADER virtual void
    CalcHessianAdd (const FiniteElement & fel,
                    const SIMD_BaseMappedIntegrationRule & mir,
                    BareSliceMatrix<SIMD<double>> wprime,
                    BareSliceVector<> elvecu,
                    SliceMatrix<> hessian,   
                    LocalHeap & lh) const;
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<float> x, 
	   FlatVector<float> flux,
	   LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<Complex> x, 
	   FlatVector<Complex> flux,
	   LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationRule & mir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<double> flux,
	   LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationRule & mir,
	   BareSliceVector<Complex> x, 
	   BareSliceMatrix<Complex> flux,
	   LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<SIMD<double>> flux) const;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<Complex> x, 
	   BareSliceMatrix<SIMD<Complex>> flux) const;

    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		FlatVector<double> flux,
		BareSliceVector<double> x, 
		LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		FlatVector<Complex> flux,
		BareSliceVector<Complex> x, 
		LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<double> flux,
		BareSliceVector<double> x, 
		LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<Complex> flux,
		BareSliceVector<Complex> x, 
		LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const;

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<Complex>> flux,
              BareSliceVector<Complex> x) const;


    NGS_DLL_HEADER virtual void
    ApplyLinearizedTrans (const FiniteElement & fel,
                          const BaseMappedIntegrationRule & mir,
                          SliceVector<double> elveclin,
                          FlatMatrix<double> flux,
                          BareSliceVector<double> x, 
                          LocalHeap & lh) const
    {
      ApplyTrans (fel, mir, flux, x, lh);
    }

    // simd version, x += B'(elveclin)^T flux; default via CalcLinearizedMatrix for nonlinear diffops
    NGS_DLL_HEADER virtual void
    AddLinearizedTrans (const FiniteElement & fel,
                          const SIMD_BaseMappedIntegrationRule & mir,
                          BareSliceVector<double> elveclin,
                          BareSliceMatrix<SIMD<double>> flux,
                          BareSliceVector<double> x,
                          LocalHeap & lh) const;


    /// calculates matrix on reference element
    
    // dimension on refelement (e.g. 2 for surface gradient)
    NGS_DLL_HEADER virtual int DimRef() const;
    
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
                const IntegrationPoint & ip,
                BareSliceMatrix<double,ColMajor> mat,
                LocalHeap & lh) const;
    
    NGS_DLL_HEADER virtual void
    CalcTransformationMatrix (const BaseMappedIntegrationPoint & mip,
                              SliceMatrix<double> trans,
                              LocalHeap & lh) const;

    NGS_DLL_HEADER
    virtual string GenerateTransformationCode (string invar, string outvar, bool trans,
                                               bool curved = false) const;
    
    NGS_DLL_HEADER virtual shared_ptr<CoefficientFunction> DiffShape (shared_ptr<CoefficientFunction> proxy,
                                                       shared_ptr<CoefficientFunction> dir,
                                                       bool Eulerian = false) const;

    NGS_DLL_HEADER virtual list<tuple<string,double>> Timing (const FiniteElement & fel, const BaseMappedIntegrationRule & mir) const;

    NGS_DLL_HEADER virtual void CheckElement (const FiniteElement& fel) const { } 
    // true if CheckElement tests the element family
    virtual bool DeclaresElementType () const { return false; }
    // acts on volume elements only at facet points (facet spaces)
    virtual bool OnlyOnFacets () const { return false; }
    // dimensions of the physical space and of the reference element, -1: unknown
    virtual int DimSpace () const { return -1; }
    virtual int DimElement () const { return -1; }
  };

  // one differential operator per co-dimension of the element.
  // level = number of traces taken, operators of lower co-dimension are dropped
  class EvaluatorTable
  {
    shared_ptr<DifferentialOperator> ops[4];
    VorB level = VOL;
  public:
    EvaluatorTable () = default;
    EvaluatorTable (shared_ptr<DifferentialOperator> vol, shared_ptr<DifferentialOperator> bnd,
                    shared_ptr<DifferentialOperator> bbnd, shared_ptr<DifferentialOperator> bbbnd,
                    VorB alevel = VOL)
      : ops{vol, bnd, bbnd, bbbnd} { SetLevel(alevel); }
    // a single operator, at its own co-dimension
    explicit EvaluatorTable (shared_ptr<DifferentialOperator> op)
    {
      if (op)
        {
          level = op->VB();
          ops[level] = op;
        }
    }

    shared_ptr<DifferentialOperator> & operator[] (int vb) { return ops[vb]; }
    // the operator for elements of co-dimension vb
    const shared_ptr<DifferentialOperator> & operator[] (int vb) const { return ops[vb]; }

    VorB Level() const { return level; }
    void SetLevel (VorB alevel)
    {
      level = alevel;
      for (int i = 0; i < level; i++)
        ops[i] = nullptr;
    }
    // the operator of the level
    const shared_ptr<DifferentialOperator> & Active() const { return ops[level]; }
    // the operator of the level, or the next one below it
    shared_ptr<DifferentialOperator> First() const
    {
      for (int i = level; i < 4; i++)
        if (ops[i]) return ops[i];
      return nullptr;
    }
    // level of the first operator
    VorB FirstLevel() const
    {
      for (int i = level; i < 4; i++)
        if (ops[i]) return VorB(i);
      return level;
    }
    bool HasTrace() const { return level < BBBND && ops[level+1]; }
    // one level down, empty if there is none
    EvaluatorTable Trace() const
    {
      if (level == BBBND) return EvaluatorTable();
      EvaluatorTable trace(*this);
      trace.SetLevel (VorB(level+1));
      return trace;
    }
  };

  

}



#endif
