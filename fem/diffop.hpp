#ifndef FILE_DIFFOP
#define FILE_DIFFOP

/*********************************************************************/
/* File:   diffop.hpp                                                */
/* Author: Joachim Schoeberl                                         */
/* Date:   24. Nov. 2009                                             */
/*********************************************************************/

namespace ngfem
{



  /**
     Differential Operator.
     Base-class for template-polymorphismus.
     Provides application and transpose-application.
     Operations can be applied for one integration point, or for the whole integration rule at once.
  */
  template<class DOP>
  class DiffOp
  {
  public:

    static string Name() { return typeid(DiffOp<DOP>()).name(); }
    static constexpr bool SUPPORT_PML = false;
    // static Array<int> GetDimensions() { return Array<int> ( { DOP::DIM_DMAT } ); };
    static INT<1> GetDimensions() { return { DOP::DIM_DMAT }; };
    static bool SupportsVB (VorB checkvb) { return DOP::DIM_SPACE-DOP::DIM_ELEMENT == int(checkvb); }


    typedef void DIFFOP_TRACE;  // 
    
    /**
       Computes the B-matrix. 
       The height is DIM_DMAT, the width is fel.GetNDof(). 
       FEL is the FiniteElement type specified in the BDB-Integrator 
       mip is the mapped integration point containing the Jacobi-Matrix 
       MAT is the resulting matrix (usually a FixedHeightMatrix)
    */

    template <typename FEL, typename MIP, typename MAT>
    static void GenerateMatrix (const FEL & fel, const MIP & mip,
				MAT & mat, LocalHeap & lh)
    { 
      cout << "bad function calling" << endl;
    }

    /// tbd
    template <typename FEL, typename MIR, typename MAT>
    static void GenerateMatrixIR (const FEL & fel, const MIR & mir,
                                  MAT & mat, LocalHeap & lh)
    {
      for (size_t i = 0; i < mir.Size(); i++)
        DOP::GenerateMatrix (fel, mir[i], mat.Rows(i*DOP::DIM_DMAT, (i+1)*DOP::DIM_DMAT), lh);
    }

    template <typename FEL, typename MIR>
    static void GenerateMatrixSIMDIR (const FEL & fel, const MIR & mir, BareSliceMatrix<SIMD<double>> mat)
    {
      throw ExceptionNOSIMD (string("generate matrix simdir not implemented for diffop ") + typeid(DOP).name());
    }
    /**
       Applies the B-matrix.
       Computes matrix-vector product with the B-matrix
    */
    template <typename FEL, typename MIP, class TVX, class TVY>
    static void Apply (const FEL & fel, const MIP & mip,
		       const TVX & x, TVY & y,
		       LocalHeap & lh)
    {
      // typedef typename TVY::TSCAL TSCAL;
      typedef typename MIP::TSCAL TSCAL;
      HeapReset hr(lh);

      FlatMatrixFixHeight<DOP::DIM_DMAT, TSCAL> mat(DOP::DIM*fel.GetNDof(), lh);
      DOP::GenerateMatrix (fel, mip, mat, lh);
      y = mat * x;
    }

    /// Computes B-matrix times element vector in many points
    template <typename FEL, class MIR, class TVX, class TVY>
    static void ApplyIR (const FEL & fel, const MIR & mir,
			 const TVX & x, TVY & y,
			 LocalHeap & lh)
    {
      for (size_t i = 0; i < mir.Size(); i++)
        {
          HeapReset hr(lh);
          DOP::Apply (fel, mir[i], x, y.Row(i), lh);
        }
    }

    /// Computes B-matrix times element vector in many points
    template <typename FEL, class MIR, class TVX, class TVY>
    static void ApplySIMDIR (const FEL & fel, const MIR & mir,
                             const TVX & x, TVY & y)
    // LocalHeap & lh)
    {
      throw ExceptionNOSIMD (string("apply simdir not implemented for diffop ") + typeid(DOP).name());
    }


    /// Computes Transpose (B-matrix) times point value
    template <typename FEL, typename MIP, class TVX, class TVY>
    static void ApplyTrans (const FEL & fel, const MIP & mip,
			    const TVX & x, TVY & y,
			    LocalHeap & lh) 
    {
      // typedef typename TVY::TSCAL TSCAL;
      typedef typename MIP::TSCAL TSCAL;

      HeapReset hr(lh);

      FlatMatrixFixHeight<DOP::DIM_DMAT, TSCAL> mat(DOP::DIM*fel.GetNDof(), lh);
      DOP::GenerateMatrix (fel, mip, mat, lh);
      y.Range(0,DOP::DIM*fel.GetNDof()) = Trans (mat) * x;
    }

    /// Computes Transpose (B-matrix) times point value
    template <typename FEL, typename MIP, class TVX, class TVY>
    static void ApplyTransAdd (const FEL & fel, const MIP & mip,
                               const TVX & x, TVY & y,
                               LocalHeap & lh) 
    {
      typedef typename MIP::TSCAL TSCAL;

      HeapReset hr(lh);

      FlatMatrixFixHeight<DOP::DIM_DMAT, TSCAL> mat(DOP::DIM*fel.GetNDof(), lh);
      DOP::GenerateMatrix (fel, mip, mat, lh);
      y.Range(DOP::DIM*fel.GetNDof()) += Trans (mat) * x;
    }


    /// Computes Transpose (B-matrix) times point value
    template <typename FEL, typename MIR, class TVX, class TVY>
    static void ApplyTransIR (const FEL & fel, const MIR & mir,
			      const TVX & x, TVY & y,
			      LocalHeap & lh) 
    {
      y.Range(0,DOP::DIM*fel.GetNDof()) = 0.0;
      for (size_t i = 0; i < mir.Size(); i++)
        {
          HeapReset hr(lh);        
          ApplyTransAdd (fel, mir[i], x.Row(i), y, lh);
        }
    }

    /// Computes Transpose (B-matrix) times point value    
    template <typename FEL, class MIR, class TVX, class TVY>
    static void AddTransSIMDIR (const FEL & fel, const MIR & mir,
                                const TVX & x, TVY & y)
    // LocalHeap & lh)
    {
      throw ExceptionNOSIMD (string("AddTrans simdir not implemented for diffop ") + typeid(DOP).name());
    }


    static shared_ptr<CoefficientFunction>
    DiffShape (shared_ptr<CoefficientFunction> proxy,
               shared_ptr<CoefficientFunction> dir,
               bool Eulerian) 
    {
      throw Exception (string("shape derivative not implemented for DifferentialOperator")+Name());
    }

  };





  

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
    NGS_DLL_HEADER virtual ~DifferentialOperator ();
    void SetDimensions (const Array<int> & adims) { dimensions = adims; }
    void SetVectorSpaceEmbedding (Matrix <> emb)
    { vsembedding = emb; vsdim = emb.Width(); }
    optional<FlatMatrix<>> GetVSEmbedding() const { return vsembedding; }
    
    ///
    virtual string Name() const { return typeid(*this).name(); }
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
    /// total polynomial degree is reduced by this order (i.e. minimal difforder)
    int DiffOrder() const { return difforder; } 

    virtual IntRange UsedDofs(const FiniteElement & fel) const { return IntRange(0, fel.GetNDof()); }

    virtual bool operator== (const DifferentialOperator & diffop2) const { return false; }

    /// calculates the matrix
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		SliceMatrix<double,ColMajor> mat,   
		LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & bfel,
		const BaseMappedIntegrationPoint & bmip,
		SliceMatrix<Complex,ColMajor> mat, 
		LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		SliceMatrix<double,ColMajor> mat,   
		LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		SliceMatrix<Complex,ColMajor> mat,   
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
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<double> x, 
	   FlatVector<double> flux,
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

    
    virtual shared_ptr<CoefficientFunction> DiffShape (shared_ptr<CoefficientFunction> proxy,
                                                       shared_ptr<CoefficientFunction> dir,
                                                       bool Eulerian = false) const
    {
      throw Exception (string("shape derivative not implemented for DifferentialOperator")+typeid(*this).name());
    }

    NGS_DLL_HEADER virtual list<tuple<string,double>> Timing (const FiniteElement & fel, const BaseMappedIntegrationRule & mir) const;
  };



  class BlockDifferentialOperator : public DifferentialOperator
  {
  protected:
    shared_ptr<DifferentialOperator> diffop;
    int dim;
    int comp;
  public:
    BlockDifferentialOperator (shared_ptr<DifferentialOperator> adiffop, 
			       int adim, int acomp = -1)
      : DifferentialOperator(adim*adiffop->Dim(), adim*adiffop->BlockDim(),
                             adiffop->VB(), adiffop->DiffOrder()),
        diffop(adiffop), dim(adim), comp(acomp)
    {
      if(adiffop->Dimensions().Size()==0)
	SetDimensions( { BlockDim() } );
    }

    NGS_DLL_HEADER virtual ~BlockDifferentialOperator ();
    
    virtual string Name() const override { return diffop->Name(); }
    shared_ptr<DifferentialOperator> BaseDiffOp() const { return diffop; }
    virtual bool SupportsVB (VorB checkvb) const override { return diffop->SupportsVB(checkvb); }
    
    virtual IntRange UsedDofs(const FiniteElement & fel) const override { return dim*diffop->UsedDofs(fel); }

    shared_ptr<DifferentialOperator> GetTrace() const override
    {
      if (auto diffoptrace = diffop->GetTrace())
        return make_shared<BlockDifferentialOperator> (diffoptrace,
                                                       dim, comp);
      else
        return nullptr;
    }
    
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		SliceMatrix<double,ColMajor> mat, 
		LocalHeap & lh) const override;    

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const SIMD_BaseMappedIntegrationRule & mir,
		BareSliceMatrix<SIMD<double>> mat) const override;
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<SIMD<double>> flux) const override;
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<double> flux,
                BareSliceVector<double> x, 
                LocalHeap & lh) const override;
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<Complex> flux,
                BareSliceVector<Complex> x, 
                LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<double> flux,
		BareSliceVector<double> x, 
		LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<Complex> flux,
		BareSliceVector<Complex> x, 
		LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const override;

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<Complex>> flux,
              BareSliceVector<Complex> x) const override;


    shared_ptr<CoefficientFunction> DiffShape (shared_ptr<CoefficientFunction> proxy,
                                               shared_ptr<CoefficientFunction> dir,
                                               bool Eulerian) const override;
  };




  class BlockDifferentialOperatorTrans : public DifferentialOperator
  {
  protected:
    shared_ptr<DifferentialOperator> diffop;
    int dim;
    int comp;
  public:
    BlockDifferentialOperatorTrans (shared_ptr<DifferentialOperator> adiffop, 
                                    int adim, int acomp = -1)
      : DifferentialOperator(adim*adiffop->Dim(), adim*adiffop->BlockDim(),
                             adiffop->VB(), adiffop->DiffOrder()),
        diffop(adiffop), dim(adim), comp(acomp)
    {
      // dimensions = Array<int> ( { adim, adiffop->Dim() });
      SetDimensions ( { adim, adiffop->Dim() } );
    }

    NGS_DLL_HEADER virtual ~BlockDifferentialOperatorTrans ();
    
    virtual string Name() const override { return diffop->Name(); }
    shared_ptr<DifferentialOperator> BaseDiffOp() const { return diffop; }
    virtual bool SupportsVB (VorB checkvb) const override { return diffop->SupportsVB(checkvb); }
    
    virtual IntRange UsedDofs(const FiniteElement & fel) const override { return dim*diffop->UsedDofs(fel); }

    shared_ptr<DifferentialOperator> GetTrace() const override
    {
      if (auto diffoptrace = diffop->GetTrace())      
        return make_shared<BlockDifferentialOperatorTrans> (diffoptrace,
                                                            dim, comp);
      else
        return nullptr;
    }
    

    
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		SliceMatrix<double,ColMajor> mat, 
		LocalHeap & lh) const override;    

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const SIMD_BaseMappedIntegrationRule & mir,
		BareSliceMatrix<SIMD<double>> mat) const override;
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<SIMD<double>> flux) const override;
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<double> flux,
                BareSliceVector<double> x, 
                LocalHeap & lh) const override;
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<Complex> flux,
                BareSliceVector<Complex> x, 
                LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const override;

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<Complex>> flux,
              BareSliceVector<Complex> x) const override;


    shared_ptr<CoefficientFunction> DiffShape (shared_ptr<CoefficientFunction> proxy,
                                               shared_ptr<CoefficientFunction> dir,
                                               bool Eulerian) const override;
  };

  
  // like BlockDifferentialOperator, but element is CompoundFE here
  class VectorDifferentialOperator : public DifferentialOperator
  {
  protected:
    shared_ptr<DifferentialOperator> diffop;
    int dim;
  public:
    VectorDifferentialOperator (shared_ptr<DifferentialOperator> adiffop, 
                                     int adim)
      : DifferentialOperator(adim*adiffop->Dim(), /* adim* */ adiffop->BlockDim(),
                             adiffop->VB(), adiffop->DiffOrder()),
        diffop(adiffop), dim(adim)
    {
      if (adiffop->Dimensions().Size() == 0)
        // dimensions = Array<int> ( { adim });
        SetDimensions ( { adim } );
      else
        // dimensions = Array<int> ( { adim, adiffop->Dim() });
        SetDimensions ( { adim, adiffop->Dim() } );
    }

    NGS_DLL_HEADER virtual ~VectorDifferentialOperator ();
    
    virtual string Name() const override { return diffop->Name(); }
    shared_ptr<DifferentialOperator> BaseDiffOp() const { return diffop; }
    virtual bool SupportsVB (VorB checkvb) const override { return diffop->SupportsVB(checkvb); }
    
    virtual IntRange UsedDofs(const FiniteElement & fel) const override { return IntRange(0, fel.GetNDof()); }

    shared_ptr<DifferentialOperator> GetTrace() const override
    {
      if (auto diffoptrace = diffop->GetTrace())      
        return make_shared<VectorDifferentialOperator> (diffoptrace, dim);
      else
        return nullptr;
    }
    
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		SliceMatrix<double,ColMajor> mat, 
		LocalHeap & lh) const override;    

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const SIMD_BaseMappedIntegrationRule & mir,
		BareSliceMatrix<SIMD<double>> mat) const override;
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   BareSliceVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<SIMD<double>> flux) const override;
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<double> flux,
                BareSliceVector<double> x, 
                LocalHeap & lh) const override;
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<Complex> flux,
                BareSliceVector<Complex> x, 
                LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const override;

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<Complex>> flux,
              BareSliceVector<Complex> x) const override;

    shared_ptr<CoefficientFunction> DiffShape (shared_ptr<CoefficientFunction> proxy,
                                               shared_ptr<CoefficientFunction> dir,
                                               bool Eulerian) const override;
  };


  class MatrixDifferentialOperator : public DifferentialOperator
  {
  protected:
    shared_ptr<DifferentialOperator> diffop;
    int vdim;
  public:
    MatrixDifferentialOperator (shared_ptr<DifferentialOperator> adiffop, 
                                int avdim)
      : DifferentialOperator(sqr(avdim)*adiffop->Dim(), adiffop->BlockDim(),
                             adiffop->VB(), adiffop->DiffOrder()),
        diffop(adiffop), vdim(avdim)
    {
      if (adiffop->Dimensions().Size() == 0)
        // dimensions = Array<int> ( { avdim, avdim });
        SetDimensions ( { avdim, avdim } );
      else
        throw Exception("no matrix-valued of vector-valued possible");
    }

    NGS_DLL_HEADER virtual ~MatrixDifferentialOperator ();
    
    virtual string Name() const override { return diffop->Name(); }
    shared_ptr<DifferentialOperator> BaseDiffOp() const { return diffop; }
    virtual bool SupportsVB (VorB checkvb) const override { return diffop->SupportsVB(checkvb); }
    
    virtual IntRange UsedDofs(const FiniteElement & fel) const override { return IntRange(0, fel.GetNDof()); }

    shared_ptr<DifferentialOperator> GetTrace() const override
    {
      if (auto diffoptrace = diffop->GetTrace())      
        return make_shared<MatrixDifferentialOperator> (diffoptrace, vdim);
      else
        return nullptr;
    }
    
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		SliceMatrix<double,ColMajor> mat, 
		LocalHeap & lh) const override;    
  };

  class SymMatrixDifferentialOperator : public DifferentialOperator
  {
  protected:
    shared_ptr<DifferentialOperator> diffop;
    int vdim;
  public:
    NGS_DLL_HEADER SymMatrixDifferentialOperator (shared_ptr<DifferentialOperator> adiffop, 
                                                  int avdim);

    NGS_DLL_HEADER virtual ~SymMatrixDifferentialOperator ();
    
    virtual string Name() const override { return diffop->Name(); }
    shared_ptr<DifferentialOperator> BaseDiffOp() const { return diffop; }
    virtual bool SupportsVB (VorB checkvb) const override { return diffop->SupportsVB(checkvb); }
    
    virtual IntRange UsedDofs(const FiniteElement & fel) const override
    { return IntRange(0, fel.GetNDof()); }

    shared_ptr<DifferentialOperator> GetTrace() const override
    {
      if (auto diffoptrace = diffop->GetTrace())      
        return make_shared<SymMatrixDifferentialOperator> (diffoptrace, vdim);
      else
        return nullptr;
    }
    
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		SliceMatrix<double,ColMajor> mat, 
		LocalHeap & lh) const override;    

    NGS_DLL_HEADER virtual void
    CalcMatrixVS (const FiniteElement & fel,
                  const BaseMappedIntegrationPoint & mip,
                  SliceMatrix<double,ColMajor> mat, 
                  LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void 
    CalcMatrix (const FiniteElement & bfel,
                const SIMD_BaseMappedIntegrationRule & mir,
                BareSliceMatrix<SIMD<double>> bmat) const override;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<SIMD<double>> flux) const override;

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const override;
  };



  
  class SymDevMatrixDifferentialOperator : public DifferentialOperator
  {
  protected:
    shared_ptr<DifferentialOperator> diffop;
    int vdim;
  public:
    NGS_DLL_HEADER SymDevMatrixDifferentialOperator (shared_ptr<DifferentialOperator> adiffop, 
                                                     int avdim);

    NGS_DLL_HEADER virtual ~SymDevMatrixDifferentialOperator ();
    
    virtual string Name() const override { return diffop->Name(); }
    shared_ptr<DifferentialOperator> BaseDiffOp() const { return diffop; }
    virtual bool SupportsVB (VorB checkvb) const override { return diffop->SupportsVB(checkvb); }
    
    virtual IntRange UsedDofs(const FiniteElement & fel) const override
    { return IntRange(0, fel.GetNDof()); }

    shared_ptr<DifferentialOperator> GetTrace() const override
    {
      if (auto diffoptrace = diffop->GetTrace())      
        return make_shared<SymMatrixDifferentialOperator> (diffoptrace, vdim);
      else
        return nullptr;
    }
    
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		SliceMatrix<double,ColMajor> mat, 
		LocalHeap & lh) const override;    

    NGS_DLL_HEADER virtual void
    CalcMatrixVS (const FiniteElement & fel,
                  const BaseMappedIntegrationPoint & mip,
                  SliceMatrix<double,ColMajor> mat, 
                  LocalHeap & lh) const override;
    /*
    NGS_DLL_HEADER virtual void 
    CalcMatrix (const FiniteElement & bfel,
                const SIMD_BaseMappedIntegrationRule & mir,
                BareSliceMatrix<SIMD<double>> bmat) const override;

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<SIMD<double>> flux) const override;

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const override;
    */
  };

  
  

  /**
     Connect compile-time polymorph DiffOp to run-time polymorph DifferentialOperator.
   */
  template <typename DIFFOP>
  class T_DifferentialOperator : public DifferentialOperator
  {
  protected:
    enum { DIM_SPACE   = DIFFOP::DIM_SPACE };
    enum { DIM_ELEMENT = DIFFOP::DIM_ELEMENT };
    enum { DIM_DMAT    = DIFFOP::DIM_DMAT };
    enum { DIM         = DIFFOP::DIM };

  public:
    T_DifferentialOperator()
      : DifferentialOperator(DIFFOP::DIM_DMAT, 1, VorB(int(DIM_SPACE)-int(DIM_ELEMENT)), DIFFOP::DIFFORDER)
    {
      Array<int> hdims;
      hdims = DIFFOP::GetDimensions();
      SetDimensions ( hdims );
    }

    virtual string Name() const override { return DIFFOP::Name(); }
    
    virtual bool operator== (const DifferentialOperator & diffop2) const override
    { return typeid(*this) == typeid(diffop2); }

    virtual bool SupportsVB (VorB checkvb) const override { return DIFFOP::SupportsVB(checkvb); }

    virtual shared_ptr<DifferentialOperator> GetTrace() const override
    {
      if constexpr (is_same_v<void,typename DIFFOP::DIFFOP_TRACE>)
                     return nullptr;
      else
        return make_shared<T_DifferentialOperator<typename DIFFOP::DIFFOP_TRACE>> ();
    }
    
    virtual void
    CalcMatrix (const FiniteElement & bfel,
		const BaseMappedIntegrationPoint & bmip,
		SliceMatrix<double,ColMajor> mat, 
		LocalHeap & lh) const override;

    virtual void
    CalcMatrix (const FiniteElement & bfel,
		const BaseMappedIntegrationPoint & bmip,
		SliceMatrix<Complex,ColMajor> mat, 
		LocalHeap & lh) const override;

    virtual void
    CalcMatrix (const FiniteElement & bfel,
		const BaseMappedIntegrationRule & bmir,
		SliceMatrix<double,ColMajor> mat, 
		LocalHeap & lh) const override;

    virtual void
    CalcMatrix (const FiniteElement & fel,
		const SIMD_BaseMappedIntegrationRule & mir,
		BareSliceMatrix<SIMD<double>> mat) const override;
    
#ifndef FASTCOMPILE
    virtual void
    Apply (const FiniteElement & bfel,
	   const BaseMappedIntegrationPoint & bmip,
	   BareSliceVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const override;

    virtual void
    Apply (const FiniteElement & bfel,
	   const BaseMappedIntegrationRule & bmir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<double> flux,
	   LocalHeap & lh) const override;

    virtual void
    Apply (const FiniteElement & bfel,
	   const BaseMappedIntegrationPoint & bmip,
	   BareSliceVector<Complex> x, 
	   FlatVector<Complex> flux,
	   LocalHeap & lh) const override;

    virtual void
    Apply (const FiniteElement & bfel,
	   const BaseMappedIntegrationRule & bmir,
	   BareSliceVector<Complex> x, 
	   BareSliceMatrix<Complex> flux,
	   LocalHeap & lh) const override;

    virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<double> x, 
	   BareSliceMatrix<SIMD<double>> flux) const override;

    virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<Complex> x, 
	   BareSliceMatrix<SIMD<Complex>> flux) const override;


    virtual void
    ApplyTrans (const FiniteElement & bfel,
		const BaseMappedIntegrationPoint & bmip,
		FlatVector<double> flux,
		BareSliceVector<double> x, 
		LocalHeap & lh) const override;

    virtual void
    ApplyTrans (const FiniteElement & bfel,
		const BaseMappedIntegrationPoint & bmip,
		FlatVector<Complex> flux,
		BareSliceVector<Complex> x, 
		LocalHeap & lh) const override;

    virtual void
    ApplyTrans (const FiniteElement & bfel,
		const BaseMappedIntegrationRule & bmir,
		FlatMatrix<double> flux,
		BareSliceVector<double> x, 
		LocalHeap & lh) const override;

    virtual void
    ApplyTrans (const FiniteElement & bfel,
		const BaseMappedIntegrationRule & bmir,
		FlatMatrix<Complex> flux,
		BareSliceVector<Complex> x, 
		LocalHeap & lh) const override;

    virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const override;

    virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<Complex>> flux,
              BareSliceVector<Complex> x) const override;

#endif
    shared_ptr<CoefficientFunction> DiffShape (shared_ptr<CoefficientFunction> proxy,
                                               shared_ptr<CoefficientFunction> dir,
                                               bool Eulerian) const override
    {
      return DIFFOP::DiffShape(proxy, dir, Eulerian);
    }
    
  };








  
  // new design, code is still experimental ...
  
  template <typename DOP, typename F>
  class T_FunctionDiffOp : public DifferentialOperator
  {

    // possible conversion from vector to scalar 
    class V2VS 
    {
      FlatVector<> v;
      FlatVector<Complex> vc;
    public:
      V2VS (FlatVector<> av) : v(av), vc(0,(Complex*)nullptr) { ; }
      V2VS (FlatVector<Complex> avc) : v(0,(double*)nullptr), vc(avc) { ; }
      
      template <int D>
      V2VS (Vec<D> av) : v(av) { ; }
      
      operator double () { return v(0); }
      
      operator FlatVector<> () { return v; }
      
      template <int D>
      operator Vec<D> () { return v; }

      template <int D>
      operator Vec<D,Complex> () { return vc; }
    };
    
    
    int dim;
    const F & func;
  public:
    
    T_FunctionDiffOp (const F & afunc, int adim) : func(afunc), dim(adim) { ; }
    
    virtual int Dim() const { return dim; }
    virtual int DiffOrder () const { return 0; }
    virtual void Apply (const FiniteElement & fel,
			const BaseMappedIntegrationPoint & mip,
			FlatVector<double> x, 
			FlatVector<double> flux,
			LocalHeap & lh) const 
    {
      Vec<DOP::DIM_DMAT> u;
      DOP::Apply (fel, static_cast<const MappedIntegrationPoint<3,3>&> (mip), x, u, lh);
      flux = func(V2VS(u));
    }

/*
    virtual void Apply (const FiniteElement & fel,
			const BaseMappedIntegrationPoint & mip,
			FlatVector<Complex> x, 
			FlatVector<Complex> flux,
			LocalHeap & lh) const 
    {
      Vec<DOP::DIM_DMAT,Complex> u;
      DOP::Apply (fel, static_cast<const MappedIntegrationPoint<3,3>&> (mip), x, u, lh);
      flux = func(V2VS(u));
    }
*/
  };
  
  
  template <typename DOP, typename F>
  shared_ptr<DifferentialOperator> CreateFunctionDiffOp (const DOP & dop, 
                                                         const F & func, int dim = 1)
  {
    return make_shared<T_FunctionDiffOp<DOP, F>> (func, dim);
  }


  /* examples:
     
  double myexp (double x)  { return exp(x); }
  Vec<1> myexpVec (FlatVector<> x)  { return Vec<1> (exp(x(0))); }

  CreateFunctionDiffOp(DiffOpId<2>(), myexpVec));
  */







  
}


#endif
