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
    NGS_DLL_HEADER virtual ~DifferentialOperator () = default;

    virtual void DoArchive(Archive & ar) { ; }

    void SetDimensions (const Array<int> & adims) { dimensions = adims; }
    void SetVectorSpaceEmbedding (Matrix <> emb)
    { vsembedding = emb; vsdim = emb.Width(); }
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
    /// total polynomial degree is reduced by this order (i.e. minimal difforder)
    int DiffOrder() const { return difforder; } 

    virtual IntRange UsedDofs(const FiniteElement & fel) const { return IntRange(0, fel.GetNDof()); }

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
    
    NGS_DLL_HEADER virtual shared_ptr<CoefficientFunction> DiffShape (shared_ptr<CoefficientFunction> proxy,
                                                       shared_ptr<CoefficientFunction> dir,
                                                       bool Eulerian = false) const;

    NGS_DLL_HEADER virtual list<tuple<string,double>> Timing (const FiniteElement & fel, const BaseMappedIntegrationRule & mir) const;
  };

  

}



#endif
