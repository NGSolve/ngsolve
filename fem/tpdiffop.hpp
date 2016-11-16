#ifndef TPDIFFOP_HPP
#define TPDIFFOP_HPP

namespace ngfem
{
  class TPDifferentialOperator : public DifferentialOperator
  {
    ArrayMem<shared_ptr<DifferentialOperator>,2> evaluators;
  public:
    NGS_DLL_HEADER TPDifferentialOperator() { ; }
    NGS_DLL_HEADER TPDifferentialOperator(FlatArray<shared_ptr<DifferentialOperator> > aevaluators)
      {
        evaluators.SetSize( aevaluators.Size() );
        evaluators = aevaluators;
      }
    /// destructor
    NGS_DLL_HEADER virtual ~TPDifferentialOperator () {}
    /// dimension of range
    NGS_DLL_HEADER virtual int Dim() const
    {
        int dim = 0;
        for (auto eval : evaluators)
          dim = max2(dim, eval->Dim());
        return dim;
    }
    /// number of copies of finite element by BlockDifferentialOperator
    NGS_DLL_HEADER virtual int BlockDim() const { return 1; }
    /// does it live on the boundary ?
    virtual bool Boundary() const { return false; }

    /// total polynomial degree is reduced by this order (i.e. minimal difforder)
    virtual int DiffOrder() const
    {
        int ord = evaluators[0]->DiffOrder();
        for (auto eval : evaluators)
          ord = min2(ord, eval->DiffOrder());
        return ord;
    }

    virtual IntRange UsedDofs(const FiniteElement & fel) const { return IntRange(0, fel.GetNDof()); }

    virtual bool operator== (const TPDifferentialOperator & diffop2) const { return false; }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationRule & mir,
	   FlatVector<double> x, 
	   FlatMatrix<double> flux,
	   LocalHeap & lh) const;
       
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<double> flux,
		FlatVector<double> x, 
		LocalHeap & lh) const;

    /// calculates the matrix
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		SliceMatrix<double,ColMajor> mat,   
		LocalHeap & lh) const { throw Exception("2 not implemented"); }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & bfel,
		const BaseMappedIntegrationPoint & bmip,
		SliceMatrix<Complex,ColMajor> mat, 
		LocalHeap & lh) const { throw Exception("3 not implemented"); }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		SliceMatrix<double,ColMajor> mat,   
		LocalHeap & lh) const { throw Exception("4 not implemented"); }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		SliceMatrix<Complex,ColMajor> mat,   
		LocalHeap & lh) const { throw Exception("5 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   FlatVector<double> x, 
	   FlatVector<double> flux,
	   LocalHeap & lh) const { throw Exception("6 not implemented"); }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationPoint & mip,
	   FlatVector<Complex> x, 
	   FlatVector<Complex> flux,
	   LocalHeap & lh) const { throw Exception("8 not implemented"); }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
	   const BaseMappedIntegrationRule & mir,
	   FlatVector<Complex> x, 
	   FlatMatrix<Complex> flux,
	   LocalHeap & lh) const { throw Exception("9 not implemented"); }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
	   const SIMD_BaseMappedIntegrationRule & bmir,
	   BareSliceVector<double> x, 
	   ABareMatrix<double> flux) const { throw ExceptionNOSIMD("nosimd"); }
    // LocalHeap & lh) const { throw Exception("not implemented"); }
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		FlatVector<double> flux,
		FlatVector<double> x, 
		LocalHeap & lh) const { throw Exception("11 not implemented"); }

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationPoint & mip,
		FlatVector<Complex> flux,
		FlatVector<Complex> x, 
		LocalHeap & lh) const { throw Exception("12 not implemented"); }

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
		const BaseMappedIntegrationRule & mir,
		FlatMatrix<Complex> flux,
		FlatVector<Complex> x, 
		LocalHeap & lh) const { throw Exception("14 not implemented"); }

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              ABareMatrix<double> flux,
              BareSliceVector<double> x) const { throw ExceptionNOSIMD("nosimd"); }
    // LocalHeap & lh) const;
  };

}
#endif // TPDIFFOP_HPP
