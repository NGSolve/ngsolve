#ifndef TPDIFFOP_HPP
#define TPDIFFOP_HPP

#include "diffop.hpp"


namespace ngfem
{
  class TPDifferentialOperator : public DifferentialOperator
  {
    ArrayMem<shared_ptr<DifferentialOperator>,2> evaluators;
  public:
    NGS_DLL_HEADER TPDifferentialOperator() : DifferentialOperator(1,1,VOL,1) { ; }
    NGS_DLL_HEADER TPDifferentialOperator(FlatArray<shared_ptr<DifferentialOperator> > aevaluators,int adim,int ablockdim,bool abnd,int adifforder) : DifferentialOperator(adim,ablockdim,abnd ? BND : VOL,adifforder)
      {
        evaluators.SetSize( aevaluators.Size() );
        evaluators = aevaluators;
      }
    /// destructor
    NGS_DLL_HEADER virtual ~TPDifferentialOperator () {}

    virtual shared_ptr<DifferentialOperator> GetTrace() const override { return nullptr; }
    

    virtual IntRange UsedDofs(const FiniteElement & fel) const override { return IntRange(0, fel.GetNDof()); }

    virtual bool operator== (const TPDifferentialOperator & diffop2) const { return false; }
    
    shared_ptr<DifferentialOperator> & GetEvaluators( int num)
    {
      return evaluators[num];
    }
      using DifferentialOperator::Apply;
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
        const BaseMappedIntegrationRule & mir,
        BareSliceVector<double> x, 
        BareSliceMatrix<double> flux,
        LocalHeap & lh) const override;
       
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
        const BaseMappedIntegrationRule & mir,
        FlatMatrix<double> flux,
        BareSliceVector<double> x, 
        LocalHeap & lh) const override;
    
    NGS_DLL_HEADER virtual void
    ApplyY(const FiniteElement &fel,
        const BaseMappedIntegrationRule & miry,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;
        
    NGS_DLL_HEADER virtual void
    ApplyYTrans(const FiniteElement &fel,
        const BaseMappedIntegrationRule & miry,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    ApplyX(const FiniteElement &fel,
        const BaseMappedIntegrationRule & mirx,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;
        
    NGS_DLL_HEADER virtual void
    ApplyXTrans(const FiniteElement &fel,
        const BaseMappedIntegrationRule & mirx,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;

    /// calculates the matrix
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                BareSliceMatrix<double,ColMajor> mat,   
                LocalHeap & lh) const override { throw Exception("2 not implemented"); }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & bfel,
                const BaseMappedIntegrationPoint & bmip,
                BareSliceMatrix<Complex,ColMajor> mat, 
                LocalHeap & lh) const override{ throw Exception("3 not implemented"); }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                BareSliceMatrix<double,ColMajor> mat,   
                LocalHeap & lh) const override { throw Exception("4 not implemented"); }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                BareSliceMatrix<Complex,ColMajor> mat,   
                LocalHeap & lh) const override { throw Exception("5 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationPoint & mip,
           BareSliceVector<double> x, 
           FlatVector<double> flux,
           LocalHeap & lh) const override { throw Exception("6 not implemented"); }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationPoint & mip,
           BareSliceVector<Complex> x, 
           FlatVector<Complex> flux,
           LocalHeap & lh) const override { throw Exception("8 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationRule & mir,
           BareSliceVector<Complex> x, 
           BareSliceMatrix<Complex> flux,
           LocalHeap & lh) const override { throw Exception("9 not implemented"); }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
        const SIMD_BaseMappedIntegrationRule & bmir,
        BareSliceVector<double> x, 
        BareSliceMatrix<SIMD<double>> flux) const override { throw ExceptionNOSIMD("nosimd"); }
    // LocalHeap & lh) const { throw Exception("not implemented"); }
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<double> flux,
                BareSliceVector<double> x, 
                LocalHeap & lh) const override { throw Exception("11 not implemented"); }

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<Complex> flux,
                BareSliceVector<Complex> x, 
                LocalHeap & lh) const override { throw Exception("12 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                FlatMatrix<Complex> flux,
                BareSliceVector<Complex> x, 
                LocalHeap & lh) const override { throw Exception("13 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const override { throw ExceptionNOSIMD("nosimd"); }
    // LocalHeap & lh) const;
  };

  
  class TPBlockDifferentialOperator : public BlockDifferentialOperator
  {
    //ArrayMem<shared_ptr<DifferentialOperator>,2> evaluators;
    //shared_ptr<DifferentialOperator> diffop;
    //int blockdim;
  public:
    //NGS_DLL_HEADER TPBlockDifferentialOperator() : DifferentialOperator(1,1,VOL,1) { ; }
    NGS_DLL_HEADER TPBlockDifferentialOperator(shared_ptr<DifferentialOperator> adiffop,int ablockdim) : /*diffop(adiffop), blockdim(ablockdim),*/ BlockDifferentialOperator(adiffop,ablockdim)
      { ; }
    /// destructor
    NGS_DLL_HEADER virtual ~TPBlockDifferentialOperator () {}

    virtual IntRange UsedDofs(const FiniteElement & fel) const override { return IntRange(0, fel.GetNDof()); }

    virtual bool operator== (const DifferentialOperator & diffop2) const override { return false; }
    
    //int BlockDim() const { return blockdim; }
    
    shared_ptr<DifferentialOperator> & GetEvaluators( int num)
    {
      return dynamic_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(num);
    }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationRule & mir,
           BareSliceVector<double> x, 
           BareSliceMatrix<double> flux,
           LocalHeap & lh) const override;
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                FlatMatrix<double> flux,
                BareSliceVector<double> x, 
                LocalHeap & lh) const override;
    
    NGS_DLL_HEADER virtual void
    ApplyY(const FiniteElement &fel,
        const BaseMappedIntegrationRule & miry,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;
        
    NGS_DLL_HEADER virtual void
    ApplyYTrans(const FiniteElement &fel,
        const BaseMappedIntegrationRule & miry,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    ApplyX(const FiniteElement &fel,
        const BaseMappedIntegrationRule & mirx,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;
        
    NGS_DLL_HEADER virtual void
    ApplyXTrans(const FiniteElement &fel,
                const BaseMappedIntegrationRule & mirx,
                FlatMatrix<double> flux,
                SliceMatrix<double> x,
                LocalHeap & lh) const;

    /// calculates the matrix
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                BareSliceMatrix<double,ColMajor> mat,   
                LocalHeap & lh) const override { throw Exception("2 not implemented"); }
    
      using BlockDifferentialOperator::CalcMatrix;
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & bfel,
                const BaseMappedIntegrationPoint & bmip,
                BareSliceMatrix<Complex,ColMajor> mat, 
                LocalHeap & lh) const override { throw Exception("3 not implemented"); }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                BareSliceMatrix<double,ColMajor> mat,   
                LocalHeap & lh) const override { throw Exception("4 not implemented"); }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                BareSliceMatrix<Complex,ColMajor> mat,   
                LocalHeap & lh) const override { throw Exception("5 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationPoint & mip,
           BareSliceVector<double> x, 
           FlatVector<double> flux,
           LocalHeap & lh) const override { throw Exception("6 not implemented"); }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationPoint & mip,
           BareSliceVector<Complex> x, 
           FlatVector<Complex> flux,
           LocalHeap & lh) const override { throw Exception("8 not implemented"); }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationRule & mir,
           BareSliceVector<Complex> x, 
           BareSliceMatrix<Complex> flux,
           LocalHeap & lh) const override { throw Exception("9 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
        const SIMD_BaseMappedIntegrationRule & bmir,
        BareSliceVector<double> x, 
        BareSliceMatrix<SIMD<double>> flux) const override { throw ExceptionNOSIMD("nosimd"); }
    // LocalHeap & lh) const { throw Exception("not implemented"); }
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<double> flux,
                BareSliceVector<double> x, 
                LocalHeap & lh) const override { throw Exception("11 not implemented"); }

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<Complex> flux,
                BareSliceVector<Complex> x, 
                LocalHeap & lh) const override { throw Exception("12 not implemented"); }

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                FlatMatrix<Complex> flux,
                BareSliceVector<Complex> x, 
                LocalHeap & lh) const override { throw Exception("13 not implemented"); }

    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const override { throw ExceptionNOSIMD("nosimd"); }
    // LocalHeap & lh) const;
  };






  class TPBlockDifferentialOperator2 : public BlockDifferentialOperator
  {
    //ArrayMem<shared_ptr<DifferentialOperator>,2> evaluators;
    //shared_ptr<DifferentialOperator> diffop;
    //int blockdim;
  public:
    //NGS_DLL_HEADER TPBlockDifferentialOperator() : DifferentialOperator(1,1,VOL,1) { ; }
    NGS_DLL_HEADER TPBlockDifferentialOperator2(shared_ptr<DifferentialOperator> adiffop,int ablockdim) : /*diffop(adiffop), blockdim(ablockdim),*/ BlockDifferentialOperator(adiffop,ablockdim)
      { ; }
    /// destructor
    NGS_DLL_HEADER virtual ~TPBlockDifferentialOperator2 () {}

    virtual IntRange UsedDofs(const FiniteElement & fel) const override { return IntRange(0, fel.GetNDof()); }

    virtual bool operator== (const DifferentialOperator & diffop2) const override { return false; }
    
    //int BlockDim() const { return blockdim; }
    
    shared_ptr<DifferentialOperator> & GetEvaluators( int num)
    {
      return dynamic_cast<TPDifferentialOperator*>(diffop.get())->GetEvaluators(num);
    }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationRule & mir,
           BareSliceVector<double> x, 
           BareSliceMatrix<double> flux,
           LocalHeap & lh) const override;
       
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                FlatMatrix<double> flux,
                BareSliceVector<double> x, 
                LocalHeap & lh) const override;
    
    NGS_DLL_HEADER virtual void
    ApplyY(const FiniteElement &fel,
        const BaseMappedIntegrationRule & miry,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;
        
    NGS_DLL_HEADER virtual void
    ApplyYTrans(const FiniteElement &fel,
        const BaseMappedIntegrationRule & miry,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void
    ApplyX(const FiniteElement &fel,
        const BaseMappedIntegrationRule & mirx,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;
        
    NGS_DLL_HEADER virtual void
    ApplyXTrans(const FiniteElement &fel,
        const BaseMappedIntegrationRule & mirx,
        FlatMatrix<double> flux,
        SliceMatrix<double> x,
        LocalHeap & lh) const;

    /// calculates the matrix
      using BlockDifferentialOperator::CalcMatrix;
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                BareSliceMatrix<double,ColMajor> mat,   
                LocalHeap & lh) const override { throw Exception("2 not implemented"); }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & bfel,
                const BaseMappedIntegrationPoint & bmip,
                BareSliceMatrix<Complex,ColMajor> mat, 
                LocalHeap & lh) const override { throw Exception("3 not implemented"); }

    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                BareSliceMatrix<double,ColMajor> mat,   
                LocalHeap & lh) const override { throw Exception("4 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    CalcMatrix (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                BareSliceMatrix<Complex,ColMajor> mat,   
                LocalHeap & lh) const override { throw Exception("5 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationPoint & mip,
           BareSliceVector<double> x, 
           FlatVector<double> flux,
           LocalHeap & lh) const override { throw Exception("6 not implemented"); }

    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationPoint & mip,
           BareSliceVector<Complex> x, 
           FlatVector<Complex> flux,
           LocalHeap & lh) const override { throw Exception("8 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & fel,
           const BaseMappedIntegrationRule & mir,
           BareSliceVector<Complex> x, 
           BareSliceMatrix<Complex> flux,
           LocalHeap & lh) const override{ throw Exception("9 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    Apply (const FiniteElement & bfel,
           const SIMD_BaseMappedIntegrationRule & bmir,
           BareSliceVector<double> x, 
           BareSliceMatrix<SIMD<double>> flux) const override { throw ExceptionNOSIMD("nosimd"); }
    // LocalHeap & lh) const { throw Exception("not implemented"); }
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<double> flux,
                BareSliceVector<double> x, 
                LocalHeap & lh) const override { throw Exception("11 not implemented"); }

    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationPoint & mip,
                FlatVector<Complex> flux,
                BareSliceVector<Complex> x, 
                LocalHeap & lh) const override { throw Exception("12 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    ApplyTrans (const FiniteElement & fel,
                const BaseMappedIntegrationRule & mir,
                FlatMatrix<Complex> flux,
                BareSliceVector<Complex> x, 
                LocalHeap & lh) const override { throw Exception("13 not implemented"); }
    
    NGS_DLL_HEADER virtual void
    AddTrans (const FiniteElement & bfel,
              const SIMD_BaseMappedIntegrationRule & bmir,
              BareSliceMatrix<SIMD<double>> flux,
              BareSliceVector<double> x) const override { throw ExceptionNOSIMD("nosimd"); }
    // LocalHeap & lh) const;
  };
}
#endif // TPDIFFOP_HPP
