#include <comp.hpp>

namespace ngcomp
{

  class NumberFiniteElement : public FiniteElement
  {
    ELEMENT_TYPE et;
  public:
    NumberFiniteElement (ELEMENT_TYPE _et)
      : FiniteElement(1, 0) { ; }
    HD virtual ELEMENT_TYPE ElementType() const { return et; }
  };


  class NumberDiffOp : public DiffOp<NumberDiffOp>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 0 };
    enum { DIM_ELEMENT = 0 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };
  
    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      mat(0,0) = 1;
    }

    using DiffOp<NumberDiffOp>::ApplySIMDIR;
    /*
    template <typename FEL, class MIR, class TVX, class TVY>
    static void ApplySIMDIR (const FEL & fel, const MIR & mir,
                             const TVX & x, TVY & y)
    */
    static void ApplySIMDIR (const FiniteElement & fel, const SIMD_BaseMappedIntegrationRule & mir,
                             BareSliceVector<double> x, BareSliceMatrix<SIMD<double>> y)
    {
      for (size_t i = 0; i < mir.Size(); i++)
        y(0,i) = x(0);
    }

    using DiffOp<NumberDiffOp>::AddTransSIMDIR;
    /// Computes Transpose (B-matrix) times point value    
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    {
      SIMD<double> sum = 0.0;
      for (size_t i = 0; i < mir.Size(); i++)
        sum += x(0,i);
      y(0) += HSum(sum);
    }
    
  };





  class NumberFESpace : public FESpace
  {
  public:
    NumberFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false)
      : FESpace (ama, flags)
    { 
      evaluator[VOL] = make_shared<T_DifferentialOperator<NumberDiffOp>>();
      evaluator[BND] = make_shared<T_DifferentialOperator<NumberDiffOp>>();
      is_atomic_dof = BitArray(1);
      is_atomic_dof = true;
    }

    virtual size_t GetNDof() const { return 1; }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const
    {
      if (DefinedOn(ei))
        return *new (lh) NumberFiniteElement(ma->GetElType(ei));
      else
        return *new (lh) DummyFE<ET_POINT>();
    }

    virtual void GetDofNrs (ElementId ei, Array<int> & dnums) const
    {
      if (DefinedOn(ei))
        {
          dnums.SetSize(1);
          dnums[0] = 0;
        }
      else
        dnums.SetSize(0);
    }

  };


  static RegisterFESpace<NumberFESpace> init ("number");

}
