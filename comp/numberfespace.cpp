#include <comp.hpp>

namespace ngcomp
{

  class NumberFiniteElement : public FiniteElement
  {
  public:
    NumberFiniteElement ()
      : FiniteElement(1, 0) { ; }
    HD virtual ELEMENT_TYPE ElementType() const { return ET_POINT; }
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
                             BareSliceVector<double> x, ABareSliceMatrix<double> y)
    {
      for (int i = 0; i < mir.IR().GetNIP(); i++)
        y(0,i) = x(0);
    }

    using DiffOp<NumberDiffOp>::AddTransSIMDIR;
    /// Computes Transpose (B-matrix) times point value    
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                ABareSliceMatrix<double> x, BareSliceVector<double> y)
    {
      double sum = 0.0;
      for (int i = 0; i < mir.IR().GetNIP(); i++)
        sum += x(0,i);
      y(0) += sum;
    }
    
  };





  class NumberFESpace : public FESpace
  {
  public:
    NumberFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags=false)
      : FESpace (ama, flags)
    { 
      evaluator = make_shared<T_DifferentialOperator<NumberDiffOp>>();
      boundary_evaluator = make_shared<T_DifferentialOperator<NumberDiffOp>>();
      is_atomic_dof = BitArray(1);
      is_atomic_dof = true;
    }

    virtual int GetNDof() const { return 1; }

    virtual FiniteElement & GetFE (ElementId ei, Allocator & lh) const
    {
      if (DefinedOn(ei))
        return *new (lh) NumberFiniteElement();
      else
        return *new (lh) DummyFE<ET_POINT>();
    }

    virtual void GetDofNrs (int elnr, Array<int> & dnums) const
    {
      if (DefinedOn(ElementId(VOL,elnr)))
        {
          dnums.SetSize(1);
          dnums[0] = 0;
        }
      else
        dnums.SetSize(0);
    }

    virtual void GetSDofNrs (int selnr, Array<int> & dnums) const
    {
      if (DefinedOn(ElementId(BND,selnr)))
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
