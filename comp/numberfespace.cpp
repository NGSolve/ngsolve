#include <comp.hpp>
#include "numberfespace.hpp"

namespace ngcomp
{

  class NumberFiniteElement : public FiniteElement
  {
    ELEMENT_TYPE et;
  public:
    NumberFiniteElement (ELEMENT_TYPE _et)
      : FiniteElement(1, 0), et(_et) { ; }
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
    
    static bool SupportsVB (VorB checkvb) { return true; }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      if(fel.GetNDof()) // if not dummyFE
        mat(0,0) = 1;
    }

    static void GenerateMatrixSIMDIR (const FiniteElement & bfel,
                                      const SIMD_BaseMappedIntegrationRule & mir,
                                      BareSliceMatrix<SIMD<double>> mat)
    {
      if(bfel.GetNDof()) // if not dummyFE
        mat.Row(0).AddSize(mir.Size()) = SIMD<double>(1);
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
      if(fel.GetNDof()) // if not dummyFE
        for (size_t i = 0; i < mir.Size(); i++)
          y(0,i) = x(0);
    }

    using DiffOp<NumberDiffOp>::AddTransSIMDIR;
    /// Computes Transpose (B-matrix) times point value    
    static void AddTransSIMDIR (const FiniteElement & bfel, const SIMD_BaseMappedIntegrationRule & mir,
                                BareSliceMatrix<SIMD<double>> x, BareSliceVector<double> y)
    {
      if(bfel.GetNDof() == 0) // if not dummyFE
        return;
      SIMD<double> sum = 0.0;
      for (size_t i = 0; i < mir.Size(); i++)
        sum += x(0,i);
      y(0) += HSum(sum);
    }
    
  };





  NumberFESpace::NumberFESpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags)
      : FESpace (ama, flags)
    {
      type = "number";
      evaluator[VOL] = make_shared<T_DifferentialOperator<NumberDiffOp>>();
      evaluator[BND] = make_shared<T_DifferentialOperator<NumberDiffOp>>();
      evaluator[BBND] = make_shared<T_DifferentialOperator<NumberDiffOp>>();
      evaluator[BBBND] = make_shared<T_DifferentialOperator<NumberDiffOp>>();


    if (dimension > 1)
      {
        for (auto vb : { VOL,BND, BBND, BBBND })
          evaluator[vb] = make_shared<BlockDifferentialOperator> (evaluator[vb], dimension);
      }
      
      is_atomic_dof = BitArray(1);
      is_atomic_dof = true;
    }

  void NumberFESpace::Update()
    {
      SetNDof( ((ma->GetCommunicator().Size() > 1) && (ma->GetCommunicator().Rank() == 0)) ? 0 : 1 );
    }

  
  FiniteElement & NumberFESpace::GetFE (ElementId ei, Allocator & lh) const
    {
      if (DefinedOn(ei))
        return *new (lh) NumberFiniteElement(ma->GetElType(ei));
      else
        return *new (lh) DummyFE<ET_POINT>();
    }

  void NumberFESpace::GetDofNrs (ElementId ei, Array<int> & dnums) const
    {
      if (DefinedOn(ei))
        {
          dnums.SetSize(1);
          dnums[0] = 0;
        }
      else
        dnums.SetSize(0);
    }


  void NumberFESpace :: GetGlobalDofNrs (int gnr, Array<int> & dnums) const
  {
    /** If parallel, there is one global DOF shared by all ranks except rank 0.
	If not parallel, there is also one DOF. **/
    if ( IsParallel() && (GetParallelDofs()->GetCommunicator().Size() > 1) && (GetParallelDofs()->GetCommunicator().Rank() == 0) )
      { dnums.SetSize0(); }
    else
      { dnums.SetSize(1); dnums[0] = 0; }
  }


  static RegisterFESpace<NumberFESpace> init ("number");

}
