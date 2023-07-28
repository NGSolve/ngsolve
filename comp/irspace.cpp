#include <comp.hpp>
#include "irspace.hpp"

namespace ngcomp
{

  class IRFiniteElement : public FiniteElement
  {
    ELEMENT_TYPE et;
    IntegrationRule ir;
  public:
    IRFiniteElement (ELEMENT_TYPE aet, int aorder)
      : FiniteElement(0, aorder), et(aet), ir(aet, 2*aorder)
    {
      ndof = ir.Size();
    }
    ELEMENT_TYPE ElementType() const override { return et; }
    
    void Interpolate (const ElementTransformation & trafo, 
                      const CoefficientFunction & func, SliceMatrix<> coefs, LocalHeap & lh) const override
    {
      HeapReset hr(lh);
      BaseMappedIntegrationRule & mir = trafo(ir, lh);
      func.Evaluate (mir, coefs);
    }

  };


  class IRDiffOp : public DiffOp<IRDiffOp>
  {
  public:
    enum { DIM = 1 };
    enum { DIM_SPACE = 0 };
    enum { DIM_ELEMENT = 0 };
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 0 };

    static INT<0> GetDimensions() { return INT<0>(); };
    static string Name() { return "IdIR"; }
    static bool SupportsVB (VorB checkvb) { return true; }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      mat.Row(0).Range(0, fel.GetNDof()) = 0;
      mat(0, mip.IP().Nr()) = 1;
    }

    using DiffOp<IRDiffOp>::ApplyIR;
    template <class MIR, class TMY>
    static void ApplyIR (const FiniteElement & fel, const MIR & mir,
                         BareSliceVector<double> x, TMY y,
			 LocalHeap & lh)
    {
      y.Col(0).Range(0, fel.GetNDof()) = x;
    }

    using DiffOp<IRDiffOp>::ApplyTransIR;    
    template <class MIR>
    static void ApplyTransIR (const FiniteElement & fel, 
			      const MIR & mir,
			      FlatMatrix<double> x, BareSliceVector<double> y,
			      LocalHeap & lh)
    {
      y.Range(0, fel.GetNDof()) = x.Col(0);
    }
    
    
    template <typename FEL, typename MIR>
    static void GenerateMatrixSIMDIR (const FEL & fel, const MIR & mir, BareSliceMatrix<SIMD<double>> mat)
    {
      int ndof = fel.GetNDof();
      mat.AddSize(ndof, mir.Size()) = SIMD<double> (0.0);
      SliceMatrix<double> hmat(ndof, ndof, SIMD<double>::Size()*mat.Dist(), (double*)mat.Data());
      hmat.Diag() = 1;
      // cout << "diffopir, genmat:" << endl << mat.AddSize(ndof, mir.Size()) << endl;
    }
    
    using DiffOp<IRDiffOp> :: ApplySIMDIR;
    template <typename FEL, class MIR, class TVX>
    static void ApplySIMDIR (const FEL & fel, const MIR & mir,
                             const TVX & x, BareSliceMatrix<SIMD<double>> y)
    {
      int i = 0;
      constexpr int SW = SIMD<double>::Size();
      int ndof = fel.GetNDof();
      // range of safe access
      for ( ; (i + 1) * SW <= ndof; i++)
        y(0, i) = SIMD<double>(&x(SW*i));

      // handle possible overhead
      if (const int mask = ndof - SW * i; mask != 0)
        y(0, i) = If(SIMD<mask64>(mask), SIMD<double>(&x(SW*i, mask)), SIMD<double>(x(ndof - 1)));

      // cout << "y (final) " << y << endl;
    }

    using DiffOp<IRDiffOp> :: AddTransSIMDIR;    
    template <typename FEL, class MIR, class TVY>
    static void AddTransSIMDIR (const FEL & fel, const MIR & mir,
                                BareSliceMatrix<SIMD<double>> x, TVY & y)
    // LocalHeap & lh)
    {
      int i = 0;
      constexpr int SW = SIMD<double>::Size();
      int ndof = fel.GetNDof();
      for ( ; i*SW < ndof; i++)
        x(0,i).Store(&y(SW*i), ndof-SW*i);
    }
    
    
  };





  IntegrationRuleSpace::IntegrationRuleSpace (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags)
      : FESpace (ama, flags)
    {
      type = "irspace";
      evaluator[VOL] = make_shared<T_DifferentialOperator<IRDiffOp>>();
      
      if (dimension > 1)
        evaluator[VOL] = make_shared<BlockDifferentialOperator> (evaluator[VOL], dimension);
    }

  void IntegrationRuleSpace::Update()
  {
    //cout << "update irspace" << endl;
    firsteldof.SetSize(ma->GetNE()+1);
    size_t ndof = 0;
    for (auto i : Range(ma->GetNE()))
      {
        firsteldof[i] = ndof;
        if (!DefinedOn(ElementId(VOL, i)))
            continue;

        IntegrationRule ir(ma->GetElType( { VOL, i } ), 2*order);
        ndof += ir.Size();
      }
    firsteldof.Last() = ndof;
    //cout << "firstel = " << firsteldof << endl;
    SetNDof(ndof);

    UpdateCouplingDofArray();
  }

  void IntegrationRuleSpace::UpdateCouplingDofArray()
  {
    // all dofs are local dofs
    ctofdof.SetSize(ndof);
    ctofdof = LOCAL_DOF;
  }
  
  
  FiniteElement & IntegrationRuleSpace::GetFE (ElementId ei, Allocator & lh) const
  {
    if (ei.VB() == VOL && DefinedOn(ei))
      return *new (lh) IRFiniteElement(ma->GetElType(ei), order);
    else
      return SwitchET (ma->GetElType(ei), [&] (auto et) -> FiniteElement&
                       {
                         return *new (lh) DummyFE<et.ElementType()> ();
                       });
  }
  
  void IntegrationRuleSpace::GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    if (ei.VB() == VOL)
      dnums = IntRange(firsteldof[ei.Nr()], firsteldof[ei.Nr()+1]);
    else
      dnums.SetSize0();
  }

  std::map<ELEMENT_TYPE, IntegrationRule> IntegrationRuleSpace::GetIntegrationRules() const
  {
    std::map<ELEMENT_TYPE, IntegrationRule> irs;
    irs[ET_TRIG] = IntegrationRule(ET_TRIG, 2*order);
    irs[ET_QUAD] = IntegrationRule(ET_QUAD, 2*order);
    irs[ET_HEX] = IntegrationRule(ET_HEX, 2*order);
    irs[ET_TET] = IntegrationRule(ET_TET, 2*order);
    irs[ET_PRISM] = IntegrationRule(ET_PRISM, 2*order);
    irs[ET_PYRAMID] = IntegrationRule(ET_PYRAMID, 2*order);
    return irs;
  }




    IntegrationRuleSpaceSurface::IntegrationRuleSpaceSurface (shared_ptr<MeshAccess> ama, const Flags & flags, bool checkflags)
      : FESpace (ama, flags)
    {
      type = "irspacesurface";
      evaluator[VOL] = make_shared<T_DifferentialOperator<IRDiffOp>>();
      evaluator[BND] = make_shared<T_DifferentialOperator<IRDiffOp>>();
      
      if (dimension > 1)
        {
          evaluator[VOL] = make_shared<BlockDifferentialOperator> (evaluator[VOL], dimension);
          evaluator[BND] = make_shared<BlockDifferentialOperator> (evaluator[BND], dimension);
        }
    }

  void IntegrationRuleSpaceSurface::Update()
  {
    //cout << "update irspace" << endl;
    firsteldof.SetSize(ma->GetNSE()+1);
    size_t ndof = 0;
    for (auto i : Range(ma->GetNSE()))
      {
        firsteldof[i] = ndof;
        if (!DefinedOn(ElementId(BND, i)))
            continue;

        IntegrationRule ir(ma->GetElType( { BND, i } ), 2*order);
        ndof += ir.Size();
      }
    firsteldof.Last() = ndof;
    //cout << "firstel = " << firsteldof << endl;
    SetNDof(ndof);

    UpdateCouplingDofArray();
  }

  void IntegrationRuleSpaceSurface::UpdateCouplingDofArray()
  {
    // all dofs are local dofs
    ctofdof.SetSize(ndof);
    ctofdof = LOCAL_DOF;
  }
  
  
  FiniteElement & IntegrationRuleSpaceSurface::GetFE (ElementId ei, Allocator & lh) const
  {
    if (ei.VB() == BND && DefinedOn(ei))
      return *new (lh) IRFiniteElement(ma->GetElType(ei), order);
    else
      return SwitchET (ma->GetElType(ei), [&] (auto et) -> FiniteElement&
                       {
                         return *new (lh) DummyFE<et.ElementType()> ();
                       });
  }
  
  void IntegrationRuleSpaceSurface::GetDofNrs (ElementId ei, Array<int> & dnums) const
  {
    if (ei.VB() == BND)
      dnums = IntRange(firsteldof[ei.Nr()], firsteldof[ei.Nr()+1]);
    else
      dnums.SetSize0();
  }

  std::map<ELEMENT_TYPE, IntegrationRule> IntegrationRuleSpaceSurface::GetIntegrationRules() const
  {
    std::map<ELEMENT_TYPE, IntegrationRule> irs;
    irs[ET_SEGM] = IntegrationRule(ET_SEGM, 2*order);
    irs[ET_TRIG] = IntegrationRule(ET_TRIG, 2*order);
    irs[ET_QUAD] = IntegrationRule(ET_QUAD, 2*order);
    return irs;
  }

  
  static RegisterFESpace<IntegrationRuleSpace> init ("irspace");
  static RegisterFESpace<IntegrationRuleSpaceSurface> initsurf ("irspacesurface");
}
