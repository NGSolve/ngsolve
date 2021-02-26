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
                      const CoefficientFunction & func, BareSliceMatrix<> coefs, LocalHeap & lh) const override
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
    static bool SupportsVB (VorB checkvb) { return true; }

    template <typename MIP, typename MAT>
    static void GenerateMatrix (const FiniteElement & fel, const MIP & mip,
				MAT && mat, LocalHeap & lh)
    {
      mat.Row(0).Range(0, fel.GetNDof()) = 0;
      mat(0, mip.IP().Nr()) = 1;
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
        IntegrationRule ir(ma->GetElType( { VOL, i } ), 2*order);
        ndof += ir.Size();
      }
    firsteldof.Last() = ndof;
    //cout << "firstel = " << firsteldof << endl;
    SetNDof(ndof);
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

  
  static RegisterFESpace<IntegrationRuleSpace> init ("irspace");
}
