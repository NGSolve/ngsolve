#include <comp.hpp>    // provides FESpace, ...
#include "../fem/tscalarfe_impl.hpp"


/*
HIGHER ORDER TRIANGULAR FINITE ELEMENTS WITH MASS LUMPING FOR THE WAVE EQUATION
G. COHEN, P. JOLY, J. E. ROBERTS, AND N. TORDJMAN
SIAM J. NUMER. ANAL. ⃝c 2001 Society for Industrial and Applied Mathematics Vol. 38, No. 6, pp. 2047–2078
*/


namespace ngcomp
{


  class H1LumpingTrig2 : public T_ScalarFiniteElementFO<H1LumpingTrig2,ET_TRIG,7,3>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
    {
      Tx l1 = ip.x;
      Tx l2 = ip.y;
      Tx l3 = 1-l1-l2;
      Tx bub = 27*l1*l2*l3;
      shape[0] = 2*l1*(l1-0.5) + 1.0/9 * bub ;
      shape[1] = 2*l2*(l2-0.5) + 1.0/9 * bub ;
      shape[2] = 2*l3*(l3-0.5) + 1.0/9 * bub ;
      shape[3] = 4*l1*l3 - 4.0/9 * bub;
      shape[4] = 4*l2*l3 - 4.0/9 * bub;
      shape[5] = 4*l1*l2 - 4.0/9 * bub;
      shape[6] = bub;
    }
  };

  // template class T_ScalarFiniteElementFO<FE_Segm0,ET_SEGM,1,0>;
  
  H1LumpingFESpace :: H1LumpingFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    type = "h1lumpingfespace";
    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
  }

  DocInfo H1LumpingFESpace :: GetDocu()
  {
    auto docu = FESpace::GetDocu();
    docu.short_docu = "H1-FESpace with nodal basis for mass lumping.";
    docu.long_docu =
      R"raw_string(at the moment only for second order + bub on trigs.
)raw_string";      
    // docu.Arg("secondorder") = "bool = False\n"
    // "  Use second order basis functions";
    return docu;
  }

  void H1LumpingFESpace :: Update()
  {
    nvert = ma->GetNV();
    nedge = ma->GetNEdges();
    SetNDof (nvert+nedge+ma->GetNE(VOL));
  }

  void H1LumpingFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize(0);
    dnums += ma->GetElement(ei).Vertices();
    for (auto e : ma->GetElement(ei).Edges())
      dnums.Append (nvert+e);
    dnums.Append(nvert+nedge+ei.Nr());
  }

  /*
    Allocate finite element class, using custom allocator alloc
  */
  FiniteElement & H1LumpingFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    switch (ma->GetElement(ei).GetType())
      {
        case ET_TRIG:
          return * new (alloc) H1LumpingTrig2;
      default:
          throw Exception("H1Lumping: Element of type "+ToString(ma->GetElement(ei).GetType()) + 
                          " not available\n");
      }
  }

  std::map<ELEMENT_TYPE, IntegrationRule> H1LumpingFESpace :: GetIntegrationRules() const
  {
    IntegrationRule ir7;
    ir7.Append ( IntegrationPoint( 1, 0, 0, 1.0/40));
    ir7.Append ( IntegrationPoint( 0, 1, 0, 1.0/40));
    ir7.Append ( IntegrationPoint( 0, 0, 0, 1.0/40));
    ir7.Append ( IntegrationPoint( 0.5, 0, 0, 1.0/15));
    ir7.Append ( IntegrationPoint( 0, 0.5, 0, 1.0/15));
    ir7.Append ( IntegrationPoint( 0.5, 0.5, 0, 1.0/15));
    ir7.Append ( IntegrationPoint( 1.0/3, 1.0/3, 0, 9.0/40));
    std::map<ELEMENT_TYPE, IntegrationRule> rules;
    rules[ET_TRIG] = move(ir7);
    return rules;
  }
  
}
