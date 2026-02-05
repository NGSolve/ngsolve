#include "JKMspace.hpp"
#include "../fem/scalarfe.hpp"
#include "../fem/hdivdivfe.hpp"
#include <bdbequations.hpp>
#include <diffop_impl.hpp>
#include <hdivdiv_equations.hpp>





namespace ngcomp
{

  // ***************************** JKM FESpace *****************************

  
  JKM_FESpace::JKM_FESpace(shared_ptr<MeshAccess> ama, const Flags &flags, bool checkflags)
    : FESpace(ama, flags)
  {
    order = int(flags.GetNumFlag("order", 1));

    evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpIdHDivDiv<2>>>();
    flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpDivHDivDiv<2>>>();
  } 

  DocInfo JKM_FESpace::GetDocu()
  {
    DocInfo docu = FESpace::GetDocu();
    docu.short_docu = "Johnson–Krizek–Mercier Finite Element Space";
    return docu;
  }


  std::map<ELEMENT_TYPE, IntegrationRule> JKM_FESpace::GetIntegrationRules(int bonus_intorder) const
  {
    IntegrationRule irtrig(ET_TRIG, 2*order+bonus_intorder);

    IntegrationRule ir;

    auto map1 = [](IntegrationPoint ip)
    {
      double x = ip(0), y = ip(1);
      return IntegrationPoint(x+1./3*y, 1./3*y, 0, ip.Weight()/3);
    };

    auto map2 = [](IntegrationPoint ip)
    {
      double x = ip(0), y = ip(1);
      return IntegrationPoint(1./3*x, y+1./3*x, 0, ip.Weight()/3);
    };
    
    auto map3 = [](IntegrationPoint ip)
    {
      double x = ip(0), y = ip(1);
      return IntegrationPoint(1.0/3+2.0/3*x-1.0/3*y, 1.0/3-1.0/3*x+2.0/3*y, 0, ip.Weight()/3);
    };

    
    for (auto ip : irtrig)
      {
        ir += map1(ip);
        ir += map2(ip);
        ir += map3(ip);
      }
    
    std::map<ELEMENT_TYPE, IntegrationRule> rules;
    rules[ET_TRIG] = std::move(ir);
    return rules;
  }
  
  void JKM_FESpace::Update()
  {
    FESpace::Update();
    first_edge_dof.SetSize(ma->GetNEdges()+1);
    first_element_dof.SetSize(ma->GetNE()+1);

    size_t ndof = 0;
    for (size_t i = 0; i < ma->GetNEdges(); i++)
      {
        first_edge_dof[i] = ndof;
        ndof += 2;
      }
    first_edge_dof[ma->GetNEdges()] = ndof;

    for (size_t i = 0; i < ma->GetNE(); i++)  
      {
        first_element_dof[i] = ndof;
        ndof += 9;
      }
    first_element_dof[ma->GetNE()] = ndof;
    
    SetNDof (ndof);
  }

  void JKM_FESpace::FinalizeUpdate()
  {
    FESpace::FinalizeUpdate();
  }

  void JKM_FESpace::GetDofNrs(ElementId ei, Array<int> &dnums) const
  {
    dnums.SetSize0(); 
    Ngs_Element ngel = ma->GetElement(ei);
    
    for (auto e : ngel.Edges())
      dnums += IntRange(first_edge_dof[e], first_edge_dof[e+1]);
    
    if (ei.VB() == VOL)
      dnums += IntRange(first_element_dof[ei.Nr()],
                        first_element_dof[ei.Nr()+1]);
  }

  FiniteElement &JKM_FESpace::GetFE(ElementId ei, Allocator &alloc) const
  {
    Ngs_Element ngel = ma->GetElement(ei);
    
    switch (ma->GetElType(ei))
      {
      case ET_TRIG:
        {
          auto el = new (alloc) JKMFE_Triangle(order);
          el -> SetVertexNumbers (ngel.Vertices());
          return *el;
        }
      default:
        throw Exception("JKMFESpace::GetFE not implemented yet");
      }
  }
  
}

