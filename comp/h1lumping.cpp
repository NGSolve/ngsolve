// #include <comp.hpp>    // provides FESpace, ...

#include "h1lumping.hpp"
#include "../fem/tscalarfe_impl.hpp"
#include "../fem/h1lofe.hpp"
#include "bdbequations.hpp"


/*
HIGHER ORDER TRIANGULAR FINITE ELEMENTS WITH MASS LUMPING FOR THE WAVE EQUATION
G. COHEN, P. JOLY, J. E. ROBERTS, AND N. TORDJMAN
SIAM J. NUMER. ANAL. ⃝c 2001 Society for Industrial and Applied Mathematics Vol. 38, No. 6, pp. 2047–2078


3D elements:
NEW HIGHER-ORDER MASS-LUMPED TETRAHEDRAL ELEMENTS FOR WAVE PROPAGATION MODELLING*
S. GEEVERS1, W.A. MULDER2,3 AND J.J.W. VAN DER VEGT1
https://arxiv.org/pdf/1803.10065.pdf
*/


namespace ngcomp
{


  class H1LumpingSegm2 : public T_ScalarFiniteElementFO<H1LumpingSegm2,ET_SEGM,3,2>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<1,Tx> ip, TFA & shape) 
    {
      Tx lam[] = { ip.x,  1-ip.x };
      for (int i = 0; i < 2; i++)
        shape[i] = 2*lam[i]*(lam[i]-0.5);
      shape[2] = 4 * lam[0] * lam[1];
    }
  };

  

  class H1LumpingTrig2 : public T_ScalarFiniteElementFO<H1LumpingTrig2,ET_TRIG,7,3>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<2,Tx> ip, TFA & shape) 
    {
      Tx lam[] = { ip.x, ip.y, 1-ip.x-ip.y };
      Tx bub = 27*lam[0]*lam[1]*lam[2];
      for (int i = 0; i < 3; i++)
        shape[i] = 2*lam[i]*(lam[i]-0.5) + 1.0/9 * bub;

      const EDGE * edges = ElementTopology::GetEdges (ET_TRIG);
      for (int i = 0; i < 3; i++)
        shape[i+3] = 4 * lam[edges[i][0]] * lam[edges[i][1]] - 4.0/9*bub;
      
      shape[6] = bub;      
    }
  };


  class H1LumpingTet2 : public T_ScalarFiniteElementFO<H1LumpingTet2,ET_TET,15,4>
  {
  public:
    template<typename Tx, typename TFA>  
    static INLINE void T_CalcShape (TIP<3,Tx> ip, TFA & shape) 
    {
      Tx lam[] = { ip.x, ip.y, ip.z, 1-ip.x-ip.y-ip.z };
      Tx bub = 256*lam[0]*lam[1]*lam[2]*lam[3];
      Tx bubf[] = { 27*lam[1]*lam[2]*lam[3] - 27./64*bub,
                    27*lam[0]*lam[2]*lam[3] - 27./64*bub,
                    27*lam[0]*lam[1]*lam[3] - 27./64*bub,
                    27*lam[0]*lam[1]*lam[2] - 27./64*bub };
      Tx sumbubf = bubf[0]+bubf[1]+bubf[2]+bubf[3];

      for (int i = 0; i < 4; i++)
        shape[i] = 2*lam[i]*(lam[i]-0.5)+1./8*bub + 1.0/9*(sumbubf-bubf[i]);

      const EDGE * edges = ElementTopology::GetEdges (ET_TET);
      for (int i = 0; i < 6; i++)
        shape[i+4] = 4 * lam[edges[i][0]] * lam[edges[i][1]] - 1.0/4 * bub
          - 4.0/9 * (sumbubf-bubf[edges[i][0]]-bubf[edges[i][1]]);
      
      for (int i = 0; i < 4; i++)
        shape[10+i] = bubf[i];
      shape[14] = bub;      
    }
  };



  
  H1LumpingFESpace :: H1LumpingFESpace (shared_ptr<MeshAccess> ama, const Flags & flags)
    : FESpace (ama, flags)
  {
    type = "h1lumpingfespace";

    if (ma->GetDimension() == 2)
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
      }
    else
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdH1<3,2>>>();        
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
      }
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
    nface = (ma->GetDimension() == 2) ? 0 : ma->GetNFaces();
    switch (order)
      {
      case 1:
        SetNDof (nvert); break;
      case 2:
        SetNDof (nvert+nedge+nface+ma->GetNE(VOL)); break;
      default:
        throw Exception("H1LumpingFESpace only supports order 1 or 2");
      }
  }

  void H1LumpingFESpace :: GetDofNrs (ElementId ei, Array<DofId> & dnums) const
  {
    dnums.SetSize(0);
    dnums += ma->GetElement(ei).Vertices();

    if (order == 1) return;
    
    for (auto e : ma->GetElement(ei).Edges())
      dnums.Append (nvert+e);

    if (ma->GetDimension()==3)
      for (auto f : ma->GetElement(ei).Faces())
        dnums.Append (nvert+nedge+f);

    if (ei.VB()==VOL)
      dnums.Append(nvert+nedge+nface+ei.Nr());
  }

  FiniteElement & H1LumpingFESpace :: GetFE (ElementId ei, Allocator & alloc) const
  {
    switch (order)
      {
      case 1:
        switch (ma->GetElement(ei).GetType())
          {
          case ET_SEGM:
            return * new (alloc) ScalarFE<ET_SEGM,1>;
          case ET_TRIG:
          return * new (alloc)  ScalarFE<ET_TRIG,1>;
          case ET_TET:
            return * new (alloc) ScalarFE<ET_TET,1>;
          default:
            throw Exception("H1Lumping: Element of type "+ToString(ma->GetElement(ei).GetType()) + 
                            " not available\n");
          }
        
      case 2:
        switch (ma->GetElement(ei).GetType())
          {
          case ET_SEGM:
            return * new (alloc) H1LumpingSegm2;
          case ET_TRIG:
          return * new (alloc) H1LumpingTrig2;
          case ET_TET:
            return * new (alloc) H1LumpingTet2;
          default:
            throw Exception("H1Lumping: Element of type "+ToString(ma->GetElement(ei).GetType()) + 
                            " not available\n");
          }
      }
    throw Exception("H1LumpingFESpace - undefined order or element");
  }

  std::map<ELEMENT_TYPE, IntegrationRule> H1LumpingFESpace :: GetIntegrationRules() const
  {
    std::map<ELEMENT_TYPE, IntegrationRule> rules;

    switch (order)
      {
      case 1:
        {
          IntegrationRule ir3;
          ir3.Append ( IntegrationPoint( 1, 0, 0, 1.0/6));
          ir3.Append ( IntegrationPoint( 0, 1, 0, 1.0/6));
          ir3.Append ( IntegrationPoint( 0, 0, 0, 1.0/6));        
          rules[ET_TRIG] = std::move(ir3);
        
          IntegrationRule ir4;  // tet
          ir4.Append ( IntegrationPoint( 1, 0, 0, 1.0/24) );
          ir4.Append ( IntegrationPoint( 0, 1, 0, 1.0/24) );
          ir4.Append ( IntegrationPoint( 0, 0, 1, 1.0/24) );
          ir4.Append ( IntegrationPoint( 0, 0, 0, 1.0/24) );
          rules[ET_TET] = std::move(ir4);
          break;
        }

      case 2:
        {
          IntegrationRule ir7;
          ir7.Append ( IntegrationPoint( 1, 0, 0, 1.0/40));
          ir7.Append ( IntegrationPoint( 0, 1, 0, 1.0/40));
          ir7.Append ( IntegrationPoint( 0, 0, 0, 1.0/40));
          ir7.Append ( IntegrationPoint( 0.5, 0, 0, 1.0/15));
          ir7.Append ( IntegrationPoint( 0, 0.5, 0, 1.0/15));
          ir7.Append ( IntegrationPoint( 0.5, 0.5, 0, 1.0/15));
          ir7.Append ( IntegrationPoint( 1.0/3, 1.0/3, 0, 9.0/40));
          
          rules[ET_TRIG] = std::move(ir7);
          
          
          IntegrationRule ir15;  // tet
          ir15.Append ( IntegrationPoint( 1, 0, 0, 17./5040) );
          ir15.Append ( IntegrationPoint( 0, 1, 0, 17./5040) );
          ir15.Append ( IntegrationPoint( 0, 0, 1, 17./5040) );
          ir15.Append ( IntegrationPoint( 0, 0, 0, 17./5040) );
          
          ir15.Append ( IntegrationPoint( 0.5, 0, 0,   2./315) );      
          ir15.Append ( IntegrationPoint( 0.5, 0.5, 0, 2./315) );     
          ir15.Append ( IntegrationPoint( 0.5, 0, 0.5, 2./315) );     
          ir15.Append ( IntegrationPoint( 0, 0.5, 0,   2./315) );      
          ir15.Append ( IntegrationPoint( 0, 0, 0.5,   2./315) );      
          ir15.Append ( IntegrationPoint( 0, 0.5, 0.5, 2./315) );
          
          ir15.Append ( IntegrationPoint( 1./3, 1./3, 1./3, 9./560) );
          ir15.Append ( IntegrationPoint( 0, 1./3, 1./3, 9./560) );
          ir15.Append ( IntegrationPoint( 1./3, 0, 1./3, 9./560) );
          ir15.Append ( IntegrationPoint( 1./3, 1./3, 0, 9./560) );    
          
          ir15.Append ( IntegrationPoint( 1./4, 1./4, 1./4, 16./315) );
          
          rules[ET_TET] = std::move(ir15);
          break;
        }
        
      default:
        throw Exception ("H1LumpingFESpace, undefined order");
      }
    
    return rules;
  }
  
}
