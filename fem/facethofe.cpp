/*********************************************************************/
/* File:   facethofe.hpp                                             */
/* Author: A. Sinwel, H. Egger, J. Schoeberl                         */
/* Date:   2008                                                      */
/*********************************************************************/

#define FILE_FACETHOFE_CPP

 
#include <fem.hpp>
#include <tscalarfe_impl.hpp>
#include <facethofe.hpp>






namespace ngfem
{ 

  
  template <> template<typename Tx, typename TFA>  
  void FacetFE<ET_SEGM> :: T_CalcShapeFNr (int fnr, Tx x[2], TFA & shape) const
  {
    shape[0] = 1.0;
  }

  
  template <> template<typename Tx, typename TFA>  
  void FacetFE<ET_TRIG> :: T_CalcShapeFNr (int fnr, Tx x[2], TFA & shape) const
  {
    Tx lam[3] = { x[0], x[1], 1-x[0]-x[1] };
    
    INT<2> e = GetVertexOrientedEdge (fnr);
    int p = facet_order[fnr];
    
    LegendrePolynomial::Eval (p, lam[e[1]]-lam[e[0]], shape);
  }



  // --------------------------------------------------------

  template<> template<typename Tx, typename TFA>  
  void FacetFE<ET_QUAD> :: T_CalcShapeFNr (int fnr, Tx hx[2], TFA & shape) const
    {
      Tx x = hx[0], y = hx[1];
      Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
      
      INT<2> e = GetVertexOrientedEdge (fnr);
      int p = facet_order[fnr];

      LegendrePolynomial::Eval (p, sigma[e[1]]-sigma[e[0]], shape);
    }


  template <> template<typename Tx, typename TFA>  
  void FacetFE<ET_TET> :: T_CalcShapeFNr (int fnr, Tx hx[3], TFA & shape) const
  {
    Tx lam[4] = { hx[0], hx[1], hx[2], 1-hx[0]-hx[1]-hx[2] };
    
    INT<4> f = GetVertexOrientedFace (fnr);
    int p = facet_order[fnr];
    
    DubinerBasis3::Eval (p, lam[f[0]], lam[f[1]], shape);
  }






  // --------------------------------------------------------
  
  template <> template<typename Tx, typename TFA>  
  void FacetFE<ET_HEX> :: T_CalcShapeFNr (int fnr, Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];
    Tx sigma[8]={(1-x)+(1-y)+(1-z),x+(1-y)+(1-z),x+y+(1-z),(1-x)+y+(1-z),
		 (1-x)+(1-y)+z,x+(1-y)+z,x+y+z,(1-x)+y+z}; 
    
    int p = facet_order[fnr];
    
    INT<4> f = GetVertexOrientedFace (fnr);	  
	    
    Tx xi  = sigma[f[0]] - sigma[f[1]]; 
    Tx eta = sigma[f[0]] - sigma[f[3]];
    
    ArrayMem<Tx,20> polx(p+1), poly(p+1);
    
    LegendrePolynomial::Eval(p, xi, polx);
    LegendrePolynomial::Eval(p, eta, poly);
    
    for (int i = 0, ii = 0; i <= p; i++)
      for (int j = 0; j <= p; j++)
	shape[ii++] = polx[i] * poly[j];
  }


  // --------------------------------------------------------
  

  template <> template<typename Tx, typename TFA>  
  void FacetFE<ET_PRISM> :: T_CalcShapeFNr (int fnr, Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];
    Tx lam[6] = { x, y, 1-x-y, x, y, 1-x-y };
    Tx muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
    
    
    INT<4> f = GetVertexOrientedFace (fnr);
    
    int p = facet_order[fnr];
    
    if (fnr < 2)
      DubinerBasis3::Eval (p, lam[f[0]], lam[f[1]], shape);
    else
      {
	Tx xi  = lam[f[0]]+muz[f[0]] - lam[f[1]]-muz[f[1]];
	Tx eta = lam[f[0]]+muz[f[0]] - lam[f[3]]-muz[f[3]];
	
	ArrayMem<Tx,20> polx(p+1), poly(p+1);
	
	LegendrePolynomial::Eval (p, xi, polx);
	LegendrePolynomial::Eval (p, eta, poly);

	for (int i = 0, ii = 0; i <= p; i++)
	  for (int j = 0; j <= p; j++)
	    shape[ii++] = polx[i] * poly[j];
      }
  }
  

  // --------------------------------------------------------
  

  template <> template<typename Tx, typename TFA>  
  void FacetFE<ET_PYRAMID> :: T_CalcShapeFNr (int fnr, Tx hx[3], TFA & shape) const
  {
    Tx x = hx[0], y = hx[1], z = hx[2];
      
    // if (z == 1.) z -= 1e-10;
    z *= (1-1e-10);
    
    Tx xt = x / (1-z);
    Tx yt = y / (1-z);
    
    Tx sigma[4]  = { (1-xt)+(1-yt), xt+(1-yt), xt+yt, (1-xt)+yt };
    Tx lambda[4] = { (1-xt)*(1-yt), xt*(1-yt), xt*yt, (1-xt)*yt };
    Tx lam[5];
      
    for (int i = 0; i < 4; i++)  
      lam[i] = lambda[i] * (1-z);
    lam[4] = z;
    
    
    INT<4> f = GetVertexOrientedFace (fnr);
    
    int p = facet_order[fnr];
    
    if (fnr < 4)
      DubinerBasis3::Eval (p, lam[f[0]], lam[f[1]], shape);
    else
      {
	Tx xi  = sigma[f[0]]-sigma[f[1]];
	Tx eta = sigma[f[0]]-sigma[f[3]];
	
	ArrayMem<Tx,20> polx(p+1), poly(p+1);
	
	LegendrePolynomial::Eval (p, xi, polx);
	LegendrePolynomial::Eval (p, eta, poly);
	
	for (int i = 0, ii = 0; i <= p; i++)
	  for (int j = 0; j <= p; j++)
	    shape[ii++] = polx[i] * poly[j];
      }
  }


  template class FacetFE<ET_SEGM>;
  template class FacetFE<ET_TRIG>;
  template class FacetFE<ET_QUAD>;
  template class FacetFE<ET_TET>;
  template class FacetFE<ET_HEX>;
  template class FacetFE<ET_PRISM>;
  template class FacetFE<ET_PYRAMID>;
  template class FacetFiniteElement_Family<ET_SEGM>;
  template class FacetFiniteElement_Family<ET_TRIG>;
  template class FacetFiniteElement_Family<ET_QUAD>;
  template class FacetFiniteElement_Family<ET_TET>;
  template class FacetFiniteElement_Family<ET_HEX>;
  template class FacetFiniteElement_Family<ET_PRISM>;
  template class FacetFiniteElement_Family<ET_PYRAMID>;
}
 
