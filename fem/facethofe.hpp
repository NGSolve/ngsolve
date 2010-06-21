#ifndef FILE_FACETHOFE
#define FILE_FACETHOFE

/*********************************************************************/
/* File:   facethofe.hpp                                             */
/* Author: A. Sinwel, H. Egger, J. Schoeberl                         */
/* Date:   2008                                                      */
/*********************************************************************/


#include "l2hofe.hpp"
#include "tscalarfe.cpp"



namespace ngfem
{
  
  template <ELEMENT_TYPE ET> class FacetFE;


  template <ELEMENT_TYPE ET> 
  class FacetFiniteElement_Family :
    public FacetVolumeFiniteElement<ET_trait<ET>::DIM>, 
    public T_ScalarFiniteElement2<FacetFE<ET>, ET>,
    public ET_trait<ET> 
  {
  public:
    FacetFiniteElement_Family ()
      : FacetVolumeFiniteElement<ET_trait<ET>::DIM> (ET_TET) 
    { ; }
  };


  //------------------------------------------------------------
  template <>
  class FacetFE<ET_TRIG> : public FacetFiniteElement_Family<ET_TRIG>
  {
  protected:
    L2HighOrderFE<ET_SEGM> facets2[3];
  public:
    FacetFE(); 
    virtual void ComputeNDof();

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[2], TFA & shape) const
    {
      Tx lam[3] = { x[0], x[1], 1-x[0]-x[1] };
      
      INT<2> e = GetEdgeSort (facetnr, vnums);
      SetZero (shape, 0, ndof);

      int p = facet_order[facetnr];
      int ii = first_facet_dof[facetnr];
      LegendrePolynomial::Eval (p, lam[e[1]]-lam[e[0]], shape.Addr(ii));
    }
  };



  // --------------------------------------------------------
  template <>
  class FacetFE<ET_QUAD> : public FacetFiniteElement_Family<ET_QUAD>
  {
  protected:
    L2HighOrderFE<ET_SEGM> facets2[4];
  public:
    FacetFE();
    virtual void ComputeNDof();

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[2], TFA & shape) const
    {
      Tx x = hx[0], y = hx[1];
      Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
      
      INT<2> e = GetEdgeSort (facetnr, vnums);
      SetZero (shape, 0, ndof);

      int p = facet_order[facetnr];
      int ii = first_facet_dof[facetnr];
      LegendrePolynomial::Eval (p, sigma[e[1]]-sigma[e[0]], shape.Addr(ii));
    }
  };



  // --------------------------------------------------------
  template <>
  class FacetFE<ET_TET> : public FacetFiniteElement_Family<ET_TET>
  {
  protected:
    L2HighOrderFE<ET_TRIG> facets2[4];
  public:
    FacetFE(); 
    virtual void ComputeNDof();

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const
    {
      Tx lam[4] = { hx[0], hx[1], hx[2], 1-hx[0]-hx[1]-hx[2] };

      INT<4> f = GetFaceSort (facetnr, vnums);

      SetZero (shape, 0, ndof);
      
      int p = facet_order[facetnr];
      int ii = first_facet_dof[facetnr];
      DubinerBasis::Eval (p, lam[f[0]], lam[f[1]], shape.Addr(ii));
    }
  };



  // --------------------------------------------------------
  template <>
  class FacetFE<ET_HEX> : public FacetFiniteElement_Family<ET_HEX>
  {
  protected:
    L2HighOrderFE<ET_QUAD> facetsq[6];
  public:
    FacetFE();
    virtual void ComputeNDof();

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[3], TFA & shape) const
    {
      ;
    }

  };

  // --------------------------------------------------------
  template <>
  class FacetFE<ET_PRISM> : public FacetFiniteElement_Family<ET_PRISM>
  {
  protected:
    L2HighOrderFE<ET_TRIG> facetst[2];
    L2HighOrderFE<ET_QUAD> facetsq[3];

  public:
    FacetFE();
    virtual void ComputeNDof();

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const
    {
      Tx x = hx[0], y = hx[1], z = hx[2];
      Tx lam[6] = { x, y, 1-x-y, x, y, 1-x-y };
      Tx muz[6]  = { 1-z, 1-z, 1-z, z, z, z };
      

      INT<4> f = GetFaceSort (facetnr, vnums);

      SetZero (shape, 0, ndof);
      
      int p = facet_order[facetnr];
      int ii = first_facet_dof[facetnr];

      if (facetnr < 2)
	DubinerBasis::Eval (p, lam[f[0]], lam[f[1]], shape.Addr(ii));
      else
	{
	  Tx xi  = lam[f[0]]+muz[f[0]] - lam[f[1]]-muz[f[1]];
	  Tx eta = lam[f[0]]+muz[f[0]] - lam[f[3]]-muz[f[3]];

	  ArrayMem<Tx,20> polx(p+1), poly(p+1);

	  LegendrePolynomial::Eval (p, xi, polx);
	  LegendrePolynomial::Eval (p, eta, poly);

	  for (int i = 0; i <= p; i++)
	    for (int j = 0; j <= p; j++)
	      shape[ii++] = polx[i] * poly[j];
	}
    }

  };


  // --------------------------------------------------------
  template <>
  class FacetFE<ET_PYRAMID> : public FacetFiniteElement_Family<ET_PYRAMID>
  {
  protected:
    int tnr, qnr; // active facets
    L2HighOrderFE<ET_TRIG> trig;
    L2HighOrderFE<ET_QUAD> quad;
  public:
    FacetFE() : FacetFiniteElement_Family<ET_PYRAMID> () {  qnr=tnr=-1; };
    
    virtual void CalcFacetShape (int afnr, const IntegrationPoint & ip, FlatVector<> shape) const; 
    virtual void SetFacet(int afnr) const;
  
    virtual void ComputeNDof();
    //     virtual const FiniteElement & GetFacetFE(int fnr, LocalHeap& lh) const;


    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const
    {
      Tx x = hx[0], y = hx[1], z = hx[2];
      
      if (z == 1.) z -= 1e-10;
      
      Tx xt = x / (1-z);
      Tx yt = y / (1-z);
      
      Tx sigma[4]  = { (1-xt)+(1-yt), xt+(1-yt), xt+yt, (1-xt)+yt };
      Tx lambda[4] = { (1-xt)*(1-yt), xt*(1-yt), xt*yt, (1-xt)*yt };
      Tx lam[5];
      
      for (int i = 0; i < 4; i++)  
	lam[i] = lambda[i] * (1-z);
      lam[4] = z;


      INT<4> f = GetFaceSort (facetnr, vnums);

      SetZero (shape, 0, ndof);
      
      int p = facet_order[facetnr];
      int ii = first_facet_dof[facetnr];

      if (facetnr < 4)
	DubinerBasis::Eval (p, lam[f[0]], lam[f[1]], shape.Addr(ii));
      else
	{
	  Tx xi  = sigma[f[0]]-sigma[f[1]];
	  Tx eta = sigma[f[0]]-sigma[f[3]];

	  ArrayMem<Tx,20> polx(p+1), poly(p+1);

	  LegendrePolynomial::Eval (p, xi, polx);
	  LegendrePolynomial::Eval (p, eta, poly);

	  for (int i = 0; i <= p; i++)
	    for (int j = 0; j <= p; j++)
	      shape[ii++] = polx[i] * poly[j];
	}
    }

  };

}

#endif
