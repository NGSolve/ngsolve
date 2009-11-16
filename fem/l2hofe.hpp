#ifndef FILE_L2HOFE
#define FILE_L2HOFE

/*********************************************************************/
/* File:   l2hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include "tscalarfe.hpp"


namespace ngfem
{

  /**
     Base class for L2 - high order finite elements, i.e., a
     discontinuous approximation
  */
  template<int D>
  class L2HighOrderFiniteElement : virtual public ScalarFiniteElement<D>
  {
  protected:

    enum { DIM = D };
    int vnums[8];  
    INT<DIM> order_inner; 

  public:
    /// global vertex numbers define ordering of vertices
    void SetVertexNumbers (FlatArray<int> & avnums);

    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }

    /// set polynomial order
    void SetOrder (int o) 
    { for (int i = 0; i < DIM; i++) order_inner[i] = o; }

    /// different orders in differnt directions
    void SetOrder (INT<DIM> oi) { order_inner = oi; }

    /// calculate number of dofs
    virtual void ComputeNDof () = 0; 

  
    virtual void GetInternalDofs (Array<int> & idofs) const; 
  };



  /**
     Template family of L2 - high order finite elements.
     The template argument is the element shape
  */
  template <ELEMENT_TYPE ET> class L2HighOrderFE;


  /**
     Barton-Nackman base class for L2 - high order finite elements
  */
  template <ELEMENT_TYPE ET>
  class T_L2HighOrderFiniteElement : 
    public L2HighOrderFiniteElement<ET_trait<ET>::DIM>,
    public T_ScalarFiniteElement2< L2HighOrderFE<ET>, ET >,
    public ET_trait<ET> 
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using ScalarFiniteElement<DIM>::eltype;
    using ScalarFiniteElement<DIM>::dimspace;

    using L2HighOrderFiniteElement<DIM>::vnums;
    using L2HighOrderFiniteElement<DIM>::order_inner;

    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_EDGE;
    using ET_trait<ET>::N_FACE;
    using ET_trait<ET>::FaceType;
    using ET_trait<ET>::GetEdgeSort;
    using ET_trait<ET>::GetFaceSort;



  public:

    T_L2HighOrderFiniteElement () 
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
	vnums[i] = i;
      dimspace = DIM;
      eltype = ET;
    }

    virtual void ComputeNDof();

    /*
      virtual void CalcShape (const IntegrationPoint & ip, 
      FlatVector<> shape) const;

      virtual void CalcDShape (const IntegrationPoint & ip, 
      FlatMatrixFixWidth<DIM> dshape) const;
    */
  };






  /**
     L2 high order 1D finite element
  */
  template <>
  class L2HighOrderFE<ET_SEGM> : public T_L2HighOrderFiniteElement<ET_SEGM>
  {
  public:
    L2HighOrderFE () { ; }
    L2HighOrderFE (int aorder);

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[1], TFA & shape) const
    {
      Tx x = hx[0];
      // orient
      if ( vnums[0] > vnums[1])
	x = 1-x;
    
      LegendrePolynomial (order, 2*x-1, shape);
    }
  };


  /**
     L2 high order triangular finite element
  */
  template <> 
  class L2HighOrderFE<ET_TRIG> : public T_L2HighOrderFiniteElement<ET_TRIG>
  {
  public:
    L2HighOrderFE () { ; }
    L2HighOrderFE (int aorder);

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[2], TFA & shape) const
    {
      Tx x = hx[0];
      Tx y = hx[1];
      
      Tx lami[3] = { x, y, 1-x-y };
      int fav[3] = { 0, 1, 2};
      if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 
      if(vnums[fav[1]] > vnums[fav[2]]) swap(fav[1],fav[2]);
      if(vnums[fav[0]] > vnums[fav[1]]) swap(fav[0],fav[1]); 	
      x = lami[fav[0]]; y=lami[fav[1]];
      
      int n = order_inner[0];
      ArrayMem<Tx, 20> polx(n+1), poly(n+1);

      ScaledLegendrePolynomial (n, 2*x+y-1, 1-y, polx);
    
      for (int i = 0, ii = 0; i <= order_inner[0]; i++)
	{
	  JacobiPolynomial (n-i, 2*y-1, 2*i+1, 0, poly);
	  for (int j = 0; j <= order_inner[0]-i; j++)
	    shape[ii++] = polx[i] * poly[j];
	}
    }
  };


  /**
     L2 high order quadrilateral finite element
  */
  template <> 
  class L2HighOrderFE<ET_QUAD> : public T_L2HighOrderFiniteElement<ET_QUAD>
  {
  public:
    L2HighOrderFE () { ; }
    L2HighOrderFE (int aorder);

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[2], TFA & shape) const
    {
      Tx x = hx[0];
      Tx y = hx[1];
      
      // orient: copied from h1
      // Tx lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
      
      Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
      int fmax = 0; 
      for (int j = 1; j < 4; j++)
	if (vnums[j] > vnums[fmax]) fmax = j;
      int f1 = (fmax+3)%4; 
      int f2 = (fmax+1)%4; 
      if(vnums[f2] > vnums[f1]) swap(f1,f2);  // fmax > f1 > f2; 
      x = sigma[fmax]-sigma[f1]; 
      y = sigma[fmax]-sigma[f2]; 
      
      int ii = 0;
      
      int n = max(order_inner[0],order_inner[1]);
      ArrayMem<Tx, 20> polx(n+1), poly(n+1);
      
      LegendrePolynomial (n, x, polx);
      LegendrePolynomial (n, y, poly);
      
      for (int i = 0; i <= order_inner[0]; i++)
	for (int j = 0; j <= order_inner[1]; j++)
	  shape[ii++] = polx[i] * poly[j];
    }
  };


  /**
     L2 high order tetrahedral finite element
  */
  template <> 
  class L2HighOrderFE<ET_TET> : public T_L2HighOrderFiniteElement<ET_TET>
  {
  public:
    L2HighOrderFE () { ; }
    L2HighOrderFE (int aorder);

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];
    
      // no orientation necessary
      int n=order;
      ArrayMem<Tx, 20> polx(n+1), poly(n+1), polz(n+1);
    
      // Polynomials orthogonal w.r.t L2 inner product
      ScaledLegendrePolynomial ( n, 2*x+y+z-1, 1-y-z, polx);
  
      for (int i = 0, ii = 0; i <= order; i++)
	{
	  ScaledJacobiPolynomial (n, 2*y+z-1, 1-z, 2*i+1, 0, poly);
	  for (int j = 0; j <= order-i; j++)
	    {
	      JacobiPolynomial(n, 2*z-1, 2*i+2*j+2, 0, polz);
	      for (int k = 0; k <= order-i-j; k++)
		shape[ii++] = polx[i] * poly[j] * polz[k];
	    }
	}
    }
  };


  /**
     L2 high order prismatic finite element
  */
  template <> 
  class L2HighOrderFE<ET_PRISM> : public T_L2HighOrderFiniteElement<ET_PRISM>
  {
  public:
    L2HighOrderFE () { ; }
    L2HighOrderFE (int aorder);

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];
    
      // no orientation necessary
      int p=order_inner[0];
      int q=order_inner[1];
    
      ArrayMem<Tx, 20> polx(p+1), poly(p+1), polz(q+1);

      LegendrePolynomial (p, 2*x-1, polx);
      ScaledLegendrePolynomial (p, 2*y+x-1, 1-x, poly);
    
      LegendrePolynomial (q, 2*z-1, polz);
  
      int ii = 0;
      for (int i = 0; i <= p; i++)
	for (int j = 0; j <= p-i; j++)
	  for (int k = 0; k <= q; k++)
	    shape[ii++] = polx[i] * poly[j] * polz[k];
    }
  };


  /**
     L2 high order pyramid finite element
  */
  template <> 
  class L2HighOrderFE<ET_PYRAMID> : public T_L2HighOrderFiniteElement<ET_PYRAMID>
  {
  public:
    L2HighOrderFE () { ; }
    L2HighOrderFE (int aorder);

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const
    {
      Tx x = hx[0];
      Tx y = hx[1];
      Tx z = hx[2];

      if (z == 1) z -= 1e-8;
      Tx xt = 2 * (x / (1-z)) - 1;
      Tx yt = 2 * (y / (1-z)) - 1;

      VectorMem<10, Tx> polsx(order+1);
      VectorMem<10, Tx> polsy(order+1);
    
      ArrayMem<Tx, 20> memz(sqr(order+1));
      FlatMatrix<Tx> polsz(order+1, &memz[0]);

      Tx fac = 1.0;
      for (int i = 0; i <= order; i++)
	{
	  JacobiPolynomial (order, 2*z-1, 2*i+2, 0, polsz.Row(i));
	  polsz.Row(i) *= fac;
	  fac *= (1-z);
	}

      LegendrePolynomial (order, xt, polsx);
      LegendrePolynomial (order, yt, polsy);
    
      int ii = 0;
      for (int iz = 0; iz <= order; iz++)
	for (int ix = 0; ix <= order-iz; ix++)
	  for (int iy = 0; iy <= order-iz; iy++, ii++)
	    shape[ii] = polsx(ix) * polsy(iy) * polsz(max(ix,iy), iz);
    }

  };


  /**
     L2 high order hexahedral finite element
  */
  template <> 
  class L2HighOrderFE<ET_HEX> : public T_L2HighOrderFiniteElement<ET_HEX>
  {
  public:
    L2HighOrderFE () { ; }
    L2HighOrderFE (int aorder);

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const
    {
      Tx x = hx[0], y = hx[1], z = hx[2];
    
      // no orientation necessary
      int p=order_inner[0];
      int q=order_inner[1];
      int r=order_inner[2];
    
      ArrayMem<Tx, 20> polx(p+1), poly(q+1), polz(r+1);

      LegendrePolynomial (p, 2*x-1, polx);
      LegendrePolynomial (q, 2*y-1, poly);
      LegendrePolynomial (r, 2*z-1, polz);
  
      for (int i = 0, ii = 0; i <= p; i++)
	for (int j = 0; j <= q; j++)
	  for (int k = 0; k <= r; k++)
	    shape[ii++] = polx[i] * poly[j] * polz[k];
    }

  };

}

#endif
