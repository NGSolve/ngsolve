#ifndef FILE_H1HOFE
#define FILE_H1HOFE

/*********************************************************************/
/* File:   h1hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include "tscalarfe.hpp"


namespace ngfem
{

  /**
     High order finite elements for H^1
  */
  template<int DIM>
  class H1HighOrderFiniteElement : virtual public ScalarFiniteElement<DIM>
  {
  public:
    int vnums[8];
    INT<3> order_cell;
    INT<2> order_face[6];
    int order_edge[12];

    using ScalarFiniteElement<DIM>::eltype;

  public:
    void SetVertexNumbers (const FlatArray<int> & avnums)
    {
      for (int i = 0; i < avnums.Size(); i++)
        vnums[i] = avnums[i];
    }

    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }

    void SetOrderCell (int oi)   { order_cell = INT<3> (oi,oi,oi); }
    void SetOrderCell (INT<3> oi)  { order_cell = oi; }

    void SetOrderFace (FlatArray<int> & of);
    void SetOrderFace (FlatArray<INT<2> > & of);
    void SetOrderFace (int nr, INT<2> order) { order_face[nr] = order; }

    void SetOrderEdge (FlatArray<int> & oe);
    void SetOrderEdge (int nr, int order) { order_edge[nr] = order; }


    /// high order elements need extra configuration. update ndof and order
    virtual void ComputeNDof () = 0;
  };



  template <ELEMENT_TYPE ET> class H1HighOrderFE;


  /**
     Barton-Nackman base class for H1 - high order finite elements
  */
  template <ELEMENT_TYPE ET>
  class T_H1HighOrderFiniteElement : 
    public H1HighOrderFiniteElement<ET_trait<ET>::DIM>, 
    public T_ScalarFiniteElement2< H1HighOrderFE<ET>, ET >,
    public ET_trait<ET> 
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using ScalarFiniteElement<DIM>::eltype;
    using ScalarFiniteElement<DIM>::dimspace;

    using H1HighOrderFiniteElement<DIM>::vnums;
    using H1HighOrderFiniteElement<DIM>::order_edge;
    using H1HighOrderFiniteElement<DIM>::order_face;
    using H1HighOrderFiniteElement<DIM>::order_cell;


    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_EDGE;
    using ET_trait<ET>::N_FACE;
    using ET_trait<ET>::FaceType;
    using ET_trait<ET>::GetEdgeSort;
    using ET_trait<ET>::GetFaceSort;

    typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

    typedef TrigShapesInnerLegendre T_TRIGSHAPES;
    // typedef TrigShapesInnerJacobi T_TRIGSHAPES;

  public:

    T_H1HighOrderFiniteElement () 
    {
      for (int i = 0; i < N_VERTEX; i++)
	vnums[i] = i;
      eltype = ET;
    }

    T_H1HighOrderFiniteElement (int aorder) 
    {
      eltype = ET;

      for (int i = 0; i < N_VERTEX; i++)
	vnums[i] = i;

      for (int i = 0; i < N_EDGE; i++)
        order_edge[i] = aorder;
      for (int i = 0; i < N_FACE; i++)
        order_face[i] = INT<2> (aorder,aorder);
      if (DIM == 3)
        order_cell = INT<3> (aorder,aorder,aorder);
      
      order = aorder;
    }


    virtual void GetDofs (Array<Dof> & dofs) const;

    virtual void ComputeNDof();
    virtual void GetInternalDofs (Array<int> & idofs) const;

    /*
      virtual void CalcShape (const IntegrationPoint & ip, 
      FlatVector<> shape) const;

      virtual double Evaluate (const IntegrationPoint & ip, FlatVector<double> x) const;

      virtual void CalcDShape (const IntegrationPoint & ip, 
      FlatMatrixFixWidth<DIM> dshape) const;

      virtual void CalcMappedDShape (const SpecificIntegrationPoint<DIM,DIM> & sip, 
      FlatMatrixFixWidth<DIM> dshape) const;
    */
  };







  /**
     High order segment finite element
  */

  template <> 
  class H1HighOrderFE<ET_SEGM> : public T_H1HighOrderFiniteElement<ET_SEGM>
  {
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement<ET_SEGM> (aorder) 
    { ndof = (order+1); }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[1], TFA & shape) const
    {
      Tx x = hx[0];
      Tx lami[2] = { x, 1-x };
      
      shape[0] = lami[0];
      shape[1] = lami[1];

      INT<2> e = GetEdgeSort (0, vnums);
      T_ORTHOPOL::Calc (order_edge[0], lami[e[1]]-lami[e[0]], shape.Addr(2));
    }
    
  };


  /**
     High order triangular finite element
  */
  template <>
  class H1HighOrderFE<ET_TRIG> : public T_H1HighOrderFiniteElement<ET_TRIG>
  {
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement<ET_TRIG> (aorder) 
    { ndof = (order+1)*(order+2)/2; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[2], TFA & shape) const
    {
      Tx lam[3] = { x[0], x[1], 1-x[0]-x[1] };

      for (int i = 0; i < 3; i++)
	shape[i] = lam[i];

      int ii = 3;
    
      // edge dofs
      for (int i = 0; i < N_EDGE; i++)
	if (order_edge[i] >= 2)
	  { 
	    INT<2> e = GetEdgeSort (i, vnums);

	    ii += T_ORTHOPOL::CalcTrigExt (order_edge[i], 
					   lam[e[1]]-lam[e[0]], 1-lam[e[0]]-lam[e[1]], 
					   shape.Addr(ii));
	  }

      int p = order_face[0][0];
      if (p >= 3)
	{
	  INT<4> f = GetFaceSort (0, vnums);

	  ArrayMem<Tx, 20> polx(p-2), poly(p-2);
	  T_TRIGSHAPES::CalcSplitted (p, lam[f[2]]-lam[f[1]],
				      lam[f[0]], polx, poly);
	  for (int i = 0; i <= p-3; i++)
	    for (int j = 0; j <= p-3-i; j++)
	      shape[ii++] = polx[i] * poly[j];
	}
    }
  };


  /**
     High order quadrilateral finite element
  */
  template <>
  class H1HighOrderFE<ET_QUAD> : public T_H1HighOrderFiniteElement<ET_QUAD>
  {
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement<ET_QUAD> (aorder) 
    { ndof = (order+1)*(order+1); }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[2], TFA & shape) const
    {
      Tx x = hx[0], y = hx[1];
      Tx lami[4] = {(1-x)*(1-y),x*(1-y),x*y,(1-x)*y};  
      Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
    

      // vertex shapes
      for(int i=0; i < N_VERTEX; i++) shape[i] = lami[i]; 
      int ii = 4;

      ArrayMem<Tx,20> polxi(order+1), poleta(order+1);
     
      // edge dofs
      for (int i = 0; i < N_EDGE; i++)
	{
	  int p = order_edge[i];
	  INT<2> e = GetEdgeSort (i, vnums);	  
        
	  Tx xi = sigma[e[1]]-sigma[e[0]]; 
	  Tx lam_e = lami[e[0]]+lami[e[1]];
        
	  ii += T_ORTHOPOL::CalcMult (p, xi, lam_e, shape.Addr(ii));
	}    
    
      INT<2> p = order_face[0];
      if (p[0] >= 2 && p[1] >= 2)
	{
	  INT<4> f = GetFaceSort (0, vnums);  // vnums[f[0]] > vnums[f[1]] > vnums[f[3]]

	  Tx xi = sigma[f[0]]-sigma[f[1]]; 
	  Tx eta = sigma[f[0]]-sigma[f[3]]; 
	
	  T_ORTHOPOL::Calc(p[0], xi,polxi);
	  T_ORTHOPOL::Calc(p[1],eta,poleta);
	
	  for (int k = 0; k <= p[0]-2; k++)
	    for (int j = 0; j <= p[1]-2; j++)
	      shape[ii++] = polxi[k] * poleta[j];
	}
    }
      
  };


  /**
     High order tetrahedral finite element
  */
  template <>
  class H1HighOrderFE<ET_TET> : public T_H1HighOrderFiniteElement<ET_TET>
  {
    typedef TetShapesInnerLegendre T_INNERSHAPES;
    typedef TetShapesFaceLegendre T_FACESHAPES;

    // typedef TetShapesInnerJacobi T_INNERSHAPES;
    // typedef TetShapesFaceJacobi T_FACESHAPES;
    // typedef TetShapesFaceOpt1 T_FACESHAPES;

  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement<ET_TET> (aorder) 
    { ndof = (order+1)*(order+2)*(order+3)/6; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const; 
  };


  /** 
      High order prismatic finite element
  */
  template <>
  class H1HighOrderFE<ET_PRISM> : public T_H1HighOrderFiniteElement<ET_PRISM>
  {
    // typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;
    // typedef TrigShapesInnerJacobi T_TRIGFACESHAPES;
    // typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement<ET_PRISM> (aorder) 
    { ndof = (order+1)*(order+2)*(order+1)/2; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const; 
  };



  /**
     High order hexahedral finite element
  */
  template <> 
  class H1HighOrderFE<ET_HEX> : public T_H1HighOrderFiniteElement<ET_HEX>
  {
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement<ET_HEX> (aorder) 
    { ndof = (order+1)*(order+1)*(order+1); }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const; 
  };


  /**
     High order pyramid finite element
  */
  template<>
  class H1HighOrderFE<ET_PYRAMID> : public T_H1HighOrderFiniteElement<ET_PYRAMID>
  {
    // typedef TrigShapesInnerLegendre T_TRIGSHAPES;

  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement<ET_PYRAMID> (aorder) 
    { ndof = (order+2)*(order+1)*(2*order+3) / 6; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const; 
  };

}


#endif
