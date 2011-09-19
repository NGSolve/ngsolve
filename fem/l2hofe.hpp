#ifndef FILE_L2HOFE
#define FILE_L2HOFE

/*********************************************************************/
/* File:   l2hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/


#include "tscalarfe.hpp"
#include "precomp.hpp"


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

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using ScalarFiniteElement<DIM>::eltype;


  public:
    /// global vertex numbers define ordering of vertices
    template <typename TA>
    void SetVertexNumbers (const TA & avnums)
    { for (int i = 0; i < avnums.Size(); i++) vnums[i] = avnums[i]; }

    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }

    /// set polynomial order
    void SetOrder (int p)  { order_inner = p; }

    /// different orders in differnt directions
    void SetOrder (INT<DIM> p) { order_inner = p; }

    /// calculate number of dofs
    virtual void ComputeNDof () = 0; 
  
    virtual void GetInternalDofs (Array<int> & idofs) const; 

    virtual void PrecomputeTrace () = 0;

    void CalcTraceMatrix (int facet, FlatMatrix<> & trace) const;

    virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const
    {
      Matrix<> trace(fcoefs.Size(), coefs.Size());
      CalcTraceMatrix(facet, trace);
      fcoefs = trace * coefs;
    }

    virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const
    {
      Matrix<> trace(fcoefs.Size(), coefs.Size());
      CalcTraceMatrix(facet, trace);
      coefs = Trans (trace) * fcoefs;
    }
  };



  /**
     Template family of L2 - high order finite elements.
     The template argument is the element shape
  */
  // template <ELEMENT_TYPE ET> class L2HighOrderFE;


  /**
     Barton-Nackman base class for L2 - high order finite elements
  */
  template <ELEMENT_TYPE ET>
  class T_L2HighOrderFiniteElement : 
    public L2HighOrderFiniteElement<ET_trait<ET>::DIM>,
    public ET_trait<ET> 
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using ScalarFiniteElement<DIM>::eltype;

    using L2HighOrderFiniteElement<DIM>::vnums;
    using L2HighOrderFiniteElement<DIM>::order_inner;

    using ET_trait<ET>::N_VERTEX;
    using ET_trait<ET>::N_EDGE;
    using ET_trait<ET>::N_FACE;
    using ET_trait<ET>::FaceType;
    using ET_trait<ET>::GetEdgeSort;
    using ET_trait<ET>::GetFaceSort;
    using ET_trait<ET>::PolDimension;


  public:

    T_L2HighOrderFiniteElement () 
    { eltype = ET; }

    T_L2HighOrderFiniteElement (int aorder) 
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++) vnums[i] = i;
      eltype = ET;

      order = aorder;
      order_inner = aorder;
      ndof = PolDimension (aorder);
    }

    virtual void ComputeNDof();
  };






  template <ELEMENT_TYPE ET> class L2HighOrderFE_Shape;

  template <int DIM>
  class PrecomputedScalShapes
  {
  public:
    Matrix<> shapes;
    Matrix<> dshapes;
    
    PrecomputedScalShapes (int nip, int ndof)
      : shapes(nip, ndof), dshapes (DIM*nip, ndof)
    { ; }

  };


  template <ELEMENT_TYPE ET, template <ELEMENT_TYPE ET> class SHAPES = L2HighOrderFE_Shape> 
  class NGS_DLL_HEADER L2HighOrderFE :  public T_L2HighOrderFiniteElement<ET>,
					public T_ScalarFiniteElement2<SHAPES<ET>, ET >
  { 
  protected:
    using T_L2HighOrderFiniteElement<ET>::ndof;
    using T_L2HighOrderFiniteElement<ET>::order;
    using T_L2HighOrderFiniteElement<ET>::vnums;

    enum { DIM = ET_trait<ET>::DIM };
    typedef PrecomputedShapesContainer<PrecomputedScalShapes<DIM> > TPRECOMP;
    static TPRECOMP precomp;

    typedef HashTable<INT<2>, Matrix<>*> TPRECOMP_TRACE;
    static TPRECOMP_TRACE precomp_trace;
  public:
    L2HighOrderFE () { ; }

    L2HighOrderFE (int aorder)
      : T_L2HighOrderFiniteElement<ET> (aorder) { ; }


    virtual void PrecomputeTrace ();
    virtual void PrecomputeShapes (const IntegrationRule & ir);

    virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const;

    virtual void EvaluateGrad (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<DIM> values) const;

    virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> values, FlatVector<> coefs) const;

    virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const;

    virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const;
  };











  /**
     L2 high order 1D finite element
  */
  template <>
  class L2HighOrderFE_Shape<ET_SEGM> : public L2HighOrderFE<ET_SEGM>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const
    {
      Tx lam[2] = { x[0], 1-x[0] };
      INT<2> e = GetEdgeSort (0, vnums);
      LegendrePolynomial (order, lam[e[1]]-lam[e[0]], shape);
    }
  };


  /**
     L2 high order triangular finite element
  */
  template <> 
  class L2HighOrderFE_Shape<ET_TRIG> : public L2HighOrderFE<ET_TRIG>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const
    {
      Tx lam[3] = { x[0], x[1], 1-x[0]-x[1] };
      INT<4> f = GetFaceSort (0, vnums);
      int p = order_inner[0];
      DubinerBasis::Eval (p, lam[f[0]], lam[f[1]], shape);
    }
  };
  
    
  /**
     L2 high order quadrilateral finite element
  */
  template <> 
  class L2HighOrderFE_Shape<ET_QUAD> : public L2HighOrderFE<ET_QUAD>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[], TFA & shape) const
    {
      Tx x = hx[0], y = hx[1];
      
      Tx sigma[4] = {(1-x)+(1-y),x+(1-y),x+y,(1-x)+y};  
      
      INT<4> f = GetFaceSort (0, vnums);  
      
      Tx xi = sigma[f[0]]-sigma[f[1]]; 
      Tx eta = sigma[f[0]]-sigma[f[3]]; 
      
      int n = max(order_inner[0],order_inner[1]);
      ArrayMem<Tx, 20> polx(n+1), poly(n+1);
      
      LegendrePolynomial (n, xi, polx);
      LegendrePolynomial (n, eta, poly);
      
      for (int i = 0, ii = 0; i <= order_inner[0]; i++)
	for (int j = 0; j <= order_inner[1]; j++)
	  shape[ii++] = polx[i] * poly[j];
    }
  };


  /**
     L2 high order tetrahedral finite element
  */
  template <> 
  class L2HighOrderFE_Shape<ET_TET> : public L2HighOrderFE<ET_TET>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const
    {
      Tx lami[4] = { x[0], x[1], x[2], 1-x[0]-x[1]-x[2] };
      
      int sort[4];
      for (int i = 0; i < 4; i++) sort[i] = i;
    
      if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
      if (vnums[sort[2]] > vnums[sort[3]]) Swap (sort[2], sort[3]);
      if (vnums[sort[0]] > vnums[sort[2]]) Swap (sort[0], sort[2]);
      if (vnums[sort[1]] > vnums[sort[3]]) Swap (sort[1], sort[3]);
      if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);

      Tx lamis[4];
      for (int i = 0; i < 4; i++)
	lamis[i] = lami[sort[i]];

      ArrayMem<Tx, 20> memx(sqr(order+1));
      ArrayMem<Tx, 20> memy(sqr(order+1));

      FlatMatrix<Tx> polsx(order+1, &memx[0]);
      FlatMatrix<Tx> polsy(order+1, &memy[0]);
      VectorMem<10, Tx> polsz(order+1);
    
      for (int i = 0; i <= order; i++)
	JacobiPolynomial (order, 2*lamis[0]-1, 2*i+2, 0, polsx.Row(i));
      for (int i = 0; i <= order; i++)
	ScaledJacobiPolynomial (order, lamis[1]-lamis[2]-lamis[3], 1-lamis[0], 2*i+1, 0, polsy.Row(i));

      ScaledLegendrePolynomial (order, lamis[2]-lamis[3], lamis[2]+lamis[3], polsz);

      for (int i = 0, ii = 0; i <= order; i++)
	for (int j = 0; j <= order-i; j++)
	  for (int k = 0; k <= order-i-j; k++, ii++)
	    shape[ii] = polsz(k) * polsy(k, j) * polsx(j+k, i);
    }

  };


  /**
     L2 high order prismatic finite element
  */
  template <> 
  class L2HighOrderFE_Shape<ET_PRISM> : public L2HighOrderFE<ET_PRISM>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[], TFA & shape) const
    {
      Tx lami[3] = { hx[0], hx[1], 1-hx[0]-hx[1] };

      int sort[3];
      for (int i = 0; i < 3; i++) sort[i] = i;
    
      if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);
      if (vnums[sort[1]] > vnums[sort[2]]) Swap (sort[1], sort[2]);
      if (vnums[sort[0]] > vnums[sort[1]]) Swap (sort[0], sort[1]);

      Tx lamis[3];
      for (int i = 0; i < 3; i++)
	lamis[i] = lami[sort[i]];

      Tx x = lamis[0];
      // Tx y = lamis[1];
      Tx z = hx[2];

      int p=order_inner[0];
      int q=order_inner[1];

      ArrayMem<Tx, 20> memx(sqr(p+1));
      FlatMatrix<Tx> polsx(p+1, &memx[0]);

      VectorMem<10, Tx> polsy(p+1);
      VectorMem<10, Tx> polsz(q+1);
    
      for (int i = 0; i <= p; i++)
	JacobiPolynomial (p, 2*x-1, 2*i+1, 0, polsx.Row(i));

      ScaledLegendrePolynomial (order, lamis[1]-lamis[2], lamis[1]+lamis[2], polsy);
      LegendrePolynomial (order, 2*z-1, polsz);

      int ii = 0;
      for (int k = 0; k <= q; k++)
	for (int i = 0; i <= p; i++)
	  for (int j = 0; j <= p-i; j++)
	    shape[ii++] = polsx(j,i) * polsy(j) * polsz(k);
    }

  };


  /**
     L2 high order pyramid finite element
  */
  template <> 
  class L2HighOrderFE_Shape<ET_PYRAMID> : public L2HighOrderFE<ET_PYRAMID>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const;
  };


  /**
     L2 high order hexahedral finite element
  */
  template <> 
  class L2HighOrderFE_Shape<ET_HEX> : public L2HighOrderFE<ET_HEX>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const;
  };
}

#endif
