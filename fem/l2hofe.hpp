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

    virtual void EvaluateGrad (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<DIM> values) const
    {
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	{
	  FlatVector<> vval(DIM*values.Height(), &values(0,0));
	  vval = pre->dshapes * coefs;
	}
      else
	T_ScalarFiniteElement2< SHAPES<ET>, ET > :: EvaluateGrad (ir, coefs, values);
    }


    virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> values, FlatVector<> coefs) const
    {
      int classnr =  ET_trait<ET>::GetClassNr (vnums);

      PrecomputedScalShapes<DIM> * pre = precomp.Get (classnr, order, ir.GetNIP());
      if (pre)
	coefs = Trans (pre->dshapes) * FlatVector<> (DIM*ndof, &values(0,0));  // values.Height !!!
      else
	T_ScalarFiniteElement2< SHAPES<ET>, ET > :: EvaluateGradTrans (ir, values, coefs);
    }

    virtual void GetTrace (int facet, FlatVector<> coefs, FlatVector<> fcoefs) const;
    /*
    {
      int classnr =  ET_trait<ET>::GetFacetClassNr (facet, vnums);
      if (precomp_trace.Used (INT<2> (order, classnr)))
	{
	  fcoefs = (*precomp_trace.Get (INT<2> (order, classnr))) * coefs;
	}
      else
	L2HighOrderFiniteElement<DIM>::GetTrace (facet, coefs, fcoefs);
    }
    */
    virtual void GetTraceTrans (int facet, FlatVector<> fcoefs, FlatVector<> coefs) const;
    /*
    {
      int classnr =  ET_trait<ET>::GetFacetClassNr (facet, vnums);
      if (precomp_trace.Used (INT<2> (order, classnr)))
	{
	  coefs = Trans(*precomp_trace.Get (INT<2> (order, classnr))) * fcoefs;
	}
      else
	L2HighOrderFiniteElement<DIM>::GetTraceTrans (facet, fcoefs, coefs);
    }
    */
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
    void T_CalcShape (Tx x[], TFA & shape) const;
  };


  /**
     L2 high order prismatic finite element
  */
  template <> 
  class L2HighOrderFE_Shape<ET_PRISM> : public L2HighOrderFE<ET_PRISM>
  {
  public:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const;
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
