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
    using ET_trait<ET>::PolDimension;



  public:

    T_L2HighOrderFiniteElement () 
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
	vnums[i] = i;
      dimspace = DIM;
      eltype = ET;
    }

    T_L2HighOrderFiniteElement (int aorder) 
    {
      for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
	vnums[i] = i;
      dimspace = DIM;
      eltype = ET;

      order = aorder;
      order_inner = aorder;
      ndof = PolDimension (aorder);
    }

    virtual void ComputeNDof();
  };






  template <ELEMENT_TYPE ET> class L2HighOrderFE_Shape;



  template <ELEMENT_TYPE ET> 
  class  NGS_DLL_HEADER L2HighOrderFE :  public T_L2HighOrderFiniteElement<ET>,
					 public T_ScalarFiniteElement2< L2HighOrderFE_Shape<ET>, ET >

  {   
  public:
    L2HighOrderFE () { ; }

    L2HighOrderFE (int aorder)
      : T_L2HighOrderFiniteElement<ET> (aorder) { ; }
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
    void T_CalcShape (Tx x[], TFA & shape) const;
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
