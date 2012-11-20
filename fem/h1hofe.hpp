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
    /// global vertex numbers for shape orientation
    int vnums[8];
    /// order of internal shapes (3d only)
    INT<3> order_cell;
    /// order of face shapes / internal in 2d
    INT<2> order_face[6];
    /// order of edge shapes
    int order_edge[12];

    using ScalarFiniteElement<DIM>::eltype;

    bool nodalp2;  // 

  public:
    H1HighOrderFiniteElement () 
      : nodalp2(false) { ; }

    /// assignes vertex numbers
    template <typename TA> 
    void SetVertexNumbers (const TA & avnums)
    { for (int i = 0; i < avnums.Size(); i++) vnums[i] = avnums[i]; }

    /// assign vertex number
    void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }

    /// set anisotropic cell order
    void SetOrderCell (INT<3> oi)  { order_cell = oi; }

    /// set isotropic or anisotropic face orders
    template <typename TA>
    void SetOrderFace (const TA & of)
    { for (int i = 0; i < of.Size(); i++) order_face[i] = of[i]; }

    /// set anisotropic face order for face nr
    void SetOrderFace (int nr, INT<2> order) { order_face[nr] = order; }

    /// set edge orders
    template <typename TA>
    void SetOrderEdge (const TA & oe)
    { for (int i = 0; i < oe.Size(); i++) order_edge[i] = oe[i]; }

    /// set edge order for edge nr
    void SetOrderEdge (int nr, int order) { order_edge[nr] = order; }

    void SetNodalP2 (bool anp2) { nodalp2 = anp2; }

    /// high order elements need extra configuration. update ndof and order
    virtual void ComputeNDof () = 0;
  };




  /**
     Barton-Nackman base class for H1 - high order finite elements
  */
  template <ELEMENT_TYPE ET>
  class T_H1HighOrderFiniteElement : 
    public H1HighOrderFiniteElement<ET_trait<ET>::DIM>, 
    public ET_trait<ET> 
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };

    using ScalarFiniteElement<DIM>::ndof;
    using ScalarFiniteElement<DIM>::order;
    using ScalarFiniteElement<DIM>::eltype;

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
    using ET_trait<ET>::PolDimension;
    using ET_trait<ET>::PolBubbleDimension;

  public:

    T_H1HighOrderFiniteElement () { ; }

    T_H1HighOrderFiniteElement (int aorder) 
    {
      ndof = PolDimension (aorder);

      for (int i = 0; i < N_VERTEX; i++) vnums[i] = i;
      for (int i = 0; i < N_EDGE; i++) order_edge[i] = aorder;
      for (int i = 0; i < N_FACE; i++) order_face[i] = aorder;   
      if (DIM == 3) order_cell = aorder; 

      order = aorder;
    }

    virtual void ComputeNDof();
  };





  /// shape function engine for high order h1-elements
  template <ELEMENT_TYPE ET> class H1HighOrderFE_Shape;


  /**
     High order finite elements for H1.
     These are the actual finite element classes to be used.
     Operations are inherited from T_ScalarFiniteElement2, and shape functions are provided by the shape template
   */
  template <ELEMENT_TYPE ET> 
  class  NGS_DLL_HEADER H1HighOrderFE :  public T_H1HighOrderFiniteElement<ET>,
					 public T_ScalarFiniteElement2< H1HighOrderFE_Shape<ET>, ET >

  {   
  public:
    /// minimal constructor, orders will be set later
    H1HighOrderFE () { ; }

    /// builds a functional element of order aorder.
    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement<ET> (aorder) { ; }
  };



  /**
     High order 0D finite element
  */

  template <> 
  class H1HighOrderFE_Shape<ET_POINT> : public H1HighOrderFE<ET_POINT>
  {
  public:
    /// generic shape function
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[], TFA & shape) const
    {
      shape[0] = 1.0;
    }
  };



  /**
     High order segment finite element
  */

  template <> 
  class H1HighOrderFE_Shape<ET_SEGM> : public H1HighOrderFE<ET_SEGM>
  {
  public:
    /// generic shape function
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const;
  };


  /**
     High order triangular finite element
  */
  template <>
  class H1HighOrderFE_Shape<ET_TRIG> : public H1HighOrderFE<ET_TRIG>
  {
  public:
    /// generic shape function
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const;
  };


  /**
     High order quadrilateral finite element
  */
  template <>
  class NGS_DLL_HEADER H1HighOrderFE_Shape<ET_QUAD> : public H1HighOrderFE<ET_QUAD>
  {
  public:
    /// generic shape function
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[], TFA & shape) const;
  };


  /**
     High order tetrahedral finite element
  */
  template <>
  class NGS_DLL_HEADER H1HighOrderFE_Shape<ET_TET> : public H1HighOrderFE<ET_TET>
  {
  public:
    /// generic shape function
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[], TFA & shape) const; 
  };


  /** 
      High order prismatic finite element
  */
  template <>
  class H1HighOrderFE_Shape<ET_PRISM> : public H1HighOrderFE<ET_PRISM>
  {
  public:
    /// generic shape function
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[], TFA & shape) const; 
  };



  /**
     High order hexahedral finite element
  */
  template <> 
  class H1HighOrderFE_Shape<ET_HEX> : public H1HighOrderFE<ET_HEX>
  {
  public:
    /// generic shape function
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[], TFA & shape) const; 
  };


  /**
     High order pyramid shape functions
  */
  template <> 
  class H1HighOrderFE_Shape<ET_PYRAMID> : public H1HighOrderFE<ET_PYRAMID>
  {
  public:
    /// generic shape function
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[], TFA & shape) const; 
  };
}


#endif
