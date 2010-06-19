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
    }

    T_H1HighOrderFiniteElement (int aorder) 
    {
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
  };




  template <ELEMENT_TYPE ET>
  class T_H1HighOrderFiniteElement2 
    : public T_H1HighOrderFiniteElement<ET>,
      public T_ScalarFiniteElement2< H1HighOrderFE<ET>, ET >
  { 
  public:    
    using T_H1HighOrderFiniteElement<ET>::GetInternalDofs;

    T_H1HighOrderFiniteElement2 () { ; }
    T_H1HighOrderFiniteElement2 (int aorder) 
      :  T_H1HighOrderFiniteElement<ET> (aorder) 
    { ; }
  };





  /**
     High order segment finite element
  */

  template <> 
  class NGS_DLL_HEADER H1HighOrderFE<ET_SEGM> : public T_H1HighOrderFiniteElement2<ET_SEGM>
  {
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement2<ET_SEGM> (aorder) 
    { ndof = (order+1); }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[1], TFA & shape) const;
  };


  /**
     High order triangular finite element
  */
  template <>
  class NGS_DLL_HEADER H1HighOrderFE<ET_TRIG> 
    : public T_H1HighOrderFiniteElement2<ET_TRIG>
  {
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement2<ET_TRIG> (aorder) 
    { ndof = (order+1)*(order+2)/2; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[2], TFA & shape) const;
  };


  /**
     High order quadrilateral finite element
  */
  template <>
  class NGS_DLL_HEADER H1HighOrderFE<ET_QUAD> : public T_H1HighOrderFiniteElement2<ET_QUAD>
  {
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement2<ET_QUAD> (aorder) 
    { ndof = (order+1)*(order+1); }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[2], TFA & shape) const;
  };


  /**
     High order tetrahedral finite element
  */
  template <>
  class NGS_DLL_HEADER H1HighOrderFE<ET_TET> : public T_H1HighOrderFiniteElement2<ET_TET>
  {
    typedef TetShapesInnerLegendre T_INNERSHAPES;
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
    : T_H1HighOrderFiniteElement2<ET_TET> (aorder) 
    { ndof = (aorder+1)*(aorder+2)*(aorder+3)/6; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const; 
  };


  /** 
      High order prismatic finite element
  */
  template <>
  class NGS_DLL_HEADER H1HighOrderFE<ET_PRISM> : public T_H1HighOrderFiniteElement2<ET_PRISM>
  {
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement2<ET_PRISM> (aorder) 
    { ndof = (order+1)*(order+2)*(order+1)/2; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const; 
  };



  /**
     High order hexahedral finite element
  */
  template <> 
  class NGS_DLL_HEADER H1HighOrderFE<ET_HEX> : public T_H1HighOrderFiniteElement2<ET_HEX>
  {
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement2<ET_HEX> (aorder) 
    { ndof = (order+1)*(order+1)*(order+1); }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const; 
  };


  /**
     High order pyramid finite element
  */
  template<>
  class NGS_DLL_HEADER H1HighOrderFE<ET_PYRAMID> : public T_H1HighOrderFiniteElement2<ET_PYRAMID>
  {
  public:
    H1HighOrderFE () { ; }

    H1HighOrderFE (int aorder)
      : T_H1HighOrderFiniteElement2<ET_PYRAMID> (aorder) 
    { ndof = (order+2)*(order+1)*(2*order+3) / 6; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[3], TFA & shape) const; 
  };
}


#endif
