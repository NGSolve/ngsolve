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
    void T_CalcShape (Tx hx[1], TFA & shape) const;
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
    void T_CalcShape (Tx hx[2], TFA & shape) const;
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
    void T_CalcShape (Tx hx[2], TFA & shape) const;
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
    void T_CalcShape (Tx hx[3], TFA & shape) const;
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
    void T_CalcShape (Tx hx[3], TFA & shape) const;
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
    void T_CalcShape (Tx hx[3], TFA & shape) const;
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
    void T_CalcShape (Tx hx[3], TFA & shape) const;
  };

}

#endif
