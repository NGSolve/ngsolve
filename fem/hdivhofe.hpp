#ifndef FILE_HDIVHOFE_
#define FILE_HDIVHOFE_ 

/*********************************************************************/
/* File:   hdivhofe.hpp                                              */
/* Author: A. Becirovic, S. Zaglmayr, J. Schoeberl                   */
/* Date:   15. Feb. 2003                                             */
/*********************************************************************/


#include "thdivfe.hpp"

namespace ngfem
{
  

#undef HDIV_OLD




///
template <int DIM>
class HDivHighOrderFiniteElement : virtual public HDivFiniteElement<DIM>
{
protected:
  int vnums[8];
 
  INT<3> order_inner;
  INT<2> order_face[6];  // 3D only
  int order_edge[12];   // 2D only

  // bool augmented;

  bool discontinuous;
  bool ho_div_free;

  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;


public:
  HDivHighOrderFiniteElement () { ; }
  HDivHighOrderFiniteElement (ELEMENT_TYPE aeltype);

  void SetVertexNumbers (FlatArray<int> & avnums);
  void SetOrderEdge(FlatArray<int> & oe);
  void SetOrderFace (FlatArray<int> & of);
  void SetOrderInner (int oi);
  void SetOrderFace (FlatArray<INT<2> > & of);
  void SetOrderInner (INT<3> oi); 

  void SetDiscontinuous (bool disc) { discontinuous = disc; };  
  void SetHODivFree (bool aho_div_free) { ho_div_free = aho_div_free; };  

  virtual void ComputeNDof () = 0;
  
  int EdgeOrientation (int enr) const
  {
    const EDGE * edges = ElementTopology::GetEdges (this->eltype);
    return (vnums[edges[enr][1]] > vnums[edges[enr][0]]) ? 1 : -1;
  }
  
  virtual void GetFacetDofs(int i, Array<int> & dnums) const
  { *testout  << " GetFacetDofs for nothing " << endl; dnums.SetSize(0);}; 
};


template <int D>
class HDivHighOrderNormalFiniteElement : public HDivNormalFiniteElement<D>
{

  //public:
  // enum { DIM = D };

protected:
  int vnums[4];
  INT<2> order_inner;

  int ned; // number of edges in element
  int nv; // number of vertices in element


  bool augmented;

public:
  ///
  HDivHighOrderNormalFiniteElement (ELEMENT_TYPE aeltype);


  void SetVertexNumbers (FlatArray<int> & avnums);

  void SetOrderInner (int oi);
  void SetOrderInner (INT<2> oi);

  virtual void ComputeNDof () = 0;

  int EdgeOrientation (int enr) const
  {
    const EDGE * edges = ElementTopology::GetEdges (this->eltype);
    return (vnums[edges[enr][1]] > vnums[edges[enr][0]]) ? 1 : -1;
  }
  
};


template <class T_ORTHOPOL = TrigExtensionMonomial>
class HDivHighOrderNormalSegm : public HDivHighOrderNormalFiniteElement<1>
{
public:

  HDivHighOrderNormalSegm (int aorder);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip,
			   FlatVector<> shape) const;

};

template <class T_ORTHOPOL = TrigExtensionMonomial>
class HDivHighOrderNormalTrig : public HDivHighOrderNormalFiniteElement<2>
{
public:

  HDivHighOrderNormalTrig (int aorder);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip,
			   FlatVector<> shape) const;

};

template <class T_ORTHOPOL = TrigExtensionMonomial>
class HDivHighOrderNormalQuad : public HDivHighOrderNormalFiniteElement<2>
{
public:

  HDivHighOrderNormalQuad (int aorder);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip,
			   FlatVector<> shape) const;

};




template <ELEMENT_TYPE ET> class HDivHighOrderFE;


template <ELEMENT_TYPE ET>
class T_HDivHighOrderFiniteElement 
  : public HDivHighOrderFiniteElement<ET_trait<ET>::DIM>,
    public T_HDivFiniteElement< HDivHighOrderFE<ET>, ET>
    
{
protected:
  enum { DIM = ET_trait<ET>::DIM };
  
  using HDivFiniteElement<DIM>::ndof;
  using HDivFiniteElement<DIM>::order;
  using HDivFiniteElement<DIM>::eltype;
  using HDivFiniteElement<DIM>::dimspace;

  using HDivHighOrderFiniteElement<DIM>::order_edge;
  using HDivHighOrderFiniteElement<DIM>::order_face;
  using HDivHighOrderFiniteElement<DIM>::order_inner;
  using HDivHighOrderFiniteElement<DIM>::ho_div_free;
  using HDivHighOrderFiniteElement<DIM>::discontinuous;


  using HDivHighOrderFiniteElement<DIM>::vnums;
  
  
public:
  T_HDivHighOrderFiniteElement () 
    : HDivHighOrderFiniteElement<DIM> (ET)
  {
    for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
      vnums[i] = i;
    dimspace = DIM;
    eltype = ET;
  }

  T_HDivHighOrderFiniteElement (int aorder) 
  {
    if (DIM == 2)
      for (int i = 0; i < ET_trait<ET>::N_EDGE; i++)
        order_edge[i] = aorder;
    else
      for (int i=0; i < ET_trait<ET>::N_FACE; i++) 
        order_face[i] = INT<2> (aorder,aorder); 
    
    order_inner = INT<3> (aorder,aorder,aorder);

    for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
      vnums[i] = i;
    dimspace = DIM;
    eltype = ET;
  }

  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;


};




template <>
class HDivHighOrderFE<ET_TRIG> : public T_HDivHighOrderFiniteElement<ET_TRIG>
{
public:
  HDivHighOrderFE (int aorder);

  /// compute shape
  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[2], TFA & shape) const; 

  virtual void GetFacetDofs(int i, Array<int> & dnums) const; 
};

///

template <>
class HDivHighOrderFE<ET_QUAD> : public T_HDivHighOrderFiniteElement<ET_QUAD>
{
public:
  HDivHighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[2], TFA & shape) const; 

  virtual void GetFacetDofs(int i, Array<int> & dnums) const; 
};



// template <class T_ORTHOPOL = TrigExtensionMonomial>
template<> 
class HDivHighOrderFE<ET_TET> : public T_HDivHighOrderFiniteElement<ET_TET>
{
   typedef TetShapesInnerLegendre T_INNERSHAPES;
   typedef TetShapesFaceLegendre T_FACESHAPES; 
public:

  HDivHighOrderFE (int aorder);
  // virtual void ComputeNDof();
  // virtual void GetInternalDofs (Array<int> & idofs) const;
  
  /// compute shape
  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[2], TFA & shape) const; 

  virtual void GetFacetDofs(int i, Array<int> & dnums) const; 
};


// template <class T_ORTHOPOL = TrigExtensionMonomial>
template<>
class HDivHighOrderFE<ET_PRISM> : public HDivHighOrderFiniteElement<3>
{
  typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;
public:

  HDivHighOrderFE (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip,
			  FlatMatrixFixWidth<3> shape) const;

  /// compute Div of shape
 // virtual void CalcDivShape (const IntegrationPoint & ip,
  //			     FlatVector<> shape) const;
  /// compute Div numerical diff
  //void CalcNumDivShape( const IntegrationPoint & ip,
  //			FlatVector<> divshape) const;

  virtual void GetFacetDofs(int i, Array<int> & dnums) const; 
};

template<> 
class HDivHighOrderFE<ET_HEX> : public HDivHighOrderFiniteElement<3>
{
public:

  HDivHighOrderFE (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;
  

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip,
			  FlatMatrixFixWidth<3> shape) const;

  /// compute Div of shape
 virtual void CalcDivShape (const IntegrationPoint & ip,
			     FlatVector<> shape) const;
  /// compute Div numerical diff
  //void CalcNumDivShape( const IntegrationPoint & ip,
  //			FlatVector<> divshape) const;
 virtual void GetFacetDofs(int i, Array<int> & dnums) const; 

};

}



#endif



