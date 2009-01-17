#ifndef FILE_HDIVHOFE_
#define FILE_HDIVHOFE_ 

/*********************************************************************/
/* File:   hdivhofe.hpp                                              */
/* Author: Almedin Becirovic                                         */
/* Date:   15. Feb. 2003                                             */
/*********************************************************************/


#undef HDIV_OLD




///
template <int D>
class HDivHighOrderFiniteElement : public HDivFiniteElement<D>
{

  //public:
  // enum { DIM = D };

protected:
public:
  int vnums[8];
 
  #ifdef HDIV_OLD
  int order_inner;
  int order_face[6];
  #else 
  INT<3> order_inner;
  INT<2> order_face[6];
  #endif 
  int order_edge[12];

  int ned; // number of edges in element
  int nv; // number of vertices in element
  int nf; // number of faces in element

  int ndof_edge;
  int ndof_face;
  int ndof_inner;

  bool augmented;
  bool discontinuous;

public:
  ///
  HDivHighOrderFiniteElement (ELEMENT_TYPE aeltype);


  void SetVertexNumbers (FlatArray<int> & avnums);
  void SetOrderEdge(FlatArray<int> & oe);
  void SetOrderFace (FlatArray<int> & of);
  void SetOrderInner (int oi);
#ifndef HDIV_OLD
  void SetOrderFace (FlatArray<INT<2> > & of);
  void SetOrderInner (INT<3> oi); 
#endif
  void SetDiscontinuous(bool disc) { discontinuous=disc; };  

  virtual void ComputeNDof () = 0;

  int EdgeOrientation (int enr) const
  {
    const EDGE * edges = ElementTopology::GetEdges (this->eltype);
    return (vnums[edges[enr][1]] > vnums[edges[enr][0]]) ? 1 : -1;
  }
  
  virtual void GetFacetDofs(int i, ARRAY<int> & dnums) const
  { *testout  << " GetFacetDofs for nothing " << endl; dnums.SetSize(0);}; 
};


template <int D>
class HDivHighOrderNormalFiniteElement : public HDivNormalFiniteElement<D>
{

  //public:
  // enum { DIM = D };

protected:
  int vnums[4];
#ifdef HDIV_OLD
  int order_inner;
#else
  INT<2> order_inner;
#endif
  int ned; // number of edges in element
  int nv; // number of vertices in element


  bool augmented;

public:
  ///
  HDivHighOrderNormalFiniteElement (ELEMENT_TYPE aeltype);


  void SetVertexNumbers (FlatArray<int> & avnums);

  void SetOrderInner (int oi);
#ifndef HDIV_OLD
  void SetOrderInner (INT<2> oi);
#endif
  

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

///
template <class T_ORTHOPOL = TrigExtensionMonomial>
class HDivHighOrderTrig : public HDivHighOrderFiniteElement<2>
{
public:

  HDivHighOrderTrig (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (ARRAY<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip,
			  FlatMatrixFixWidth<2> shape) const;

  /// compute Div of shape
  virtual void CalcDivShape (const IntegrationPoint & ip,
			     FlatVector<> shape) const;
  /// compute Div numerical diff
  //void CalcNumDivShape( const IntegrationPoint & ip,
  //			FlatVector<> divshape) const;

  virtual void GetFacetDofs(int i, ARRAY<int> & dnums) const; 
};

///

template <class T_ORTHOPOL = TrigExtensionMonomial>
class HDivHighOrderQuad : public HDivHighOrderFiniteElement<2>
{
public:

  HDivHighOrderQuad (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (ARRAY<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip,
			  FlatMatrixFixWidth<2> shape) const;


  /// compute Div of shape
  virtual void CalcDivShape (const IntegrationPoint & ip,
			     FlatVector<> shape) const;
  /// compute Div numerical diff
 void CalcNumDivShape( const IntegrationPoint & ip,
  			FlatVector<> divshape) const;

 virtual void GetFacetDofs(int i, ARRAY<int> & dnums) const; 
};

template <class T_ORTHOPOL = TrigExtensionMonomial>
class HDivHighOrderTet : public HDivHighOrderFiniteElement<3>
{
   typedef TetShapesInnerLegendre T_INNERSHAPES;
   typedef TetShapesFaceLegendre T_FACESHAPES; 
public:

  HDivHighOrderTet (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (ARRAY<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip,
			  FlatMatrixFixWidth<3> shape) const;

  /// compute Div of shape
  virtual void CalcDivShape (const IntegrationPoint & ip,
			     FlatVector<> shape) const;
  /// compute Div numerical diff
  //void CalcNumDivShape( const IntegrationPoint & ip,
  //			FlatVector<> divshape) const;

  virtual void GetFacetDofs(int i, ARRAY<int> & dnums) const; 
};

template <class T_ORTHOPOL = TrigExtensionMonomial>
class HDivHighOrderPrism : public HDivHighOrderFiniteElement<3>
{
  typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;
public:

  HDivHighOrderPrism (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (ARRAY<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip,
			  FlatMatrixFixWidth<3> shape) const;

  /// compute Div of shape
 // virtual void CalcDivShape (const IntegrationPoint & ip,
  //			     FlatVector<> shape) const;
  /// compute Div numerical diff
  //void CalcNumDivShape( const IntegrationPoint & ip,
  //			FlatVector<> divshape) const;

  virtual void GetFacetDofs(int i, ARRAY<int> & dnums) const; 
};

template <class T_ORTHOPOL = TrigExtensionMonomial>
class HDivHighOrderHex : public HDivHighOrderFiniteElement<3>
{
public:

  HDivHighOrderHex (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (ARRAY<int> & idofs) const;
  

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip,
			  FlatMatrixFixWidth<3> shape) const;

  /// compute Div of shape
 virtual void CalcDivShape (const IntegrationPoint & ip,
			     FlatVector<> shape) const;
  /// compute Div numerical diff
  //void CalcNumDivShape( const IntegrationPoint & ip,
  //			FlatVector<> divshape) const;
 virtual void GetFacetDofs(int i, ARRAY<int> & dnums) const; 

};





#endif



