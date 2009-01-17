#ifndef FILE_L2HOFE
#define FILE_L2HOFE

/*********************************************************************/
/* File:   l2hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/



/**
  High order finite elements for L_2
*/
template<int DIM>
class L2HighOrderFiniteElement : public NodalFiniteElement<DIM>
{
public:
  
  int vnums[8]; //he: make oriented L2 elements => can be used for facets
  INT<3> order_inner; 


  L2HighOrderFiniteElement (int dim, ELEMENT_TYPE aeltype);

  virtual void SetVertexNumbers (FlatArray<int> & avnums);
  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);

  void SetOrder (int o);
  
  void SetOrderInner (INT<3> oi);

  virtual void ComputeNDof () = 0; 

  // const int * GetVNums() const { return vnums; }

  virtual void GetInternalDofs (ARRAY<int> & idofs) const; 
  
private:
  
  virtual const ARRAY<typename NodalFiniteElement<DIM>::IPData> & GetIPData () const
  {
    throw Exception ("GetIPData not available for L2HighOrderFE");
  }
};

/**
  High order 1D finite element
*/
class L2HighOrderSegm : public L2HighOrderFiniteElement<1>
{
public:
  L2HighOrderSegm (int aorder=0);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;

  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;
};


/**
  High order triangular finite element
*/
class L2HighOrderTrig : public L2HighOrderFiniteElement<2>
{
public:
  L2HighOrderTrig (int aorder=0);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;

  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;
};


/**
  High order quadrilateral finite element
*/
class L2HighOrderQuad : public L2HighOrderFiniteElement<2>
{
  //typedef TrigShapesInnerLegendre T_INNERSHAPES;
  // typedef TrigShapesInnerJacobi T_INNERSHAPES;
public:
  L2HighOrderQuad (int aorder=0);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;

  /// compute gradient of shape 
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const; 
};


/**
  High order tetrahedral finite element
 */
class L2HighOrderTet : public L2HighOrderFiniteElement<3>
{
  public:
    L2HighOrderTet (int aorder=0);
    virtual void ComputeNDof();

  /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
                            FlatVector<> shape) const;

  /// compute gradient of shape
    virtual void CalcDShape (const IntegrationPoint & ip, FlatMatrix<> dshape) const;
};



/**
  High order prismatic finite element
 */
class L2HighOrderPrism : public L2HighOrderFiniteElement<3>
{
  public:
    L2HighOrderPrism (int aorder=0);
    virtual void ComputeNDof();

  /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
                            FlatVector<> shape) const;

  /// compute gradient of shape
    virtual void CalcDShape (const IntegrationPoint & ip, FlatMatrix<> dshape) const;
};


/**
  High order hexahedral finite element
 */
class L2HighOrderHex : public L2HighOrderFiniteElement<3>
{
  public:
    L2HighOrderHex (int aorder=0);
    virtual void ComputeNDof();

  /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
                            FlatVector<> shape) const;

  /// compute gradient of shape
    virtual void CalcDShape (const IntegrationPoint & ip, FlatMatrix<> dshape) const;
};




#ifdef abc
/**
  High order quadrilateral finite element
*/
template <class T_EXT>
class L2HighOrderQuad : public L2HighOrderFiniteElement<2>
{
  typedef T_EXT T_EDGESHAPES;
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;

public:
  L2HighOrderQuad (int aorder);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
 
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

private:
  void CalcShapeDShape (const IntegrationPoint & ip, 
			ARRAY<AutoDiff<2> > & shape) const;
};


/**
  High order tetrahedral finite element
*/
template <class T_EXT>
class L2HighOrderTet : public L2HighOrderFiniteElement<3>
{
  typedef T_EXT T_EDGESHAPES;

  typedef TetShapesInnerLegendre T_INNERSHAPES;
  typedef TetShapesFaceLegendre T_FACESHAPES;

  // typedef TetShapesInnerJacobi T_INNERSHAPES;
  // typedef TetShapesFaceJacobi T_FACESHAPES;

  // typedef TetShapesFaceOpt1 T_FACESHAPES;

  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
  // typedef VertexStandard T_VERTEXSHAPES;
public:
  L2HighOrderTet (int aorder);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

private:
  void CalcShapeDShape (const IntegrationPoint & ip, 
			ARRAY<AutoDiff<3> > & shape) const;
};


/**
  High order prismatic finite element
*/
template <class T_EXT>
class L2HighOrderPrism : public L2HighOrderFiniteElement<3>
{
  typedef T_EXT T_EDGESHAPES;
  typedef TrigShapesInnerLegendre T_TRIGSHAPES;
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;

public:
  L2HighOrderPrism (int aorder);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

  void CalcShapeDShape (const IntegrationPoint & ip, 
			ARRAY<AutoDiff<3> > & shape) const;
};



/**
  High order hexahedral finite element
*/
template <class T_EXT>
class L2HighOrderHex : public L2HighOrderFiniteElement<3>
{
  typedef T_EXT T_EDGESHAPES;
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;

public:
  L2HighOrderHex (int aorder);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

  void CalcShapeDShape (const IntegrationPoint & ip, 
			ARRAY<AutoDiff<3> > & shape) const;
};


/**
  High order pyramid finite element
*/
template <class T_EXT>
class L2HighOrderPyramid : public L2HighOrderFiniteElement<3>
{
  typedef T_EXT T_EDGESHAPES;
  typedef TrigShapesInnerLegendre T_TRIGSHAPES;
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;

public:
  L2HighOrderPyramid (int aorder);
  virtual void ComputeNDof();

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /*
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;
  */
};
#endif


#endif
