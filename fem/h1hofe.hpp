#ifndef FILE_H1HOFE
#define FILE_H1HOFE

/*********************************************************************/
/* File:   h1hofe.hpp                                                */
/* Author: Start                                                     */
/* Date:   6. Feb. 2003                                              */
/*********************************************************************/

/**
  High order finite elements for H^1
*/
template<int DIM>
class H1HighOrderFiniteElement : public NodalFiniteElement<DIM>
{
public:
  int vnums[8];
  INT<3> order_inner;
  INT<2> order_face[6];
  int order_edge[12];
  int order_vertex[8];   // for augmented

  int augmented; // 0..l.o, 1..l.o + h.o,  2..all monomials

public:
  ///
  H1HighOrderFiniteElement (ELEMENT_TYPE aeltype);
  
  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);
  virtual void GetDofs (Array<Dof> & dofs) const;

  void SetOrderInner (int oi);
  void SetOrderInner (INT<3> oi);
  void SetOrderFace (FlatArray<int> & of);
  void SetOrderFace (FlatArray<INT<2> > & of);
  void SetOrderEdge (FlatArray<int> & oe);
  void SetOrderVertex (FlatArray<int> & ov);
  void SetAugmented (int aa);

  virtual void ComputeNDof () = 0;
  const int * GetVNums() const { return vnums; }
};



/**
  High order segment finite element
*/
template <class T_ORTHOPOL = IntegratedLegendreMonomialExt>
class H1HighOrderSegm : public H1HighOrderFiniteElement<1>
{
private:
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;

public:
  H1HighOrderSegm (int aorder);
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
template <class T_ORTHOPOL = IntegratedLegendreMonomialExt>
class H1HighOrderTrig : public H1HighOrderFiniteElement<2>
{
  typedef TrigShapesInnerLegendre T_INNERSHAPES;
  // typedef TrigShapesInnerJacobi T_INNERSHAPES;

  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
  // typedef VertexStandard T_VERTEXSHAPES;
  

public:

  H1HighOrderTrig (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

private: 

  template<typename Tx, typename Ty, typename TFA>  
  void T_CalcShape (Tx x, Ty y, TFA & shape) const; 
};


/**
  High order quadrilateral finite element
*/
template <class T_ORTHOPOL = IntegratedLegendreMonomialExt>
class H1HighOrderQuad : public H1HighOrderFiniteElement<2>
{
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;

public:
  H1HighOrderQuad (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

private:
  template<typename Tx, typename Ty, typename TFA>  
  void T_CalcShape (Tx x, Ty y, TFA & shape) const; 

  //  void CalcShapeDShape (const IntegrationPoint & ip, 
  //			Array<AutoDiff<2> > & shape) const;
};


/**
  High order tetrahedral finite element
*/


template <class T_ORTHOPOL = IntegratedLegendreMonomialExt>
class H1HighOrderTet : public H1HighOrderFiniteElement<3>
{
  typedef TetShapesInnerLegendre T_INNERSHAPES;
  typedef TetShapesFaceLegendre T_FACESHAPES;

  // typedef TetShapesInnerJacobi T_INNERSHAPES;
  // typedef TetShapesFaceJacobi T_FACESHAPES;

  // typedef TetShapesFaceOpt1 T_FACESHAPES;

  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
  // typedef VertexStandard T_VERTEXSHAPES;
public:
  H1HighOrderTet (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

  /// compute dshape, matrix: ndof x spacedim
  virtual void CalcMappedDShape (const BaseSpecificIntegrationPoint & sip, 
                                 FlatMatrix<> dshape) const;


private:
  // void CalcShapeDShape (const IntegrationPoint & ip, 
  // 			Array<AutoDiff<3> > & shape) const;
  
  
  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void T_CalcShape (Tx x, Ty y, Tz z, TFA & shape) const; 
};


/** 
  High order prismatic finite element
*/
template <class T_ORTHOPOL = IntegratedLegendreMonomialExt>
class H1HighOrderPrism : public H1HighOrderFiniteElement<3>
{
  typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;
  //typedef TrigShapesInnerJacobi T_TRIGFACESHAPES;
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
  bool plate;   // thin plate: top and bottom are internal dofs

public:
  H1HighOrderPrism (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

  void SetPlate (int aplate) { plate = aplate; }

private:
  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void T_CalcShape (Tx x, Ty y, Tz z, TFA & shape) const; 
  
  //  void CalcShapeDShape (const IntegrationPoint & ip, 
  //			Array<AutoDiff<3> > & shape) const;
};



/**
  High order hexahedral finite element
*/
template <class T_ORTHOPOL = IntegratedLegendreMonomialExt>
class H1HighOrderHex : public H1HighOrderFiniteElement<3>
{
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;

public:
  H1HighOrderHex (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

private:
  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void T_CalcShape (Tx x, Ty y, Tz z, TFA & shape) const; 

  //  void CalcShapeDShape (const IntegrationPoint & ip, 
  //  			Array<AutoDiff<3> > & shape) const;
};


/**
  High order pyramid finite element
*/
template <class T_ORTHOPOL = IntegratedLegendreMonomialExt>
class H1HighOrderPyramid : public H1HighOrderFiniteElement<3>
{
  typedef TrigShapesInnerLegendre T_TRIGSHAPES;
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;

public:
  H1HighOrderPyramid (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;
  
  /// compute gradient of shape
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   FlatMatrix<> dshape) const;

private:
  template<typename Tx, typename Ty, typename Tz, typename TFA>  
  void T_CalcShape (Tx x, Ty y, Tz z, TFA & shape) const; 

  // void CalcShapeDShape (const IntegrationPoint & ip, 
  // Array<AutoDiff<3> > & shape) const;
};




#endif
