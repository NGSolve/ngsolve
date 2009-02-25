#ifndef FILE_HCURLHOFE_
#define FILE_HCURLHOFE_  

// zur Numerischen Berechnung von curl(phi) bei TETS: 
// #define NUMCURLTET 
// #define NUMCURLPRISM 
// #define NUMCURLHEX
#define NUMCURLPYR
#undef newinner // test for hcurl-pyr-inner 
 
/*********************************************************************/
/* File:   hcurlhofe.hpp                                             */
/* Author: Sabine Zaglmayr                                           */
/* Date:   20. Maerz 2003                                            */
/*********************************************************************/
   

/**
   High order H(curl) finite element
 */
template <int D>
class HCurlHighOrderFiniteElement : public HCurlFiniteElement<D> 
{
 
public:
  //  enum { DIM = D };

protected:
  // int ndof; 
  int vnums[8]; 
  INT<3> order_inner;
  INT<2> order_face[6];
  int order_edge[12];
  //int ndof_edge; 
  int ned; // number of edges in element  
  int nv; // number of vertices in element 
  int nf; // number of faces in element
  int usegrad_face[6]; 
  int usegrad_cell; 
  int usegrad_edge[12]; 
  int order_vertex[8];   // for augmented
  int augmented; 
  bool discontinuous;
  
public:
  ///
  HCurlHighOrderFiniteElement (ELEMENT_TYPE aeltype);

  virtual void SetVertexNumbers (FlatArray<int> & avnums, LocalHeap & lh);
  void SetOrderInner (int oi);
  void SetOrderInner (INT<3> oi);
  void SetOrderFace (FlatArray<int> & of);
  void SetOrderFace (FlatArray<INT<2> > & of); 
  void SetOrderEdge (FlatArray<int> & oen);
  void SetUsegradEdge(FlatArray<int> & uge); 
  void SetUsegradFace(FlatArray<int> & ugf); 
  void SetUsegradCell(int ugc); 
  void SetOrderVertex (FlatArray<int> & ov);
  void SetAugmented (int aa);
  
  virtual void ComputeNDof () = 0;

  void PrintInfo() const;

  virtual void SetDiscontinuous ( bool adiscont ) { discontinuous = adiscont; }
};
 
/// 
template <class T_ORTHOPOL> 
class HCurlHighOrderSegm:  public HCurlHighOrderFiniteElement<1>
{
private: 
   typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
 
public:
  
  HCurlHighOrderSegm (int aorder);
  virtual void ComputeNDof();
 
  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<1> shape) const;
  
};

///
template <class T_ORTHOPOL> 
class HCurlHighOrderTrig : public HCurlHighOrderFiniteElement<2>
{
private:
  typedef TrigShapesInnerLegendre T_INNERSHAPES; 
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
public:

  HCurlHighOrderTrig (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;
 
  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<2> shape) const;
 
  /// compute Curl of shape 
  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<1> curlshape) const; 
 
};


///
template <class T_ORTHOPOL>
class HCurlHighOrderQuad : public HCurlHighOrderFiniteElement<2>
{
private: 
   typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
public:
  HCurlHighOrderQuad (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;
 
  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<2> shape) const;

   /// compute Curl of shape 
  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<1> curlshape) const; 
  
};


///
template <class T_ORTHOPOL = IntegratedLegendreMonomialExt> 
class HCurlHighOrderTet : public HCurlHighOrderFiniteElement<3>
{
private:
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
  typedef TetShapesFaceLegendre T_FACESHAPES; 
  typedef TetShapesInnerLegendre T_INNERSHAPES; 
public:
  HCurlHighOrderTet (int aorder);

  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;
 #ifndef NUMCURLTET
  ///compute Curl of shape 
  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<3> shape) const;
#endif 

  virtual Vec<3> EvaluateCurlShape (const IntegrationPoint & ip, 
                                    FlatVector<double> x,
                                    LocalHeap & lh) const;
};

///
template <class T_ORTHOPOL>
class HCurlHighOrderHex : public HCurlHighOrderFiniteElement<3>
{
private: 
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
public:
  HCurlHighOrderHex (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;
 
  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;
#ifndef NUMCURLHEX
  /// compute Curl of shape 
  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<3> curlshape) const;
#endif
};

///
template <class T_ORTHOPOL = IntegratedLegendreMonomialExt> 
class HCurlHighOrderPrism : public HCurlHighOrderFiniteElement<3>
{
private:
  typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
public:
  HCurlHighOrderPrism (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

/// compute curl of shape
#ifndef NUMCURLPRISM
  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<3> shape) const;
#endif
 
};

///
template <class T_ORTHOPOL>
class HCurlHighOrderPyr : public HCurlHighOrderFiniteElement<3>
{
  typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;  
public:

  HCurlHighOrderPyr (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;
 
  /// compute shape
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<3> shape) const;

  /// compute Curl of shape 
#ifndef NUMCURLPYR
  virtual void CalcCurlShape (const IntegrationPoint & ip, 
			      FlatMatrixFixWidth<3> curlshape) const;
#endif 
};



#endif

