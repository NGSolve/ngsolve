#ifndef FILE_HCURLHOFE_
#define FILE_HCURLHOFE_  

/*********************************************************************/
/* File:   hcurlhofe.hpp                                             */
/* Author: Sabine Zaglmayr                                           */
/* Date:   20. Maerz 2003                                            */
/*                                                                   */
/* AutoCurl - revision: J. Schoeberl, March 2009                     */
/*********************************************************************/
   
namespace ngfem
{

/**
   High order H(curl) finite element of dimension D
 */
template <int D>
class HCurlHighOrderFiniteElement : public HCurlFiniteElement<D> 
{
protected:
  int vnums[8]; 
  int order_edge[12];
  INT<2> order_face[6];
  INT<3> order_cell;

  int usegrad_edge[12]; 
  int usegrad_face[6]; 
  int usegrad_cell; 

  bool discontinuous;
  
public:
  HCurlHighOrderFiniteElement (ELEMENT_TYPE aeltype);
  HCurlHighOrderFiniteElement () { discontinuous = false; }

  void SetVertexNumber (int nr, int vnum) { vnums[nr] = vnum; }
  void SetOrderEdge (int nr, int order) { order_edge[nr] = order; }
  void SetOrderFace (int nr, INT<2> order) { order_face[nr] = order; }

  void SetUseGradEdge(int nr, int uge) { usegrad_edge[nr] = uge; }
  void SetUseGradFace(int nr, int ugf) { usegrad_face[nr] = ugf; }

  void SetVertexNumbers (FlatArray<int> & avnums);
  void SetOrderCell (int oi);
  void SetOrderCell (INT<3> oi);
  void SetOrderFace (FlatArray<int> & of);
  void SetOrderFace (FlatArray<INT<2> > & of); 
  void SetOrderEdge (FlatArray<int> & oen);
  void SetUsegradEdge(FlatArray<int> & uge); 
  void SetUsegradFace(FlatArray<int> & ugf); 
  void SetUsegradCell(int ugc); 
  void SetDiscontinuous ( bool adiscont ) { discontinuous = adiscont; }
  
  virtual void ComputeNDof () = 0;
  void PrintInfo() const;
};


  /** 
      HCurlHighOrderFE of shape ET.
      The template specialization provides the shape functions.
   */
template <ELEMENT_TYPE ET> class HCurlHighOrderFE;



  /**
     HCurlHighOrderFE of shape ET.
     provides access functions, shape funcitons are provided by CalcShape template
  */
template <ELEMENT_TYPE ET>
class T_HCurlHighOrderFiniteElement 
  : public HCurlHighOrderFiniteElement<ET_trait<ET>::DIM>, public ET_trait<ET> 

{
protected:
  enum { DIM = ET_trait<ET>::DIM };
  
  using HCurlFiniteElement<DIM>::DIM_CURL;
  using HCurlFiniteElement<DIM>::ndof;
  using HCurlFiniteElement<DIM>::order;
  using HCurlFiniteElement<DIM>::eltype;
  // using HCurlFiniteElement<DIM>::dimspace;

  using HCurlHighOrderFiniteElement<DIM>::vnums;
  using HCurlHighOrderFiniteElement<DIM>::order_edge;
  using HCurlHighOrderFiniteElement<DIM>::order_face;
  using HCurlHighOrderFiniteElement<DIM>::order_cell;

  using HCurlHighOrderFiniteElement<DIM>::usegrad_edge;
  using HCurlHighOrderFiniteElement<DIM>::usegrad_face;
  using HCurlHighOrderFiniteElement<DIM>::usegrad_cell;

  using HCurlHighOrderFiniteElement<DIM>::discontinuous;


  using ET_trait<ET>::N_VERTEX;
  using ET_trait<ET>::N_EDGE;
  using ET_trait<ET>::N_FACE;
  using ET_trait<ET>::FaceType;
  using ET_trait<ET>::GetEdgeSort;
  using ET_trait<ET>::GetFaceSort;
  


  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

  // typedef TrigShapesInnerLegendre T_TRIGSHAPES;
  // typedef TrigShapesInnerJacobi T_TRIGSHAPES;

public:

  T_HCurlHighOrderFiniteElement () 
  {
    for (int i = 0; i < N_VERTEX; i++)
      vnums[i] = i;
    // dimspace = DIM;
    eltype = ET;
  }

  T_HCurlHighOrderFiniteElement (int aorder);

  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  virtual void CalcShape (const IntegrationPoint & ip, 
                          FlatMatrixFixWidth<DIM> shape) const;

  virtual void CalcCurlShape (const IntegrationPoint & ip, 
                              FlatMatrixFixWidth<DIM_CURL> curlshape) const;

  virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                                FlatMatrixFixWidth<DIM> shape) const;

  virtual void CalcMappedCurlShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                                    FlatMatrixFixWidth<DIM_CURL> curlshape) const;

  /*
  virtual Vec <DIM_CURL_TRAIT<ET_trait<ET>::DIM>::DIM>
  EvaluateCurlShape (const IntegrationPoint & ip, 
                     FlatVector<double> x,
                     LocalHeap & lh) const;
  */
};

 
/// A segment high order H(curl) element
template <>
class HCurlHighOrderFE<ET_SEGM>:  public HCurlHighOrderFiniteElement<1>
{
private: 
  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
 
public:
  HCurlHighOrderFE ();
  HCurlHighOrderFE (int aorder);
  virtual void ComputeNDof();
 
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatMatrixFixWidth<1> shape) const;
};

/// A triangular high order H(curl) element
template <>
class HCurlHighOrderFE<ET_TRIG> : public T_HCurlHighOrderFiniteElement<ET_TRIG>
{
private:
  typedef TrigShapesInnerLegendre T_INNERSHAPES; 
  // typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
public:
  HCurlHighOrderFE () { ; }
  HCurlHighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[2], TFA & shape) const; 
};

/// A quadrilateral high order H(curl) element
template <>
class HCurlHighOrderFE<ET_QUAD> : public T_HCurlHighOrderFiniteElement<ET_QUAD>
{
private: 
  // typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
public:
  HCurlHighOrderFE () { ; }
  HCurlHighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[2], TFA & shape) const; 
};

/// A tetrahedral high order H(curl) element
template <>
class HCurlHighOrderFE<ET_TET> : public T_HCurlHighOrderFiniteElement<ET_TET>
{
private:
  typedef TetShapesFaceLegendre T_FACESHAPES; 
  typedef TetShapesInnerLegendre T_INNERSHAPES; 
public:
  HCurlHighOrderFE () { ; }
  HCurlHighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[3], TFA & shape) const; 
};

/// A hexahedral high order H(curl) element
template <>
class HCurlHighOrderFE<ET_HEX> : public T_HCurlHighOrderFiniteElement<ET_HEX>
{
public:
  HCurlHighOrderFE () { ; }
  HCurlHighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[3], TFA & shape) const; 
};

/// A prismatic high order H(curl) element
template <>
class HCurlHighOrderFE<ET_PRISM> : public T_HCurlHighOrderFiniteElement<ET_PRISM>
{
private:
  typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;

public:
  HCurlHighOrderFE () { ; }
  HCurlHighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[3], TFA & shape) const; 
};

/// A pyramidal high order H(curl) element
template <>
class HCurlHighOrderFE<ET_PYRAMID> : public T_HCurlHighOrderFiniteElement<ET_PYRAMID>
{
  typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;  
public:
  HCurlHighOrderFE () { ; }
  HCurlHighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[3], TFA & shape) const; 
};

}

#endif

