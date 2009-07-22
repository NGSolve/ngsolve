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
   High order H(curl) finite element
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



template <ELEMENT_TYPE ET> class HCurlHighOrderFE;

template <ELEMENT_TYPE ET>
class T_HCurlHighOrderFiniteElement 
  : public HCurlHighOrderFiniteElement<ET_trait<ET>::DIM>
{
protected:
  enum { DIM = ET_trait<ET>::DIM };
  
  using HCurlFiniteElement<DIM>::DIM_CURL;
  using HCurlFiniteElement<DIM>::ndof;
  using HCurlFiniteElement<DIM>::order;
  using HCurlFiniteElement<DIM>::eltype;
  using HCurlFiniteElement<DIM>::dimspace;

  using HCurlHighOrderFiniteElement<DIM>::vnums;
  using HCurlHighOrderFiniteElement<DIM>::order_edge;
  using HCurlHighOrderFiniteElement<DIM>::order_face;
  using HCurlHighOrderFiniteElement<DIM>::order_cell;

  using HCurlHighOrderFiniteElement<DIM>::usegrad_edge;
  using HCurlHighOrderFiniteElement<DIM>::usegrad_face;
  using HCurlHighOrderFiniteElement<DIM>::usegrad_cell;

  using HCurlHighOrderFiniteElement<DIM>::discontinuous;


  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

  // typedef TrigShapesInnerLegendre T_TRIGSHAPES;
  // typedef TrigShapesInnerJacobi T_TRIGSHAPES;

public:

  T_HCurlHighOrderFiniteElement () 
  {
    for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
      vnums[i] = i;
    dimspace = DIM;
    eltype = ET;
  }

  T_HCurlHighOrderFiniteElement (int aorder) 
  {
    for (int i = 0; i < ET_trait<ET>::N_EDGE; i++)
      order_edge[i] = aorder;
    for (int i=0; i < ET_trait<ET>::N_FACE; i++) 
      order_face[i] = INT<2> (aorder,aorder); 
    if (DIM == 3)
      order_cell = INT<3> (aorder,aorder,aorder);

    for(int i = 0; i < ET_trait<ET>::N_EDGE; i++)
      usegrad_edge[i] = 1;
    for(int i=0; i < ET_trait<ET>::N_FACE; i++)
      usegrad_face[i] = 1;
    if (DIM == 3)
      usegrad_cell = 1;

    for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
      vnums[i] = i;
    dimspace = DIM;
    eltype = ET;
  }


  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  virtual void CalcShape (const IntegrationPoint & ip, 
                          FlatMatrixFixWidth<DIM> shape) const;

  virtual void CalcCurlShape (const IntegrationPoint & ip, 
                              FlatMatrixFixWidth<DIM_CURL> curlshape) const;

  virtual void CalcMappedShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
                                FlatMatrixFixWidth<DIM> shape) const;

  virtual void CalcMappedCurlShape (const SpecificIntegrationPoint<DIM,DIM> & sip,
                                    FlatMatrixFixWidth<DIM_CURL> curlshape) const;

  virtual Vec <DIM_CURL_TRAIT<ET_trait<ET>::DIM>::DIM>
  EvaluateCurlShape (const IntegrationPoint & ip, 
                     FlatVector<double> x,
                     LocalHeap & lh) const;
};

 
/// 
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

///
template <>
class HCurlHighOrderFE<ET_TRIG> : public T_HCurlHighOrderFiniteElement<ET_TRIG>
{
private:
  typedef TrigShapesInnerLegendre T_INNERSHAPES; 
  typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
public:
  HCurlHighOrderFE () { ; }
  HCurlHighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[2], TFA & shape) const; 
};

///
template <>
class HCurlHighOrderFE<ET_QUAD> : public T_HCurlHighOrderFiniteElement<ET_QUAD>
{
private: 
   typedef VertexExtensionOptimal<3> T_VERTEXSHAPES;
public:
  HCurlHighOrderFE () { ; }
  HCurlHighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[2], TFA & shape) const; 
};

///
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

///
template <>
class HCurlHighOrderFE<ET_HEX> : public T_HCurlHighOrderFiniteElement<ET_HEX>
{
public:
  HCurlHighOrderFE () { ; }
  HCurlHighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[3], TFA & shape) const; 
};

///
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

///
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

