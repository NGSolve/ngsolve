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
class H1HighOrderFiniteElement : virtual public ScalarFiniteElement<DIM>
{
public:
  int vnums[8];
  INT<3> order_inner;
  INT<2> order_face[6];
  int order_edge[12];

  using ScalarFiniteElement<DIM>::eltype;

public:
  
  void SetVertexNumbers (FlatArray<int> & avnums);
  virtual void GetDofs (Array<Dof> & dofs) const;

  void SetOrderInner (int oi);
  void SetOrderInner (INT<3> oi);
  void SetOrderFace (FlatArray<int> & of);
  void SetOrderFace (FlatArray<INT<2> > & of);
  void SetOrderEdge (FlatArray<int> & oe);

  /// high order elements need extra configuration. update ndof and order
  virtual void ComputeNDof () = 0;
  // const int * GetVNums() const { return vnums; }
};


template <ELEMENT_TYPE ET> class H1HighOrderFE;


template <ELEMENT_TYPE ET>
class T_H1HighOrderFiniteElement : public H1HighOrderFiniteElement<ET_trait<ET>::DIM>
{
  enum { DIM = ET_trait<ET>::DIM };

  using ScalarFiniteElement<DIM>::ndof;
  using ScalarFiniteElement<DIM>::order;
  using ScalarFiniteElement<DIM>::eltype;
  using ScalarFiniteElement<DIM>::dimspace;

  using H1HighOrderFiniteElement<DIM>::vnums;
  using H1HighOrderFiniteElement<DIM>::order_edge;
  using H1HighOrderFiniteElement<DIM>::order_face;
  using H1HighOrderFiniteElement<DIM>::order_inner;

public:

  T_H1HighOrderFiniteElement () 
  {
    for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
      vnums[i] = i;
    dimspace = DIM;
    eltype = ET;
  }

  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  virtual void CalcShape (const IntegrationPoint & ip, 
                          FlatVector<> shape) const;

  virtual void CalcDShape (const IntegrationPoint & ip, 
                           FlatMatrix<> dshape) const;

  /*
  /// compute dshape, matrix: ndof x spacedim
  virtual void CalcMappedDShape (const BaseSpecificIntegrationPoint & sip, 
                                 FlatMatrix<> dshape) const;

  {
    const SpecificIntegrationPoint<3,3> & sip = 
      static_cast<const SpecificIntegrationPoint<3,3> &> (bsip);
    
    const IntegrationPoint & ip = sip.IP();
    AutoDiff<3> x(ip(0)), y(ip(1)), z(ip(2));
    for (int j = 0; j < 3; j++)
      {
        x.DValue(j) = sip.GetJacobianInverse()(0,j);
        y.DValue(j) = sip.GetJacobianInverse()(1,j);
        z.DValue(j) = sip.GetJacobianInverse()(2,j);
      }

    ArrayMem<AutoDiff<3>,40> sds(ndof);
    T_CalcShape<AutoDiff<3>, AutoDiff<3>, AutoDiff<3>, FlatArray<AutoDiff<3> > >
      (x,y,z,sds);
    for (int i = 0; i < ndof; i++)
      for (int j = 0; j < 3; j++)
	dshape(i, j) = sds[i].DValue(j);
  }

  */

};








/**
  High order segment finite element
*/

template <> 
class H1HighOrderFE<ET_SEGM> : public T_H1HighOrderFiniteElement<ET_SEGM>
{
private:
  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
public:
  H1HighOrderFE () { ; }
  H1HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx x[1], TFA & shape) const; 
};


/**
  High order triangular finite element
*/
template <>
class H1HighOrderFE<ET_TRIG> : public T_H1HighOrderFiniteElement<ET_TRIG>
{
  typedef TrigShapesInnerLegendre T_INNERSHAPES;
  // typedef TrigShapesInnerJacobi T_INNERSHAPES;

  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
public:
  H1HighOrderFE () { ; }
  H1HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx x[2], TFA & shape) const; 
};


/**
  High order quadrilateral finite element
*/
template <>
class H1HighOrderFE<ET_QUAD> : public T_H1HighOrderFiniteElement<ET_QUAD>
{
  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
public:
  H1HighOrderFE () { ; }
  H1HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx x[2], TFA & shape) const; 
};


/**
  High order tetrahedral finite element
*/
template <>
class H1HighOrderFE<ET_TET> : public T_H1HighOrderFiniteElement<ET_TET>
{
  typedef TetShapesInnerLegendre T_INNERSHAPES;
  typedef TetShapesFaceLegendre T_FACESHAPES;
  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

  // typedef TetShapesInnerJacobi T_INNERSHAPES;
  // typedef TetShapesFaceJacobi T_FACESHAPES;

  // typedef TetShapesFaceOpt1 T_FACESHAPES;

public:
  H1HighOrderFE () { ; }
  H1HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[3], TFA & shape) const; 
};


/** 
  High order prismatic finite element
*/
template <>
class H1HighOrderFE<ET_PRISM> : public T_H1HighOrderFiniteElement<ET_PRISM>
{
  typedef TrigShapesInnerLegendre T_TRIGFACESHAPES;
  //typedef TrigShapesInnerJacobi T_TRIGFACESHAPES;
  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

public:
  H1HighOrderFE () { ; }
  H1HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[3], TFA & shape) const; 
};



/**
  High order hexahedral finite element
*/
template <> 
class H1HighOrderFE<ET_HEX> : public T_H1HighOrderFiniteElement<ET_HEX>
{
  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;
public:
  H1HighOrderFE () { ; }
  H1HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[3], TFA & shape) const; 
};


/**
  High order pyramid finite element
*/
template<>
class H1HighOrderFE<ET_PYRAMID> : public T_H1HighOrderFiniteElement<ET_PYRAMID>
{
  typedef TrigShapesInnerLegendre T_TRIGSHAPES;
  typedef IntegratedLegendreMonomialExt T_ORTHOPOL;

public:
  H1HighOrderFE () { ; }
  H1HighOrderFE (int aorder);
  virtual void ComputeNDof();
  virtual void GetInternalDofs (Array<int> & idofs) const;

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx hx[3], TFA & shape) const; 
};




#endif
