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
template<int D>
class L2HighOrderFiniteElement : virtual public ScalarFiniteElement<D>
{
  enum { DIM = D };
public:
  int vnums[8];  
  INT<DIM> order_inner; 

  virtual void SetVertexNumbers (FlatArray<int> & avnums);

  void SetOrder (int o);
  void SetOrder (INT<DIM> oi);

  virtual void ComputeNDof () = 0; 
  virtual void GetInternalDofs (Array<int> & idofs) const; 
};



template <ELEMENT_TYPE ET> class L2HighOrderFE;

template <ELEMENT_TYPE ET>
class T_L2HighOrderFiniteElement : public L2HighOrderFiniteElement<ET_trait<ET>::DIM>
{
  enum { DIM = ET_trait<ET>::DIM };

  using ScalarFiniteElement<DIM>::ndof;
  using ScalarFiniteElement<DIM>::order;
  using ScalarFiniteElement<DIM>::eltype;
  using ScalarFiniteElement<DIM>::dimspace;

  using L2HighOrderFiniteElement<DIM>::vnums;
  using L2HighOrderFiniteElement<DIM>::order_inner;

public:

  T_L2HighOrderFiniteElement () 
  {
    for (int i = 0; i < ET_trait<ET>::N_VERTEX; i++)
      vnums[i] = i;
    dimspace = DIM;
    eltype = ET;
  }

  virtual void ComputeNDof();

  virtual void CalcShape (const IntegrationPoint & ip, 
                          FlatVector<> shape) const;

  virtual void CalcDShape (const IntegrationPoint & ip, 
                           FlatMatrixFixWidth<DIM> dshape) const;

};






/**
  High order 1D finite element
*/
template <>
class L2HighOrderFE<ET_SEGM> : public T_L2HighOrderFiniteElement<ET_SEGM>
{
public:
  L2HighOrderFE () { ; }
  L2HighOrderFE (int aorder);
  // virtual void ComputeNDof();

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx x[1], TFA & shape) const; 
};


/**
  High order triangular finite element
*/
template <> 
class L2HighOrderFE<ET_TRIG> : public T_L2HighOrderFiniteElement<ET_TRIG>
{
public:
  L2HighOrderFE () { ; }
  L2HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx x[2], TFA & shape) const; 
};

template <> 
class L2HighOrderFE<ET_QUAD> : public T_L2HighOrderFiniteElement<ET_QUAD>
{
public:
  L2HighOrderFE () { ; }
  L2HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx x[2], TFA & shape) const; 
};



template <> 
class L2HighOrderFE<ET_TET> : public T_L2HighOrderFiniteElement<ET_TET>
{
public:
  L2HighOrderFE () { ; }
  L2HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx x[3], TFA & shape) const; 
};


template <> 
class L2HighOrderFE<ET_PRISM> : public T_L2HighOrderFiniteElement<ET_PRISM>
{
public:
  L2HighOrderFE () { ; }
  L2HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx x[3], TFA & shape) const; 
};


template <> 
class L2HighOrderFE<ET_HEX> : public T_L2HighOrderFiniteElement<ET_HEX>
{
public:
  L2HighOrderFE () { ; }
  L2HighOrderFE (int aorder);

  template<typename Tx, typename TFA>  
  void T_CalcShape (Tx x[3], TFA & shape) const; 
};


#endif
