#ifndef FILE_HIGHORDERFE
#define FILE_HIGHORDERFE

/*********************************************************************/
/* File:   highorderfe.hpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   28. Oct. 2000                                             */
/*********************************************************************/

/*
  Old style high order Finite Elements
  nodal value shape functions (non - hierarchical)
*/


/**
   high order finite element
 */
class FE_SegmP : public ScalarFiniteElement<1>
{
  ///
  IPDataArray ipdata;
  ///
  int p;
  ///
  bool hierarchical;
public:
  ///
  FE_SegmP(int ap, bool ahierarchical = 1);
  ///
  virtual ~FE_SegmP();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;

  virtual const IPDataArray & GetIPData () const { return ipdata; }
};





///
class FE_TrigP : public ScalarFiniteElement<2>
{
  ///
  IPDataArray ipdata;
  ///
  int p;
  ///
  bool hierarchical;
public:
  ///
  FE_TrigP(int ap);
  ///
  virtual ~FE_TrigP();

  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  ngbla::FlatVector<> shape) const;

  virtual const IPDataArray & GetIPData () const { return ipdata; }
};








#ifdef NONE

///
class FE_QuadP : public ScalarFiniteElement<2>
{
  ///
  IPDataArray ipdata;
  ///
  FE_SegmP segm;
  ///
  FE_SegmP *segm2;
  ///
  int p, p2;
  ///
  bool hierarchical;
public:
  ///
  FE_QuadP(int ap, int ap2 = -1);
  ///
  virtual ~FE_QuadP();

  ///
  virtual int SpatialDim () const { return 2; }
  ///
  virtual int GetNDof () const;
  ///
  virtual int Order () const { return max2 (p, p2); }
  ///
  virtual ELEMENT_TYPE ElementType() const { return ET_QUAD; }
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  ngbla::FlatVector<> shape) const;
  /**
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   DenseMatrix & dshape,
			   int comp = 1) const;
  */
  virtual const IPDataArray & GetIPData () const { return ipdata; }
};
#endif




///
class FE_TetP : public ScalarFiniteElement<3>
{
  ///
  IPDataArray ipdata;
  ///
  int p;
  ///
  bool hierarchical;
public:
  ///
  FE_TetP(int ap);
  ///
  virtual ~FE_TetP();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  ngbla::FlatVector<> shape) const;

  ///
  virtual const IPDataArray & GetIPData () const { return ipdata; }
};


#ifdef NONE

///
class FE_PrismP : public ScalarFiniteElement<3>
{
  ///
  IPDataArray ipdata;
  ///
  int p;
public:
  ///
  FE_PrismP(int ap);
  ///
  virtual ~FE_PrismP();

  ///
  virtual int SpatialDim () const { return 3; }
  ///
  virtual int GetNDof () const;
  ///
  virtual int Order () const { return p+2; }
  ///
  virtual ELEMENT_TYPE ElementType() const { return ET_PRISM; }
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  ngbla::FlatVector<> shape) const;
  /**
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   DenseMatrix & dshape,
			   int comp = 1) const;
  */
  virtual const IPDataArray & GetIPData () const { return ipdata; }
};




///
class FE_HexP : public ScalarFiniteElement<3>
{
  ///
  IPDataArray ipdata;
  Array<IPData> ipdata;
  ///
  int p;
public:
  ///
  FE_HexP(int ap);
  ///
  virtual ~FE_HexP();

  ///
  virtual int SpatialDim () const { return 3; }
  ///
  virtual int GetNDof () const;
  ///
  virtual int Order () const { return p+1; }
  ///
  virtual ELEMENT_TYPE ElementType() const { return ET_HEX; }
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  ngbla::FlatVector<> shape) const;
  /**
  virtual void CalcDShape (const IntegrationPoint & ip, 
			   DenseMatrix & dshape,
			   int comp = 1) const;
  */
  virtual const IPDataArray & GetIPData () const { return ipdata; }
};
#endif








#ifdef ABC

///
class FE_Augmented_SegmP : public ScalarFiniteElement<1>
{
  ///
  Array<IPData> ipdata;
  ///
  int p;
public:
  ///
  FE_Augmented_SegmP(int ap);
  ///
  virtual ~FE_Augmented_SegmP();
  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;

  virtual const Array<IPData> & GetIPData () const { return ipdata; }
};



///
class FE_Augmented_TrigP : public ScalarFiniteElement<2>
{
  ///
  Array<IPData> ipdata;
  ///
  Matrix<> mat_unify;
  ///
  int p;
public:
  ///
  FE_Augmented_TrigP(int ap);
  ///
  virtual ~FE_Augmented_TrigP();

  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  FlatVector<> shape) const;

  ///
  virtual void Unify (FlatVector<> & uelem) const;
  ///
  virtual void UnifyTrans (FlatVector<> & uelem) const;

  virtual const Array<IPData> & GetIPData () const { return ipdata; }
};



///
class FE_Augmented_TetP : public ScalarFiniteElement<3>
{
  ///
  Array<IPData> ipdata;
  ///
  int p;
  ///
  Matrix<> mat_unify;
public:
  ///
  FE_Augmented_TetP(int ap);
  ///
  virtual ~FE_Augmented_TetP();

  ///
  virtual void CalcShape (const IntegrationPoint & ip, 
			  ngbla::FlatVector<> shape) const;

  ///
  virtual void Unify (FlatVector<> & uelem) const;
  ///
  virtual void UnifyTrans (FlatVector<> & uelem) const;

  ///
  virtual const Array<IPData> & GetIPData () const { return ipdata; }
};
#endif



#endif
