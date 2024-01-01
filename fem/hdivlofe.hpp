#ifndef FILE_HDIVLOFE
#define FILE_HDIVLOFE

/*********************************************************************/
/* File:   hdivlofe.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   5. Jul. 2001                                              */
/*********************************************************************/


#include "hdivfe.hpp"

namespace ngfem
{

  template <ELEMENT_TYPE ET>
  class HDivDummyFE : public HDivFiniteElement<ET_trait<ET>::DIM>
  {
  public:
    HDivDummyFE() : HDivFiniteElement<ET_trait<ET>::DIM> (0,0) { ; }
    virtual ELEMENT_TYPE ElementType() const override { return ET; }
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceMatrix<> shape) const override { ; }
  };

  template <ELEMENT_TYPE ET>
  class HDivNormalDummyFE : public HDivNormalFiniteElement<ET_trait<ET>::DIM>
  {
  public:
    HDivNormalDummyFE() : HDivNormalFiniteElement<ET_trait<ET>::DIM> (0,0) { ; }
    virtual ELEMENT_TYPE ElementType() const override { return ET; }
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const override { ; }
  };


  ///
  class FE_RTTrig0 : public HDivFiniteElement<2>
  {
    ///
  public:
    ///
    FE_RTTrig0();
    ///
    virtual ~FE_RTTrig0();
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceMatrix<> shape) const override;
  };



  /// RT0 + Curl B3
  class FE_RTTrig0plus : public HDivFiniteElement<2>
  {
    ///
  public:
    ///
    FE_RTTrig0plus();
    ///
    virtual ~FE_RTTrig0plus();
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceMatrix<> shape) const override;
  };




  ///
  class FE_BDMTrig1 : public HDivFiniteElement<2>
  {
    ///
    static Matrix<> trans;

  public:
    ///
    FE_BDMTrig1();
    ///
    virtual ~FE_BDMTrig1();
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceMatrix<> shape) const override;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<2> shape) const override;
  
    ///
    void Orthogonalize();
  };









  class HDivNormalSegm0 : public HDivNormalFiniteElement<1>
  {
  public:

    HDivNormalSegm0 ();
    virtual ELEMENT_TYPE ElementType() const override { return ET_SEGM; }
    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip,
			    FlatVector<> shape) const override;
  };










#ifdef ABC

  /// BDM1 + H(div)-bubble  curl B3
  class FE_BDMTrig1plus : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;
    ///
    static FlatMatrix<> trans;

  public:

    ///
    FE_BDMTrig1plus();
    ///
    virtual ~FE_BDMTrig1plus();

    ///
    virtual int SpatialDim () const { return 2; }
    ///
    virtual int GetNDof () const override { return 7; }
    ///
    virtual int Order () const override { return 2; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatVector<> & shape,
			     int comp = 1) const override;

    ///
    void Orthogonalize();

    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }
  };






  ///
  class FE_BDFMTrig2 : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;
    ///
    static FlatMatrix<> trans;

  public:

    ///
    FE_BDFMTrig2();
    ///
    virtual ~FE_BDFMTrig2();

    ///
    virtual int SpatialDim () const { return 2; }
    ///
    virtual int GetNDof () const override { return 9; }
    ///
    virtual int Order () const override{ return 2; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;

    /// full P2
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatVector<> & shape,
			     int comp = 1) const override;

    ///
    void Orthogonalize();

    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }
  };








  ///
  class FE_BDMTrig2 : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;
    ///
    static FlatMatrix<> trans;

  public:

    ///
    FE_BDMTrig2();
    ///
    virtual ~FE_BDMTrig2();

    ///
    virtual int SpatialDim () const { return 2; }
    ///
    virtual int GetNDof () const override { return 12; }
    ///
    virtual int Order () const override { return 3; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }
  
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatVector<> & shape,
			     int comp = 1) const override;

    ///
    void Orthogonalize();

    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }

  };







  ///
  class FE_BDMTrig2plus : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;
    ///
    static FlatMatrix<> trans;

  public:

    ///
    FE_BDMTrig2plus();
    ///
    virtual ~FE_BDMTrig2plus();

    ///
    virtual int SpatialDim () const { return 2; }
    ///
    virtual int GetNDof () const override { return 14; }
    ///
    virtual int Order () const override { return 2; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatVector<> & shape,
			     int comp = 1) const override;

    ///
    void Orthogonalize();

    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }

  };


#endif

  ///
  class FE_RTQuad0 : public HDivFiniteElement<2>
  {
  protected:
    ///
    // static Array<IPData> ipdata;
  
  public:
  
    ///
    FE_RTQuad0();
    ///
    virtual ~FE_RTQuad0();
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_QUAD; }

    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceMatrix<> shape) const override;

  };



#ifdef OLD


  ///
  class FE_BDMQuad1 : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;
    ///
    static FlatMatrix<> trans;

  public:

    ///
    FE_BDMQuad1();
    ///
    virtual ~FE_BDMQuad1();

    ///
    virtual int SpatialDim () const { return 2; }
    ///
    virtual int GetNDof () const override { return 8; }
    ///
    virtual int Order () const override { return 2; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_QUAD; }
  
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatVector<> & shape,
			     int comp = 1) const override;

    ///
    void Orthogonalize();

    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }

  };








  ///
  class FE_RTSegm0 : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;

  public:

    ///
    FE_RTSegm0();
    ///
    virtual ~FE_RTSegm0();

    ///
    virtual int SpatialDim () const { return 1; }
    ///
    virtual int GetNDof () const override { return 1; }
    ///
    virtual int Order () const override { return 1; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_SEGM; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;

    ///
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrix<> & dshape,
			     int comp = 1) const;

    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }

  };









  ///
  class FE_RTSegm1 : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;

  public:

    ///
    FE_RTSegm1();
    ///
    virtual ~FE_RTSegm1();

    ///
    virtual int SpatialDim () const { return 1; }
    ///
    virtual int GetNDof () const override { return 2; }
    ///
    virtual int Order () const override { return 1; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_SEGM; }


    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;

    ///
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrix<> & dshape,
			     int comp = 1) const;


    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }
  };




  ///
  class FE_RTSegm2 : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;

  public:

    ///
    FE_RTSegm2();
    ///
    virtual ~FE_RTSegm2();

    ///
    virtual int SpatialDim () const { return 1; }
    ///
    virtual int GetNDof () const override { return 3; }
    ///
    virtual int Order () const override { return 3; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_SEGM; }


    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;

    ///
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrix<> & dshape,
			     int comp = 1) const;


    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }
  };


#endif




  ///
  class FE_BDMTet1 : public HDivFiniteElement<3>
  {
    ///
// static Array<IPData> ipdata;
    ///
    static Matrix<> trans;

  public:

    ///
    FE_BDMTet1();
    ///
    virtual ~FE_BDMTet1();

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceMatrix<> shape) const override;

    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> shape) const override;

    ///
    void Orthogonalize();

  };



#ifdef ABC

  ///
  class FE_BDFMTet2 : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;
    ///
    static FlatMatrix<> trans;
    ///
    static FlatMatrix<> trans2;


  public:

    ///
    FE_BDFMTet2();
    ///
    virtual ~FE_BDFMTet2();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const override { return 18; }
    ///
    virtual int Order () const override { return 2; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_TET; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;

    /// full p2
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatVector<> & shape,
			     int comp = 1) const override;

    /// full P1
    virtual void CalcShape2 (const IntegrationPoint & ip, 
			     FlatVector<> & shape,
			     int comp = 1) const override;


    ///
    void Orthogonalize();

    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }
  };










  /**
     Space: Q1,1 + z-bubbles
  */
  class FE_BDMPrism1 : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;
    ///
    static FlatMatrix<> trans;

  public:

    ///
    FE_BDMPrism1();
    ///
    virtual ~FE_BDMPrism1();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const override { return 18; }
    ///
    virtual int Order () const override { return 1; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_PRISM; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;
  
    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatVector<> & shape,
			     int comp = 1) const override;
  
    ///
    void Orthogonalize();

    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }
  };




  /**
     Space: Q1,1 + z-bubbles
  */
  class FE_BDMPrism1p : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;
    ///
    static FlatMatrix<> trans;

  public:

    ///
    FE_BDMPrism1p();
    ///
    virtual ~FE_BDMPrism1p();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const override { return 21; }
    ///
    virtual int Order () const override { return 1; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_PRISM; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;
  
    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatVector<> & shape,
			     int comp = 1) const override;
  
    ///
    void Orthogonalize();

    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }
  };






  /**
     Space: BDFMTrig2 x P1  +  P1 x P2
     total: 2*9+3*3 = 27 dofs
     extern:
     3*4 quad dofs
     2*3 trig dofs
  */
  class FE_BDFMPrism2 : public HDivFiniteElement
  {
    ///
// static Array<IPData*> ipdata;
    ///
    static FlatMatrix<> trans;
    ///
    static FlatMatrix<> trans2;

  public:

    ///
    FE_BDFMPrism2();
    ///
    virtual ~FE_BDFMPrism2();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const override { return 27; }
    ///
    virtual int Order () const override { return 3; }
    ///
    virtual ELEMENT_TYPE ElementType() const override { return ET_PRISM; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> & shape,
			    int comp = 1) const override;
  
    ///
    virtual void CalcShape1 (const IntegrationPoint & ip, 
			     FlatVector<> & shape,
			     int comp = 1) const override;
  
    ///
    void Orthogonalize();

    ///
    // virtual const Array<IPData*> & GetIPData () const 
    // { return ipdata; }
  };



#endif

  

}


#endif
