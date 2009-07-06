#ifndef FILE_SCALARFE
#define FILE_SCALARFE

/*********************************************************************/
/* File:   scalarfe.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

namespace ngfem
{


  /**
     Scalar finite element.
     Provides shape functions and derivaties.
  */
  template <int D>
  class ScalarFiniteElement : public FiniteElement
  {
  public:
    virtual string ClassName(void) const {return "ScalarFiniteElement";}
    ///
    ScalarFiniteElement () { dimspace = D; }
    ///
    ScalarFiniteElement (ELEMENT_TYPE aeltype, 
			 int andof = 0, int aorder = 0)
      : FiniteElement (D, aeltype, andof, aorder) { ; }
    ///
    virtual ~ScalarFiniteElement () { ; }
    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;

    /**
       returns shape functions in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatVector<> GetShape (const IntegrationPoint & ip, 
				 LocalHeap & lh) const
    {
      FlatVector<> shape(ndof, lh);
      CalcShape (ip, shape);
      return shape;
    }

    /**
       returns derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatMatrixFixWidth<D> 
    GetDShape (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      FlatMatrixFixWidth<D> dshape(ndof, lh);
      CalcDShape (ip, dshape);
      return dshape;
    }


    /// compute shape
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const = 0;
  
    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<D> dshape) const;

    /// compute dshape, matrix: ndof x spacedim
    virtual void CalcMappedDShape (const SpecificIntegrationPoint<D,D> & sip, 
				   FlatMatrixFixWidth<D> dshape) const;

    virtual double
    Evaluate (const IntegrationPoint & ip, 
	      FlatVector<double> x, LocalHeap & lh) const
    {
      return InnerProduct (GetShape(ip, lh), x);
    }  


    /**
       returns second derivatives in point ip.
       returns stored values for valid ip.IPNr(), else computes values
    */
    const FlatMatrix<> GetDDShape (const IntegrationPoint & ip, LocalHeap & lh) const
    {
      FlatMatrix<> ddshape(ndof, dimspace*dimspace, lh);
      CalcDDShape (ip, ddshape);
      return ddshape;
    }

    /// compute dshape, matrix: ndof x (spacedim spacedim)
    virtual void CalcDDShape (const IntegrationPoint & ip, 
			      FlatMatrix<> ddshape) const;



    virtual void EvaluateShapeGrid (const IntegrationRuleTP<D> & ir,
				    const FlatVector<double> coefs,
				    FlatVector<double> gridvalues,
				    LocalHeap & lh) const;
				  
    virtual void EvaluateShapeGridTrans (const IntegrationRuleTP<D> & ir,
					 const FlatVector<double> gridvalues,
					 FlatVector<double> coefs,
					 LocalHeap & lh) const;
				  
    virtual void EvaluateDShapeGrid (const IntegrationRuleTP<D> & ir,
				     const FlatVector<double> coefs,
				     FlatMatrixFixWidth<D> gridvalues,
				     LocalHeap & lh) const;
				  
    virtual void EvaluateDShapeGridTrans (const IntegrationRuleTP<D> & ir,
					  const FlatMatrixFixWidth<D> gridvalues,
					  FlatVector<double> coefs,
					  LocalHeap & lh) const;
  };











  /* ********************************** Segm ********************************* */




  /* ********************************* Trigs ******************************* */




  /* ***************************** Tet *************************************** */



  /* ***************************** Quads ********************************* */



  /* **************************** Pyramid Elements *********************** */

  /*
  /// pyramid of order 0
  class FE_Pyramid0 : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;

  public:

  ///
    FE_Pyramid0();
    ///
    virtual ~FE_Pyramid0();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 1; }
    ///
    virtual int Order () const { return 0; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_PYRAMID; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
    ///
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;

    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;

  };
  */


  /*
  /// pyramid of order 1
  class FE_Pyramid1 : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;

  public:

  ///
    FE_Pyramid1();
    ///
    virtual ~FE_Pyramid1();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 5; }
    ///
    virtual int Order () const { return 3; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_PYRAMID; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;
			  

    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;
  };
  */

  /*
  /// pyramid of order 2
  class FE_Pyramid2 : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;

  public:

  ///
    FE_Pyramid2();
    ///
    virtual ~FE_Pyramid2();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 13; }
    ///
    virtual int Order () const { return 4; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_PYRAMID; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
  };

  */



  /* ******************************* Prism Elements ********************* */


  /*
  /// prism of order 0
  class FE_Prism0 : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;

  public:

  ///
    FE_Prism0();
    ///
    virtual ~FE_Prism0();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;

    ///
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;

    ///
    virtual const IntegrationRule & NodalIntegrationRule() const;
  };
  */

  /*
  /// prism of order 1
  class FE_Prism1 : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;

  public:

  ///
    FE_Prism1();
    ///
    virtual ~FE_Prism1();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
  
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;
			  
    virtual const IntegrationRule & NodalIntegrationRule() const;
  };
  */



  /// prism of order 2
  class FE_Prism2 : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;

  public:

  ///
    FE_Prism2();
    ///
    virtual ~FE_Prism2();
    /*
   ///
   virtual int SpatialDim () const { return 3; }
   ///
   virtual int GetNDof () const { return 18; }
   ///
   virtual int Order () const { return 3; }
   ///
   virtual ELEMENT_TYPE ElementType() const { return ET_PRISM; }
    */
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;
  };





  /// in plane second order
  class FE_Prism2aniso : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;

  public:

  ///
    FE_Prism2aniso();
    ///
    virtual ~FE_Prism2aniso();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<3> dshape) const;
  };



  /// in plane second order
  class FE_Prism2HBaniso : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;

  public:

  ///
    FE_Prism2HBaniso();
    ///
    virtual ~FE_Prism2HBaniso();
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  

  };






  /// in plane third order
  class FE_Prism3aniso : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;

  public:

  ///
    FE_Prism3aniso();
    ///
    virtual ~FE_Prism3aniso();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 20; }
    ///
    virtual int Order () const { return 3; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_PRISM; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
  };





  /* **************************** Hex elements ************************* */


  /*
  ///
  class FE_Hex0 : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;
  public:
  ///
    FE_Hex0();
    ///
    virtual ~FE_Hex0();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 1; }
    ///
    virtual int Order () const { return 1; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_HEX; }

    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
  };
  */

  /*
  ///
  class FE_Hex1 : public ScalarFiniteElement<3>
  {
    ///
// static IPDataArray ipdata;
  public:
  ///
    FE_Hex1();
    ///
    virtual ~FE_Hex1();

    ///
    virtual int SpatialDim () const { return 3; }
    ///
    virtual int GetNDof () const { return 8; }
    ///
    virtual int Order () const { return 1; }
    ///
    virtual ELEMENT_TYPE ElementType() const { return ET_HEX; }
    ///
    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;
			  
  };
  */

}

#endif
