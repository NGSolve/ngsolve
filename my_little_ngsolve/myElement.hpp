#ifndef FILE_MYELEMENT_HPP
#define FILE_MYELEMENT_HPP

/*********************************************************************/
/* File:   myElement.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*
  
My own simple first and second order triangular finite elements

*/


namespace ngfem
{

  /*
    A linear triangular finite element is a scalar element in two dimension, 
    thus we derive it from the ScalarFiniteElement<2> base class:
   */

  class MyLinearTrig : public ScalarFiniteElement<2>
  {
  public:
    // constructor
    MyLinearTrig ();

    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }
    
    /*
      Calculate the vector of shape functions in the point ip.
      ip is given in the reference element.
     */
    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceVector<> shape) const;
  
    /*
      Calculate the matrix of derivatives of the shape functions in the point ip.
      dshape is an 3 by 2 matrix in our case.
     */
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             SliceMatrix<> dshape) const;
  };


  /*
    A triangular finite element with second order basis functions
   */
  class MyQuadraticTrig : public ScalarFiniteElement<2>
  {
  public:
    MyQuadraticTrig ();
    virtual ELEMENT_TYPE ElementType() const { return ET_TRIG; }
    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceVector<> shape) const;
  
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             SliceMatrix<> dshape) const;
  };

}

#endif

