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
     Our finite element base classes. The elements have to implement
     the methods CalcShape and CalcDShape to compute shape functions
     and gradients.
  */

  class MyBaseElement : public FiniteElement
  {
  public:
    MyBaseElement(int ndof, int order) : FiniteElement(ndof, order) {}

    virtual void CalcShape(const IntegrationPoint& ip,
                           BareSliceVector<> shape) const = 0;
    virtual void CalcDShape(const IntegrationPoint& ip,
                            BareSliceMatrix<> dshape) const = 0;
  };

  /*
    A linear triangular finite element.
  */
  class MyLinearTrig : public MyBaseElement
  {
  public:
    MyLinearTrig () : MyBaseElement(3, 1) {}
    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    /*
      Calculate the vector of shape functions in the point ip.
      ip is given in the reference element.
     */
    void CalcShape (const IntegrationPoint & ip, 
                    BareSliceVector<> shape) const override;
  
    /*
      Calculate the matrix of derivatives of the shape functions in the point ip.
      dshape is an 3 by 2 matrix in our case.
     */
    void CalcDShape (const IntegrationPoint & ip, 
                     BareSliceMatrix<> dshape) const override;
  };


  /*
    A triangular finite element with second order basis functions
   */
  class MyQuadraticTrig : public MyBaseElement
  {
  public:
    MyQuadraticTrig () : MyBaseElement(6,2) {}
    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    void CalcShape (const IntegrationPoint & ip, 
                    BareSliceVector<> shape) const override;
  
    void CalcDShape (const IntegrationPoint & ip, 
                     BareSliceMatrix<> dshape) const override;
  };

}

#endif // FILE_MYELEMENT_HPP

