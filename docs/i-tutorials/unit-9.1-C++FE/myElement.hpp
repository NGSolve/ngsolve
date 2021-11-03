#ifndef FILE_MYELEMENT_HPP
#define FILE_MYELEMENT_HPP


/*
  
My own simple first and second order triangular finite elements

*/


namespace ngfem
{

  /* 
     Our finite element base class. 
     The elements have to implement shape functions and their derivatives.
  */

  class MyBaseElement : public FiniteElement
  {
  public:
    /*
      Initialize with number of degrees of freedom,
      and polynomial order of basis functions.
    */
    MyBaseElement(int ndof, int order) : FiniteElement(ndof, order) {}

    /*
      Calculate the vector of shape functions in the point ip.
      ip is given on the reference element.
     */
    virtual void CalcShape(const IntegrationPoint & ip,
                           SliceVector<> shape) const = 0;
                           
    /*
      Calculate the matrix of derivatives of the shape functions in the point ip.
      dshape is a matrix of height ndof, and width of space dimension.
     */
    virtual void CalcDShape(const IntegrationPoint & ip,
                            SliceMatrix<> dshape) const = 0;
  };

  /*
    A linear triangular finite element.
  */
  class MyLinearTrig : public MyBaseElement
  {
  public:
    MyLinearTrig () : MyBaseElement(3, 1) {}
    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    void CalcShape (const IntegrationPoint & ip, 
                    SliceVector<> shape) const override;
  
    void CalcDShape (const IntegrationPoint & ip, 
                     SliceMatrix<> dshape) const override;
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
                    SliceVector<> shape) const override;
  
    void CalcDShape (const IntegrationPoint & ip, 
                     SliceMatrix<> dshape) const override;
  };

}

#endif // FILE_MYELEMENT_HPP

