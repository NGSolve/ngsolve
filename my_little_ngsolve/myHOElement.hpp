#ifndef FILE_MYHOELEMENT_HPP
#define FILE_MYHOELEMENT_HPP

/*********************************************************************/
/* File:   myHOElement.hpp                                           */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*
  
It is also simple to implement high order elements ...

*/


namespace ngfem
{


  /*
    A Segment finite element with arbitrary order hierarchic basis
    functions
   */
  class MyHighOrderSegm : public ScalarFiniteElement<1>
  {
    int vnums[2];
  public:
    MyHighOrderSegm (int order);
    void SetVertexNumber (int i, int v) { vnums[i] = v; }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            FlatVector<> shape) const;
  
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             FlatMatrixFixWidth<1> dshape) const;

  private:
    template <class T>
    void T_CalcShape (const T & x, FlatArray<T> & shape) const;
  };


  /*
    A triangular finite element with arbitrary order hierarchic basis
    functions
   */
  class MyHighOrderTrig : public ScalarFiniteElement<2>
  {
    int vnums[3];
  public:
    MyHighOrderTrig (int order);
    void SetVertexNumber (int i, int v) { vnums[i] = v; }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            FlatVector<> shape) const;
  
    virtual void CalcDShape (const IntegrationPoint & ip, 
                             FlatMatrixFixWidth<2> dshape) const;

  private:
    template <class T>
    void T_CalcShape (const T & x, const T & y, FlatArray<T> & shape) const;
  };
}

#endif

