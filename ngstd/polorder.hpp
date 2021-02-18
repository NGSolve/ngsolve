#ifndef FILE_POLORDER
#define FILE_POLORDER

/**************************************************************************/
/* File:   polorder.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   10. July. 13                                                   */
/**************************************************************************/

namespace ngstd
{

  // Automatic differentiation datatype

  
  template <int DIM>
  class PolOrder
  {
    INT<DIM> order;

  public:
    /*
    typedef int TELEM;
    typedef int TSCAL;
    typedef int TV_COL;
    typedef int TV_ROW;
    enum { HEIGHT = 1 };
    enum { WIDTH = 1 };
    enum { IS_COMPLEX = 0 };
    */
  public:
    PolOrder (INT<DIM> ao) : order(ao) { ; }
    PolOrder (double /* d */) : order(0) { ; }
    PolOrder (int /* d */) : order(0) { ; }
    PolOrder () : order(0) { ; }

    PolOrder & operator= (INT<DIM> ao) { order = ao; return *this; }
    PolOrder & operator= (int /* i */) { order = 0; return *this; }
    PolOrder & operator= (double /* d */) { order = 0; return *this; }

    int operator() (int i) const { return order[i]; }
    int & operator() (int i) { return order[i]; }
    INT<DIM> GetOrder() const { return order; }
  };

  template <int DIM>  
  inline PolOrder<DIM> operator+ (PolOrder<DIM> a, PolOrder<DIM> b)
  {
    PolOrder<DIM> result;
    for (int i = 0; i < DIM; i++) 
      result(i) = max2(a(i), b(i));
    return result;
  }

  template <int DIM>  
  inline PolOrder<DIM> operator- (PolOrder<DIM> a, PolOrder<DIM> b)  { return a+b; }
  template <int DIM>  
  inline PolOrder<DIM> operator* (PolOrder<DIM> a, PolOrder<DIM> b)
  {
    PolOrder<DIM> result;
    for (int i = 0; i < DIM; i++) 
      result(i) = a(i)+b(i);
    return result;
  }
  template <int DIM>  
  inline PolOrder<DIM> operator- (PolOrder<DIM> a)  { return a; }

  // currently the best opinion ...
  template <int DIM>  
  inline PolOrder<DIM> operator/ (PolOrder<DIM> a, PolOrder<DIM> b) { return a*b; }


  template <int DIM>  
  inline PolOrder<DIM> operator- (double /* a */, PolOrder<DIM> b)  { return b; }
  template <int DIM>  
  inline PolOrder<DIM> operator+ (double /* a */, PolOrder<DIM> b)  { return b; }
  template <int DIM>  
  inline PolOrder<DIM> operator* (double /* a */, PolOrder<DIM> b)  { return b; }

  template <int DIM>  
  inline PolOrder<DIM> operator- (PolOrder<DIM> b, double /* a */)  { return b; }
  template <int DIM>  
  inline PolOrder<DIM> operator+ (PolOrder<DIM> b, double /* a */)  { return b; }
  template <int DIM>  
  inline PolOrder<DIM> operator* (PolOrder<DIM> b, double /* a */)  { return b; }
  template <int DIM>  
  inline PolOrder<DIM> operator/ (PolOrder<DIM> b, double /* a */)  { return b; }



  template <int DIM, typename T>  
  inline PolOrder<DIM> operator+= (PolOrder<DIM> & a, const T & b)  { a = a + b; return a; }
  template <int DIM, typename T>  
  inline PolOrder<DIM> operator-= (PolOrder<DIM> & a, const T & b)  { a = a - b; return a; }
  template <int DIM, typename T>  
  inline PolOrder<DIM> operator*= (PolOrder<DIM> & a, const T & b)  { a = a * b; return a; }
  template <int DIM, typename T>  
  inline PolOrder<DIM> operator/= (PolOrder<DIM> & a, const T & b)  { a = a / b; return a; }

  template <int DIM>
  inline bool operator== (PolOrder<DIM> /* a */, double /* b */) { return false; }

  template<int DIM>
  inline ostream & operator<< (ostream & ost, const PolOrder<DIM> & x)
  {
    ost << x.GetOrder();
    return ost;
  }
}


#endif
