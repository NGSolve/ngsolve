#ifndef FILE_AUTODIFF
#define FILE_AUTODIFF

/**************************************************************************/
/* File:   autodiff.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   24. Oct. 02                                                    */
/**************************************************************************/

namespace ngstd
{

// Automatic differentiation datatype


/**
   Datatype for automatic differentiation.
   Contains function value and D derivatives. Algebraic
   operations are overloaded by using product-rule etc. etc. 
**/
template <int D, typename SCAL = double>
class AutoDiff
{
  SCAL val;
  SCAL dval[D?D:1];
public:

  typedef AutoDiff<D,SCAL> TELEM;
  typedef SCAL TSCAL;


  /// elements are undefined
  INLINE AutoDiff  () throw() { }; 
  // { val = 0; for (int i = 0; i < D; i++) dval[i] = 0; }  // !

  /// copy constructor
  INLINE AutoDiff  (const AutoDiff & ad2) throw()
  {
    val = ad2.val;
    for (int i = 0; i < D; i++)
      dval[i] = ad2.dval[i];
  }

  /// initial object with constant value
  INLINE AutoDiff  (SCAL aval) throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
  }

  /// init object with (val, e_diffindex)
  INLINE AutoDiff  (SCAL aval, int diffindex)  throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
    dval[diffindex] = 1;
  }

  INLINE AutoDiff (SCAL aval, const SCAL * grad)
  {
    val = aval;
    LoadGradient (grad);
  }

  /// assign constant value
  INLINE AutoDiff & operator= (SCAL aval) throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
    return *this;
  }

  /// returns value
  INLINE SCAL Value() const throw() { return val; }
  
  /// returns partial derivative
  INLINE SCAL DValue (int i) const throw() { return dval[i]; }

  ///
  INLINE void StoreGradient (SCAL * p) const 
  {
    for (int i = 0; i < D; i++)
      p[i] = dval[i];
  }

  INLINE void LoadGradient (const SCAL * p) 
  {
    for (int i = 0; i < D; i++)
      dval[i] = p[i];
  }

  /// access value
  INLINE SCAL & Value() throw() { return val; }

  /// accesses partial derivative 
  INLINE SCAL & DValue (int i) throw() { return dval[i]; }
};


//@{  AutoDiff helper functions.

/// prints AutoDiff
template<int D, typename SCAL>
inline ostream & operator<< (ostream & ost, const AutoDiff<D,SCAL> & x)
{
  ost << x.Value() << ", D = ";
  for (int i = 0; i < D; i++)
    ost << x.DValue(i) << " ";
  return ost;
}

/// AutoDiff plus AutoDiff
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator+ (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value () = x.Value()+y.Value();
  // AutoDiff<D,SCAL> res(x.Value()+y.Value());
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i) + y.DValue(i);
  return res;
}


/// AutoDiff minus AutoDiff
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator- (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x.Value()-y.Value();
  // AutoDiff<D,SCAL> res (x.Value()-y.Value());
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i) - y.DValue(i);
  return res;
}

/// double plus AutoDiff
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator+ (double x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = y.DValue(i);
  return res;
}

/// AutoDiff plus double
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator+ (const AutoDiff<D,SCAL> & y, double x) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = y.DValue(i);
  return res;
}


/// minus AutoDiff
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator- (const AutoDiff<D,SCAL> & x) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = -x.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = -x.DValue(i);
  return res;
}

/// AutoDiff minus double
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator- (const AutoDiff<D,SCAL> & x, double y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x.Value()-y;
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i);
  return res;
}

///
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator- (double x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x-y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = -y.DValue(i);
  return res;
}


/// double times AutoDiff
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator* (double x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x*y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x*y.DValue(i);
  return res;
}

/// AutoDiff times double
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator* (const AutoDiff<D,SCAL> & y, double x) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x*y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x*y.DValue(i);
  return res;
}

/// AutoDiff times AutoDiff
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator* (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  SCAL hx = x.Value();
  SCAL hy = y.Value();

  res.Value() = hx*hy;
  for (int i = 0; i < D; i++)
    res.DValue(i) = hx*y.DValue(i) + hy*x.DValue(i);

  return res;
}

/// AutoDiff times AutoDiff
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> sqr (const AutoDiff<D,SCAL> & x) throw()
{
  AutoDiff<D,SCAL> res;
  SCAL hx = x.Value();
  res.Value() = hx*hx;
  hx *= 2;
  for (int i = 0; i < D; i++)
    res.DValue(i) = hx*x.DValue(i);
  return res;
}

/// Inverse of AutoDiff
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> Inv (const AutoDiff<D,SCAL> & x)
{
  AutoDiff<D,SCAL> res(1.0 / x.Value());
  for (int i = 0; i < D; i++)
    res.DValue(i) = -x.DValue(i) / (x.Value() * x.Value());
  return res;
}


/// AutoDiff div AutoDiff
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator/ (const AutoDiff<D,SCAL> & x, const AutoDiff<D,SCAL> & y)
{
  return x * Inv (y);
}

/// AutoDiff div double
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> operator/ (const AutoDiff<D,SCAL> & x, double y)
{
  return (1/y) * x;
}

  /// double div AutoDiff
  template<int D, typename SCAL>
  INLINE AutoDiff<D,SCAL> operator/ (double x, const AutoDiff<D,SCAL> & y)
  {
    return x * Inv(y);
  }
  


  
  template <int D, typename SCAL, typename SCAL2>
  INLINE AutoDiff<D,SCAL> & operator+= (AutoDiff<D,SCAL> & x, SCAL2 y) throw()
  {
    x.Value() += y;
    return x;
  }


  /// 
  template <int D, typename SCAL>
  INLINE AutoDiff<D,SCAL> & operator+= (AutoDiff<D,SCAL> & x, AutoDiff<D,SCAL> y)
  {
    x.Value() += y.Value();
    for (int i = 0; i < D; i++)
      x.DValue(i) += y.DValue(i);
    return x;
  }

  ///
  template <int D, typename SCAL>
  INLINE AutoDiff<D,SCAL> & operator-= (AutoDiff<D,SCAL> & x, AutoDiff<D,SCAL> y)
  {
    x.Value() -= y.Value();
    for (int i = 0; i < D; i++)
      x.DValue(i) -= y.DValue(i);
    return x;

  }

  template <int D, typename SCAL, typename SCAL2>
  INLINE AutoDiff<D,SCAL> & operator-= (AutoDiff<D,SCAL> & x, SCAL2 y)
  {
    x.Value() -= y;
    return x;
  }

  ///
  template <int D, typename SCAL>
  INLINE AutoDiff<D,SCAL> & operator*= (AutoDiff<D,SCAL> & x, AutoDiff<D,SCAL> y) 
  {
    for (int i = 0; i < D; i++)
      x.DValue(i) = x.DValue(i)*y.Value() + x.Value() * y.DValue(i);
    x.Value() *= y.Value();
    return x;
  }

  ///
  template <int D, typename SCAL, typename SCAL2>
  INLINE AutoDiff<D,SCAL> & operator*= (AutoDiff<D,SCAL> & x, SCAL2 y) 
  {
    x.Value() *= y;
    for (int i = 0; i < D; i++)
      x.DValue(i) *= y;
    return x;
  }

  ///
  template <int D, typename SCAL>
  INLINE AutoDiff<D,SCAL> & operator/= (AutoDiff<D,SCAL> & x, SCAL y) 
  {
    SCAL iy = 1.0 / y;
    x.Value() *= iy;
    for (int i = 0; i < D; i++)
      x.DValue(i) *= iy;
    return x;
  }




  /// 
  template <int D, typename SCAL>
  INLINE bool operator== (AutoDiff<D,SCAL> x, SCAL val2) 
  {
    return x.Value() == val2;
  }

  ///
  template <int D, typename SCAL>
  INLINE bool operator!= (AutoDiff<D,SCAL> x, SCAL val2) throw()
  {
    return x.Value() != val2;
  }

  ///
  template <int D, typename SCAL>
  INLINE bool operator< (AutoDiff<D,SCAL> x, SCAL val2) throw()
  {
    return x.Value() < val2;
  }
  
  ///
  template <int D, typename SCAL>
  INLINE bool operator> (AutoDiff<D,SCAL> x, SCAL val2) throw()
  {
    return x.Value() > val2;
  }




template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> fabs (const AutoDiff<D,SCAL> & x)
{
  double abs = std::fabs (x.Value());
  AutoDiff<D,SCAL> res( abs );
  if (abs != 0.0)
    for (int i = 0; i < D; i++)
      res.DValue(i) = x.DValue(i) / abs;
  else
    for (int i = 0; i < D; i++)
      res.DValue(i) = 0.0;
  return res;
}

template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> sqrt (const AutoDiff<D,SCAL> & x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = std::sqrt(x.Value());
  for (int j = 0; j < D; j++)
    res.DValue(j) = 0.5 / res.Value() * x.DValue(j);
  return res;
}

using std::log;
template <int D, typename SCAL>
AutoDiff<D,SCAL> log (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = std::log(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) / x.Value();
  return res;
}



//@}

}

#endif


