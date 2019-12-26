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

  template <int D, typename SCAL = double> class AutoDiffRec;


/**
   Datatype for automatic differentiation.
   Contains function value and D derivatives. Algebraic
   operations are overloaded by using product-rule etc. etc. 
**/
template <int D, typename SCAL = double>
class AutoDiff // : public AlignedAlloc<AutoDiff<D,SCAL>>
{
  SCAL val;
  SCAL dval[D?D:1];
public:

  typedef AutoDiff<D,SCAL> TELEM;
  typedef SCAL TSCAL;


  /// elements are undefined
  // INLINE AutoDiff  () throw() { };
  AutoDiff() = default;
  // { val = 0; for (int i = 0; i < D; i++) dval[i] = 0; }  // !

  /// copy constructor
  AutoDiff  (const AutoDiff & ad2) = default;
  /*
  INLINE AutoDiff  (const AutoDiff & ad2) throw()
  {
    val = ad2.val;
    for (int i = 0; i < D; i++)
      dval[i] = ad2.dval[i];
  }
  */
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

  AutoDiff & operator= (const AutoDiff & ad2) = default;
  
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
  template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
INLINE AutoDiff<D,SCAL> operator+ (SCAL2 x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = y.DValue(i);
  return res;
}

/// AutoDiff plus double
  template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
INLINE AutoDiff<D,SCAL> operator+ (const AutoDiff<D,SCAL> & y, SCAL2 x) throw()
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
  template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
INLINE AutoDiff<D,SCAL> operator- (const AutoDiff<D,SCAL> & x, SCAL2 y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x.Value()-y;
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i);
  return res;
}

///
  template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
INLINE AutoDiff<D,SCAL> operator- (SCAL2 x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x-y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = -y.DValue(i);
  return res;
}


/// double times AutoDiff
  template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
INLINE AutoDiff<D,SCAL> operator* (SCAL2 x, const AutoDiff<D,SCAL> & y) throw()
{
  AutoDiff<D,SCAL> res;
  res.Value() = x*y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x*y.DValue(i);
  return res;
}

/// AutoDiff times double
  template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>

  INLINE AutoDiff<D,SCAL> operator* (const AutoDiff<D,SCAL> & y, SCAL2 x) throw()
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
template<int D, typename SCAL, typename SCAL2,
         typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
INLINE AutoDiff<D,SCAL> operator/ (const AutoDiff<D,SCAL> & x, SCAL2 y)
{
  return (1.0/y) * x;
}

  /// double div AutoDiff
template<int D, typename SCAL, typename SCAL2,
         typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
  INLINE AutoDiff<D,SCAL> operator/ (SCAL2 x, const AutoDiff<D,SCAL> & y)
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
  double abs = fabs (x.Value());
  AutoDiff<D,SCAL> res( abs );
  if (abs != 0.0)
    for (int i = 0; i < D; i++)
      res.DValue(i) = x.Value()*x.DValue(i) / abs;
  else
    for (int i = 0; i < D; i++)
      res.DValue(i) = 0.0;
  return res;
}

using std::sqrt;
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> sqrt (const AutoDiff<D,SCAL> & x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = sqrt(x.Value());
  for (int j = 0; j < D; j++)
    res.DValue(j) = 0.5 / res.Value() * x.DValue(j);
  return res;
}

using std::log;
template <int D, typename SCAL>
AutoDiff<D,SCAL> log (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = log(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) / x.Value();
  return res;
}

using std::exp;
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> exp (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = exp(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * res.Value();
  return res;
}

using std::pow;
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> pow (AutoDiff<D,SCAL> x, AutoDiff<D,SCAL> y )
{
  return exp(log(x)*y);
}

using std::sin;
/*
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> sin (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = sin(x.Value());
  SCAL c = cos(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * c;
  return res;
}
*/

template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> sin (AutoDiff<D,SCAL> x)
{
  return sin(AutoDiffRec<D,SCAL>(x));
}
  
using std::cos;
/*
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> cos (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = cos(x.Value());
  SCAL ms = -sin(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * ms;
  return res;
}
*/
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> cos (AutoDiff<D,SCAL> x)
{
  return cos(AutoDiffRec<D,SCAL>(x));
}

using std::tan;
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> tan (AutoDiff<D,SCAL> x)
{ return sin(x) / cos(x); }

using std::sinh;
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> sinh (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = sinh(x.Value());
  SCAL ch = cosh(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * ch;
  return res;
}

using std::cosh;
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> cosh (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = cosh(x.Value());
  SCAL sh = sinh(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * sh;
  return res;
}

using std::floor;
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> floor (const AutoDiff<D,SCAL> & x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = floor(x.Value());
  for (int j = 0; j < D; j++)
    res.DValue(j) = 0.0;
  return res;
}

using std::ceil;
template<int D, typename SCAL>
INLINE AutoDiff<D,SCAL> ceil (const AutoDiff<D,SCAL> & x)
{
  AutoDiff<D,SCAL> res;
  res.Value() = ceil(x.Value());
  for (int j = 0; j < D; j++)
    res.DValue(j) = 0.0;
  return res;
}


using std::atan;
/*
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> atan (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  SCAL a = atan(x.Value());
  res.Value() = a;
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k)/(1+x.Value()*x.Value()) ;
  return res;
}
*/
template <int D, typename SCAL>
AutoDiff<D,SCAL> atan (AutoDiff<D,SCAL> x)
{
  return atan (AutoDiffRec<D,SCAL> (x));
}

using std::atan2;
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> atan2 (AutoDiff<D,SCAL> x, AutoDiff<D,SCAL> y)
{
  AutoDiff<D,SCAL> res;
  SCAL a = atan2(x.Value(), y.Value());
  res.Value() = a;
  for (int k = 0; k < D; k++)
    res.DValue(k) = (x.Value()*y.DValue(k)-y.Value()*x.DValue(k))/(y.Value()*y.Value()+x.Value()*x.Value());
  return res;
}


using std::acos;
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> acos (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  SCAL a = acos(x.Value());
  res.Value() = a;
  SCAL da = -1 / sqrt(1-x.Value()*x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k)*da;
  return res;
}


using std::asin;
template <int D, typename SCAL>
INLINE AutoDiff<D,SCAL> asin (AutoDiff<D,SCAL> x)
{
  AutoDiff<D,SCAL> res;
  SCAL a = asin(x.Value());
  res.Value() = a;
  SCAL da = 1 / sqrt(1-x.Value()*x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k)*da;
  return res;
}




  template <int D, typename SCAL, typename TB, typename TC>
  auto IfPos (AutoDiff<D,SCAL> a, TB b, TC c) // -> decltype(IfPos (a.Value(), b, c))
  {
    return IfPos (a.Value(), b, c);
  }

  template <int D, typename SCAL>
  INLINE AutoDiff<D,SCAL> IfPos (SCAL /* SIMD<double> */ a, AutoDiff<D,SCAL> b, AutoDiff<D,SCAL> c)
  {
    AutoDiff<D,SCAL> res;
    res.Value() = IfPos (a, b.Value(), c.Value());
    for (int j = 0; j < D; j++)
      res.DValue(j) = IfPos (a, b.DValue(j), c.DValue(j));
    return res;
  }

  template <int D, typename SCAL, typename TC>
  INLINE AutoDiff<D,SCAL> IfPos (SCAL /* SIMD<double> */ a, AutoDiff<D,SCAL> b, TC c)
  {
    return IfPos (a, b, AutoDiff<D,SCAL> (c));
  }

//@}


  
  template <int D, typename SCAL>
  class AutoDiffRec //  : public AlignedAlloc<AutoDiffRec<D,SCAL>>
  {
    AutoDiffRec<D-1, SCAL> rec;
    SCAL last;

  public:
    INLINE AutoDiffRec () = default;
    INLINE AutoDiffRec (const AutoDiffRec &) = default;
    INLINE AutoDiffRec (AutoDiffRec<D-1,SCAL> _rec, SCAL _last) : rec(_rec), last(_last) { ; }
    INLINE AutoDiffRec & operator= (const AutoDiffRec &) = default;

    INLINE AutoDiffRec (SCAL aval) : rec(aval), last(0.0) { ; }
    INLINE AutoDiffRec (SCAL aval, int diffindex) : rec(aval, diffindex), last((diffindex==D-1) ? 1.0 : 0.0) { ; }
    INLINE AutoDiffRec (const AutoDiff<D,SCAL> & ad)
    {
      Value() = ad.Value();
      for (int i = 0; i < D; i++)
        DValue(i) = ad.DValue(i);
    }
    
    INLINE AutoDiffRec & operator= (SCAL aval) { rec = aval; last = 0.0; return *this; }
    INLINE SCAL Value() const { return rec.Value(); }
    INLINE SCAL DValue(int i) const { return (i == D-1) ? last : rec.DValue(i); }
    INLINE SCAL & Value() { return rec.Value(); }
    INLINE SCAL & DValue(int i) { return (i == D-1) ? last : rec.DValue(i); }
    INLINE auto Rec() const { return rec; }
    INLINE auto Last() const { return last; }
    INLINE auto & Rec() { return rec; }
    INLINE auto & Last() { return last; }
    INLINE operator AutoDiff<D,SCAL> () const
    {
      AutoDiff<D,SCAL> res(Value());
      for (int i = 0; i < D; i++)
        res.DValue(i) = DValue(i);
      return res;
    }
  };

  template<int D, typename SCAL>
  ostream & operator<< (ostream & ost, AutoDiffRec<D,SCAL> ad)
  {
    return ost << AutoDiff<D,SCAL> (ad);
  }

  template <typename SCAL>
  class AutoDiffRec<0,SCAL> : public AlignedAlloc<AutoDiffRec<0,SCAL>>
  {
    SCAL val;
  public:
    INLINE AutoDiffRec () = default;
    INLINE AutoDiffRec (const AutoDiffRec &) = default;
    INLINE AutoDiffRec (SCAL _val) : val(_val) { ; }
    INLINE AutoDiffRec (SCAL _val, SCAL /* _dummylast */) : val(_val) { ; }

    INLINE AutoDiffRec & operator= (const AutoDiffRec &) = default;
    INLINE AutoDiffRec & operator= (SCAL aval) { val = aval; return *this; }

    INLINE SCAL Value() const { return val; }
    INLINE SCAL DValue(int i) const { return SCAL(0); }
    INLINE SCAL & Value() { return val; }
    // SCAL & DValue(int i) { return val; }
    INLINE auto Rec() const { return val; }
    INLINE auto Last() const { return SCAL(0); }
    INLINE auto & Rec() { return val; }
    INLINE auto & Last() { return val; }
    INLINE operator AutoDiff<0,SCAL> () const { return AutoDiff<0,SCAL>(); }
  };


  template <typename SCAL>
  class AutoDiffRec<1,SCAL> : public AlignedAlloc<AutoDiffRec<1,SCAL>>
  {
    SCAL val;
    SCAL last;
  public:
    INLINE AutoDiffRec () = default;
    INLINE AutoDiffRec (const AutoDiffRec &) = default;
    INLINE AutoDiffRec (SCAL _val) : val(_val), last(0.0) { ; }
    INLINE AutoDiffRec (SCAL _val, SCAL _last) : val(_val), last(_last) { ; }
    INLINE AutoDiffRec (SCAL aval, int diffindex) : val(aval), last((diffindex==0) ? 1.0 : 0.0) { ; }    
    INLINE AutoDiffRec (const AutoDiff<1,SCAL> & ad)
    {
      Value() = ad.Value();
      DValue(0) = ad.DValue(0);
    }

    INLINE AutoDiffRec & operator= (const AutoDiffRec &) = default;
    INLINE AutoDiffRec & operator= (SCAL aval) { val = aval; last = 0.0; return *this; }

    INLINE SCAL Value() const { return val; }
    INLINE SCAL DValue(int i) const { return last; }
    INLINE SCAL & Value() { return val; }
    INLINE SCAL & DValue(int i) { return last; }
    INLINE auto Rec() const { return val; }
    INLINE auto Last() const { return last; }
    INLINE auto & Rec() { return val; }
    INLINE auto & Last() { return last; }

    INLINE operator AutoDiff<1,SCAL> () const
    {
      AutoDiff<1,SCAL> res(Value());
      res.DValue(0) = DValue(0);
      return res;
    }
  };

  template <int D, typename SCAL>
  INLINE AutoDiffRec<D,SCAL> operator+ (double a, AutoDiffRec<D,SCAL> b)
  {
    return AutoDiffRec<D,SCAL> (a+b.Rec(), b.Last());
  }
  
  template <int D, typename SCAL>
  INLINE AutoDiffRec<D,SCAL> operator+ (AutoDiffRec<D,SCAL> a, double b)
  {
    return AutoDiffRec<D,SCAL> (a.Rec()+b, a.Last());
  }
  
  template <int D, typename SCAL>  
  INLINE AutoDiffRec<D,SCAL> operator+ (AutoDiffRec<D,SCAL> a, AutoDiffRec<D,SCAL> b)
  {
    return AutoDiffRec<D,SCAL> (a.Rec()+b.Rec(), a.Last()+b.Last());
  }

  template <int D, typename SCAL>
  INLINE AutoDiffRec<D,SCAL> operator- (double b, AutoDiffRec<D,SCAL> a)
  {
    return AutoDiffRec<D,SCAL> (b-a.Rec(), -a.Last());
  }
  
  template <int D, typename SCAL>
  INLINE AutoDiffRec<D,SCAL> operator- (AutoDiffRec<D,SCAL> a, double b)
  {
    return AutoDiffRec<D,SCAL> (a.Rec()-b, a.Last());
  }
  
  template <int D, typename SCAL>  
  INLINE AutoDiffRec<D,SCAL> operator- (AutoDiffRec<D,SCAL> a, AutoDiffRec<D,SCAL> b)
  {
    return AutoDiffRec<D,SCAL> (a.Rec()-b.Rec(), a.Last()-b.Last());
  }

  /// minus AutoDiff
  template<int D, typename SCAL>
  INLINE AutoDiffRec<D,SCAL> operator- (AutoDiffRec<D,SCAL> a)
  {
    return AutoDiffRec<D,SCAL> (-a.Rec(), -a.Last());
  }

  template <int D, typename SCAL>  
  INLINE AutoDiffRec<D,SCAL> operator* (AutoDiffRec<D,SCAL> a, AutoDiffRec<D,SCAL> b)
  {
    return AutoDiffRec<D,SCAL> (a.Rec()*b.Rec(), a.Value()*b.Last()+b.Value()*a.Last());
  }

  template <int D, typename SCAL, typename SCAL1>  
  INLINE AutoDiffRec<D,SCAL> operator* (AutoDiffRec<D,SCAL> b, SCAL1 a)
  {
    return AutoDiffRec<D,SCAL> (a*b.Rec(), a*b.Last());
  }

  template <int D, typename SCAL, typename SCAL1>  
  INLINE AutoDiffRec<D,SCAL> operator* (SCAL1 a, AutoDiffRec<D,SCAL> b)
  {
    return AutoDiffRec<D,SCAL> (a*b.Rec(), a*b.Last());
  }

  template <int D, typename SCAL>
  INLINE AutoDiffRec<D,SCAL> & operator+= (AutoDiffRec<D,SCAL> & a, AutoDiffRec<D,SCAL> b)
  {
    a.Rec() += b.Rec();
    a.Last() += b.Last();
    return a;
  }

  template <int D, typename SCAL>
  INLINE AutoDiffRec<D,SCAL> & operator-= (AutoDiffRec<D,SCAL> & a, double b)
  {
    a.Rec() -= b;
    return a;
  }

  template <int D, typename SCAL>
  INLINE AutoDiffRec<D,SCAL> & operator-= (AutoDiffRec<D,SCAL> & a, AutoDiffRec<D,SCAL> b)
  {
    a.Rec() -= b.Rec();
    a.Last() -= b.Last();
    return a;
  }


  template <int D, typename SCAL>  
  INLINE AutoDiffRec<D,SCAL> & operator*= (AutoDiffRec<D,SCAL> & a, AutoDiffRec<D,SCAL> b)
  {
    a = a*b;
    return a;
  }

  
  template <int D, typename SCAL, typename SCAL2>  
  INLINE AutoDiffRec<D,SCAL> & operator*= (AutoDiffRec<D,SCAL> & b, SCAL2 a)
  {
    b.Rec() *= a;
    b.Last() *= a;
    return b;
  }

  /// Inverse of AutoDiffRec

  template <typename SCAL>
  auto Inv1 (SCAL x)  { return 1.0/x; }

  template<int D, typename SCAL>
  INLINE AutoDiffRec<D,SCAL> Inv1 (AutoDiffRec<D,SCAL> x)
  {
    return AutoDiffRec<D,SCAL> (Inv1(x.Rec()), (-sqr(1.0/x.Value())) * x.Last());
  }

  /// AutoDiffRec div AutoDiffRec
  template<int D, typename SCAL>
  INLINE AutoDiffRec<D,SCAL> operator/ (const AutoDiffRec<D,SCAL> & x, const AutoDiffRec<D,SCAL> & y)
  {
    return x * Inv1 (y);
  }


  template <int D, typename SCAL>
  auto sin (AutoDiffRec<D,SCAL> x)
  {
    return AutoDiffRec<D,SCAL> (sin(x.Rec()), cos(x.Value())*x.Last());
  }

  template <int D, typename SCAL>
  auto cos (AutoDiffRec<D,SCAL> x)
  {
    return AutoDiffRec<D,SCAL> (cos(x.Rec()), -sin(x.Value())*x.Last());
  }

  template <int D, typename SCAL>
  auto tan (AutoDiffRec<D,SCAL> x)
  {
    return sin(x) / cos(x);
  }

  template <int D, typename SCAL>
  auto atan (AutoDiffRec<D,SCAL> x)
  {
    return AutoDiffRec<D,SCAL> (atan(x.Rec()), (1./(1.+x.Value()*x.Value()))*x.Last());
  }

  
}

#endif


