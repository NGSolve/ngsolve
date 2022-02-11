#ifndef FILE_AUTODIFFDIFF
#define FILE_AUTODIFFDIFF

/**************************************************************************/
/* File:   autodiffdiff.hpp                                               */
/* Author: Joachim Schoeberl                                              */
/* Date:   13. June. 05                                                   */
/**************************************************************************/

namespace ngstd
{
  using ngcore::IfPos;

// Automatic second differentiation datatype


/**
   Datatype for automatic differentiation.  Contains function value,
   D first derivatives, and D*D second derivatives. Algebraic operations are
   overloaded by using product-rule etc. etc.
**/
template <int D, typename SCAL = double>
class AutoDiffDiff
{
  SCAL val;
  SCAL dval[D?D:1];
  SCAL ddval[D?D*D:1];
public:

  typedef AutoDiffDiff<D, SCAL> TELEM;


  /// elements are undefined
  AutoDiffDiff  () throw() { ; }

  /// copy constructor
  AutoDiffDiff  (const AutoDiffDiff & ad2) throw()
  {
    val = ad2.val;
    for (int i = 0; i < D; i++)
      dval[i] = ad2.dval[i];
    for (int i = 0; i < D*D; i++)
      ddval[i] = ad2.ddval[i];
  }

  /// initial object with constant value
  AutoDiffDiff  (SCAL aval) throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
    for (int i = 0; i < D*D; i++)
      ddval[i] = 0;
  }

  /// initial object with value and derivative
  AutoDiffDiff  (const AutoDiff<D, SCAL> & ad2) throw()
  {
    val = ad2.Value();
    for (int i = 0; i < D; i++)
      dval[i] = ad2.DValue(i);
    for (int i = 0; i < D*D; i++)
      ddval[i] = 0;
  }

  /// init object with (val, e_diffindex)
  AutoDiffDiff  (SCAL aval, int diffindex)  throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
    for (int i = 0; i < D*D; i++)
      ddval[i] = 0;
    dval[diffindex] = 1;
  }

  INLINE AutoDiffDiff (SCAL aval, const SCAL * grad)
  {
    val = aval;
    LoadGradient (grad);
    for (int i = 0; i < D*D; i++)
      ddval[i] = 0;
  }

  INLINE AutoDiffDiff (SCAL aval, const SCAL * grad, const SCAL * hesse)
  {
    val = aval;
    LoadGradient (grad);
    LoadHessian (hesse);
  }

  /// assign constant value
  AutoDiffDiff & operator= (SCAL aval) throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
    for (int i = 0; i < D*D; i++)
      ddval[i] = 0;
    return *this;
  }

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

  INLINE void StoreHessian (SCAL * p) const 
  {
    for (int i = 0; i < D*D; i++)
      p[i] = ddval[i];
  }

  INLINE void LoadHessian (const SCAL * p) 
  {
    for (int i = 0; i < D*D; i++)
      ddval[i] = p[i];
  }

  /// returns value
  SCAL Value() const throw() { return val; }

  /// returns partial derivative
  SCAL DValue (int i) const throw() { return dval[i]; }

  AutoDiff<D,SCAL> DValueAD (int i) const
  {
    AutoDiff<D,SCAL> r(dval[i]);
    for (int j = 0; j < D; j++)
      r.DValue(j) = ddval[i*D+j];
    return r;
  }
  
  /// returns partial derivative
  SCAL DDValue (int i) const throw() { return ddval[i]; }

  /// returns partial derivative
  SCAL DDValue (int i, int j) const throw() { return ddval[i*D+j]; }

  /// access value
  SCAL & Value() throw() { return val; }

  /// accesses partial derivative 
  SCAL & DValue (int i) throw() { return dval[i]; }

  /// accesses partial derivative 
  SCAL & DDValue (int i) throw() { return ddval[i]; }

  /// accesses partial derivative 
  SCAL & DDValue (int i, int j) throw() { return ddval[i*D+j]; }

  explicit operator AutoDiff<D,SCAL> () const
  { return AutoDiff<D,SCAL> (val, &dval[0]); }
  
  /// add autodiffdiff object
  AutoDiffDiff<D, SCAL> & operator+= (const AutoDiffDiff<D, SCAL> & y) throw()
  {
    val += y.val;
    for (int i = 0; i < D; i++)
      dval[i] += y.dval[i];
    for (int i = 0; i < D*D; i++)
      ddval[i] += y.ddval[i];
    return *this;
  }

  /// subtract autodiffdiff object
  AutoDiffDiff<D, SCAL> & operator-= (const AutoDiffDiff<D, SCAL> & y) throw()
  {
    val -= y.val;
    for (int i = 0; i < D; i++)
      dval[i] -= y.dval[i];
    for (int i = 0; i < D*D; i++)
      ddval[i] -= y.ddval[i];
    return *this;
  }

  /// multiply with autodiffdiff object
  AutoDiffDiff<D, SCAL> & operator*= (const AutoDiffDiff<D, SCAL> & y) throw()
  {
    for (int i = 0; i < D*D; i++)
      ddval[i] = val * y.ddval[i] + y.val * ddval[i];

    for (int i = 0; i < D; i++)
      for (int j = 0; j < D; j++)
	ddval[i*D+j] += dval[i] * y.dval[j] + dval[j] * y.dval[i];

    for (int i = 0; i < D; i++)
      {
	dval[i] *= y.val;
	dval[i] += val * y.dval[i];
      }
    val *= y.val;
    return *this;
  }

  /// multiply with scalar
  AutoDiffDiff<D, SCAL> & operator*= (const SCAL & y) throw()
  {
    for ( int i = 0; i < D*D; i++ )
      ddval[i] *= y;
    for (int i = 0; i < D; i++)
      dval[i] *= y;
    val *= y;
    return *this;
  }

  /// divide by scalar
  AutoDiffDiff<D, SCAL> & operator/= (const SCAL & y) throw()
  {
    SCAL iy = 1.0 / y;
    for ( int i = 0; i < D*D; i++ )
      ddval[i] *= iy;
    for (int i = 0; i < D; i++)
      dval[i] *= iy;
    val *= iy;
    return *this;
  }

  /// same value
  bool operator== (SCAL val2) throw()
  {
    return val == val2;
  }

  /// different values 
  bool operator!= (SCAL val2) throw()
  {
    return val != val2;
  }

  /// less 
  bool operator< (SCAL val2) throw()
  {
    return val < val2;
  }
  
  /// greater
  bool operator> (SCAL val2) throw()
  {
    return val > val2;
  }
};


//@{  AutoDiff helper functions.

/// Prints AudoDiffDiff
template<int D, typename SCAL>
inline ostream & operator<< (ostream & ost, const AutoDiffDiff<D, SCAL> & x)
{
  ost << x.Value() << ", D = ";
  for (int i = 0; i < D; i++)
    ost << x.DValue(i) << " ";
  ost << ", DD = ";
  for (int i = 0; i < D*D; i++)
    ost << x.DDValue(i) << " ";
  return ost;
}

///
template<int D, typename SCAL>
inline AutoDiffDiff<D, SCAL> operator+ (const AutoDiffDiff<D, SCAL> & x, const AutoDiffDiff<D, SCAL> & y) throw()
{
  AutoDiffDiff<D, SCAL> res;
  res.Value () = x.Value()+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i) + y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = x.DDValue(i) + y.DDValue(i);
  return res;
}


/// 
template<int D, typename SCAL>
inline AutoDiffDiff<D, SCAL> operator- (const AutoDiffDiff<D, SCAL> & x, const AutoDiffDiff<D, SCAL> & y) throw()
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = x.Value()-y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i) - y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = x.DDValue(i) - y.DDValue(i);
  return res;
}


///
template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
inline AutoDiffDiff<D, SCAL> operator+ (SCAL2 x, const AutoDiffDiff<D, SCAL> & y) throw()
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = x+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = y.DDValue(i);
  return res;
}

///
template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
inline AutoDiffDiff<D, SCAL> operator+ (const AutoDiffDiff<D, SCAL> & y, SCAL2 x) throw()
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = x+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = y.DDValue(i);
  return res;
}


///
template<int D, typename SCAL>
inline AutoDiffDiff<D, SCAL> operator- (const AutoDiffDiff<D, SCAL> & x) throw()
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = -x.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = -x.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = -x.DDValue(i);
  return res;
}

///
template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
inline AutoDiffDiff<D, SCAL> operator- (const AutoDiffDiff<D, SCAL> & x, SCAL2 y) throw()
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = x.Value()-y;
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = x.DDValue(i);
  return res;
}

///
template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
inline AutoDiffDiff<D, SCAL> operator- (SCAL2 x, const AutoDiffDiff<D, SCAL> & y) throw()
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = x-y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = -y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = -y.DDValue(i);
  return res;
}


///
template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
inline AutoDiffDiff<D, SCAL> operator* (SCAL2 x, const AutoDiffDiff<D, SCAL> & y) throw()
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = x*y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x*y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = x*y.DDValue(i);
  return res;
}

///
template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
inline AutoDiffDiff<D, SCAL> operator* (const AutoDiffDiff<D, SCAL> & y, SCAL2 x) throw()
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = x*y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x*y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = x*y.DDValue(i);
  return res;
}

///
template<int D, typename SCAL>
inline AutoDiffDiff<D, SCAL> operator* (const AutoDiffDiff<D, SCAL> & x, const AutoDiffDiff<D, SCAL> & y) throw()
{
  AutoDiffDiff<D, SCAL> res;
  SCAL hx = x.Value();
  SCAL hy = y.Value();

  res.Value() = hx*hy;
  for (int i = 0; i < D; i++)
    res.DValue(i) = hx*y.DValue(i) + hy*x.DValue(i);

  for (int i = 0; i < D; i++)
    for (int j = 0; j < D; j++)
      res.DDValue(i,j) = hx * y.DDValue(i,j) + hy * x.DDValue(i,j)
	+ x.DValue(i) * y.DValue(j) + x.DValue(j) * y.DValue(i);

  return res;
}



template<int D, typename SCAL>
inline AutoDiffDiff<D, SCAL> Inv (const AutoDiffDiff<D, SCAL> & x)
{
  AutoDiffDiff<D, SCAL> res(1.0 / x.Value());
  for (int i = 0; i < D; i++)
    res.DValue(i) = -x.DValue(i) / (x.Value() * x.Value());

  SCAL fac1 = 2/(x.Value()*x.Value()*x.Value());
  SCAL fac2 = 1/sqr(x.Value());
  for (int i = 0; i < D; i++)
    for (int j = 0; j < D; j++)
      res.DDValue(i,j) = fac1*x.DValue(i)*x.DValue(j) - fac2*x.DDValue(i,j);
  return res;
}


template<int D, typename SCAL>
inline AutoDiffDiff<D, SCAL> operator/ (const AutoDiffDiff<D, SCAL> & x, const AutoDiffDiff<D, SCAL> & y)
{
  return x * Inv (y);
}

template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
inline AutoDiffDiff<D, SCAL> operator/ (const AutoDiffDiff<D, SCAL> & x, SCAL2 y)
{
  return (1/y) * x;
}

template<int D, typename SCAL, typename SCAL2,
           typename std::enable_if<std::is_convertible<SCAL2,SCAL>::value, int>::type = 0>
inline AutoDiffDiff<D, SCAL> operator/ (SCAL2 x, const AutoDiffDiff<D, SCAL> & y)
{
  return x * Inv(y);
}


template<int D, typename SCAL>
inline AutoDiffDiff<D, SCAL> sqrt (const AutoDiffDiff<D, SCAL> & x)
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = sqrt(x.Value());
  for (int j = 0; j < D; j++)
    res.DValue(j) = IfZero(x.DValue(j),SCAL{0.},0.5 / res.Value() * x.DValue(j));

  
  for (int i = 0; i < D; i++)
    for (int j = 0; j < D; j++)
      res.DDValue(i,j) = IfZero(x.DDValue(i,j)+x.DValue(i) * x.DValue(j),SCAL{0.},0.5/res.Value() * x.DDValue(i,j) - 0.25 / (x.Value()*res.Value()) * x.DValue(i) * x.DValue(j));

  return res;
}

// df(u)/dx  = exp(x) * du/dx
// d^2 f(u) / dx^2 = exp(x) * (du/dx)^2 + exp(x) * d^2u /dx^2
template <int D, typename SCAL>
INLINE AutoDiffDiff<D, SCAL> exp (AutoDiffDiff<D, SCAL> x)
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = exp(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * res.Value();
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = (x.DValue(k) * x.DValue(l)+x.DDValue(k,l)) * res.Value();
  return res;
}

using std::pow;
template <int D, typename SCAL>
INLINE AutoDiffDiff<D,SCAL> pow (AutoDiffDiff<D,SCAL> x, AutoDiffDiff<D,SCAL> y )
{
  return exp(log(x)*y);
}

template <int D, typename SCAL>
INLINE AutoDiffDiff<D, SCAL> log (AutoDiffDiff<D, SCAL> x)
{
  AutoDiffDiff<D, SCAL> res;
  res.Value() = log(x.Value());
  SCAL xinv = 1.0/x.Value();
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * xinv;
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = -xinv*xinv*x.DValue(k) * x.DValue(l) + xinv * x.DDValue(k,l);
  return res;
}



template <int D, typename SCAL>
INLINE AutoDiffDiff<D, SCAL> sin (AutoDiffDiff<D, SCAL> x)
{
  AutoDiffDiff<D, SCAL> res;
  SCAL s = sin(x.Value());
  SCAL c = cos(x.Value());
  
  res.Value() = s;
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * c;
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = -s * x.DValue(k) * x.DValue(l) + c * x.DDValue(k,l);
  return res;
}


template <int D, typename SCAL>
INLINE AutoDiffDiff<D, SCAL> cos (AutoDiffDiff<D, SCAL> x)
{
  AutoDiffDiff<D, SCAL> res;
  SCAL s = sin(x.Value());
  SCAL c = cos(x.Value());
  
  res.Value() = c;
  for (int k = 0; k < D; k++)
    res.DValue(k) = -s * x.DValue(k);
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = -c * x.DValue(k) * x.DValue(l) - s * x.DDValue(k,l);
  return res;
}

template <int D, typename SCAL>
INLINE AutoDiffDiff<D, SCAL> tan (AutoDiffDiff<D, SCAL> x)
{ return sin(x) / cos(x); }


template <int D, typename SCAL>
INLINE AutoDiffDiff<D, SCAL> atan (AutoDiffDiff<D, SCAL> x)
{
  AutoDiffDiff<D, SCAL> res;
  SCAL a = atan(x.Value());
  res.Value() = a;
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k)/(1+x.Value()*x.Value()) ;
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = -2*x.Value()/((1+x.Value()*x.Value())*(1+x.Value()*x.Value())) * x.DValue(k) * x.DValue(l) + x.DDValue(k,l)/(1+x.Value()*x.Value());
  return res;
}

template <int D, typename SCAL>
INLINE AutoDiffDiff<D, SCAL> atan2 (AutoDiffDiff<D, SCAL> x,AutoDiffDiff<D, SCAL> y)
{
  AutoDiffDiff<D, SCAL> res;
  SCAL a = atan2(x.Value(), y.Value());
  res.Value() = a;
  for (int k = 0; k < D; k++)
    res.DValue(k) = (x.Value()*y.DValue(k)-y.Value()*x.DValue(k))/(y.Value()*y.Value()+x.Value()*x.Value());

  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = (x.DValue(k)*y.DValue(l)+x.Value()*y.DDValue(l,k) - y.DValue(k)*x.DValue(l) - y.Value()*x.DDValue(l,k))/(y.Value()*y.Value()+x.Value()*x.Value()) - 2 * (x.Value()*y.DValue(k)-y.Value()*x.DValue(k)) * (x.Value()*x.DValue(k) + y.Value()*y.DValue(k))/( (y.Value()*y.Value()+x.Value()*x.Value()) * (y.Value()*y.Value()+x.Value()*x.Value()) );
  return res;
}



using std::acos;
template <int D, typename SCAL>
INLINE AutoDiffDiff<D,SCAL> acos (AutoDiffDiff<D,SCAL> x)
{
  AutoDiffDiff<D,SCAL> res;
  SCAL a = acos(x.Value());
  res.Value() = a;
  auto omaa = 1-x.Value()*x.Value();
  auto s = sqrt(omaa);
  SCAL da = -1 / s;
  SCAL dda = -x.Value() / (s*omaa);
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k)*da;
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = dda * x.DValue(k) * x.DValue(l) + da * x.DDValue(k,l);
  
  return res;
}


using std::acos;
template <int D, typename SCAL>
INLINE AutoDiffDiff<D,SCAL> asin (AutoDiffDiff<D,SCAL> x)
{
  AutoDiffDiff<D,SCAL> res;
  SCAL a = asin(x.Value());
  res.Value() = a;
  auto omaa = 1-x.Value()*x.Value();
  auto s = sqrt(omaa);
  SCAL da = 1 / s;
  SCAL dda = x.Value() / (s*omaa);
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k)*da;
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = dda * x.DValue(k) * x.DValue(l) + da * x.DDValue(k,l);
  
  return res;
}


template <int D, typename SCAL>
INLINE AutoDiffDiff<D, SCAL> sinh (AutoDiffDiff<D, SCAL> x)
{
  AutoDiffDiff<D, SCAL> res;
  SCAL sh = sinh(x.Value());
  SCAL ch = cosh(x.Value());
  
  res.Value() = sh;
  for (int k = 0; k < D; k++)
    res.DValue(k) = x.DValue(k) * ch;
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = sh * x.DValue(k) * x.DValue(l) + ch * x.DDValue(k,l);
  return res;
}


template <int D, typename SCAL>
INLINE AutoDiffDiff<D, SCAL> cosh (AutoDiffDiff<D, SCAL> x)
{
  AutoDiffDiff<D, SCAL> res;
  SCAL sh = sinh(x.Value());
  SCAL ch = cosh(x.Value());
  
  res.Value() = ch;
  for (int k = 0; k < D; k++)
    res.DValue(k) = sh * x.DValue(k);
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = ch * x.DValue(k) * x.DValue(l) + sh * x.DDValue(k,l);
  return res;
}

template <int D, typename SCAL>
INLINE AutoDiffDiff<D, SCAL> erf (AutoDiffDiff<D, SCAL> x)
{
  AutoDiffDiff<D, SCAL> res;
  SCAL derf = 2. / sqrt(M_PI) * exp(- x.Value() * x.Value());
  
  res.Value() = erf(x.Value());
  for (int k = 0; k < D; k++)
    res.DValue(k) = - derf * x.DValue(k);
  for (int k = 0; k < D; k++)
    for (int l = 0; l < D; l++)
      res.DDValue(k,l) = derf * (x.DDValue(k, l) - 2 * x.Value() * x.DValue(k) * x.DValue(l));
  return res;
}

using std::floor;
template<int D, typename SCAL>
INLINE AutoDiffDiff<D,SCAL> floor (const AutoDiffDiff<D,SCAL> & x)
{
  return floor(x.Value());
}

using std::ceil;
template<int D, typename SCAL>
INLINE AutoDiffDiff<D,SCAL> ceil (const AutoDiffDiff<D,SCAL> & x)
{
  return ceil(x.Value());  
}


template <int D, typename SCAL, typename TB, typename TC>
auto IfPos (AutoDiffDiff<D,SCAL> a, TB b, TC c) -> decltype(IfPos (a.Value(), b, c))
{
  return IfPos (a.Value(), b, c);
}

template <int D, typename SCAL>
INLINE AutoDiffDiff<D,SCAL> IfPos (SCAL /* SIMD<double> */ a, AutoDiffDiff<D,SCAL> b, AutoDiffDiff<D,SCAL> c)
{
  AutoDiffDiff<D,SCAL> res;
  res.Value() = IfPos (a, b.Value(), c.Value());
  for (int j = 0; j < D; j++)
  {
    res.DValue(j) = IfPos (a, b.DValue(j), c.DValue(j));
    res.DDValue(j) = IfPos (a, b.DDValue(j), c.DDValue(j));
  }
  return res;
}

template <int D, typename SCAL, typename TC>
INLINE AutoDiffDiff<D,SCAL> IfPos (SCAL /* SIMD<double> */ a, AutoDiffDiff<D,SCAL> b, TC c)
{
  return IfPos (a, b, AutoDiffDiff<D,SCAL> (c));
}




//@}

}


#endif
