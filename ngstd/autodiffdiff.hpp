#ifndef FILE_AUTODIFFDIFF
#define FILE_AUTODIFFDIFF

/**************************************************************************/
/* File:   autodiffdiff.hpp                                               */
/* Author: Joachim Schoeberl                                              */
/* Date:   13. June. 05                                                   */
/**************************************************************************/

namespace ngstd
{

// Automatic second differentiation datatype


/**
   Datatype for automatic differentiation.  Contains function value,
   D first derivatives, and D*D second derivatives. Algebraic operations are
   overloaded by using product-rule etc. etc.
**/
template <int D>
class AutoDiffDiff
{
  double val;
  double dval[D];
  double ddval[D*D];
public:

  typedef AutoDiffDiff<D> TELEM;
  typedef double TSCAL;


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
  AutoDiffDiff  (double aval) throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
    for (int i = 0; i < D*D; i++)
      ddval[i] = 0;
  }

  /// initial object with value and derivative
  AutoDiffDiff  (const AutoDiff<D> & ad2) throw()
  {
    val = ad2.Value();
    for (int i = 0; i < D; i++)
      dval[i] = ad2.DValue(i);
    for (int i = 0; i < D*D; i++)
      ddval[i] = 0;
  }

  /// init object with (val, e_diffindex)
  AutoDiffDiff  (double aval, int diffindex)  throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
    for (int i = 0; i < D*D; i++)
      ddval[i] = 0;
    dval[diffindex] = 1;
  }

  /// assign constant value
  AutoDiffDiff & operator= (double aval) throw()
  {
    val = aval;
    for (int i = 0; i < D; i++)
      dval[i] = 0;
    for (int i = 0; i < D*D; i++)
      ddval[i] = 0;
    return *this;
  }

  /// returns value
  double Value() const throw() { return val; }

  /// returns partial derivative
  double DValue (int i) const throw() { return dval[i]; }

  /// returns partial derivative
  double DDValue (int i) const throw() { return ddval[i]; }

  /// returns partial derivative
  double DDValue (int i, int j) const throw() { return ddval[i*D+j]; }

  /// access value
  double & Value() throw() { return val; }

  /// accesses partial derivative 
  double & DValue (int i) throw() { return dval[i]; }

  /// accesses partial derivative 
  double & DDValue (int i) throw() { return ddval[i]; }

  /// accesses partial derivative 
  double & DDValue (int i, int j) throw() { return ddval[i*D+j]; }

  /// add autodiffdiff object
  AutoDiffDiff<D> & operator+= (const AutoDiffDiff<D> & y) throw()
  {
    val += y.val;
    for (int i = 0; i < D; i++)
      dval[i] += y.dval[i];
    for (int i = 0; i < D*D; i++)
      ddval[i] += y.ddval[i];
    return *this;
  }

  /// subtract autodiffdiff object
  AutoDiffDiff<D> & operator-= (const AutoDiffDiff<D> & y) throw()
  {
    val -= y.val;
    for (int i = 0; i < D; i++)
      dval[i] -= y.dval[i];
    for (int i = 0; i < D*D; i++)
      ddval[i] -= y.ddval[i];
    return *this;
  }

  /// multiply with autodiffdiff object
  AutoDiffDiff<D> & operator*= (const AutoDiffDiff<D> & y) throw()
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
  AutoDiffDiff<D> & operator*= (const double & y) throw()
  {
    for ( int i = 0; i < D*D; i++ )
      ddval[i] *= y;
    for (int i = 0; i < D; i++)
      dval[i] *= y;
    val *= y;
    return *this;
  }

  /// divide by scalar
  AutoDiffDiff<D> & operator/= (const double & y) throw()
  {
    double iy = 1.0 / y;
    for ( int i = 0; i < D*D; i++ )
      ddval[i] *= iy;
    for (int i = 0; i < D; i++)
      dval[i] *= iy;
    val *= iy;
    return *this;
  }

  /// same value
  bool operator== (double val2) throw()
  {
    return val == val2;
  }

  /// different values 
  bool operator!= (double val2) throw()
  {
    return val != val2;
  }

  /// less 
  bool operator< (double val2) throw()
  {
    return val < val2;
  }
  
  /// greater
  bool operator> (double val2) throw()
  {
    return val > val2;
  }
};


//@{  AutoDiff helper functions.

/// Prints AudoDiffDiff
template<int D>
inline ostream & operator<< (ostream & ost, const AutoDiffDiff<D> & x)
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
template<int D>
inline AutoDiffDiff<D> operator+ (const AutoDiffDiff<D> & x, const AutoDiffDiff<D> & y) throw()
{
  AutoDiffDiff<D> res;
  res.Value () = x.Value()+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i) + y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = x.DDValue(i) + y.DDValue(i);
  return res;
}


/// 
template<int D>
inline AutoDiffDiff<D> operator- (const AutoDiffDiff<D> & x, const AutoDiffDiff<D> & y) throw()
{
  AutoDiffDiff<D> res;
  res.Value() = x.Value()-y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i) - y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = x.DDValue(i) - y.DDValue(i);
  return res;
}


///
template<int D>
inline AutoDiffDiff<D> operator+ (double x, const AutoDiffDiff<D> & y) throw()
{
  AutoDiffDiff<D> res;
  res.Value() = x+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = y.DDValue(i);
  return res;
}

///
template<int D>
inline AutoDiffDiff<D> operator+ (const AutoDiffDiff<D> & y, double x) throw()
{
  AutoDiffDiff<D> res;
  res.Value() = x+y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = y.DDValue(i);
  return res;
}


///
template<int D>
inline AutoDiffDiff<D> operator- (const AutoDiffDiff<D> & x) throw()
{
  AutoDiffDiff<D> res;
  res.Value() = -x.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = -x.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = -x.DDValue(i);
  return res;
}

///
template<int D>
inline AutoDiffDiff<D> operator- (const AutoDiffDiff<D> & x, double y) throw()
{
  AutoDiffDiff<D> res;
  res.Value() = x.Value()-y;
  for (int i = 0; i < D; i++)
    res.DValue(i) = x.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = x.DDValue(i);
  return res;
}

///
template<int D>
inline AutoDiffDiff<D> operator- (double x, const AutoDiffDiff<D> & y) throw()
{
  AutoDiffDiff<D> res;
  res.Value() = x-y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = -y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = -y.DDValue(i);
  return res;
}


///
template<int D>
inline AutoDiffDiff<D> operator* (double x, const AutoDiffDiff<D> & y) throw()
{
  AutoDiffDiff<D> res;
  res.Value() = x*y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x*y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = x*y.DDValue(i);
  return res;
}

///
template<int D>
inline AutoDiffDiff<D> operator* (const AutoDiffDiff<D> & y, double x) throw()
{
  AutoDiffDiff<D> res;
  res.Value() = x*y.Value();
  for (int i = 0; i < D; i++)
    res.DValue(i) = x*y.DValue(i);
  for (int i = 0; i < D*D; i++)
    res.DDValue(i) = x*y.DDValue(i);
  return res;
}

///
template<int D>
inline AutoDiffDiff<D> operator* (const AutoDiffDiff<D> & x, const AutoDiffDiff<D> & y) throw()
{
  AutoDiffDiff<D> res;
  double hx = x.Value();
  double hy = y.Value();

  res.Value() = hx*hy;
  for (int i = 0; i < D; i++)
    res.DValue(i) = hx*y.DValue(i) + hy*x.DValue(i);

  for (int i = 0; i < D; i++)
    for (int j = 0; j < D; j++)
      res.DDValue(i,j) = hx * y.DDValue(i,j) + hy * x.DDValue(i,j)
	+ x.DValue(i) * y.DValue(j) + x.DValue(j) * y.DValue(i);

  return res;
}



template<int D>
inline AutoDiffDiff<D> Inv (const AutoDiffDiff<D> & x)
{
  AutoDiffDiff<D> res(1.0 / x.Value());
  for (int i = 0; i < D; i++)
    res.DValue(i) = -x.DValue(i) / (x.Value() * x.Value());
  cout << "ADD::Inv not implemented" << endl;
  return res;
}


template<int D>
inline AutoDiffDiff<D> operator/ (const AutoDiffDiff<D> & x, const AutoDiffDiff<D> & y)
{
  return x * Inv (y);
}

template<int D>
inline AutoDiffDiff<D> operator/ (const AutoDiffDiff<D> & x, double y)
{
  return (1/y) * x;
}

template<int D>
inline AutoDiffDiff<D> operator/ (double x, const AutoDiffDiff<D> & y)
{
  return x * Inv(y);
}


template<int D>
inline AutoDiffDiff<D> sqrt (const AutoDiffDiff<D> & x)
{
  AutoDiffDiff<D> res;
  res.Value() = std::sqrt(x.Value());
  for (int j = 0; j < D; j++)
    res.DValue(j) = 0.5 / res.Value() * x.DValue(j);

  
  for (int i = 0; i < D; i++)
    for (int j = 0; j < D; j++)
      res.DDValue(i,j) = 0.5/res.Value() * x.DDValue(i,j) - 0.25 / (x.Value()*res.Value()) * x.DValue(i) * x.DValue(j);

  return res;
}





//@}

}


#endif
