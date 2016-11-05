#ifndef COMPLEX_H_INCLUDED
#define COMPLEX_H_INCLUDED


namespace ngstd
{
template <class T>
class MyComplex {
public:
  enum { D = 2 };

  INLINE MyComplex() : a(0), b(0) { }
  INLINE MyComplex(const T &sa) : a(sa), b(0) { }
  INLINE MyComplex(const MyComplex<T> &S) : a(S.a), b(S.b) { }
  // MyComplex(const T U[D]) : a(U[0]), b(U[1]) { }
  INLINE MyComplex(const T &sa, const T &sb) : a(sa), b(sb) { }


  INLINE const MyComplex<T> & set(const T &, const T &);

  INLINE const MyComplex<T> operator + (const MyComplex<T> &) const;
  INLINE const MyComplex<T> operator - (const MyComplex<T> &) const;
  INLINE const MyComplex<T> operator * (const MyComplex<T> &) const;
  INLINE const MyComplex<T> operator / (const MyComplex<T> &) const;
  INLINE const MyComplex<T> operator + (const T &) const;
  INLINE const MyComplex<T> operator - (const T &) const;
  INLINE const MyComplex<T> operator * (const T &) const;
  INLINE const MyComplex<T> operator / (const T &) const;

  INLINE const MyComplex<T> & operator = (const MyComplex<T> &);
  INLINE const MyComplex<T> & operator += (const MyComplex<T> &);
  INLINE const MyComplex<T> & operator -= (const MyComplex<T> &);
  INLINE const MyComplex<T> & operator *= (const MyComplex<T> &);
  INLINE const MyComplex<T> & operator /= (const MyComplex<T> &);
  INLINE const MyComplex<T> & operator = (const T &);
  INLINE const MyComplex<T> & operator += (const T &);
  INLINE const MyComplex<T> & operator -= (const T &);
  INLINE const MyComplex<T> & operator *= (const T &);
  INLINE const MyComplex<T> & operator /= (const T &);

  INLINE T real() const { return a; }
  INLINE T imag() const { return b; }

  INLINE const MyComplex<T> operator - () const;
  INLINE const T & operator [] (const int) const;
  INLINE T & operator [] (const int);

  std::ofstream & saveBin(std::ofstream &) const;
  std::ifstream & loadBin(std::ifstream &);

  INLINE T normalize(void);

  INLINE ~MyComplex(void) { }

  union
  {
    T V[D];
    struct
    {
      T a, b;
    };
  };
};


template <class T>
INLINE const T absSq(const MyComplex<T> &A)
{
  return A.a*A.a + A.b*A.b;
}

template <class T>
INLINE const T arg(const MyComplex<T> &A)
{
  return atan2(A.b, A.a);
}

template <class T>
const MyComplex<T> conj(const MyComplex<T> &A)
{
  return MyComplex<T>(A.a, -A.b);
}

template <class T>
INLINE const T re(const MyComplex<T> &A)
{
  return A.a;
}

template <class T>
INLINE const T im(const MyComplex<T> &A)
{
  return A.b;
}

template <class T>
INLINE const T real(const MyComplex<T> &A)
{
  return A.a;
}

template <class T>
INLINE const T imag(const MyComplex<T> &A)
{
  return A.b;
}


template <class T>
INLINE MyComplex<T> sqrt(const MyComplex<T> &C)
{
  T l = abs(C);
  T p = 0.5*arg(C);

  return MyComplex<T>(l*cos(p), l*sin(p));
}

template <class T>
INLINE const T abs(const MyComplex<T> &A)
{
  return std::sqrt(A.a*A.a + A.b*A.b);
}

  using std::exp;

  /*
INLINE double exp (double a)
{
  return std::exp(a);
}
  */

template <class T>
INLINE MyComplex<T> exp(const MyComplex<T> &C)
{
  T k = exp(C.a);
  return MyComplex<T>(cos(C.b)*k, sin(C.b)*k);
}

template <class T>
INLINE MyComplex<T> pow(const MyComplex<T> &C, const T &m)
{
  T l = pow(abs(C), m);
  T p = m*arg(C);

  return MyComplex<T>(l*cos(p), l*sin(p));
}

template <class T>
INLINE MyComplex<T> pow(const T &a, const MyComplex<T> &C)
{
  T p = pow(a, C.a);
  T l = log(C.a);

  return MyComplex<T>(p*cos(C.b*l), p*sin(C.b*l));
}

template <class T>
std::ostream & operator << (std::ostream &vout, const MyComplex<T> &Q)
{
  return vout << "" << std::setw(14) << Q.a << " + " << std::setw(14) << Q.b << "i";
}

template <class T>
INLINE const MyComplex<T> operator + (const T &l, const MyComplex<T> &R)
{
  return MyComplex<T>(l + R.a, R.b);
}

template <class T>
INLINE const MyComplex<T> operator - (const T &l, const MyComplex<T> &R)
{
  return MyComplex<T>(l - R.a, -R.b);
}

template <class T>
INLINE const MyComplex<T> operator * (const T &l, const MyComplex<T> &R)
{
  return MyComplex<T>(l*R.a, l*R.b);
}

template <class T>
INLINE const MyComplex<T> operator / (const T &l, const MyComplex<T> &R)
{
  T z = absSq(R);
  return MyComplex<T>(l*R.a/z, -l*R.b/z);
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::set(const T &sa, const T &sb)
{
  a = sa; b = sb;
  return *this;
}

template <class T>
INLINE const MyComplex<T> MyComplex<T>::operator + (const MyComplex<T> &R) const
{
  return MyComplex<T>(a + R.a, b + R.b);
}

template <class T>
INLINE const MyComplex<T> MyComplex<T>::operator - (const MyComplex<T> &R) const
{
  return MyComplex<T>(a - R.a, b - R.b);
}

template <class T>
INLINE const MyComplex<T> MyComplex<T>::operator * (const MyComplex<T> &R) const
{
  return MyComplex<T>(a*R.a - b*R.b, a*R.b + b*R.a);
}

template <class T>
INLINE const MyComplex<T> MyComplex<T>::operator / (const MyComplex<T> &R) const
{
  T z = abs(R);
  return MyComplex<T>((a*R.a + b*R.b)/z, (b*R.a - a*R.b)/z);
}

template <class T>
INLINE const MyComplex<T> MyComplex<T>::operator + (const T &r) const
{
  return MyComplex<T>(a + r, b);
}

template <class T>
INLINE const MyComplex<T> MyComplex<T>::operator - (const T &r) const
{
  return MyComplex<T>(a - r, b);
}

template <class T>
INLINE const MyComplex<T> MyComplex<T>::operator * (const T &r) const
{
  return MyComplex<T>(a*r, b*r);
}

template <class T>
INLINE const MyComplex<T> MyComplex<T>::operator / (const T &r) const
{
  return MyComplex<T>(a/r, b/r);
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::operator = (const MyComplex<T> &R)
{
  return set(R.a, R.b);
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::operator += (const MyComplex<T> &R)
{
  return set(a + R.a, b + R.b);
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::operator -= (const MyComplex<T> &R)
{
  return set(a - R.a, b - R.b);
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::operator *= (const MyComplex<T> &R)
{
  return set(a*R.a - a*R.b, a*R.b + b*R.a);
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::operator /= (const MyComplex<T> &R)
{
  T z = abs(R);
  return set((a*R.a + a*R.b)/z, (b*R.a - a*R.b)/z);
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::operator = (const T &r)
{
  return set(r, static_cast<T>(0));
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::operator += (const T &r)
{
  return set(a + r, b);
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::operator -= (const T &r)
{
  return set(a - r, b);
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::operator *= (const T &r)
{
  return set(a*r, b*r);
}

template <class T>
INLINE const MyComplex<T> & MyComplex<T>::operator /= (const T &r)
{
  return set(a/r, b/r);
}

template <class T>
INLINE const MyComplex <T> MyComplex<T>::operator - () const
{
  return MyComplex(-a, -b);
}

template <class T>
INLINE const T & MyComplex<T>::operator [] (const int en) const
{
  return V[en];
}

template <class T>
INLINE T & MyComplex<T>::operator [] (const int en)
{
  return V[en];
}

template <class T>
inline std::ofstream & MyComplex<T>::saveBin(std::ofstream &savf) const
{
  savf.write((char *)V, D*sizeof(T));
  return savf;
}

template <class T>
inline std::ifstream & MyComplex<T>::loadBin(std::ifstream &loaf)
{
  loaf.read((char *)V, D*sizeof(T));
  return loaf;
}

INLINE bool operator== (MyComplex<double> a, MyComplex<double> b)
{
  return (real(a) == real(b)) && (imag(a) == imag(b));
}



template <class T>
INLINE T MyComplex<T>::normalize(void)
{
  T l = sqrt(a*a + b*b);

  if(fabs(l) > (T)0)
    {
      a /= l;
      b /= l;
    }

  return l;
}

}

#endif 
