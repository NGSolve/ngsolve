#include "expr.hpp"

// arithmetics with error control (better say error monitoring)


namespace ngbla
{
  
  
  class EC_double
  {
    double val;
    double err;
  public:
    EC_double () : val(0), err(0) { }
    EC_double (double a)
      : val(a), err(fabs(a)*1e-16) { }
    EC_double (double a, double aerr)
      : val(a), err(aerr)
    {
      Check();
    }
    EC_double (EC_double&) = default;    

    EC_double & operator= (EC_double&) = default;
    EC_double & operator= (double a)
    {
      val = a; err = fabs(a)*1e-16;
      return *this;
    }
    
    double Val() const { return val; }
    double Err() const { return err; }    
    operator double() const { return val; }

    EC_double & operator*= (EC_double b)
    {
      *this = EC_double(*this)*b;
      return *this;
    }

    EC_double & operator+= (EC_double b)
    {
      val += b.Val();
      err += b.Err();
      Check();
      return *this;
    }
    
    EC_double & operator-= (EC_double b)
    {
      *this = EC_double(*this) + (-b);
      return *this;
    }

    EC_double operator- () const
    {
      return EC_double(-val, err);
    }

    void Check() const
    {
      if (err > 1e-4 * abs(val))
        {
          // throw Exception("roundoff is big, val = "+ToString(val)+", err = ", ToString(err));
          cout << "roundoff is big, val = "+ToString(val)+", err = " + ToString(err) << endl;
        }
    }
  };


  inline EC_double operator* (EC_double a, EC_double b)
  {
    return EC_double(a.Val()*b.Val(), abs(a.Val())*b.Err()+abs(b.Val())*a.Err());
  }

  inline EC_double operator* (EC_double a, double b)
  {
    return a * EC_double(b);
  }

  inline EC_double operator* (double b, EC_double a)
  {
    return a * EC_double(b);
  }
  
  inline EC_double operator+ (EC_double a, EC_double b)
  {
    return EC_double(a.Val()+b.Val(), a.Err()+b.Err());
  }

  typedef std::complex<EC_double> EC_Complex;

  inline auto operator* (double a, EC_Complex b)
  {
    return EC_double(a) * b;
  }

  inline auto operator* (EC_Complex a, const Complex b)
  {
    return a*EC_Complex(b);
  }

  inline double L2Norm2 (EC_Complex val)
  {
    return sqr(val.real().Val()) + sqr(val.imag().Val());
  }


  // typedef std::complex<__128float> Complex128;
  // typedef std::complex<double> Complex128;  



  
  template<> struct is_scalar_type<ngfem::EC_Complex> { static constexpr bool value = true; };
  template<> struct is_scalar_type<ngfem::EC_double> { static constexpr bool value = true; };
  template <> inline constexpr bool IsComplex<ngfem::EC_Complex> () { return true; }

  
}
