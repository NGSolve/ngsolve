#ifndef FILE_BASESCALAR
#define FILE_BASESCALAR

/*********************************************************************/
/* File:   basescalar.hpp                                            */
/* Author: Natalia Tylek, Joachim Schoeberl                          */
/* Date:   Mar. 2026                                                 */
/*********************************************************************/

#include <bla.hpp>

namespace ngla
{
  using namespace ngbla;
  
  class Scalar : public std::variant<double, Complex>
  {
  public:
    using std::variant<double, Complex>::variant;
    Scalar operator- () {
      return std::visit([](auto val) {
        return Scalar(-val);
      }, *this);
    }
  };


  class BaseScalar
  {
    Scalar scal;

  public:
    BaseScalar() : scal{double(0)} { }
    virtual ~BaseScalar() = default;
    virtual bool IsComplex() const { return std::holds_alternative<Complex>(scal); } 
    virtual void Set (double d) { scal = d; }
    virtual void Set (Complex c) { scal = c; }
    virtual double GetD () const { return get<double>(scal); }
    virtual Complex GetC () const { return get<Complex>(scal); }
  };


  INLINE ostream & operator<< (ostream & ost, const BaseScalar & scal)
  {
    if (scal.IsComplex())
      ost << scal.GetC();
    else
      ost << scal.GetD();
    return ost;
  }
  
}


#endif
