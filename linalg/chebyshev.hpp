#ifndef FILE_CHEBYSHEV
#define FILE_CHEBYSHEV

/**************************************************************************/
/* File:   chebyshev.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   30. Jun. 01                                                    */
/**************************************************************************/


#include "basematrix.hpp"

namespace ngla
{

  /**
     Chebyshev iteraion
  */ 
  class NGS_DLL_HEADER ChebyshevIteration : public BaseMatrix
  {
  protected:
    ///
    const BaseMatrix *a, *c;
    ///
    int steps;
    ///
    double lmin, lmax;
  public:
    ///
    ChebyshevIteration (const BaseMatrix & aa, const BaseMatrix & ac, int steps);
    
    bool IsComplex() const override { return a->IsComplex(); } 
    ///
    void SetBounds (double almin, double almax);
    ///
    void Mult (const BaseVector & v, BaseVector & prod) const override;
    ///
    AutoVector CreateRowVector () const override { return a->CreateColVector(); }
    AutoVector CreateColVector () const override { return a->CreateRowVector(); }
  };

}

#endif
