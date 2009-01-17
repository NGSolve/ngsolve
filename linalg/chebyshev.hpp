#ifndef FILE_CHEBYSHEV
#define FILE_CHEBYSHEV

/**************************************************************************/
/* File:   chebyshev.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   30. Jun. 01                                                    */
/**************************************************************************/

/**
   Chebyshev iteraion
  */ 
class ChebyshevIteration : public BaseMatrix
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
  ///
  void SetBounds (double almin, double almax);
  ///
  virtual void Mult (const BaseVector & v, BaseVector & prod) const;
  ///
  virtual BaseVector * CreateVector () const;
};

#endif
