#ifndef FILE_EIGEN
#define FILE_EIGEN

/**************************************************************************/
/* File:   eigen.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   5. Jul. 96                                                     */
/**************************************************************************/


#include "basematrix.hpp"

namespace ngla
{

  /**
     Lanczos Eigen system calculation
  */ 

  class NGS_DLL_HEADER EigenSystem
  {
    ///
    const BaseMatrix *a, *c;
    ///
    Array<double> ai, bi;
    ///
    double prec;
    ///
    int maxsteps;
  
  public:
    ///
    EigenSystem (const BaseMatrix & aa);
    ///
    EigenSystem (const BaseMatrix & aa, const BaseMatrix & ac);
    ///
    void SetMatrix (const BaseMatrix & aa);
    ///
    void SetPrecond (const BaseMatrix & ac);
    ///
    void SetMaxSteps (int amaxsteps);
    ///
    void SetPrecision (double aprec);

    ///
    int Calc();

    ///
    double EigenValue (int nr) const;
    ///
    double MaxEigenValue () const;
    ///
    int NumEigenValues () const;
    ///
    void PrintEigenValues (ostream & ost) const;
  };

}

#endif
