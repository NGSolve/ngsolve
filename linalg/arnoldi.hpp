#ifndef FILE_ARNOLDI
#define FILE_ARNOLDI


/**************************************************************************/
/* File:   arnoldi.hpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   5. Jul. 96                                                     */
/**************************************************************************/


#include "basematrix.hpp"

namespace ngla
{
  /**
     Arnoldi Eigenvalue Solver.

     Solve the generalized evp
     
     A x = lam B x

     B must by symmetric and (in theory) positive definite
     A can be non-symmetric

     It uses a shift-and-invert Arnoldi method 
   */

  template <typename SCAL>
  class NGS_DLL_HEADER Arnoldi
  {
    shared_ptr<BaseMatrix> a;
    shared_ptr<BaseMatrix> b;
    shared_ptr<BitArray> freedofs;
    SCAL shift;
    optional<string> inversetype;
  public:
    Arnoldi (shared_ptr<BaseMatrix> aa, shared_ptr<BaseMatrix> ab, shared_ptr<BitArray> afreedofs = nullptr)
      : a(aa), b(ab), freedofs(afreedofs)
    { 
      shift = 1.0;
    }

    void SetShift (SCAL ashift)
    { shift = ashift; }
    void SetInverseType (optional<string> ainv)
    { inversetype = ainv; }

    void Calc (int numval, Array<Complex> & lam, int nev, 
               Array<shared_ptr<BaseVector>> & evecs, 
               shared_ptr<BaseMatrix> pre = nullptr) const;
  };
}

#endif
