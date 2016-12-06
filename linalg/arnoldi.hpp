#ifndef FILE_ARNOLDI
#define FILE_ARNOLDI


/**************************************************************************/
/* File:   arnoldi.hpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   5. Jul. 96                                                     */
/**************************************************************************/

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
    const BaseMatrix & a;
    const BaseMatrix & b;
    shared_ptr<BitArray> freedofs;
    SCAL shift;

  public:
    Arnoldi (const BaseMatrix & aa, const BaseMatrix & ab, shared_ptr<BitArray> afreedofs = nullptr)
      : a(aa), b(ab), freedofs(afreedofs)
    { 
      shift = 1.0;
    }

    void SetShift (SCAL ashift)
    { shift = ashift; }

    void Calc (int numval, Array<Complex> & lam, int nev, 
               Array<shared_ptr<BaseVector>> & evecs, 
               const BaseMatrix * pre = NULL) const;
  };
}

#endif
