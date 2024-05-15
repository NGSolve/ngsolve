#ifndef FILE_NGS_INTERPOLATE
#define FILE_NGS_INTERPOLATE

/**********************************************************************/
/* File:   interpolate.cpp                                            */
/* Author: L Kogler, M Neunteufel, J Schoeberl                        */
/* Date:   June 2020                                                  */
/**********************************************************************/

/* 
   Interpolation of CoefficientFunctions using
   dual shapes
*/

#include "fespace.hpp"

namespace ngcomp
{
  
  class InterpolateProxy : public ProxyFunction
  {
  protected:
    shared_ptr<CoefficientFunction> func;
    shared_ptr<FESpace> space;
    bool testfunction;
    shared_ptr<DifferentialOperator> final_diffop;
    int bonus_intorder;
  public:
    InterpolateProxy (shared_ptr<CoefficientFunction> func,
                      shared_ptr<FESpace> aspace,
                      bool testfunction,
                      shared_ptr<DifferentialOperator> diffop,                      
                      int bonus_intorder=0, VorB vb=VOL);

    shared_ptr<ProxyFunction> GetAdditionalProxy (string name) const override;

    shared_ptr<CoefficientFunction> Diff (const CoefficientFunction * var, shared_ptr<CoefficientFunction> dir) const override;    
  };
    
    
  NGS_DLL_HEADER
    shared_ptr<CoefficientFunction> InterpolateCF (shared_ptr<CoefficientFunction> func, shared_ptr<FESpace> space,
                                                   int bonus_intorder=0,
                                                   optional<string> opname=nullopt);
  
  
}

#endif // FILE_NGS_INTERPOLATE
