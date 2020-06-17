/**********************************************************************/
/* File:   interpolate.cpp                                            */
/* Author: L Kogler, M Neunteufel, J Schoeberl                        */
/* Date:   June 2020                                                  */
/**********************************************************************/

/* 
   Interpolation of CoefficientFunctions using
   dual shapes
*/


namespace ngcomp
{
  
  class InterpolateProxy : public ProxyFunction
  {
  protected:
    shared_ptr<CoefficientFunction> func;
    shared_ptr<FESpace> space;
    int bonus_intorder;
  public:
    InterpolateProxy (shared_ptr<CoefficientFunction> func,
                      shared_ptr<FESpace> aspace,
                      bool testfunction,
                      int bonus_intorder = 0);
  };
    
    
  NGS_DLL_HEADER
    shared_ptr<CoefficientFunction> InterpolateCF (shared_ptr<CoefficientFunction> func, shared_ptr<FESpace> space,
                                                   int bonus_intorder = 0);
  
  
}
