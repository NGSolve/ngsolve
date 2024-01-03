#ifndef FILE_THDIVFE
#define FILE_THDIVFE

/*********************************************************************/
/* File:   thdivfe.hpp                                               */
/* Author: Joachim Schoeberl                                         */
/* Date:   5. Jul. 2001                                              */
/*********************************************************************/

#include "hdivfe.hpp"
#include "thcurlfe.hpp"  // for Du

namespace ngfem
{






  template <class FEL, ELEMENT_TYPE ET>
  class T_HDivFiniteElement 
    : public HDivFiniteElement<ET_trait<ET>::DIM>
    
  {
    enum { DIM = ET_trait<ET>::DIM };

  public:

    virtual void CalcShape (const IntegrationPoint & ip, 
			    BareSliceMatrix<> shape) const override;
    
    virtual void CalcDivShape (const IntegrationPoint & ip, 
			       BareSliceVector<> divshape) const override;

#ifndef FASTCOMPILE

    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & bmip,
                                  BareSliceMatrix<> shape) const override;

    virtual void CalcMappedShape (const BaseMappedIntegrationRule & bmir, BareSliceMatrix<> shapes) const override;

    virtual void CalcMappedShape (const SIMD<MappedIntegrationPoint<DIM,DIM>> & mip,
				  BareSliceMatrix<SIMD<double>> shape) const override;

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                  BareSliceMatrix<SIMD<double>> shapes) const override;

    virtual void CalcMappedNormalShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                        BareSliceMatrix<SIMD<double>> shapes) const override;
    
    using HDivFiniteElement<ET_trait<ET>::DIM>::CalcMappedDivShape;
    virtual void CalcMappedDivShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                     BareSliceMatrix<SIMD<double>> divshapes) const override;

    virtual void Evaluate (const IntegrationRule & ir, 
			   FlatVector<double> coefs, 
			   BareSliceMatrix<> vals) const override;

    virtual void EvaluateTrans (const IntegrationRule & ir, 
                                BareSliceMatrix<> vals,
                                FlatVector<double> coefs) const override;


    virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const override;
    virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values, BareSliceVector<> coefs) const override;

    virtual void EvaluateDiv (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareVector<SIMD<double>> values) const override;
    
    virtual void AddDivTrans (const SIMD_BaseMappedIntegrationRule & ir, BareVector<SIMD<double>> values,
                              BareSliceVector<> coefs) const override;
    
    
#endif
  };

}


#endif
