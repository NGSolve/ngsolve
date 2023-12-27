#ifndef FILE_THCURLFE
#define FILE_THCURLFE

/*********************************************************************/
/* File:   thcurlfe.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   5. Sep. 2013                                              */
/*********************************************************************/



#include "hcurlfe.hpp"

namespace ngfem
{


  
  /**
     HCurlHighOrderFE of shape ET.
     provides access functions, shape functions are provided by CalcShape template
  */
  template <ELEMENT_TYPE ET, typename SHAPES,
            typename BASE = HCurlFiniteElement<ET_trait<ET>::DIM>>
  class T_HCurlHighOrderFiniteElement : public BASE 
  {
  protected:
    enum { DIM = ET_trait<ET>::DIM };
  
    using BASE::DIM_CURL;
    using BASE::ndof;
    using BASE::order;

  public:

    NGS_DLL_HEADER T_HCurlHighOrderFiniteElement () { ; }
    HD virtual ELEMENT_TYPE ElementType() const override { return ET; }
    
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (TIP<DIM,Tx> tip, TFA & shape) const
    { 
      static_cast<const SHAPES*> (this) -> T_CalcShape (tip, shape);
    }

    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            BareSliceMatrix<> shape) const override;

    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                                BareSliceMatrix<> curlshape) const override;
#ifndef FASTCOMPILE
    
    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & mip,
                                  BareSliceMatrix<> shape) const override;

    virtual void CalcMappedShape (const BaseMappedIntegrationRule & bmir, BareSliceMatrix<> shapes) const override;

    virtual void CalcMappedShape (const SIMD<BaseMappedIntegrationPoint> & bmip,
                                  BareSliceMatrix<SIMD<double>> shape) const override;

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                  BareSliceMatrix<SIMD<double>> shapes) const override;

    virtual void CalcMappedCurlShape (const BaseMappedIntegrationPoint & mip,
                                      BareSliceMatrix<> curlshape) const override;

    virtual void CalcMappedCurlShape (const BaseMappedIntegrationRule & mir, 
                                      BareSliceMatrix<> curlshape) const override;

    virtual void CalcMappedCurlShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                      BareSliceMatrix<SIMD<double>> curlshapes) const override;
    
    virtual Vec <DIM_CURL_(DIM)>
    EvaluateCurlShape (const IntegrationPoint & ip, 
                       BareSliceVector<double> x,
                       LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void 
    EvaluateCurl (const IntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<> curl) const override;

    using BASE::Evaluate;
    NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const override;
    
    NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const override;        
    NGS_DLL_HEADER virtual void EvaluateCurl (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const override;     
    NGS_DLL_HEADER virtual void EvaluateCurl (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const;  // actually not in base-class !!!!

    NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                             BareSliceVector<> coefs) const override;
    NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values,
                                             BareSliceVector<Complex> coefs) const override;
    NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                                 BareSliceVector<> coefs) const override;
    NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values,
                                                 BareSliceVector<Complex> coefs) const override;
#endif
    
  };







  /**
     Base-element for template polymorphism.
     Barton and Nackman Trick
  */
  
  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  class T_HCurlFiniteElementFO 
    : public T_HCurlHighOrderFiniteElement<ET, FEL>
  {
  public:

    T_HCurlFiniteElementFO ()
    {
      this -> ndof = NDOF;
      this -> order = ORDER;
    }
  };

}



#endif
