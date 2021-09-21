#ifndef FILE_THCURLFE
#define FILE_THCURLFE

/*********************************************************************/
/* File:   thcurlfe.hpp                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   5. Sep. 2013                                              */
/*********************************************************************/

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

    /*
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (const Tx & hx, TFA & shape) const
    { 
      static_cast<const SHAPES*> (this) -> T_CalcShape (hx, shape);
    }
    */
    
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (TIP<DIM,Tx> tip, TFA & shape) const
    { 
      static_cast<const SHAPES*> (this) -> T_CalcShape (tip, shape);
    }

    // const SHAPES & Cast() const { return static_cast<const SHAPES&> (*this); }
    
    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceMatrix<> shape) const override;

    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                                SliceMatrix<> curlshape) const override;
#ifndef FASTCOMPILE
    
    virtual void CalcMappedShape (const BaseMappedIntegrationPoint & mip,
                                  SliceMatrix<> shape) const override;

    virtual void CalcMappedShape (const BaseMappedIntegrationRule & bmir, SliceMatrix<> shapes) const override;

    virtual void CalcMappedShape (const SIMD<BaseMappedIntegrationPoint> & bmip,
                                  BareSliceMatrix<SIMD<double>> shape) const override;

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                  BareSliceMatrix<SIMD<double>> shapes) const override;

    virtual void CalcMappedCurlShape (const BaseMappedIntegrationPoint & mip,
                                      SliceMatrix<> curlshape) const override;

    virtual void CalcMappedCurlShape (const MappedIntegrationRule<DIM,DIM> & mir, 
                                      SliceMatrix<> curlshape) const override;

    virtual void CalcMappedCurlShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                      BareSliceMatrix<SIMD<double>> curlshapes) const override;
    
    virtual Vec <DIM_CURL_(DIM)>
    EvaluateCurlShape (const IntegrationPoint & ip, 
                       BareSliceVector<double> x,
                       LocalHeap & lh) const override;

    NGS_DLL_HEADER virtual void 
    EvaluateCurl (const IntegrationRule & ir, BareSliceVector<> coefs, FlatMatrixFixWidth<DIM_CURL_(DIM)> curl) const override;
#endif

    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const override;
    
    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const override;        
    HD NGS_DLL_HEADER virtual void EvaluateCurl (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const override;     
    HD NGS_DLL_HEADER virtual void EvaluateCurl (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const;  // actually not in base-class !!!!

    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                             BareSliceVector<> coefs) const override;
    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values,
                                             BareSliceVector<Complex> coefs) const override;
    HD NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                                 BareSliceVector<> coefs) const override;
    HD NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values,
                                                 BareSliceVector<Complex> coefs) const override; 
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
