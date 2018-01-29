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
    HD virtual ELEMENT_TYPE ElementType() const { return ET; }

    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx hx[2], TFA & shape) const
    { 
      static_cast<const SHAPES*> (this) -> T_CalcShape (hx, shape);
    }

    virtual void CalcShape (const IntegrationPoint & ip, 
                            SliceMatrix<> shape) const;

    virtual void CalcCurlShape (const IntegrationPoint & ip, 
                                SliceMatrix<> curlshape) const;
#ifndef FASTCOMPILE
    virtual void CalcMappedShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                                  SliceMatrix<> shape) const;

    virtual void CalcMappedShape (const MappedIntegrationRule<DIM,DIM> & mir, 
                                  SliceMatrix<> shape) const;

    virtual void CalcMappedShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                  BareSliceMatrix<SIMD<double>> shapes) const;

    virtual void CalcMappedCurlShape (const MappedIntegrationPoint<DIM,DIM> & mip,
                                      SliceMatrix<> curlshape) const;

    virtual void CalcMappedCurlShape (const MappedIntegrationRule<DIM,DIM> & mir, 
                                      SliceMatrix<> curlshape) const;

    virtual void CalcMappedCurlShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                      BareSliceMatrix<SIMD<double>> curlshapes) const;
    
    virtual Vec <DIM_CURL_(DIM)>
    EvaluateCurlShape (const IntegrationPoint & ip, 
                       FlatVector<double> x,
                       LocalHeap & lh) const;

    NGS_DLL_HEADER virtual void 
    EvaluateCurl (const IntegrationRule & ir, FlatVector<> coefs, FlatMatrixFixWidth<DIM_CURL_(DIM)> curl) const;
#endif

    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const;
    
    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const;        
    HD NGS_DLL_HEADER virtual void EvaluateCurl (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<> coefs, BareSliceMatrix<SIMD<double>> values) const;     
    HD NGS_DLL_HEADER virtual void EvaluateCurl (const SIMD_BaseMappedIntegrationRule & ir, BareSliceVector<Complex> coefs, BareSliceMatrix<SIMD<Complex>> values) const;

    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                             BareSliceVector<> coefs) const;
    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values,
                                             BareSliceVector<Complex> coefs) const;
    HD NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<double>> values,
                                                 BareSliceVector<> coefs) const;
    HD NGS_DLL_HEADER virtual void AddCurlTrans (const SIMD_BaseMappedIntegrationRule & ir, BareSliceMatrix<SIMD<Complex>> values,
                                                 BareSliceVector<Complex> coefs) const;


    
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
