#ifndef FILE_TSCALARFE
#define FILE_TSCALARFE

/*********************************************************************/
/* File:   tscalarfe.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#include "scalarfe.hpp"


namespace ngfem
{

  
  /**
     Base-element for template polymorphism.
     Barton and Nackman Trick for elements with non-static CalcShape method
  */

  template <class FEL, ELEMENT_TYPE ET, 
            class BASE = ScalarFiniteElement<ET_trait<ET>::DIM> >

  class T_ScalarFiniteElement : public BASE   
  {
  public:
    // enum { DIM = ET_trait<ET>::DIM };
    static constexpr int DIM = ET_trait<ET>::DIM;

    // using BASE::eltype;
    using BASE::ndof;
    using BASE::order;

    INLINE T_ScalarFiniteElement () { ; }
    // virtual ~T_ScalarFiniteElement() { ; }

    HD virtual ELEMENT_TYPE ElementType() const final { return ET; }
    // HD NGS_DLL_HEADER virtual int Dim () const override { return DIM; } 

    
    HD NGS_DLL_HEADER virtual void CalcShape (const IntegrationPoint & ip, 
					      BareSliceVector<> shape) const override;

    HD NGS_DLL_HEADER virtual void CalcDShape (const IntegrationPoint & ip, 
					       BareSliceMatrix<> dshape) const override;

    
#ifndef FASTCOMPILE
    HD NGS_DLL_HEADER virtual void CalcShape (const IntegrationRule & ir, 
                                              BareSliceMatrix<> shape) const override;
    /// compute shape, row is shape nr, col is ip nr
    HD NGS_DLL_HEADER 
    virtual void CalcShape (const SIMD_IntegrationRule & ir, 
                            BareSliceMatrix<SIMD<double>> shape) const override;
    
    HD NGS_DLL_HEADER virtual double Evaluate (const IntegrationPoint & ip, 
					       BareSliceVector<double> x) const override;
    
    HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, 
					     BareSliceVector<double> coefs, 
					     BareSliceVector<double> vals) const override;

    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_IntegrationRule & ir,
                                             BareSliceVector<> coefs,
                                             BareVector<SIMD<double>> values) const override;

    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_IntegrationRule & ir,
                                             SliceMatrix<> coefs,
                                             BareSliceMatrix<SIMD<double>> values) const override;
    
    HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, SliceMatrix<> coefs, BareSliceMatrix<> values) const override;

    HD NGS_DLL_HEADER virtual void EvaluateTrans (const IntegrationRule & ir, 
                                                  BareSliceVector<> vals, 
                                                  BareSliceVector<double> coefs) const override;
    
    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_IntegrationRule & ir,
                                             BareVector<SIMD<double>> values,
                                             BareSliceVector<> coefs) const override;

    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_IntegrationRule & ir,
                                             BareSliceMatrix<SIMD<double>> values,
                                             SliceMatrix<> coefs) const override; 

    HD NGS_DLL_HEADER virtual Vec<DIM> EvaluateGrad (const IntegrationPoint & ip, 
                                                     BareSliceVector<> x) const override;

    HD NGS_DLL_HEADER virtual void EvaluateGrad (const IntegrationRule & ir, 
                                                 BareSliceVector<double> coefs, 
                                                 BareSliceMatrix<> vals) const override;

    HD NGS_DLL_HEADER virtual void EvaluateGrad (const SIMD_BaseMappedIntegrationRule & ir,
                                                 BareSliceVector<> coefs,
                                                 BareSliceMatrix<SIMD<double>> values) const override;

    HD NGS_DLL_HEADER virtual void EvaluateGrad (const SIMD_IntegrationRule & ir,
                                                 BareSliceVector<> coefs,
                                                 BareSliceMatrix<SIMD<double>> values) const override;

    HD NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, 
                                                      BareSliceMatrix<> vals, 
                                                      BareSliceVector<double> coefs) const override;

    HD NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, 
                                                      SliceMatrix<> values, 
                                                      SliceMatrix<> coefs) const override;

    HD NGS_DLL_HEADER virtual void AddGradTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                                 BareSliceMatrix<SIMD<double>> values,
                                                 BareSliceVector<> coefs) const override;

    HD NGS_DLL_HEADER virtual void AddGradTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                                 BareSliceMatrix<SIMD<double>> values,
                                                 SliceMatrix<> coefs) const override;

    HD NGS_DLL_HEADER virtual void CalcMappedDShape (const BaseMappedIntegrationPoint & mip, 
                                                     BareSliceMatrix<> dshape) const override;

    HD NGS_DLL_HEADER virtual void CalcMappedDShape (const BaseMappedIntegrationRule & mip, 
                                                     BareSliceMatrix<> dshape) const override;

    HD NGS_DLL_HEADER 
    virtual void CalcMappedDShape (const SIMD_BaseMappedIntegrationRule & mir, 
                                   BareSliceMatrix<SIMD<double>> dshapes) const override;
    
    /// compute dshape, matrix: ndof x (spacedim spacedim)
    NGS_DLL_HEADER virtual void CalcDDShape (const IntegrationPoint & ip, 
                                             BareSliceMatrix<> ddshape) const override;

    NGS_DLL_HEADER virtual void CalcMappedDDShape (const BaseMappedIntegrationPoint & mip, 
                                                   BareSliceMatrix<> ddshape) const override;
    
    // NGS_DLL_HEADER virtual void GetPolOrders (FlatArray<PolOrder<DIM> > orders) const;
    
    NGS_DLL_HEADER virtual void AddDualTrans (const IntegrationRule & ir, BareSliceVector<double> values, BareSliceVector<> coefs) const override;
    NGS_DLL_HEADER virtual void AddDualTrans (const SIMD_IntegrationRule & ir, BareVector<SIMD<double>> values, BareSliceVector<> coefs) const override;    
#endif

    NGS_DLL_HEADER 
    virtual void CalcDualShape (const BaseMappedIntegrationPoint & mip, BareSliceVector<> shape) const override;

    
    NGS_DLL_HEADER virtual bool GetDiagDualityMassInverse (FlatVector<> diag) const override;
    
  protected:
    /*
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (Tx x[], TFA & shape) const
    {
      static_cast<const FEL*> (this) -> T_CalcShape (x, shape);
    }
    */
    
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (const TIP<DIM,Tx> & ip, TFA && shape) const
    {
      static_cast<const FEL*> (this) -> T_CalcShape (ip, shape);
    }

    template<typename Tx, typename TFA>  
    INLINE void T_CalcDualShape (const TIP<DIM,Tx> & ip, TFA & shape) const
    {
      throw Exception (string("T_CalcDualShape not implemented for element ")+typeid(*this).name());       
      // static_cast<const FEL*> (this) -> T_CalcDualShape (ip, shape);
    }

    bool GetDiagDualityMassInverse2 (FlatVector<> diag) const { return false; }

    /*
    void CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
    {
      // throw Exception (string("dual shape not implemented for element ")+typeid(*this).name());
      double imeas = 1.0/mip.GetMeasure();
      shape = 0.0;
      static_cast<const FEL*> (this)->        
        T_CalcDualShape (GetTIP<DIM>(mip.IP()), SBLambda ( [&](int j, double val) { shape(j) = imeas * val; }));
    }
    */
    
  };


  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  class T_ScalarFiniteElementFO : public T_ScalarFiniteElement<FEL,ET>
  {
  public:
    INLINE T_ScalarFiniteElementFO ()
    {
      this->ndof = NDOF; 
      this->order = ORDER; 
    }
    // virtual ~T_ScalarFiniteElementFO () { ; }
  };




  template <ELEMENT_TYPE ET>
  class ScalarDummyFE : public T_ScalarFiniteElement<ScalarDummyFE<ET>,ET>
  {
  public:
    INLINE NGS_DLL_HEADER ScalarDummyFE()
    {
      this->ndof= 0;
      this->order = 0;
    }
    HD NGS_DLL_HEADER virtual ~ScalarDummyFE() { ; }
    HD NGS_DLL_HEADER virtual void GetDiagMassMatrix (FlatVector<> mass) const { ; } 
    template<typename Tx, typename TFA>  
    INLINE static void T_CalcShape (TIP<ngfem::Dim(ET),Tx> ip, TFA & shape) 
    { ; }

    INLINE void CalcDualShape2 (const BaseMappedIntegrationPoint & mip, SliceVector<> shape) const
    { ; }
  };

  extern template class  T_ScalarFiniteElement<ScalarDummyFE<ET_POINT>,ET_POINT>;
  extern template class  T_ScalarFiniteElement<ScalarDummyFE<ET_SEGM>,ET_SEGM>;
  extern template class  T_ScalarFiniteElement<ScalarDummyFE<ET_TRIG>,ET_TRIG>;
  extern template class  T_ScalarFiniteElement<ScalarDummyFE<ET_QUAD>,ET_QUAD>;
  extern template class  T_ScalarFiniteElement<ScalarDummyFE<ET_TET>,ET_TET>;
  extern template class  T_ScalarFiniteElement<ScalarDummyFE<ET_PRISM>,ET_PRISM>;
  extern template class  T_ScalarFiniteElement<ScalarDummyFE<ET_PYRAMID>,ET_PYRAMID>;
  extern template class  T_ScalarFiniteElement<ScalarDummyFE<ET_HEX>,ET_HEXAMID>;
  extern template class  T_ScalarFiniteElement<ScalarDummyFE<ET_HEX>,ET_HEX>;

  extern template class  ScalarDummyFE<ET_POINT>;
  extern template class  ScalarDummyFE<ET_SEGM>;
  extern template class  ScalarDummyFE<ET_TRIG>;
  extern template class  ScalarDummyFE<ET_QUAD>;
  extern template class  ScalarDummyFE<ET_TET>;
  extern template class  ScalarDummyFE<ET_PRISM>;
  extern template class  ScalarDummyFE<ET_PYRAMID>;
  extern template class  ScalarDummyFE<ET_HEXAMID>;
  extern template class  ScalarDummyFE<ET_HEX>;

}


#endif
