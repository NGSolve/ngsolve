#ifndef FILE_TSCALARFE
#define FILE_TSCALARFE

/*********************************************************************/
/* File:   tscalarfe.hpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

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
    enum { DIM = ET_trait<ET>::DIM };

    // using BASE::eltype;
    using BASE::ndof;
    using BASE::order;

    INLINE T_ScalarFiniteElement () { ; }
    // virtual ~T_ScalarFiniteElement() { ; }

    HD virtual ELEMENT_TYPE ElementType() const { return ET; }
    
    HD NGS_DLL_HEADER virtual void CalcShape (const IntegrationPoint & ip, 
					      SliceVector<> shape) const;

    HD NGS_DLL_HEADER virtual void CalcDShape (const IntegrationPoint & ip, 
					       SliceMatrix<> dshape) const;

    
#ifndef FASTCOMPILE
    HD NGS_DLL_HEADER virtual void CalcShape (const IntegrationRule & ir, 
                                              SliceMatrix<> shape) const;
    
    HD NGS_DLL_HEADER virtual double Evaluate (const IntegrationPoint & ip, 
					       SliceVector<double> x) const;
    
    HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, 
					     SliceVector<double> coefs, 
					     FlatVector<double> vals) const;

    HD NGS_DLL_HEADER virtual void Evaluate (const SIMD_IntegrationRule & ir,
                                             BareSliceVector<> coefs,
                                             ABareVector<double> values) const;

    HD NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, SliceMatrix<> coefs, SliceMatrix<> values) const;

    HD NGS_DLL_HEADER virtual void EvaluateTrans (const IntegrationRule & ir, 
                                                  FlatVector<> vals, 
                                                  SliceVector<double> coefs) const;
    
    HD NGS_DLL_HEADER virtual void AddTrans (const SIMD_IntegrationRule & ir,
                                             ABareVector<double> values,
                                             BareSliceVector<> coefs) const;
    
    HD NGS_DLL_HEADER virtual Vec<DIM> EvaluateGrad (const IntegrationPoint & ip, 
                                                     SliceVector<> x) const;

    HD NGS_DLL_HEADER virtual void EvaluateGrad (const IntegrationRule & ir, 
                                                 SliceVector<double> coefs, 
                                                 FlatMatrixFixWidth<DIM> vals) const;

    HD NGS_DLL_HEADER virtual void EvaluateGrad (const SIMD_BaseMappedIntegrationRule & ir,
                                                 BareSliceVector<> coefs,
                                                 ABareMatrix<double> values) const;

    HD NGS_DLL_HEADER virtual void EvaluateGrad (const SIMD_IntegrationRule & ir,
                                                 BareSliceVector<> coefs,
                                                 ABareMatrix<double> values) const;

    HD NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, 
                                                      FlatMatrixFixWidth<DIM> vals, 
                                                      SliceVector<double> coefs) const;

    HD NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, 
                                                      SliceMatrix<> values, 
                                                      SliceMatrix<> coefs) const;

    HD NGS_DLL_HEADER virtual void AddGradTrans (const SIMD_BaseMappedIntegrationRule & ir,
                                                 ABareMatrix<double> values,
                                                 BareSliceVector<> coefs) const;

/*    
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     const std::function<void(int,Vec<DIM>)> & callback) const;
                       */

    HD NGS_DLL_HEADER virtual void CalcMappedDShape (const MappedIntegrationPoint<DIM,DIM> & mip, 
                                                  SliceMatrix<> dshape) const;

    HD NGS_DLL_HEADER virtual void CalcMappedDShape (const MappedIntegrationRule<DIM,DIM> & mip, 
				   SliceMatrix<> dshape) const;

#endif

    // NGS_DLL_HEADER virtual void GetPolOrders (FlatArray<PolOrder<DIM> > orders) const;
    
  protected:
    template<typename Tx, typename TFA>  
    INLINE void T_CalcShape (Tx x[], TFA & shape) const
    {
      static_cast<const FEL*> (this) -> T_CalcShape (x, shape);
    }
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
  class ScalarDummyFE : public T_ScalarFiniteElementFO<ScalarDummyFE<ET>,ET,0,0>
  {
  public:
    INLINE NGS_DLL_HEADER ScalarDummyFE() { ; }
    HD NGS_DLL_HEADER virtual ~ScalarDummyFE() { ; }
    template<typename Tx, typename TFA>  
    INLINE static void T_CalcShape (Tx x[1], TFA & shape) 
    { ; }
  };
}


namespace ngbla
{

  template <int DIM, typename SCAL = double>
  class AD2Vec : public MatExpr<AD2Vec<DIM,SCAL> >
  {
    AutoDiff<DIM,SCAL> ad;
  public:
    INLINE AD2Vec (double d) : ad(d) { ; }
    INLINE AD2Vec (AutoDiff<DIM,SCAL> aad) : ad(aad) { ; }
    INLINE SCAL operator() (int i) const { return ad.DValue(i); }
    INLINE SCAL operator() (int i, int j) const { return ad.DValue(i); }
    INLINE AutoDiff<DIM,SCAL> Data() const { return ad; }

    INLINE int Size () const { return DIM; }
    INLINE int Height () const { return DIM; }
    INLINE int Width () const { return 1; }
  };

}



#endif
