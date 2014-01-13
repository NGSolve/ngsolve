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

    T_ScalarFiniteElement () { ; }
    // virtual ~T_ScalarFiniteElement() { ; }

    virtual ELEMENT_TYPE ElementType() const { return ET; }
    
    NGS_DLL_HEADER virtual void CalcShape (const IntegrationPoint & ip, 
			    SliceVector<> shape) const;
    
    NGS_DLL_HEADER virtual double Evaluate (const IntegrationPoint & ip, FlatVector<double> x) const;
    
    NGS_DLL_HEADER virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const;
    
    NGS_DLL_HEADER virtual void EvaluateTrans (const IntegrationRule & ir, FlatVector<> vals, FlatVector<double> coefs) const;
    NGS_DLL_HEADER virtual void EvaluateGrad (const IntegrationRule & ir, FlatVector<double> coefs, FlatMatrixFixWidth<DIM> vals) const;

    NGS_DLL_HEADER virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> vals, FlatVector<double> coefs) const;

    NGS_DLL_HEADER virtual void CalcDShape (const IntegrationPoint & ip, 
			     SliceMatrix<> dshape) const;
/*    
    virtual void CalcDShape (const IntegrationPoint & ip, 
			     const std::function<void(int,Vec<DIM>)> & callback) const;
                       */
    /*
    virtual void CalcMappedDShape (const MappedIntegrationPoint<DIM,DIM> & mip, 
				   FlatMatrixFixWidth<DIM> dshape) const;
    */
    NGS_DLL_HEADER virtual void CalcMappedDShape (const MappedIntegrationPoint<DIM,DIM> & mip, 
				   SliceMatrix<> dshape) const;

    NGS_DLL_HEADER virtual void CalcMappedDShape (const MappedIntegrationRule<DIM,DIM> & mip, 
				   SliceMatrix<> dshape) const;

    NGS_DLL_HEADER virtual void GetPolOrders (FlatArray<PolOrder<DIM> > orders) const;
    
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
    T_ScalarFiniteElementFO ()
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
    NGS_DLL_HEADER ScalarDummyFE() { ; }
    NGS_DLL_HEADER virtual ~ScalarDummyFE() { ; }
    template<typename Tx, typename TFA>  
    static void T_CalcShape (Tx x[1], TFA & shape) 
    { ; }
  };


}


#endif
