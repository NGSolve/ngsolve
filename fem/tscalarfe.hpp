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

  template <class FEL, ELEMENT_TYPE ET, class BASE = ScalarFiniteElement<ET_trait<ET>::DIM> >
  class T_ScalarFiniteElement : public BASE   
  {
  public:
    enum { DIM = ET_trait<ET>::DIM };

    // using BASE::eltype;
    using BASE::ndof;
    using BASE::order;

    T_ScalarFiniteElement () { /* eltype = ET; */ }
    virtual ~T_ScalarFiniteElement() = 0;

    virtual ELEMENT_TYPE ElementType() const { return ET; }

    virtual void CalcShape (const IntegrationPoint & ip, 
			    FlatVector<> shape) const;

    virtual double Evaluate (const IntegrationPoint & ip, FlatVector<double> x) const;
    
    virtual void Evaluate (const IntegrationRule & ir, FlatVector<double> coefs, FlatVector<double> vals) const;

    virtual void EvaluateTrans (const IntegrationRule & ir, FlatVector<> vals, FlatVector<double> coefs) const;

    virtual void EvaluateGrad (const IntegrationRule & ir, FlatVector<double> coefs, FlatMatrixFixWidth<DIM> vals) const;

    virtual void EvaluateGradTrans (const IntegrationRule & ir, FlatMatrixFixWidth<DIM> vals, FlatVector<double> coefs) const;

    virtual void CalcDShape (const IntegrationPoint & ip, 
			     FlatMatrixFixWidth<DIM> dshape) const;

    virtual void 
    CalcMappedDShape (const MappedIntegrationPoint<DIM,DIM> & sip, 
                      FlatMatrixFixWidth<DIM> dshape) const;


    virtual void GetPolOrders (FlatArray<PolOrder<DIM> > orders) const;

  protected:
    template<typename Tx, typename TFA>  
    void T_CalcShape (Tx x[], TFA & shape) const
    {
      static_cast<const FEL*> (this) -> T_CalcShape (x, shape);
    }
  };


  template <class FEL, ELEMENT_TYPE ET, class BASE>
  inline T_ScalarFiniteElement<FEL,ET,BASE> :: 
  ~T_ScalarFiniteElement() { ; }



  template <class FEL, ELEMENT_TYPE ET, int NDOF, int ORDER>
  class T_ScalarFiniteElementFO : public T_ScalarFiniteElement<FEL,ET>
  {
  public:
    T_ScalarFiniteElementFO () 
    {
      this->ndof = NDOF; 
      this->order = ORDER; 
    }

    virtual ~T_ScalarFiniteElementFO ();
  };

}


#endif
