#ifndef FILE_MYINTEGRATOR_HPP
#define FILE_MYINTEGRATOR_HPP

/*********************************************************************/
/* File:   myIntegrator.hpp                                          */
/* Author: Joachim Schoeberl                                         */
/* Date:   26. Apr. 2009                                             */
/*********************************************************************/

/*
  
My own simple integrators for the Poisson Equation

*/


namespace ngfem
{
  
  // integrator for \int \lambda(x) \nabla u \nabla v dx
  class MyLaplaceIntegrator : public BilinearFormIntegrator
  {
    CoefficientFunction * coef_lambda;
  public:
    MyLaplaceIntegrator (CoefficientFunction * acoef) 
      : coef_lambda(acoef)
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new MyLaplaceIntegrator (coeffs[0]);
    }

    virtual ~MyLaplaceIntegrator ()  { ; }

    virtual string Name () const { return "MyLaplace"; }

    // components of flux
    virtual int DimFlux () const { return 2; }

    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }

    // Calculates the element matrix
    virtual void
    AssembleElementMatrix (const FiniteElement & fel,
			   const ElementTransformation & eltrans, 
			   FlatMatrix<double> & elmat,
			   LocalHeap & lh) const;

  };



  // integrator for \int f v dx
  class MySourceIntegrator : public LinearFormIntegrator
  {
    CoefficientFunction * coef_f;
  public:
    MySourceIntegrator (CoefficientFunction * acoef) 
      : coef_f(acoef)
    { ; }

    static Integrator * Create (Array<CoefficientFunction*> & coeffs)
    {
      return new MySourceIntegrator (coeffs[0]);
    }

    virtual ~MySourceIntegrator ()  { ; }

    virtual string Name () const { return "MySource"; }

    virtual bool BoundaryForm () const { return false; }

    // Calculates the right hand side element vector
    virtual void
    AssembleElementVector (const FiniteElement & fel,
			   const ElementTransformation & eltrans, 
			   FlatVector<double> & elvec,
			   LocalHeap & lh) const;
  };


}
#endif
