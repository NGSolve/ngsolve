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
    shared_ptr<CoefficientFunction> coef_lambda;
  public:
    MyLaplaceIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : coef_lambda(coeffs[0])
    { ; }

    virtual string Name () const { return "MyLaplace"; }

    virtual int DimElement () const { return 2; }
    virtual int DimSpace () const { return 2; }


    // it is not a boundary integral (but a domain integral)
    virtual bool BoundaryForm () const { return false; }

    // Calculates the element matrix
    virtual void
    CalcElementMatrix (const FiniteElement & fel,
                       const ElementTransformation & eltrans, 
                       FlatMatrix<double> elmat,
                       LocalHeap & lh) const;


    // flux postprocessing 
    virtual int DimFlux () const { return 2; }

    virtual void
    CalcFlux (const FiniteElement & fel,
	      const BaseMappedIntegrationPoint & bsip,
              FlatVector<double> elx, 
	      FlatVector<double> flux,
	      bool applyd,
	      LocalHeap & lh) const;
  };



  // integrator for \int f v dx
  class MySourceIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_f;
  public:
    MySourceIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : coef_f(coeffs[0])
    { ; }

    virtual string Name () const { return "MySource"; }

    virtual bool BoundaryForm () const { return false; }

    // Calculates the right hand side element vector
    virtual void
    CalcElementVector (const FiniteElement & fel,
		       const ElementTransformation & eltrans, 
		       FlatVector<double> elvec,
		       LocalHeap & lh) const;
  };


}
#endif
