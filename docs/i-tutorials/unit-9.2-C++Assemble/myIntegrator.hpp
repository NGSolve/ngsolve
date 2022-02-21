#ifndef FILE_MYINTEGRATOR_HPP
#define FILE_MYINTEGRATOR_HPP


/*
  
My own simple integrators for the Poisson Equation

*/


namespace ngfem
{


  // integrator for \int f v dx
  class MySourceIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_f;
  public:
    MySourceIntegrator (shared_ptr<CoefficientFunction> coef)
      : coef_f(coef) { ; }

    string Name () const override { return "MySource"; }

    VorB VB() const override { return VOL; }
    
    // Calculates the element source vector
    void CalcElementVector (const FiniteElement & fel,
                            const ElementTransformation & eltrans, 
                            FlatVector<double> elvec,
                            LocalHeap & lh) const override;
  };




  
  // integrator for \int \lambda(x) \nabla u \nabla v dx
  class MyLaplaceIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef_lambda;
  public:
    MyLaplaceIntegrator (shared_ptr<CoefficientFunction> coef)
      : coef_lambda(coef) { ; }

    string Name () const override { return "MyLaplace"; }

    xbool IsSymmetric () const override { return true; }
    VorB VB() const override { return VOL; }

    // Calculates the element matrix
    void CalcElementMatrix (const FiniteElement & fel,
                            const ElementTransformation & eltrans, 
                            FlatMatrix<double> elmat,
                            LocalHeap & lh) const override;
  };

}

#endif
