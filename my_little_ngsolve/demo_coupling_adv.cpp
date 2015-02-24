/*
  

In this demo, we take the solution of one PDE, and use it as right
hand side of a second one.

For simplicity, we want to solve the first PDE,

-\Delta u = f

the second one is

-\Delta w = u


Similar problems occur when coupling different physical fields.

Please include this file to the src files given in netgen/ngsolve/Makefile
*/



#include <solve.hpp>
using namespace ngsolve;


/*
  DifferentialOperator can compute some function of u and its derivatives ...
*/
class MyDifferentialOperator : public DifferentialOperator
{
public:
  virtual int Dim() const { return 1; }
  virtual bool Boundary() const { return false; }
  virtual int DiffOrder() const { return 0; }

  virtual void
  CalcMatrix (const FiniteElement & fel,
	      const BaseMappedIntegrationPoint & mip,
	      FlatMatrix<double> mat, 
	      LocalHeap & lh) const 
  { cout << "never used ... uuuups" << endl; }
  
  virtual void Apply (const FiniteElement & fel,
		      const BaseMappedIntegrationPoint & mip,
		      FlatVector<double> x, 
		      FlatVector<double> flux,
		      LocalHeap & lh) const
    {
      double u = static_cast<const ScalarFiniteElement<2>&> (fel).Evaluate (mip.IP(), x);
      flux(0) = u*u*u;
    }

};


class NumProcCouplingDemoAdv : public NumProc
{
protected:
  shared_ptr<GridFunction> gfu;
  shared_ptr<LinearForm> lff;

public:
  /*
    In the constructor, the solver class gets the flags from the pde - input file.
    the PDE class apde constains all bilinear-forms, etc...
  */
  NumProcCouplingDemoAdv (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    lff = apde->GetLinearForm (flags.GetStringFlag ("linearform", "f"));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));

    auto coefu = 
      make_shared<GridFunctionCoefficientFunction> (gfu, make_shared<MyDifferentialOperator>());

    // add a source integrator to the linearform
    // the coefficient is  DiffOp (gfu)
    lff -> AddIntegrator (CreateLFI("source", 2, coefu));
  }
  
  virtual void Do (LocalHeap & lh) 
  {
    // re-assemble f with latest u
    lff -> Assemble (lh);
  }
};

static RegisterNumProc<NumProcCouplingDemoAdv> npinit1("democoupling_adv");

