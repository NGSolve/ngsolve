/*

Solve the nonlinear problem

min_v   1/2 \| \nabla v \|^2  + \int phi(v) - f v 


with a convex function phi(v)

Example: phi(v) = v^4

*/


#include <solve.hpp>


using namespace ngfem;
using namespace ngsolve;


class MyNonlinearIntegrator : public BilinearFormIntegrator
{
public:
  
  MyNonlinearIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs)
  { ; }
  
  virtual string Name () const { return "MyNonlinear"; }

  virtual bool BoundaryForm () const { return 0; }
  
  virtual void
  CalcElementMatrix (const FiniteElement & fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> elmat,
		     LocalHeap & lh) const
  {
    FlatVector<> elveclin(fel.GetNDof(), lh);
    elveclin = 0;
    CalcLinearizedElementMatrix (fel, eltrans, elveclin, elmat, lh);
  }


  // compute the functional
  virtual double Energy (const FiniteElement & bfel, 
			 const ElementTransformation & eltrans, 
			 FlatVector<double> elx, 
			 LocalHeap & lh) const
  {
    HeapReset hr(lh);
    const ScalarFiniteElement<2> & fel = static_cast<const ScalarFiniteElement<2>&> (bfel);

    IntegrationRule ir(fel.ElementType(), 2*fel.Order());
    MappedIntegrationRule<2,2> mir(ir, eltrans, lh);

    FlatVector<> vals(ir.GetNIP(), lh);
    fel.Evaluate (ir, elx, vals);

    double energy = 0;
    for (int i = 0 ; i < ir.GetNIP(); i++)
      energy += mir[i].GetWeight() * pow(vals(i),4);
    return energy;
  }


  // compute the gradient at a given point
  virtual void 
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
                      FlatVector<double> elx,      // element vector
		      FlatVector<double> ely,      // element gradient
		      void * precomputed,
		      LocalHeap & lh) const
  {
    const ScalarFiniteElement<2> & fel = static_cast<const ScalarFiniteElement<2>&> (bfel);

    IntegrationRule ir(fel.ElementType(), 2*fel.Order());
    MappedIntegrationRule<2,2> mir(ir, eltrans, lh);

    FlatVector<> vals(ir.GetNIP(), lh);
    fel.Evaluate (ir, elx, vals);

    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        double phiprime = 4 * pow(vals(i),3);
        vals(i) = mir[i].GetWeight() * phiprime;
      }

    fel.EvaluateTrans (ir, vals, ely);
  }

  // compute the Hesse Matrix at point elveclin
  virtual void
  CalcLinearizedElementMatrix (const FiniteElement & bfel,
			       const ElementTransformation & eltrans,
			       FlatVector<double> elveclin,
			       FlatMatrix<double> elmat,
			       LocalHeap & lh) const
  {
    const ScalarFiniteElement<2> & fel = static_cast<const ScalarFiniteElement<2>&> (bfel);
    int ndof = fel.GetNDof();

    elmat = 0;
    
    FlatVector<> shape(ndof, lh);
    IntegrationRule ir(fel.ElementType(), 2*fel.Order());

    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        MappedIntegrationPoint<2,2> mip(ir[i], eltrans); 

        fel.CalcShape (ir[i], shape);

        double uilin = InnerProduct(shape, elveclin);
        double phiprimeprime = 12 * pow(uilin,2);

        double fac = mip.GetWeight();
        elmat += (fac * phiprimeprime) * shape * Trans(shape);
      }
  }
};





class NumProcNonlinearSolve : public NumProc
{
protected:
  shared_ptr<BilinearForm> bfa;
  shared_ptr<LinearForm> lff;
  shared_ptr<GridFunction> gfu;

  int maxit;

public:
  NumProcNonlinearSolve (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  { 
    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearforma", "a"));
    lff = apde->GetLinearForm (flags.GetStringFlag ("linearform", "f"));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));

    maxit = int ( flags.GetNumFlag ("maxit", 30));
  }

  virtual ~NumProcNonlinearSolve()
  { ; }


  virtual void Do(LocalHeap & lh)
  {
    cout << "nonlinmag solver called" << endl;

    BaseVector & vecu = gfu->GetVector();
    const BaseVector & vecf = lff->GetVector();

    auto uold = vecu.CreateVector();
    auto d = vecu.CreateVector();
    auto w = vecu.CreateVector();

    BilinearFormApplication applya(bfa);

    double err, err0;
    double energy, energyold;

    d = vecf - applya * vecu;
    err0 = L2Norm(*d);

    for (int i = 1; i <= maxit; i++)
      {
	cout << "newton it " << i << endl;
	
	bfa -> AssembleLinearization (vecu, lh);

	// bfa->GetMatrix().SetInverseType (SPARSECHOLESKY);
	auto inva = bfa -> GetMatrix().InverseMatrix();

	d = vecf - applya * vecu;
	err = L2Norm(*d);
	energy = bfa->Energy(vecu) - InnerProduct (vecf, vecu);

	cout << " err = " << err/err0;
	cout << " energy = " << energy << endl;

	energyold = energy;

	w = *inva * d;
	uold = vecu;
	int lin_its = 0;
	double tau = 1;

	do
	  {
	    vecu = uold + tau * w;
	    energy = bfa->Energy(vecu) - InnerProduct (vecf, vecu);

	    cout << "tau = " << tau
		 << " energy = " << energy << endl;

	    tau *= 0.5;
	  }
	while (energy > energyold && lin_its++ < 30 && err > 1e-7*err0);

	if (err < 1e-7*err0) break;
      }
  }
};



static RegisterNumProc<NumProcNonlinearSolve> npinit("nonlinearsolve");
static RegisterBilinearFormIntegrator<MyNonlinearIntegrator> initnl ("mynonlinear", 2, 0);
