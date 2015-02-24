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
  
  MyNonlinearIntegrator (Array<CoefficientFunction*> & coeffs)
  { ; }
  
  virtual string Name () const { return "MyNonlinear"; }

  virtual bool BoundaryForm () const { return 0; }
  
  virtual void
  CalcElementMatrix (const FiniteElement & fel,
		     const ElementTransformation & eltrans, 
		     FlatMatrix<double> & elmat,
		     LocalHeap & lh) const
  {
    FlatVector<> elveclin(fel.GetNDof(), lh);
    elveclin = 0;
    CalcLinearizedElementMatrix (fel, eltrans, elveclin, elmat, lh);
  }


  // compute the functional
  virtual double Energy (const FiniteElement & bfel, 
			 const ElementTransformation & eltrans, 
			 const FlatVector<double> & elx, 
			 LocalHeap & lh) const
  {
    const ScalarFiniteElement<2> & fel = static_cast<const ScalarFiniteElement<2>&> (bfel);
    int ndof = fel.GetNDof();

    FlatVector<> shape(ndof, lh);
    IntegrationRule ir(fel.ElementType(), 2*fel.Order());

    double energy = 0;

    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        HeapReset hr(lh);
        
        MappedIntegrationPoint<2,2> mip(ir[i], eltrans); 

        fel.CalcShape (ir[i], shape);

        double ui = InnerProduct (shape, elx);
        double phi = pow(ui,4);
        double fac = fabs (mip.GetJacobiDet()) * mip.IP().Weight();
        
        energy += fac * phi;
      }
    return energy;
  }


  // compute the gradient at a given point
  virtual void 
  ApplyElementMatrix (const FiniteElement & bfel, 
		      const ElementTransformation & eltrans, 
		      const FlatVector<double> & elx,    // element vector
		      FlatVector<double> & ely,          // element gradient
		      void * precomputed,
		      LocalHeap & lh) const
  {
    const ScalarFiniteElement<2> & fel = static_cast<const ScalarFiniteElement<2>&> (bfel);
    int ndof = fel.GetNDof();

    FlatVector<> shape(ndof, lh);
    IntegrationRule ir(fel.ElementType(), 2*fel.Order());

    ely = 0;

    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        HeapReset hr(lh);
        
        MappedIntegrationPoint<2,2> mip(ir[i], eltrans); 

        fel.CalcShape (ir[i], shape);

        double ui = InnerProduct (shape, elx);
        double phiprime = 4 * pow(ui,3);

        double fac = fabs (mip.GetJacobiDet()) * mip.IP().Weight();
        
        ely += fac * phiprime * shape;
      }
  }

  // compute the Hesse Matrix at point elveclin
  virtual void
  CalcLinearizedElementMatrix (const FiniteElement & bfel,
			       const ElementTransformation & eltrans,
			       FlatVector<double> & elveclin,
			       FlatMatrix<double> & elmat,
			       LocalHeap & lh) const
  {
    const ScalarFiniteElement<2> & fel = static_cast<const ScalarFiniteElement<2>&> (bfel);
    int ndof = fel.GetNDof();

    elmat.AssignMemory (ndof, ndof, lh);
    elmat = 0;
    

    FlatVector<> shape(ndof, lh);
    IntegrationRule ir(fel.ElementType(), 2*fel.Order());

    for (int i = 0 ; i < ir.GetNIP(); i++)
      {
        HeapReset hr(lh);
        MappedIntegrationPoint<2,2> mip(ir[i], eltrans); 

        fel.CalcShape (ir[i], shape);

        double uilin = InnerProduct(shape, elveclin);
        double phiprimeprime = 12 * pow(uilin,2);

        double fac = fabs (mip.GetJacobiDet()) * mip.IP().Weight();

        elmat += (fac * phiprimeprime) * shape * Trans(shape);
      }
  }
};





class NumProcNonlinearSolve : public NumProc
{
protected:
  BilinearForm * bfa;
  LinearForm * lff;
  GridFunction * gfu;

  int maxit;

public:
  NumProcNonlinearSolve (PDE & apde, const Flags & flags)
    : NumProc (apde)
  { 
    bfa = pde.GetBilinearForm (flags.GetStringFlag ("bilinearforma", "a"));
    lff = pde.GetLinearForm (flags.GetStringFlag ("linearform", "f"));
    gfu = pde.GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));

    maxit = int ( flags.GetNumFlag ("maxit", 30));
  }

  virtual ~NumProcNonlinearSolve()
  { ; }


  virtual void Do(LocalHeap & lh)
  {
    cout << "nonlinmag solver called" << endl;

    BaseVector & vecu = gfu->GetVector();
    const BaseVector & vecf = lff->GetVector();

    BaseVector & uold = *vecu.CreateVector();
    BaseVector & d = *vecu.CreateVector();
    BaseVector & w = *vecu.CreateVector();

    BilinearFormApplication applya(bfa);

    double err, errold, err0;
    double energy, energyold;

    d = vecf - applya * vecu;
    err0 = L2Norm(d);

    for (int i = 1; i <= maxit; i++)
      {
	cout << "newton it " << i << endl;
	
	bfa -> AssembleLinearization (vecu, lh);

	// bfa->GetMatrix().SetInverseType (SPARSECHOLESKY);
	BaseMatrix & inva = *bfa -> GetMatrix().InverseMatrix();

	d = vecf - applya * vecu;
	err = L2Norm(d);
	energy = bfa->Energy(vecu) - InnerProduct (vecf, vecu);

	cout << " err = " << err/err0;
	cout << " energy = " << energy << endl;

	errold = err;
	energyold = energy;

	w = inva * d;
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

	delete &inva;
	if (err < 1e-7*err0) break;
      }
  }
};



static RegisterNumProc<NumProcNonlinearSolve> npinit("nonlinearsolve");
static RegisterBilinearFormIntegrator<MyNonlinearIntegrator> initnl ("mynonlinear", 2, 0);
