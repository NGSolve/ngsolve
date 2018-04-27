/*
  
Solver for a hyperbolic pde

Solves

M d^2u/dt^2  +  A u = f

by the Newmark method


Please include this file to the src files given in netgen/ngsolve/Makefile
*/



#include <solve.hpp>

using namespace ngsolve;



/*
  Every solver is a class derived from the class NumProc.
  It collects objects (such as bilinear-forms, gridfunctions) and parameters.
*/

class NumProcHyperbolic : public NumProc
{
protected:
  // bilinear-form for the stiffness matrix
  shared_ptr<BilinearForm> bfa;
  // bilinear-form for the mass-matrix
  shared_ptr<BilinearForm> bfm;
  // linear-form providing the right hand side
  shared_ptr<LinearForm> lff;
  // solution vector
  shared_ptr<GridFunction> gfu;

  // time step
  double dt;
  // total time
  double tend;

public:
    
  /*
    In the constructor, the solver class gets the flags from the pde - input file.
    the PDE class apde contains all bilinear-forms, etc...
  */
  NumProcHyperbolic (shared_ptr<PDE> apde, const Flags & flags)
    : NumProc (apde)
  {
    // in the input-file, you specify the bilinear-forms for the stiffness and for the mass-term
    // like  "-bilinearforma=k". Default arguments are 'a' and 'm'

    bfa = apde->GetBilinearForm (flags.GetStringFlag ("bilinearforma", "a"));
    bfm = apde->GetBilinearForm (flags.GetStringFlag ("bilinearformm", "m"));
    lff = apde->GetLinearForm (flags.GetStringFlag ("linearform", "f"));
    gfu = apde->GetGridFunction (flags.GetStringFlag ("gridfunction", "u"));

    dt = flags.GetNumFlag ("dt", 0.001);
    tend = flags.GetNumFlag ("tend", 1);
  }

  virtual ~NumProcHyperbolic() 
  { ; }


  // solve at one level
  virtual void Do(LocalHeap & lh)
  {
    cout << "solve hyperbolic pde" << endl;
      
    // reference to the matrices provided by the bi-forms.
    // will be of type SparseSymmetricMatrix<double> for scalar problems

    const BaseMatrix & mata = bfa->GetMatrix();
    const BaseMatrix & matm = bfm->GetMatrix();
    const BaseVector & vecf = lff->GetVector();
    BaseVector & vecu = gfu->GetVector();

    // creates a matrix of same type:
    auto summat = matm.CreateMatrix();

    // create additional vectors:
    auto d = vecu.CreateVector();
    auto w = vecu.CreateVector();
    auto v = vecu.CreateVector();
    auto a = vecu.CreateVector();
    auto anew = vecu.CreateVector();


    // matrices matm and mata have the same memory layout. The arrays of values 
    // can be accessed and manipulated as vectors:

    summat->AsVector() = matm.AsVector() + (dt*dt/4) * mata.AsVector();

    // A sparse matrix can compute a sparse factorization. One has to cast to a sparse matrix:
    // auto invmat = * dynamic_cast<BaseSparseMatrix&> (*summat) . InverseMatrix();
    auto & invmat = * summat -> InverseMatrix();


    // Newmark method
    vecu = 0;
    v = 0;
    a = 0;

    for (double t = 0; t <= tend; t += dt)
      {
	cout << "t = " << t << endl;

	*w = vecu + dt * *v + (dt*dt/4) * *a;

	double fac = (t < 1) ? 1 : 0;
	*d = fac * vecf - mata * *w;
	*anew = invmat * *d;
	
	vecu += dt * *v + (dt*dt/4) * *a + (dt*dt/4) * *anew;
	*v += (dt/2) * *a + (dt/2) * *anew;

	*a = *anew;

	// update visualization
	Ng_Redraw ();	  
      }
  }






  virtual string GetClassName () const
  {
    return "Hyperbolic Solver (Demo)";
  }

  virtual void PrintReport (ostream & ost) const
  {
    ost << GetClassName() << endl
	<< "Bilinear-form A = " << bfa->GetName() << endl
	<< "Bilinear-form M = " << bfm->GetName() << endl
	<< "Linear-form     = " << lff->GetName() << endl
	<< "Gridfunction    = " << gfu->GetName() << endl
	<< "dt              = " << dt << endl
	<< "tend            = " << tend << endl;
  }

  ///
  static void PrintDoc (ostream & ost)
  {
    ost << 
      "\n\nNumproc Hyperbolic:\n" \
      "------------------\n" \
      "Solves a hyperbolic partial differential equation by an implicit Euler method\n\n" \
      "Required flags:\n" 
      "-bilinearforma=<bfname>\n" 
      "    bilinear-form providing the stiffness matrix\n" \
      "-bilinearformm=<bfname>\n" 
      "    bilinear-form providing the mass matrix\n" \
      "-linearform=<lfname>\n" \
      "    linear-form providing the right hand side\n" \
      "-gridfunction=<gfname>\n" \
      "    grid-function to store the solution vector\n" 
      "-dt=<value>\n"
      "    time step\n"
      "-tend=<value>\n"
      "    total time\n"
	<< endl;
  }
};


static RegisterNumProc<NumProcHyperbolic> nphyper("hyperbolic");

  
