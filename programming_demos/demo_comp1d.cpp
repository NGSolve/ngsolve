/*

  NGSolve 1D finite element demo

*/

// ng-soft header files
#include <comp.hpp>
#include <ng_mesh1d.hpp>
using namespace ngcomp;

 

int main (int argc, char **argv)
{
  cout << "1D demo" << endl;

  MyMPI mympi (argc, argv);

  
  LocalHeap lh(100000, "main heap");

  Mesh1D ma(20);  // a mesh with 4 1D - elements

  cout << "elements: " << ma.GetNE() << endl;
  cout << "vertices: " << ma.GetNV() << endl;

  Array<double> dirichlet(1);
  dirichlet[0] = 1;
  H1HighOrderFESpace fes(ma, 
			 Flags() 
                         .SetFlag("order",3)
			 // .SetFlag("print")
			 .SetFlag("dirichlet", dirichlet)
			 );

  T_GridFunction<double> gfu (fes, "gfu", Flags());

  T_BilinearFormSymmetric<double> bfa(fes, "bfa", 
				      Flags()
				      .SetFlag("symmetric")
				      // .SetFlag("print")
				      // .SetFlag("printelmat")
				      );


  bfa.AddIntegrator (new LaplaceIntegrator<1> (new ConstantCoefficientFunction(1)));
  bfa.AddIntegrator (new RobinIntegrator<1> (new ConstantCoefficientFunction(1)));



  T_LinearForm<double> lff(fes, "lff", Flags());

  Array<EvalFunction*> asource(1);
  asource[0] = new EvalFunction ("sin(x)*y");

  lff.AddIntegrator (new SourceIntegrator<1> (new DomainVariableCoefficientFunction<3>(asource)));


  fes.Update(lh);
  fes.FinalizeUpdate(lh);
  gfu.Update();
  bfa.Assemble(lh);
  lff.Assemble(lh);


  BaseMatrix & mata = bfa.GetMatrix();
  BaseVector & vecf = lff.GetVector();
  BaseVector & vecu = gfu.GetVector();
  
  BaseMatrix * jacobi = dynamic_cast<const BaseSparseMatrix&> (mata).
    CreateJacobiPrecond(fes.GetFreeDofs());
  
  CGSolver<double> inva (mata, *jacobi);
  inva.SetPrintRates();
  inva.SetMaxSteps(100);

  cout << "CG solver:" << endl;
  vecu = inva * vecf;


  ofstream out("solution.out");
  GridFunctionCoefficientFunction cfu (gfu);
  for (int i = 0; i < ma.GetNE(); i++)
    {
      HeapReset hr(lh);
      const ElementTransformation & trafo = ma.GetTrafo(i, false, lh);
      for (double xi = 0; xi <= 1; xi += 0.1)
	{
	  IntegrationPoint ip(xi);
	  MappedIntegrationPoint<1,1> mip(ip, trafo);
	  out << mip(0) << " " << cfu.Evaluate (mip) << endl;
	}
    }
}
