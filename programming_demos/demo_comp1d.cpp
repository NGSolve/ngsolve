/*

  NGSolve finite element demo

*/

// ng-soft header files
#include <comp.hpp>
using namespace ngcomp;

 

int main (int argc, char **argv)
{
#ifdef PARALLEL
  MPI_Init (&argc, &argv);
  ngs_comm = MPI_COMM_WORLD;
#endif 

  // Ng_LoadGeometry ("cube.geo");
  // Ng_LoadMesh ("cube.vol");
  
  LocalHeap lh(10000000, "main heap");
  MeshAccess ma;

  Array<double> dirichlet(1);
  dirichlet[0] = 1;
  H1HighOrderFESpace fes(ma, 
			 Flags() .SetFlag("order",2)
			 .SetFlag("print")
			 .SetFlag("dirichlet", dirichlet));

  // NodalFESpace fes(ma, Flags().SetFlag("order",1).SetFlag("print") );
  
  T_GridFunction<double> gfu (fes, "gfu", Flags());

  T_BilinearFormSymmetric<double> bfa(fes, "bfa", 
				      Flags().SetFlag("symmetric").SetFlag("print").SetFlag("printelmat"));

  bfa.AddIntegrator (new LaplaceIntegrator<1> (new ConstantCoefficientFunction(1)));
  // bfa.AddIntegrator (new RobinIntegrator<1> (new ConstantCoefficientFunction(1)));



  T_LinearForm<double> lff(fes, "lff", Flags());

  Array<EvalFunction*> asource(1);
  asource[0] = new EvalFunction ("sin(x)*y");
  asource[0]->Print(cout);

  lff.AddIntegrator (new SourceIntegrator<1> (new DomainVariableCoefficientFunction<3>(asource)));


  fes.Update(lh);
  fes.FinalizeUpdate(lh);
  gfu.Update();
  bfa.Assemble(lh);
  lff.Assemble(lh);


  BaseMatrix & mata = bfa.GetMatrix();
  BaseVector & vecf = lff.GetVector();
  BaseVector & vecu = gfu.GetVector();
  
  cout << "freedofs = " << *fes.GetFreeDofs() << endl;
  BaseMatrix * jacobi = dynamic_cast<const BaseSparseMatrix&> (mata).CreateJacobiPrecond(fes.GetFreeDofs());

  CGSolver<double> inva (mata, *jacobi);
  inva.SetPrintRates();
  inva.SetMaxSteps(1000);

  vecu = inva * vecf;


#ifdef PARALLEL
  MPI_Finalize ();
#endif

}
