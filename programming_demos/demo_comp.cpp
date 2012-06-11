/*

  NGSolve finite element demo

*/

// ng-soft header files
#include <comp.hpp>
using namespace ngcomp;



int main ()
{
  Ng_LoadGeometry ("cube.geo");
  Ng_LoadMesh ("cube.vol");
  
  LocalHeap lh(10000000, "main heap");
  MeshAccess ma;

  H1HighOrderFESpace fes(ma, Flags().SetFlag("order",2) );
  
  T_GridFunction<double> gfu (fes, "gfu", Flags());

  T_BilinearFormSymmetric<double> bfa(fes, "bfa", 
				      Flags().SetFlag("symmetric"));

  bfa.AddIntegrator (new LaplaceIntegrator<3> (new ConstantCoefficientFunction(1)));
  bfa.AddIntegrator (new RobinIntegrator<3> (new ConstantCoefficientFunction(1)));



  T_LinearForm<double> lff(fes, "lff", Flags());

  Array<EvalFunction*> asource(1);
  asource[0] = new EvalFunction ("sin(x)*y");
  asource[0]->Print(cout);

  lff.AddIntegrator (new SourceIntegrator<3> (new DomainVariableCoefficientFunction<3>(asource)));


  fes.Update(lh);
  gfu.Update();
  bfa.Assemble(lh);
  lff.Assemble(lh);


  BaseMatrix & mata = bfa.GetMatrix();
  BaseVector & vecf = lff.GetVector();
  BaseVector & vecu = gfu.GetVector();
  
  BaseMatrix * jacobi = dynamic_cast<const BaseSparseMatrix&> (mata).CreateJacobiPrecond();

  CGSolver<double> inva (mata, *jacobi);
  inva.SetPrintRates();
  inva.SetMaxSteps(1000);

  vecu = inva * vecf;
}
