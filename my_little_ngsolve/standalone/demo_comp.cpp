/*

  NGSolve finite element demo

*/

// ng-soft header files
#include <comp.hpp>
using namespace ngcomp;


int main (int argc, char **argv)
{
  MyMPI mympi(argc, argv);

  Ng_LoadGeometry ("cube.geo");

  MeshAccess ma;
  ma.LoadMesh ("cube.vol");
  
  LocalHeap lh(10000000, "main heap");

  H1HighOrderFESpace fes(ma, { "order=2" });
  
  T_GridFunction<double> gfu (fes);

  T_BilinearFormSymmetric<double> bfa(fes, "bfa", { "symmetric", "printelmat" });

  bfa.AddIntegrator (new LaplaceIntegrator<3> (new ConstantCoefficientFunction(1)));
  bfa.AddIntegrator (new RobinIntegrator<3> (new ConstantCoefficientFunction(1)));

  T_LinearForm<double> lff(fes, "lff", { } );

  Array<EvalFunction*> asource(1);
  asource[0] = new EvalFunction ("sin(x)*y");
  // asource[0]->Print(cout);

  lff.AddIntegrator (new SourceIntegrator<3> (new DomainVariableCoefficientFunction<3>(asource)));

  fes.Update(lh);
  fes.FinalizeUpdate(lh);

  gfu.Update();
  bfa.Assemble(lh);
  lff.Assemble(lh);

  BaseMatrix & mata = bfa.GetMatrix();
  BaseVector & vecf = lff.GetVector();
  BaseVector & vecu = gfu.GetVector();

  
  BaseMatrix * mat = &mata;

#ifdef PARALLEL
  const ParallelMatrix * pmat = dynamic_cast<const ParallelMatrix*> (mat);
  if (pmat) mat = &pmat->GetMatrix();
#endif

  BaseMatrix * jacobi = 
    dynamic_cast<const BaseSparseMatrix&> (*mat).CreateJacobiPrecond();

    // BaseMatrix * jacobi = mata.InverseMatrix();

  
  CGSolver<double> inva (mata, *jacobi);
  inva.SetPrintRates();
  inva.SetMaxSteps(1000);

  vecu = inva * vecf;
}
