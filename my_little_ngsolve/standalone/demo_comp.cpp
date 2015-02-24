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
  LocalHeap lh(10000000, "main heap");


  auto ma  = make_shared<MeshAccess> ("cube.vol");

  auto fes = make_shared<H1HighOrderFESpace> (ma, Flags({ "order=2" }));
  auto gfu = make_shared<T_GridFunction<double>> (fes);
  auto bfa = make_shared<T_BilinearFormSymmetric<double>> (fes, "bfa", Flags({ "symmetric", "printelmat" }));

  bfa->AddIntegrator (make_shared<LaplaceIntegrator<3>> (make_shared<ConstantCoefficientFunction>(1)));
  bfa->AddIntegrator (make_shared<RobinIntegrator<3>> (make_shared<ConstantCoefficientFunction>(1)));

  auto lff = make_shared<T_LinearForm<double>>(fes, "lff", Flags());

  Array<shared_ptr<EvalFunction>> asource(1);
  asource[0] = make_shared<EvalFunction> ("sin(x)*y");
  // asource[0]->Print(cout);

  lff->AddIntegrator (make_shared<SourceIntegrator<3>> (make_shared<DomainVariableCoefficientFunction>(asource)));

  fes->Update(lh);
  fes->FinalizeUpdate(lh);

  gfu->Update();
  bfa->Assemble(lh);
  lff->Assemble(lh);

  BaseMatrix & mata = bfa->GetMatrix();
  BaseVector & vecf = lff->GetVector();
  BaseVector & vecu = gfu->GetVector();

  
  BaseMatrix * mat = &mata;

#ifdef PARALLEL
  const ParallelMatrix * pmat = dynamic_cast<const ParallelMatrix*> (mat);
  if (pmat) mat = &pmat->GetMatrix();
#endif

  shared_ptr<BaseMatrix> jacobi = 
    dynamic_cast<const BaseSparseMatrix&> (*mat).CreateJacobiPrecond();

    // BaseMatrix * jacobi = mata.InverseMatrix();

  
  CGSolver<double> inva (mata, *jacobi);
  inva.SetPrintRates();
  inva.SetMaxSteps(1000);

  vecu = inva * vecf;
}
