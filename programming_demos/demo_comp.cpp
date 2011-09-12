/*

  NGSolve finite element demo

*/


#include <iostream>
#include <sstream>
#include <iomanip>

#include <time.h>
#include <math.h>
#include <complex>



// ng-soft header files
#include <comp.hpp>


using namespace std;
using namespace ngstd;
using namespace ngbla;
using namespace ngfem;
using namespace ngcomp;


namespace netgen {
  ostream * testout = &cout;
}



int main ()
{
  Ng_LoadGeometry ("cube.geo");
  Ng_LoadMesh ("cube.vol");
  
  LocalHeap lh(10000000, "main heap");
  MeshAccess ma;

  Flags fesflags;
  fesflags.SetFlag ("order", 2);
  H1HighOrderFESpace fes(ma, fesflags);
  
  Flags uflags;
  T_GridFunction<double> gfu (fes, "gfu", uflags);

  Flags aflags;
  T_BilinearFormSymmetric<double> bfa(fes, "bfa", aflags);

  bfa.AddIntegrator (new LaplaceIntegrator<3> (new ConstantCoefficientFunction(1)));
  bfa.AddIntegrator (new RobinIntegrator<3> (new ConstantCoefficientFunction(1)));




  Flags fflags;
  T_LinearForm<double> lff(fes, "lff", fflags);

  Array<EvalFunction*> asource(1);
  asource[0] = new EvalFunction ("sin(x)*y");
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
