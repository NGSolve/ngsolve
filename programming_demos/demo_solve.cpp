#include <solve.hpp>
using namespace ngsolve;

int main(int argc, char ** argv)
{
  if (argc < 2)
    {
      cout << "Usage:  demo_solve filename" << endl;
      exit(1);
    }

  MyMPI mympi(argc, argv);

  MeshAccess ma;
  PDE pde(ma);

  pde.LoadPDE (argv[1]); // "d1_square.pde");

  pde.Solve();

  return 0;
}
