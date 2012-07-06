#include <solve.hpp>
using namespace ngsolve;

extern void ParallelRun();

int main(int argc, char ** argv)
{
#ifdef PARALLEL
  MPI_Init (&argc, &argv);
  ngs_comm = MPI_COMM_WORLD;
  NGSOStream::SetGlobalActive (MyMPI_GetId() == 0);
#endif 

  MeshAccess ma;
  PDE pde(ma);
  pde.LoadPDE ("d1_square.pde");

  pde.Solve();

#ifdef PARALLEL
  MPI_Finalize ();
#endif

  return 0;
}
