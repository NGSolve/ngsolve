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

  cout << "id = " << MyMPI_GetId() << ", ma.getne = " << ma.GetNE() << endl;
  int ntasks = MyMPI_GetNTasks();

  Array<int> nsend(ntasks);
  nsend = 3;

  Table<int> dist_data(nsend), recv_data(nsend);
  for (int i = 0; i < ntasks; i++)
    dist_data[i] = MyMPI_GetId();

  MyMPI_Barrier();
  cout << "try something" << endl;

  {
    MPI_Status status;
    int flag;
    MPI_Iprobe (0, MPI_TAG_SOLVE, MPI_COMM_WORLD, &flag, &status);
    cout << "iprove gives " << flag << endl;
  }

  char ch;
  cin >> ch;
  MyMPI_Barrier();

  MPI_Comm comm = ngs_comm;
  Array<MPI_Request> requests;

  for (int i = 0; i < ntasks; i++)
    {
      if (i == MyMPI_GetId()) continue;

      if (nsend[i])
	requests.Append (MyMPI_ISend (dist_data[i], i, MPI_TAG_SOLVE, comm));
      if (nsend[i])
	requests.Append (MyMPI_IRecv (recv_data[i], i, MPI_TAG_SOLVE, comm));
    }

  if (requests.Size())
    MyMPI_WaitAll (requests);


  cout << "success" << endl;

  pde.Solve();

#ifdef PARALLEL
  MPI_Finalize ();
#endif

  return 0;
}
