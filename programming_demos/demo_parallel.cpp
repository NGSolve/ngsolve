// ng-soft header files
#include <solve.hpp>
using namespace ngstd;

int main (int argc, char ** argv)
{
  MyMPI mympi(argc, argv);
  
  if (MyMPI_GetNTasks() < 2)
    {
      cout << "run with 'mpirun -np 2 demo_parallel'" << endl;
      exit (1);
    }

  int id = MyMPI_GetId();
  cout << "hi from process " << id << endl;

  if (id == 0)
    {
      Array<bool> sa(5);
      sa = 1;
      MyMPI_Send (sa, 1);
    }
  if (id == 1)
    {
      Array<bool> ra(10);
      ra = 0;
      MyMPI_Recv (ra.Range(2,7), 0);
      cout << "got array" << endl << ra << endl;
    }
  
  return 0;
}
