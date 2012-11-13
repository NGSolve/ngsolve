#include <solve.hpp>


int main(int argc, char ** argv)
{
  if (argc < 2)
    {
      std::cout << "Usage:  ngs filename" << std::endl;
      exit(1);
    }

  ngsolve::MyMPI mympi(argc, argv);

  ngsolve::PDE pde; 

  pde.LoadPDE (argv[1]); 
  pde.Solve();

  return 0;
}
