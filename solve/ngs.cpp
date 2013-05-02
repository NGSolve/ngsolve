#include <solve.hpp>

namespace netgen
{
  int h_argc;
  char ** h_argv;
}

int main(int argc, char ** argv)
{
  if (argc < 2)
    {
      std::cout << "Usage:  ngs filename" << std::endl;
      exit(1);
    }

  netgen::h_argc = argc;
  netgen::h_argv = argv;


  ngsolve::MyMPI mympi(argc, argv);

  ngsolve::PDE pde; 

  try
    {
      pde.LoadPDE (argv[argc-1]); 
      pde.Solve();
    }

  catch(ngstd::Exception & e)
    {
      std::cout << "Caught exception: " << std::endl
                << e.What() << std::endl;
    };

  return 0;
}
