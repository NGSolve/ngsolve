#include <solve.hpp>
#include <cstdlib>

namespace netgen
{
  int h_argc;
  char ** h_argv;
}

extern int dummy_bvp;

void printMeasurement( string name ) {
    int nr = NgProfiler::GetNr(name);
    cout << name << ": " << nr << endl;
    if(nr < 0) return;
    double time =  NgProfiler::GetTime(nr);
    string timer_name = name;

    for (int i=0; i< timer_name.length(); i++) {
        char & c = timer_name[i];
        if(c==' ' || c==':' || c=='-') {
            timer_name.erase(i,1);
            i--;
        }
    }
    if(time == 0) return;
    cout << "<DartMeasurement name=" << '"' << timer_name << '"' << endl;
    cout << "type=\"numeric/double\">" << time << "</DartMeasurement>" << endl;

    double mflops = NgProfiler::GetFlops(nr) / NgProfiler::GetTime(nr) * 1e-6; 
    if(mflops == 0) return;
    cout << "<DartMeasurement name=" << '"' << timer_name << "MFlops" << '"' << endl;
    cout << "type=\"numeric/double\">" << mflops << "</DartMeasurement>" << endl;
}

int main(int argc, char ** argv)
{
  int retcode = EXIT_SUCCESS;
  if (argc < 2)
    {
      std::cout << "Usage:  ngs filename" << std::endl;
      exit(EXIT_FAILURE);
    }

  netgen::h_argc = argc;
  netgen::h_argv = argv;
  dummy_bvp = 17;

  ngsolve::MyMPI mympi(argc, argv);
  cout << IM(1) << "CTEST_FULL_OUTPUT" << endl;

  try
    {
      string filename = argv[argc-1];
      auto pde = ngcomp::LoadPDE (filename);
      pde->Solve();
      printMeasurement( "Matrix assembling" );
      printMeasurement( "Solver - Total" );
      printMeasurement( "CG solver" );
      printMeasurement( "SparseMatrixSymmetric::MultAdd" );
      printMeasurement( "SparseMatrixSymmetric::MultAdd1" );
      printMeasurement( "SparseMatrixSymmetric::MultAdd2" );
      NgProfiler::Print (stdout);
    }

  catch(ngstd::Exception & e)
    {
      std::cout << "Caught exception: " << std::endl
                << e.What() << std::endl;
      retcode = EXIT_FAILURE;
    };


#ifdef PARALLEL
  int id;
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  char filename[100];
  sprintf (filename, "ngs.prof.%d", id);
  FILE *prof = fopen(filename,"w");
  NgProfiler::Print (prof);
  fclose(prof);
#endif


  return retcode;
}

