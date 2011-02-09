#ifdef PARALLEL

#include <solve.hpp>
#include <comp.hpp>
#include <parallelngs.hpp>
#include <ngstd.hpp>

using namespace ngsolve;
using namespace ngparallel;


extern AutoPtr<PDE>  pde;
extern MeshAccess * ma;
extern ParallelMeshAccess * ngparallel::parallelma;




void NGS_ParallelRun ( const string & message )
{
  MPI_Status status;
  
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  
  if ( message == "ngs_pdefile" )
    {
      string pdefilename;
#ifdef SCALASCA
#pragma pomp inst begin (recvpdefile)
#endif
      ngparallel::MyMPI_Recv ( pdefilename, 0);

#ifdef SCALASCA
#pragma pomp inst end (recvpdefile)
#endif
      if ( ma ) delete ma;
      if ( parallelma ) delete parallelma;
      ma = new MeshAccess;
      pde.Reset(new PDE ( *ma ));
      parallelma = new ParallelMeshAccess ( *ma );
      pde -> SetParallelMeshAccess ( parallelma );
      pde -> LoadPDE (pdefilename, 1, 0);
    } 

  else if ( message == "ngs_solvepde" )
    {
      try
	{
	  cout << "parallel solve bvp" << endl;
	  pde -> SolveBVP();
	}
      catch (exception & e)
	{
	  cerr << "\n\ncaught exception in SolveBVP:\n " 
	       << typeid(e).name() << endl;
	  pde->SetGood (false);
	}
#ifdef WIN32
      catch (CException * e)
	{
	  TCHAR msg[255];
	  e->GetErrorMessage(msg, 255);
	  cerr << "\n\ncaught Exception in SolveBVP:\n"
	       << msg << "\n\n";
	  pde->SetGood (false);
	  // got_exception = true;
	}
#endif
      catch (ngstd::Exception & e)
	{
	  cerr << "\n\ncaught Exception in SolveBVP:\n"
	       << e.What() << "\n\n";
	  pde->SetGood (false);
	  // got_exception = true;
	}
    }


  (*testout) << message << " done! " << endl;

  return;
}



void Parallel_Exit ()
{
  delete parallelma;
}


#endif
