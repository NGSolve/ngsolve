#ifdef PARALLEL

#include <solve.hpp>
#include <parallelngs.hpp>


using namespace ngsolve;

extern AutoPtr<PDE>  pde;
extern MeshAccess * ma;

extern "C" void NGS_ParallelRun ( const string & message );

namespace ngparallel
{
  MPI_Comm ngs_comm;
}


void * SolveBVP2(void *)
{
  if (pde && pde->IsGood())
    pde->SolveBVP();

  Ng_SetRunning (0); 
  return NULL;
}


void NGS_ParallelRun ( const string & message )
{

#ifdef _OPENMP
  omp_set_num_threads (1);
#endif


  if ( message == "ngs_pdefile" )
    {
      MPI_Comm_dup ( MPI_COMM_WORLD, &ngs_comm);
      // ngs_comm = MPI_COMM_WORLD;


      delete ma;
      ma = new MeshAccess;
      pde.Reset(new PDE ( *ma ));


      /*
      string pdefilename;
      MyMPI_Recv (pdefilename, 0);
      pde -> LoadPDE (pdefilename, 1, 0);
      */

      // transfer file contents, not filename
      string pdefiledata;
      MyMPI_Recv (pdefiledata, 0);

      istringstream pdefile (pdefiledata);
      pde -> LoadPDE (pdefile, 1, 0);
    } 

  else if ( message == "ngs_solvepde" )
    {
      RunParallel (SolveBVP2, NULL);
    }

  return;
}


void Parallel_Exit ()
{
  ;
}


#endif
