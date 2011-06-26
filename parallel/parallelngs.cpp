#include <solve.hpp>
#include <parallelngs.hpp>

using namespace ngsolve;


#ifdef PARALLEL

extern AutoPtr<PDE>  pde;
extern MeshAccess * ma;


extern "C" void NGS_ParallelRun ( const string & message );

void NGS_ParallelRun ( const string & message )
{
  if ( message == "ngs_pdefile" )
    {
      string pdefilename;

      MyMPI_Recv ( pdefilename, 0);
      
      if ( ma ) delete ma;
      ma = new MeshAccess;

      pde.Reset(new PDE ( *ma ));
      pde -> LoadPDE (pdefilename, 1, 0);
    } 

  else if ( message == "ngs_solvepde" )
    {
      try
	{
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
	}
#endif
      catch (ngstd::Exception & e)
	{
	  cerr << "\n\ncaught Exception in SolveBVP:\n"
	       << e.What() << "\n\n";
	  pde->SetGood (false);
	}
    }

  return;
}


void Parallel_Exit ()
{
  ;
}


#endif
