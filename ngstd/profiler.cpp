/**************************************************************************/
/* File:   profiler.cpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2002                                                  */
/**************************************************************************/

#include <ngstd.hpp>
#ifdef PARALLEL
#include <mpi.h>
#endif

namespace ngstd
{
  using namespace ngstd;

  double NgProfiler::tottimes[SIZE];
  double NgProfiler::starttimes[SIZE];
  long int NgProfiler::counts[SIZE];
  double NgProfiler::flops[SIZE];
  double NgProfiler::loads[SIZE];
  double NgProfiler::stores[SIZE];
  string NgProfiler::names[SIZE];
  int NgProfiler::usedcounter[SIZE];
  

  NgProfiler :: NgProfiler()
  {
    for (int i = 0; i < SIZE; i++)
      {
	tottimes[i] = 0;
	usedcounter[i] = 0;
	flops[i] = 0;
      }

    total_timer = CreateTimer ("total CPU time");
    StartTimer (total_timer);
  }

  NgProfiler :: ~NgProfiler()
  {
    // StopTimer (total_timer);

    //ofstream prof;
    //prof.open("ng.prof");


    // ofstream-constructor may be called after STL-stuff is destructed,
    // which leads to an "order of destruction"-problem,
    // thus we use the C-variant:

    FILE *prof = fopen("ng.prof","w");
    Print (prof); 
    fclose(prof);
    if (getenv ("NGSPROFILE"))
      {
	char filename[100];

#ifdef PARALLEL
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sprintf (filename, "ngs.prof.%d", id);
#else
	sprintf (filename, "ngs.prof");
#endif
	
	if (id == 0) printf ("write profile to file ngs.prof\n"); 
	FILE *prof = fopen(filename,"w");
	Print (prof);
	fclose(prof);
      }
  }


  //   void NgProfiler :: Print (ostream & prof)
  //   {
  //     for (int i = 0; i < SIZE; i++)
  //       if (counts[i] != 0 || usedcounter[i] != 0)
  // 	{
  // 	  prof.setf (ios::fixed, ios::floatfield);
  // 	  prof.setf (ios::showpoint);

  // 	  prof << "job " << setw(3) << i 
  // 	       << " calls " << setw(8) << counts[i] 
  // 	       << ", time " << setprecision(2) << setw(6) << double(tottimes[i]) / CLOCKS_PER_SEC << " sec";
  // 	  if (flops[i]) prof << ", MFlops = " << flops[i] / (double(tottimes[i]) / CLOCKS_PER_SEC) * 1e-6;
  // 	  if (usedcounter[i]) prof << " " << names[i];
  // 	  prof << endl;
  // 	}
  //   }


  void NgProfiler :: Print (FILE * prof)
  {
#ifdef USE_TIMEOFDAY
    double fac = 1;
#else
    double fac = 1.0/CLOCKS_PER_SEC;
#endif

    for (int i = 0; i < SIZE; i++)
      if (counts[i] != 0 || usedcounter[i] != 0)
	{
	  // fprintf(prof,"job %3i calls %8i, time %6.2f sec",i,counts[i],double(tottimes[i]) / CLOCKS_PER_SEC);
          fprintf(prof,"job %3i calls %8li, time %6.4f sec",i,counts[i],tottimes[i]*fac);
	  if(flops[i])
	    fprintf(prof,", MFlops = %6.2f",flops[i] / (double(tottimes[i])*fac) * 1e-6);
	  if(loads[i])
	    fprintf(prof,", MLoads = %6.2f",loads[i] / (double(tottimes[i])*fac) * 1e-6);
	  if(stores[i])
	    fprintf(prof,", MStores = %6.2f",stores[i] / (double(tottimes[i])*fac) * 1e-6);
	  if(usedcounter[i])
	    fprintf(prof," %s",names[i].c_str());
	  fprintf(prof,"\n");
	}
  }


  int NgProfiler :: CreateTimer (const string & name)
  {
    for (int i = SIZE-1; i > 0; i--)
      if (!usedcounter[i])
	{
	  usedcounter[i] = 1;
	  names[i] = name;
	  return i;
	}
    throw Exception ("no more timer available");
    return -1;
  }


  NgProfiler prof;



#ifdef  VTRACE
  Timer * Timer::stack_top = NULL;
#endif


}
