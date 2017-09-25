/**************************************************************************/
/* File:   profiler.cpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2002                                                  */
/**************************************************************************/

#include <ngstd.hpp>
/*
#ifdef PARALLEL
#include <mpi.h>
#endif
*/
namespace ngstd
{
  using namespace ngstd;

  std::chrono::time_point<std::chrono::system_clock> wall_time_start = std::chrono::system_clock::now();

  double NgProfiler::tottimes[SIZE];
  double NgProfiler::starttimes[SIZE];
  long int NgProfiler::counts[SIZE];
  double NgProfiler::flops[SIZE];
  double NgProfiler::loads[SIZE];
  double NgProfiler::stores[SIZE];
  string NgProfiler::names[SIZE];
  int NgProfiler::usedcounter[SIZE];
  string NgProfiler::filename;

  size_t dummy_thread_times[NgProfiler::SIZE];
  size_t * NgProfiler::thread_times = dummy_thread_times;
  size_t dummy_thread_flops[NgProfiler::SIZE];
  size_t * NgProfiler::thread_flops = dummy_thread_flops;

  NgProfiler :: NgProfiler()
  {
    for (int i = 0; i < SIZE; i++)
      {
	tottimes[i] = 0;
	usedcounter[i] = 0;
	flops[i] = 0;
      }

    // total_timer = CreateTimer ("total CPU time");
    // StartTimer (total_timer);
  }

  NgProfiler :: ~NgProfiler()
  {
    // StopTimer (total_timer);

    //ofstream prof;
    //prof.open("ng.prof");


    // ofstream-constructor may be called after STL-stuff is destructed,
    // which leads to an "order of destruction"-problem,
    // thus we use the C-variant:

    /*
    FILE *prof = fopen("ng.prof","w");
    Print (prof); 
    fclose(prof);
    */
    
    if (filename.length())
      {
	// printf ("write profile to file %s\n", filename.c_str());
	// cout << "write profile to file " << filename << endl;
	FILE *prof = fopen(filename.c_str(),"w");
	Print (prof);
	fclose(prof);
      }

    /*
    if (getenv ("NGSPROFILE"))
      {
	char filename[100];

#ifdef PARALLEL
	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	sprintf (filename, "ngs.prof.%d", id);
	printf ("write profiles to files ngs.prof.<id>\n"); 
#else
	sprintf (filename, "ngs.prof");
	printf ("write profile to file ngs.prof\n"); 
#endif
	
	FILE *prof = fopen(filename,"w");
	Print (prof);
	fclose(prof);
      }
    */
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
    static mutex createtimer_mutex;
    int nr = -1;
    {
      lock_guard<mutex> guard(createtimer_mutex);
      for (int i = SIZE-1; i > 0; i--)
	if (!usedcounter[i])
	  {
	    usedcounter[i] = 1;
	    names[i] = name;
	    nr = i;
	    break;
	  }
    }
    if (nr > -1) return nr;
    static bool first_overflow = true;
    if (first_overflow)
      {
        first_overflow = false;
        cerr << "no more timer available, reusing last one" << endl;
      }
    return 0;
    // throw Exception ("no more timer available");
  }

  void NgProfiler :: Reset () 
  {
      for(int i=0; i<SIZE; i++) {
          tottimes[i] = 0;
          counts[i] = 0;
          flops[i] = 0;
          loads[i] = 0;
          stores[i] = 0;
      }
  }

  NgProfiler prof;



#ifdef  VTRACE
#ifdef PARALLEL
  Timer * Timer::stack_top = NULL;
#endif
#endif


}
