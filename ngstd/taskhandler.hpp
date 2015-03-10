#ifndef FILE_TASKHANDLER
#define FILE_TASKHANDLER

/*********************************************************************/
/* File:   taskhandler.hpp                                           */
/* Author: M. Hochsterger, J. Schoeberl                              */
/* Date:   10. Mar. 2015                                             */
/*********************************************************************/


#ifdef USE_NUMA
#include <numa.h>
#endif


namespace ngstd
{
  
  class TaskHandler
  {
    
    function<void(int)> func;
    volatile int jobnr = 0;
    volatile int cnt;
    volatile int done = false;
    int num_nodes;
    
  public:
    
    TaskHandler()
    {
#ifdef USE_NUMA
      numa_available();
      num_nodes = numa_max_node() + 1;
#else
      num_nodes = 1;
#endif
    }


    void CreateTask (function<void(int)> afunc)
    {
      func = afunc;
      
      cnt = 0;
      jobnr++;
      
      int thd = omp_get_thread_num();
      int thds = omp_get_num_threads();

#ifdef USE_NUMA
      numa_run_on_node (num_nodes*thd/thds);
#endif
      
      func(thd);

#pragma omp atomic
      cnt++;
      
      while (cnt < thds)
	; 
    }

    void Done()
    {
      done = true;
    }

    void Loop()
    {
      int thd = omp_get_thread_num();
      int thds = omp_get_num_threads();

#ifdef USE_NUMA
      numa_run_on_node (num_nodes*thd/thds);
#endif

      int jobdone = 0;
      
      while (!done)
	{
	  // usleep (1);
	  // _mm_pause();
	  if (jobnr > jobdone)
	    {
	      func(thd);
	      jobdone++;
#pragma omp atomic
	      cnt++;
	    }
	}
    }
    
  };
  
}



#endif
