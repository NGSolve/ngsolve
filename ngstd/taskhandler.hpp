#ifndef FILE_TASKHANDLER
#define FILE_TASKHANDLER

/*********************************************************************/
/* File:   taskhandler.hpp                                           */
/* Author: M. Hochsterger, J. Schoeberl                              */
/* Date:   10. Mar. 2015                                             */
/*********************************************************************/


#ifdef USE_NUMA
#include <numa.h>
#include <sched.h>
#endif


namespace ngstd
{
  
  class TaskHandler
  {
    
    function<void(int)> func;
    volatile int jobnr = 0;

    volatile atomic<int> start_cnt[8]; // max nodes
    volatile atomic<int> complete_cnt;
    volatile int done = false;
    volatile atomic<int> inside;
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

      inside = 0;
      
    }


    void CreateTask (function<void(int)> afunc)
    {
      func = afunc;
      
      for (int j = 0; j < num_nodes; j++)
	start_cnt[j] = 0;
      complete_cnt = 0;

      jobnr++;
      
      int thd = omp_get_thread_num();
      int thds = omp_get_num_threads();

      int tasks_per_node = thds / num_nodes;
      int mynode = num_nodes * thd/thds;

#ifdef USE_NUMA
      numa_run_on_node (mynode);
#endif

      int mytask;
      // #pragma omp atomic capture
      mytask = start_cnt[mynode]++;
      
      if (mytask < tasks_per_node)
	{
	  func(mytask + mynode*tasks_per_node);
	  // #pragma omp atomic
	  complete_cnt++;
	}
      
      // cout << "A" << flush;
      while (complete_cnt < thds)
	{
	  // cout << "comp cnt = " << complete_cnt << endl;
	  // cout << "start cnt = " << start_cnt[0] << " " << start_cnt[1] << endl;
	}
      // cout << "after: start cnt = " << start_cnt[0] << " " << start_cnt[1] << endl;
      // cout << "B" << flush;
      while (inside) ; 
      // cout << "C" << flush;

    }


    void Done()
    {
      done = true;
    }

    void Loop()
    {
      int thd = omp_get_thread_num();
      int thds = omp_get_num_threads();

      int tasks_per_node = thds / num_nodes;
      int mynode = num_nodes * thd/thds;

#ifdef USE_NUMA
      numa_run_on_node (mynode);
#endif

      int jobdone = 0;
      
      while (!done)
	{
	  
	  int curjob;
	  curjob = jobnr;

	  if (curjob > jobdone)
	    {
	      inside++;

	      if (curjob == jobnr) // still the same job
		
		while (1)
		  {
		    int mytask;
		  // #pragma omp atomic capture
		      mytask = start_cnt[mynode]++;
		      
		      if (mytask >= tasks_per_node) break;
		      
		      func(mytask + mynode*tasks_per_node);
		      
		      // #pragma omp atomic
			complete_cnt++;
		  }
	      
	      jobdone = jobnr;

	      inside--;
	    }

	  // sched_yield();
	}
    }
    
  };
  
}



#endif
