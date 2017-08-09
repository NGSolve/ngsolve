/********************************************************************/
/* File:   taskmanager.cpp                                          */
/* Author: M. Hochsterger, J. Schoeberl                             */
/* Date:   10. Mar. 2015                                            */
/********************************************************************/

#include <ngstd.hpp>
#include <thread>

#include "taskmanager.hpp"
#include "paje_interface.hpp"

#ifdef USE_MKL
#include <mkl.h>
#endif

namespace ngstd
{
  TaskManager * task_manager = nullptr;
  bool TaskManager :: use_paje_trace = false;
  int TaskManager :: max_threads = getenv("NGS_NUM_THREADS") ? atoi(getenv("NGS_NUM_THREADS")) : std::thread::hardware_concurrency();
  int TaskManager :: num_threads = 1;
  // #ifndef __clang__      
  thread_local int TaskManager :: thread_id = 0;
  // #else
  // __thread int TaskManager :: thread_id;
  // #endif

  const function<void(TaskInfo&)> * TaskManager::func;
  const function<void()> * TaskManager::startup_function = nullptr;
  const function<void()> * TaskManager::cleanup_function = nullptr;

  
  static mutex copyex_mutex;

  int EnterTaskManager ()
  {
    if (task_manager)
      {
        // no task manager started
        return 0;
      }

    task_manager = new TaskManager();

    cout << IM(3) << "task-based parallelization (C++11 threads) using "<< task_manager->GetNumThreads() << " threads" << endl;

#ifdef USE_NUMA
    numa_run_on_node (0);
#endif

#ifndef WIN32
    // master has maximal priority !
    int policy;
    struct sched_param param;
    pthread_getschedparam(pthread_self(), &policy, &param);
    param.sched_priority = sched_get_priority_max(policy);
    pthread_setschedparam(pthread_self(), policy, &param);
#endif // WIN32

    
    task_manager->StartWorkers();

    ParallelFor (Range(100), [&] (int i) { ; });    // startup
    return task_manager->GetNumThreads();
  }


  void ExitTaskManager (int num_threads)
  {
    if(num_threads > 0)
      {
        task_manager->StopWorkers();
        delete task_manager;
        task_manager = nullptr;
      }
  }

  void RunWithTaskManager (function<void()> alg)
  {
    int num_threads = EnterTaskManager();
    alg();
    ExitTaskManager(num_threads);
  }




  void TaskManager :: SetNumThreads(int amax_threads)
    { 
      if(task_manager && task_manager->active_workers>0)
        {
          cerr << "Warning: can't change number of threads while TaskManager active!" << endl;
          return;
        }
      max_threads = amax_threads;
    }


  TaskManager :: TaskManager()
    {
      num_threads = GetMaxThreads();
  
      if (MyMPI_GetNTasks() > 1) num_threads = 1;

#ifdef USE_NUMA
      numa_available();
      num_nodes = numa_max_node() + 1;
      if (num_nodes > num_threads) num_nodes = num_threads;

      for (int j = 0; j < num_nodes; j++)
        {
          void * mem = numa_alloc_onnode (sizeof(NodeData), j);
          nodedata[j] = new (mem) NodeData;
	  complete[j] = -1;
          workers_on_node[j] = 0;          
          // nodedata[j]->participate = -1;
        }
#else
      num_nodes = 1;
      nodedata[0] = new NodeData;
      complete[0] = -1;
      workers_on_node[0] = 0;
      // nodedata[0]->participate = -1;
#endif

      jobnr = 0;
      done = 0;
      sleep = false;
      sleep_usecs = 1000;
      active_workers = 0;

      static int cnt = 0;
      char buf[100];
      if (use_paje_trace)
	sprintf(buf, "ng%d.trace", cnt++);
      else
        buf[0] = 0;
      //sprintf(buf, "");
      trace = new PajeTrace(num_threads, buf);
    }


  TaskManager :: ~TaskManager ()
  {
    delete trace;
    trace = nullptr;
  }


  void TaskManager :: StartWorkers()
  {
    done = false;
    sync.SetSize(num_threads);
    sync[0] = new atomic<int>(0);
    completed_tasks = 0;
    nodedata[0]->completed_tasks = 0;
    for (int i = 1; i < num_threads; i++)
      {
        std::thread([this,i]() { this->Loop(i); }).detach();
      }

    size_t alloc_size = num_threads*NgProfiler::SIZE;
    NgProfiler::thread_times = new size_t[alloc_size];
    for (size_t i = 0; i < alloc_size; i++)
      NgProfiler::thread_times[i] = 0;
    NgProfiler::thread_flops = new size_t[alloc_size];
    for (size_t i = 0; i < alloc_size; i++)
      NgProfiler::thread_flops[i] = 0;

    while (active_workers < num_threads-1)
      ;
  }
  extern size_t dummy_thread_times[NgProfiler::SIZE];
  extern size_t dummy_thread_flops[NgProfiler::SIZE];
  void TaskManager :: StopWorkers()
  {
    done = true;

    // collect timings
    for (size_t i = 0; i < num_threads; i++)
      for (size_t j = 0; j < NgProfiler::SIZE; j++)
        {
          NgProfiler::tottimes[j] += 1.0/3.1e9 * NgProfiler::thread_times[i*NgProfiler::SIZE+j];
          NgProfiler::flops[j] += NgProfiler::thread_flops[i*NgProfiler::SIZE+j];
        }
    delete [] NgProfiler::thread_times;
    NgProfiler::thread_times = dummy_thread_times;
    delete [] NgProfiler::thread_flops;
    NgProfiler::thread_flops = dummy_thread_flops;    
    while (active_workers)
      ;
    delete sync[0];
    // cout << "workers all stopped !!!!!!!!!!!!!!!!!!!" << endl;
  }



  void TaskManager :: CreateJob (const function<void(TaskInfo&)> & afunc,
                                 int antasks)
  {
    trace->StartJob(jobnr, afunc.target_type());
    /*
    for (int j = 0; j < num_nodes; j++)
      {
        while (nodedata[j]->participate > 0);
        int oldval = 0;
        while (!nodedata[j]->participate.compare_exchange_weak (oldval, -1))
          oldval = 0;
      }
    */
    func = &afunc;
    sync[0]->store(1); // , memory_order_release);

    ntasks.store (antasks); // , memory_order_relaxed);
    ex = nullptr;

    // atomic_thread_fence (memory_order_release);

    /*
    for (int j = 0; j < num_nodes; j++)
      {
        // nodedata[j]->start_cnt.store (0, memory_order_relaxed);
        // nodedata[j]->complete_cnt.store (0, memory_order_relaxed);

        // complete[j].store (0, memory_order_relaxed);
      }
    */

    nodedata[0]->start_cnt.store (0, memory_order_relaxed);

    // complete_cnt = 0;
    jobnr++;
    
    for (int j = 0; j < num_nodes; j++)
      {
        // nodedata[j]->participate.store (0, memory_order_release);
        // nodedata[j]->start_cnt = 0;
        nodedata[j]->participate |= 1;
        // nodedata[j]->participate.store (1, memory_order_release);
      }
    if (startup_function) (*startup_function)();
    
    int thd = 0;

    int thds = GetNumThreads();

    // int tasks_per_node = thds / num_nodes;
    int mynode = num_nodes * thd/thds;

    IntRange mytasks = Range(int(ntasks)).Split (mynode, num_nodes);
      
    // #ifdef USE_NUMA
    // numa_run_on_node (mynode);
    // #endif

    NodeData & mynode_data = *(nodedata[mynode]);

    TaskInfo ti;
    ti.nthreads = thds;
    ti.thread_nr = thd;
    ti.nnodes = num_nodes;
    ti.node_nr = mynode;

    try
      {
        while (1)
          {
            int mytask = mynode_data.start_cnt++;
            if (mytask >= mytasks.Size()) break;
            
            ti.task_nr = mytasks.First()+mytask;
            ti.ntasks = ntasks;

              {
		RegionTracer t(ti.thread_nr, jobnr, RegionTracer::ID_JOB, ti.task_nr);
                (*func)(ti); 
                mynode_data.completed_tasks++;
              }

	      // if (++mynode_data.complete_cnt == mytasks.Size())
              // complete[mynode] = true;
          }

      }
    catch (Exception e)
      {
        {
          // cout << "got exception in TM:" << ex->What() << endl;
          lock_guard<mutex> guard(copyex_mutex);
          delete ex;
          ex = new Exception (e);
          mynode_data.start_cnt = mytasks.Size();
          // complete[mynode] = true;
        }
      }

    if (cleanup_function) (*cleanup_function)();
    
    for (int j = 0; j < num_nodes; j++)
      if (workers_on_node[j])
        {
          while (complete[j] != jobnr)
            ;
          /*
          while (completed_tasks+ntasks/num_nodes != nodedata[j]->completed_tasks)
            cout << "master, check node " << j << " node complete = " << nodedata[j]->completed_tasks << " should be " << completed_tasks+ntasks/num_nodes << endl;
            ;
          */
        }


    completed_tasks += ntasks / num_nodes;
//     for (int j = 0; j < num_nodes; j++)
//       if (completed_tasks != nodedata[j]->completed_tasks)
//         cout << "tasks missing: node = " << j 
//              << "global tasks = " << completed_tasks
//              << "node tasks = " << nodedata[j]->completed_tasks << endl;

    if (ex)
      throw Exception (*ex);

    trace->StopJob();
    for (auto ap : sync)
      ap->load(); // memory_order_acquire);
  }
    
  void TaskManager :: Loop(int thd)
  {
    static Timer tADD("add entry counter");
    static Timer tCASready1("spin-CAS ready tick1");
    static Timer tCASready2("spin-CAS ready tick2");
    static Timer tCASyield("spin-CAS yield");
    static Timer tCAS1("spin-CAS wait");
    static Timer texit("exit zone");
    static Timer tdec("decrement");
    thread_id = thd;

    sync[thread_id] = new atomic<int>(0);

    int thds = GetNumThreads();

    int mynode = num_nodes * thd/thds;

    NodeData & mynode_data = *(nodedata[mynode]);



    TaskInfo ti;
    ti.nthreads = thds;
    ti.thread_nr = thd;
    ti.nnodes = num_nodes;
    ti.node_nr = mynode;

      
#ifdef USE_NUMA
    numa_run_on_node (mynode);
#endif
    active_workers++;
    workers_on_node[mynode]++;
    int jobdone = 0;


#ifdef USE_MKL
    auto mkl_max = mkl_get_max_threads();
    mkl_set_num_threads_local(1);
#endif

    
    while (!done)
      {
        if (complete[mynode] > jobdone)
          jobdone = complete[mynode];

        if (jobnr == jobdone)
          {
            // RegionTracer t(ti.thread_nr, tCASyield, ti.task_nr);            
            if(sleep)
              this_thread::sleep_for(chrono::microseconds(sleep_usecs));
            else
              {
#ifdef WIN32
                this_thread::yield();
#else  // WIN32
                sched_yield();
#endif // WIN32
              }
            continue;
          }
        
        /*
        while (mynode_data.participate.load(memory_order_relaxed) == -1)
          {
            RegionTracer t(ti.thread_nr, tCAS1, ti.task_nr);            
            if (done.load(memory_order_relaxed))
              {
                active_workers--;
                return;
              }
          }
        */


        /*
        int oldpart = 0;
        while (! mynode_data.participate.compare_exchange_weak (oldpart, oldpart+1))
	  {
            RegionTracer t(ti.thread_nr, tCAS, ti.task_nr);
	    if (oldpart == -1) oldpart = 0;
            if (done.load(memory_order_relaxed))
	      {
		active_workers--;
		return;
	      }
	  }
        */

        /*
        int oldpart = 0;
        if (! mynode_data.participate.compare_exchange_weak (oldpart, oldpart+1))
          continue;
        */

        // non-atomic fast check ...
        {
          // RegionTracer t(ti.thread_nr, tADD, ti.task_nr);
        if ( (mynode_data.participate & 1) == 0) continue;

        {
          int oldval = mynode_data.participate += 2;
          if ( (oldval & 1) == 0)
            { // job not active, going out again
              mynode_data.participate -= 2;
              continue;
            }
        }
        }

        // for (auto ap : sync)
        // ap->load(); // memory_order_acquire);
        
        for (int j = 0; j < num_nodes; j++)
          while (completed_tasks > nodedata[j]->completed_tasks)
            ;
        /*
          cout << "mynode = " << mynode << ", j = " << completed_tasks << " node comp = " << nodedata[j]->completed_tasks << endl;
        */
        

        // atomic_thread_fence (memory_order_acquire);
        if (startup_function) (*startup_function)();
        
        IntRange mytasks = Range(int(ntasks)).Split (mynode, num_nodes);
          
        try
          {
            
            while (1)
              {
                // int mytask = mynode_data.start_cnt++;
                if (mynode_data.start_cnt >= mytasks.Size()) break;
		int mytask = mynode_data.start_cnt.fetch_add(1, memory_order_relaxed);
                if (mytask >= mytasks.Size()) break;
                
                ti.task_nr = mytasks.First()+mytask;
                ti.ntasks = ntasks;
                
                  {
		    RegionTracer t(ti.thread_nr, jobnr, RegionTracer::ID_JOB, ti.task_nr);
                    (*func)(ti);
                    mynode_data.completed_tasks++;
                  }
		  // if (++mynode_data.complete_cnt == mytasks.Size())
                  // complete[mynode] = true;
              }

          }
        catch (Exception e)
          {
            {
              // cout << "got exception in TM" << endl; 
              lock_guard<mutex> guard(copyex_mutex);
              delete ex;
              ex = new Exception (e);
              mynode_data.start_cnt = mytasks.Size();
              // complete[mynode] = true;
            }
          }

#ifndef __MIC__
        atomic_thread_fence (memory_order_release);     
#endif // __MIC__

        if (cleanup_function) (*cleanup_function)();
        sync[thread_id]->store(1); // , memory_order_release);

        jobdone = jobnr;

        mynode_data.participate-=2;

	{
	  int oldpart = 1;
	  if (mynode_data.participate.compare_exchange_strong (oldpart, 0))
	    {
              if (jobdone < jobnr.load())
                { // reopen gate
                  mynode_data.participate |= 1;                  
                }
              else
                {
                  if (mynode != 0)
                    mynode_data.start_cnt = 0;
                  complete[mynode] = jobnr.load(); 
                }
	    }	      
	}
      }
    

#ifdef USE_MKL
    mkl_set_num_threads_local(mkl_max);
#endif

    delete sync[thread_id];
    workers_on_node[mynode]--;
    active_workers--;
  }


  list<tuple<string,double>> TaskManager :: Timing ()
  {
    /*
    list<tuple<string,double>>timings;
    double time =
      RunTiming
      ( [&] ()
        {
          ParallelJob ( [] (TaskInfo ti) { ; } ,
                        TasksPerThread(1) );
        });
    timings.push_back (make_tuple("parallel job with 1 task per thread", time*1e9));
    
    time =
      RunTiming
      ( [&] ()
        {
          ParallelJob ( [] (TaskInfo ti) { ; } ,
                        TasksPerThread(10) );
        });
    timings.push_back (make_tuple("parallel job with 10 tasks per thread", time*1e9));

    time =
      RunTiming
      ( [&] ()
        {
          ParallelJob ( [] (TaskInfo ti) { ; } ,
                        TasksPerThread(100) );
        });
    timings.push_back (make_tuple("parallel job with 100 tasks per thread", time*1e9));

    return timings;
    */


    
    // this is the old function moved from the py-interface:
    list<tuple<string,double>>timings;           
    double starttime, time;
    double maxtime = 0.5;
    size_t steps;
    
    starttime = WallTime();
    steps = 0;
    do
      {
        for (size_t i = 0; i < 1000; i++)
          ParallelJob ( [] (TaskInfo ti) { ; },
                        TasksPerThread(1));
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("ParallelJob 1 task/thread", time/steps*1e9));


    starttime = WallTime();
    steps = 0;
    do
      {
        for (size_t i = 0; i < 1000; i++)
          ParallelJob ( [] (TaskInfo ti) { ; },
                        TasksPerThread(100));
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("ParallelJob 100 task/thread", time/steps*1e9));

    
    starttime = WallTime();
    steps = 0;
    do
      {
        for (int k = 0; k < 10000; k++)
          {
            SharedLoop2 sl(1000);
            steps += 1;
          }
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("SharedLoop init", time/steps*1e9));
    
    starttime = WallTime();
    steps = 0;
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            SharedLoop2 sl(1000);
            ParallelJob ( [&sl] (TaskInfo ti)
                          {
                            // auto beg = sl.begin();
                            // auto end = sl.end();
                          } );
            steps += 1;
          }
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("SharedLoop begin/end", time/steps*1e9));
    
    
    starttime = WallTime();
    steps = 0;
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            SharedLoop2 sl(1000);
            ParallelJob ( [&sl] (TaskInfo ti)
                          {
                            for (auto i : sl)
                              ; 
                          } );
            steps += 1000;
          }
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("SharedLoop 1000", time/steps*1e9));
    
    starttime = WallTime();
    steps = 0;
    do
      {
        SharedLoop2 sl(1000000);
        ParallelJob ( [&sl] (TaskInfo ti)
                      {
                        for (auto i : sl)
                          ; 
                      } );
        steps += 1000000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("SharedLoop 1000000", time/steps*1e9));
    
    return timings;
  }
  
}
