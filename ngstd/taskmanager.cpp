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


#include "../linalg/concurrentqueue.h"



namespace ngstd
{
  TaskManager * task_manager = nullptr;
  bool TaskManager :: use_paje_trace = false;
  int TaskManager :: max_threads = getenv("NGS_NUM_THREADS") ? atoi(getenv("NGS_NUM_THREADS")) : std::thread::hardware_concurrency();
  int TaskManager :: num_threads = 1;

  
#ifndef __clang__      
  thread_local int TaskManager :: thread_id = 0;
#else
  __thread int TaskManager :: thread_id;
#endif
  
  thread_local int thread_id = 0;

  const function<void(TaskInfo&)> * TaskManager::func;
  const function<void()> * TaskManager::startup_function = nullptr;
  const function<void()> * TaskManager::cleanup_function = nullptr;

  atomic<int> TaskManager::ntasks;
  Exception * TaskManager::ex;
  
  atomic<int> TaskManager::jobnr;
  
  atomic<int> TaskManager::complete[8];   // max nodes
  atomic<int> TaskManager::done;
  atomic<int> TaskManager::active_workers;
  atomic<int> TaskManager::workers_on_node[8];   // max nodes

  
  int TaskManager::sleep_usecs = 1000;
  bool TaskManager::sleep = false;

  TaskManager::NodeData *TaskManager::nodedata[8];
  int TaskManager::num_nodes;
  
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
      // if (MyMPI_GetNTasks() > 1) num_threads = 1;

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
        }
#else
      num_nodes = 1;
      nodedata[0] = new NodeData;
      complete[0] = -1;
      workers_on_node[0] = 0;
#endif

      jobnr = 0;
      done = 0;
      sleep = false;
      sleep_usecs = 1000;
      active_workers = 0;

      static int cnt = 0;
      char buf[100];
      if (use_paje_trace)
        {
#ifdef PARALLEL
          sprintf(buf, "ng%d_rank%d.trace", cnt++, MyMPI_GetId());
#else
          sprintf(buf, "ng%d.trace", cnt++);
#endif
        }
      else
        buf[0] = 0;
      //sprintf(buf, "");
      trace = new PajeTrace(num_threads, buf);
    }


  TaskManager :: ~TaskManager ()
  {
    delete trace;
    trace = nullptr;
    num_threads = 1;
  }

  /*
  int TaskManager :: GetThreadId()
  {
    return thread_id;
  }
  */
  
  void TaskManager :: StartWorkers()
  {
    done = false;

    for (int i = 1; i < num_threads; i++)
      {
        std::thread([this,i]() { this->Loop(i); }).detach();
      }
    thread_id = 0;
    
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

  static size_t calibrate_init_tsc = __rdtsc();
  typedef std::chrono::system_clock TClock;
  static TClock::time_point calibrate_init_clock = TClock::now();
  
  void TaskManager :: StopWorkers()
  {
    done = true;
    double delta_tsc = __rdtsc()-calibrate_init_tsc;
    double delta_sec = std::chrono::duration<double>(TClock::now()-calibrate_init_clock).count();
    double frequ = (delta_sec != 0) ? delta_tsc/delta_sec : 2.7e9;
    
    // cout << "cpu frequ = " << frequ << endl;
    // collect timings
    for (size_t i = 0; i < num_threads; i++)
      for (size_t j = NgProfiler::SIZE; j-- > 0; )
        {
          if (!NgProfiler::usedcounter[j]) break;
          NgProfiler::tottimes[j] += 1.0/frequ * NgProfiler::thread_times[i*NgProfiler::SIZE+j];
          NgProfiler::flops[j] += NgProfiler::thread_flops[i*NgProfiler::SIZE+j];
        }
    delete [] NgProfiler::thread_times;
    NgProfiler::thread_times = dummy_thread_times;
    delete [] NgProfiler::thread_flops;
    NgProfiler::thread_flops = dummy_thread_flops;
    
    while (active_workers)
      ;
  }

  /////////////////////// NEW: nested tasks using concurrent queue

  struct TNestedTask
  {
    const function<void(TaskInfo&)> * func;
    atomic<int> * endcnt;
    int mynr;
    int total;

    TNestedTask () { ; }
    TNestedTask (const function<void(TaskInfo&)> & _func,
                 int _mynr, int _total,
                 atomic<int> & _endcnt)
      : func(&_func), mynr(_mynr), total(_total), endcnt(&_endcnt)
    {
      ;
    }
  };

  typedef moodycamel::ConcurrentQueue<TNestedTask> TQueue; 
  typedef moodycamel::ProducerToken TPToken; 
  typedef moodycamel::ConsumerToken TCToken; 
  
  static TQueue taskqueue;

  void AddTask (const function<void(TaskInfo&)> & afunc,
                atomic<int> & endcnt)
                
  {
    TPToken ptoken(taskqueue); 

    int num = endcnt;
    for (int i = 0; i < num; i++)
      taskqueue.enqueue (ptoken, { afunc, i, num, endcnt });
  }

  mutex m;
  bool ProcessTask()
  {
    TNestedTask task;
    TCToken ctoken(taskqueue); 
    
    if (taskqueue.try_dequeue(ctoken, task))
      {
        TaskInfo ti;
        ti.task_nr = task.mynr;
        ti.ntasks = task.total;
        ti.thread_nr = thread_id;
        ti.nthreads = TaskManager::GetNumThreads();
        /*
        {
          lock_guard<mutex> guard(m);
          cout << "process nested, nr = " << ti.task_nr << "/" << ti.ntasks << endl;
        }
        */
        (*task.func)(ti);
        --*task.endcnt;
        return true;
      }
    return false;
  }


  void TaskManager :: CreateJob (const function<void(TaskInfo&)> & afunc,
                                 int antasks)
  {
    if (num_threads == 1 || !task_manager) //  || func)
      {
        if (startup_function) (*startup_function)();
        
        TaskInfo ti;
        ti.ntasks = antasks;
        ti.thread_nr = 0; ti.nthreads = 1;
        // ti.node_nr = 0; ti.nnodes = 1;
        for (ti.task_nr = 0; ti.task_nr < antasks; ti.task_nr++)
          afunc(ti);

        if (cleanup_function) (*cleanup_function)();        
        return;
      }


    if (func)
      { // we are already parallel, use nested tasks
        // startup for inner function not supported ...
        // if (startup_function) (*startup_function)();

        if (antasks == 1)
          {
            TaskInfo ti;
            ti.task_nr = 0;
            ti.ntasks = 1;
            ti.thread_nr = 0; ti.nthreads = 1;
            afunc(ti);
            return;
          }
        
        atomic<int> endcnt(antasks);
        AddTask (afunc, endcnt);
        while (endcnt > 0)
          {
            ProcessTask();
          }
        
        // if (cleanup_function) (*cleanup_function)();
        return;
      }
    
    
    trace->StartJob(jobnr, afunc.target_type());

    func = &afunc;

    ntasks.store (antasks); // , memory_order_relaxed);
    ex = nullptr;


    nodedata[0]->start_cnt.store (0, memory_order_relaxed);

    jobnr++;
    
    for (int j = 0; j < num_nodes; j++)
      nodedata[j]->participate |= 1;

    if (startup_function) (*startup_function)();
    
    int thd = 0;
    int thds = GetNumThreads();
    int mynode = num_nodes * thd/thds;

    IntRange mytasks = Range(int(ntasks)).Split (mynode, num_nodes);
    NodeData & mynode_data = *(nodedata[mynode]);

    TaskInfo ti;
    ti.nthreads = thds;
    ti.thread_nr = thd;
    // ti.nnodes = num_nodes;
    // ti.node_nr = mynode;

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
            }
          }

      }
    catch (Exception e)
      {
        {
          lock_guard<mutex> guard(copyex_mutex);
          delete ex;
          ex = new Exception (e);
          mynode_data.start_cnt = mytasks.Size();
        }
      }

    if (cleanup_function) (*cleanup_function)();
    
    for (int j = 0; j < num_nodes; j++)
      if (workers_on_node[j])
        {
          while (complete[j] != jobnr)
            _mm_pause();
        }

    func = nullptr;
    if (ex)
      throw Exception (*ex);

    trace->StopJob();
  }
    
  void TaskManager :: Loop(int thd)
  {
    /*
    static Timer tADD("add entry counter");
    static Timer tCASready1("spin-CAS ready tick1");
    static Timer tCASready2("spin-CAS ready tick2");
    static Timer tCASyield("spin-CAS yield");
    static Timer tCAS1("spin-CAS wait");
    static Timer texit("exit zone");
    static Timer tdec("decrement");
    */
    thread_id = thd;

    int thds = GetNumThreads();

    int mynode = num_nodes * thd/thds;

    NodeData & mynode_data = *(nodedata[mynode]);



    TaskInfo ti;
    ti.nthreads = thds;
    ti.thread_nr = thd;
    // ti.nnodes = num_nodes;
    // ti.node_nr = mynode;

      
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
            while (ProcessTask()); // do the nested tasks
                   
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

        {
          // RegionTracer t(ti.thread_nr, tADD, ti.task_nr);

          // non-atomic fast check ...
          if ( (mynode_data.participate & 1) == 0) continue;

          int oldval = mynode_data.participate += 2;
          if ( (oldval & 1) == 0)
            { // job not active, going out again
              mynode_data.participate -= 2;
              continue;
            }
        }

        if (startup_function) (*startup_function)();
        
        IntRange mytasks = Range(int(ntasks)).Split (mynode, num_nodes);
          
        try
          {
            
            while (1)
              {
                if (mynode_data.start_cnt >= mytasks.Size()) break;
		int mytask = mynode_data.start_cnt.fetch_add(1, memory_order_relaxed);
                if (mytask >= mytasks.Size()) break;
                
                ti.task_nr = mytasks.First()+mytask;
                ti.ntasks = ntasks;
                
                {
                  RegionTracer t(ti.thread_nr, jobnr, RegionTracer::ID_JOB, ti.task_nr);
                  (*func)(ti);
                }
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
            }
          }

#ifndef __MIC__
        atomic_thread_fence (memory_order_release);     
#endif // __MIC__

        if (cleanup_function) (*cleanup_function)();

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
            SharedLoop sl(5);
            ParallelJob ( [&sl] (TaskInfo ti)
                          {
                            for (auto i : sl)
                              (void)i;  // silence warning
                          } );
          }
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("short SharedLoop", time/steps*1e9));
    

    starttime = WallTime();
    steps = 0;
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            SharedLoop sl1(5), sl2(5), sl3(5), sl4(5), sl5(5);
            ParallelJob ( [&sl1, &sl2, &sl3, &sl4, &sl5] (TaskInfo ti)
                          {
                            for (auto i : sl1)
                              (void)i;  // silence warning
                            for (auto i : sl2)
                              (void)i;  // silence warning
                            for (auto i : sl3)
                              (void)i;  // silence warning
                            for (auto i : sl4)
                              (void)i;  // silence warning
                            for (auto i : sl5)
                              (void)i;  // silence warning
                          } );
          }
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("5 short SharedLoops", time/steps*1e9));
    

    starttime = WallTime();
    steps = 0;
    SharedLoop2 sl2(5);
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            sl2.Reset(5);
            ParallelJob ( [&sl2] (TaskInfo ti)
                          {
                            for (auto i : sl2)
                              (void)i;  // silence warning                              
                          } );
          }
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("short SharedLoop2", time/steps*1e9));

    {
    starttime = WallTime();
    steps = 0;
    SharedLoop2 sl1(5), sl2(5), sl3(5), sl4(5), sl5(5);
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            sl1.Reset(5);
            sl2.Reset(5);
            sl3.Reset(5);
            sl4.Reset(5);
            sl5.Reset(5);
            ParallelJob ( [&sl1,&sl2,&sl3,&sl4,&sl5] (TaskInfo ti)
                          {
                            for (auto i : sl1)
                              (void)i;  // silence warning                              
                            for (auto i : sl2)
                              (void)i;  // silence warning                              
                            for (auto i : sl3)
                              (void)i;  // silence warning                              
                            for (auto i : sl4)
                              (void)i;  // silence warning                              
                            for (auto i : sl5)
                              (void)i;  // silence warning                              
                          } );
          }
        steps += 1000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("5 short SharedLoop2", time/steps*1e9));
    }

    
    starttime = WallTime();
    steps = 0;
    {
    SharedLoop2 sl(1000);
    do
      {
        for (int k = 0; k < 1000; k++)
          {
            sl.Reset(1000);
            ParallelJob ( [&sl] (TaskInfo ti)
                          {
                            for (auto i : sl)
                              (void)i;  // silence warning                               
                          } );
            steps += 1000;
          }
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("SharedLoop2 1000, time per iteration", time/steps*1e9));
    }

    {
    starttime = WallTime();
    steps = 0;
    SharedLoop2 sl(1000000);
    do
      {
        sl.Reset(1000000);
        ParallelJob ( [&sl] (TaskInfo ti)
                      {
                        for (auto i : sl)
                          (void)i;  // silence warning
                      } );
        steps += 1000000;
        time = WallTime()-starttime;
      }
    while (time < maxtime);
    timings.push_back(make_tuple("SharedLoop2 1000000, time per iteration", time/steps*1e9));
    }
    
    return timings;
  }
  
}
