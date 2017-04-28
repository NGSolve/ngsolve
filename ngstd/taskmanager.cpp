/********************************************************************/
/* File:   taskmanager.cpp                                          */
/* Author: M. Hochsterger, J. Schoeberl                             */
/* Date:   10. Mar. 2015                                            */
/********************************************************************/

#include <ngstd.hpp>
#include <thread>

#include "taskmanager.hpp"
#include "paje_interface.hpp"


namespace ngstd
{
  TaskManager * task_manager = nullptr;
  bool TaskManager :: use_paje_trace = false;
  int TaskManager :: max_threads = getenv("NGS_NUM_THREADS") ? atoi(getenv("NGS_NUM_THREADS")) : std::thread::hardware_concurrency();
  int TaskManager :: num_threads = 1;
#ifndef __clang__      
  thread_local int TaskManager :: thread_id;
#else
  __thread int TaskManager :: thread_id;
#endif

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
        }
#else
      num_nodes = 1;
      nodedata[0] = new NodeData;
      nodedata[0]->participate = 0;
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

    for (int i = 1; i < num_threads; i++)
      {
        std::thread([this,i]() { this->Loop(i); }).detach();
      }

    while (active_workers < num_threads-1)
      ;
  }
 
  void TaskManager :: StopWorkers()
  {
    done = true;
    while (active_workers)
      ;
    // cout << "workers all stopped !!!!!!!!!!!!!!!!!!!" << endl;
  }



  void TaskManager :: CreateJob (const function<void(TaskInfo&)> & afunc,
                                 int antasks)
  {

    trace->StartJob(jobnr, afunc.target_type());
    for (int j = 0; j < num_nodes; j++)
      {
        while (nodedata[j]->participate > 0);
        int oldval = 0;
        while (!nodedata[j]->participate.compare_exchange_weak (oldval, -1))
          oldval = 0;
      }

    func = &afunc;


    ntasks.store (antasks, memory_order_relaxed);
    ex = nullptr;

    // atomic_thread_fence (memory_order_release);

    for (int j = 0; j < num_nodes; j++)
      {
        nodedata[j]->start_cnt.store (0, memory_order_relaxed);
        nodedata[j]->complete_cnt.store (0, memory_order_relaxed);

        complete[j].store (0, memory_order_relaxed);
      }

    // complete_cnt = 0;
    jobnr++;
      
    for (int j = 0; j < num_nodes; j++)
      nodedata[j]->participate.store (0, memory_order_release);

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
              }

            if (++mynode_data.complete_cnt == mytasks.Size())
              complete[mynode] = true;
          }

      }
    catch (Exception e)
      {
        {
          lock_guard<mutex> guard(copyex_mutex);
          delete ex;
          ex = new Exception (e);
          mynode_data.start_cnt = mytasks.Size();
          complete[mynode] = true;
        }
      }

    if (cleanup_function) (*cleanup_function)();
    
    for (int j = 0; j < num_nodes; j++)
      while (!complete[j])
        ;

    if (ex)
      throw Exception (*ex);

    trace->StopJob();
    // atomic_thread_fence (memory_order_acquire);
  }
    
  void TaskManager :: Loop(int thd)
  {
    thread_id = thd;
    active_workers++;

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
      
    int jobdone = 0;
      
    while (!done)
      {
        if (jobnr == jobdone)
          {
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
          
        while (mynode_data.participate == -1);
          
        int oldpart = 0;
        while (! mynode_data.participate.compare_exchange_weak (oldpart, oldpart+1))
          if (oldpart == -1) oldpart = 0;

        // atomic_thread_fence (memory_order_acquire);
        if (startup_function) (*startup_function)();
        
        IntRange mytasks = Range(int(ntasks)).Split (mynode, num_nodes);
          
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
                if (++mynode_data.complete_cnt == mytasks.Size())
                  complete[mynode] = true;
              }

          }
        catch (Exception e)
          {
            {
              lock_guard<mutex> guard(copyex_mutex);
              delete ex;
              ex = new Exception (e);
              mynode_data.start_cnt = mytasks.Size();
              complete[mynode] = true;
            }
          }

#ifndef __MIC__
        atomic_thread_fence (memory_order_release);     
#endif // __MIC__

        if (cleanup_function) (*cleanup_function)();
        
        jobdone = jobnr;
        mynode_data.participate--;
      }

    active_workers--;
  }

}
