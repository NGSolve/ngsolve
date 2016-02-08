/********************************************************************/
/* File:   taskmanager.cpp                                          */
/* Author: M. Hochsterger, J. Schoeberl                             */
/* Date:   10. Mar. 2015                                            */
/********************************************************************/

#include <ngstd.hpp>
#include <thread>

#include "taskmanager.hpp"
#include "paje_interface.hpp"

/*
#ifdef USE_NUMA
#include <numa.h>
#include <sched.h>
#endif
*/


// #define CPP11_THREADS 

namespace ngstd
{
  TaskManager * task_manager = nullptr;
  bool TaskManager :: use_paje_trace = false;

  int EnterTaskManager ()
  {
    if (task_manager)
      {
        // no task manager started
        return 0;
      }

    task_manager = new TaskManager();

    cout << IM(3) << "new task-based parallelization (C++11 threads)" << endl;

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


    task_manager->StartWorkersWithCppThreads();

    ParallelFor (100, [&] (int i) { ; });    // startup
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

#ifdef CPP11_THREADS
  void RunWithTaskManager (function<void()> alg)
  {
    int num_threads = EnterTaskManager();
    alg();
    ExitTaskManager(num_threads);
  }

#else // openmp-threads
  void RunWithTaskManager (function<void()> alg)
  {
    if (task_manager)
      {
        alg();
        return;
      }

    task_manager = new TaskManager();

    
    Exception * ex = nullptr;

#pragma omp parallel
    {
#pragma omp single 
      {
        try
          {
            cout << IM(3) << "new task-based parallelization" << endl;
            
#ifdef USE_NUMA
            int num_nodes = task_manager->GetNumNodes();
            
            int thd = omp_get_thread_num();
            int thds = omp_get_num_threads();
            int mynode = num_nodes * thd/thds;
            numa_run_on_node (mynode);
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
            ParallelFor (100, [&] (int i) { ; });    // startup

            alg();
            
            task_manager->StopWorkers();
          }
        catch (Exception e)
          {
            ex = new Exception (e);
          }
      }
    }

    if (ex) throw Exception (*ex);


    
    delete task_manager;
    task_manager = nullptr;
  }

#endif





  TaskManager :: TaskManager()
    {
      num_threads = omp_get_max_threads();
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
      active_workers = 0;

      static int cnt = 0;
      char buf[100];
      if (use_paje_trace)
	sprintf(buf, "ng%d.trace", cnt++);
      else
	sprintf(buf, "");
      trace = new PajeTrace(num_threads, buf);
    }


  TaskManager :: ~TaskManager ()
  {
    delete trace;
    trace = nullptr;
  }

  void TaskManager :: StartWorkersWithCppThreads()
  {
    done = false;

    for (int i = 1; i < num_threads; i++)
      {
        std::thread([this,i]() { this->Loop(i); }).detach();
      }

    while (active_workers < num_threads-1)
      ;
  }


  void TaskManager :: StartWorkers()
  {
#ifdef CPP11_THREADS
    StartWorkersWithCppThreads();
#else
    done = false;

    for (int i = 0; i < num_threads-1; i++)
#pragma omp task
      {
        Loop(omp_get_thread_num());
      }

    while (active_workers < num_threads-1)
      ;
    // cout << "workers are all active !!!!!!!!!!!" << endl;
#endif
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
    

#ifdef CPP11_THREADS
    int thd = 0;
#else
    int thd = omp_get_thread_num();
#endif

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
#pragma omp critical(copyex)
        {
          delete ex;
          ex = new Exception (e);
          mynode_data.start_cnt = mytasks.Size();
          complete[mynode] = true;
        }
      }
    
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
#ifdef WIN32
            this_thread::yield();
#else  // WIN32
            sched_yield();
#endif // WIN32
            continue;
          }
          
        while (mynode_data.participate == -1);
          
        int oldpart = 0;
        while (! mynode_data.participate.compare_exchange_weak (oldpart, oldpart+1))
          if (oldpart == -1) oldpart = 0;

        // atomic_thread_fence (memory_order_acquire);

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
#pragma omp critical(copyex)
            {
              delete ex;
              ex = new Exception (e);
              mynode_data.start_cnt = mytasks.Size();
              complete[mynode] = true;
            }
          }

#ifndef __MIC__
        atomic_thread_fence (memory_order_release);     
#endif // __MIC__

        jobdone = jobnr;
        mynode_data.participate--;
      }

    active_workers--;
  }

}
