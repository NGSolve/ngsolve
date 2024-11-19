/*********************************************************************/
/* File:   cuda_linalg.cpp                                           */
/* Author: Joachim Schoeberl, Matthias Hochsteger                    */
/* Date:   11. Aug. 2014                                             */
/*********************************************************************/

#include <la.hpp>
#include "cuda_linalg.hpp"
#include <nvfunctional>

using namespace std;



#include "concurrentqueue.h"

class Task
{
public:
  int nr, size;
  const nvstd::function<void(int nr, int size)> * pfunc;
  // std::atomic<int> * cnt;
  int * cnt;
  
  Task & operator++(int)
  {
    nr++;
    return *this;
  }
  Task & operator*() { return *this; }
};

  
  typedef moodycamel::ConcurrentQueue<Task> TQueue; 
  typedef moodycamel::ProducerToken TPToken; 
  typedef moodycamel::ConsumerToken TCToken; 
  
  
  // __device__ static std::atomic<bool> stop{false};
  __device__ static bool stop{false};
  // __device__ static std::vector<std::thread> threads;
  __device__ static TQueue * queue;
  
  __device__ void Worker()
  {
    stop = false;
          
    TPToken ptoken(*queue); 
    TCToken ctoken(*queue); 
    
    while(true)
      {
        if (stop) break;
        
        Task task;
        if(!queue->try_dequeue_from_producer(ptoken, task)) 
          if(!queue->try_dequeue(ctoken, task))  
            continue; 
        
        (*task.pfunc)(task.nr, task.size);
        // (*task.cnt)++;
        atomicAdd(task.cnt, 1);
      }
  }

__global__ void StopWorkers()
{
  stop = true;
  }

  
__global__ void RunParallel (int num,
                             const nvstd::function<void(int nr, int size)> & func)
  {
    TPToken ptoken(*queue);
    TCToken ctoken(*queue);
    
    // std::atomic<int> cnt{0};
    int cnt = 0;

    for (size_t i = 0; i < num; i++)
      {
        Task task;
        task.nr = i;
        task.size = num;
        task.pfunc = &func;
        task.cnt = &cnt;
        queue->enqueue(ptoken, task);
      }

    /*
    // faster with bulk enqueue (error with gcc-Release)
    Task firsttask;
    firsttask.nr = 0;
    firsttask.size = num;
    firsttask.pfunc=&func;
    firsttask.cnt = &cnt;
    queue.enqueue_bulk (ptoken, firsttask, num);    
    */
    
    while (cnt < num)
      {
        Task task;
        if(!queue->try_dequeue_from_producer(ptoken, task)) 
          if(!queue->try_dequeue(ctoken, task))
            continue; 
        
        (*task.pfunc)(task.nr, task.size);
        // (*task.cnt)++;
        atomicAdd(task.cnt, 1);

      }
  }


// #ifdef __CUDACC__

__global__ void RunAll() // int n, F f)
{
  int bid = blockIdx.x;
  if (bid == 0)
        ; // RunParallel (10, [] __device__ (int nr, int size) { ; });
  else
    Worker();
}

__host__ void StartDeviceJob()
{
  RunAll<<<10, 64>>>();
}

// #endif


namespace ngla
{
  
  ;
  
}
