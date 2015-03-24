#ifndef FILE_TASKMANAGER
#define FILE_TASKMANAGER

/*********************************************************************/
/* File:   taskmanager.hpp                                           */
/* Author: M. Hochsterger, J. Schoeberl                              */
/* Date:   10. Mar. 2015                                             */
/*********************************************************************/




namespace ngstd
{

  class TaskInfo
  {
  public:
    int task_nr;
    int ntasks;

    int thread_nr;
    int nthreads;

    int node_nr;
    int nnodes;
  };

  
  class TaskManager
  {

    class NodeData
    {
    public:
      atomic<int> start_cnt;
      atomic<int> complete_cnt;
      atomic<int> participate;

      NodeData() : start_cnt(0), participate(0) { ; }
    };
    
    function<void(TaskInfo&)> func;
    volatile int ntasks;
    Exception * ex;

    atomic<int> jobnr;

    atomic<int> complete[8];   // max nodes
    atomic<int> done;

    NodeData *nodedata[8];

    int num_nodes;
    
  public:
    
    TaskManager();

    void CreateJob (function<void(TaskInfo&)> afunc, 
                    int antasks = omp_get_max_threads());



    template <typename TFUNC>
    INLINE void ParallelFor (IntRange r, TFUNC f)
    {
      CreateJob 
        ([r, f] (TaskInfo & ti) 
         {
           auto myrange = r.Split (ti.task_nr, ti.ntasks);
           for (auto i : myrange) f(i);
         });
    }
    


    void Done() { done = true; }


    void Loop();
  };








  extern TaskManager * task_manager;
  
  void RunWithTaskManager (function<void()> alg);



  template <typename TFUNC>
  INLINE void ParallelFor (IntRange r, TFUNC f)
  {
    if (task_manager)
      task_manager -> ParallelFor (r, f);
    else
      for (auto i : r) f(i);
  }



  
  
  
  /*
    Usage example:

    ShareLoop myloop(100);
    task_manager->CreateJob ([]()
    {
      for (int i : myloop)
        cout << "i = " << i << endl;
    });

  */
  
  class SharedLoop
  {
    atomic<int> cnt;
    IntRange r;

    
    class SharedIterator
    {
      atomic<int> & cnt;
      int myval;
      int endval;
    public:
      SharedIterator (atomic<int> & acnt, int aendval, bool begin_iterator) 
        : cnt (acnt)
      {
        endval = aendval;
        myval = begin_iterator ? cnt++ : endval;
        if (myval > endval) myval = endval;
      }
      
      SharedIterator & operator++ () 
      {
        myval = cnt++; 
        if (myval > endval) myval = endval;
        return *this; 
      }
      
      int operator* () const { return myval; }
      bool operator!= (const SharedIterator & it2) const { return myval != it2.myval; }
    };
    
    
  public:
    SharedLoop (IntRange ar) : r(ar) { cnt = r.begin(); }
    SharedIterator begin() { return SharedIterator (cnt, r.end(), true); }
    SharedIterator end()   { return SharedIterator (cnt, r.end(), false); }
  };






}



#endif
