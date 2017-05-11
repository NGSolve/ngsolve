#ifndef FILE_TASKMANAGER
#define FILE_TASKMANAGER

/*********************************************************************/
/* File:   taskmanager.hpp                                           */
/* Author: M. Hochsterger, J. Schoeberl                              */
/* Date:   10. Mar. 2015                                             */
/*********************************************************************/

#include <atomic>
#include <thread>



namespace ngstd
{
  class PajeTrace;

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

  NGS_DLL_HEADER extern class TaskManager * task_manager;
  
  class TaskManager
  {
//     PajeTrace *trace;

    class NodeData
    {
    public:
      atomic<int> start_cnt{0};
      // atomic<int> complete_cnt;
      atomic<int> participate{0};
      // atomic<int> participate_exit;

      // NodeData() : start_cnt(0), participate(0), participate_exit(0) { ; }
    };
    
    static const function<void(TaskInfo&)> * func;
    static const function<void()> * startup_function;
    static const function<void()> * cleanup_function;
    atomic<int> ntasks;
    Exception * ex;

    int jobnr;

    atomic<int> complete[8];   // max nodes
    atomic<int> done;
    atomic<int> active_workers;
    atomic<int> workers_on_node[8];   // max nodes

    int sleep_usecs;
    bool sleep;

    NodeData *nodedata[8];

    int num_nodes;
    static int num_threads;
    static int max_threads;
#ifndef __clang__    
    static thread_local int thread_id;
#else
    static __thread int thread_id;
#endif
    
    static bool use_paje_trace;
  public:
    
    TaskManager();
    ~TaskManager();


    void StartWorkers();
    void StopWorkers();

    void SuspendWorkers(int asleep_usecs = 1000 )
      {
        sleep_usecs = asleep_usecs;
        sleep = true;
      }
    void ResumeWorkers() { sleep = false; }

    static void SetNumThreads(int amax_threads);
    static int GetMaxThreads() { return max_threads; }
    // static int GetNumThreads() { return task_manager ? task_manager->num_threads : 1; }
    static int GetNumThreads() { return num_threads; }
    static int GetThreadId() { return task_manager ? task_manager->thread_id : 0; }
    int GetNumNodes() const { return num_nodes; }

    static void SetPajeTrace (bool use)  { use_paje_trace = use; }
    
    NGS_DLL_HEADER void CreateJob (const function<void(TaskInfo&)> & afunc, 
                    int antasks = task_manager->GetNumThreads());

    static void SetStartupFunction (const function<void()> & func) { startup_function = &func; }
    static void SetStartupFunction () { startup_function = nullptr; }
    static void SetCleanupFunction (const function<void()> & func) { cleanup_function = &func; }
    static void SetCleanupFunction () { cleanup_function = nullptr; }    

    void Done() { done = true; }
    void Loop(int thread_num);

    static list<tuple<string,double>> Timing ();
  };








  
  void RunWithTaskManager (function<void()> alg);

  // For Python context manager
  int  EnterTaskManager ();
  void ExitTaskManager (int num_threads);

  INLINE int TasksPerThread (int tpt)
  {
    // return task_manager ? tpt*task_manager->GetNumThreads() : 1;
    return tpt*TaskManager::GetNumThreads();
  }
  

  class TotalCosts
  {
    size_t cost;
  public:
    TotalCosts (size_t _cost) : cost(_cost) { ; }
    size_t operator ()() { return cost; }
  };

  template <typename TR, typename TFUNC>
  INLINE void ParallelFor (T_Range<TR> r, TFUNC f, 
                           int antasks = task_manager ? task_manager->GetNumThreads() : 0,
                           TotalCosts costs = 1000)
  {
    if (task_manager && costs() >= 1000)

      task_manager -> CreateJob 
        ([r, f] (TaskInfo & ti) 
         {
           auto myrange = r.Split (ti.task_nr, ti.ntasks);
           for (auto i : myrange) f(i);
         }, 
         antasks);

    else

      for (auto i : r) f(i);
  }

  /*
  template <typename TFUNC>
  INLINE void ParallelFor (size_t n, TFUNC f, 
                           int antasks = task_manager ? task_manager->GetNumThreads() : 0)
  {
    ParallelFor (IntRange (n), f, antasks);
  }
  */
  template <typename ...Args>
  INLINE void ParallelFor (size_t n, Args...args)
  {
    ParallelFor (IntRange (n), args...);
  }
  
  template <typename TR, typename TFUNC>
  INLINE void ParallelForRange (T_Range<TR> r, TFUNC f, 
                                int antasks = task_manager ? task_manager->GetNumThreads() : 0,
                                TotalCosts costs = 1000)
  {
    if (task_manager && costs() >= 1000)

      task_manager -> CreateJob 
        ([r, f] (TaskInfo & ti) 
         {
           auto myrange = r.Split (ti.task_nr, ti.ntasks);
           f(myrange);
         }, 
         antasks);

    else

      f(r);
  }

  /*
  template <typename TFUNC>
  INLINE void ParallelForRange (size_t n, TFUNC f, 
                                int antasks = task_manager ? task_manager->GetNumThreads() : 0)
  {
    ParallelForRange (IntRange(n), f, antasks);
  }
  */
  template <typename ...Args>
  INLINE void ParallelForRange (size_t n, Args...args)
  {
    ParallelForRange (IntRange(n), args...);
  }
  
  template <typename TFUNC>
  INLINE void ParallelJob (TFUNC f, 
                           int antasks = task_manager ? task_manager->GetNumThreads() : 1)
  {
    if (task_manager)

      task_manager -> CreateJob (f, antasks);

    else
      
      {
        TaskInfo ti;
        // ti.task_nr = 0; ti.ntasks = 1;
        ti.ntasks = antasks;
        ti.thread_nr = 0; ti.nthreads = 1;
        ti.node_nr = 0; ti.nnodes = 1;
        for (ti.task_nr = 0; ti.task_nr < antasks; ti.task_nr++)
          f(ti);
      }
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


  /*
class alignas(4096) AtomicRange
{
  mutex lock;
  int begin;
  int end;
public:
  
  void Set (IntRange r)
  {
    lock_guard<mutex> guard(lock);
    begin = r.begin();
    end = r.end();
  }

  IntRange Get()
  {
    lock_guard<mutex> guard(lock);
    return IntRange(begin, end);
  }

  bool PopFirst (int & first)
  {
    lock_guard<mutex> guard(lock);
    bool non_empty = end > begin;
    first = begin;
    if (non_empty) begin++;
    return non_empty;
  }

  bool PopHalf (IntRange & r)
  {
    lock_guard<mutex> guard(lock);
    bool non_empty = end > begin;
    if (non_empty)
      {
        int mid = (begin+end+1)/2;
        r = IntRange(begin, mid);
        begin = mid;
      }
    return non_empty;
  }
};
*/



  // lock free popfirst
  // faster for large loops, bug slower for small loops (~1000) ????

  class alignas(4096) AtomicRange
{
  mutex lock;
  atomic<int> begin;
  int end;
public:
  
  void Set (IntRange r)
  {
    lock_guard<mutex> guard(lock);
    // begin = r.begin();
    begin.store(r.begin(), std::memory_order_relaxed);
    end = r.end();
  }
  
  void SetNoLock (IntRange r)
  {
    begin.store(r.begin(), std::memory_order_relaxed);
    end = r.end();
  }

  IntRange Get()
  {
    lock_guard<mutex> guard(lock);
    return IntRange(begin, end);
  }

  bool PopFirst (int & first)
  {
    // int oldbegin = begin;
    int oldbegin = begin.load(std::memory_order_relaxed);
    if (oldbegin >= end) return false;
    while (!begin.compare_exchange_weak (oldbegin, oldbegin+1,
                                         std::memory_order_relaxed, std::memory_order_relaxed))
      if (oldbegin >= end) return false;        

    first = oldbegin;
    return true;
  }

  bool PopHalf (IntRange & r)
  {
    // int oldbegin = begin;
    int oldbegin = begin.load(std::memory_order_relaxed);    
    if (oldbegin >= end) return false;
    
    lock_guard<mutex> guard(lock);    
    while (!begin.compare_exchange_weak (oldbegin, (oldbegin+end+1)/2,
                                         std::memory_order_relaxed, std::memory_order_relaxed))
      if (oldbegin >= end) return false;        

    r = IntRange(oldbegin, (oldbegin+end+1)/2);
    return true;
  }
};

  


  inline ostream & operator<< (ostream & ost, AtomicRange & r)
  {
    ost << r.Get();
    return ost;
  }




  class SharedLoop2
  {
    Array<AtomicRange> ranges;
    atomic<int> processed;
    int total;
    
    class SharedIterator
    {
      FlatArray<AtomicRange> ranges;
      atomic<int> & processed;
      int total;
      int myval;
      int processed_by_me = 0;
      int me;
      int steal_from;
    public:
      SharedIterator (FlatArray<AtomicRange> _ranges, atomic<int> & _processed, int _total, bool begin_it)
        : ranges(_ranges), processed(_processed), total(_total)
      {
        me = TaskManager::GetThreadId();
        steal_from = me;
        if (begin_it)
          GetNext();
      }
      
      SharedIterator & operator++ () { GetNext(); return *this;}

      /*
      void GetNext()
      {
        while (1)
          {
            int nr;
            if (ranges[me].PopFirst(nr))
              {
                processed_by_me++;
                myval = nr;
                return;
              }

            processed += processed_by_me;
            processed_by_me = 0;

            // done with my work, going to steal ...
            while (1)
              {
                if (processed >= total) return;
                // steal_from = (steal_from + 1) % ranges.Size();
                steal_from++;
                if (steal_from == ranges.Size()) steal_from = 0;
            
                // steal half of the work reserved for 'from':
                IntRange steal;
                if (ranges[steal_from].PopHalf(steal))
                  {
                    // ranges[me].Set(steal);
                    // break;
                    myval = steal.First();
                    processed_by_me++;                    
                    if (myval+1 < steal.Next())
                      ranges[me].Set (IntRange(myval+1, steal.Next()));
                    return;
                  }
              }
          }
      }
      */
      
      void GetNext()
      {
        int nr;
        if (ranges[me].PopFirst(nr))
          {
            processed_by_me++;
            myval = nr;
            return;
          }
        GetNext2();
      }

      void GetNext2()
      {
        processed += processed_by_me;
        processed_by_me = 0;
        
        // done with my work, going to steal ...
        while (1)
          {
            if (processed >= total) return;
            // steal_from = (steal_from + 1) % ranges.Size();
            steal_from++;
            if (steal_from == ranges.Size()) steal_from = 0;
            
            // steal half of the work reserved for 'from':
            IntRange steal;
            if (ranges[steal_from].PopHalf(steal))
              {
                myval = steal.First();
                processed_by_me++;                    
                if (myval+1 < steal.Next())
                  ranges[me].Set (IntRange(myval+1, steal.Next()));
                return;
              }
          }
      }
      
      int operator* () const { return myval; }
      bool operator!= (const SharedIterator & it2) const { return processed < total; }
    };
    
    
  public:
    SharedLoop2 (IntRange r)
      : ranges(TaskManager::GetMaxThreads()), processed(0)
    {
      total = r.Size();
      for (size_t i = 0; i < ranges.Size(); i++)
        ranges[i].SetNoLock (r.Split(i,ranges.Size()));
    }
    
    SharedIterator begin() { return SharedIterator (ranges, processed, total, true); }
    SharedIterator end()   { return SharedIterator (ranges, processed, total, false); }
  };





  class Partitioning
  {
    Array<size_t> part;
    size_t total_costs;
  public:
    Partitioning () { ; }

    template <typename T>
    Partitioning (const Array<T> & apart) { part = apart; }

    template <typename T>
    Partitioning & operator= (const Array<T> & apart) { part = apart; return *this; }

    size_t GetTotalCosts() const { return total_costs; }

    template <typename TFUNC>
    void Calc (size_t n, TFUNC costs, int size = task_manager ? task_manager->GetNumThreads() : 1)
    {
      Array<size_t> prefix (n);

      size_t sum = 0;
      for (auto i : ngstd::Range(n))
        {
          sum += costs(i);
          prefix[i] = sum;
        }
      total_costs = sum;
      part.SetSize (size+1);
      part[0] = 0;

      for (int i = 1; i <= size; i++)
        part[i] = BinSearch (prefix, sum*i/size);      
    }
    
    size_t Size() const { return part.Size()-1; }
    IntRange operator[] (size_t i) const { return ngstd::Range(part[i], part[i+1]); }
    IntRange Range() const { return ngstd::Range(part[0], part[Size()]); }




  private:
    template <typename Tarray>
    int BinSearch(const Tarray & v, size_t i) {
      int n = v.Size();
      if (n == 0) return 0;
      
      int first = 0;
      int last = n-1;
      if(v[0]>i) return 0;
      if(v[n-1] <= i) return n;
      while(last-first>1) {
        int m = (first+last)/2;
        if(v[m]<i)
          first = m;
        else
          last = m;
      }
      return first;
    }
  };

  
  inline ostream & operator<< (ostream & ost, const Partitioning & part)
  {
    for (int i : Range(part.Size()))
      ost << part[i] << " ";
    return ost;
  }
  

  // tasks must be a multiple of part.size
  template <typename TFUNC>
  INLINE void ParallelFor (const Partitioning & part, TFUNC f, int tasks_per_thread = 1)
  {
    if (task_manager)
      {
        int ntasks = tasks_per_thread * task_manager->GetNumThreads();
        if (ntasks % part.Size() != 0)
          throw Exception ("tasks must be a multiple of part.size");

        task_manager -> CreateJob 
          ([&] (TaskInfo & ti) 
           {
             int tasks_per_part = ti.ntasks / part.Size();
             int mypart = ti.task_nr / tasks_per_part;
             int num_in_part = ti.task_nr % tasks_per_part;
             
             auto myrange = part[mypart].Split (num_in_part, tasks_per_part);
             for (auto i : myrange) f(i);
           }, ntasks);
      }
    else
      {
        for (auto i : part.Range())
          f(i);
      }
  }





  template <typename TFUNC>
  INLINE void ParallelForRange (const Partitioning & part, TFUNC f,
                                int tasks_per_thread = 1, TotalCosts costs = 1000)
  {
    if (task_manager && costs() >= 1000)
      {
        int ntasks = tasks_per_thread * task_manager->GetNumThreads();
        if (ntasks % part.Size() != 0)
          throw Exception ("tasks must be a multiple of part.size");

        task_manager -> CreateJob 
          ([&] (TaskInfo & ti) 
           {
             int tasks_per_part = ti.ntasks / part.Size();
             int mypart = ti.task_nr / tasks_per_part;
             int num_in_part = ti.task_nr % tasks_per_part;
             
             auto myrange = part[mypart].Split (num_in_part, tasks_per_part);
             f(myrange);
           }, ntasks);
      }
    else
      {
        f(part.Range());
      }
  }





}



#endif
