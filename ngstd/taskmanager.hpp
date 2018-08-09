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

    // int node_nr;
    // int nnodes;
  };

  NGS_DLL_HEADER extern class TaskManager * task_manager;
  
  class TaskManager
  {
//     PajeTrace *trace;

    class alignas(64) NodeData : public AlignedAlloc<NodeData>
    {
    public:
      atomic<int> start_cnt{0};
      atomic<int> participate{0};
    };
    
    static const function<void(TaskInfo&)> * func;
    static const function<void()> * startup_function;
    static const function<void()> * cleanup_function;
    static atomic<int> ntasks;
    static Exception * ex;

    static atomic<int> jobnr;

    static atomic<int> complete[8];   // max nodes
    static atomic<int> done;
    static atomic<int> active_workers;
    static atomic<int> workers_on_node[8];   // max nodes
    // Array<atomic<int>*> sync;
    static int sleep_usecs;
    static bool sleep;

    static NodeData *nodedata[8];

    static int num_nodes;
    NGS_DLL_HEADER static int num_threads;
    NGS_DLL_HEADER static int max_threads;



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
    static int GetThreadId() { return thread_id; } 
    int GetNumNodes() const { return num_nodes; }

    static void SetPajeTrace (bool use)  { use_paje_trace = use; }
    
    NGS_DLL_HEADER static void CreateJob (const function<void(TaskInfo&)> & afunc, 
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
                           int antasks = TaskManager::GetNumThreads(),
                           TotalCosts costs = 1000)
  {
    // if (task_manager && costs() >= 1000)

    TaskManager::CreateJob 
        ([r, f] (TaskInfo & ti) 
         {
           auto myrange = r.Split (ti.task_nr, ti.ntasks);
           for (auto i : myrange) f(i);
         }, 
         antasks);

      /*
    else
      for (auto i : r) f(i);
      */
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
                                int antasks = TaskManager::GetNumThreads(),
                                TotalCosts costs = 1000)
  {
    // if (task_manager && costs() >= 1000)

    TaskManager::CreateJob 
        ([r, f] (TaskInfo & ti) 
         {
           auto myrange = r.Split (ti.task_nr, ti.ntasks);
           f(myrange);
         }, 
         antasks);
    /*
    else
      f(r);
    */
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
                           int antasks = TaskManager::GetNumThreads())
  {
    TaskManager::CreateJob (f, antasks);
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
   /*
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

  // IntRange Get()
  // {
  //   lock_guard<mutex> guard(lock);
  //   return IntRange(begin, end);
  // }

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


  // inline ostream & operator<< (ostream & ost, AtomicRange & r)
  // {
  //   ost << r.Get();
  //   return ost;
  // }
  */


   
   class alignas(4096) AtomicRange : public AlignedAlloc<AtomicRange>
  {
    atomic<size_t> begin;
    atomic<size_t> end;
  public:
    
    void Set (IntRange r)
    {
      begin.store(numeric_limits<size_t>::max(), memory_order_release);
      end.store(r.end(), memory_order_release);
      begin.store(r.begin(), std::memory_order_release);
    }
  
    void SetNoLock (IntRange r)
    {
      end.store(r.end(), memory_order_release);
      begin.store(r.begin(), std::memory_order_release);
    }

    // IntRange Get()
    // {
    //   lock_guard<mutex> guard(lock);
    //   return IntRange(begin, end);
    // }
    
    bool PopFirst (size_t & first)
    {
      first = begin++;
      return first < end;
      /*
      // int oldbegin = begin;
      size_t oldbegin = begin.load(std::memory_order_acquire);
      if (oldbegin >= end) return false;
      while (!begin.compare_exchange_weak (oldbegin, oldbegin+1,
                                           std::memory_order_relaxed, std::memory_order_relaxed))
        if (oldbegin >= end) return false;
      
      first = oldbegin;
      return true;
      */
    }
    
    bool PopHalf (IntRange & r)
    {
      // int oldbegin = begin;
      size_t oldbegin = begin.load(std::memory_order_acquire);
      size_t oldend = end.load(std::memory_order_acquire);
      if (oldbegin >= oldend) return false;
      
      // lock_guard<mutex> guard(lock);    
      while (!begin.compare_exchange_weak (oldbegin, (oldbegin+oldend+1)/2,
                                           std::memory_order_relaxed, std::memory_order_relaxed))
        {
          oldend = end.load(std::memory_order_acquire);
          if (oldbegin >= oldend) return false;
        }
      
      r = IntRange(oldbegin, (oldbegin+oldend+1)/2);
      return true;
    }
  };
  



  class SharedLoop2
  {
    Array<AtomicRange> ranges;
    atomic<size_t> processed;
    atomic<size_t> total;
    atomic<int> participants;
    
    class SharedIterator
    {
      FlatArray<AtomicRange> ranges;
      atomic<size_t> & processed;
      size_t total;
      size_t myval;
      size_t processed_by_me = 0;
      int me;
      int steal_from;
    public:
      SharedIterator (FlatArray<AtomicRange> _ranges, atomic<size_t> & _processed, size_t _total,
                      int _me, bool begin_it)
        : ranges(_ranges), processed(_processed), total(_total)
      {
        if (begin_it)
          {
            // me = TaskManager::GetThreadId();
            me = _me;
            steal_from = me;
            GetNext();
          }
      }
      ~SharedIterator()
      {
        if (processed_by_me)
          processed += processed_by_me;
      }
      
      SharedIterator & operator++ () { GetNext(); return *this;}

      void GetNext()
      {
        size_t nr;
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
      
      size_t operator* () const { return myval; }
      bool operator!= (const SharedIterator & it2) const { return processed < total; }
    };
    
    
  public:
    SharedLoop2 ()
      : ranges(TaskManager::GetNumThreads())
    { ; }
    
    SharedLoop2 (IntRange r)
      : ranges(TaskManager::GetNumThreads())
    {
      Reset (r);
    }

    void Reset (IntRange r)
    {
      for (size_t i = 0; i < ranges.Size(); i++)
        ranges[i].SetNoLock (r.Split(i,ranges.Size()));
      
      total.store(r.Size(), std::memory_order_relaxed);
      participants.store(0, std::memory_order_relaxed);
      processed.store(0, std::memory_order_release);
    }
    
    SharedIterator begin()
    {
      /*
      int me = participants++;
      if (me < ranges.Size())
        return SharedIterator (ranges, processed, total, me, true);
      else
        // more participants than buckets. set processed to total, and the loop is terminated immediately
        return SharedIterator (ranges, total, total, me, true);
      */
      return SharedIterator (ranges, processed, total, TaskManager::GetThreadId(), true);      
    }
    
    SharedIterator end()   { return SharedIterator (ranges, processed, total, -1, false); }
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

      /*
      size_t sum = 0;
      for (auto i : ngstd::Range(n))
        {
          sum += costs(i);
          prefix[i] = sum;
        }
      total_costs = sum;
      */
      
      Array<size_t> partial_sums(TaskManager::GetNumThreads()+1);
      partial_sums[0] = 0;
      ParallelJob
        ([&] (TaskInfo ti)
         {
           IntRange r = IntRange(n).Split(ti.task_nr, ti.ntasks);
           size_t mysum = 0;
           for (size_t i : r)
             {
               size_t c = costs(i);
               mysum += c;
               prefix[i] = c;
             }
           partial_sums[ti.task_nr+1] = mysum;
         });
      
      for (size_t i = 1; i < partial_sums.Size(); i++)
        partial_sums[i] += partial_sums[i-1];
      total_costs = partial_sums.Last();
      
      ParallelJob
        ([&] (TaskInfo ti)
         {
           IntRange r = IntRange(n).Split(ti.task_nr, ti.ntasks);
           size_t mysum = partial_sums[ti.task_nr];
           for (size_t i : r)
             {
               mysum += prefix[i];
               prefix[i] = mysum;
             }
         });
      

      part.SetSize (size+1);
      part[0] = 0;

      for (int i = 1; i <= size; i++)
        part[i] = BinSearch (prefix, total_costs*i/size);      
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





  template <typename FUNC, typename OP, typename T>
  auto ParallelReduce (size_t n, FUNC f, OP op, T initial1)
  {
    typedef decltype (op(initial1,initial1)) TRES;
    TRES initial(initial1);
    /*
    for (size_t i = 0; i < n; i++)
      initial = op(initial, f(i));
    */
    Array<TRES> part_reduce(TaskManager::GetNumThreads());
    ParallelJob ([&] (TaskInfo ti)
                 {
                   auto r = Range(n).Split(ti.task_nr, ti.ntasks);
                   auto var = initial;
                   for (auto i : r)
                     var = op(var, f(i));
                   part_reduce[ti.task_nr] = var;
                 });
    for (auto v : part_reduce)
      initial = op(initial, v);
    return initial;
  }





  

  //  some suggar for working with arrays 

  template <typename T> template <typename T2>
  const FlatArray<T> FlatArray<T>::operator= (ParallelValue<T2> val)
  {
    ParallelForRange (Size(),
                      [this, val] (IntRange r)
                      {
                        for (auto i : r)
                          (*this)[i] = val;
                      });
    return *this;
  }

  template <typename T> template <typename T2>
  const FlatArray<T> FlatArray<T>::operator= (ParallelFunction<T2> func)
  {
    ParallelForRange (Size(),
                      [this, func] (IntRange r)
                      {
                        for (auto i : r)
                          (*this)[i] = func(i);
                      });
    return *this;
  }

class Tasks
{
  size_t num;
public:
  explicit Tasks (size_t _num = TaskManager::GetNumThreads()) : num(_num) { ; }
  auto GetNum() const { return num; } 
};

/* currently not used, plus causing problems on MSVC 2017
template <typename T, typename std::enable_if<ngstd::has_call_operator<T>::value, int>::type = 0>                                  
inline ParallelFunction<T> operator| (const T & func, Tasks tasks)
{
  return func;
}

template <typename T, typename std::enable_if<!ngstd::has_call_operator<T>::value, int>::type = 0>                                  
inline ParallelValue<T> operator| (const T & obj, Tasks tasks)
{
  return obj;
}

inline Tasks operator "" _tasks_per_thread (unsigned long long n)
{
  return Tasks(n * TaskManager::GetNumThreads());
}
*/

/*
  thought to be used as:   array = 1 | tasks
class DefaultTasks
{
public:
  operator Tasks () const { return TaskManager::GetNumThreads(); }
};
static DefaultTasks tasks;
*/







#ifdef USE_NUMA

template <typename T>
class NumaInterleavedArray : public Array<T>
{
  T * numa_ptr;
  size_t numa_size;
public:
  NumaInterleavedArray () { numa_size = 0; numa_ptr = nullptr; }
  NumaInterleavedArray (size_t s)
    : Array<T> (s, (T*)numa_alloc_interleaved(s*sizeof(T)))
  {
    numa_ptr = this->data;
    numa_size = s;
  }

  ~NumaInterleavedArray ()
  {
    numa_free (numa_ptr, numa_size*sizeof(T));
  }

  NumaInterleavedArray & operator= (T val)
  {
    Array<T>::operator= (val);      
    return *this;
  }

  NumaInterleavedArray & operator= (NumaInterleavedArray && a2)
  {
    Array<T>::operator= ((Array<T>&&)a2);  
    ngstd::Swap (numa_ptr, a2.numa_ptr);
    ngstd::Swap (numa_size, a2.numa_size);
    return *this;
  }

  void Swap (NumaInterleavedArray & b)
  {
    Array<T>::Swap(b);    
    ngstd::Swap (numa_ptr, b.numa_ptr);
    ngstd::Swap (numa_size, b.numa_size);
  }

  void SetSize (size_t size)
  {
    cerr << "************************* NumaDistArray::SetSize not overloaded" << endl;
    Array<T>::SetSize(size);
  }
};

template <typename T>
class NumaDistributedArray : public Array<T>
{
  T * numa_ptr;
  size_t numa_size;
public:
  NumaDistributedArray () { numa_size = 0; numa_ptr = nullptr; }
  NumaDistributedArray (size_t s)
    : Array<T> (s, (T*)numa_alloc_local(s*sizeof(T)))
  {
    numa_ptr = this->data;
    numa_size = s;

    /* int avail = */ numa_available();   // initialize libnuma
    int num_nodes = numa_num_configured_nodes();
    size_t pagesize = numa_pagesize();
    
    int npages = ceil ( double(s)*sizeof(T) / pagesize );

    // cout << "size = " << numa_size << endl;
    // cout << "npages = " << npages << endl;

    for (int i = 0; i < num_nodes; i++)
      {
        int beg = (i * npages) / num_nodes;
        int end = ( (i+1) * npages) / num_nodes;
        // cout << "node " << i << " : [" << beg << "-" << end << ")" << endl;
        numa_tonode_memory(numa_ptr+beg*pagesize/sizeof(T), (end-beg)*pagesize, i);
      }
  }

  ~NumaDistributedArray ()
  {
    numa_free (numa_ptr, numa_size*sizeof(T));
  }

  NumaDistributedArray & operator= (NumaDistributedArray && a2)
  {
    Array<T>::operator= ((Array<T>&&)a2);  
    ngstd::Swap (numa_ptr, a2.numa_ptr);
    ngstd::Swap (numa_size, a2.numa_size);
    return *this;
  }

  void Swap (NumaDistributedArray & b)
  {
    Array<T>::Swap(b);    
    ngstd::Swap (numa_ptr, b.numa_ptr);
    ngstd::Swap (numa_size, b.numa_size);
  }

  void SetSize (size_t size)
  {
    cerr << "************************* NumaDistArray::SetSize not overloaded" << endl;
    Array<T>::SetSize(size);
  }
};



template <typename T>
class NumaLocalArray : public Array<T>
{
  T * numa_ptr;
  size_t numa_size;
public:
  NumaLocalArray () { numa_size = 0; numa_ptr = nullptr; }
  NumaLocalArray (size_t s)
    : Array<T> (s, (T*)numa_alloc_local(s*sizeof(T)))
  {
    numa_ptr = this->data;
    numa_size = s;
  }

  ~NumaLocalArray ()
  {
    numa_free (numa_ptr, numa_size*sizeof(T));
  }

  NumaLocalArray & operator= (T val)
  {
    Array<T>::operator= (val);      
    return *this;
  }
  
  NumaLocalArray & operator= (NumaLocalArray && a2)
  {
    Array<T>::operator= ((Array<T>&&)a2);  
    ngstd::Swap (numa_ptr, a2.numa_ptr);
    ngstd::Swap (numa_size, a2.numa_size);
    return *this;
  }

  void Swap (NumaLocalArray & b)
  {
    Array<T>::Swap(b);    
    ngstd::Swap (numa_ptr, b.numa_ptr);
    ngstd::Swap (numa_size, b.numa_size);
  }

  void SetSize (size_t size)
  {
    cerr << "************************* NumaDistArray::SetSize not overloaded" << endl;
    Array<T>::SetSize(size);
  }
};


#else

  template <typename T>
  using NumaDistributedArray = Array<T>;

  template <typename T> 
  using NumaInterleavedArray = Array<T>;
  
  template <typename T>
  using NumaLocalArray = Array<T>;
  
#endif

}



#endif
