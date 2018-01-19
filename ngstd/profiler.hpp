#ifndef FILE_PROFILER
#define FILE_PROFILER

/**************************************************************************/
/* File:   profiler.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   5. Jan. 2005                                                  */
/**************************************************************************/


// Philippose - 27 January 2010
// Windows does not provide a "sys/time.h" include, 
// and neither does it define the function "gettimeofday" 
// anywhere..... This is a workaround...
#ifdef _WIN32
#include <sys/timeb.h>
#include <sys/types.h>
#include <winsock.h>
#include <intrin.h>   // for __rdtsc()  CPU time step counter

inline void gettimeofday(struct timeval* t,void* timezone)
{       struct _timeb timebuffer;
        _ftime( &timebuffer );
        t->tv_sec=timebuffer.time;
        t->tv_usec=1000*timebuffer.millitm;
}

#else
#include <sys/time.h>
#include <x86intrin.h>   // for __rdtsc()  CPU time step counter
#endif

#include <chrono>

namespace ngstd {
  extern NGS_DLL_HEADER std::chrono::time_point<std::chrono::system_clock> wall_time_start;

  // Time in seconds since program start
  inline double WallTime ()
    {
      std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = now-wall_time_start;
      return elapsed_seconds.count();
    }
}

#ifdef VTRACE
#include "vt_user.h"
#else
#define VT_USER_START(n)
#define VT_USER_END(n)
#define VT_TRACER(n)
#define VT_ON()
#define VT_OFF()
#endif

#include "paje_interface.hpp"

namespace ngstd
{


  /**
     A built-in profile
  */
  class NgProfiler
  {
  public:
    /// maximal number of timers
    enum { SIZE = 8*1024 };
    //  static long int tottimes[SIZE];
    // static long int starttimes[SIZE];

    NGS_DLL_HEADER static double tottimes[SIZE];
    NGS_DLL_HEADER static double starttimes[SIZE];

    NGS_DLL_HEADER static long int counts[SIZE];
    NGS_DLL_HEADER static double flops[SIZE];
    NGS_DLL_HEADER static double loads[SIZE];
    NGS_DLL_HEADER static double stores[SIZE];
    NGS_DLL_HEADER static string names[SIZE];
    NGS_DLL_HEADER static int usedcounter[SIZE];

    NGS_DLL_HEADER static size_t * thread_times;
    NGS_DLL_HEADER static size_t * thread_flops;
  private:

    // int total_timer;
    static string filename;
  public: 
    /// create new profile
    NgProfiler();
    /// delete profiler
    ~NgProfiler();

    static void SetFileName (const string & afilename) { filename = afilename; }

    /// create new timer, use integer index
    NGS_DLL_HEADER static int CreateTimer (const string & name);

    NGS_DLL_HEADER static void Reset ();


#ifndef NOPROFILE


#ifdef USE_TIMEOFDAY
    static void StartTimer (int nr) 
    { 
      timeval time;
      gettimeofday (&time, 0);
      // starttimes[nr] = time.tv_sec + 1e-6 * time.tv_usec;
      // #pragma omp atomic
      AsAtomic(tottimes[nr]) += -(time.tv_sec + 1e-6 * time.tv_usec);
      // #pragma omp atomic
      AsAtomic(counts[nr])++; 
      VT_USER_START (const_cast<char*> (names[nr].c_str())); 
    }

    static void StopTimer (int nr) 
    { 
      timeval time;
      gettimeofday (&time, 0);
      // tottimes[nr] += time.tv_sec + 1e-6 * time.tv_usec - starttimes[nr];
      // #pragma omp atomic
      AsAtomic(tottimes[nr]) += time.tv_sec + 1e-6 * time.tv_usec;
      VT_USER_END (const_cast<char*> (names[nr].c_str())); 
    }

    static void StartThreadTimer (size_t nr, size_t tid)
    {
      thread_times[tid*SIZE+nr] -= __rdtsc();
    }

    static void StopThreadTimer (size_t nr, size_t tid)
    {
      thread_times[tid*SIZE+nr] += __rdtsc();
    }

    static void AddThreadFlops (size_t nr, size_t tid, size_t flops)
    {
      thread_flops[tid*SIZE+nr] += flops;
    }
      
#else
  
    /// start timer of index nr
    static void StartTimer (int nr) 
    {
      starttimes[nr] = clock(); counts[nr]++; 
      VT_USER_START (const_cast<char*> (names[nr].c_str())); 
    }

    /// stop timer of index nr
    static void StopTimer (int nr) 
    { 
      tottimes[nr] += clock()-starttimes[nr]; 
      VT_USER_END (const_cast<char*> (names[nr].c_str())); 
    }

    static void StartThreadTimer (size_t nr, size_t tid)
    {
      thread_times[tid*SIZE+nr] -= __rdtsc();
    }

    static void StopThreadTimer (size_t nr, size_t tid)
    {
      thread_times[tid*SIZE+nr] += __rdtsc();
    }

    static void AddThreadFlops (size_t nr, size_t tid, size_t flops)
    {
      thread_flops[tid*SIZE+nr] += flops;
    }
    
#endif


    /// if you know number of flops, provide them to obtain the MFlop - rate
    static void AddFlops (int nr, double aflops) { flops[nr] += aflops; }
    static void AddLoads (int nr, double aloads) { loads[nr] += aloads; }
    static void AddStores (int nr, double astores) { stores[nr] += astores; }
#else

    static void StartTimer (int nr) { ; }
    static void StopTimer (int nr) { ; }
    static void AddFlops (int nr, double aflops) { ; };
    static void AddLoads (int nr, double aflops) { ; };
    static void AddStores (int nr, double aflops) { ; };
#endif

    static int GetNr (const string & name)
    {
      for (int i = SIZE-1; i >= 0; i--)
        if (names[i] == name)
          return i;
      return -1;
    }

    static double GetTime (int nr)
    {
#ifdef USE_TIMEOFDAY
      return tottimes[nr];
#else
      return tottimes[nr]/CLOCKS_PER_SEC;
#endif
    }

    static double GetTime (const string & name)
    {
      for (int i = SIZE-1; i >= 0; i--)
        if (names[i] == name)
          return GetTime (i);
      return 0;
    }


    static long int GetCounts (int nr)
    {
      return counts[nr];
    }

    static long int GetFlops (int nr)
    {
      return flops[nr];
    }

    /// change name
    static void SetName (int nr, const string & name) { names[nr] = name; }
    static string GetName (int nr) { return names[nr]; }
    /// print profile
    NGS_DLL_HEADER static void Print (FILE * ost);
    //static void Print (ostream & ost);

    /**
       Timer object.
       Start / stop timer at constructor / destructor.
    */
    class RegionTimer
    {
      int nr;
    public:
      /// start timer
      RegionTimer (int anr) : nr(anr) { StartTimer (nr); }
      /// stop timer
      ~RegionTimer () { StopTimer (nr); }
    };
  };



  
#ifndef VTRACE
  class Timer
  {
    int timernr;
    int priority;
  public:
    Timer (const string & name, int apriority = 1)
      : priority(apriority)
    {
      timernr = NgProfiler::CreateTimer (name);
    }
    void SetName (const string & name)
    {
      NgProfiler::SetName (timernr, name);
    }
    void Start () 
    {
      if (priority <= 2) 
	NgProfiler::StartTimer (timernr);
      if (priority <= 1)
        if(trace) trace->StartTimer(timernr);
    }
    void Stop () 
    {
      if (priority <= 2) 
	NgProfiler::StopTimer (timernr);
      if (priority <= 1)
        if(trace) trace->StopTimer(timernr);
    }
    void AddFlops (double aflops)
    {
      if (priority <= 2)
	NgProfiler::AddFlops (timernr, aflops);
    }

    double GetTime () { return NgProfiler::GetTime(timernr); }
    long int GetCounts () { return NgProfiler::GetCounts(timernr); }
    double GetMFlops () 
    { return NgProfiler::GetFlops(timernr) 
        / NgProfiler::GetTime(timernr) * 1e-6; }
    operator int () { return timernr; }
  };
  
  
  /**
     Timer object.
       Start / stop timer at constructor / destructor.
  */
  class RegionTimer
  {
    Timer & timer;
  public:
    /// start timer
    RegionTimer (Timer & atimer) : timer(atimer) { timer.Start(); }
    /// stop timer
    ~RegionTimer () { timer.Stop(); }
  };

  class ThreadRegionTimer
  {
    size_t nr;
    size_t tid;
  public:
    /// start timer
    ThreadRegionTimer (size_t _nr, size_t _tid) : nr(_nr), tid(_tid)
    { NgProfiler::StartThreadTimer(nr, tid); }
    /// stop timer
    ~ThreadRegionTimer ()
    { NgProfiler::StopThreadTimer(nr, tid); }
  };

  class RegionTracer
    {
      int nr;
      int thread_id;
    public:
      static constexpr int ID_JOB = PajeTrace::Task::ID_JOB;
      static constexpr int ID_NONE = PajeTrace::Task::ID_NONE;
      static constexpr int ID_TIMER = PajeTrace::Task::ID_TIMER;

      /// start trace
      RegionTracer (int athread_id, int region_id, int id_type = ID_NONE, int additional_value = -1 )
        : thread_id(athread_id)
        {
          nr = trace->StartTask (athread_id, region_id, id_type, additional_value);
        }
      /// start trace with timer
      RegionTracer (int athread_id, Timer & timer, int additional_value = -1 )
        : thread_id(athread_id)
        {
          nr = trace->StartTask (athread_id, (int)timer, ID_TIMER, additional_value);
        }

      /// set user defined value
      void SetValue( int additional_value )
      {
        trace->SetTask( thread_id, nr, additional_value );
      }

      /// stop trace
      ~RegionTracer ()
        {
          trace->StopTask (thread_id, nr);
        }
    };
#else


  
#ifdef PARALLEL
  class Timer
  {
    static Timer * stack_top;

    int timer_id;
    Timer * prev;
    int priority;
    string name;
  public:
    Timer (const string & aname, int apriority = 1)
      : name(aname)
    {
      priority = apriority;
      // timer_id = VT_USER_DEF(name.c_str());
      timer_id = VT_USER_DEF(name.c_str(), NULL, NULL, -1);
    }
    void Start () 
    {
      if (priority == 1)
	{
	  prev = stack_top;
	  stack_top = this;
	  if (prev)
	    // VT_USER_END_ID (prev -> timer_id);
	    // VT_USER_END2 (prev -> timer_id);
	    VT_USER_END (prev -> name.c_str());
	  // VT_USER_START_ID(timer_id);
	  // VT_USER_START2 (timer_id);
	  VT_USER_START (name.c_str());
	}
    }
    void Stop () 
    {
      if (priority == 1)
	{
	  // VT_USER_END_ID(timer_id);
	  // VT_USER_END2(timer_id);
	  VT_USER_END(name.c_str());
	  if (prev != NULL)
	    // VT_USER_START_ID(prev -> timer_id);
	    // VT_USER_START2(prev -> timer_id);
	    VT_USER_START(prev -> name.c_str());
	  stack_top = prev;
	}
    }
    void SetName (const string & /* st */) { ; }
    void AddFlops (double aflops)  { ; }
    double GetTime () { return 0; }
    long int GetCounts () { return 0; }
    operator int () { return timer_id; }
  };

  class RegionTracer
    {
    public:
      static constexpr int ID_JOB = -1;

      /// start trace
      RegionTracer (int athread_id, int region_id, int id_type = -1, int additional_value = -1 ){}
      /// start trace with timer
      RegionTracer (int athread_id, Timer & timer, int additional_value = -1 ){}
      /// set user defined value
      void SetValue( int additional_value ){}
      /// stop trace
      ~RegionTracer (){}
    };

#else

  class Timer
  {
    int timer_id;

  public:
    Timer (const string & name, int priority = 1)
    {
      timer_id = VT_USER_DEF(name.c_str());
    }
    void Start () 
    {
      VT_USER_START_ID(timer_id);
    }
    void Stop () 
    {
      VT_USER_END_ID(timer_id);
    }

    void SetName (const string & st) { ; }
    void AddFlops (double aflops)  { ; }
    double GetTime () { return 0; }
    long int GetCounts () { return 0; }
  };

#endif



  
  /**
     Timer object.
       Start / stop timer at constructor / destructor.
  */
  class RegionTimer
  {
    Timer & timer;
  public:
    /// start timer
    RegionTimer (Timer & atimer) : timer(atimer) { timer.Start(); }
    /// stop timer
    ~RegionTimer () { timer.Stop(); }
  };
#endif

  // Helper function for timings
  // Run f() at least min_iterations times until max_time seconds elapsed
  // returns minimum runtime for a call of f()
  template<typename TFunc>
  double RunTiming( TFunc f, double max_time = 0.5, int min_iterations = 10 )
  {
      // Make sure the whole test run does not exceed maxtime
      double tend = WallTime()+max_time;

      // warmup
      f();

      double tres = std::numeric_limits<double>::max();
      int iteration = 0;
      while(WallTime()<tend || iteration++ < min_iterations)
      {
          double t = -WallTime();
          f();
          t += WallTime();
          tres = min2(tres, t);
      }

      return tres;
  }


}


#endif
