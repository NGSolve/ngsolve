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


inline void gettimeofday(struct timeval* t,void* timezone)
{       struct _timeb timebuffer;
        _ftime( &timebuffer );
        t->tv_sec=timebuffer.time;
        t->tv_usec=1000*timebuffer.millitm;
}

inline double WallTime ()
{
  struct _timeb timebuffer;
  _ftime( &timebuffer );
  return timebuffer.time+1e-3*timebuffer.millitm;
}


#else

#include <sys/time.h>

inline double WallTime ()
{
  timeval time;
  gettimeofday (&time, 0);
  return time.tv_sec + 1e-6 * time.tv_usec;
}

#endif

#ifdef VTRACE
#include "vt_user.h"
#else
#define VT_USER_START(n)
#define VT_USER_END(n)
#define VT_TRACER(n)
#define VT_ON()
#define VT_OFF()
#endif

namespace ngstd
{


  /**
     A built-in profile
  */
  class NgProfiler
  {
    /// maximal number of timers
    enum { SIZE = 1000 };

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

    int total_timer;
  public: 
    /// create new profile
    NgProfiler();
    /// delete profiler
    ~NgProfiler();

    /// create new timer, use integer index
    NGS_DLL_HEADER static int CreateTimer (const string & name);


#ifndef NOPROFILE


#ifdef USE_TIMEOFDAY
    static void StartTimer (int nr) 
    { 
      timeval time;
      gettimeofday (&time, 0);
      starttimes[nr] = time.tv_sec + 1e-6 * time.tv_usec;
      counts[nr]++; 
      VT_USER_START (const_cast<char*> (names[nr].c_str())); 
    }

    static void StopTimer (int nr) 
    { 
      timeval time;
      gettimeofday (&time, 0);
      tottimes[nr] += time.tv_sec + 1e-6 * time.tv_usec - starttimes[nr];
      VT_USER_END (const_cast<char*> (names[nr].c_str())); 
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

    static double GetTime (int nr)
    {
#ifdef USE_TIMEOFDAY
      return tottimes[nr];
#else
      return tottimes[nr]/CLOCKS_PER_SEC;
#endif
    }

    static long int GetCounts (int nr)
    {
      return counts[nr];
    }

    /// change name
    static void SetName (int nr, const string & name) { names[nr] = name; }
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
  public:
    Timer (const string & name)
    {
      timernr = NgProfiler::CreateTimer (name);
    }
    void Start () 
    {
      NgProfiler::StartTimer (timernr);
    }
    void Stop () 
    {
      NgProfiler::StopTimer (timernr);
    }
    void AddFlops (double aflops)
    {
      NgProfiler::AddFlops (timernr, aflops);
    }

    double GetTime () { return NgProfiler::GetTime(timernr); }
    long int GetCounts () { return NgProfiler::GetCounts(timernr); }

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
#else


  

  class Timer
  {
    static Timer * stack_top;

    int timer_id;
    Timer * prev;

  public:
    Timer (const string & name)
    {
      timer_id = VT_USER_DEF(name.c_str());
    }
    void Start () 
    {
      prev = stack_top;
      stack_top = this;
      if (prev)
	VT_USER_END_ID (prev -> timer_id);
      VT_USER_START_ID(timer_id);
    }
    void Stop () 
    {
      VT_USER_END_ID(timer_id);
      if (prev != NULL)
	VT_USER_START_ID(prev -> timer_id);
      stack_top = prev;
    }
    void AddFlops (double aflops)
    {
      ; // NgProfiler::AddFlops (timernr, aflops);
    }

    double GetTime () { return 0; }
    long int GetCounts () { return 0; }
    // operator int () { return timernr; }
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
#endif

}


#endif
