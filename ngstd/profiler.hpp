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
#else

#include <sys/time.h>

#endif

  /*
    #ifdef USE_TIMEOFDAY
    #include <sys/time.h>
    #endif
  */


#ifdef VTRACE
#include "vt_user.h"
#else
#define VT_USER_START(n)
#define VT_USER_END(n)
#define VT_TRACER(n)
#endif

namespace ngstd
{


  /**
     A built-in profile
  */
  class NGS_DLL_HEADER NgProfiler
  {
    /// maximal number of timers
    enum { SIZE = 1000 };

    //  static long int tottimes[SIZE];
    // static long int starttimes[SIZE];

    static double tottimes[SIZE];
    static double starttimes[SIZE];

    static long int counts[SIZE];
    static double flops[SIZE];
    static string names[SIZE];
    static int usedcounter[SIZE];

    int total_timer;
  public: 
    /// create new profile
    NgProfiler();
    /// delete profiler
    ~NgProfiler();

    /// create new timer, use integer index
    static int CreateTimer (const string & name);


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
#else

    static void StartTimer (int nr) { ; }
    static void StopTimer (int nr) { ; }
    static void AddFlops (int nr, double aflops) { ; };
#endif

    static double GetTime (int nr)
    {
#ifdef USE_TIMEOFDAY
      return tottimes[nr];
#else
      return tottimes[nr]/CLOCKS_PER_SEC;
#endif
    }

    /// change name
    static void SetName (int nr, const string & name) { names[nr] = name; }
    /// print profile
    static void Print (FILE * ost);
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

}


#endif
