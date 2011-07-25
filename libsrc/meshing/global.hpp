#ifndef FILE_GLOBAL
#define FILE_GLOBAL


/**************************************************************************/
/* File:   global.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Okt. 95                                                    */
/**************************************************************************/

/*
  global functions and variables
*/

namespace netgen
{

  ///
  DLL_HEADER extern double GetTime ();
  extern void ResetTime ();

  ///
  extern int testmode;

  /// calling parameters
  // extern Flags parameters;

  // extern DLL_HEADER MeshingParameters mparam;

  extern Array<int> tets_in_qualclass;

  class multithreadt
  {
  public:
    int pause;
    int testmode;
    int redraw;
    int drawing;
    int terminate;
    int running;
    double percent;
    const char * task;
    bool demorunning;
    multithreadt();
  };

  extern volatile multithreadt multithread;

  extern string ngdir;
  extern DebugParameters debugparam;
  extern bool verbose;  
}

#endif
