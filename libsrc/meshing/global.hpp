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
  DLL_HEADER extern void ResetTime ();

  ///
  DLL_HEADER extern int testmode;

  /// calling parameters
  // extern Flags parameters;

  // extern DLL_HEADER MeshingParameters mparam;

  DLL_HEADER extern Array<int> tets_in_qualclass;

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
    string * tcl_todo = new string("");  // tcl commands set from parallel thread
    multithreadt();
  };

  DLL_HEADER extern volatile multithreadt multithread;

  DLL_HEADER extern string ngdir;
  DLL_HEADER extern DebugParameters debugparam;
  DLL_HEADER extern bool verbose;

  DLL_HEADER extern int h_argc;
  DLL_HEADER extern char ** h_argv;


  DLL_HEADER extern weak_ptr<Mesh> global_mesh;
  DLL_HEADER void SetGlobalMesh (shared_ptr<Mesh> m);
}

#endif
