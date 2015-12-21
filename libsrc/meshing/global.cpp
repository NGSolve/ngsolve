#include <mystdlib.h>
#include "meshing.hpp"


namespace netgen
{
  // stringstream emptystr;
  // ostream * testout = &emptystr;
  // testout -> clear(ios::failbit);

  // ostream * testout = &cout;
  ostream * testout = new ostream(0);

  // NetgenOutStream * testout = new NetgenOutStream;

  ostream * mycout = &cout;
  ostream * myerr = &cerr;

  // some functions (visualization) still need a global mesh 
  DLL_HEADER shared_ptr<Mesh> mesh;
  DLL_HEADER shared_ptr<NetgenGeometry> ng_geometry;

  weak_ptr<Mesh> global_mesh;
  void SetGlobalMesh (shared_ptr<Mesh> m)
  {
    PrintMessage(5, "set global mesh");
    global_mesh = m;
  }
  

  
  //  Flags parameters;
  int silentflag = 0;
  int testmode = 0;

  volatile multithreadt multithread;

  string ngdir = ".";

  // parallel netgen
  int id = 0, ntasks = 1;


  void Ng_PrintDest(const char * s)
  {
    if (id == 0)
      (*mycout) << s << flush;
  }

  DLL_HEADER void MyError(const char * ch)
  {
    cout << ch;
    (*testout) << "Error !!! " << ch << endl << flush;
  }

  static clock_t starttimea;
  void ResetTime ()
  {
    starttimea = clock();
  }

  double GetTime ()
  {
    return double(clock() - starttimea) / CLOCKS_PER_SEC;
  }



  Array<int> tets_in_qualclass;

  int h_argc = 0;
  char ** h_argv = NULL;

  multithreadt :: multithreadt()
  {
    pause =0;
    testmode = 0;
    redraw = 0;
    drawing = 0;
    terminate = 0;
    running = 0;
    percent = 0;
    task = "";
  }

  DebugParameters debugparam;
  bool verbose = 0;

  int timestamp = 0;
  int GetTimeStamp() 
  { 
    return timestamp; 
  }

  int NextTimeStamp()
  {
    timestamp++;
    return timestamp;
  }
}
