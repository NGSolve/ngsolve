/**************************************************************************/
/* File:   ngsolve.cpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 2000                                                  */
/**************************************************************************/

/*
   Finite element package NGSolve
*/

#include <solve.hpp>
#include <parallelngs.hpp>

#include <thread>
#include <meshing/visual_interface.hpp>

using namespace std;
using namespace ngsolve;

typedef void *ClientData;

#include <nginterface.h>

namespace netgen
{
    DLL_HEADER extern string ngdir;
    DLL_HEADER extern bool netgen_executable_started;
}

using netgen::Ng_Tcl_SetResult;
using netgen::Ng_Tcl_CreateCommand;
using netgen::NG_TCL_OK;
using netgen::NG_TCL_ERROR;
using netgen::h_argc;
using netgen::h_argv;
/*
using netgen::NG_TCL_STATIC;
using netgen::NG_TCL_VOLATILE;
*/
#define NG_TCL_VOLATILE		((Tcl_FreeProc *) 1)
#define NG_TCL_STATIC		((Tcl_FreeProc *) 0)
#define NG_TCL_DYNAMIC		((Tcl_FreeProc *) 3)


#ifdef SOCKETS
#include "markus/jobmanager.hpp"
#endif


#ifdef NGS_PYTHON
// #include <dlfcn.h>  // dlopen of python lib ???
#include "../ngstd/python_ngstd.hpp"

extern PythonEnvironment pyenv;

void * ptr = (void*)PyOS_InputHook;



std::thread::id pythread_id = std::this_thread::get_id();
std::thread::id mainthread_id = std::this_thread::get_id();


void SpawnPython ()
{
  if(pythread_id != mainthread_id) 
    {
      cout << "Python thread already running!" << endl;
    }
  else
    {
      std::thread([]()
                  {
                    AcquireGIL gil_lock;
                    try{
                      Ng_SetRunning (1); 
                      pythread_id = std::this_thread::get_id();
                      pyenv.exec(
                                 "import ngsolve.__console;"
                                 "_vars2 = globals();"
                                 "_vars2.update(locals());"
                                 "ngsolve.__console.startConsole(_vars2)"
                      );
                      Ng_SetRunning (0); 
                    }
                    catch (py::error_already_set const &) {
                      PyErr_Print();
                    }
                    cout << "Python shell finished." << endl;
                    
                    pythread_id = mainthread_id;
                    
                  }).detach();
      /*
      cout << IM(1)
          << "ngs-python objects are available in ngstd, bla, ...\n"
          << "to import the whole bunch of objets enter\n\n"
          << "from ngsolve.ngstd import *\n"
          << "from ngsolve.bla import *\n"
          << "from ngsolve.fem import *\n"
          << "from ngsolve.la import *\n"
          << "from ngsolve.comp import *\n"
          << "from ngsolve.solve import *\n"
          //              << "from ngsolve.ngmpi import *\n"
          << "dir()\n"
          << endl << endl;
      */
#ifdef PARALLEL
      cout << IM(1) << "To start the mpi shell call" << endl << "MpiShell()" << endl << endl;
#endif
    }
}

#endif


// #include "/home/joachim/netgen-mesher/netgen/libsrc/include/meshing.hpp"
// volatile int & running = netgen::multithread.running;
// static bool got_exception = false;


int NGS_PrintRegistered (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, const char *argv[])
{
  ngfem::GetIntegrators().Print (cout);
  ngsolve::GetFESpaceClasses().Print (cout);
  ngsolve::GetPreconditionerClasses().Print (cout);

  return NG_TCL_OK;
}

int NGS_Help (ClientData clientData,
	      Tcl_Interp * interp,
	      int argc, const char *argv[])
{
  if (argc >= 2)
    {
      string topics = argv[1];

      stringstream str;

      if (topics == "constant")
	{
	  str << "heapsize = <num bytes>\n"
	      << "   size for optimized memory handler\n\n"
	      << "testout = <filename>\n"
	      << "   filename for testoutput\n\n"
	      << "numthreads = <num>\n"
	      << "   threads for openmp parallelization\n\n"
	      << "geometryorder = <num>\n"
	      << "   curved elements of this polynomial order\n\n"
	      << "refinep = 0|1\n"
	      << "   increase p instead of mesh refinement\n\n"
	      << "refinehp = 0|1\n"
	      << "   h-refinement only for singular elements, otherwise p\n\n"
	      << endl;
	}

      if (topics == "coefficient")
	{
	  str << "define coefficient <name>" << endl;
	  str << "val1, val2, val3, ...." << endl;
	}

      if (topics == "bilinearform")
	{
	  ;
	}


      cout << str.str();
      Ng_Tcl_SetResult (interp, (char*)str.str().c_str(), NG_TCL_VOLATILE);
    }
  return NG_TCL_OK;
}



#ifdef _NGSOLVE_SOCKETS_HPP
  class SocketOutArchive : public Archive
  {
    Socket & sock;
  public:
    SocketOutArchive (Socket & asock) : Archive(true), sock(asock)
    {
      auto versions = GetLibraryVersions();
      (*this) & versions;
    }

    const VersionInfo& GetVersion(const std::string& library)
    { return GetLibraryVersions().at(library); }

    using Archive::operator&;
    virtual Archive & operator & (std::byte & d)
    {
      sock.Tsend(uint8_t(d));
      return *this;
    }

    virtual Archive & operator & (float & f)
    {
      sock.Tsend(f);
      return *this;
    }
    virtual Archive & operator & (double & d) 
    {
      sock.Tsend(d);
      return *this; 
    }
    virtual Archive & operator & (int & i) 
    {
      sock.Tsend(i);
      return *this; 
    }
    virtual Archive & operator & (short int & i) 
    {
      sock.Tsend(i);
      return *this; 
    }
    virtual Archive & operator & (long & i) 
    {
      sock.Tsend(i);
      return *this; 
    }
    virtual Archive & operator & (size_t & i) 
    {
      sock.Tsend(i);
      return *this; 
    }
    virtual Archive & operator & (unsigned char & i) 
    { 
      sock.Tsend (i);
      return *this; 
    }
    virtual Archive & operator & (bool & b) 
    { 
      sock.Tsend (b);
      return *this; 
    }
    virtual Archive & operator & (string & str) 
    {
      sock.send (str);
      return *this; 
    }
    virtual Archive & operator & (char *& str) 
    {
      sock.send (string (str));
      return *this; 
    }
  };


  class SocketInArchive : public Archive
  {
    std::map<std::string, VersionInfo> vinfo{};
    Socket & sock;
  public:
    SocketInArchive (Socket & asock) : Archive(false), sock(asock)
    { (*this) & vinfo; }
    const VersionInfo& GetVersion(const std::string& library)
    { return vinfo[library]; }

    using Archive::operator&;
    virtual Archive & operator & (std::byte & d)
    {
      uint8_t tmp;
      sock.Trecv (tmp);
      d = byte(tmp);
      return *this;
    }

    virtual Archive & operator & (float & f)
    {
      sock.Trecv (f);
      return *this;
    }
    virtual Archive & operator & (double & d) 
    {
      sock.Trecv (d);
      return *this; 
    }
    virtual Archive & operator & (int & i) 
    {
      sock.Trecv (i);
      return *this; 
    }
    virtual Archive & operator & (short int & i) 
    {
      sock.Trecv (i);
      return *this; 
    }
    virtual Archive & operator & (unsigned char & i) 
    {
      sock.Trecv(i);
      return *this; 
    }
    virtual Archive & operator & (long & i) 
    {
      sock.Trecv (i);
      return *this; 
    }
    virtual Archive & operator & (size_t & i) 
    {
      sock.Trecv (i);
      return *this; 
    }
    virtual Archive & operator & (bool & b) 
    {
      sock.Trecv(b);
      return *this; 
    }

    virtual Archive & operator & (string & str) 
    { 
      sock.recv (str);
      return *this;
    }

    virtual Archive & operator & (char *& str) 
    { 
      string hstr;
      sock.recv (hstr);
      str = new char [hstr.length()+1];
      strcpy (str, hstr.c_str());
      return *this;
    }
  };


  class WorkerOutArchive : public Archive
  {
    // Socket & sock;
  public:
    WorkerOutArchive () : Archive(true) { ; }

    virtual Archive & operator & (std::byte & d)
    {
      return *this;
    }
    virtual Archive & operator & (float & f)
    {
      return *this;
    }
    virtual Archive & operator & (double & d) 
    {
      return *this; 
    }
    virtual Archive & operator & (int & i) 
    {
      return *this; 
    }
    virtual Archive & operator & (short int & i) 
    {
      return *this; 
    }
    virtual Archive & operator & (long int & i) 
    {
      return *this; 
    }
    virtual Archive & operator & (size_t & i) 
    {
      return *this; 
    }
    virtual Archive & operator & (unsigned char & i) 
    { 
      return *this; 
    }
    virtual Archive & operator & (bool & b) 
    { 
      return *this; 
    }
    virtual Archive & operator & (string & str) 
    {
      return *this; 
    }
    virtual Archive & operator & (char *& str) 
    {
      return *this; 
    }
  };


#endif

#ifndef WIN32
static pthread_t socket_thread;

void MyRunParallel ( void * (*fun)(void *), void * in)
{
  pthread_attr_t attr;
  pthread_attr_init (&attr);
  pthread_attr_setstacksize(&attr, 1000000);
  pthread_create (&socket_thread, &attr, fun, in);
}



#endif





int NGS_LoadPy (ClientData clientData,
		Tcl_Interp * interp,
		int argc, const char *argv[])
{
  if(!netgen::netgen_executable_started) {
      Ng_Tcl_SetResult (interp, (char*)"This feature is not available when running from Python", NG_TCL_STATIC);
      return NG_TCL_ERROR;
  }

  if (Ng_IsRunning())
    {
      Ng_Tcl_SetResult (interp, (char*)"Thread already running", NG_TCL_STATIC);
      return NG_TCL_ERROR;
    }

  if (argc >= 2)
    {
      try
	{
	  string filename = argv[1];
	  cout << IM(3) << "(should) load python file '" << filename << "'" << endl;

#ifdef NGS_PYTHON
	  {
        std::thread([](string init_file_)
          {
            AcquireGIL gil_lock;
            try{
              Ng_SetRunning (1); 
              pythread_id = std::this_thread::get_id();

              // change working directory to the python file
              stringstream s;
              s << "import os" << endl
                << "os.chdir(os.path.dirname(os.path.abspath('" << init_file_ << "')))" << endl
                << "del os" << endl;

              // set sys.argv
              s << "import sys" << endl;
              s << "sys.argv = [";
              if(h_argc > 1)
                for(auto i : Range(2, h_argc))
                  {
                    if(i>0) s << ',' << endl;
                    s << '"' << h_argv[i] << '"';
                  }
              s << ']' << endl;
              pyenv.exec(s.str());

              pyenv.exec_file(init_file_.c_str());
              Ng_SetRunning (0); 
            }
            catch (py::error_already_set const &) {
              PyErr_Print();
            }
            cout << IM(3) << "Finished executing " << init_file_  << endl;
            
            pythread_id = mainthread_id;
            
          }, filename).detach();
// 	    AcquireGIL gil_lock;
// 	    string command = string("import ") + filename;
// 	    pyenv.exec (command.c_str());
	    // pyenv.exec("from module1 import *");
	  }
#else
	  cout << "want to load python file, but python not enabled" << endl;	  
#endif
	  
	  return NG_TCL_OK;
	}
      catch (Exception & e)
	{
	  cerr << "\n\nCaught Exception in NGS_LoadPy:\n"
	       << e.What() << endl;
	  
	  ostringstream ost;
	  ost << "Exception in NGS_LoadPDE: \n " << e.What() << endl;
	  Ng_Tcl_SetResult (interp, (char*)ost.str().c_str(), NG_TCL_VOLATILE);
	  return NG_TCL_ERROR;
	}

    }
  Ng_Tcl_SetResult (interp, (char*)"no filename", NG_TCL_STATIC);
  return NG_TCL_ERROR;
}


int NGS_EnterCommand (ClientData clientData,
                      Tcl_Interp * interp,
                      int argc, const char *argv[])
{
  cout << "Enter command: ";
  string st;
  char ch;
  do
    {
      cin.get(ch);
      st += ch;
    }
  while (ch != '\n');
  cout << "command = " << st << endl;

  return NG_TCL_OK;
}




int NGS_SaveSolution (ClientData clientData,
		      Tcl_Interp * interp,
		      int argc, const char *argv[])
{
  Ng_Tcl_SetResult (interp, (char*)"Cannot save solution", NG_TCL_STATIC);
  return NG_TCL_ERROR;
}


int NGS_LoadSolution (ClientData clientData,
		      Tcl_Interp * interp,
		      int argc, const char *argv[])
{
  Ng_Tcl_SetResult (interp, (char*)"Cannot load solution", NG_TCL_STATIC);
  return NG_TCL_ERROR;
}




int NGS_PythonShell (ClientData clientData,
                    Tcl_Interp * interp,
                    int argc, const char *argv[])
{
  if(!netgen::netgen_executable_started) {
      Ng_Tcl_SetResult (interp, (char*)"This feature is not available when running from Python", NG_TCL_STATIC);
      return NG_TCL_ERROR;
  }
#ifdef NGS_PYTHON
  SpawnPython();
  return NG_TCL_OK;
#else
  cerr << "Sorry, you have to compile ngsolve with Python" << endl;
  return NG_TCL_ERROR;
#endif
}



int NGS_PrintMemoryUsage (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, const char *argv[])
{
  // netgen::BaseMoveableMem::Print ();
  
  // netgen::BaseDynamicMem::Print ();

  return NG_TCL_OK;
}



int NGS_PrintTiming (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, const char *argv[])
{
  ngstd::NgProfiler::Print (stdout);
  return NG_TCL_OK;
}




int NGS_GetData (ClientData clientData,
		 Tcl_Interp * interp,
		 int argc, const char *argv[])
{
  constexpr int BS = 1000;
  static char buf[BS];
  buf[0] = 0;
  stringstream str;

  snprintf (buf, BS, "0");
  str << buf;
  Ng_Tcl_SetResult (interp, (char*)str.str().c_str(), NG_TCL_VOLATILE);
  return NG_TCL_OK;
}


#ifdef NOT_WORKING
// from multibody ?
char playanimfile[256];
extern void PlayAnimFile(const char* name, int speed, int maxcnt);
/*
void * PlayAnim(void *)
{
  PlayAnimFile(playanimfile);

  running = 0;
  return NULL;
}
*/


int NGS_PlayAnim (ClientData clientData,
		  Tcl_Interp * interp,
		  int argc, const char *argv[])
{
  /*
  if (running)
    {
      Ng_Tcl_SetResult (interp, "Thread already running", NG_TCL_STATIC);
      return NG_TCL_ERROR;
    }
  running = 1;
  */
  int speed = 1;
  int maxcnt = 10;
  const char* name = "geom";
  if (argc >= 4)
    {
      speed = atoi(argv[1]);
      maxcnt = atoi(argv[2]);
      name = argv[3];
    }
  PlayAnimFile(name, speed, maxcnt);
  //  RunParallel (PlayAnim, NULL);
  return NG_TCL_OK;
}
#endif





int NGS_Set (ClientData clientData,
	     Tcl_Interp * interp,
	     int argc, const char *argv[])
{
  if (argc >= 3 && strcmp (argv[1], "time") == 0)
    {
      double time = double (atof (argv[2])) * 1e-6;
      cout << "NGS time = " << time << endl;
    }
  return NG_TCL_OK;
}



extern "C" int NGSolve_Init (Tcl_Interp * interp);
extern "C" void NGSolve_Exit ();

// tcl package dynamic load
extern "C" int NGS_DLL_HEADER Ngsolve_Init (Tcl_Interp * interp)
{
  return NGSolve_Init(interp);
}


// tcl package dynamic load
extern "C" int NGS_DLL_HEADER Ngsolve_Unload (Tcl_Interp * interp)
{
  return NG_TCL_OK;
}



namespace ngsolve { 
  // namespace bvp_cpp { extern int link_it; }
}

/*
namespace ngfem {
  namespace bdbequations_cpp { extern int link_it; }
  // extern int link_it_h1hofefo;
}
*/
namespace ngcomp {
  extern int link_it_hdivhofes;
}



int NGSolve_Init (Tcl_Interp * interp)
{
  cout << "NGSolve-" << ngsolve_version << endl;

#ifdef LAPACK
  cout << "Using Lapack" << endl;
#else
  cout << "sorry, no lapack" << endl; 
#endif


if(is_pardiso_available)
  cout << "Including sparse direct solver Pardiso" << endl;

#ifdef USE_UMFPACK
  cout << "Including sparse direct solver UMFPACK" << endl;
#endif

#ifdef USE_SUPERLU
  cout << "Including sparse direct solver SuperLU by Lawrence Berkeley National Laboratory" << endl;
#endif

  if (getenv ("NGSPROFILE"))
    NgProfiler::SetFileName (string("ngs.prof"));
  
#ifdef PARALLEL
  // in MyMPI
  // if (MyMPI_GetNTasks(MPI_COMM_WORLD) > 1)
  // TaskManager::SetNumThreads (1);
#endif
  cout << "Running parallel using " << TaskManager::GetMaxThreads() << " thread(s)" << endl;



  
#ifdef SOCKETS
  if(netgen::serversocketmanager.Good())
    serverjobmanager.Reset(new ServerJobManager(netgen::serversocketmanager,pde,ma,interp)); //!

#endif


  

#ifdef NGS_PYTHON

  // Only initialize a python environment if Netgen was started using the executable
  // In case we start the GUI from python, no further initialization is needed here
  if(netgen::netgen_executable_started)
  {
    Py_Initialize();

#if PY_VERSION_HEX < 0x03090000
    PyEval_InitThreads();
#endif

    py::module::import("__main__");

    // only relevant for gui+mpi -> undefined
    // if (MyMPI_GetId(MPI_COMM_WORLD) == 0)
      {
        pyenv.exec("from ngsolve import *");
        // Release GIL on this thread and reset thread state
        // to enable python operations on other threads
        PyEval_SaveThread();

        // Dummy SpawnPython (to avoid nasty Python GIL error)
        std::thread([]()
                    {
                      AcquireGIL gil_lock;
                      try{
                      pyenv.exec("from ngsolve import *");
                      pyenv.exec("from netgen import *");
                      }
                      catch (py::error_already_set const &) {
                        PyErr_Print();
                      }
                    }).detach();
      }
  }
#endif


  Ng_Tcl_CreateCommand (interp, "NGS_PrintRegistered", NGS_PrintRegistered);
  Ng_Tcl_CreateCommand (interp, "NGS_Help", NGS_Help);
  Ng_Tcl_CreateCommand (interp, "NGS_LoadPy", NGS_LoadPy);
  Ng_Tcl_CreateCommand (interp, "NGS_EnterCommand", NGS_EnterCommand);
  Ng_Tcl_CreateCommand (interp, "NGS_PythonShell", NGS_PythonShell);
  Ng_Tcl_CreateCommand (interp, "NGS_PrintMemoryUsage", NGS_PrintMemoryUsage);
  Ng_Tcl_CreateCommand (interp, "NGS_PrintTiming", NGS_PrintTiming);

  Ng_Tcl_CreateCommand (interp, "NGS_GetData", NGS_GetData);
  Ng_Tcl_CreateCommand (interp, "NGS_Set", NGS_Set);

  /*
  const Array<NumProcs::NumProcInfo*> & npi = GetNumProcs().GetNumProcs();

  Array<int> sort(npi.Size());
  Array<string> names(npi.Size());
  for (int i = 0; i < npi.Size(); i++)
    {
      sort[i] = i;
      names[i] = npi[i]->name;
    }
  BubbleSort (names, sort);

  Tcl_Eval (interp, "menu .ngmenusolvehelpnp");
  for (int ii = 0; ii < npi.Size(); ii++)
    {
      int i = sort[ii];
      string command =
	".ngmenusolvehelpnp add command -label \"numproc " + npi[i]->name + "\" \
	-command { tk_messageBox -title \"Help\" -message  [ NGS_Help  numproc " + npi[i]->name + "] -type ok }" ;

      Tcl_Eval (interp, (char *) command.c_str());
    }
  */

  // trick for forcing linking static libs
  // ngsolve::bvp_cpp::link_it = 0;
  // ngfem::bdbequations_cpp::link_it = 0;
  // ngfem::link_it_h1hofefo = 0;



  return NG_TCL_OK;
}


/*
void NGSolve_Exit ()
{
  cout << "NGSolve says good bye" << endl;
  // delete pde;
  // delete ma;
  // ma = 0;

  // Parallel_Exit();
}
*/






#ifdef PARALLEL
				    

extern "C" void NGS_ParallelRun (const string & message);

#ifdef NGS_PYTHON
void Parallel_InitPython ()
{
  static bool python_initialized = false;
  if (!python_initialized)
    {
      cout << "ini (parallel) python" << endl;
      Py_Initialize();
#if PY_VERSION_HEX < 0x03090000
      PyEval_InitThreads();
#endif

      
      py::module m = py::module::import("__main__");
      // pyenv = PythonEnvironment (m);
      {
	m.def ("Redraw", 
	       []() {Ng_Redraw();});
      }
      
      // cout << "ini python complete" << endl;	  

      pyenv.exec("from ngsolve import *");
      //PyEval_ReleaseLock();
      PyEval_SaveThread();

      python_initialized = true;
    }
}
#endif




void NGS_ParallelRun (const string & message)
{
  TaskManager::SetNumThreads (1);

  NGSOStream::SetGlobalActive (false);

  // make sure to have lapack loaded (important for newer MKL !!!)
  Matrix<> a(100), b(100), c(100);
  a = 1; b = 2; c = a*b | Lapack;


  if (getenv ("NGSPROFILE"))
    {
      stringstream filename;
      filename << "ngs.prof." << NgMPI_Comm(NG_MPI_COMM_WORLD).Rank();
      NgProfiler::SetFileName (filename.str());
    }

  if ( message == "ngs_loadngs" )
    {
      // MPI_Comm_dup ( MPI_COMM_WORLD, &ngs_comm);
    }


#ifdef NGS_PYTHON
  else if (message.substr(0,7) == "ngs_py " ) 
    {

      string command = message.substr(7);
      std::thread( [](string a_command){
	  Parallel_InitPython ();
	  AcquireGIL gil_lock;
	  pythread_id = std::this_thread::get_id();
	  // PythonEnvironment & py_env = PythonEnvironment::getInstance();
	  pyenv.exec(a_command);
	  pythread_id = mainthread_id;
	  
	}, command).detach();
    }
#endif
  return;
}


void Parallel_Exit ()
{
  ;
}


#endif





#ifdef WIN32
void * __cdecl my_operator_new_replacement(size_t _count)
{
	void * p = operator new(_count);
//	cout << "alloc, cnt = " << _count << ", ptr = " << p << endl;
  return p;
  }

void __cdecl my_operator_delete_replacement(void * _ptr)
{
//  cout << "delete, ptr = " << _ptr << endl;
  delete (_ptr);
}

void * __cdecl my_operator_new_array_replacement(size_t _count)
{
		void * p = operator new[] (_count);
//	cout << "alloc [], cnt = " << _count << ", ptr = " << p << endl;
  return p;

//   return operator new[](_count);
}

void __cdecl my_operator_delete_array_replacement(void * _ptr)
{
//  cout << "delete[], ptr = " << _ptr << endl;
  delete [] (_ptr);
}
#endif
