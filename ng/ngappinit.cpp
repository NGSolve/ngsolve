/*
  The main function of netgen.
  This file is a modification of tkAppInit.c from the Tcl/Tk package
*/

#include <mystdlib.h> 
#include <inctcl.hpp>
#include <meshing.hpp>

#ifdef NGPYTHON
#include <boost/python.hpp>
#endif

/*
namespace netgen
{
  int id = 0, ntasks = 1;
}
*/

#ifdef PARALLEL
#include <mpi.h>

extern void ParallelRun();
namespace netgen
{
  MPI_Comm mesh_comm;
}
#endif


#include "../libsrc/interface/writeuser.hpp"

namespace netgen
{
  Flags parameters;
}
 


using netgen::parameters;
using netgen::ngdir;
using netgen::verbose;
using netgen::Array;
using netgen::RegisterUserFormats;




/*
 * The following variable is a special hack that is needed in order for
 * Sun shared libraries to be used for Tcl.
 */

// extern "C" int matherr();
// int *tclDummyMathPtr = (int *) matherr;

extern "C" int Ng_ServerSocketManagerInit (int port);
extern "C" int Ng_ServerSocketManagerRun (void);
 
bool nodisplay = false;
bool shellmode = false;


/*
 *
 *     The Netgen main function
 *
 */

int main(int argc, char ** argv)
{

#ifdef PARALLEL
  int mpi_required = MPI_THREAD_MULTIPLE;
#ifdef VTRACE
  mpi_required = MPI_THREAD_SINGLE;
#endif
  int mpi_provided;
  MPI_Init_thread(&argc, &argv, mpi_required, &mpi_provided);          

  MPI_Comm_size(MPI_COMM_WORLD, &netgen::ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &netgen::id);
  

  MPI_Comm_dup ( MPI_COMM_WORLD, &netgen::mesh_comm);
#endif

  if ( netgen::id == 0 )
    {
      cout << "NETGEN-" << PACKAGE_VERSION << endl;
      
      cout << "Developed by Joachim Schoeberl at" << endl
	   << "2010-xxxx Vienna University of Technology" << endl
	   << "2006-2010 RWTH Aachen University" << endl
           << "1996-2006 Johannes Kepler University Linz" << endl;
      
#ifdef OCCGEOMETRY
      cout << "Including OpenCascade geometry kernel" << endl;
#endif
      
#ifdef ACIS
      cout << "Including ACIS geometry kernel" << endl;
#endif

#ifdef LINUX
      //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
      //feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
      //cout << "Handle Floating Point Exceptions: " << fegetexcept() << endl; 
#endif  

#ifdef DEBUG
      cout << "You are running the debug version !" << endl;
#endif


#ifdef PARALLEL
      if (netgen::ntasks == 1)
	{
	  cout << "Run parallel Netgen with 'mpirun -np xy netgen'" << endl;
	}
      else
	{
	  cout << "Running MPI - parallel using " 
	       << netgen::ntasks << " processor" 
	       << ((netgen::ntasks > 1) ? "s " : " ") << endl;
	  
	  cout << "MPI-version = " << MPI_VERSION << '.' << MPI_SUBVERSION << endl;

	  if (mpi_provided == MPI_THREAD_MULTIPLE)
	    cout << "multithreaded MPI is supported" << endl;
	}
#endif
    }


  netgen::h_argc = argc;
  netgen::h_argv = argv;

  // command line arguments:
  for (int i = 1; i < argc; i++)
    {
      if (argv[i][0] == '-')
	parameters.SetCommandLineFlag (argv[i]);
      else
        {
          if (strstr(argv[i], ".py"))
            parameters.SetFlag ("py", argv[i]);
          else
            parameters.SetFlag ("geofile", argv[i]);
        }
    }


  if (getenv ("NETGENDIR") && strlen (getenv ("NETGENDIR")))
    ngdir = getenv ("NETGENDIR");
  else
    ngdir = ".";
  
  verbose = parameters.GetDefineFlag ("V");

  if (verbose)
    cout << "NETGENDIR = " << ngdir << endl;
  

  if ( netgen::id == 0 )
    {
      if (parameters.StringFlagDefined ("testout"))      
        netgen::testout = new ofstream (parameters.GetStringFlag ("testout", "test.out"));


#ifdef SOCKETS
      Ng_ServerSocketManagerInit(static_cast<int>(parameters.GetNumFlag("serversocket",-1)));
      if(parameters.GetNumFlag("serversocket",-1) > 0 && !parameters.GetDefineFlag("display")) 
        nodisplay = true;
#endif
  
      if(parameters.GetDefineFlag("batchmode"))
        nodisplay = true;

    
      if(parameters.GetDefineFlag("shellmode"))
        {
          nodisplay = true;
          shellmode = true;
        }

      Tcl_FindExecutable(NULL);

      // initialize application
      Tcl_Interp * myinterp = Tcl_CreateInterp ();
      if (Tcl_AppInit (myinterp) == TCL_ERROR)
        {
          cerr << "Exit Netgen due to initialization problem" << endl;
          exit (1);
        }



      // parse tcl-script
      int errcode;

      bool internaltcl = false;
      if (shellmode)
        internaltcl = false;
  
      if (verbose)
        {
          cout << "Tcl header version = " << TCL_PATCH_LEVEL << endl;
          Tcl_Eval (myinterp, "puts \"Tcl runtime version = [info patchlevel] \";");
        }

      if (parameters.GetDefineFlag ("internaltcl"))
        internaltcl=true;
      if (parameters.GetDefineFlag ("externaltcl"))
        internaltcl=false;

 

      if (internaltcl)
        {
          if (verbose)
            cout << "using internal Tcl-script" << endl;
      
          // connect to one string 
          extern const char * ngscript[];
          const char ** hcp = ngscript;
          int len = 0;
          while (*hcp)
            len += strlen (*hcp++); 

          char * tr1 = new char[len+1];
          *tr1 = 0;
          hcp = ngscript;
      
          char * tt1 = tr1;
          while (*hcp)
            {
              strcat (tt1, *hcp); 
              tt1 += strlen (*hcp++);
            }
      
          errcode = Tcl_Eval (myinterp, tr1);
          delete [] tr1;
        }

      else

        {
          string startfile = ngdir + "/ng.tcl";
      
          if (verbose)
            cout << "Load Tcl-script from " << startfile << endl;
      
          errcode = Tcl_EvalFile (myinterp, (char*)startfile.c_str());
        }

      if (errcode)
        {
          cout << "Error in Tcl-Script:" << endl;
          // cout << "result = " << myinterp->result << endl;
          cout << "result = " << Tcl_GetStringResult (myinterp) << endl;
          // cout << "in line " << myinterp->errorLine << endl;

          cout << "\nMake sure to set environment variable NETGENDIR to directory containing ng.tcl" << endl;
          exit (1);
        }


      // lookup user file formats and insert into format list:
      Array<const char*> userformats;
      Array<const char*> extensions;
      RegisterUserFormats (userformats, extensions);

      ostringstream fstr;
      
      tcl_const char * exportft = Tcl_GetVar (myinterp, "exportfiletype", 0);
      for (int i = 1; i <= userformats.Size(); i++)
	{
	  fstr << ".ngmenu.file.filetype add radio -label \"" 
	       << userformats.Get(i) << "\" -variable exportfiletype\n";
	  fstr << "lappend meshexportformats { {" << userformats.Get(i) << "} {" << extensions.Get(i) << "} }\n";
	}

      Tcl_Eval (myinterp, (char*)fstr.str().c_str());
      Tcl_SetVar (myinterp, "exportfiletype", exportft, 0);


#ifdef SOCKETS
      Ng_ServerSocketManagerRun();
#endif

      // start event-loop
      Tk_MainLoop();
      Tcl_DeleteInterp (myinterp); 
      Tcl_Exit(0);
    }

#ifdef PARALLEL
  else
    {
      ParallelRun();
      MPI_Finalize();
    }  

#endif
  
  return 0;		
}



/*
extern "C" int Tix_Init (Tcl_Interp * interp);
extern "C" int Itcl_Init (Tcl_Interp * interp);
extern "C" int Itk_Init (Tcl_Interp * interp);
*/
extern "C" int Ng_Init (Tcl_Interp * interp);
extern "C" int Ng_Vis_Init (Tcl_Interp * interp);



// extern Tcl_PackageInitProc * Tk_SafeInit;

/*
 *
 * Initialize packages
 *
 */

// extern "C" int NGSolve_Init (Tcl_Interp * interp);


int Tcl_AppInit(Tcl_Interp * interp)
{

  if (Tcl_Init(interp) == TCL_ERROR) { 
    cerr << "Problem in Tcl_Init: " << endl;
    cout << "result = " << Tcl_GetStringResult (interp) << endl;
    // cerr << interp->result << endl;
    // return TCL_ERROR;
  }
  
  if (!nodisplay && Tk_Init(interp) == TCL_ERROR) {
    cerr << "Problem in Tk_Init: " << endl;
    cout << "result = " << Tcl_GetStringResult (interp) << endl;
    // cerr << interp->result << endl;
    // return TCL_ERROR;
  }


  if (Ng_Init(interp) == TCL_ERROR) {
    cerr << "Problem in Ng_Init: " << endl;
    cout << "result = " << Tcl_GetStringResult (interp) << endl;
    // cerr << interp->result << endl;
    // return TCL_ERROR;
  }
 
  if (!nodisplay && Ng_Vis_Init(interp) == TCL_ERROR) {
    cerr << "Problem in Ng_Vis_Init: " << endl;
    cout << "result = " << Tcl_GetStringResult (interp) << endl;
    // cerr << interp->result << endl;
    // return TCL_ERROR;
  }






#ifdef TRAFO
  //   extern int Trafo_Init (Tcl_Interp * interp);
  //   if (Trafo_Init(interp) == TCL_ERROR) 
  //     {
  //       cerr << "Problem in Trafo_Init: " << endl;
  //       cerr << interp->result << endl;
  //       return TCL_ERROR;
  //     }
#endif

#ifdef EBGELAST
  extern int EBGElast_Init (Tcl_Interp * interp);
  if(EBGElast_Init(interp) == TCL_ERROR)
    {
      cerr << "Problem in EBGElast_Init: " << endl;
      cerr << interp->result << endl;
      return TCL_ERROR;
    }

#endif

#ifdef SMALLTRAFO
  extern int SmallModels_Init (Tcl_Interp * interp);
  if(SmallModels_Init(interp) == TCL_ERROR)
    {
      cerr << "Problem in SmallModel_Init: " << endl;
      cerr << interp->result << endl;
      return TCL_ERROR;
    }

#endif

#ifdef SOCKETS
  extern int Ng_Socket_Init (Tcl_Interp * interp);
  if ( Ng_Socket_Init(interp) == TCL_ERROR)
    {
      cerr << "Problem in Ng_Socket_Init: " << endl;
      cerr << interp->result << endl;
      return TCL_ERROR;
    }

#endif


#ifdef ZUGSTANGE
  extern int Zugstange_Init (Tcl_Interp * interp);
  if (Zugstange_Init(interp) == TCL_ERROR)
    {
      cerr << "Problem in Zugstange_Init: " << endl;
      cerr << interp->result << endl;
      return TCL_ERROR;
    }
#endif


  Tcl_StaticPackage(interp, "Tk", Tk_Init, 0);
  return TCL_OK;
}





// link MKL with netgen
// necessary for MKL 11.x, since MKL complains if started
// from the ngsolve shared library

#ifdef LINKMKL
extern "C" double ddot_(int *n, double *dx, int *incx, double *dy, 
                        int *incy);
void mkldummy()
{
  int n = 1, one = 1;
  double a = 1, b = 1; 
  ddot_(&n, &a, &one, &b, &one);
}
#endif
