/*
  The main function of netgen.
  This file is a modification of tkAppInit.c from the Tcl/Tk package
*/

#include <mystdlib.h>
#include "incvis.hpp"
#include <meshing.hpp>

#ifdef LINUX
// #include <fenv.h>
#endif

/*
#ifndef WIN32
#include <dlfcn.h>
#endif
*/

#ifdef PARALLEL
#include <mpi.h>

namespace netgen
{
  int id, ntasks;
  MPI_Group MPI_HIGHORDER_WORLD;
  MPI_Comm MPI_HIGHORDER_COMM;
}
#endif


#include "parallelfunc.hpp"


namespace netgen
{
#include "writeuser.hpp"
  extern string ngdir;
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
  // parallel profiling
#pragma pomp inst init

  MPI_Init(&argc, &argv);          

  MPI_Comm_size(MPI_COMM_WORLD, &netgen::ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &netgen::id);
  
  MPI_Group MPI_GROUP_WORLD;

  int n_ho = netgen::ntasks - 1;
  int * process_ranks = new int[netgen::ntasks-1];
  for ( int i = 0; i < netgen::ntasks-1; i++ )
    process_ranks[i] = i+1;

  MPI_Comm_group ( MPI_COMM_WORLD, &MPI_GROUP_WORLD);
  MPI_Group_incl ( MPI_GROUP_WORLD, n_ho, process_ranks, & netgen::MPI_HIGHORDER_WORLD);
  MPI_Comm_create ( MPI_COMM_WORLD, netgen::MPI_HIGHORDER_WORLD, & netgen::MPI_HIGHORDER_COMM);

#pragma pomp inst begin(main)
#endif

  if ( netgen::id == 0 )
    {
      cout << "NETGEN-" << PACKAGE_VERSION << endl;
      
      cout << "Developed at RWTH Aachen University, Germany" << endl
           << "and Johannes Kepler University Linz, Austria" << endl;

      
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

#ifdef USE_SUPERLU
      cout << "Including sparse direct solver SuperLU by Lawrence Berkeley National Laboratory" << endl;
#endif


#ifdef PARALLEL
      cout << "Including MPI " << endl;
      cout << "Using " <<  netgen::ntasks << " processor" 
           << ((netgen::ntasks > 1) ? "s " : " ") << endl;
#endif
    }
  else
    {
      ;// nodisplay = true;
    }


  // command line arguments:
  for (int i = 1; i < argc; i++)
    {
      if (argv[i][0] == '-')
	parameters.SetCommandLineFlag (argv[i]);
      else
	parameters.SetFlag ("geofile", argv[i]);
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
  
#ifdef PARALLEL
      internaltcl = false;
#endif

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
      // Tcl_SetVar (myinterp, "exportfiletype", "Neutral Format", 0);
      Tcl_SetVar (myinterp, "exportfiletype", exportft, 0);


      // For adding an application, parse the file here,
      // and call the init-procedure below
      // #define DEMOAPP
#ifdef DEMOAPP  
      Tcl_EvalFile (myinterp, "demoapp/demoapp.tcl");
#endif

#ifdef ADDON
      Tcl_EvalFile (myinterp, "addon/addon.tcl");
#endif

#ifdef SOCKETS
      Ng_ServerSocketManagerRun();
#endif

      // start event-loop
      Tk_MainLoop();
  
      Tcl_DeleteInterp (myinterp); 

#ifdef PARALLEL
#pragma pomp inst altend(main)

      // MPI beenden
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
#endif
      
      Tcl_Exit(0);
    }
#ifdef PARALLEL
  else
    {
      // main for parallel processors    
      ParallelRun();

#pragma pomp inst end(main)

      // MPI beenden
      MPI_Barrier(MPI_COMM_WORLD);
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

  // if ITcl and ITk are installed on the system, then
  // they must also be initialized
  /*
    if (Itcl_Init(interp) == TCL_ERROR) {
    cerr << "Problem in Itcl_Init: " << endl;
    cerr << interp->result << endl;
    // return TCL_ERROR;
    }

    if (Itk_Init(interp) == TCL_ERROR) {
    cerr << "Problem in Itk_Init: " << endl;
    cerr << interp->result << endl;
    // return TCL_ERROR;
    }
  */
    /*
  if (!nodisplay && Tix_Init(interp) == TCL_ERROR) {
    cerr << "Problem in Tix_Init: " << endl;
    cerr << interp->result << endl;
    // return TCL_ERROR;
  }
    */

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


  /*
  if (NGSolve_Init(interp) == TCL_ERROR) 
    return TCL_ERROR;
  */

#ifdef DEMOAPP
  extern int DemoApp_Init (Tcl_Interp * interp);
  if (DemoApp_Init(interp) == TCL_ERROR) 
    {
      return TCL_ERROR;
    }
#endif 
#ifdef ADDON
  extern int AddOn_Init (Tcl_Interp * interp);
  if (AddOn_Init(interp) == TCL_ERROR) 
    {
      return TCL_ERROR;
    }
#endif
#ifdef METIS_OLD
  extern int NgMetis_Init (Tcl_Interp * interp);
  if (NgMetis_Init(interp) == TCL_ERROR) 
    {
      return TCL_ERROR;
    }
#endif    

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
