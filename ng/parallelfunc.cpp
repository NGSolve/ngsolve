#ifdef PARALLEL

#ifdef OCCGEOMETRY
#include <occgeom.hpp>
#endif

// #include <mystdlib.h>

#include <meshing.hpp>

// #include "incvis.hpp"
#include <visual.hpp>


#ifdef PARALLELGL

#ifndef WIN32

#define GLX_GLXEXT_LEGACY
#define GLX_GLXEXT_PROTOTYPES

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>  /* for XA_RGB_DEFAULT_MAP atom */
#include <GL/glx.h>  
#include <GL/glxext.h>
#endif

#endif



#include <geometry2d.hpp>
#include <stlgeom.hpp>

// #include <incvis.hpp>
// #include <visual.hpp>
// #include <mystdlib.h>
// #include <myadt.hpp>
// #include <linalg.hpp>
// #include <csg.hpp>



#include <meshing.hpp>

// #include "parallel.hpp"
#include "parallelfunc.hpp"










namespace netgen
{
#include "../interface/writeuser.hpp"
  extern string ngdir;

#ifdef OPENGL
  extern VisualSceneMesh vsmesh;
#endif
}

void Parallel_Exit();


namespace netgen {
  extern AutoPtr<Mesh>  mesh;

  // geometry: either CSG, or, if an other is non-null, 
  // then the other
  extern AutoPtr<CSGeometry> geometry ;

  extern STLGeometry * stlgeometry;
  extern AutoPtr<SplineGeometry2d> geometry2d ;

#ifdef OCCGEOMETRY
  extern OCCGeometry * occgeometry;
#endif

#ifdef ACIS
  extern ACISGeometry * acisgeometry;
#endif

}

using namespace netgen;
using netgen::RegisterUserFormats;

#ifdef PARALLEL
void Ng_Exit ()
{
#ifdef NGSOLVE
  Parallel_Exit();
#endif

  delete stlgeometry;
  stlgeometry = NULL;
   

#ifdef OCCGEOMETRY 
  delete occgeometry; 
  occgeometry = 0;
#endif


  geometry.Reset (NULL);
  geometry2d.Reset (NULL);
    
  //    delete testout;
  return;
}
#endif



void ParallelRun()
{   
  string message;
      
  MPI_Status status;
      
  //       int id, rc, ntasks;

  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);



  bool test = true;



  testout = new ofstream (string("testout_proc") + id  );      

  while ( test )
    {
#ifdef SCALASCA
#pragma pomp inst begin (message)
#endif

      (*testout) << "wait for mess " << endl;
      MyMPI_Recv ( message, 0 );
      (*testout) << "message " << message << endl;

#ifdef SCALASCA
#pragma pomp inst end (message)
#endif

#ifdef NGSOLVE

      if ( message.compare(0, 3, "ngs") == 0 ) 
        {
          PrintMessage ( 2, "Starting NgSolve routine ", message ) ;
          NGS_ParallelRun (message);
        }
      else
#endif
	
        if ( message == "mesh" )
          {
  VT_USER_START ("Mesh::ReceiveParallelMesh");
            mesh.Reset( new netgen::Mesh);
            mesh->ReceiveParallelMesh();
  VT_USER_END ("Mesh::ReceiveParallelMesh");
          }

        else if ( message == "overlap++" )
          {
	    PrintMessage (1, "overlap++++++");
            mesh -> UpdateOverlap();
          }

        else if ( message == "visualize" )
          {
            cout << "p" << id << ": ACHTUNG - alles wieder zumachen, sonst geht nix mehr :)" << endl;
            cout << "Tcl-disabled" << endl;
            /*	      

            // initialize application
            Tcl_Interp * myinterp = Tcl_CreateInterp ();
            if (Tcl_AppInit (myinterp) == TCL_ERROR)
              {
                cerr << "Exit Netgen due to initialization problem" << endl;
                exit (1);
              }

            string startfile = ngdir + "/libsrc/parallel/ng_parallel.tcl";
	      
            if (verbose)
              cout << "Load Tcl-script from " << startfile << endl;
	      
            int errcode = Tcl_EvalFile (myinterp, (char*)startfile.c_str());

            if (errcode)
              {
                cout << "Error in Tcl-Script:" << endl;
                cout << "result = " << myinterp->result << endl;
                cout << "in line " << myinterp->errorLine << endl;
		  
                if (myinterp->errorLine == 1)
                  cout << "\nMake sure to set environment variable NETGENDIR" << endl;
		  
                exit (1);
              }

            // lookup user file formats and insert into format list:
            ARRAY<const char*> userformats;
            RegisterUserFormats (userformats);
	      
            ostringstream fstr;
            for (int i = 1; i <= userformats.Size(); i++)
              fstr << ".ngmenu.file.filetype add radio -label \"" 
                   << userformats.Get(i) << "\" -variable exportfiletype\n";
	      
	      
            Tcl_Eval (myinterp, (char*)fstr.str().c_str());
            Tcl_SetVar (myinterp, "exportfiletype", "Neutral Format", 0);



            Tk_MainLoop();
	      
	      
            Tcl_DeleteInterp (myinterp); 
*/	      

          }


#ifdef PARALLELGL

        else if ( message == "redraw" )
          {
            // draw into the same GLX - drawing context
	    // works on parallel machine, but 
            // did not manage to get glXImportContextEXT working on Laptop (JS)

	    string redraw_cmd;
	    MyMPI_Recv (redraw_cmd, 0);

	    // PrintMessage (1, "Redraw - ", redraw_cmd);
                  
	    static string displname;
	    static int curDrawable, contextid;
	    
	    static Display * display = NULL;
	    static GLXContext context;
	    static XVisualInfo * visinfo = 0;
	    
	    // if (!display)
	    if (redraw_cmd == "init")
	      {
		MyMPI_Recv (displname, 0);
		MyMPI_Recv (curDrawable, 0);
		MyMPI_Recv (contextid, 0);
		
		
		display = XOpenDisplay (displname.c_str());
		
		/*
		  PrintMessage (3, "displ - name = ", displname);
		  PrintMessage (3, "display = ", display,
		  " display props: ",
		  " screen w = ", XDisplayWidth (display, 0),
		  " , h = ", XDisplayHeight (display, 0));
		*/

		Window win;
		int wx, wy;
		unsigned int ww, wh, bw, depth;
		cout << "got drawable: " << curDrawable << endl;
		
		XGetGeometry(display, curDrawable, &win,
			     &wx, &wy, &ww, &wh,
			     &bw, &depth);
		
                cout << "P" << id << ": window-props:  x = " << wx << ", y = " << wy 
                     << ", w = " << ww << ", h = " << wh << ", depth = " << depth << endl;
		
#define VISUAL
#ifdef VISUAL
		
                // make a new GLXContext
                // first, generate a visual (copied from togl)
		
                int     attrib_list[1000];
		
#  define MAX_ATTEMPTS 12
                static int ci_depths[MAX_ATTEMPTS] = {
                  8, 4, 2, 1, 12, 16, 8, 4, 2, 1, 12, 16
                };
                static int dbl_flags[MAX_ATTEMPTS] = {
                  0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1
                };
		
                /* It may take a few tries to get a visual */
		
		for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++) 
		  {
		    int attrib_count = 0;
		    attrib_list[attrib_count++] = GLX_USE_GL;
		      
                      
		    /* RGB[A] mode */
		    attrib_list[attrib_count++] = GLX_RGBA;
		    attrib_list[attrib_count++] = GLX_RED_SIZE;
		    attrib_list[attrib_count++] = 1;
		    attrib_list[attrib_count++] = GLX_GREEN_SIZE;
		    attrib_list[attrib_count++] = 1;
		    attrib_list[attrib_count++] = GLX_BLUE_SIZE;
		    attrib_list[attrib_count++] = 1;
		    // attrib_list[attrib_count++] = GLX_ALPHA_SIZE;
		    // attrib_list[attrib_count++] = 1;
                      
		    attrib_list[attrib_count++] = GLX_DEPTH_SIZE;
		    attrib_list[attrib_count++] = 1;

		    attrib_list[attrib_count++] = GLX_DOUBLEBUFFER;

		    attrib_list[attrib_count++] = None;
                    
		    visinfo = glXChooseVisual(display, 0, 
					      attrib_list);
		    if (visinfo) {
		      /* found a GLX visual! */
		      // cout << "found VISINFO !!!" << endl;

		      /*
		      int hi = 0;
		      std::cout << "attribs = ";
		      while (attrib_list[hi] != None)
			std::cout << attrib_list[hi++] << " ";
		      std::cout << std::endl;
		      */

		      break;
		    }
		  }

		if (!visinfo)
		  cerr << "no VISINFO found" << endl;

                // context = glXCreateContext( display, visinfo, 0, /* curContext, */ False );
		context = glXCreateContext( display, visinfo, glXImportContextEXT ( display, contextid ), False );
                cout << "context = " << context << endl;

                glXMakeCurrent (display, curDrawable, context);


#else
                // try to get GLXcontext from the master. 
                // this needs an indirect context (BUT DOES NOT WORK ????)
		
		context = glXImportContextEXT ( display, contextid );


		PrintMessage (1, "GLX-contextid = " , contextid,
			      " imported context ", context);


		// glXMakeCurrent (display, curDrawable, context);
#endif

		// PrintMessage (1, "redraw - init complete");
		int hi = id;
		MyMPI_Send (hi, 0);
	      }
	    
	    if (redraw_cmd == "broadcast")
	      {
		vsmesh.Broadcast ();
	      }

	    if (redraw_cmd == "linelist")
	      {
		glXMakeCurrent (display, curDrawable, context);
		vsmesh.BuildLineList();
		glXMakeCurrent (display, None, NULL);
	      }

	    if (redraw_cmd == "filledlist")
	      {
		glXMakeCurrent (display, curDrawable, context);
		vsmesh.BuildFilledList();
		glXMakeCurrent (display, None, NULL);
	      }


	    if (redraw_cmd == "solsurfellist")
	      {
		glXMakeCurrent (display, curDrawable, context);
		vssolution.DrawSurfaceElements();
		glXMakeCurrent (display, None, NULL);
	      }

	    
	  }
#endif






	else if ( message ==  "end" )
	      {
		PrintMessage (1, "EXIT");
		test = false;
		// end netgen
		Ng_Exit();
	      }


	    else
	      {
		PrintMessage ( 1, "received unidentified message '" + message + "'\n");
	     
		test = false;
	      }
	  
	  }
      
      return;
    }




#endif
