/*

The interface between the GUI and the netgen library

*/


#include <mystdlib.h>
#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>

#include <geometry2d.hpp>
#include <stlgeom.hpp>
#include <meshing.hpp>


#ifdef OCCGEOMETRY
#include <occgeom.hpp>
#endif

#include "../libsrc/meshing/bcfunctions.hpp"

#include <incvis.hpp>
#include <visual.hpp>


#ifdef SOCKETS
#include "../libsrc/sockets/sockets.hpp"
#include "../libsrc/sockets/socketmanager.hpp"
#endif

// #include <parallel.hpp>

// to be sure to include the 'right' togl-version
#include "togl_1_7.h"
// #include "Togl2/togl.h"

extern bool nodisplay;


#include <nginterface.h>



#ifdef _MSC_VER
// Philippose - 30/01/2009
// MSVC Express Edition Support
#ifdef MSVC_EXPRESS

// #include <pthread.h>

  static pthread_t meshingthread;
  void RunParallel ( void * (*fun)(void *), void * in)
  {
    if (netgen::mparam.parthread)
     {
	     pthread_attr_t attr;
	     pthread_attr_init (&attr);
	     // the following call can be removed if not available:
	     pthread_attr_setstacksize(&attr, 1000000);
	     //pthread_create (&meshingthread, &attr, fun, NULL);
	     pthread_create (&meshingthread, &attr, fun, in);
     }
     else
     fun (in);
  }

#else // Using MS VC++ Standard / Enterprise / Professional edition

  // Afx - Threads need different return - value:

  static void* (*sfun)(void *);
  unsigned int fun2 (void * val)
  {
    sfun (val);
    return 0;
  }

  void RunParallel ( void* (*fun)(void *), void * in)
  {
    sfun = fun;
    if (netgen::mparam.parthread)
      AfxBeginThread (fun2, in);
    //AfxBeginThread (fun2, NULL);
    else
      fun (in);
  }

#endif // #ifdef MSVC_EXPRESS

#else  // For #ifdef _MSC_VER

// #include <pthread.h>

  static pthread_t meshingthread;
  void RunParallel ( void * (*fun)(void *), void * in)
  {
    if (netgen::mparam.parthread)
      {
	pthread_attr_t attr;
	pthread_attr_init (&attr);
	// the following call can be removed if not available:
	pthread_attr_setstacksize(&attr, 1000000);
	//pthread_create (&meshingthread, &attr, fun, NULL);
	pthread_create (&meshingthread, &attr, fun, in);
      }
    else
      fun (in);
  }

#endif // #ifdef _MSC_VER




namespace netgen
{
#include "../libsrc/interface/writeuser.hpp"
#include "demoview.hpp"
}


#ifdef ACIS
#include "ng_acis.hpp"
#endif


#ifdef JPEGLIB
#include <jpeglib.h>
#endif

#ifdef FFMPEG
extern "C" {
  /*
#include <ffmpeg/avcodec.h>
#include <ffmpeg/avformat.h>
#include <ffmpeg/swscale.h>
  */
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libswscale/swscale.h>
}
#endif

#ifdef NGSOLVE
extern "C" void NGSolve_Exit();
#endif

// extern void * ngsolve_handle;


namespace netgen
{

  NetgenOutStream operator<< ( ostream & ost, Imp  imp )
  {
    return ( NetgenOutStream ( &ost, imp ) );
  }

  NetgenOutStream operator<< ( ostream & ost, Proc proc )
  {
    return ( NetgenOutStream ( &ost, proc ) );
  }


  NetgenOutStream operator<< ( ostream & ost, Procs & procs )
  {
    return ( NetgenOutStream ( &ost, procs ) );
  }


  // global variable mesh (should not be used in libraries)
  AutoPtr<Mesh> mesh;

  // geometry: either CSG, or, if an other is non-null,
  // then the other
  AutoPtr<CSGeometry> geometry (new CSGeometry(""));
  STLGeometry * stlgeometry = NULL;
  AutoPtr<SplineGeometry2d> geometry2d (0);

#ifdef OCCGEOMETRY
  OCCGeometry * occgeometry = NULL;
#endif

  NetgenGeometry * ng_geometry = new NetgenGeometry;

  Tcl_Interp * tcl_interp;

#ifdef SOCKETS
  AutoPtr<ClientSocket> clientsocket;
  ServerSocketManager serversocketmanager;
  //Array< AutoPtr < ServerInfo > > servers;
  Array< ServerInfo* > servers;
  AutoPtr<ServerSocketUserNetgen> serversocketusernetgen;
#endif

  // OpenGL near and far clipping planes
  extern double pnear;
  extern double pfar;


  // visualization scenes, pointer vs selects which one is drawn:

  static VisualScene vscross;
  static VisualSceneGeometry vsgeom;
  static VisualSceneGeometry2d vsgeom2d;
#ifdef OPENGL
  static VisualSceneSTLGeometry vsstlgeom;
  extern VisualSceneSTLMeshing vsstlmeshing;
#endif

#ifdef OCCGEOMETRY
  static VisualSceneOCCGeometry vsoccgeom;
#endif



#ifdef OPENGL
  extern VisualSceneSurfaceMeshing vssurfacemeshing;
#endif
  extern VisualSceneMesh vsmesh;
  extern VisualSceneMeshDoctor vsmeshdoc;
  static VisualSceneSpecPoints vsspecpoints;

  VisualSceneSolution vssolution;


#ifdef STEP
  static VisualSceneSTEPGeometry vsstepgeom;
#endif



  VisualScene *vs = &vscross;



  static char * err_needsmesh = (char*) "This operation needs a mesh";
  static char * err_needsstlgeometry = (char*) "This operation needs an STL geometry";
  static char * err_jobrunning = (char*) "Meshing Job already running";





#ifndef SMALLLIB
  // Destination for messages, errors, ...
  DLL_HEADER void Ng_PrintDest(const char * s)
  {
#ifdef PARALLEL
    int id, ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
#else
    int id = 0; int ntasks = 1;
#endif


    if ( ntasks == 1 )
      (*mycout) << s << flush;
    else
      (*mycout) << "p" << id << ": " << s << flush ;
  }

  DLL_HEADER void MyError(const char * ch)
  {
    cout << ch;
    (*testout) << "Error !!! " << ch << endl << flush;
  }
#endif

  static clock_t starttimea;
  void ResetTime ()
  {
    starttimea = clock();
  }

#ifndef SMALLLIB
  DLL_HEADER double GetTime ()
  {
    return double(clock() - starttimea) / CLOCKS_PER_SEC;
  }
#endif





  // file handling ..
  int Ng_New (ClientData clientData,
	      Tcl_Interp * interp,
	      int argc, tcl_const char *argv[])
  {
    if (strcmp (argv[1], "mesh") == 0)
      mesh.Reset();


    if (strcmp (argv[1], "geom") == 0)
      {
	geometry.Reset (new CSGeometry (""));

	if (stlgeometry)
	  {
	    delete stlgeometry;
	    stlgeometry = NULL;
	  }
#ifdef OCCGEOMETRY
	if (occgeometry)
	  {
	    delete occgeometry;
	    occgeometry = NULL;
	  }
#endif
      }

    return TCL_OK;
  }




  int Ng_ImportMesh (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[]);

  int Ng_LoadMesh (ClientData clientData,
		   Tcl_Interp * interp,
		   int argc, tcl_const char *argv[])
  {
    string filename (argv[1]);

    if ( (strlen (filename.c_str()) > 4) &&
	 strcmp (&filename[strlen (filename.c_str())-4], ".vol") != 0 )
      {
	return Ng_ImportMesh(clientData,interp,argc,argv);
      }

    PrintMessage (1, "load mesh from file ", filename);

    mesh.Reset (new Mesh());

    try
      {
	//mesh -> Load (filename);
	ifstream infile(filename.c_str());
	mesh -> Load(infile);
	string auxstring;
	if(infile.good())
	  {
	    infile >> auxstring;
	    if(auxstring == "csgsurfaces")
	      {
		if (geometry)
		  geometry.Reset (new CSGeometry (""));
		if (stlgeometry)
		  {
		    delete stlgeometry;
		    stlgeometry = NULL;
		  }
#ifdef OCCGEOMETRY
		if (occgeometry)
		  {
		    delete occgeometry;
		    occgeometry = NULL;
		  }
#endif
		geometry2d.Reset (0);


		geometry -> LoadSurfaces(infile);
	      }
	  }
      }
    catch (NgException e)
      {
	PrintMessage (3, e.What());
	return TCL_ERROR;
      }

    PrintMessage (2,  mesh->GetNP(), " Points, ",
		  mesh->GetNE(), " Elements.");

    return TCL_OK;
  }





  int Ng_SaveMesh (ClientData clientData,
		   Tcl_Interp * interp,
		   int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    const string filename (argv[1]);
    PrintMessage (1, "Save mesh to file ", filename, ".... Please Wait!");

    ofstream outfile(filename.c_str());
    mesh -> Save (outfile);

    outfile << endl << endl << "endmesh" << endl << endl;

    if (geometry && geometry->GetNSurf()) geometry->SaveSurfaces(outfile);

    PrintMessage (1, "Save mesh to file .... DONE!");
    return TCL_OK;
  }





  int Ng_MergeMesh (ClientData clientData,
		    Tcl_Interp * interp,
		    int argc, tcl_const char *argv[])
  {
    string filename (argv[1]);

    PrintMessage (1, "merge with mesh from file ", filename);

    try
      {
	//mesh -> Merge (filename);
	ifstream infile(filename.c_str());
	const int offset = (geometry) ? geometry->GetNSurf() : 0;
	mesh -> Merge(infile,offset);
	string auxstring;
	if(infile.good())
	  {
	    infile >> auxstring;
	    if(auxstring == "csgsurfaces")
	      geometry -> LoadSurfaces(infile);
	  }
      }
    catch (NgException e)
      {
	PrintMessage (3, e.What());
	return TCL_ERROR;
      }

    PrintMessage (2,  mesh->GetNP(), " Points, ",
		  mesh->GetNSE(), " Surface Elements.");

    return TCL_OK;
  }





  int Ng_ExportMesh (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    string filename (argv[1]);
    string filetype (argv[2]);
    PrintMessage (1, "Export mesh to file ", filename, ".... Please Wait!");

    if (WriteUserFormat (filetype, *mesh, *geometry, filename))
      {
	ostringstream ost;
	ost << "Sorry, nothing known about file format " << filetype << endl;
	Tcl_SetResult (interp, (char*)ost.str().c_str(), TCL_VOLATILE);
	return TCL_ERROR;
      }

    PrintMessage (1, "Export mesh to file .... DONE!");
    return TCL_OK;
  }



  int Ng_ImportMesh (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
    const string filename (argv[1]);
    PrintMessage (1, "import mesh from ", filename);

    mesh.Reset (new Mesh());

    ReadFile (*mesh, filename);
    PrintMessage (2, mesh->GetNP(), " Points, ",
		  mesh->GetNE(), " Elements.");

    mesh->SetGlobalH (mparam.maxh);
    mesh->CalcLocalH();

    return TCL_OK;
  }



  int Ng_ImportSolution (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    const char * filename = argv[1];
    PrintMessage (1, "Import solution from file ", filename);

    ImportSolution (filename);
    return TCL_OK;
  }



  static DemoView * demoview = 0;
  int Ng_ShowDemo (ClientData clientData,
		   Tcl_Interp * interp,
		   int argc, tcl_const char *argv[])
  {
    const char * filename = argv[1];
    PrintMessage (1, "Show demo ", filename);
    demoview = new DemoView (filename);

    return TCL_OK;
  }



  int Ng_DemoSetTime (ClientData clientData,
		      Tcl_Interp * interp,
		      int argc, tcl_const char *argv[])
  {
    cout << "demosettime, time = " << argv[1] << endl;
    int result = -1;

    static char strminusone[] = "-1";
    static char str0[] = "0";

    if (demoview)
      result = demoview->SetTime (atof (argv[1]));

    if (result == -1)
      Tcl_SetResult (interp, strminusone, TCL_STATIC);
    else
      Tcl_SetResult (interp, str0, TCL_STATIC);

    return TCL_OK;
  }







  int Ng_SaveSolution (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    const char * filename = argv[1];
    PrintMessage (1, "Save solution to file ", filename);

    vssolution.SaveSolutionData (filename);
    return TCL_OK;
  }







  int Ng_LoadGeometry (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    tcl_const char * lgfilename = argv[1];

#ifdef LOG_STREAM
    (*logout) << "Load geometry file: " << lgfilename << endl;
#endif

#ifdef STAT_STREAM
    (*statout) << lgfilename << " & " << endl;
#endif

    if (geometry)
      geometry.Reset (new CSGeometry (""));
    if (stlgeometry)
      {
	delete stlgeometry;
	stlgeometry = NULL;
      }
#ifdef OCCGEOMETRY
    if (occgeometry)
      {
	delete occgeometry;
	occgeometry = NULL;
      }
#endif
    geometry2d.Reset (0);


    try
      {

	ifstream infile(lgfilename);

	if (strlen(lgfilename) < 4)
	  {
	    cout << "ERROR: cannot recognise file format!" << endl;
	  }
	else
	  {
	    if (strcmp (&lgfilename[strlen(lgfilename)-3], "geo") == 0)
	      {
		// strcpy (geomfilename, lgfilename);
		PrintMessage (1, "Load geometry file ", lgfilename);


		extern CSGeometry * ParseCSG (istream & istr);
		// ifstream infile(lgfilename);
		CSGeometry * hgeom = ParseCSG (infile);
		ng_geometry = hgeom;
		if (hgeom)
		  geometry.Reset (hgeom);
		else
		  {
		    geometry.Reset (new CSGeometry (""));
		    Tcl_SetResult (interp, (char*)"geo-file should start with 'algebraic3d'", TCL_STATIC);
		    return TCL_ERROR;
		  }

		//geometry -> FindIdenticSurfaces(geometry->GetIdEps() * geometry->MaxSize()); // 1e-8*geometry->MaxSize()
		geometry -> FindIdenticSurfaces(1e-8 * geometry->MaxSize()); // 1e-8*geometry->MaxSize()
	      }
	    else if (strcmp (&lgfilename[strlen(lgfilename)-3], "ngg") == 0)
	      {
		//		strcpy (geomfilename, lgfilename);

		PrintMessage (1, "Load new geometry file ", lgfilename);
		geometry.Reset (new CSGeometry(""));
		geometry -> Load (infile);
	      }

	    // strcpy (geomfilename, lgfilename);
	    // (*mycout) << "Load geometry file " << lgfilename << endl;

	    else if (strcmp (&lgfilename[strlen(lgfilename)-3], "stl") == 0)
	      {
		// strcpy (geomfilename, lgfilename);
		PrintMessage (1, "Load stl geometry file ", lgfilename);
		stlgeometry = STLGeometry :: Load (infile);
		ng_geometry = stlgeometry;
		stlgeometry->edgesfound = 0;
	      }
	    else if ((strcmp (&lgfilename[strlen(lgfilename)-4], "iges") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-3], "igs") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-3], "IGS") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-4], "IGES") == 0))
	      {
#ifdef OCCGEOMETRY
		// strcpy (geomfilename, lgfilename);
		PrintMessage (1, "Load IGES geometry file ", lgfilename);
		occgeometry = LoadOCC_IGES (lgfilename);
		ng_geometry = occgeometry;
#else
		Tcl_SetResult (interp, (char*)"IGES import requires the OpenCascade geometry kernel. "
			       "Please install OpenCascade as described in the Netgen-website",
			       TCL_STATIC);
		return TCL_ERROR;
#endif
	      }

	    else if (strcmp (&lgfilename[strlen(lgfilename)-3], "sat") == 0)
	      {
#ifdef ACIS
		PrintMessage (1, "Load ACIS geometry file ", lgfilename);
		acisgeometry = netgen::LoadACIS_SAT (lgfilename);
#endif
	      }
	    else if ((strcmp (&lgfilename[strlen(lgfilename)-4], "step") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-3], "stp") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-3], "STP") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-4], "STEP") == 0))
	      {
#ifdef ACISxxx
		PrintMessage (1, "Load STEP geometry file ", lgfilename);
		acisgeometry = netgen::LoadACIS_STEP (lgfilename);
#else
#ifdef OCCGEOMETRY
		// strcpy (geomfilename, lgfilename);
		PrintMessage (1, "Load STEP geometry file ", lgfilename);
		occgeometry = LoadOCC_STEP (lgfilename);
		ng_geometry = occgeometry;
#else
		Tcl_SetResult (interp, (char*)"IGES import requires the OpenCascade geometry kernel. "
			       "Please install OpenCascade as described in the Netgen-website",
			       TCL_STATIC);
		return TCL_ERROR;
#endif
#endif
	      }
	    else if ((strcmp (&lgfilename[strlen(lgfilename)-4], "brep") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-4], "Brep") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-4], "BREP") == 0))
	      {
#ifdef OCCGEOMETRY
		// strcpy (geomfilename, lgfilename);
		PrintMessage (1, "Load BREP geometry file ", lgfilename);
		occgeometry = LoadOCC_BREP (lgfilename);
		ng_geometry = occgeometry;
#else
		Tcl_SetResult (interp, (char*)"BREP import requires the OpenCascade geometry kernel. "
			       "Please install OpenCascade as described in the Netgen-website",
			       TCL_STATIC);
		return TCL_ERROR;
#endif
	      }

	    else if (strcmp (&lgfilename[strlen(lgfilename)-4], "stlb") == 0)
	      {
		// strcpy (geomfilename, lgfilename);

		PrintMessage (1, "Load stl geometry file ", lgfilename, " in binary format");
		stlgeometry = STLGeometry :: LoadBinary (infile);
		ng_geometry = stlgeometry;
		stlgeometry->edgesfound = 0;
	      }

	    else if (strcmp (&lgfilename[strlen(lgfilename)-3], "nao") == 0)
	      {
		// strcpy (geomfilename, lgfilename);

		PrintMessage (1, "Load naomi (F. Kickinger) geometry file ", lgfilename);
		stlgeometry = STLGeometry :: LoadNaomi (infile);
		ng_geometry = stlgeometry;
		stlgeometry->edgesfound = 0;
	      }

	    else if (strcmp (&lgfilename[strlen(lgfilename)-4], "in2d") == 0)
	      {
		// strcpy (geomfilename, lgfilename);
		geometry2d.Reset (new SplineGeometry2d());
		geometry2d -> Load (lgfilename);
		ng_geometry = geometry2d.Ptr();
	      }
	  }
      }
    catch (NgException e)
      {
	Tcl_SetResult (interp, const_cast<char*> (e.What().c_str()), TCL_VOLATILE);
	return TCL_ERROR;
      }

    mesh.Reset();

    return TCL_OK;
  }










  int Ng_SaveGeometry (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    if (argc == 2)
      {
	const char * cfilename = argv[1];
	PrintMessage (1, "Save geometry to file ", cfilename);

	if (strlen(cfilename) < 4) {cout << "ERROR: can not recognise file format!!!" << endl;}
	else
	  {
#ifdef ACIS
	    if (acisgeometry)
	      {
		char * filename = const_cast<char*> (argv[1]);
		if (strcmp (&filename[strlen(filename)-3], "sat") == 0)
		  {
		    acisgeometry -> SaveSATFile (filename);
		  }
	      }
#endif
#ifdef OCCGEOMETRY
	    if (occgeometry)
	      {
		char * filename = const_cast<char*> (argv[1]);
		if (strcmp (&filename[strlen(filename)-3], "igs") == 0)
		  {
		    IGESControl_Writer writer("millimeters", 1);
		    writer.AddShape (occgeometry->shape);
		    writer.Write (filename);
		  }
		else if (strcmp (&filename[strlen(filename)-3], "stp") == 0)
		  {
		    STEPControl_Writer writer;
		    writer.Transfer (occgeometry->shape, STEPControl_AsIs);
		    writer.Write (filename);
		  }
		else if (strcmp (&filename[strlen(filename)-3], "stl") == 0)
		  {
		    StlAPI_Writer writer;
		    writer.ASCIIMode() = Standard_True;
		    writer.Write (occgeometry->shape, filename);
		  }
		else if (strcmp (&filename[strlen(filename)-4], "stlb") == 0)
		  {
		    StlAPI_Writer writer;
		    writer.ASCIIMode() = Standard_False;
		    writer.Write (occgeometry->shape, filename);
		  }
	      }
	    else
#endif
	      if (strcmp (&cfilename[strlen(cfilename)-3], "ngg") == 0)
		{
		  if (geometry)
		    {
		      ofstream of(cfilename);
		      geometry->Save (of);
		    }
		}
	      else if (strlen(cfilename) > 3 &&
		       strcmp (&cfilename[strlen(cfilename)-3], "stl") == 0)
		{
		  if (stlgeometry)
		    stlgeometry->Save (cfilename);
		}
	      else if (strlen(cfilename) > 4 &&
		       strcmp (&cfilename[strlen(cfilename)-4], "stlb") == 0)
		{
		  if (stlgeometry)
		    stlgeometry->SaveBinary (cfilename,"Binary STL Geometry");
		}
	      else if (strlen(cfilename) > 4 &&
		       strcmp (&cfilename[strlen(cfilename)-4], "stle") == 0)
		{
		  if (stlgeometry)
		    stlgeometry->SaveSTLE (cfilename);
		}
	  }
      }

    else
      {
	if (geometry)
	  geometry->Save (cout);
      }
    return TCL_OK;
  }



  int Ng_ParseGeometry (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
#ifdef OCCGEOMETRY
    if (!stlgeometry && !geometry2d && !occgeometry)
#else
      if (!stlgeometry && !geometry2d)
#endif
	{
          double detail = atof (Tcl_GetVar (interp, "::geooptions.detail", 0));
          double facets = atof (Tcl_GetVar (interp, "::geooptions.facets", 0));
	  Box<3> box (geometry->BoundingBox());

          if (atoi (Tcl_GetVar (interp, "::geooptions.drawcsg", 0)))
	    geometry->CalcTriangleApproximation(box, detail, facets);
	}
    return TCL_OK;
  }




  /*
    NgLock * ngpkg_lock = NULL;


    int Ng_Lock (ClientData clientData,
    Tcl_Interp * interp,
    int argc, tcl_const char *argv[])
    {
    delete ngpkg_lock;
    ngpkg_lock = NULL;

    if(strcmp (argv[1], "mesh") == 0)
    {
    ngpkg_lock = new NgLock(mesh->Mutex());
    }
    else if(strcmp (argv[1], "unlock") == 0)
    {
    ;
    }

    return TCL_OK;
    }
  */

  int Ng_GeometryOptions (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
    const char * command = argv[1];

    if (strcmp (command, "get") == 0)
      {
	char buf[20];
	Point3d pmin = geometry->BoundingBox ().PMin();
	Point3d pmax = geometry->BoundingBox ().PMax();

	sprintf (buf, "%5.1lf", pmin.X());
        Tcl_SetVar (interp, "::geooptions.minx", buf, 0);
	sprintf (buf, "%5.1lf", pmin.Y());
        Tcl_SetVar (interp, "::geooptions.miny", buf, 0);
	sprintf (buf, "%5.1lf", pmin.Z());
        Tcl_SetVar (interp, "::geooptions.minz", buf, 0);

	sprintf (buf, "%5.1lf", pmax.X());
        Tcl_SetVar (interp, "::geooptions.maxx", buf, 0);
	sprintf (buf, "%5.1lf", pmax.Y());
        Tcl_SetVar (interp, "::geooptions.maxy", buf, 0);
	sprintf (buf, "%5.1lf", pmax.Z());
        Tcl_SetVar (interp, "::geooptions.maxz", buf, 0);
      }
    else if (strcmp (command, "set") == 0)
      {
        Point<3> pmin (atof (Tcl_GetVar (interp, "::geooptions.minx", 0)),
                       atof (Tcl_GetVar (interp, "::geooptions.miny", 0)),
                       atof (Tcl_GetVar (interp, "::geooptions.minz", 0)));
        Point<3> pmax (atof (Tcl_GetVar (interp, "::geooptions.maxx", 0)),
                       atof (Tcl_GetVar (interp, "::geooptions.maxy", 0)),
                       atof (Tcl_GetVar (interp, "::geooptions.maxz", 0)));
	Box<3> box (pmin, pmax);
	geometry -> SetBoundingBox (box);
	CSGeometry::SetDefaultBoundingBox (box);
      }

    return TCL_OK;
  }





  // attempt of a simple modeller

  int Ng_CreatePrimitive (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
    tcl_const char * classname = argv[1];
    tcl_const char * name = argv[2];

    cout << "Create primitive, class = " << classname
	 << ", name = " << name << endl;

    Primitive * nprim = Primitive::CreatePrimitive (classname);
    Solid * nsol = new Solid (nprim);

    char sname[100];
    for (int j = 1; j <= nprim->GetNSurfaces(); j++)
      {
	sprintf (sname, "%s,%d", name, j);
	geometry -> AddSurface (sname, &nprim->GetSurface(j));
	nprim -> SetSurfaceId (j, geometry->GetNSurf());
      }

    geometry->SetSolid (name, nsol);

    return TCL_OK;
  }


  int Ng_SetPrimitiveData (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    tcl_const char * name = argv[1];
    tcl_const char * value = argv[2];

    Array<double> coeffs;


    cout << "Set primitive data, name = " << name
	 << ", value = " << value  << endl;


    istringstream vst (value);
    double val;
    while (!vst.eof())
      {
	vst >> val;
	coeffs.Append (val);
      }

    ((Primitive*)
     geometry->GetSolid (name)->GetPrimitive())->SetPrimitiveData (coeffs);

    return TCL_OK;
  }



  int Ng_SetSolidData (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    tcl_const char * name = argv[1];
    tcl_const char * val = argv[2];

    cout << "Set Solid Data, name = " << name
	 << ", value = " << val << endl;

    istringstream vst (val);

    Solid * nsol = Solid::CreateSolid (vst, geometry->GetSolids());
    geometry->SetSolid (name, nsol);

    return TCL_OK;
  }


  int Ng_GetPrimitiveData (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    tcl_const char * name = argv[1];
    tcl_const char * classnamevar = argv[2];
    tcl_const char * valuevar = argv[3];

    const char * classname;

    Array<double> coeffs;

    geometry->GetSolid (name)->GetPrimitive()->GetPrimitiveData (classname, coeffs);

    ostringstream vst;
    for (int i = 1; i <= coeffs.Size(); i++)
      vst << coeffs.Get(i) << " ";

    cout << "GetPrimitiveData, name = " << name
	 << ", classnamevar = " << classnamevar
	 << ", classname = " << classname << endl
	 << " valuevar = " << valuevar
	 << ", values = " << vst.str() << endl;

    Tcl_SetVar  (interp, classnamevar, (char*)classname, 0);
    Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);

    return TCL_OK;
  }

  int Ng_GetSolidData (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    tcl_const char * name = argv[1];
    tcl_const char * valuevar = argv[2];

    ostringstream vst;

    const Solid * sol = geometry->GetSolid (name);
    sol->GetSolidData (vst);

    cout << "GetSolidData, name = " << name << ", data = " << vst.str() << endl;

    Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);

    return TCL_OK;
  }


  int Ng_GetPrimitiveList (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    tcl_const char * valuevar = argv[1];
    int i;

    stringstream vst;

    for (i = 1; i <= geometry->GetNSolids(); i++)
      {
	const Solid * sol = geometry->GetSolid(i);
	if (sol->GetPrimitive())
	  vst << sol->Name() << " ";
      }

    cout << "primnames = " << vst.str() << endl;

    Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);

    return TCL_OK;
  }



  int Ng_GetSurfaceList (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[])
  {
    tcl_const char * valuevar = argv[1];
    int i;

    stringstream vst;

    for (i = 1; i <= geometry->GetNSurf(); i++)
      {
	const Surface * surf = geometry->GetSurface(i);
	vst << surf->Name() << " ";
      }

    cout << "surfnames = " << vst.str() << endl;

    Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);

    return TCL_OK;
  }


  int Ng_GetSolidList (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    tcl_const char * valuevar = argv[1];
    int i;

    stringstream vst;

    for (i = 1; i <= geometry->GetNSolids(); i++)
      {
	const Solid * sol = geometry->GetSolid(i);
	if (!sol->GetPrimitive())
	  vst << sol->Name() << " ";
      }

    cout << "solnames = " << vst.str() << endl;

    Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);

    return TCL_OK;
  }


  int Ng_TopLevel (ClientData clientData,
		   Tcl_Interp * interp,
		   int argc, tcl_const char *argv[])
  {
    int i;
    /*
      for (i = 0; i < argc; i++)
      cout << argv[i] << ", ";
      cout << endl;
    */

    if (strcmp (argv[1], "getlist") == 0)
      {
	stringstream vst;

	for (i = 0; i < geometry->GetNTopLevelObjects(); i++)
	  {
	    const Solid * sol;
	    const Surface * surf;
	    geometry->GetTopLevelObject (i, sol, surf);

	    if (!surf)
	      vst << "{ " << sol->Name() << " } ";
	    else
	      vst << "{ " << sol->Name() << " " << surf->Name() << " } ";
	  }

	tcl_const char * valuevar = argv[2];
	Tcl_SetVar  (interp, valuevar, (char*)vst.str().c_str(), 0);
      }

    if (strcmp (argv[1], "set") == 0)
      {
	tcl_const char * solname = argv[2];
	tcl_const char * surfname = argv[3];
	Solid * sol = (Solid*)geometry->GetSolid (solname);
	Surface * surf = (Surface*)geometry->GetSurface (surfname);
	geometry->SetTopLevelObject (sol, surf);
      }

    if (strcmp (argv[1], "remove") == 0)
      {
	tcl_const char * solname = argv[2];
	tcl_const char * surfname = argv[3];
	Solid * sol = (Solid*)geometry->GetSolid (solname);
	Surface * surf = (Surface*)geometry->GetSurface (surfname);
	geometry->RemoveTopLevelObject (sol, surf);
      }

    if (strcmp (argv[1], "setprop") == 0)
      {
	tcl_const char * solname = argv[2];
	tcl_const char * surfname = argv[3];
	tcl_const char * propvar = argv[4];
	Solid * sol = (Solid*)geometry->GetSolid (solname);
	Surface * surf = (Surface*)geometry->GetSurface (surfname);
	TopLevelObject * tlo = geometry->GetTopLevelObject (sol, surf);

	if (!tlo) return TCL_OK;

	char varname[50];
	sprintf (varname, "%s(red)", propvar);
	double red = atof (Tcl_GetVar (interp, varname, 0));
	sprintf (varname, "%s(blue)", propvar);
	double blue = atof (Tcl_GetVar (interp, varname, 0));
	sprintf (varname, "%s(green)", propvar);
	double green = atof (Tcl_GetVar (interp, varname, 0));
	tlo -> SetRGB (red, green, blue);

	sprintf (varname, "%s(visible)", propvar);
	tlo -> SetVisible (bool(atoi (Tcl_GetVar (interp, varname, 0))));
	sprintf (varname, "%s(transp)", propvar);
	tlo -> SetTransparent (bool(atoi (Tcl_GetVar (interp, varname, 0))));
      }

    if (strcmp (argv[1], "getprop") == 0)
      {
	tcl_const char * solname = argv[2];
	tcl_const char * surfname = argv[3];
	tcl_const char * propvar = argv[4];

	Solid * sol = (Solid*)geometry->GetSolid (solname);
	Surface * surf = (Surface*)geometry->GetSurface (surfname);
	TopLevelObject * tlo = geometry->GetTopLevelObject (sol, surf);

	if (!tlo) return TCL_OK;

	char varname[50], varval[10];

	sprintf (varname, "%s(red)", propvar);
	sprintf (varval, "%lf", tlo->GetRed());
	Tcl_SetVar (interp, varname, varval, 0);

	sprintf (varname, "%s(green)", propvar);
	sprintf (varval, "%lf", tlo->GetGreen());
	Tcl_SetVar (interp, varname, varval, 0);

	sprintf (varname, "%s(blue)", propvar);
	sprintf (varval, "%lf", tlo->GetBlue());
	Tcl_SetVar (interp, varname, varval, 0);

	sprintf (varname, "%s(visible)", propvar);
	sprintf (varval, "%d", tlo->GetVisible());
	Tcl_SetVar (interp, varname, varval, 0);

	sprintf (varname, "%s(transp)", propvar);
	sprintf (varval, "%d", tlo->GetTransparent());
	Tcl_SetVar (interp, varname, varval, 0);
      }


    return TCL_OK;
  }






  int Ng_ReadStatus (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
    char buf[20], lstring[200];
    if (mesh.Ptr())
      {
	sprintf (buf, "%d", mesh->GetNP());
        Tcl_SetVar  (interp, "::status_np", buf, 0);
	sprintf (buf, "%d", mesh->GetNE());
        Tcl_SetVar  (interp, "::status_ne", buf, 0);
	sprintf (buf, "%d", mesh->GetNSE());
        Tcl_SetVar  (interp, "::status_nse", buf, 0);
      }
    else
      {
        Tcl_SetVar  (interp, "::status_np", "0", 0);
        Tcl_SetVar  (interp, "::status_ne", "0", 0);
        Tcl_SetVar  (interp, "::status_nse", "0", 0);
      }

    if (multithread.running)
      Tcl_SetVar (interp, "::status_working", "working", 0);
    else
      Tcl_SetVar (interp, "::status_working", "       ", 0);

    Tcl_SetVar (interp, "::status_task", const_cast<char *>(multithread.task), 0);
    sprintf (buf, "%lf", multithread.percent);
    Tcl_SetVar  (interp, "::status_percent", buf, 0);

    int i;
    lstring[0] = 0;
    for (i = 1; i <= tets_in_qualclass.Size(); i++)
      {
	sprintf (buf, " %d", tets_in_qualclass.Get(i));
	strcat (lstring, buf);
      }
    for (i = tets_in_qualclass.Size()+1; i <= 20; i++)
      strcat (lstring, " 0");
    Tcl_SetVar  (interp, "::status_tetqualclasses", lstring, 0);

    return TCL_OK;
  }


  int Ng_MemInfo (ClientData clientData,
		  Tcl_Interp * interp,
		  int argc, tcl_const char *argv[])
  {
    if (argc < 2) return TCL_ERROR;

    if (strcmp (argv[1], "usedmb") == 0)
      { // returns string of 512 '0' or '1'

	static char usedmb[513];
	for (int i = 0; i < 512; i++)
	  usedmb[i] = (i % 7 == 0) ? '1' : '0';

	usedmb[512] = 0;
	BaseDynamicMem::GetUsed (512, usedmb);
	Tcl_SetResult (interp, usedmb, TCL_STATIC);
	return TCL_OK;
      }

    return TCL_ERROR;
  }



  int Ng_BCProp (ClientData clientData,
		 Tcl_Interp * interp,
		 int argc, tcl_const char *argv[])
  {
    static char buf[100];

    if (argc < 2)
      {
	Tcl_SetResult (interp, (char*)"Ng_BCProp needs arguments", TCL_STATIC);
	return TCL_ERROR;
      }

    if (strcmp (argv[1], "setbc") == 0)
      {
	int facenr = atoi (argv[2]);
	int bcnr = atoi (argv[3]);
	if (mesh.Ptr() && facenr >= 1 && facenr <= mesh->GetNFD())
	  mesh->GetFaceDescriptor (facenr).SetBCProperty (bcnr);
      }

    if (strcmp (argv[1], "setall") == 0)
      {
	int bcnr = atoi (argv[2]);
	if (mesh.Ptr())
	  {
	    int nfd = mesh->GetNFD();
	    for (int i = 1; i <= nfd; i++)
	      mesh->GetFaceDescriptor (i).SetBCProperty (bcnr);
	  }
      }

    if (strcmp (argv[1], "getbc") == 0)
      {
	int facenr = atoi (argv[2]);
	if (mesh.Ptr() && facenr >= 1 && facenr <= mesh->GetNFD())
	  {
	    sprintf (buf, "%d", mesh->GetFaceDescriptor(facenr).BCProperty());
	  }
	else
	  {
	    strcpy (buf, "0");
	  }
	Tcl_SetResult (interp, buf, TCL_STATIC);
      }

    if (strcmp (argv[1], "getbcname") == 0)
      {
	int facenr = atoi (argv[2]);
	if (mesh.Ptr() && facenr >= 1 && facenr <= mesh->GetNFD())
	  {
	    sprintf (buf, "%s", mesh->GetFaceDescriptor(facenr).GetBCName().c_str());
	  }
	else
	  {
	    strcpy (buf, "-");
	  }
	Tcl_SetResult (interp, buf, TCL_STATIC);
      }


    if (strcmp (argv[1], "getactive") == 0)
      {
	sprintf (buf, "%d", vsmesh.SelectedFace());
	Tcl_SetResult (interp, buf, TCL_STATIC);
      }

    if (strcmp (argv[1], "setactive") == 0)
      {
	int facenr = atoi (argv[2]);
	if (mesh.Ptr() && facenr >= 1 && facenr <= mesh->GetNFD())
	  {
	    vsmesh.SetSelectedFace (facenr);
	  }
      }

    if (strcmp (argv[1], "getnfd") == 0)
      {
	if (mesh.Ptr())
	  sprintf (buf, "%d", mesh->GetNFD());
	else
	  sprintf (buf, "0");
	Tcl_SetResult (interp, buf, TCL_STATIC);
      }

    return TCL_OK;
  }



  // Philippose - 30/01/2009
  // TCL interface function for the Local Face Mesh size
  // definition functionality
  int Ng_SurfaceMeshSize (ClientData clientData,
		                    Tcl_Interp * interp,
		                    int argc, tcl_const char *argv[])
  {
#ifdef OCCGEOMETRY

    static char buf[100];

    if (argc < 2)
    {
	   Tcl_SetResult (interp, (char *)"Ng_SurfaceMeshSize needs arguments", TCL_STATIC);
	   return TCL_ERROR;
    }

    if (!occgeometry)
    {
      Tcl_SetResult (interp, (char *)"Ng_SurfaceMeshSize currently supports only OCC (STEP/IGES) Files", TCL_STATIC);
	   return TCL_ERROR;
    }

    // Update the face mesh sizes to reflect the global maximum mesh size
    for(int i = 1; i <= occgeometry->NrFaces(); i++)
    {
           if(!occgeometry->GetFaceMaxhModified(i))
           {
              occgeometry->SetFaceMaxH(i, mparam.maxh);
           }   
    }

    if (strcmp (argv[1], "setsurfms") == 0)
    {
	   int facenr = atoi (argv[2]);
	   double surfms = atof (argv[3]);
	   if (occgeometry && facenr >= 1 && facenr <= occgeometry->NrFaces())
	     occgeometry->SetFaceMaxH(facenr, surfms);

    }

    if (strcmp (argv[1], "setall") == 0)
    {
	   double surfms = atof (argv[2]);
	   if (occgeometry)
	   {
	     int nrFaces = occgeometry->NrFaces();
	     for (int i = 1; i <= nrFaces; i++)
	      occgeometry->SetFaceMaxH(i, surfms);
	   }
    }

    if (strcmp (argv[1], "getsurfms") == 0)
    {
	   int facenr = atoi (argv[2]);
	   if (occgeometry && facenr >= 1 && facenr <= occgeometry->NrFaces())
	   {
	     sprintf (buf, "%5.2f", occgeometry->GetFaceMaxH(facenr));
	   }
	   else
	   {
	     sprintf (buf, "%5.2f", mparam.maxh);
	   }
	   Tcl_SetResult (interp, buf, TCL_STATIC);
    }

    if (strcmp (argv[1], "getactive") == 0)
    {
	   sprintf (buf, "%d", occgeometry->SelectedFace());
	   Tcl_SetResult (interp, buf, TCL_STATIC);
    }

    if (strcmp (argv[1], "setactive") == 0)
    {
	   int facenr = atoi (argv[2]);
	   if (occgeometry && facenr >= 1 && facenr <= occgeometry->NrFaces())
	   {
	     occgeometry->SetSelectedFace (facenr);

        occgeometry->LowLightAll();
        occgeometry->fvispar[facenr-1].Highlight();
        occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
	   }
    }

    if (strcmp (argv[1], "getnfd") == 0)
    {
	   if (occgeometry)
	     sprintf (buf, "%d", occgeometry->NrFaces());
	   else
	     sprintf (buf, "0");
	   Tcl_SetResult (interp, buf, TCL_STATIC);
    }
    return TCL_OK;
#else // No OCCGEOMETRY 

    Tcl_SetResult (interp, (char *)"Ng_SurfaceMeshSize currently supports only OCC (STEP/IGES) Files", TCL_STATIC);
    return TCL_ERROR;
    
#endif // OCCGEOMETRY
  }




  // Philippose - 10/03/2009
  // TCL interface function for the Automatic Colour-based
  // definition of boundary conditions for OCC Geometry
  int Ng_AutoColourBcProps (ClientData clientData,
		                      Tcl_Interp * interp,
		                      int argc, tcl_const char *argv[])
  {
     if(argc < 1)
     {
        Tcl_SetResult (interp, (char *)"Ng_AutoColourBcProps needs arguments", TCL_STATIC);
        return TCL_ERROR;
     }

     if(!mesh.Ptr())
     {
        Tcl_SetResult (interp, (char *)"Ng_AutoColourBcProps: Valid netgen mesh required...please mesh the Geometry first", TCL_STATIC);
	     return TCL_ERROR;
     }

     if(strcmp(argv[1], "auto") == 0)
     {
        AutoColourBcProps(*mesh, 0);
     }

     if(strcmp(argv[1], "profile") == 0)
     {
        AutoColourBcProps(*mesh, argv[2]);
     }

     return TCL_OK;
  }
                          
                          
  
  int Ng_SetNextTimeStamp  (ClientData clientData,
			    Tcl_Interp * interp,
			    int argqc, tcl_const char *argv[])
  {
    if (mesh.Ptr())
      mesh -> SetNextTimeStamp();
    return TCL_OK;
  }



  int Ng_Refine  (ClientData clientData,
		  Tcl_Interp * interp,
		  int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

#ifdef ACIS
    if (acisgeometry)
      {
	ACISRefinementSurfaces ref (*acisgeometry);
	ACISMeshOptimize2dSurfaces opt(*acisgeometry);
	ref.Set2dOptimizer(&opt);
	ref.Refine (*mesh);
      }
#endif
    else
      {
	ng_geometry -> GetRefinement().Refine(*mesh);
      }

    return TCL_OK;
  }

  int Ng_SecondOrder  (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }
    
    const_cast<Refinement&> (ng_geometry -> GetRefinement()).MakeSecondOrder (*mesh);

    return TCL_OK;
  }


  void * HighOrderDummy (void *)
  {
    //  mparam.elementorder = atoi (Tcl_GetVar (interp, "options.elementorder", 0));
    const char * savetask = multithread.task;

    Refinement & ref = const_cast<Refinement&> (ng_geometry -> GetRefinement());
    mesh -> GetCurvedElements().BuildCurvedElements (&ref, mparam.elementorder);

    multithread.task = savetask;
    multithread.running = 0;
    multithread.terminate = 1;

    mesh -> SetNextMajorTimeStamp();
    return 0;
  }

  int Ng_HighOrder  (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    multithread.running = 1;
    multithread.terminate = 0;

    mparam.elementorder = atoi(argv[1]);


    HighOrderDummy(NULL);

    return TCL_OK;
  }



  void * ValidateDummy (void *)
  {
    RefinementSTLGeometry ref (*stlgeometry);
    ref.ValidateSecondOrder (*mesh);

    multithread.running = 0;
    return NULL;
  }



  int Ng_ValidateSecondOrder  (ClientData clientData,
			       Tcl_Interp * interp,
			       int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    multithread.running = 1;
    RunParallel (ValidateDummy, NULL);

    return TCL_OK;
  }


  int Ng_ZRefinement  (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    ZRefinementOptions opt;
    opt.minref = 5;

    if (argc >= 2) opt.minref = atoi (argv[1]);

    ZRefinement (*mesh, geometry.Ptr(), opt);

    return TCL_OK;
  }

  int Ng_HPRefinement  (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    int levels = atoi (argv[1]);


    Refinement & ref = const_cast<Refinement&> (ng_geometry -> GetRefinement());
    HPRefinement (*mesh, &ref, levels);
    return TCL_OK;
  }


  int Ng_LoadMeshSize  (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    mesh->LoadLocalMeshSize(argv[1]);
    return TCL_OK;
  }


  int Ng_MeshSizeFromSurfaceMesh  (ClientData clientData,
				   Tcl_Interp * interp,
				   int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    mesh->SetGlobalH (mparam.maxh);
    mesh->CalcLocalH();

    return TCL_OK;
  }

  int Ng_SingularPointMS (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
    double globh = mparam.maxh;
    for (int i = 1; i <= geometry->singpoints.Size(); i++)
      geometry->singpoints.Get(i)->SetMeshSize (*mesh, globh);
    return TCL_OK;
  }



  int Ng_SingularEdgeMS (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    double globh = mparam.maxh;
    for (int i = 1; i <= geometry->singedges.Size(); i++)
      geometry->singedges.Get(i)->SetMeshSize (*mesh, globh);
    return TCL_OK;
  }


  // Philippose Rajan - 13 June 2009
  // Added a new TCL function call for the generation 
  // of prismatic boundary layers
  int Ng_GenerateBoundaryLayer (ClientData clientData,
           Tcl_Interp * interp,
           int argc, tcl_const char *argv[])
  {
     if (!mesh.Ptr())
     {
        Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
        return TCL_ERROR;
     }

     if(multithread.running)
     {
        Tcl_SetResult(interp, err_jobrunning, TCL_STATIC);
        return TCL_ERROR;
     }

     GenerateBoundaryLayer (*mesh, mparam);
     return TCL_OK;
  }


  int Ng_InsertVirtualBL (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    InsertVirtualBoundaryLayer (*mesh);
    return TCL_OK;
  }

  int Ng_CutOffAndCombine (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    Mesh othermesh;
    othermesh.Load (argv[1]);
    othermesh.SetGlobalH (mparam.maxh);
    othermesh.CalcLocalH();

    CutOffAndCombine (*mesh, othermesh);
    return TCL_OK;
  }


  int Ng_HelmholtzMesh (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    HelmholtzMesh (*mesh);
    return TCL_OK;
  }




  int Ng_SetMeshingParameters  (ClientData clientData,
				Tcl_Interp * interp,
				int argc, tcl_const char *argv[])
  {
    mparam.maxh = atof (Tcl_GetVar (interp, "::options.meshsize", 0));
    mparam.minh = atof (Tcl_GetVar (interp, "::options.minmeshsize", 0));
    mparam.meshsizefilename = Tcl_GetVar (interp, "::options.meshsizefilename", 0);
    if (!strlen (mparam.meshsizefilename))
      mparam.meshsizefilename = NULL;

    mparam.curvaturesafety = atof (Tcl_GetVar (interp, "::options.curvaturesafety", 0));
    mparam.segmentsperedge = atof (Tcl_GetVar (interp, "::options.segmentsperedge", 0));
    mparam.badellimit = atof (Tcl_GetVar (interp, "::options.badellimit", 0));
    mparam.secondorder = atoi (Tcl_GetVar (interp, "::options.secondorder", 0));
    mparam.elementorder = atoi (Tcl_GetVar (interp, "::options.elementorder", 0));
    mparam.quad = atoi (Tcl_GetVar (interp, "::options.quad", 0));

    mparam.inverttets = atoi (Tcl_GetVar (interp, "::options.inverttets", 0));
    mparam.inverttrigs = atoi (Tcl_GetVar (interp, "::options.inverttrigs", 0));
    mparam.uselocalh = atoi (Tcl_GetVar (interp, "::options.localh", 0));
    mparam.grading = atof (Tcl_GetVar (interp, "::options.grading", 0));
    mparam.delaunay = atoi (Tcl_GetVar (interp, "::options.delaunay", 0));
    mparam.checkoverlap = atoi (Tcl_GetVar (interp, "::options.checkoverlap", 0));
    mparam.checkoverlappingboundary = atoi (Tcl_GetVar (interp, "::options.checkoverlappingboundary", 0));
    mparam.checkchartboundary = atoi (Tcl_GetVar (interp, "::options.checkchartboundary", 0));
    mparam.optsteps3d = atoi (Tcl_GetVar (interp, "::options.optsteps3d", 0));
    mparam.optsteps2d = atoi (Tcl_GetVar (interp, "::options.optsteps2d", 0));
    mparam.opterrpow = atof (Tcl_GetVar (interp, "::options.opterrpow", 0));

    mparam.parthread = atoi (Tcl_GetVar (interp, "::options.parthread", 0));
    mparam.elsizeweight = atof (Tcl_GetVar (interp, "::options.elsizeweight", 0));

    mparam.autozrefine = atoi (Tcl_GetVar (interp, "::options.autozrefine", 0));

    extern int printmessage_importance;
    extern int printdots;
    printmessage_importance = atoi (Tcl_GetVar (interp, "::options.printmsg", 0));
    printdots = (printmessage_importance >= 4);

    //BaseMoveableMem::totalsize = 0;
    // 1048576 * atoi (Tcl_GetVar (interp, "::options.memory", 0));
    if (mesh.Ptr())
      {
	mesh->SetGlobalH (mparam.maxh);
	mesh->SetMinimalH (mparam.minh);
      }

    return TCL_OK;
  }



  int Ng_SetDebugParameters  (ClientData clientData,
			      Tcl_Interp * interp,
			      int argc, tcl_const char *argv[])
  {
    debugparam.slowchecks = atoi (Tcl_GetVar (interp, "::debug.slowchecks", 0));
    debugparam.debugoutput = atoi (Tcl_GetVar (interp, "::debug.debugoutput", 0));
    debugparam.haltexistingline = atoi (Tcl_GetVar (interp, "::debug.haltexistingline", 0));
    debugparam.haltoverlap = atoi (Tcl_GetVar (interp, "::debug.haltoverlap", 0));
    debugparam.haltsuccess = atoi (Tcl_GetVar (interp, "::debug.haltsuccess", 0));
    debugparam.haltnosuccess = atoi (Tcl_GetVar (interp, "::debug.haltnosuccess", 0));
    debugparam.haltlargequalclass = atoi (Tcl_GetVar (interp, "::debug.haltlargequalclass", 0));
    debugparam.haltsegment = atoi (Tcl_GetVar (interp, "::debug.haltsegment", 0));
    debugparam.haltnode = atoi (Tcl_GetVar (interp, "::debug.haltnode", 0));
    debugparam.haltface = atoi (Tcl_GetVar (interp, "::debug.haltface", 0));
    debugparam.haltsegmentp1 = atoi (Tcl_GetVar (interp, "::debug.haltsegmentp1", 0));
    debugparam.haltsegmentp2 = atoi (Tcl_GetVar (interp, "::debug.haltsegmentp2", 0));
    debugparam.haltfacenr = atoi (Tcl_GetVar (interp, "::debug.haltfacenr", 0));
    return TCL_OK;
  }



  int Ng_SetSTLParameters  (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[])
  {

    stlparam.yangle =
      atof (Tcl_GetVar (interp, "::stloptions.yangle", 0));
    stlparam.contyangle =
      atof (Tcl_GetVar (interp, "::stloptions.contyangle", 0));
    stlparam.edgecornerangle =
      atof (Tcl_GetVar (interp, "::stloptions.edgecornerangle", 0));
    stlparam.chartangle =
      atof (Tcl_GetVar (interp, "::stloptions.chartangle", 0));
    stlparam.outerchartangle =
      atof (Tcl_GetVar (interp, "::stloptions.outerchartangle", 0));

    stlparam.usesearchtree =
      atoi (Tcl_GetVar (interp, "::stloptions.usesearchtree", 0));


    stlparam.atlasminh =
      atof (Tcl_GetVar (interp, "::stloptions.atlasminh", 0));

    stlparam.resthsurfcurvfac =
      atof (Tcl_GetVar (interp, "::stloptions.resthsurfcurvfac", 0));
    stlparam.resthsurfcurvenable =
      atoi (Tcl_GetVar (interp, "::stloptions.resthsurfcurvenable", 0));

    stlparam.resthatlasfac =
      atof (Tcl_GetVar (interp, "::stloptions.resthatlasfac", 0));
    stlparam.resthatlasenable =
      atoi (Tcl_GetVar (interp, "::stloptions.resthatlasenable", 0));

    stlparam.resthchartdistfac =
      atof (Tcl_GetVar (interp, "::stloptions.resthchartdistfac", 0));
    stlparam.resthchartdistenable =
      atoi (Tcl_GetVar (interp, "::stloptions.resthchartdistenable", 0));

    stlparam.resthlinelengthfac =
      atof (Tcl_GetVar (interp, "::stloptions.resthlinelengthfac", 0));
    stlparam.resthlinelengthenable =
      atoi (Tcl_GetVar (interp, "::stloptions.resthlinelengthenable", 0));

    stlparam.resthcloseedgefac =
      atof (Tcl_GetVar (interp, "::stloptions.resthcloseedgefac", 0));
    stlparam.resthcloseedgeenable =
      atoi (Tcl_GetVar (interp, "::stloptions.resthcloseedgeenable", 0));

    stlparam.resthedgeanglefac =
      atof (Tcl_GetVar (interp, "::stloptions.resthedgeanglefac", 0));
    stlparam.resthedgeangleenable =
      atoi (Tcl_GetVar (interp, "::stloptions.resthedgeangleenable", 0));

    stlparam.resthsurfmeshcurvfac =
      atof (Tcl_GetVar (interp, "::stloptions.resthsurfmeshcurvfac", 0));
    stlparam.resthsurfmeshcurvenable =
      atoi (Tcl_GetVar (interp, "::stloptions.resthsurfmeshcurvenable", 0));

    stlparam.recalc_h_opt =
      atoi (Tcl_GetVar (interp, "::stloptions.recalchopt", 0));
    //  stlparam.Print (cout);
    return TCL_OK;
  }



#ifdef OCCGEOMETRY
  int Ng_SetOCCParameters  (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[])
  {
    occparam.resthcloseedgefac =
      atof (Tcl_GetVar (interp, "::stloptions.resthcloseedgefac", 0));

    occparam.resthcloseedgeenable =
      atoi (Tcl_GetVar (interp, "::stloptions.resthcloseedgeenable", 0));

    return TCL_OK;
  }
#endif // OCCGEOMETRY



  int Ng_GetCommandLineParameter  (ClientData clientData,
				   Tcl_Interp * interp,
				   int argc, tcl_const char *argv[])
  {
    if (argc != 2)
      {
	Tcl_SetResult (interp, (char*)"Ng_GetCommandLineParameter needs 1 parameter",
                       TCL_STATIC);
	return TCL_ERROR;
      }

    static char buf[10];

    if (parameters.StringFlagDefined (argv[1]))
      Tcl_SetResult (interp,
		     (char*)parameters.GetStringFlag (argv[1], NULL), TCL_STATIC);
    else if (parameters.NumFlagDefined (argv[1]))
      {
	sprintf (buf, "%lf", parameters.GetNumFlag (argv[1], 0));
	Tcl_SetResult (interp, buf, TCL_STATIC);
      }
    else if (parameters.GetDefineFlag (argv[1]))
      Tcl_SetResult (interp, (char*)"defined", TCL_STATIC);
    else
      Tcl_SetResult (interp, (char*)"undefined", TCL_STATIC);

    return TCL_OK;
  }


  static int perfstepsstart;
  static int perfstepsend;
  static char* optstring = NULL;
  static char* optstringcsg = NULL;

  void * MeshingDummy (void *)
  {

    const char * savetask = multithread.task;
    multithread.task = "Generate Mesh";

    ResetTime();


    try
      {

#ifdef LOG_STREAM
	(*logout) << "Start meshing" << endl;
	(*logout) << "Meshing parameters:" << endl;
	mparam.Print (*logout);
#endif

#ifdef ACIS
	if (acisgeometry)
	  {
	    ACISGenerateMesh(*acisgeometry, mesh.Ptr(), perfstepsstart, perfstepsend, optstring);
	  }
	else
#endif
	  if (geometry2d)
	    {
	      extern void MeshFromSpline2D (SplineGeometry2d & geometry2d,
					    Mesh *& mesh, MeshingParameters & mp);
	      MeshFromSpline2D (*geometry2d, mesh.Ptr(), mparam);
	    }
	  else
	    {
	      int res = ng_geometry -> GenerateMesh (mesh.Ptr(), perfstepsstart, perfstepsend, optstringcsg);
	      if (res != MESHING3_OK) 
		{
		  multithread.task = savetask;
		  multithread.running = 0;
		  return 0;
		}
	    }

	if (mparam.autozrefine && ( (NetgenGeometry*)geometry.Ptr() == ng_geometry))
	  {
	    ZRefinementOptions opt;
	    opt.minref = 5;
	    ZRefinement (*mesh, geometry.Ptr(), opt);
	    mesh -> SetNextMajorTimeStamp();
	  }
	
	if (mparam.secondorder)
	  {
	    const_cast<Refinement&> (ng_geometry -> GetRefinement()).MakeSecondOrder (*mesh);
	    mesh -> SetNextMajorTimeStamp();
	  }

	if (mparam.elementorder > 1)
	  {
	    mesh -> GetCurvedElements().BuildCurvedElements (&const_cast<Refinement&> (ng_geometry -> GetRefinement()),
							     mparam.elementorder);

	    mesh -> SetNextMajorTimeStamp();
	  }


	PrintMessage (1, "Meshing done, time = ", GetTime(), " sec");
      }

    catch (NgException e)
      {
	cout << e.What() << endl;
      }

    multithread.task = savetask;
    multithread.running = 0;

#ifdef OCCGEOMETRY
    if (occgeometry)
      {
	if (occgeometry->ErrorInSurfaceMeshing())
	  {
	    char script[] = "rebuildoccdialog";
	    // int errcode = 
	      Tcl_GlobalEval (tcl_interp, script);
	  }
      }
#endif


    return NULL;
  }


  int MeshingVal(tcl_const char* str)
  {
    if (strcmp(str, "ag") == 0) {return MESHCONST_ANALYSE;}
    if (strcmp(str, "me") == 0) {return MESHCONST_MESHEDGES;}
    if (strcmp(str, "ms") == 0) {return MESHCONST_MESHSURFACE;}
    if (strcmp(str, "os") == 0) {return MESHCONST_OPTSURFACE;}
    if (strcmp(str, "mv") == 0) {return MESHCONST_MESHVOLUME;}
    if (strcmp(str, "ov") == 0) {return MESHCONST_OPTVOLUME;}

    cout << "TCL TK ERROR, wrong meshing value, return='" << str << "'" << endl;
    return 0;
  }

  int Ng_GenerateMesh  (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    multithread.running = 1;
    multithread.terminate = 0;

    Ng_SetSTLParameters (clientData, interp, 0, argv);

#ifdef OCCGEOMETRY
    Ng_SetOCCParameters (clientData, interp, 0, argv);
#endif // OCCGEOMETRY

    Ng_SetMeshingParameters (clientData, interp, 0, argv);

    perfstepsstart = 1;
    perfstepsend = 6;

    if (optstringcsg) delete optstringcsg;
    optstringcsg = NULL;
    if (optstring) delete optstring;
    optstring = NULL;

    if (argc == 2)
      {
	perfstepsstart = 1;
	perfstepsend = MeshingVal(argv[1]);
      }
    else if (argc == 3)
      {
	perfstepsstart = MeshingVal(argv[1]);
	perfstepsend = MeshingVal(argv[2]);
      }
    else if (argc == 4)
      {
	perfstepsstart = MeshingVal(argv[1]);
	perfstepsend = MeshingVal(argv[2]);
	optstring = new char[strlen(argv[3])+1];
	strcpy(optstring, argv[3]);
	optstringcsg = new char[strlen(argv[3])+1];
	strcpy(optstringcsg, argv[3]);
      }

    RunParallel (MeshingDummy, NULL);

    return TCL_OK;
  }


  int Ng_StopMeshing  (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    multithread.terminate = 1;
    return TCL_OK;
  }

  int Ng_MeshInfo  (ClientData clientData,
		    Tcl_Interp * interp,
		    int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    ostringstream str;

    if (argc >= 2 && strcmp (argv[1], "dim") == 0)
      str << mesh->GetDimension();
    else if (argc >= 2 && strcmp (argv[1], "np") == 0)
      str << mesh->GetNP();
    else if (argc >= 2 && strcmp (argv[1], "ne") == 0)
      str << mesh->GetNE();
    else if (argc >= 2 && strcmp (argv[1], "nse") == 0)
      str << mesh->GetNSE();
    else if (argc >= 2 && strcmp (argv[1], "nseg") == 0)
      str << mesh->GetNSeg();
    else if (argc >= 2 && strcmp (argv[1], "bbox") == 0)
      {
	Point3d pmin, pmax;
	mesh->GetBox (pmin, pmax);
	str << pmin.X() << " " << pmax.X() << " "
	    << pmin.Y() << " " << pmax.Y() << " "
	    << pmin.Z() << " " << pmax.Z() << endl;
      }
    else
      {
	cout << "argv[1] = " << argv[1] << endl;
	Tcl_SetResult (interp, (char*)"Ng_MeshInfo requires an argument out of \n dim np ne", TCL_STATIC);
	return TCL_ERROR;
      }

    Tcl_SetResult  (interp, (char*)str.str().c_str(), TCL_VOLATILE);
    return TCL_OK;
  }

  int Ng_MeshQuality  (ClientData clientData,
		       Tcl_Interp * interp,
		       int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    double angles[4];
    char buf[10];

    if (mesh.Ptr())
      mesh->CalcMinMaxAngle(mparam.badellimit, angles);
    else
      {
	angles[0] = angles[1] = angles[2] = angles[3] = 0;
      }
    sprintf (buf, "%5.1lf", angles[0]);
    Tcl_SetVar (interp, argv[1], buf, 0);
    sprintf (buf, "%5.1lf", angles[1]);
    Tcl_SetVar (interp, argv[2], buf, 0);
    sprintf (buf, "%5.1lf", angles[2]);
    Tcl_SetVar (interp, argv[3], buf, 0);
    sprintf (buf, "%5.1lf", angles[3]);
    Tcl_SetVar (interp, argv[4], buf, 0);

    return TCL_OK;
  }

  int Ng_CheckSurfaceMesh  (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    mesh->FindOpenElements();
    if (mesh->CheckConsistentBoundary())
      {
	PrintMessage (1, "surface mesh not consistent, trying orientation");
	mesh->SurfaceMeshOrientation();
      }
    else
      {
	PrintMessage (1, "surface mesh consistent");
      }

    mesh->CheckOverlappingBoundary();
    return TCL_OK;
  }

  int Ng_CheckVolumeMesh  (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    mesh->CheckVolumeMesh();
    return TCL_OK;
  }


  int Ng_DeleteVolMesh  (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[])
  {
    if (mesh.Ptr())
      mesh->ClearVolumeElements();

    return TCL_OK;
  }


  int Ng_SplitSeparatedFaces (ClientData clientData,
			      Tcl_Interp * interp,
			      int argc, tcl_const char *argv[])
  {
    if (mesh.Ptr())
      mesh->SplitSeparatedFaces ();
    return TCL_OK;
  }



  int Ng_RestrictH  (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    if (argc != 3)
      return TCL_OK;
    if (!mesh.Ptr())
      return TCL_OK;


    double loch = atof (argv[2]);
    if (strcmp (argv[1], "face") == 0)
      {
	cout << "Restrict h at face to " << loch << endl;
	mesh -> RestrictLocalH  (RESTRICTH_FACE, vsmesh.SelectedFace(), loch);
      }
    if (strcmp (argv[1], "edge") == 0)
      {
	cout << "Restrict h at edge to " << loch << endl;
	mesh -> RestrictLocalH  (RESTRICTH_EDGE, vsmesh.SelectedEdge(), loch);
      }
    if (strcmp (argv[1], "point") == 0)
      {
	cout << "Restrict h at point to " << loch << endl;
	mesh -> RestrictLocalH  (RESTRICTH_POINT, vsmesh.SelectedPoint(), loch);
      }

    return TCL_OK;
  }





  int Ng_Anisotropy  (ClientData clientData,
		      Tcl_Interp * interp,
		      int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    if (argc != 2)
      return TCL_OK;
    if (!mesh.Ptr())
      return TCL_OK;

    if (strcmp (argv[1], "edge") == 0)
      {
	int edgenr = vsmesh.SelectedEdge();
	for (int i = 1; i <= mesh->GetNSeg(); i++)
	  {
	    Segment & seg = mesh->LineSegment(i);
	    if (seg.edgenr == edgenr)
	      {
		seg.singedge_left = 1 - seg.singedge_left;
		seg.singedge_right = 1 - seg.singedge_right;
	      }
	  }
      }

    return TCL_OK;
  }




  BisectionOptions biopt;

  void * BisectDummy (void *)
  {
    Refinement * ref;
    MeshOptimize2d * opt = NULL;
    if (stlgeometry)
      ref = new RefinementSTLGeometry(*stlgeometry);
#ifdef OCCGEOMETRY
    else if (occgeometry)
      ref = new OCCRefinementSurfaces (*occgeometry);
#endif
#ifdef ACIS
    else if (acisgeometry)
      {
	ref = new ACISRefinementSurfaces(*acisgeometry);
	opt = new ACISMeshOptimize2dSurfaces(*acisgeometry);
	ref->Set2dOptimizer(opt);
      }
#endif
    else
      {
	ref = new RefinementSurfaces(*geometry);
	opt = new MeshOptimize2dSurfaces(*geometry);
	ref->Set2dOptimizer(opt);
      }

    if(!mesh->LocalHFunctionGenerated())
      mesh->CalcLocalH();

    mesh->LocalHFunction().SetGrading (mparam.grading);
    ref -> Bisect (*mesh, biopt);
    mesh -> UpdateTopology();
    mesh -> GetCurvedElements().BuildCurvedElements (ref, mparam.elementorder);

    multithread.running = 0;
    delete ref;
    delete opt;
    return NULL;
  }


  int Ng_Bisect  (ClientData clientData,
		  Tcl_Interp * interp,
		  int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }


    if (multithread.running)
      {
	cout << "Thread alrad running" << endl;
	return TCL_OK;
      }
    multithread.running = 1;

    biopt.outfilename = NULL; // "ngfepp.vol";
    biopt.femcode = "fepp";
    biopt.refinementfilename = NULL;

    if (argc >= 2)
      biopt.refinementfilename = argv[1];


    //  pthread_create (&meshingthread, NULL, &BisectDummy, NULL);
    BisectDummy (0);

    /*
      extern void BisectTets (Mesh &, const CSGeometry *);
      BisectTets (*mesh, geometry);
    */
    return TCL_OK;
  }



  //   int Ng_BisectCopyMesh  (ClientData clientData,
  // 			  Tcl_Interp * interp,
  // 			  int argc, tcl_const char *argv[])
  //   {
  //     if (!mesh.Ptr())
  //       {
  // 	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
  // 	return TCL_ERROR;
  //       }
  //     if (multithread.running)
  //       {
  // 	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
  // 	return TCL_ERROR;
  //       }

  //     BisectTetsCopyMesh (*mesh, geometry.Ptr(), biopt);
  //     return TCL_OK;
  //   }






  int Ng_Split2Tets  (ClientData clientData,
		      Tcl_Interp * interp,
		      int argc, tcl_const char *argv[])
  {
    if (!mesh.Ptr())
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }
    if (multithread.running)
      {
	Tcl_SetResult (interp, err_jobrunning, TCL_STATIC);
	return TCL_ERROR;
      }

    mesh->Split2Tets ();

    return TCL_OK;
  }










  int Ng_STLDoctor (ClientData clientData,
		    Tcl_Interp * interp,
		    int argc, tcl_const char *argv[])
  {
    //cout << "STL doctor" << endl;


    stldoctor.drawmeshededges =
      atoi (Tcl_GetVar (interp, "::stldoctor.drawmeshededges", 0));

    stldoctor.geom_tol_fact =
      atof (Tcl_GetVar (interp, "::stldoctor.geom_tol_fact", 0));


    stldoctor.useexternaledges =
      atoi (Tcl_GetVar (interp, "::stldoctor.useexternaledges", 0));

    stldoctor.showfaces =
      atoi (Tcl_GetVar (interp, "::stldoctor.showfaces", 0));

    stldoctor.conecheck =
      atoi (Tcl_GetVar (interp, "::stldoctor.conecheck", 0));

    stldoctor.spiralcheck =
      atoi (Tcl_GetVar (interp, "::stldoctor.spiralcheck", 0));

    stldoctor.selectwithmouse =
      atoi (Tcl_GetVar (interp, "::stldoctor.selectwithmouse", 0));

    stldoctor.showedgecornerpoints =
      atoi (Tcl_GetVar (interp, "::stldoctor.showedgecornerpoints", 0));

    stldoctor.showmarkedtrigs =
      atoi (Tcl_GetVar (interp, "::stldoctor.showmarkedtrigs", 0));

    stldoctor.showtouchedtrigchart =
      atoi (Tcl_GetVar (interp, "::stldoctor.showtouchedtrigchart", 0));

    //cout << "smt=" << stldoctor.showmarkedtrigs << endl;

    stldoctor.dirtytrigfact =
      atof (Tcl_GetVar (interp, "::stldoctor.dirtytrigfact", 0));

    stldoctor.smoothnormalsweight =
      atof (Tcl_GetVar (interp, "::stldoctor.smoothnormalsweight", 0));

    stldoctor.smoothangle =
      atof (Tcl_GetVar (interp, "::stldoctor.smoothangle", 0));

    stldoctor.selectmode =
      atoi (Tcl_GetVar (interp, "::stldoctor.selectmode", 0));

    stldoctor.edgeselectmode =
      atoi (Tcl_GetVar (interp, "::stldoctor.edgeselectmode", 0));

    stldoctor.longlinefact =
      atoi (Tcl_GetVar (interp, "::stldoctor.longlinefact", 0));

    stldoctor.showexcluded =
      atoi (Tcl_GetVar (interp, "::stldoctor.showexcluded", 0));



    if (!stldoctor.selectwithmouse)
      {
	stldoctor.selecttrig =
	  atoi (Tcl_GetVar (interp, "::stldoctor.selecttrig", 0));

	stldoctor.nodeofseltrig =
	  atoi (Tcl_GetVar (interp, "::stldoctor.nodeofseltrig", 0));
      }

    stldoctor.showvicinity =
      atoi (Tcl_GetVar (interp, "::stldoctor.showvicinity", 0));

    stldoctor.vicinity =
      atoi (Tcl_GetVar (interp, "::stldoctor.vicinity", 0));


    if (argc >= 2)
      {
	if (!stlgeometry)
	  {
	    Tcl_SetResult (interp, err_needsstlgeometry, TCL_STATIC);
	    return TCL_ERROR;
	  }

	if (strcmp (argv[1], "destroy0trigs") == 0)
	  {
	    stlgeometry->DestroyDirtyTrigs();
	  }
	else if (strcmp (argv[1], "movepointtomiddle") == 0)
	  {
	    stlgeometry->MoveSelectedPointToMiddle();
	  }
	else if (strcmp (argv[1], "calcnormals") == 0)
	  {
	    stlgeometry->CalcNormalsFromGeometry();
	  }
	else if (strcmp (argv[1], "showchartnum") == 0)
	  {
	    stlgeometry->ShowSelectedTrigChartnum();
	  }
	else if (strcmp (argv[1], "showcoords") == 0)
	  {
	    stlgeometry->ShowSelectedTrigCoords();
	  }
	else if (strcmp (argv[1], "loadmarkedtrigs") == 0)
	  {
	    stlgeometry->LoadMarkedTrigs();
	  }
	else if (strcmp (argv[1], "savemarkedtrigs") == 0)
	  {
	    stlgeometry->SaveMarkedTrigs();
	  }
	else if (strcmp (argv[1], "neighbourangles") == 0)
	  {
	    stlgeometry->NeighbourAnglesOfSelectedTrig();
	  }
	else if (strcmp (argv[1], "vicinity") == 0)
	  {
	    stlgeometry->CalcVicinity(stldoctor.selecttrig);
	  }
	else if (strcmp (argv[1], "markdirtytrigs") == 0)
	  {
	    stlgeometry->MarkDirtyTrigs();
	  }
	else if (strcmp (argv[1], "smoothdirtytrigs") == 0)
	  {
	    stlgeometry->SmoothDirtyTrigs();
	  }
	else if (strcmp (argv[1], "smoothrevertedtrigs") == 0)
	  {
	    stlgeometry->GeomSmoothRevertedTrigs();
	  }
	else if (strcmp (argv[1], "invertselectedtrig") == 0)
	  {
	    stlgeometry->InvertTrig(stlgeometry->GetSelectTrig());
	  }
	else if (strcmp (argv[1], "deleteselectedtrig") == 0)
	  {
	    stlgeometry->DeleteTrig(stlgeometry->GetSelectTrig());
	  }
	else if (strcmp (argv[1], "smoothgeometry") == 0)
	  {
	    stlgeometry->SmoothGeometry();
	  }
	else if (strcmp (argv[1], "orientafterselectedtrig") == 0)
	  {
	    stlgeometry->OrientAfterTrig(stlgeometry->GetSelectTrig());
	  }
	else if (strcmp (argv[1], "marktoperrortrigs") == 0)
	  {
	    stlgeometry->MarkTopErrorTrigs();
	  }
	else if (strcmp (argv[1], "exportedges") == 0)
	  {
	    stlgeometry->ExportEdges();
	  }
	else if (strcmp (argv[1], "importedges") == 0)
	  {
	    stlgeometry->ImportEdges();
	  }
	else if (strcmp (argv[1], "importexternaledges") == 0)
	  {
	    stlgeometry->ImportExternalEdges(argv[2]);
	  }
	else if (strcmp (argv[1], "loadedgedata") == 0)
	  {
	    if (argc >= 3)
	      {
		stlgeometry->LoadEdgeData(argv[2]);
	      }
	  }
	else if (strcmp (argv[1], "saveedgedata") == 0)
	  {
	    if (argc >= 3)
	      {
		stlgeometry->SaveEdgeData(argv[2]);
	      }
	  }

	else if (strcmp (argv[1], "buildexternaledges") == 0)
	  {
	    stlgeometry->BuildExternalEdgesFromEdges();
	  }
	else if (strcmp (argv[1], "smoothnormals") == 0)
	  {
	    stlgeometry->SmoothNormals();
	  }
	else if (strcmp (argv[1], "marknonsmoothnormals") == 0)
	  {
	    stlgeometry->MarkNonSmoothNormals();
	  }
	else if (strcmp (argv[1], "addexternaledge") == 0)
	  {
	    stlgeometry->AddExternalEdgeAtSelected();
	  }
	else if (strcmp (argv[1], "addgeomline") == 0)
	  {
	    stlgeometry->AddExternalEdgesFromGeomLine();
	  }
	else if (strcmp (argv[1], "addlonglines") == 0)
	  {
	    stlgeometry->AddLongLinesToExternalEdges();
	  }
	else if (strcmp (argv[1], "addclosedlines") == 0)
	  {
	    stlgeometry->AddClosedLinesToExternalEdges();
	  }
	else if (strcmp (argv[1], "addnotsinglelines") == 0)
	  {
	    stlgeometry->AddAllNotSingleLinesToExternalEdges();
	  }
	else if (strcmp (argv[1], "deletedirtyexternaledges") == 0)
	  {
	    stlgeometry->DeleteDirtyExternalEdges();
	  }
	else if (strcmp (argv[1], "deleteexternaledge") == 0)
	  {
	    stlgeometry->DeleteExternalEdgeAtSelected();
	  }
	else if (strcmp (argv[1], "deletevicexternaledge") == 0)
	  {
	    stlgeometry->DeleteExternalEdgeInVicinity();
	  }

	else if (strcmp (argv[1], "addlonglines") == 0)
	  {
	    stlgeometry->STLDoctorLongLinesToCandidates();
	  }
	else if (strcmp (argv[1], "deletedirtyedges") == 0)
	  {
	    stlgeometry->STLDoctorDirtyEdgesToCandidates();
	  }
	else if (strcmp (argv[1], "undoedgechange") == 0)
	  {
	    stlgeometry->UndoEdgeChange();
	  }
	else if (strcmp (argv[1], "buildedges") == 0)
	  {
	    stlgeometry->STLDoctorBuildEdges();
	  }
	else if (strcmp (argv[1], "confirmedge") == 0)
	  {
	    stlgeometry->STLDoctorConfirmEdge();
	  }
	else if (strcmp (argv[1], "candidateedge") == 0)
	  {
	    stlgeometry->STLDoctorCandidateEdge();
	  }
	else if (strcmp (argv[1], "excludeedge") == 0)
	  {
	    stlgeometry->STLDoctorExcludeEdge();
	  }
	else if (strcmp (argv[1], "undefinededge") == 0)
	  {
	    stlgeometry->STLDoctorUndefinedEdge();
	  }
	else if (strcmp (argv[1], "setallundefinededges") == 0)
	  {
	    stlgeometry->STLDoctorSetAllUndefinedEdges();
	  }
	else if (strcmp (argv[1], "erasecandidateedges") == 0)
	  {
	    stlgeometry->STLDoctorEraseCandidateEdges();
	  }
	else if (strcmp (argv[1], "confirmcandidateedges") == 0)
	  {
	    stlgeometry->STLDoctorConfirmCandidateEdges();
	  }
	else if (strcmp (argv[1], "confirmedtocandidateedges") == 0)
	  {
	    stlgeometry->STLDoctorConfirmedToCandidateEdges();
	  }
      }

    return TCL_OK;
  }





  extern int Ng_MeshDoctor (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[]);


  int Ng_STLInfo  (ClientData clientData,
		   Tcl_Interp * interp,
		   int argc, tcl_const char *argv[])
  {
    double data[10];
    static char buf[20];

    if (!stlgeometry)
      {
	Tcl_SetResult (interp, err_needsstlgeometry, TCL_STATIC);
	return TCL_ERROR;
      }



    if (stlgeometry)
      {
	stlgeometry->STLInfo(data);
	//      cout << "NT=" << data[0] << endl;

	if (argc == 2)
	  {
	    if (strcmp (argv[1], "status") == 0)
	      {
		switch (stlgeometry->GetStatus())
		  {
		  case STLGeometry::STL_GOOD:
		    strcpy (buf, "GOOD"); break;
		  case STLGeometry::STL_WARNING:
		    strcpy (buf, "WARNING"); break;
		  case STLGeometry::STL_ERROR:
		    strcpy (buf, "ERROR"); break;
		  }
		Tcl_SetResult (interp, buf, TCL_STATIC);
		return TCL_OK;
	      }
	    if (strcmp (argv[1], "statustext") == 0)
	      {
		Tcl_SetResult (interp, (char*)stlgeometry->GetStatusText().c_str(), TCL_STATIC);
		return TCL_OK;
	      }
	    if (strcmp (argv[1], "topology_ok") == 0)
	      {
		sprintf (buf, "%d", stlgeometry->Topology_Ok());
		Tcl_SetResult (interp, buf, TCL_STATIC);
	      }
	    if (strcmp (argv[1], "orientation_ok") == 0)
	      {
		sprintf (buf, "%d", stlgeometry->Orientation_Ok());
		Tcl_SetResult (interp, buf, TCL_STATIC);
	      }
	  }
      }
    else
      {
	data[0] = 0;
	data[1] = 0;
	data[2] = 0;
	data[3] = 0;
	data[4] = 0;
	data[5] = 0;
	data[6] = 0;
	data[7] = 0;
      }




    sprintf (buf, "%i", (int)data[0]);
    Tcl_SetVar (interp, argv[1], buf, 0);

    sprintf (buf, "%5.3g", data[1]);
    Tcl_SetVar (interp, argv[2], buf, 0);
    sprintf (buf, "%5.3g", data[2]);
    Tcl_SetVar (interp, argv[3], buf, 0);
    sprintf (buf, "%5.3g", data[3]);
    Tcl_SetVar (interp, argv[4], buf, 0);

    sprintf (buf, "%5.3g", data[4]);
    Tcl_SetVar (interp, argv[5], buf, 0);
    sprintf (buf, "%5.3g", data[5]);
    Tcl_SetVar (interp, argv[6], buf, 0);
    sprintf (buf, "%5.3g", data[6]);
    Tcl_SetVar (interp, argv[7], buf, 0);

    sprintf (buf, "%i", (int)data[7]);
    Tcl_SetVar (interp, argv[8], buf, 0);

    return TCL_OK;
  }


  int Ng_STLCalcLocalH  (ClientData clientData,
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[])
  {

    Ng_SetSTLParameters (clientData, interp, argc, argv);

#ifdef OCCGEOMETRY
    Ng_SetOCCParameters (clientData, interp, argc, argv);
#endif // OCCGEOMETRY

    Ng_SetMeshingParameters (clientData, interp, argc, argv);

    if (mesh.Ptr() && stlgeometry)
      {
	mesh -> SetLocalH (stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
			   stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
			   mparam.grading);
	stlgeometry -> RestrictLocalH(*mesh, mparam.maxh);

	if (stlparam.resthsurfmeshcurvenable)
	  mesh -> CalcLocalHFromSurfaceCurvature (stlparam.resthsurfmeshcurvfac);
      }

    return TCL_OK;
  }




  SYMBOLTABLE<VisualScene*> & GetVisualizationScenes ()
  {
    static SYMBOLTABLE<VisualScene*> vss;
    return vss;
  }

  void AddVisualizationScene (const string & name,
			      VisualScene * avs)
  {
    GetVisualizationScenes().Set (name.c_str(), avs);
  }


  void SetVisualScene (Tcl_Interp * interp)
  {
    const char * vismode = vispar.selectvisual;
    // Tcl_GetVar (interp, "selectvisual", 0);
    vs = &vscross;
    if (GetVisualizationScenes().Used(vismode))
      {
	vs = GetVisualizationScenes().Get(vismode);
      }
    else if (vismode)
      {
	if (strcmp (vismode, "geometry") == 0)
	  {
	    if (stlgeometry != NULL)
	      vs = &vsstlmeshing;
	    else if (geometry2d)
	      vs = &vsgeom2d;

#ifdef OCCGEOMETRY
	    else if (occgeometry)
	      vs = &vsoccgeom;
#endif // OCCGEOMETRY
#ifdef ACIS
	    else if (acisgeometry)
	      vs = &vsacisgeom;
#endif // ACIS
	    else
	      vs = &vsgeom;

	    // vs = &vsstlgeom;
	  }
	if (strcmp (vismode, "mesh") == 0)
	  {
	    if (!meshdoctor.active)
	      vs = &vsmesh;
	    else
	      vs = &vsmeshdoc;
	  }
          // if (strcmp (vismode, "surfmeshing") == 0) vs = &vssurfacemeshing;
	if (strcmp (vismode, "specpoints") == 0) vs = &vsspecpoints;
	//      if (strcmp (vismode, "solution") == 0) vs = &vssolution;
      }
  }






#if TOGL_MAJOR_VERSION==1


  // Togl
  
  static int fontbase = 0;

  void MyOpenGLText (const char * text)
  {
    glListBase (fontbase);
    glCallLists (GLsizei(strlen(text)), GL_UNSIGNED_BYTE, text);
  }



  static void init( struct Togl *togl )
  {
    fontbase = Togl_LoadBitmapFont( togl, TOGL_BITMAP_8_BY_13 );

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);

    SetVisualScene (Togl_Interp(togl));
    vs->DrawScene();
  }

  static void zap( struct Togl *togl )
  {
    ;
  }

  static void draw( struct Togl *togl )
  {
    Tcl_Interp * interp = Togl_Interp(togl);

    SetVisualScene (interp);

#ifdef STEREO
    if (1) // vispar.stereo)
      {
	glMatrixMode (GL_MODELVIEW);
	glPushMatrix();

	glLoadIdentity ();

	//  glTranslatef (0.1, 0, 0);
	gluLookAt (0.3, 0, 6, 0, 0, 0, 0, 1, 0);
	Togl_StereoDrawBuffer(GL_BACK_RIGHT);
	vs->DrawScene();

	glLoadIdentity ();
	// glTranslatef (-0.1, 0, 0);
	gluLookAt (-0.3, 0, 6, 0, 0, 0, 0, 1, 0);
	Togl_StereoDrawBuffer(GL_BACK_LEFT);
	vs->DrawScene();
	glPopMatrix();

	Togl_SwapBuffers(togl);
      }
    else
#endif
      {
	glPushMatrix();
	glLoadIdentity();
	// gluLookAt (0, 0, 6, 0, 0, 0, 0, 1, 0);
	vs->DrawScene();

	Togl_SwapBuffers(togl);
	glPopMatrix();
      }
  }

  static void reshape( struct Togl *togl)
  {
    int w = Togl_Width (togl);
    int h = Togl_Height (togl);

    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(20.0f, double(w) / h, pnear, pfar);
    glMatrixMode(GL_MODELVIEW);

    draw (togl);
  }



#else


  // Sorry, Togl 2.0 not supported

  Tcl_Obj * togl_font;
 
  void MyOpenGLText (const char * text)
  {
    cout << "togl - text" << endl;
  }
  


  static int
  init(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
  {
    cout << "call init" << endl;

    Togl * togl = NULL;

    if (Togl_GetToglFromObj(interp, objv[1], &togl) != TCL_OK) 
      return TCL_ERROR;

    cout << "call Togl - load font (crash on my Linux64)" << endl;
    // togl_font = Togl_LoadBitmapFont( togl, "Times"); // TOGL_BITMAP_8_BY_13 );
    // togl_font = Togl_LoadBitmapFont( togl, TOGL_BITMAP_8_BY_13 );
    // togl_font = Togl_LoadBitmapFont( togl, NULL );
    // cout << "success" << endl;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);

    SetVisualScene (Togl_Interp(togl));
    vs->DrawScene();
    return TCL_OK;
  }

  static int
  zap(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
  {
    return TCL_OK;
  }


  static int
  draw(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
  {
    Togl * togl;
    if (Togl_GetToglFromObj(interp, objv[1], &togl) != TCL_OK) 
      return TCL_ERROR;
    

    SetVisualScene (interp);

    glPushMatrix();
    glLoadIdentity();
    // gluLookAt (0, 0, 6, 0, 0, 0, 0, 1, 0);
    vs->DrawScene();
    
    Togl_SwapBuffers(togl);
    glPopMatrix();

    return TCL_OK;
  }

  static int
  reshape(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
  {
    Togl * togl;
    if (Togl_GetToglFromObj(interp, objv[1], &togl) != TCL_OK) 
      return TCL_ERROR;

    int w = Togl_Width (togl);
    int h = Togl_Height (togl);

    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(20.0f, double(w) / h, pnear, pfar);
    glMatrixMode(GL_MODELVIEW);

    // draw (togl);

    return TCL_OK;
  }





#endif





#if TOGL_MAJOR_VERSION==1

#ifndef JPEGLIB
  static int Ng_SnapShot (struct Togl * togl,
                          int argc, tcl_const char *argv[])
  {
    const char * filename = argv[2];

    char str[250];
    char filename2[250];
    int len = strlen(filename);
    strcpy (filename2, filename);

    filename2[len-3] = 'p';
    filename2[len-2] = 'p';
    filename2[len-1] = 'm';
    filename2[len] = 0;

    cout << "Snapshot to file '" << filename << endl;

    int w = Togl_Width (togl);
    w = int((w + 1) / 4) * 4 + 4;
    int h = Togl_Height (togl);

    // unsigned char * buffer = new unsigned char[w*h*4];
    unsigned char * buffer = new unsigned char[w*h*3];
    glReadPixels (0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, buffer);

    ofstream outfile(filename2);
    outfile << "P6" << endl
    << "# CREATOR: Netgen" << endl
    << w << " " << h << endl
    << "255" << endl;
    for (int i = 0; i < h; i++)
    for (int j = 0; j < w; j++)
    for (int k = 0; k < 3; k++)
    outfile.put (buffer[k+3*j+3*w*(h-i-1)]);
    outfile << flush;

    delete[] buffer;

    // convert image file (Unix/Linux only):
    sprintf(str,"convert -quality 100 %s %s", filename2, filename);

    int err = system(str);
    if (err != 0)
    {
      Tcl_SetResult (Togl_Interp(togl), (char*)"Cannot convert image file", TCL_VOLATILE);
      return TCL_ERROR;
    }
    sprintf(str,"rm %s", filename2);
    system(str);

    return TCL_OK;
  }



#else


  static int Ng_SnapShot (struct Togl * togl,
			  int argc, tcl_const char *argv[])
  {
    const char * filename = argv[2];
    int len = strlen(filename);

    if (strcmp ("jpg", filename+len-3) == 0)
      {
        cout << "Snapshot to file '" << filename << "'" << endl;

        int w = Togl_Width (togl);
        w = int((w + 1) / 4) * 4 + 4;
        int h = Togl_Height (togl);

        // unsigned char * buffer = new unsigned char[w*h*4];
        unsigned char * buffer = new unsigned char[w*h*3];
        glReadPixels (0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, buffer);

        struct jpeg_compress_struct cinfo;
        struct jpeg_error_mgr jerr;
        FILE *outfile = fopen(filename,"wb");
        JSAMPROW row_pointer[1];
        int row_stride, quality = 85; // 1...100

        cinfo.err = jpeg_std_error( &jerr );
        jpeg_create_compress( &cinfo );
        jpeg_stdio_dest( &cinfo, outfile );


        cinfo.image_width = w;
        cinfo.image_height = h;
        cinfo.input_components = 3;
        cinfo.in_color_space = JCS_RGB;

        jpeg_set_defaults( &cinfo );
        jpeg_set_quality( &cinfo, quality, TRUE );
        jpeg_start_compress( &cinfo, TRUE );

        row_stride = 3*w;
        while( cinfo.next_scanline < cinfo.image_height ) {
          row_pointer[0] = &buffer[ (h-1-cinfo.next_scanline) * row_stride ];
          (void)jpeg_write_scanlines( &cinfo, row_pointer, 1 );
        }

        jpeg_finish_compress( &cinfo );
        fclose( outfile );

        jpeg_destroy_compress( &cinfo );
        fprintf( stdout, "done [ok]\n" );
        fflush( stdout );

        free( buffer );
        return TCL_OK;
      }
    else
      {
        cout << "Snapshot to " << filename << " not supported" << endl;
        return TCL_ERROR;
      }
  }


#endif





#ifdef FFMPEG

  // thanks to Mikko Lyly @ CSC, Helsinki
  
#define STATE_READY 0
#define STATE_STARTED 1


#define INBUF_SIZE 4096
#define DEFAULT_B_FRAMES 3
  // #define DEFAULT_B_FRAMES 0
#define DEFAULT_GOP_SIZE 200
  // #define DEFAULT_GOP_SIZE 10
  // #define DEFAULT_BITRATE 500000
#define DEFAULT_BITRATE 5000000
  // #define DEFAULT_MPG_BUFSIZE 500000
#define DEFAULT_MPG_BUFSIZE 500000

  typedef struct buffer_s {
    uint8_t *MPG;
    uint8_t *YUV;
    uint8_t *RGB;
    uint8_t *ROW;
  } buffer_t;

  void free_buffers( buffer_t *buff ) {
    free( buff->MPG );
    free( buff->YUV );
    free( buff->RGB );
    free( buff->ROW );
  }

  static double psnr( double d ) {
    if( d==0 )
      return INFINITY;
    return -10.0*log( d )/log( 10.0 );
  }

  void print_info( int count_frames, AVCodecContext *context, int bytes ) {
    double tmp = context->width * context->height * 255.0 * 255.0;
    double Ypsnr = psnr( context->coded_frame->error[0] / tmp );
    double quality = context->coded_frame->quality/(double)FF_QP2LAMBDA;
    char pict_type = av_get_pict_type_char(context->coded_frame->pict_type);
    cout << "video: frame=" << count_frames << " type=" << pict_type;
    cout << " size=" << bytes << " PSNR(Y)=" << Ypsnr << " dB q=" << (float)quality << endl;
  }



  static int Ng_VideoClip (struct Togl * togl,
                           int argc, tcl_const char *argv[])
  {
    static AVCodec *codec = NULL;
    static AVCodecContext *context = NULL;
    static AVFrame *YUVpicture = NULL;
    static AVFrame *RGBpicture = NULL;
    static int bytes, PIXsize, stride;
    static int y, nx, ny, ox, oy, viewp[4];
    static int i_state = STATE_READY;
    static int initialized = 0;
    static int count_frames = 0;
    static int bitrate = DEFAULT_BITRATE;
    static int gopsize = DEFAULT_GOP_SIZE;
    static int bframes = DEFAULT_B_FRAMES;
    static int MPGbufsize = DEFAULT_MPG_BUFSIZE;
    static CodecID codec_id = CODEC_ID_MPEG1VIDEO;
    static FILE *MPGfile;
    static buffer_t buff;
    static struct SwsContext *img_convert_ctx;


    if (strcmp (argv[2], "init") == 0)
      {
        // Can't initialize when running:
        //-------------------------------
        if( i_state != STATE_READY ) {
          cout << "cannot initialize: already running" << endl;
          return TCL_ERROR;
        }



        // Open output file:
        //-------------------
        const char * filename = argv[3];
        cout << "Saving videoclip to file '" << filename << "'" << endl;
        MPGfile = fopen(filename, "wb");

        // Determine picture size:
        //------------------------
        nx = Togl_Width (togl);
        nx = int((nx + 1) / 4) * 4 + 4;
        ny = Togl_Height (togl);
        ny = 2 * (ny/2);
        cout << "Width=" << nx << ", height=" << ny << endl;

        // Allocate buffers:
        //------------------
        PIXsize = nx*ny;
        stride = 3*nx;
        buff.RGB = (uint8_t*)malloc(stride*ny);
        buff.ROW = (uint8_t*)malloc(stride);
        buff.YUV = (uint8_t*)malloc(3*(PIXsize/2));
        buff.MPG = (uint8_t*)malloc(MPGbufsize);


        // Initialize libavcodec:
        //-----------------------
        if( !initialized ) {
          av_register_all();
          initialized = 1;
        }

        // Choose codec:
        //--------------
        codec = avcodec_find_encoder( codec_id );
        if( !codec ) {
          free_buffers( &buff );
          fclose( MPGfile );
          cout << "can't find codec" << endl;
          return TCL_ERROR;
        }

        // Init codec context etc.:
        //--------------------------
        context = avcodec_alloc_context();

        context->bit_rate = bitrate;
        context->width = nx;
        context->height = ny;
	context->time_base = (AVRational){ 1, 25 };
        context->gop_size = gopsize;
        context->max_b_frames = bframes;
        context->pix_fmt = PIX_FMT_YUV420P;
        context->flags |= CODEC_FLAG_PSNR;

        if( avcodec_open( context, codec ) < 0 ) {
          cout << "can't open codec" << endl;
          avcodec_close( context );
          av_free( context );
          free_buffers( &buff );
          fclose( MPGfile );
          return TCL_ERROR;
        }

        YUVpicture = avcodec_alloc_frame();

        YUVpicture->data[0] = buff.YUV;
        YUVpicture->data[1] = buff.YUV + PIXsize;
        YUVpicture->data[2] = buff.YUV + PIXsize + PIXsize / 4;
        YUVpicture->linesize[0] = nx;
        YUVpicture->linesize[1] = nx / 2;
        YUVpicture->linesize[2] = nx / 2;

        RGBpicture = avcodec_alloc_frame();

        RGBpicture->data[0] = buff.RGB;
        RGBpicture->data[1] = buff.RGB;
        RGBpicture->data[2] = buff.RGB;
        RGBpicture->linesize[0] = stride;
        RGBpicture->linesize[1] = stride;
        RGBpicture->linesize[2] = stride;

        // Set state "started":
        //----------------------
        i_state = STATE_STARTED;
        cout << "savempg: state: started" << endl;

        return TCL_OK;
      }



    else if (strcmp (argv[2], "addframe") == 0)
      {
        // Can't compress if status != started:
        //-------------------------------------
        if( i_state != STATE_STARTED ) {
          cout << "cannot add frame: codec not initialized" << endl;
          return TCL_ERROR;
        }

        // Read RGB data:
        //---------------
        glReadPixels (0, 0, nx, ny, GL_RGB, GL_UNSIGNED_BYTE, buff.RGB );

        // The picture is upside down - flip it:
        //---------------------------------------
        for( y=0; y<ny/2; y++ ) {
          uint8_t *r1 = buff.RGB + stride*y;
          uint8_t *r2 = buff.RGB + stride*(ny-1-y);
          memcpy( buff.ROW, r1, stride );
          memcpy( r1, r2, stride );
          memcpy( r2, buff.ROW, stride );
        }

        // Convert to YUV:
        //----------------
        if( img_convert_ctx == NULL )
          img_convert_ctx = sws_getContext( nx, ny, PIX_FMT_RGB24,
                                            nx, ny, PIX_FMT_YUV420P,
                                            SWS_BICUBIC, NULL, NULL, NULL );
        
        if( img_convert_ctx == NULL ) {
          cout << "can't initialize scaler context" << endl;
          return TCL_ERROR;
        }
        
        sws_scale( img_convert_ctx, RGBpicture->data, RGBpicture->linesize,
                   0, ny, YUVpicture->data, YUVpicture->linesize );


        // Encode frame:
        //--------------
        bytes = avcodec_encode_video( context, buff.MPG,
                                      MPGbufsize, YUVpicture );
        count_frames++;
        print_info( count_frames, context, bytes );
        fwrite( buff.MPG, 1, bytes, MPGfile );

        return TCL_OK;
      }



    else if (strcmp (argv[2], "finalize") == 0)
      {
        // Can't stop if status != started:
        //---------------------------------
        if( i_state != STATE_STARTED ) {
          cout << "cannot finalize: codec not initialized" << endl;
          return TCL_ERROR;
        }

        // Get the delayed frames, if any:
        //--------------------------------
        for( ; bytes; ) {
          bytes = avcodec_encode_video( context, buff.MPG, MPGbufsize, NULL );
          count_frames++;
          print_info( count_frames, context, bytes );
          fwrite( buff.MPG, 1, bytes, MPGfile );
        }

        // Add sequence end code:
        //-----------------------
        if( codec_id == CODEC_ID_MPEG1VIDEO ) {
          buff.MPG[0] = 0x00;
          buff.MPG[1] = 0x00;
          buff.MPG[2] = 0x01;
          buff.MPG[3] = 0xb7;
          fwrite( buff.MPG, 1, 4, MPGfile );
        }

        // Finalize:
        //-----------
        avcodec_close( context );
        av_free( context );
        av_free( YUVpicture );
        av_free( RGBpicture );
        free_buffers( &buff );
        fclose( MPGfile );

        i_state = STATE_READY;
        cout <<  "finalized" << endl;

        return TCL_OK;
      }
  }



#else
  static int Ng_VideoClip (struct Togl * togl,
                           int argc, tcl_const char *argv[])
  {
    Tcl_SetResult (Togl_Interp(togl), (char*)"Video not available, Netgen was not compiled with FFMPEG library", TCL_STATIC);
    return TCL_ERROR;
  }
#endif

#endif












  int Ng_MouseMove (ClientData clientData,
		    Tcl_Interp * interp,
		    int argc, tcl_const char *argv[])
  {
    int oldx, oldy;
    int newx, newy;

    oldx = atoi (argv[1]);
    oldy = atoi (argv[2]);
    newx = atoi (argv[3]);
    newy = atoi (argv[4]);

    SetVisualScene(interp);
    vs->MouseMove (oldx, oldy, newx, newy, argv[5][0]);

    return TCL_OK;
  }


  int Ng_MouseDblClick (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    int px, py;

    px = atoi (argv[1]);
    py = atoi (argv[2]);

    SetVisualScene(interp);
    vs->MouseDblClick (px, py);

    return TCL_OK;
  }


  int Ng_ZoomAll (ClientData clientData,
		  Tcl_Interp * interp,
		  int argc, tcl_const char *argv[])
  {
    SetVisualScene(interp);
    vs->BuildScene (1);

    return TCL_OK;
  }


  int Ng_Center (ClientData clientData,
		 Tcl_Interp * interp,
		 int argc, tcl_const char *argv[])
  {
    SetVisualScene(interp);
    vs->BuildScene (2);

    return TCL_OK;
  }


  int Ng_StandardRotation (ClientData clientData,
			   Tcl_Interp * interp,
			   int argc, tcl_const char *argv[])
  {
    SetVisualScene(interp);
    vs->StandardRotation (argv[1]);

    return TCL_OK;
  }

  int Ng_ArbitraryRotation (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[])
  {
    SetVisualScene(interp);
    Array<double> alpha;
    Array<Vec3d> vec;

    for(int i=1; i<argc; i+=4)
      {
	alpha.Append(atof(argv[i]));
	vec.Append(Vec3d(atof(argv[i+1]),atof(argv[i+2]),atof(argv[i+3])));
      }

    vs->ArbitraryRotation (alpha,vec);

    return TCL_OK;
  }


  int Ng_Metis (ClientData clientData,
		Tcl_Interp * interp,
		int argc, tcl_const char *argv[])
  {
#ifdef METISold

    if (!mesh)
      {
	Tcl_SetResult (interp, err_needsmesh, TCL_STATIC);
	return TCL_ERROR;
      }

    // METIS Partitioning

    if (mesh->GetDimension() == 3)
      {
	using namespace metis;

	int ne = mesh->GetNE();
	if (ne < 3)
	  {
	    Tcl_SetResult (interp, "This operation needs a volume mesh", TCL_STATIC);
	    return TCL_ERROR;
	  }

	int nn = mesh->GetNP();

	ELEMENT_TYPE elementtype = mesh->VolumeElement(1).GetType();
	int npe = mesh->VolumeElement(1).GetNP();

	for (int i = 2; i<=ne; i++)
	  if (mesh->VolumeElement(i).GetType() != elementtype)
	    {
	      Tcl_SetResult (interp, "Works in 3D only uniformal tet or hex meshes", TCL_STATIC);
	      return TCL_ERROR;
	    }

	idxtype *elmnts;
	elmnts = new idxtype[ne*npe];

	int etype;
	if (elementtype == TET)
	  etype = 2;
	else if (elementtype == HEX)
	  etype = 3;
	else
	  {
	    Tcl_SetResult (interp, "Works in 3D only uniformal tet or hex meshes", TCL_STATIC);
	    return TCL_ERROR;
	  }

	for (int i=1; i<=ne; i++)
	  for (int j=1; j<=npe; j++)
	    elmnts[(i-1)*npe+(j-1)] = mesh->VolumeElement(i).PNum(j)-1;


	int numflag = 0;
	int nparts = atoi (argv[1]);
	int edgecut;
	idxtype *epart, *npart;
	epart = new idxtype[ne];
	npart = new idxtype[nn];

	cout << "Starting Metis (" << ne << " Elements, " << nn << " Nodes, " << nparts << " Partitions) ... " << flush;

	METIS_PartMeshNodal (&ne, &nn, elmnts, &etype, &numflag, &nparts,
			     &edgecut, epart, npart);

	cout << "done" << endl;

	cout << "edge-cut: " << edgecut << ", balance: " << ComputeElementBalance(ne, nparts, epart) << endl;

	for (int i=1; i<=ne; i++)
	  mesh->VolumeElement(i).SetPartition(epart[i-1]);

	mesh->SetNextTimeStamp();
      }


#endif
    return TCL_OK;
  }






  int Ng_SetOCCVisParameters  (ClientData clientData,
			       Tcl_Interp * interp,
			       int argc, tcl_const char *argv[])
  {
#ifdef OCCGEOMETRY
    int showvolume;

    showvolume = atoi (Tcl_GetVar (interp, "::occoptions.showvolumenr", 0));

    if (showvolume != vispar.occshowvolumenr)
      {
	if (showvolume < 0 || showvolume > occgeometry->NrSolids())
	  {
	    char buf[20];
	    sprintf (buf, "%5i", vispar.occshowvolumenr);
            Tcl_SetVar (interp, "::occoptions.showvolumenr", buf, 0);
	  }
	else
	  {
	    vispar.occshowvolumenr = showvolume;
	    if (occgeometry)
	      occgeometry -> changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
	  }
      }

    int temp;

    temp = atoi (Tcl_GetVar (interp, "::occoptions.visproblemfaces", 0));

    if ((bool) temp != vispar.occvisproblemfaces)
      {
	vispar.occvisproblemfaces = temp;
	if (occgeometry)
	  occgeometry -> changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
      }

    vispar.occshowsurfaces = atoi (Tcl_GetVar (interp, "::occoptions.showsurfaces", 0));
    vispar.occshowedges = atoi (Tcl_GetVar (interp, "::occoptions.showedges", 0));
    vispar.occzoomtohighlightedentity = atoi (Tcl_GetVar (interp, "::occoptions.zoomtohighlightedentity", 0));
    vispar.occdeflection = pow(10.0,-1-atof (Tcl_GetVar (interp, "::occoptions.deflection", 0)));

#endif





#ifdef ACIS
    vispar.ACISshowfaces = atoi (Tcl_GetVar (interp, "::occoptions.showsurfaces", 0));
    vispar.ACISshowedges = atoi (Tcl_GetVar (interp, "::occoptions.showedges", 0));
    vispar.ACISshowsolidnr = atoi (Tcl_GetVar (interp, "::occoptions.showsolidnr", 0));
    vispar.ACISshowsolidnr2 = atoi (Tcl_GetVar (interp, "::occoptions.showsolidnr2", 0));

#endif



    return TCL_OK;
  }

  void SelectFaceInOCCDialogTree (int facenr)
  {
    char script[50];
    sprintf (script, "selectentity {Face %i}", facenr);
    Tcl_GlobalEval (tcl_interp, script);
  }



  int Ng_GetOCCData (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
#ifdef OCCGEOMETRY

    static char buf[1000];
    buf[0] = 0;
    stringstream str;

    if (argc >= 2)
      {
	if (strcmp (argv[1], "getentities") == 0)
	  {
	    if (occgeometry)
	      {
		occgeometry->GetTopologyTree(str);
	      }
	  }
      }

    Tcl_SetResult (interp, (char*)str.str().c_str(), TCL_VOLATILE);

#endif
    return TCL_OK;
  }

  int Ng_OCCCommand (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
#ifdef OCCGEOMETRY

    stringstream str;
    if (argc >= 2)
      {
	if (strcmp (argv[1], "isoccgeometryloaded") == 0)
	  {
	    if (occgeometry)
	      str << "1 " << flush;
	    else str << "0 " << flush;

	    Tcl_SetResult (interp, (char*)str.str().c_str(), TCL_VOLATILE);
	  }
	if (occgeometry)
	  {
	    if (strcmp (argv[1], "buildvisualizationmesh") == 0)
	      {
		occgeometry->BuildVisualizationMesh(vispar.occdeflection);
		occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
	      }
	    if (strcmp (argv[1], "mesherror") == 0)
	      {
		if (occgeometry->ErrorInSurfaceMeshing())
		  str << 1;
		else
		  str << 0;
	      }
	    if (strcmp (argv[1], "sewfaces") == 0)
	      {
		cout << "Before operation:" << endl;
		occgeometry->PrintNrShapes();
		occgeometry->SewFaces();
		occgeometry->BuildFMap();
		cout << endl << "After operation:" << endl;
		occgeometry->PrintNrShapes();
		occgeometry->BuildVisualizationMesh(vispar.occdeflection);
		occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
	      }
	    if (strcmp (argv[1], "makesolid") == 0)
	      {
		cout << "Before operation:" << endl;
		occgeometry->PrintNrShapes();
		occgeometry->MakeSolid();
		occgeometry->BuildFMap();
		cout << endl << "After operation:" << endl;
		occgeometry->PrintNrShapes();
		occgeometry->BuildVisualizationMesh(vispar.occdeflection);
		occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
	      }
	    if (strcmp (argv[1], "upgradetopology") == 0)
	      {
		cout << "Before operation:" << endl;
		occgeometry->PrintNrShapes();
		occgeometry->SewFaces();
		occgeometry->MakeSolid();
		occgeometry->BuildFMap();
		cout << endl << "After operation:" << endl;
		occgeometry->PrintNrShapes();
		occgeometry->BuildVisualizationMesh(vispar.occdeflection);
		occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
	      }
	    if (strcmp (argv[1], "shapehealing") == 0)
	      {
		occgeometry->tolerance =
		  atof (Tcl_GetVar (interp, "::occoptions.tolerance", 0));
		occgeometry->fixsmalledges =
		  atoi (Tcl_GetVar (interp, "::occoptions.fixsmalledges", 0));
		occgeometry->fixspotstripfaces =
		  atoi (Tcl_GetVar (interp, "::occoptions.fixspotstripfaces", 0));
		occgeometry->sewfaces =
		  atoi (Tcl_GetVar (interp, "::occoptions.sewfaces", 0));
		occgeometry->makesolids =
		  atoi (Tcl_GetVar (interp, "::occoptions.makesolids", 0));
		occgeometry->splitpartitions =
		  atoi (Tcl_GetVar (interp, "::occoptions.splitpartitions", 0));

		//	      cout << "Before operation:" << endl;
		//	      occgeometry->PrintNrShapes();
		occgeometry->HealGeometry();
		occgeometry->BuildFMap();
		//	      cout << endl << "After operation:" << endl;
		//	      occgeometry->PrintNrShapes();
		occgeometry->BuildVisualizationMesh(vispar.occdeflection);
		occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
	      }


	    if (strcmp (argv[1], "highlightentity") == 0)
	      {
		if (strcmp (argv[2], "Face") == 0)
		  {
		    int nr = atoi (argv[3]);
		    occgeometry->LowLightAll();

		    occgeometry->fvispar[nr-1].Highlight();
		    if (vispar.occzoomtohighlightedentity)
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONFULLCHANGE;
		    else
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
		  }
		if (strcmp (argv[2], "Shell") == 0)
		  {
		    int nr = atoi (argv[3]);
		    occgeometry->LowLightAll();

		    TopExp_Explorer exp;
		    for (exp.Init (occgeometry->shmap(nr), TopAbs_FACE);
			 exp.More(); exp.Next())
		      {
			int i = occgeometry->fmap.FindIndex (TopoDS::Face(exp.Current()));
			occgeometry->fvispar[i-1].Highlight();
		      }
		    if (vispar.occzoomtohighlightedentity)
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONFULLCHANGE;
		    else
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
		  }
		if (strcmp (argv[2], "Solid") == 0)
		  {
		    int nr = atoi (argv[3]);
		    occgeometry->LowLightAll();

		    TopExp_Explorer exp;
		    for (exp.Init (occgeometry->somap(nr), TopAbs_FACE);
			 exp.More(); exp.Next())
		      {
			int i = occgeometry->fmap.FindIndex (TopoDS::Face(exp.Current()));
			occgeometry->fvispar[i-1].Highlight();
		      }
		    if (vispar.occzoomtohighlightedentity)
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONFULLCHANGE;
		    else
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
		  }
		/*
		  if (strcmp (argv[2], "CompSolid") == 0)
		  {
		  int nr = atoi (argv[3]);
		  occgeometry->LowLightAll();

		  TopExp_Explorer exp;
		  for (exp.Init (occgeometry->cmap(nr), TopAbs_FACE);
		  exp.More(); exp.Next())
		  {
		  int i = occgeometry->fmap.FindIndex (TopoDS::Face(exp.Current()));
		  occgeometry->fvispar[i-1].Highlight();
		  }
		  occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
		  }
		*/

		if (strcmp (argv[2], "Edge") == 0)
		  {
		    int nr = atoi (argv[3]);
		    occgeometry->LowLightAll();

		    occgeometry->evispar[nr-1].Highlight();
		    if (vispar.occzoomtohighlightedentity)
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONFULLCHANGE;
		    else
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
		  }
		if (strcmp (argv[2], "Wire") == 0)
		  {
		    int nr = atoi (argv[3]);
		    occgeometry->LowLightAll();

		    TopExp_Explorer exp;
		    for (exp.Init (occgeometry->wmap(nr), TopAbs_EDGE);
			 exp.More(); exp.Next())
		      {
			int i = occgeometry->emap.FindIndex (TopoDS::Edge(exp.Current()));
			occgeometry->evispar[i-1].Highlight();
		      }
		    if (vispar.occzoomtohighlightedentity)
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONFULLCHANGE;
		    else
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
		  }

		if (strcmp (argv[2], "Vertex") == 0)
		  {
		    int nr = atoi (argv[3]);
		    occgeometry->LowLightAll();

		    occgeometry->vvispar[nr-1].Highlight();
		    if (vispar.occzoomtohighlightedentity)
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONFULLCHANGE;
		    else
		      occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
		  }

	      }



	    if (strcmp (argv[1], "show") == 0)
	      {
		int nr = atoi (argv[3]);
		occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;

		if (strcmp (argv[2], "Face") == 0)
		  {
		    occgeometry->fvispar[nr-1].Show();
		  }
		if (strcmp (argv[2], "Shell") == 0)
		  {
		    TopExp_Explorer exp;
		    for (exp.Init (occgeometry->shmap(nr), TopAbs_FACE);
			 exp.More(); exp.Next())
		      {
			int i = occgeometry->fmap.FindIndex (TopoDS::Face(exp.Current()));
			occgeometry->fvispar[i-1].Show();
		      }
		  }
		if (strcmp (argv[2], "Solid") == 0)
		  {
		    TopExp_Explorer exp;
		    for (exp.Init (occgeometry->somap(nr), TopAbs_FACE);
			 exp.More(); exp.Next())
		      {
			int i = occgeometry->fmap.FindIndex (TopoDS::Face(exp.Current()));
			occgeometry->fvispar[i-1].Show();
		      }
		  }
		if (strcmp (argv[2], "Edge") == 0)
		  {
		    occgeometry->evispar[nr-1].Show();
		  }
		if (strcmp (argv[2], "Wire") == 0)
		  {
		    TopExp_Explorer exp;
		    for (exp.Init (occgeometry->wmap(nr), TopAbs_EDGE);
			 exp.More(); exp.Next())
		      {
			int i = occgeometry->emap.FindIndex (TopoDS::Edge(exp.Current()));
			occgeometry->evispar[i-1].Show();
		      }
		  }
	      }


	    if (strcmp (argv[1], "hide") == 0)
	      {
		int nr = atoi (argv[3]);
		occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;

		if (strcmp (argv[2], "Face") == 0)
		  {
		    occgeometry->fvispar[nr-1].Hide();
		  }
		if (strcmp (argv[2], "Shell") == 0)
		  {
		    TopExp_Explorer exp;
		    for (exp.Init (occgeometry->shmap(nr), TopAbs_FACE);
			 exp.More(); exp.Next())
		      {
			int i = occgeometry->fmap.FindIndex (TopoDS::Face(exp.Current()));
			occgeometry->fvispar[i-1].Hide();
		      }
		  }
		if (strcmp (argv[2], "Solid") == 0)
		  {
		    TopExp_Explorer exp;
		    for (exp.Init (occgeometry->somap(nr), TopAbs_FACE);
			 exp.More(); exp.Next())
		      {
			int i = occgeometry->fmap.FindIndex (TopoDS::Face(exp.Current()));
			occgeometry->fvispar[i-1].Hide();
		      }
		  }
		if (strcmp (argv[2], "Edge") == 0)
		  {
		    occgeometry->evispar[nr-1].Hide();
		  }
		if (strcmp (argv[2], "Wire") == 0)
		  {
		    TopExp_Explorer exp;
		    for (exp.Init (occgeometry->wmap(nr), TopAbs_EDGE);
			 exp.More(); exp.Next())
		      {
			int i = occgeometry->emap.FindIndex (TopoDS::Edge(exp.Current()));
			occgeometry->evispar[i-1].Hide();
		      }
		  }
	      }



	    if (strcmp (argv[1], "findsmallentities") == 0)
	      {
		stringstream str("");
		occgeometry->CheckIrregularEntities(str);
		Tcl_SetResult (interp, (char*)str.str().c_str(), TCL_VOLATILE);
	      }
	    if (strcmp (argv[1], "getunmeshedfaceinfo") == 0)
	      {
		occgeometry->GetUnmeshedFaceInfo(str);
		Tcl_SetResult (interp, (char*)str.str().c_str(), TCL_VOLATILE);
	      }
	    if (strcmp (argv[1], "getnotdrawablefaces") == 0)
	      {
		occgeometry->GetNotDrawableFaces(str);
		Tcl_SetResult (interp, (char*)str.str().c_str(), TCL_VOLATILE);
	      }
	    if (strcmp (argv[1], "redrawstatus") == 0)
	      {
		int i = atoi (argv[2]);
		occgeometry->changed = i;
	      }
	    if (strcmp (argv[1], "swaporientation") == 0)
	      {
		IGESControl_Writer writer("millimeters", 1);
		writer.AddShape (occgeometry->shape);
		writer.Write ("1.igs");
		/*
		  int nr = atoi (argv[3]);

		  //	      const_cast<TopoDS_Shape&> (occgeometry->fmap(nr)).Reverse();

		  Handle_ShapeBuild_ReShape rebuild = new ShapeBuild_ReShape;
		  rebuild->Apply(occgeometry->shape);

		  TopoDS_Shape sh;

		  //	      if (strcmp (argv[2], "CompSolid") == 0) sh = occgeometry->cmap(nr);
		  if (strcmp (argv[2], "Solid") == 0) sh = occgeometry->somap(nr);
		  if (strcmp (argv[2], "Shell") == 0) sh = occgeometry->shmap(nr);
		  if (strcmp (argv[2], "Face") == 0) sh = occgeometry->fmap(nr);
		  if (strcmp (argv[2], "Wire") == 0) sh = occgeometry->wmap(nr);
		  if (strcmp (argv[2], "Edge") == 0) sh = occgeometry->emap(nr);

		  rebuild->Replace(sh, sh.Reversed(), Standard_False);

		  TopoDS_Shape newshape = rebuild->Apply(occgeometry->shape, TopAbs_SHELL, 1);
		  occgeometry->shape = newshape;

		  occgeometry->BuildFMap();
		  occgeometry->BuildVisualizationMesh();
		  occgeometry->changed = OCCGEOMETRYVISUALIZATIONHALFCHANGE;
		*/
	      }
	    if (strcmp (argv[1], "marksingular") == 0)
	      {
		int nr = atoi (argv[3]);
		cout << "marking " << argv[2] << " " << nr << endl;
		char buf[2]; buf[0] = '0'; buf[1] = 0;
		bool sing = false;
		if (strcmp (argv[2], "Face") == 0)
		  sing = occgeometry->fsingular[nr-1] = !occgeometry->fsingular[nr-1];
		if (strcmp (argv[2], "Edge") == 0)
		  sing = occgeometry->esingular[nr-1] = !occgeometry->esingular[nr-1];
		if (strcmp (argv[2], "Vertex") == 0)
		  sing = occgeometry->vsingular[nr-1] = !occgeometry->vsingular[nr-1];

		if (sing) buf[0] = '1';

                Tcl_SetVar (interp, "::ismarkedsingular", buf, 0);

		stringstream str;
		occgeometry->GetTopologyTree (str);

		char* cstr = (char*)str.str().c_str();

		(*testout) << cstr << endl;

		char helpstr[1000];

		while (strchr (cstr, '}'))
		  {
		    strncpy (helpstr, cstr+2, strlen(strchr(cstr+2, '}')));
		    (*testout) << "***" << cstr << "***" << endl;
		    cstr = strchr (cstr, '}');
		  }
	      }
	  }
      }

#endif
    return TCL_OK;
  }



#ifdef OCCGEOMETRY

  void OCCConstructGeometry (OCCGeometry & geom);

  int Ng_OCCConstruction (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
    if (occgeometry)
      OCCConstructGeometry (*occgeometry);
    return TCL_OK;
  }
#endif



#ifndef ACIS
  int Ng_ACISCommand (ClientData clientData,
		      Tcl_Interp * interp,
		      int argc, tcl_const char *argv[])
  {
    if (argc >= 2)
      {
	if (strcmp (argv[1], "isACISavailable") == 0)
	  {
            Tcl_SetResult (interp, (char*)"no", TCL_STATIC);
	    return TCL_OK;
	  }
      }
    Tcl_SetResult (interp, (char*)"undefined ACiS command", TCL_STATIC);
    return TCL_ERROR;
  }
#endif




  int Ng_SetVisParameters  (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[])
  {
    if (!Tcl_GetVar (interp, "::viewoptions.light.amb", TCL_GLOBAL_ONLY))
      return TCL_ERROR;

    vispar.lightamb = atof (Tcl_GetVar (interp, "::viewoptions.light.amb", TCL_GLOBAL_ONLY));
    vispar.lightdiff = atof (Tcl_GetVar (interp, "::viewoptions.light.diff", TCL_GLOBAL_ONLY));
    vispar.lightspec = atof (Tcl_GetVar (interp, "::viewoptions.light.spec", TCL_GLOBAL_ONLY));
    vispar.shininess = atof (Tcl_GetVar (interp, "::viewoptions.mat.shininess", TCL_GLOBAL_ONLY));
    vispar.locviewer = atoi (Tcl_GetVar (interp, "::viewoptions.light.locviewer", TCL_GLOBAL_ONLY));
    vispar.transp = atof (Tcl_GetVar (interp, "::viewoptions.mat.transp", TCL_GLOBAL_ONLY));

    vispar.clipnormal.X() = atof (Tcl_GetVar (interp, "::viewoptions.clipping.nx", TCL_GLOBAL_ONLY));
    vispar.clipnormal.Y() = atof (Tcl_GetVar (interp, "::viewoptions.clipping.ny", TCL_GLOBAL_ONLY));
    vispar.clipnormal.Z() = atof (Tcl_GetVar (interp, "::viewoptions.clipping.nz", TCL_GLOBAL_ONLY));
    vispar.clipdist = atof (Tcl_GetVar (interp, "::viewoptions.clipping.dist", TCL_GLOBAL_ONLY));
    vispar.clipenable = atoi (Tcl_GetVar (interp, "::viewoptions.clipping.enable", TCL_GLOBAL_ONLY));
    vispar.clipdomain =
      atoi (Tcl_GetVar (interp, "::viewoptions.clipping.onlydomain", TCL_GLOBAL_ONLY));
    vispar.donotclipdomain =
      atoi (Tcl_GetVar (interp, "::viewoptions.clipping.notdomain", TCL_GLOBAL_ONLY));

    vispar.clipplanetimestamp = NextTimeStamp();


    vispar.whitebackground = atoi (Tcl_GetVar (interp, "::viewoptions.whitebackground", TCL_GLOBAL_ONLY));
    vispar.drawcoordinatecross = atoi (Tcl_GetVar (interp, "::viewoptions.drawcoordinatecross", TCL_GLOBAL_ONLY));
    vispar.drawcolorbar = atoi (Tcl_GetVar (interp, "::viewoptions.drawcolorbar", TCL_GLOBAL_ONLY));
    vispar.drawnetgenlogo = atoi (Tcl_GetVar (interp, "::viewoptions.drawnetgenlogo", TCL_GLOBAL_ONLY));
    vispar.stereo = atoi (Tcl_GetVar (interp, "::viewoptions.stereo", TCL_GLOBAL_ONLY));
    vispar.colormeshsize = atoi (Tcl_GetVar (interp, "::viewoptions.colormeshsize", TCL_GLOBAL_ONLY));
    VisualScene :: SetBackGroundColor (vispar.whitebackground ? 1 : 0);

    strcpy (vispar.selectvisual, Tcl_GetVar (interp, "::selectvisual", TCL_GLOBAL_ONLY));

    //  vispar.showstltrias = atoi (Tcl_GetVar (interp, "::viewoptions.stl.showtrias", TCL_GLOBAL_ONLY));

    vispar.stlshowtrias =
      atoi (Tcl_GetVar (interp, "::stloptions.showtrias", TCL_GLOBAL_ONLY));
    vispar.stlshowfilledtrias =
      atoi (Tcl_GetVar (interp, "::stloptions.showfilledtrias", TCL_GLOBAL_ONLY));
    vispar.stlshowedges =
      atoi (Tcl_GetVar (interp, "::stloptions.showedges", TCL_GLOBAL_ONLY));
    vispar.stlshowmarktrias =
      atoi (Tcl_GetVar (interp, "::stloptions.showmarktrias", TCL_GLOBAL_ONLY));
    vispar.stlshowactivechart =
      atoi (Tcl_GetVar (interp, "::stloptions.showactivechart", TCL_GLOBAL_ONLY));
    vispar.stlchartnumber =
      atoi (Tcl_GetVar (interp, "::stloptions.chartnumber", TCL_GLOBAL_ONLY));
    vispar.stlchartnumberoffset =
      atoi (Tcl_GetVar (interp, "::stloptions.chartnumberoffset", TCL_GLOBAL_ONLY));

    vispar.occshowsurfaces =
      atoi (Tcl_GetVar (interp, "::occoptions.showsurfaces", TCL_GLOBAL_ONLY));
    vispar.occshowedges =
      atoi (Tcl_GetVar (interp, "::occoptions.showedges", TCL_GLOBAL_ONLY));

    vispar.drawoutline =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawoutline", TCL_GLOBAL_ONLY));
    vispar.drawfilledtrigs =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawfilledtrigs", TCL_GLOBAL_ONLY));
    vispar.subdivisions =
      atoi (Tcl_GetVar (interp, "::visoptions.subdivisions", TCL_GLOBAL_ONLY));
    vispar.drawbadels =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawbadels", TCL_GLOBAL_ONLY));
    vispar.drawedges =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawedges", TCL_GLOBAL_ONLY));

    vispar.drawtetsdomain =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawtetsdomain", TCL_GLOBAL_ONLY));
    vispar.drawtets =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawtets", TCL_GLOBAL_ONLY));
    vispar.drawprisms =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawprisms", TCL_GLOBAL_ONLY));
    vispar.drawpyramids =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawpyramids", TCL_GLOBAL_ONLY));
    vispar.drawhexes =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawhexes", TCL_GLOBAL_ONLY));
    vispar.shrink =
      atof (Tcl_GetVar (interp, "::viewoptions.shrink", TCL_GLOBAL_ONLY));
    vispar.drawidentified =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawidentified", TCL_GLOBAL_ONLY));
    vispar.drawpointnumbers =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawpointnumbers", TCL_GLOBAL_ONLY));
    vispar.drawedgenumbers =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawedgenumbers", TCL_GLOBAL_ONLY));
    vispar.drawfacenumbers =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawfacenumbers", TCL_GLOBAL_ONLY));
    vispar.drawelementnumbers =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawelementnumbers", TCL_GLOBAL_ONLY));
    vispar.drawdomainsurf =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawdomainsurf", TCL_GLOBAL_ONLY));

    vispar.drawededges =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawededges", TCL_GLOBAL_ONLY));
    vispar.drawedpoints =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawedpoints", TCL_GLOBAL_ONLY));
    vispar.drawedpointnrs =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawedpointnrs", TCL_GLOBAL_ONLY));
    vispar.drawedtangents =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawedtangents", TCL_GLOBAL_ONLY));
    vispar.drawededgenrs =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawededgenrs", TCL_GLOBAL_ONLY));

    vispar.drawcurveproj =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawcurveproj", TCL_GLOBAL_ONLY));
    vispar.drawcurveprojedge =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawcurveprojedge", TCL_GLOBAL_ONLY));

    vispar.centerpoint =
      atoi (Tcl_GetVar (interp, "::viewoptions.centerpoint", TCL_GLOBAL_ONLY));
    vispar.use_center_coords =
      atoi (Tcl_GetVar (interp, "::viewoptions.usecentercoords", TCL_GLOBAL_ONLY)) > 0;
    vispar.centerx =
      atof (Tcl_GetVar (interp, "::viewoptions.centerx", TCL_GLOBAL_ONLY));
    vispar.centery =
      atof (Tcl_GetVar (interp, "::viewoptions.centery", TCL_GLOBAL_ONLY));
    vispar.centerz =
      atof (Tcl_GetVar (interp, "::viewoptions.centerz", TCL_GLOBAL_ONLY));
    vispar.drawelement =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawelement", TCL_GLOBAL_ONLY));
    vispar.drawmetispartition =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawmetispartition", TCL_GLOBAL_ONLY));

    vispar.drawspecpoint =
      atoi (Tcl_GetVar (interp, "::viewoptions.drawspecpoint", TCL_GLOBAL_ONLY));
    vispar.specpointx =
      atof (Tcl_GetVar (interp, "::viewoptions.specpointx", TCL_GLOBAL_ONLY));
    vispar.specpointy =
      atof (Tcl_GetVar (interp, "::viewoptions.specpointy", TCL_GLOBAL_ONLY));
    vispar.specpointz =
      atof (Tcl_GetVar (interp, "::viewoptions.specpointz", TCL_GLOBAL_ONLY));


    vsspecpoints.len =
      atof (Tcl_GetVar (interp, "::viewoptions.specpointvlen", TCL_GLOBAL_ONLY));




#ifdef OCCGEOMETRY
    vispar.occdeflection = pow(10.0,-1-atof (Tcl_GetVar (interp, "::occoptions.deflection", TCL_GLOBAL_ONLY)));
#endif

#ifdef PARALLELGL
    vsmesh.Broadcast ();
#endif

    return TCL_OK;
  }



  int Ng_SelectSurface (ClientData clientData,
			Tcl_Interp * interp,
			int argc, tcl_const char *argv[])
  {
    int surfnr = atoi (argv[1]);
    vsgeom.SelectSurface (surfnr);
    return TCL_OK;
  }


  int Ng_BuildFieldLines (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
    vssolution.BuildFieldLinesPlot();
    return TCL_OK;
  }

#ifdef PARALLEL
  int Ng_VisualizeAll (ClientData clientData,
                       Tcl_Interp * interp,
                       int argc, tcl_const char *argv[])
  {
    int id, rc, ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    string visualizationmode = Tcl_GetVar (interp, "::selectvisual", 0);
    string scalfun = Tcl_GetVar (interp, "::visoptions.scalfunction", 0);
    for ( int dest = 1; dest < ntasks; dest++)
      {
	MyMPI_Send ( "visualize", dest );
	MyMPI_Send ( visualizationmode, dest);
	if ( visualizationmode == "solution" )
	  MyMPI_Send ( scalfun, dest);
      }
    return TCL_OK;
  }

  int Ng_VisualizeOne (ClientData clientData,
                       Tcl_Interp * interp,
                       int argc, tcl_const char *argv[])
  {
    int id, rc, ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    string visualizationmode = Tcl_GetVar (interp, "::selectvisual", 0);
    string scalfun = Tcl_GetVar (interp, "::visoptions.scalfunction", 0);

    MyMPI_Send ( "visualize", 1 );
    MyMPI_Send ( visualizationmode, 1);

    if ( visualizationmode == "solution" )
      MyMPI_Send ( scalfun, 1);
    return TCL_OK;
  }

  int Ng_IncrOverlap ( ClientData clientDate,
		       Tcl_Interp * interp,
		       int argc, tcl_const char * argv[] )
  {
    int id, rc, ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);

    for ( int dest = 1; dest < ntasks; dest++)
      {
	MyMPI_Send ( "overlap++", dest );
      }
    mesh->UpdateOverlap();
    return TCL_OK;

  }

  int Ng_SetSelectVisual (  ClientData clientDate,
                            Tcl_Interp * interp,
			    int argc, tcl_const char * argv[] )
  {
    string visualizationmode;
    MyMPI_Recv ( visualizationmode, 0);
    Tcl_SetVar (interp, "::selectvisual", visualizationmode.c_str(), 0);
    return TCL_OK;
  }

  int Ng_SetScalarFunction (  ClientData clientDate,
                              Tcl_Interp * interp,
                              int argc, tcl_const char * argv[] )
  {
    string visualizationmode;
    string scalarfun;
    visualizationmode = Tcl_GetVar (interp, "::selectvisual", 0);

    if ( visualizationmode == "solution" )
      {
	MyMPI_Recv ( scalarfun, 0);
        Tcl_SetVar (interp, "::visoptions.scalfunction", scalarfun.c_str(), 0);
      }
    return TCL_OK;
  }


#endif

  int Ng_IsParallel (ClientData clientData,
                     Tcl_Interp * interp,
                     int argc, tcl_const char *argv[])
  {
#ifdef PARALLEL
    int id, rc, ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if ( ntasks > 1 )
      Tcl_SetVar (interp, "::parallel_netgen", "1", 0);
    else
      Tcl_SetVar (interp, "::parallel_netgen", "0", 0);
#else
    Tcl_SetVar (interp, "::parallel_netgen", "0", 0);
#endif
    return TCL_OK;
  }

  int Ng_Exit (ClientData clientData,
	       Tcl_Interp * interp,
	       int argc, tcl_const char *argv[])
  {
#ifdef PARALLEL
    int id, rc, ntasks;
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if ( id != 0 )
      return TCL_OK;
#endif


    /*
    if (ngsolve_handle)
      {
        void (*ngs_exit)();
        ngs_exit = ( void (*)() ) dlsym (ngsolve_handle, "NGSolve_Exit");
        if (ngs_exit) (*ngs_exit)();
      }
    */

#ifdef NGSOLVE
    NGSolve_Exit ();
#endif




    delete stlgeometry;
    stlgeometry = NULL;

    geometry.Reset (0);
    geometry2d.Reset (0);


#ifdef ACIS
    outcome res;
    res = api_terminate_faceter();
    if(!res.ok())
      cerr << "problem with terminating acis faceter" << endl;
    res = api_terminate_constructors();
    if(!res.ok())
      cerr << "problem with terminating acis constructors" << endl;
    res = api_terminate_kernel();
    if(!res.ok())
      cerr << "problem with terminating acis kernel" << endl;
    res = api_stop_modeller();
    if(!res.ok())
      cerr << "problem with terminating acis modeller" << endl;

    //cout << "stopped acis, outcome = " << res.ok() << endl;
#endif

#ifdef PARALLEL

    for ( int dest = 1; dest < ntasks; dest++)
      MyMPI_Send ( "end", dest );
#endif

    if (testout != &cout)
      delete testout;

    return TCL_OK;
  }


#ifdef SOCKETS
  void * ServerSocketManagerRunDummy ( void * nix )
  {
    serversocketmanager.Run();
    return NULL;
  }

  extern "C" int Ng_ServerSocketManagerRun( void );
  int Ng_ServerSocketManagerRun( void )
  {
    if(mparam.parthread)
      RunParallel(ServerSocketManagerRunDummy,NULL);
    else
      serversocketmanager.Run();

    return TCL_OK;
  }

  extern "C" int Ng_ServerSocketManagerInit(int port);
  int Ng_ServerSocketManagerInit(int port)
  {
    serversocketmanager.Init(port);
    return TCL_OK;
  }
#endif //SOCKETS


  extern "C" int Ng_Init (Tcl_Interp * interp);

  //   int main_Eero (ClientData clientData,
  // 	       Tcl_Interp * interp,
  // 		 int argc, tcl_const char *argv[]);


  int Ng_Init (Tcl_Interp * interp)
  {
#ifdef SOCKETS
    if(serversocketmanager.Good())
      serversocketusernetgen.Reset(new ServerSocketUserNetgen (serversocketmanager, mesh, geometry));
#endif

    tcl_interp = interp;

    //     Tcl_CreateCommand (interp, "Ng_Eero", main_Eero,
    // 		       (ClientData)NULL,
    // 		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_New", Ng_New,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    //     Tcl_CreateCommand (interp, "Ng_Lock", Ng_Lock,
    // 		       (ClientData)NULL,
    // 		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_LoadGeometry", Ng_LoadGeometry,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SaveGeometry", Ng_SaveGeometry,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ParseGeometry", Ng_ParseGeometry,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_LoadMesh", Ng_LoadMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SaveMesh", Ng_SaveMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MergeMesh", Ng_MergeMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ExportMesh", Ng_ExportMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ImportMesh", Ng_ImportMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ImportSolution", Ng_ImportSolution,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ShowDemo", Ng_ShowDemo,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_DemoSetTime", Ng_DemoSetTime,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SaveSolution", Ng_SaveSolution,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GeometryOptions", Ng_GeometryOptions,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    // geometry
    Tcl_CreateCommand (interp, "Ng_CreatePrimitive", Ng_CreatePrimitive,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetPrimitiveData", Ng_SetPrimitiveData,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GetPrimitiveData", Ng_GetPrimitiveData,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GetPrimitiveList", Ng_GetPrimitiveList,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_GetSurfaceList", Ng_GetSurfaceList,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);



    Tcl_CreateCommand (interp, "Ng_SetSolidData", Ng_SetSolidData,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GetSolidData", Ng_GetSolidData,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GetSolidList", Ng_GetSolidList,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_TopLevel", Ng_TopLevel,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    // Philippose - 30/01/2009
    // Register the TCL Interface Command for local face mesh size
    // definition
    Tcl_CreateCommand (interp, "Ng_SurfaceMeshSize", Ng_SurfaceMeshSize,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_AutoColourBcProps", Ng_AutoColourBcProps,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);
    

    // meshing
    Tcl_CreateCommand (interp, "Ng_GenerateMesh", Ng_GenerateMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_StopMeshing", Ng_StopMeshing,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MeshInfo", Ng_MeshInfo,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MeshQuality", Ng_MeshQuality,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_CheckSurfaceMesh", Ng_CheckSurfaceMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_CheckVolumeMesh", Ng_CheckVolumeMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_DeleteVolMesh", Ng_DeleteVolMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SplitSeparatedFaces", Ng_SplitSeparatedFaces,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetNextTimeStamp", Ng_SetNextTimeStamp,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Refine", Ng_Refine,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SecondOrder", Ng_SecondOrder,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_HighOrder", Ng_HighOrder,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ValidateSecondOrder", Ng_ValidateSecondOrder,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_RestrictH", Ng_RestrictH,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Anisotropy", Ng_Anisotropy,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Bisect", Ng_Bisect,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    //     Tcl_CreateCommand (interp, "Ng_BisectCopyMesh", Ng_BisectCopyMesh,
    // 		       (ClientData)NULL,
    // 		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Split2Tets", Ng_Split2Tets,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ZRefinement", Ng_ZRefinement,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_HPRefinement", Ng_HPRefinement,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_LoadMeshSize", Ng_LoadMeshSize,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MeshSizeFromSurfaceMesh", Ng_MeshSizeFromSurfaceMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SingularEdgeMS", Ng_SingularEdgeMS,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SingularPointMS", Ng_SingularPointMS,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GenerateBoundaryLayer", Ng_GenerateBoundaryLayer,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_InsertVirtualBL", Ng_InsertVirtualBL,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_CutOffAndCombine", Ng_CutOffAndCombine,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_HelmholtzMesh", Ng_HelmholtzMesh,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ReadStatus", Ng_ReadStatus,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MemInfo", Ng_MemInfo,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_STLDoctor", Ng_STLDoctor,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MeshDoctor", Ng_MeshDoctor,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_BCProp", Ng_BCProp,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_STLInfo", Ng_STLInfo,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_STLCalcLocalH",
		       Ng_STLCalcLocalH,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetOCCVisParameters",
		       Ng_SetOCCVisParameters,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GetOCCData",
		       Ng_GetOCCData,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

#ifdef OCCGEOMETRY
    Tcl_CreateCommand (interp, "Ng_OCCConstruction",
		       Ng_OCCConstruction,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);
#endif

    Tcl_CreateCommand (interp, "Ng_OCCCommand",
		       Ng_OCCCommand,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ACISCommand",
		       Ng_ACISCommand,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_MouseMove", Ng_MouseMove,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_MouseDblClick", Ng_MouseDblClick,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_ZoomAll", Ng_ZoomAll,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Center", Ng_Center,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_StandardRotation", Ng_StandardRotation,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_ArbitraryRotation", Ng_ArbitraryRotation,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_SetVisParameters", Ng_SetVisParameters,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetMeshingParameters", Ng_SetMeshingParameters,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetDebugParameters", Ng_SetDebugParameters,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetSTLParameters", Ng_SetSTLParameters,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

#ifdef OCCGEOMETRY
    Tcl_CreateCommand (interp, "Ng_SetOCCParameters", Ng_SetOCCParameters,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);
#endif //OCCGEOMETRY


    Tcl_CreateCommand (interp, "Ng_SelectSurface", Ng_SelectSurface,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GetCommandLineParameter",
		       Ng_GetCommandLineParameter,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Exit",
		       Ng_Exit,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Metis",
		       Ng_Metis,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_BuildFieldLines",
		       Ng_BuildFieldLines,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


#ifdef PARALLEL
    Tcl_CreateCommand (interp, "Ng_VisualizeAll", Ng_VisualizeAll,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_VisualizeOne", Ng_VisualizeOne,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_IncrOverlap", Ng_IncrOverlap,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetSelectVisual", Ng_SetSelectVisual,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_SetScalarFunction",  Ng_SetScalarFunction,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

#endif
    Tcl_CreateCommand (interp, "Ng_IsParallel", Ng_IsParallel,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    /*
     * Specify the C callback functions for widget creation, display,
     * and reshape.
     */
#if TOGL_MAJOR_VERSION==1
    if (!nodisplay)
      {
	if (Togl_Init(interp) == TCL_ERROR) 
	  return TCL_ERROR;
	
	Togl_CreateFunc( init );
	Togl_DestroyFunc( zap );
	Togl_DisplayFunc( draw );
	Togl_ReshapeFunc( reshape );
	//   Togl_TimerFunc(  idle );
	Togl_CreateCommand( (char*)"Ng_SnapShot", Ng_SnapShot);
	Togl_CreateCommand( (char*)"Ng_VideoClip", Ng_VideoClip);
	//   Togl_CreateCommand("position",position);
      }
#else
    if (!nodisplay)
      {
	if (Togl_Init(interp) == TCL_ERROR) 
	  return TCL_ERROR;
	
	
	Tcl_CreateObjCommand(interp, "init", init, NULL, NULL);
	Tcl_CreateObjCommand(interp, "zap", zap, NULL, NULL);
	Tcl_CreateObjCommand(interp, "draw", draw, NULL, NULL);
	Tcl_CreateObjCommand(interp, "reshape", reshape, NULL, NULL);
	
	//   Togl_TimerFunc(  idle );
	// Togl_CreateCommand( (char*)"Ng_SnapShot", Ng_SnapShot);
	// Togl_CreateCommand( (char*)"Ng_VideoClip", Ng_VideoClip);
	//   Togl_CreateCommand("position",position);
      }

#endif


    multithread.pause = 0;
    multithread.testmode = 0;
    multithread.redraw = 0;
    multithread.drawing = 1;
    multithread.terminate = 0;
    multithread.running = 0;

    multithread.task = "";
    multithread.percent = 20;


    Tcl_LinkVar (interp, "multithread_pause",
		 (char*)&multithread.pause, TCL_LINK_INT);
    Tcl_LinkVar (interp, "multithread_testmode",
		 (char*)&multithread.testmode, TCL_LINK_INT);
    Tcl_LinkVar (interp, "multithread_redraw",
		 (char*)&multithread.redraw, TCL_LINK_INT);
    Tcl_LinkVar (interp, "multithread_drawing",
		 (char*)&multithread.drawing, TCL_LINK_INT);
    Tcl_LinkVar (interp, "multithread_terminate",
		 (char*)&multithread.terminate, TCL_LINK_INT);
    Tcl_LinkVar (interp, "multithread_running",
		 (char*)&multithread.running, TCL_LINK_INT);


    //testout->setstate(ios_base::badbit);

    myerr = &cerr;
    extern ostream * mycout;
    mycout = &cout;

    testmode = 0;

#ifdef ACIS
    outcome res;
    res = api_start_modeller (0);
    if(!res.ok())
      cerr << "problem with starting acis modeller" << endl;

#ifdef ACIS_R17
    unlock_spatial_products_661();
#endif
    res = api_initialize_kernel();
    if(!res.ok())
      cerr << "problem with starting acis kernel" << endl;
    res = api_initialize_constructors();
    if(!res.ok())
      cerr << "problem with starting acis constructors" << endl;
    res = api_initialize_faceter();
    if(!res.ok())
      cerr << "problem with starting acis faceter" << endl;
#endif


    return TCL_OK;
  }


}


