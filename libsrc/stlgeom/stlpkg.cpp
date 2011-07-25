#include <mystdlib.h>
#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>

#include <meshing.hpp>


#include <incvis.hpp>
#include <visual.hpp>

#include <stlgeom.hpp>

#include "vsstl.hpp"

extern "C" int Ng_STL_Init (Tcl_Interp * interp);



namespace netgen
{
  extern NetgenGeometry * ng_geometry;
  extern AutoPtr<Mesh> mesh;

  static VisualSceneSTLGeometry vsstlgeom;
  static VisualSceneSTLMeshing vsstlmeshing;

  char * err_needsstlgeometry = (char*) "This operation needs an STL geometry";





  class STLGeometryRegister : public GeometryRegister
  {
  public:
    virtual NetgenGeometry * Load (string filename) const;
    virtual VisualScene * GetVisualScene (const NetgenGeometry * geom) const;
    virtual void SetParameters (Tcl_Interp * interp) 
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
    }
  };



  int Ng_SetSTLParameters  (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[])
  {
    STLGeometryRegister reg;
    reg.SetParameters (interp);

    return TCL_OK;
  }

  






  int Ng_STLDoctor (ClientData clientData,
		    Tcl_Interp * interp,
		    int argc, tcl_const char *argv[])
  {
    //cout << "STL doctor" << endl;
    STLGeometry * stlgeometry = 
      dynamic_cast<STLGeometry*> (ng_geometry);
      

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









  NetgenGeometry *  STLGeometryRegister :: Load (string filename) const
  {
    const char * cfilename = filename.c_str();

    if (strcmp (&cfilename[strlen(cfilename)-3], "stl") == 0)
      {
	PrintMessage (1, "Load STL geometry file ", cfilename);

	ifstream infile(cfilename);

	STLGeometry * hgeom = STLGeometry :: Load (infile);
	hgeom -> edgesfound = 0;
	return hgeom;
      }
    else if (strcmp (&cfilename[strlen(cfilename)-4], "stlb") == 0)
      {
	PrintMessage (1, "Load STL binary geometry file ", cfilename);

	ifstream infile(cfilename);

	STLGeometry * hgeom = STLGeometry :: LoadBinary (infile);
	hgeom -> edgesfound = 0;
	return hgeom;
      }
    else if (strcmp (&cfilename[strlen(cfilename)-3], "nao") == 0)
      {
	PrintMessage (1, "Load naomi (F. Kickinger) geometry file ", cfilename);

	ifstream infile(cfilename);

	STLGeometry * hgeom = STLGeometry :: LoadNaomi (infile);
	hgeom -> edgesfound = 0;
	return hgeom;
      }

    
    return NULL;
  }









  int Ng_STLInfo  (ClientData clientData,
		   Tcl_Interp * interp,
		   int argc, tcl_const char *argv[])
  {
    double data[10];
    static char buf[20];

    STLGeometry * stlgeometry = dynamic_cast<STLGeometry*> (ng_geometry);

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



  extern int Ng_SetMeshingParameters  (ClientData clientData,
				       Tcl_Interp * interp,
				       int argc, tcl_const char *argv[]);

  int Ng_STLCalcLocalH  (ClientData clientData,    
			 Tcl_Interp * interp,
			 int argc, tcl_const char *argv[])
  {
    for (int i = 0; i < geometryregister.Size(); i++)
      geometryregister[i] -> SetParameters (interp);


    Ng_SetMeshingParameters (clientData, interp, argc, argv);

    STLGeometry * stlgeometry = dynamic_cast<STLGeometry*> (ng_geometry);
    if (mesh.Ptr() && stlgeometry)
      {
	mesh -> SetLocalH (stlgeometry->GetBoundingBox().PMin() - Vec3d(10, 10, 10),
			   stlgeometry->GetBoundingBox().PMax() + Vec3d(10, 10, 10),
			   mparam.grading);
	stlgeometry -> RestrictLocalH(*mesh, mparam.maxh);

	if (stlparam.resthsurfmeshcurvenable)
	  mesh -> CalcLocalHFromSurfaceCurvature (mparam.grading, 
						  stlparam.resthsurfmeshcurvfac);
      }

    return TCL_OK;
  }




  VisualScene * STLGeometryRegister :: GetVisualScene (const NetgenGeometry * geom) const
  {
    STLGeometry * geometry = dynamic_cast<STLGeometry*> (ng_geometry);
    if (geometry)
      {
	vsstlmeshing.SetGeometry (geometry);
	return &vsstlmeshing;
      }
    return NULL;
  }
}


using namespace netgen;

extern "C" int Ng_stl_Init (Tcl_Interp * interp);
int Ng_stl_Init (Tcl_Interp * interp)
{
  geometryregister.Append (new STLGeometryRegister);

  Tcl_CreateCommand (interp, "Ng_SetSTLParameters", Ng_SetSTLParameters,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);
  
  Tcl_CreateCommand (interp, "Ng_STLDoctor", Ng_STLDoctor,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

  Tcl_CreateCommand (interp, "Ng_STLInfo", Ng_STLInfo,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);
  
  Tcl_CreateCommand (interp, "Ng_STLCalcLocalH", Ng_STLCalcLocalH,
		     (ClientData)NULL,
		     (Tcl_CmdDeleteProc*) NULL);

  
  return TCL_OK;
}
