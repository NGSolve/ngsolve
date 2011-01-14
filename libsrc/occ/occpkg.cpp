#ifdef OCCGEOMETRY

#include <mystdlib.h>
#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <occgeom.hpp>


#include <incvis.hpp>
#include <visual.hpp>

#include "../meshing/bcfunctions.hpp"

#include "vsocc.hpp"


extern "C" int Ng_occ_Init (Tcl_Interp * interp);



namespace netgen
{
  extern NetgenGeometry * ng_geometry;
  extern AutoPtr<Mesh> mesh;
 
  char * err_needsoccgeometry = (char*) "This operation needs an OCC geometry";
  extern char * err_needsmesh;
  extern char * err_jobrunning;



                          
  class OCCGeometryRegister : public GeometryRegister
  {
  public:
    virtual NetgenGeometry * Load (string filename) const;
    virtual VisualScene * GetVisualScene (const NetgenGeometry * geom) const;

    virtual void SetParameters (Tcl_Interp * interp) 
    {
      occparam.resthcloseedgefac =
	atof (Tcl_GetVar (interp, "::stloptions.resthcloseedgefac", 0));
      occparam.resthcloseedgeenable =
	atoi (Tcl_GetVar (interp, "::stloptions.resthcloseedgeenable", 0));
    }
  };




  int Ng_SetOCCVisParameters  (ClientData clientData,
			       Tcl_Interp * interp,
			       int argc, tcl_const char *argv[])
  {
#ifdef OCCGEOMETRY
    int showvolume;
    OCCGeometry * occgeometry = dynamic_cast<OCCGeometry*> (ng_geometry);

    showvolume = atoi (Tcl_GetVar (interp, "::occoptions.showvolumenr", 0));

    if (occgeometry)
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




  int Ng_GetOCCData (ClientData clientData,
		     Tcl_Interp * interp,
		     int argc, tcl_const char *argv[])
  {
#ifdef OCCGEOMETRY
    OCCGeometry * occgeometry = dynamic_cast<OCCGeometry*> (ng_geometry);

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
    OCCGeometry * occgeometry = dynamic_cast<OCCGeometry*> (ng_geometry);

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
  /*
  void OCCConstructGeometry (OCCGeometry & geom);

  int Ng_OCCConstruction (ClientData clientData,
			  Tcl_Interp * interp,
			  int argc, tcl_const char *argv[])
  {
    if (occgeometry)
      OCCConstructGeometry (*occgeometry);
    return TCL_OK;
  }
  */
#endif




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

    OCCGeometry * occgeometry = dynamic_cast<OCCGeometry*> (ng_geometry);
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



  // Philippose - 25/07/2010
  // TCL interface function for extracting and eventually 
  // setting or editing the current colours present in the mesh
  int Ng_CurrentFaceColours (ClientData clientData,
                             Tcl_Interp * interp,
                             int argc, tcl_const char *argv[])
  {
     if(argc < 1)
     {
        Tcl_SetResult (interp, (char *)"Ng_GetCurrentFaceColours needs arguments", TCL_STATIC);
        return TCL_ERROR;
     }

     if(!mesh.Ptr())
     {
        Tcl_SetResult (interp, (char *)"Ng_GetCurrentFaceColours: Valid netgen mesh required...please mesh the Geometry first", TCL_STATIC);
	     return TCL_ERROR;
     }

     if(strcmp(argv[1], "getcolours") == 0)
     {
        stringstream outVar;
        Array<Vec3d> face_colours;
        GetFaceColours(*mesh, face_colours);

        for(int i = 0; i < face_colours.Size();i++)
        {
           outVar << "{ " << face_colours[i].X(1)
                  << " "  << face_colours[i].X(2)
                  << " "  << face_colours[i].X(3)
                  << " } ";
        }

        tcl_const char * valuevar = argv[2];
        Tcl_SetVar  (interp, valuevar, (char*)outVar.str().c_str(), 0);
     }

     if(strcmp(argv[1], "showalso") == 0)
     {
        Array<Vec3d> face_colours;
        GetFaceColours(*mesh,face_colours);

        int colourind = atoi (argv[2]);

        for(int i = 1; i <= mesh->GetNFD(); i++)
        {
           Array<SurfaceElementIndex> surfElems;
           mesh->GetSurfaceElementsOfFace(i,surfElems);

           if(ColourMatch(face_colours[colourind],mesh->GetFaceDescriptor(i).SurfColour()))
           {
              for(int j = 0; j < surfElems.Size(); j++)
              {
                 mesh->SurfaceElement(surfElems[j]).Visible(1);
              }
           }
        }

        mesh->SetNextTimeStamp();
     }

     if(strcmp(argv[1], "hidealso") == 0)
     {
        Array<Vec3d> face_colours;
        GetFaceColours(*mesh,face_colours);

        int colourind = atoi (argv[2]);

        for(int i = 1; i <= mesh->GetNFD(); i++)
        {
           Array<SurfaceElementIndex> surfElems;
           mesh->GetSurfaceElementsOfFace(i,surfElems);

           if(ColourMatch(face_colours[colourind],mesh->GetFaceDescriptor(i).SurfColour()))
           {
              for(int j = 0; j < surfElems.Size(); j++)
              {
                 mesh->SurfaceElement(surfElems[j]).Visible(0);
              }
           }
        }

        mesh->SetNextTimeStamp();
     }

     if(strcmp(argv[1], "showonly") == 0)
     {
        Array<Vec3d> face_colours;
        GetFaceColours(*mesh,face_colours);

        int colourind = atoi (argv[2]);

        for(int i = 1; i <= mesh->GetNFD(); i++)
        {
           Array<SurfaceElementIndex> surfElems;
           mesh->GetSurfaceElementsOfFace(i,surfElems);

           if(ColourMatch(face_colours[colourind],mesh->GetFaceDescriptor(i).SurfColour()))
           {
              for(int j = 0; j < surfElems.Size(); j++)
              {
                 mesh->SurfaceElement(surfElems[j]).Visible(1);
              }
           }
           else
           {
              for(int j = 0; j < surfElems.Size(); j++)
              {
                 mesh->SurfaceElement(surfElems[j]).Visible(0);
              }
           }
        }

        mesh->SetNextTimeStamp();
     }

     if(strcmp(argv[1], "hideonly") == 0)
     {
        Array<Vec3d> face_colours;
        GetFaceColours(*mesh,face_colours);

        int colourind = atoi (argv[2]);

        for(int i = 1; i <= mesh->GetNFD(); i++)
        {
           Array<SurfaceElementIndex> surfElems;
           mesh->GetSurfaceElementsOfFace(i,surfElems);

           if(ColourMatch(face_colours[colourind],mesh->GetFaceDescriptor(i).SurfColour()))
           {
              for(int j = 0; j < surfElems.Size(); j++)
              {
                 mesh->SurfaceElement(surfElems[j]).Visible(0);
              }
           }
           else
           {
              for(int j = 0; j < surfElems.Size(); j++)
              {
                 mesh->SurfaceElement(surfElems[j]).Visible(1);
              }
           }
        }

        mesh->SetNextTimeStamp();
     }

     if(strcmp(argv[1], "showall") == 0)
     {
        for(int i = 1; i <= mesh->GetNSE(); i++)
        {
           mesh->SurfaceElement(i).Visible(1);
        }

        mesh->SetNextTimeStamp();
     }

     if(strcmp(argv[1], "hideall") == 0)
     {
        for(int i = 1; i <= mesh->GetNSE(); i++)
        {
           mesh->SurfaceElement(i).Visible(0);
        }

        mesh->SetNextTimeStamp();
     }

     return TCL_OK;
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


  int Ng_SetOCCParameters  (ClientData clientData,
			    Tcl_Interp * interp,
			    int argc, tcl_const char *argv[])
  {
    OCCGeometryRegister reg;
    reg.SetParameters (interp);
    /*
    occparam.resthcloseedgefac =
      atof (Tcl_GetVar (interp, "::stloptions.resthcloseedgefac", 0));

    occparam.resthcloseedgeenable =
      atoi (Tcl_GetVar (interp, "::stloptions.resthcloseedgeenable", 0));
    */
    return TCL_OK;
  }




  NetgenGeometry *  OCCGeometryRegister :: Load (string filename) const
  {
    const char * lgfilename = filename.c_str();


    /*
    if (strcmp (&cfilename[strlen(cfilename)-3], "geo") == 0)
      {
	PrintMessage (1, "Load OCCG geometry file ", cfilename);
	
	extern OCCGeometry * ParseOCCG (istream & istr);

	ifstream infile(cfilename);

	OCCGeometry * hgeom = ParseOCCG (infile);
	if (!hgeom)
	  throw NgException ("geo-file should start with 'algebraic3d'");

	hgeom -> FindIdenticSurfaces(1e-8 * hgeom->MaxSize()); 
	return hgeom;
      }
    */


    if ((strcmp (&lgfilename[strlen(lgfilename)-4], "iges") == 0) ||
	(strcmp (&lgfilename[strlen(lgfilename)-3], "igs") == 0) ||
	(strcmp (&lgfilename[strlen(lgfilename)-3], "IGS") == 0) ||
	(strcmp (&lgfilename[strlen(lgfilename)-4], "IGES") == 0))
      {
	PrintMessage (1, "Load IGES geometry file ", lgfilename);
	OCCGeometry * occgeometry = LoadOCC_IGES (lgfilename);
	return occgeometry;
      }

    else if ((strcmp (&lgfilename[strlen(lgfilename)-4], "step") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-3], "stp") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-3], "STP") == 0) ||
		     (strcmp (&lgfilename[strlen(lgfilename)-4], "STEP") == 0))
      {
	PrintMessage (1, "Load STEP geometry file ", lgfilename);
	OCCGeometry * occgeometry = LoadOCC_STEP (lgfilename);
	return occgeometry;    
      }
    else if ((strcmp (&lgfilename[strlen(lgfilename)-4], "brep") == 0) ||
	     (strcmp (&lgfilename[strlen(lgfilename)-4], "Brep") == 0) ||
	     (strcmp (&lgfilename[strlen(lgfilename)-4], "BREP") == 0))
      {
	PrintMessage (1, "Load BREP geometry file ", lgfilename);
	OCCGeometry * occgeometry = LoadOCC_BREP (lgfilename);
	return occgeometry;
      }
    
    return NULL;
  }


  static VisualSceneOCCGeometry vsoccgeom;

  VisualScene * OCCGeometryRegister :: GetVisualScene (const NetgenGeometry * geom) const
  {
    OCCGeometry * geometry = dynamic_cast<OCCGeometry*> (ng_geometry);
    if (geometry)
      {
	vsoccgeom.SetGeometry (geometry);
	return &vsoccgeom;
      }
    return NULL;
  }



}



using namespace netgen;

int Ng_occ_Init (Tcl_Interp * interp)
{
  geometryregister.Append (new OCCGeometryRegister);


    Tcl_CreateCommand (interp, "Ng_SetOCCVisParameters",
		       Ng_SetOCCVisParameters,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_GetOCCData",
		       Ng_GetOCCData,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);

    /*
#ifdef OCCGEOMETRY
    Tcl_CreateCommand (interp, "Ng_OCCConstruction",
		       Ng_OCCConstruction,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);
#endif
    */

    Tcl_CreateCommand (interp, "Ng_OCCCommand",
		       Ng_OCCCommand,
		       (ClientData)NULL,
		       (Tcl_CmdDeleteProc*) NULL);


    Tcl_CreateCommand (interp, "Ng_SetOCCParameters", Ng_SetOCCParameters,
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

    // Philippose - 25/07/2010
    // Register the TCL Interface Command for handling the face colours 
    // present in the mesh
    Tcl_CreateCommand(interp, "Ng_CurrentFaceColours", Ng_CurrentFaceColours,
                      (ClientData)NULL,
                      (Tcl_CmdDeleteProc*) NULL);


  return TCL_OK;
}

#endif

