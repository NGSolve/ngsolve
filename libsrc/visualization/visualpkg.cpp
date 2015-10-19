#include <mystdlib.h>
// #include <incopengl.hpp>
#include <inctcl.hpp>


#include <myadt.hpp>
#include <meshing.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>


// #include <parallel.hpp>
#include <visual.hpp>

#include <limits>






namespace netgen

{

  /*
  */








  int Ng_Vis_Set (ClientData clientData,
                  Tcl_Interp * interp,
                  int argc, tcl_const char *argv[])

  {
    if (argc >= 2)
      {
        if (strcmp (argv[1], "parameters") == 0)
          {
            vssolution.imag_part = 
              atoi (Tcl_GetVar (interp, "::visoptions.imaginary", TCL_GLOBAL_ONLY));      
            vssolution.usetexture = 
              atoi (Tcl_GetVar (interp, "::visoptions.usetexture", TCL_GLOBAL_ONLY));
            if (atoi (Tcl_GetVar (interp, "::visoptions.redrawperiodic", TCL_GLOBAL_ONLY)))
              vssolution.usetexture = 2;
                
            vssolution.invcolor = 
              atoi (Tcl_GetVar (interp, "::visoptions.invcolor", TCL_GLOBAL_ONLY));       

            vssolution.clipsolution = 0;

            if (strcmp (Tcl_GetVar (interp, "::visoptions.clipsolution", TCL_GLOBAL_ONLY), 
                        "scal") == 0)
              vssolution.clipsolution = 1;
            if (strcmp (Tcl_GetVar (interp, "::visoptions.clipsolution", TCL_GLOBAL_ONLY), 
                        "vec") == 0)
              vssolution.clipsolution = 2;
            
            tcl_const char * scalname =  
              Tcl_GetVar (interp, "::visoptions.scalfunction", TCL_GLOBAL_ONLY);
            tcl_const char * vecname = 
              Tcl_GetVar (interp, "::visoptions.vecfunction", TCL_GLOBAL_ONLY);
            tcl_const char * fieldlines_vecname = 
              Tcl_GetVar (interp, "::visoptions.fieldlinesvecfunction", TCL_GLOBAL_ONLY);
                
          
            vssolution.scalfunction = -1;
            vssolution.vecfunction = -1;
            vssolution.fieldlines_vecfunction = -1;

            int pointpos; // SZ 
            const char * pch = strchr(scalname,':');
            pointpos = int(pch-scalname+1);
	    
            for (int i = 0; i < vssolution.soldata.Size(); i++)
              {
                if ( (strlen (vssolution.soldata[i]->name) == size_t(pointpos-1)) &&
                     (strncmp (vssolution.soldata[i]->name, scalname, pointpos-1) == 0) )
                  {
                    vssolution.scalfunction = i;
                    vssolution.scalcomp = atoi (scalname + pointpos);
		    if ( vssolution.scalcomp > vssolution.soldata[i]->components )
                      vssolution.scalcomp = 1;
		    char newscalname[100];
		    for ( int ii = 0; ii < pointpos; ii++ )
		      newscalname[ii] = scalname[ii];
		    newscalname[pointpos] = ':';
		    sprintf (newscalname+pointpos, "%i", vssolution.scalcomp);

                    if (strcmp (scalname, newscalname) != 0)
                      Tcl_SetVar ( interp, "::visoptions.scalfunction", newscalname, TCL_GLOBAL_ONLY );
		    scalname = Tcl_GetVar (interp, "::visoptions.scalfunction", TCL_GLOBAL_ONLY);
                  }
                if (strcmp (vssolution.soldata[i]->name, vecname) == 0)
		  vssolution.vecfunction = i;

                if (strcmp (vssolution.soldata[i]->name, fieldlines_vecname) == 0)
		  vssolution.fieldlines_vecfunction = i;
              }

            if(vssolution.fieldlines_vecfunction != -1 &&
               vssolution.vecfunction == -1)
              {
                //cout << "WARNING: Setting vector function in Visualization toolbox to value from Fieldlines toolbox!" << endl;
                vssolution.vecfunction = vssolution.fieldlines_vecfunction;
              }
               
	    // reset visoptions.scalfunction and visoptions.vecfunction if not avialable 
	    if ( vssolution.scalfunction == -1 && strcmp (scalname, "none") != 0)
              Tcl_SetVar ( interp, "::visoptions.scalfunction", "none", TCL_GLOBAL_ONLY );
	    if ( vssolution.vecfunction == -1  && strcmp (vecname, "none") != 0)
              Tcl_SetVar ( interp, "::visoptions.vecfunction", "none", TCL_GLOBAL_ONLY );

            tcl_const char * evalname = 
              Tcl_GetVar (interp, "::visoptions.evaluate", TCL_GLOBAL_ONLY);
          
            if (strcmp(evalname, "abs") == 0) vssolution.evalfunc = VisualSceneSolution::FUNC_ABS;
            if (strcmp(evalname, "abstens") == 0) vssolution.evalfunc = VisualSceneSolution::FUNC_ABS_TENSOR;
            if (strcmp(evalname, "mises") == 0) vssolution.evalfunc = VisualSceneSolution::FUNC_MISES;
            if (strcmp(evalname, "main") == 0) vssolution.evalfunc = VisualSceneSolution::FUNC_MAIN;

            vssolution.gridsize = 
              atoi (Tcl_GetVar (interp, "::visoptions.gridsize", TCL_GLOBAL_ONLY));

            vssolution.xoffset = 
              atof (Tcl_GetVar (interp, "::visoptions.xoffset", TCL_GLOBAL_ONLY));

            //    cout << "x-offset:" << vssolution.xoffset << endl;

            vssolution.yoffset = 
              atof (Tcl_GetVar (interp, "::visoptions.yoffset", TCL_GLOBAL_ONLY));

            vssolution.autoscale = 
              atoi (Tcl_GetVar (interp, "::visoptions.autoscale", TCL_GLOBAL_ONLY));


            /*
              vssolution.linear_colors = 
              atoi (Tcl_GetVar (interp, "::visoptions.lineartexture", TCL_GLOBAL_ONLY));
            */
            vssolution.logscale = 
              atoi (Tcl_GetVar (interp, "::visoptions.logscale", TCL_GLOBAL_ONLY));

            vssolution.mminval = 
              atof (Tcl_GetVar (interp, "::visoptions.mminval", TCL_GLOBAL_ONLY));
            vssolution.mmaxval = 
              atof (Tcl_GetVar (interp, "::visoptions.mmaxval", TCL_GLOBAL_ONLY));

            vssolution.showclipsolution = 
              atoi (Tcl_GetVar (interp, "::visoptions.showclipsolution", TCL_GLOBAL_ONLY));
            vssolution.showsurfacesolution = 
              atoi (Tcl_GetVar (interp, "::visoptions.showsurfacesolution", TCL_GLOBAL_ONLY));
            vssolution.lineartexture = 
              atoi (Tcl_GetVar (interp, "::visoptions.lineartexture", TCL_GLOBAL_ONLY));
            vssolution.numtexturecols = 
              atoi (Tcl_GetVar (interp, "::visoptions.numtexturecols", TCL_GLOBAL_ONLY));

            vssolution.multidimcomponent = 
              atoi (Tcl_GetVar (interp, "::visoptions.multidimcomponent", TCL_GLOBAL_ONLY));

	    vssolution.drawpointcurves = 
	      atoi (Tcl_GetVar (interp, "::visoptions.drawpointcurves", TCL_GLOBAL_ONLY));	      

            vssolution.draw_fieldlines = 
	      atoi (Tcl_GetVar (interp, "::visoptions.drawfieldlines", TCL_GLOBAL_ONLY));
            vssolution.num_fieldlines = 
              atoi (Tcl_GetVar (interp, "::visoptions.numfieldlines", TCL_GLOBAL_ONLY));
            vssolution.fieldlines_randomstart =
              atoi (Tcl_GetVar (interp, "::visoptions.fieldlinesrandomstart", TCL_GLOBAL_ONLY));

            vssolution.fieldlines_reltolerance =
              atof (Tcl_GetVar (interp, "::visoptions.fieldlinestolerance", TCL_GLOBAL_ONLY));

            if (strcmp (Tcl_GetVar (interp, "::visoptions.fieldlinesrktype", TCL_GLOBAL_ONLY), 
                        "euler") == 0)
              vssolution.fieldlines_rktype = 0;
            else if (strcmp (Tcl_GetVar (interp, "::visoptions.fieldlinesrktype", TCL_GLOBAL_ONLY), 
                             "eulercauchy") == 0)
              vssolution.fieldlines_rktype = 1;
            else if (strcmp (Tcl_GetVar (interp, "::visoptions.fieldlinesrktype", TCL_GLOBAL_ONLY), 
                             "simpson") == 0)
              vssolution.fieldlines_rktype = 2;
            else if (strcmp (Tcl_GetVar (interp, "::visoptions.fieldlinesrktype", TCL_GLOBAL_ONLY), 
                             "crungekutta") == 0)
              vssolution.fieldlines_rktype = 3;


            vssolution.fieldlines_rellength =
              atof (Tcl_GetVar (interp, "::visoptions.fieldlineslength", TCL_GLOBAL_ONLY));

            vssolution.fieldlines_maxpoints = 
              atoi (Tcl_GetVar (interp, "::visoptions.fieldlinesmaxpoints", TCL_GLOBAL_ONLY));

            vssolution.fieldlines_relthickness =
              atof (Tcl_GetVar (interp, "::visoptions.fieldlinesthickness", TCL_GLOBAL_ONLY));


            vssolution.fieldlines_fixedphase = 
              (atoi (Tcl_GetVar (interp, "::visoptions.fieldlinesonlyonephase", TCL_GLOBAL_ONLY)) != 0);

            if(vssolution.fieldlines_fixedphase)
              vssolution.fieldlines_phase =
                atof (Tcl_GetVar (interp, "::visoptions.fieldlinesphase", TCL_GLOBAL_ONLY));


            if (strcmp (Tcl_GetVar (interp, "::visoptions.fieldlinesstartarea", TCL_GLOBAL_ONLY), 
                        "box") == 0)
              vssolution.fieldlines_startarea  = 0;
            else if (strcmp (Tcl_GetVar (interp, "::visoptions.fieldlinesstartarea", TCL_GLOBAL_ONLY), 
                             "file") == 0)
              vssolution.fieldlines_startarea  = 1;
            else if (strcmp (Tcl_GetVar (interp, "::visoptions.fieldlinesstartarea", TCL_GLOBAL_ONLY), 
                             "face") == 0)
              vssolution.fieldlines_startarea  = 2;

                
            if (vssolution.fieldlines_startarea == 0)
              {
                vssolution.fieldlines_startarea_parameter.SetSize(6);
                vssolution.fieldlines_startarea_parameter[0] = atof (Tcl_GetVar (interp, "::visoptions.fieldlinesstartareap1x", TCL_GLOBAL_ONLY));
                vssolution.fieldlines_startarea_parameter[1] = atof (Tcl_GetVar (interp, "::visoptions.fieldlinesstartareap1y", TCL_GLOBAL_ONLY));
                vssolution.fieldlines_startarea_parameter[2] = atof (Tcl_GetVar (interp, "::visoptions.fieldlinesstartareap1z", TCL_GLOBAL_ONLY));
                vssolution.fieldlines_startarea_parameter[3] = atof (Tcl_GetVar (interp, "::visoptions.fieldlinesstartareap2x", TCL_GLOBAL_ONLY));
                vssolution.fieldlines_startarea_parameter[4] = atof (Tcl_GetVar (interp, "::visoptions.fieldlinesstartareap2y", TCL_GLOBAL_ONLY));
                vssolution.fieldlines_startarea_parameter[5] = atof (Tcl_GetVar (interp, "::visoptions.fieldlinesstartareap2z", TCL_GLOBAL_ONLY));
              }
            else if (vssolution.fieldlines_startarea == 1)
              {
                vssolution.fieldlines_filename = Tcl_GetVar (interp, "::visoptions.fieldlinesfilename", TCL_GLOBAL_ONLY);
              }
            else if (vssolution.fieldlines_startarea == 2)
              {
                vssolution.fieldlines_startface = atoi (Tcl_GetVar (interp, "::visoptions.fieldlinesstartface", TCL_GLOBAL_ONLY));
              }

          
            vssolution.deform =
              atoi (Tcl_GetVar (interp, "::visoptions.deformation", TCL_GLOBAL_ONLY));
            vssolution.scaledeform =
              atof (Tcl_GetVar (interp, "::visoptions.scaledeform1", TCL_GLOBAL_ONLY)) *
              atof (Tcl_GetVar (interp, "::visoptions.scaledeform2", TCL_GLOBAL_ONLY));


            if (atoi (Tcl_GetVar (interp, "::visoptions.isolines", TCL_GLOBAL_ONLY)))
              vssolution.numisolines = atoi (Tcl_GetVar (interp, "::visoptions.numiso", TCL_GLOBAL_ONLY));
            else
              vssolution.numisolines = 0;
            vssolution.draw_isosurface = 
              atoi (Tcl_GetVar (interp, "::visoptions.isosurf", TCL_GLOBAL_ONLY));

            vssolution.SetSubdivision(atoi (Tcl_GetVar (interp, "::visoptions.subdivisions", TCL_GLOBAL_ONLY)));

            vssolution.UpdateSolutionTimeStamp();
          }
      
        if (strcmp (argv[1], "parametersrange") == 0)
          {
            vssolution.invcolor = 
              atoi (Tcl_GetVar (interp, "::visoptions.invcolor", TCL_GLOBAL_ONLY));       
            vssolution.mminval = 
              atof (Tcl_GetVar (interp, "::visoptions.mminval", TCL_GLOBAL_ONLY));
            vssolution.mmaxval = 
              atof (Tcl_GetVar (interp, "::visoptions.mmaxval", TCL_GLOBAL_ONLY));
            vssolution.lineartexture = 
              atoi (Tcl_GetVar (interp, "::visoptions.lineartexture", TCL_GLOBAL_ONLY));
            vssolution.numtexturecols = 
              atoi (Tcl_GetVar (interp, "::visoptions.numtexturecols", TCL_GLOBAL_ONLY));

            if (vssolution.usetexture == 0 || vssolution.logscale)
              vssolution.UpdateSolutionTimeStamp();
          }


        if (argc >= 3 && strcmp (argv[1], "time") == 0)
          {
            vssolution.time = double (atoi (argv[2])) / 1000;
         
            vssolution.timetimestamp = NextTimeStamp();
            cout << "\rtime = " << vssolution.time << "    " << flush;
          }

      }

    
    vsmesh.SetClippingPlane ();  // for computing parameters
    vssolution.SetClippingPlane ();  // for computing parameters
    glDisable(GL_CLIP_PLANE0);

#ifdef PARALLELGL
    vsmesh.Broadcast ();
#endif    


    return TCL_OK;
  }

  int Ng_Vis_Field (ClientData clientData,
                    Tcl_Interp * interp,
                    int argc, tcl_const char *argv[])
  {
    int i;
    char buf[1000];
    buf[0] = 0;

    if (argc >= 2)
      {
        if (strcmp (argv[1], "setfield") == 0)
          {
            if (argc < 3)
              return TCL_ERROR;

            for (i = 0; i < vssolution.GetNSolData(); i++)
              if (strcmp (vssolution.GetSolData(i)->name, argv[2]) == 0)
                {
                  cout << "found soldata " << i << endl;
                }
          }

        if (strcmp (argv[1], "getnfieldnames") == 0)
          {
            sprintf (buf, "%d", vssolution.GetNSolData());
          }
      
        if (strcmp (argv[1], "getfieldname") == 0)
          {
            sprintf (buf, "%s", vssolution.GetSolData(atoi(argv[2])-1)->name);
          }

        if (strcmp (argv[1], "iscomplex") == 0)
          {
            sprintf (buf, "%d", vssolution.GetSolData(atoi(argv[2])-1)->iscomplex);
          }

        if (strcmp (argv[1], "getfieldcomponents") == 0)
          {
            sprintf (buf, "%d", vssolution.GetSolData(atoi(argv[2])-1)->components);
          }

      
        if (strcmp (argv[1], "getfieldnames") == 0)
          {
            for (i = 0; i < vssolution.GetNSolData(); i++)
              {
                strcat (buf, vssolution.GetSolData(i)->name);
                strcat (buf, " ");
              }
            strcat (buf, "var1 var2 var3");
            Tcl_SetResult (interp, buf, TCL_STATIC);
          }

        if (strcmp (argv[1], "setcomponent") == 0)
          {
            cout << "set component " << argv[2] << endl;
          }

        if (strcmp (argv[1], "getactivefield") == 0)
          {
            sprintf (buf, "1");
          }

        if (strcmp (argv[1], "getdimension") == 0)
          {
            auto mesh = vssolution.GetMesh();
            sprintf (buf, "%d", mesh->GetDimension());
          }
      }

    Tcl_SetResult (interp, buf, TCL_STATIC);
    return TCL_OK;
  }

  VisualSceneMeshDoctor vsmeshdoc;
  DLL_HEADER extern shared_ptr<Mesh> mesh;

  int Ng_MeshDoctor(ClientData clientData,
	  Tcl_Interp * interp,
	  int argc, tcl_const char *argv[])
  {
	  cout << "Mesh Doctor:" << endl;
	  int i;
	  for (i = 0; i < argc; i++)
		  cout << argv[i] << " ";
	  cout << endl;

	  meshdoctor.active =
		  atoi(Tcl_GetVar(interp, "::meshdoctor.active", 0));


	  if (argc >= 2)
	  {
		  if (strcmp(argv[1], "markedgedist") == 0)
		  {
			  vsmeshdoc.SetMarkEdgeDist(atoi(argv[2]));
		  }

		  if (strcmp(argv[1], "deletemarkedsegments") == 0)
		  {
			  for (i = 1; i <= mesh->GetNSeg(); i++)
				  if (vsmeshdoc.IsSegmentMarked(i))
					  mesh->DeleteSegment(i);

			  //	  for (i = 1; i <= mesh->GetNSE(); i++)
			  //	    mesh->SurfaceElement(i).SetIndex (1);
			  mesh->Compress();
		  }
	  }


	  vsmeshdoc.UpdateTables();
	  vsmeshdoc.BuildScene();
	  return TCL_OK;
  }


  extern "C" int Ng_Vis_Init (Tcl_Interp * interp);

  int Ng_Vis_Init (Tcl_Interp * interp)
  {
    Tcl_CreateCommand (interp, "Ng_Vis_Set", Ng_Vis_Set,
                       (ClientData)NULL,
                       (Tcl_CmdDeleteProc*) NULL);

    Tcl_CreateCommand (interp, "Ng_Vis_Field", Ng_Vis_Field,
                       (ClientData)NULL,
                       (Tcl_CmdDeleteProc*) NULL);


    return TCL_OK;
  }




}

