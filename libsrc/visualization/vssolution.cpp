#include <mystdlib.h>

#include <myadt.hpp>
#include <meshing.hpp>
#include <csg.hpp>
#include <stlgeom.hpp>


// #include <parallel.hpp>
#include <visual.hpp>

#include <limits>


namespace netgen
{

  VisualSceneSolution vssolution;
  // extern shared_ptr<Mesh> mesh;
  extern VisualSceneMesh vsmesh;


  void AddUserVisualizationObject (UserVisualizationObject * vis)
  {
    vssolution.AddUserVisualizationObject (vis);
  }


  VisualSceneSolution :: SolData :: SolData ()
    : name (0), data (0), solclass(0)
  { ; }

  VisualSceneSolution :: SolData :: ~SolData ()
  {
    delete [] name;
    delete data;
    delete solclass;
  }

  
  VisualSceneSolution :: VisualSceneSolution ()
    : VisualScene()
  {
    surfellist = 0;
    linelist = 0;
    element1dlist = 0;
    clipplanelist_scal = 0;
    clipplanelist_vec = 0;
    isolinelist = 0;
    clipplane_isolinelist = 0;
    surface_vector_list = 0;
    isosurface_list = 0;

    fieldlineslist = 0;
    pointcurvelist = 0;

    num_fieldlineslists = 0;


    surfeltimestamp = GetTimeStamp();
    surfellinetimestamp = GetTimeStamp();
    clipplanetimestamp = GetTimeStamp();
    solutiontimestamp = GetTimeStamp();
    fieldlinestimestamp = GetTimeStamp();
    pointcurve_timestamp = GetTimeStamp();
    surface_vector_timestamp = GetTimeStamp();
    isosurface_timestamp = GetTimeStamp();
    timetimestamp = GetTimeStamp();
    // AddVisualizationScene ("solution", &vssolution);
  }
  
  VisualSceneSolution :: ~VisualSceneSolution ()
  {
    ClearSolutionData();
  }

  /*
  void VisualSceneSolution :: SetMesh (shared_ptr<Mesh> amesh)
  { 
    wp_mesh = amesh; 
  }
  */

  void VisualSceneSolution :: AddSolutionData (SolData * sd)
  {
    shared_ptr<Mesh> mesh = GetMesh();

    NgLock meshlock1 (mesh->MajorMutex(), 1);
    int funcnr = -1;
    for (int i = 0; i < soldata.Size(); i++)
      {
        if (strcmp (soldata[i]->name, sd->name) == 0)
          {
            delete soldata[i];
            soldata[i] = sd;
            funcnr = i;
            break;
          }
      }

    if (funcnr == -1)
      {
        soldata.Append (sd);
        funcnr = soldata.Size()-1;
      }
    
    SolData * nsd = soldata[funcnr]; 

    nsd->size = 0;
    if (mesh)
      {
        switch (nsd->soltype)
          {
          case SOL_NODAL: nsd->size = mesh->GetNV(); break;
          case SOL_ELEMENT: nsd->size = mesh->GetNE(); break;
          case SOL_SURFACE_ELEMENT: nsd->size = mesh->GetNSE(); break;
          case SOL_NONCONTINUOUS: 
            {
              switch (nsd->order)
                {
                case 0: nsd->size =      mesh->GetNE(); break;
                case 1: nsd->size =  6 * mesh->GetNE(); break;
                case 2: nsd->size = 18 * mesh->GetNE(); break;
                }
              break;
            }
          case SOL_SURFACE_NONCONTINUOUS: 
            {
              switch (nsd->order)
                {
                case 0: nsd->size =     mesh->GetNSE(); break;
                case 1: nsd->size = 4 * mesh->GetNSE(); break;
                case 2: nsd->size = 9 * mesh->GetNSE(); break;
                }
              break;
            }
          default:
            nsd->size = 0;
          }
        solutiontimestamp = NextTimeStamp();
      }
  }

  
  void VisualSceneSolution :: ClearSolutionData ()
  {
    for (int i = 0; i < soldata.Size(); i++)
      delete soldata[i];
    soldata.SetSize (0);
  }

  void VisualSceneSolution :: UpdateSolutionTimeStamp ()
  {
    solutiontimestamp = NextTimeStamp();
  }
    
  VisualSceneSolution::SolData * VisualSceneSolution :: GetSolData (int i)
  { 
    if (i >= 0 && i < soldata.Size())
      return soldata[i];
    else 
      return NULL;
  }
  



  void VisualSceneSolution :: SaveSolutionData (const char * filename) 
  {
    shared_ptr<Mesh> mesh = GetMesh();

    PrintMessage (1, "Write solution data to file ", filename);


    if (strcmp (&filename[strlen(filename)-3], "sol") == 0)
      {
        ofstream ost(filename);
        for (int i = 0; i < soldata.Size(); i++)
          {
            const SolData & sol = *soldata[i];
      
            ost << "solution " 
                << sol.name
                << " -size=" << sol.size 
                << " -components=" << sol.components
                << " -order=" << sol.order;
            if (sol.iscomplex)
              ost << " -complex";
      
            switch (sol.soltype)
              {
              case SOL_NODAL:
                ost << " -type=nodal"; break;
              case SOL_ELEMENT:
                ost << " -type=element"; break;
              case SOL_SURFACE_ELEMENT:
                ost << " -type=surfaceelement"; break;
              case SOL_NONCONTINUOUS:
                ost << " -type=noncontinuous"; break;
              case SOL_SURFACE_NONCONTINUOUS:
                ost << " -type=surfacenoncontinuous"; break;
              default:
                cerr << "save solution data, case not handeld" << endl;
              }
      
            ost << endl;
            for (int j = 0; j < sol.size; j++)
              {
                for (int k = 0; k < sol.components; k++)
                  ost << sol.data[j*sol.dist+k] << " ";
                ost << "\n";
              }
          }
      }


    if (strcmp (&filename[strlen(filename)-3], "vtk") == 0)
      {
        string surf_fn = filename;
        surf_fn.erase (strlen(filename)-4);
        surf_fn += "_surf.vtk";

        cout << "surface mesh = " << surf_fn << endl;
        
        ofstream surf_ost(surf_fn.c_str());

        surf_ost << "# vtk DataFile Version 1.0\n"
		 << "NGSolve surface mesh\n"
		 << "ASCII\n"
		 << "DATASET UNSTRUCTURED_GRID\n\n";

        surf_ost << "POINTS " << mesh->GetNP() << " float\n";
        for (PointIndex pi = PointIndex::BASE; pi < mesh->GetNP()+PointIndex::BASE; pi++)
          {
            const MeshPoint & mp = (*mesh)[pi];
            surf_ost << mp(0) << " " << mp(1) << " " << mp(2) << "\n";
          }

        int cntverts = 0;
        for (SurfaceElementIndex sei = 0; sei < mesh->GetNSE(); sei++)
          cntverts += 1 + (*mesh)[sei].GetNP();

        surf_ost << "\nCELLS " << mesh->GetNSE() << " " << cntverts << "\n";
        for (SurfaceElementIndex sei = 0; sei < mesh->GetNSE(); sei++)
          {
            const Element2d & el = (*mesh)[sei];
            surf_ost << el.GetNP();
            for (int j = 0; j < el.GetNP(); j++)
              surf_ost << " " << el[j] - PointIndex::BASE;
            surf_ost << "\n";
          }
        surf_ost << "\nCELL_TYPES " << mesh->GetNSE() << "\n";
        for (SurfaceElementIndex sei = 0; sei < mesh->GetNSE(); sei++)
          {
            const Element2d & el = (*mesh)[sei];
            switch (el.GetType())
              {
              case QUAD: surf_ost << 9; break;
              case TRIG: surf_ost << 5; break;
              default:
                cerr << "not implemented 2378" << endl;
              }
            surf_ost << "\n";
          }


       
        ofstream ost(filename);

        ost << "# vtk DataFile Version 1.0\n"
            << "NGSolve solution\n"
            << "ASCII\n"
            << "DATASET UNSTRUCTURED_GRID\n\n";

        ost << "POINTS " << mesh->GetNP() << " float\n";
        for (PointIndex pi = PointIndex::BASE; pi < mesh->GetNP()+PointIndex::BASE; pi++)
          {
            const MeshPoint & mp = (*mesh)[pi];
            ost << mp(0) << " " << mp(1) << " " << mp(2) << "\n";
          }

        cntverts = 0;
        for (ElementIndex ei = 0; ei < mesh->GetNE(); ei++)
          cntverts += 1 + (*mesh)[ei].GetNP();

        ost << "\nCELLS " << mesh->GetNE() << " " << cntverts << "\n";
        for (ElementIndex ei = 0; ei < mesh->GetNE(); ei++)
          {
            const Element & el = (*mesh)[ei];
            ost << el.GetNP();
            for (int j = 0; j < el.GetNP(); j++)
              ost << " " << el[j] - PointIndex::BASE;
            ost << "\n";
          }
        ost << "\nCELL_TYPES " << mesh->GetNE() << "\n";
        for (ElementIndex ei = 0; ei < mesh->GetNE(); ei++)
          {
            const Element & el = (*mesh)[ei];
            switch (el.GetType())
              {
              case TET: ost << 10; break;
              default:
                cerr << "not implemented 67324" << endl;
              }
            ost << "\n";
          }


        ost << "CELL_DATA " << mesh->GetNE() << "\n";
        for (int i = 0; i < soldata.Size(); i++)
          {
            ost << "VECTORS bfield float\n";
            SolutionData & sol = *(soldata[i] -> solclass);
            double values[3];

            for (int elnr = 0; elnr < mesh->GetNE(); elnr++)
              {
                sol.GetValue (elnr, 0.25, 0.25, 0.25, values);
                ost << values[0] << " "  << values[1] << " "  << values[2] << "\n";
              }
          }

        /*
	  ost << "POINT_DATA " << mesh->GetNP() << "\n";
	  for (int i = 0; i < soldata.Size(); i++)
          {
	  ost << "VECTORS bfield float\n";
	  SolutionData & sol = *(soldata[i] -> solclass);
            
	  for (PointIndex pi = PointIndex::BASE; 
	  pi < mesh->GetNP()+PointIndex::BASE; pi++)
	  {
	  double values[3], sumvalues[3] = { 0, 0, 0 };

	  FlatArray<int> els = mesh->GetTopology().GetVertexElements(pi);

	  for (int j = 0; j < els.Size(); j++)
	  {
	  sol.GetValue (els[j]-1, 0.25, 0.25, 0.25, values);
	  for (int k = 0; k < 3; k++)
	  sumvalues[k] += values[k];
	  }
	  for (int k = 0; k < 3; k++)
	  sumvalues[k] /= els.Size();
                
	  ost << sumvalues[0] << " "  << sumvalues[1] << " "  << sumvalues[2] << "\n";
	  }
          }
        */
      } 
    
  }
  



  void VisualSceneSolution :: DrawScene ()
  {
    try
      {
    shared_ptr<Mesh> mesh = GetMesh();

    if (!mesh) 
      {
        VisualScene::DrawScene();      
        return;
      }

    // static NgLock mem_lock(mem_mutex);
    // mem_lock.Lock();

    NgLock meshlock1 (mesh->MajorMutex(), true);
    NgLock meshlock (mesh->Mutex(), true);

    BuildScene();

    CreateTexture (numtexturecols, lineartexture, 0.5, GL_MODULATE);

    glClearColor(backcolor, backcolor, backcolor, 1);
    // glClearColor(backcolor, backcolor, backcolor, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    SetLight();
    
    glPushMatrix();
    glMultMatrixd (transformationmat);

    glMatrixMode (GL_MODELVIEW); 
    
    glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
    
    glPolygonOffset (1, 1);

    glEnable (GL_POLYGON_OFFSET_FILL);

    glEnable (GL_COLOR_MATERIAL);

    if (usetexture)
      {
        SetTextureMode (usetexture);

        glMatrixMode (GL_TEXTURE);
        glLoadIdentity();
        
        if (usetexture == 1)
          {
            double hmax = maxval;
            double hmin = minval;
            if (invcolor) Swap (hmax, hmin);

            if (fabs (hmax - hmin) > 1e-30) 
              glScaled (1.0 / (hmin - hmax), 0, 0);
            else
              glScaled (1e30, 0, 0);
            
            glTranslatef (-hmax, 0, 0);
          }
        else
          {
            glTranslatef (0.5, 0, 0);
            glRotatef(360 * vssolution.time, 0, 0, -1);
            if (fabs (maxval) > 1e-10)
              glScalef(0.5/maxval, 0.5/maxval, 0.5/maxval);
            else
              glScalef (1e10, 1e10, 1e10);
          }
        glMatrixMode (GL_MODELVIEW);
      }




    if (vispar.drawfilledtrigs || vispar.drawtetsdomain > 0 || vispar.drawdomainsurf > 0)
      {
	// Change for Martin:

	// orig:
	SetClippingPlane ();  

	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);   
	// glEnable(GL_BLEND); 
	glDisable(GL_BLEND); 
	glCallList (surfellist);
	glDisable(GL_BLEND); 
	/*
	// transparent test ...
	glColor4f (1, 0, 0, 0.1);
	glEnable (GL_COLOR_MATERIAL);

	glDepthFunc(GL_GREATER); 
	glDepthMask(GL_FALSE); 
	// glBlendFunc(GL_ONE_MINUS_DST_ALPHA,GL_DST_ALPHA); 
	glBlendFunc(GL_ONE_MINUS_SRC_ALPHA,GL_SRC_ALPHA); 

	glCallList (surfellist);

	glDisable(GL_BLEND);
	glDepthFunc(GL_LEQUAL); 
	glDepthMask(GL_TRUE); 

	glCallList (surfellist);
	// end test ...
	*/
	

        glCallList (surface_vector_list);
        glDisable(GL_CLIP_PLANE0);
      }


    if (showclipsolution)
      {
	if (clipsolution == 1)
	  {
	    // Martin 
	    // orig:
	    glCallList (clipplanelist_scal);

	    // transparent experiments
	    // see http://wiki.delphigl.com/index.php/Blenden

	    /*
	    glColor4f (1, 1, 1, 0.5);
	    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);   
	    glEnable(GL_BLEND); 
	    glEnable(GL_COLOR);
	    glDepthFunc(GL_GREATER); 
	    glDepthMask(GL_FALSE); 

	    glCallList (clipplanelist_scal); 
	    glDepthFunc(GL_LEQUAL); 
	    glDepthMask(GL_TRUE); 

	    glCallList (clipplanelist_scal);
	    glDisable(GL_BLEND); 
	    */


	    /*
	      // latest transparent version ...
	    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);   
	    glEnable(GL_BLEND); 
	    glEnable(GL_DEPTH_TEST);

	    // CreateTexture (numtexturecols, lineartexture, 0.25, GL_MODULATE);
	    // glCallList (clipplanelist_scal); 

	    glEnable(GL_BLEND); 
	    // glDisable(GL_DEPTH_TEST);
	    
	    // CreateTexture (numtexturecols, lineartexture, 0.25, GL_MODULATE);
	    glCallList (clipplanelist_scal); 


	    // glDepthFunc(GL_LEQUAL); 
	    // glDepthMask(GL_TRUE); 
	    // glCallList (clipplanelist_scal);
	    glEnable(GL_DEPTH_TEST);
	    glDisable(GL_BLEND); 
	    */
	    // end test
	  } 
	if (clipsolution == 2)
	  {
	    // glDisable(GL_DEPTH_TEST);
	    glCallList (clipplanelist_vec);
	    // glEnable(GL_DEPTH_TEST);
	  }
      }



    if (draw_fieldlines)
      {
	SetClippingPlane();
        if (num_fieldlineslists <= 1)
          glCallList (fieldlineslist);
        else
          {  // animated
            int start = int (time / 10 * num_fieldlineslists);
            for (int ln = 0; ln < 10; ln++)
              {
                int nr = fieldlineslist + (start + ln) % num_fieldlineslists;
                glCallList (nr);
              }
          }
        glDisable(GL_CLIP_PLANE0);
      }

    if(drawpointcurves)
      {
	glCallList(pointcurvelist);
      }


    glMatrixMode (GL_TEXTURE);
    glLoadIdentity();
    glMatrixMode (GL_MODELVIEW);

    glDisable (GL_TEXTURE_1D);
    glDisable (GL_TEXTURE_2D);

    glDisable (GL_POLYGON_OFFSET_FILL);
    glDisable (GL_COLOR_MATERIAL);

    if (draw_isosurface)
      glCallList (isosurface_list);
    
    
    GLfloat matcol0[] = { 0, 0, 0, 1 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matcol0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, matcol0);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matcol0);
    
    glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth (1.0f);
    glColor3f (0.0f, 0.0f, 0.0f);
    glDisable (GL_LINE_SMOOTH);

    if (vispar.drawedges)
      {
        glCallList (element1dlist);
      }

    if (vispar.drawoutline && !numisolines)
      {
	SetClippingPlane ();
	glDepthMask(GL_FALSE); 
        glCallList (linelist);
	glDepthMask(GL_TRUE); 

        glDisable(GL_CLIP_PLANE0);
      }

    if (numisolines)
      {
        SetClippingPlane ();
        glCallList (isolinelist);

        glDisable(GL_CLIP_PLANE0);
        glCallList (clipplane_isolinelist);
      }


    // user visualization
    
    for (int i = 0; i < user_vis.Size(); i++)
      user_vis[i] -> Draw();

    glPopMatrix();
    
    glDisable(GL_CLIP_PLANE0);
    DrawColorBar (minval, maxval, logscale, lineartexture);
    
    if (vispar.drawcoordinatecross)
      DrawCoordinateCross ();
    DrawNetgenLogo ();
    
    glFinish();  

    
    // delete lock;
    // mem_lock.UnLock();
      }
    catch (bad_weak_ptr e)
      {
        cout << "don't have a mesh to visualize" << endl;
        VisualScene::DrawScene();      
      }
  }
  

  
  /*
  void VisualSceneSolution :: RealVec3d (const double * values, Vec3d & v, 
                                         bool iscomplex, bool imag)
  {
    if (!iscomplex)
      {
        v.X() = values[0];
        v.Y() = values[1];
        v.Z() = values[2];
      }
    else
      {
        if (!imag)
          {
            v.X() = values[0];
            v.Y() = values[2];
            v.Z() = values[4];
          }
        else
          {
            v.X() = values[1];
            v.Y() = values[3];
            v.Z() = values[5];
          }
      }
  }
  */
  Vec<3>  VisualSceneSolution :: RealVec3d (const double * values, 
					    bool iscomplex, bool imag)
  {
    Vec<3> v;
    if (!iscomplex)
      {
	for (int j = 0; j < 3; j++)
	  v(j) = values[j];
      }
    else
      {
        if (!imag)
          {
	    for (int j = 0; j < 3; j++)
	      v(j) = values[2*j];
          }
        else
          {
	    for (int j = 0; j < 3; j++)
	      v(j) = values[2*j+1];
          }
      }
    return v;
  }
  

  void VisualSceneSolution :: RealVec3d (const double * values, Vec3d & v, 
                                         bool iscomplex, double phaser, double phasei)
  {
    if (!iscomplex)
      {
        v.X() = values[0];
        v.Y() = values[1];
        v.Z() = values[2];
      }
    else
      {
        for (int i = 0; i < 3; i++)
          v.X(i+1) = phaser * values[2*i] + phasei * values[2*i+1];
      }
  }


  

  void VisualSceneSolution :: BuildScene (int zoomall)
  {
    try
      {
    shared_ptr<Mesh> mesh = GetMesh();

    if (!mesh)
      {
        VisualScene::BuildScene (zoomall);
        return;
      }

    /*
      if (!cone_list)
      {
      cone_list = glGenLists (1);
      glNewList (cone_list, GL_COMPILE);
      DrawCone (Point<3> (0,0,0), Point<3> (0,0,1), 0.4);
      glEndList();
      }
    */
    
    // vispar.colormeshsize = 1;
    
    // recalc clipping plane
    SetClippingPlane ();
    glDisable(GL_CLIP_PLANE0);
    
    
    SolData * sol = NULL;
    SolData * vsol = NULL;
  
    if (scalfunction != -1) 
      sol = soldata[scalfunction];
    if (vecfunction != -1)
      vsol = soldata[vecfunction];

    if (mesh->GetTimeStamp () > solutiontimestamp)
      {
        sol = NULL;
        vsol = NULL;
      }
 

    if (sol && sol->solclass) sol->solclass->SetMultiDimComponent (multidimcomponent);
    if (vsol && vsol->solclass) vsol->solclass->SetMultiDimComponent (multidimcomponent);

    if (!autoscale || (!sol && !vsol) )
      {
        minval = mminval;
        maxval = mmaxval;
      }
    else
      {
        if (mesh->GetTimeStamp () > surfeltimestamp ||
            vispar.clipping.timestamp > clipplanetimestamp ||
            solutiontimestamp > surfeltimestamp)
          {
            GetMinMax (scalfunction, scalcomp, minval, maxval);
          }
      }

    if (mesh->GetTimeStamp() > surfeltimestamp ||
        solutiontimestamp > surfeltimestamp || 
        zoomall)
      {
        if (mesh->GetTimeStamp() > surfeltimestamp || zoomall)
          {
            // mesh has changed
          
            Point3d pmin, pmax;
            static double oldrad = 0;
          
            mesh->GetBox (pmin, pmax, -1);
            center = Center (pmin, pmax);
            rad = 0.5 * Dist (pmin, pmax);
          
            glEnable (GL_NORMALIZE);
          
            if (rad > 1.5 * oldrad ||
                mesh->GetMajorTimeStamp() > surfeltimestamp ||
                zoomall)
              {
                CalcTransformationMatrices();
                oldrad = rad;
              }
          }
      
        DrawSurfaceElements();
      
        surfeltimestamp = max2 (solutiontimestamp, mesh->GetTimeStamp());
      }

    if (mesh->GetTimeStamp() > surfellinetimestamp ||
        subdivision_timestamp > surfellinetimestamp ||
        (deform && solutiontimestamp > surfellinetimestamp) || 
        zoomall)
      {      
        DrawSurfaceElementLines();      
        surfellinetimestamp = max2 (solutiontimestamp, mesh->GetTimeStamp());
      }


    if (vispar.drawedges)
      Draw1DElements();



    if (mesh->GetTimeStamp() > surface_vector_timestamp ||
        solutiontimestamp > surface_vector_timestamp ||
        zoomall)
      {
        if (surface_vector_list)
          glDeleteLists (surface_vector_list, 1);
      
        surface_vector_list = glGenLists (1);
        glNewList (surface_vector_list, GL_COMPILE);

        glEnable (GL_NORMALIZE);
        DrawSurfaceVectors();

        glEndList ();

        surface_vector_timestamp = 
          max2 (mesh->GetTimeStamp(), solutiontimestamp);
      }

    if (clipplanetimestamp < vispar.clipping.timestamp ||
        clipplanetimestamp < solutiontimestamp)
      {

        //      cout << "clipsolution = " << clipsolution << endl;
        if (vispar.clipping.enable && clipsolution == 2)      
          {
            // lock->UnLock();
            NgLock mlock (mesh->Mutex(), 0);
            mlock.UnLock(); 
            mesh->BuildElementSearchTree();
            mlock.Lock();

            // lock->Lock();
          }

      
        if (vispar.clipping.enable && clipsolution == 1 && sol)
	  DrawClipPlaneTrigs (); 

        if (clipplanelist_vec)
          glDeleteLists (clipplanelist_vec, 1);
      
        clipplanelist_vec = glGenLists (1);
        glNewList (clipplanelist_vec, GL_COMPILE);

        if (vispar.clipping.enable && clipsolution == 2 && vsol)
          {
            SetTextureMode (usetexture);

            if (autoscale)
              GetMinMax (vecfunction, 0, minval, maxval);

            Array<ClipPlanePoint> cpp;
            GetClippingPlaneGrid (cpp);

            for (int i = 0; i < cpp.Size(); i++)
              {
                const ClipPlanePoint & p = cpp[i];
                double values[6];
                Vec3d v;

                bool drawelem = 
                  GetValues (vsol, p.elnr, p.lami(0), p.lami(1), p.lami(2), values);
                // RealVec3d (values, v, vsol->iscomplex, imag_part);
		v = RealVec3d (values, vsol->iscomplex, imag_part);

                double val = v.Length();

                if (drawelem && val > 1e-10 * maxval)
                  {
                    v *= (rad / val / gridsize * 0.5);
                  
                    SetOpenGlColor  (val);
                    DrawCone (p.p, p.p+v, rad / gridsize * 0.2);
                  }
              }
          }

        glEndList ();
      }


    if (mesh->GetTimeStamp() > isosurface_timestamp ||
        solutiontimestamp > isosurface_timestamp ||
        zoomall)
      {
        if (isosurface_list)
          glDeleteLists (isosurface_list, 1);
      
        isosurface_list = glGenLists (1);
        glNewList (isosurface_list, GL_COMPILE);

        glEnable (GL_NORMALIZE);
        DrawIsoSurface(sol, vsol, scalcomp);

        glEndList ();

        isosurface_timestamp = 
          max2 (mesh->GetTimeStamp(), solutiontimestamp);
      }

    if(mesh->GetTimeStamp() > pointcurve_timestamp ||
       solutiontimestamp > pointcurve_timestamp)
      {
	if(pointcurvelist)
	  glDeleteLists(pointcurvelist,1);
	
		
	if(mesh->GetNumPointCurves() > 0)
	  {
	    pointcurvelist = glGenLists(1);
	    glNewList(pointcurvelist,GL_COMPILE);
	    //glColor3f (1.0f, 0.f, 0.f);
	    
	    for(int i=0; i<mesh->GetNumPointCurves(); i++)
	      {
		Box3d box;
		box.SetPoint(mesh->GetPointCurvePoint(i,0));
		for(int j=1; j<mesh->GetNumPointsOfPointCurve(i); j++)
		  box.AddPoint(mesh->GetPointCurvePoint(i,j));
		double diam = box.CalcDiam();
			     
		double thick = min2(0.1*diam, 0.001*rad);

		double red,green,blue;
		mesh->GetPointCurveColor(i,red,green,blue);
		glColor3f (red, green, blue);
		for(int j=0; j<mesh->GetNumPointsOfPointCurve(i)-1; j++)
		  {
		    DrawCylinder(mesh->GetPointCurvePoint(i,j),
				 mesh->GetPointCurvePoint(i,j+1),
				 thick);
		  }
	      }
	    glEndList();
	  }
	
      }


    if (
        numisolines && 
        (clipplanetimestamp < vispar.clipping.timestamp ||
         clipplanetimestamp < solutiontimestamp) 
        )
      {
        if (isolinelist) glDeleteLists (isolinelist, 1);
      
        isolinelist = glGenLists (1);
        glNewList (isolinelist, GL_COMPILE);

        Point<3> points[1100];
        double values[1100];
      
        int nse = mesh->GetNSE();

        CurvedElements & curv = mesh->GetCurvedElements();

        if (sol)
          {
            glBegin (GL_LINES);
          
            for (SurfaceElementIndex sei = 0; sei < nse; sei++)
              {
                const Element2d & el = (*mesh)[sei];

                bool curved = curv.IsHighOrder(); //  && curv.IsSurfaceElementCurved(sei);
              
                if (el.GetType() == TRIG || el.GetType() == TRIG6)
                  {
                    Point<3> lp1, lp2, lp3;
                    if (!curved)
                      {
                        GetPointDeformation (el[0]-1, lp1);
                        GetPointDeformation (el[1]-1, lp2);
                        GetPointDeformation (el[2]-1, lp3);
                      }
                  
                    int n = 1 << subdivisions;
                    int ii = 0;
                    int ix, iy;
                    for (iy = 0; iy <= n; iy++)
                      for (ix = 0; ix <= n-iy; ix++)
                        {
                          double x = double(ix) / n;
                          double y = double(iy) / n;
                        
                          // TODO: consider return value (bool: draw/don't draw element)
                          GetSurfValue (sol, sei, -1, x, y, scalcomp, values[ii]);
                          Point<2> xref(x,y);
                        
                          if (curved)
                            mesh->GetCurvedElements().
                              CalcSurfaceTransformation (xref, sei, points[ii]);
                          else
                            points[ii] = lp3 + x * (lp1-lp3) + y * (lp2-lp3);
                        
                          if (deform)
                            {
                              points[ii] += GetSurfDeformation (sei, -1, x, y);
                            }
                          ii++;
                        }
                  
                    ii = 0;
                    for (iy = 0; iy < n; iy++, ii++)
                      for (ix = 0; ix < n-iy; ix++, ii++)
                        {
                          int index[] = { ii, ii+1, ii+n-iy+1,
                                          ii+1, ii+n-iy+2, ii+n-iy+1 };
                        
                          DrawIsoLines (points[index[0]], points[index[1]], points[index[2]],
                                        values[index[0]], values[index[1]], values[index[2]]);

                          if (ix < n-iy-1) 
                            DrawIsoLines (points[index[3]], points[index[4]], points[index[5]],
                                          values[index[3]], values[index[4]], values[index[5]]);
                        }    
                  }
              
              
                if (el.GetType() == QUAD || el.GetType() == QUAD6 || el.GetType() == QUAD8 )
                  {
                    Point<3> lpi[4];
                    Vec<3> vx = 0.0, vy = 0.0, vtwist = 0.0, def;
                    if (!curved)
                      {
                        for (int j = 0; j < 4; j++)
                          GetPointDeformation (el[j]-1, lpi[j]);
                        vx = lpi[1]-lpi[0];
                        vy = lpi[3]-lpi[0];
                        vtwist = (lpi[0]-lpi[1]) + (lpi[2]-lpi[3]);
                      }

                    int n = 1 << subdivisions;
                    int ix, iy, ii = 0;
                    for (iy = 0; iy <= n; iy++)
                      for (ix = 0; ix <= n; ix++, ii++)
                        {
                          double x = double(ix) / n;
                          double y = double(iy) / n;
                        
                          // TODO: consider return value (bool: draw/don't draw element)
                          GetSurfValue (sol, sei, -1, x, y, scalcomp, values[ii]);
                          Point<2> xref(x,y);
                        
                          if (curved)
                            mesh->GetCurvedElements().
                              CalcSurfaceTransformation (xref, sei, points[ii]);
                          else
                            points[ii] = lpi[0] + x * vx + y * vy + x*y * vtwist;
                        
                          if (deform)
                            points[ii] += GetSurfDeformation (sei, -1, x, y);
                        }
                  
                    ii = 0;
                    for (iy = 0; iy < n; iy++, ii++)
                      for (ix = 0; ix < n; ix++, ii++)
                        {
                          DrawIsoLines (points[ii], points[ii+1], points[ii+n+1],
                                        values[ii], values[ii+1], values[ii+n+1]);
                          DrawIsoLines (points[ii+1], points[ii+n+2], points[ii+n+1],
                                        values[ii+1], values[ii+n+2], values[ii+n+1]);
                        }       
                  }
              }
            glEnd();
          }
        glEndList ();

        if (clipplane_isolinelist) glDeleteLists (clipplane_isolinelist, 1);
            
        if (vispar.clipping.enable && clipsolution == 1 && sol)
          {
            clipplane_isolinelist = glGenLists (1);
            glNewList (clipplane_isolinelist, GL_COMPILE);

            Array<ClipPlaneTrig> cpt;
            Array<ClipPlanePoint> pts;
            GetClippingPlaneTrigs (cpt, pts);  
            bool drawelem;
          
            glNormal3d (-clipplane[0], -clipplane[1], -clipplane[2]);
          
            if (numisolines)
              for (int i = 0; i < cpt.Size(); i++)
                {
                  const ClipPlaneTrig & trig = cpt[i];
                  double vali[3];
                  for (int j = 0; j < 3; j++)
                    {
                      Point<3> lami = pts[trig.points[j].pnr].lami;
                      drawelem = GetValue (sol, trig.elnr, lami(0), lami(1), lami(2),
                                           scalcomp, vali[j]);
                    }
                  if ( drawelem )
                    DrawIsoLines (pts[trig.points[0].pnr].p,
                                  pts[trig.points[1].pnr].p,
                                  pts[trig.points[2].pnr].p,
                                  // trig.points[1].p,
                                  // trig.points[2].p,
                                  vali[0], vali[1], vali[2]);
                }
            glEndList ();
          }
        glEnd();
      }
  
    clipplanetimestamp = max2 (vispar.clipping.timestamp, solutiontimestamp);
      }
    catch (bad_weak_ptr e)
      {
        cout << "vssolution::buildscene: don't have a mesh to visualize" << endl;
        VisualScene::BuildScene (zoomall);
      }
  }
  
  void  VisualSceneSolution :: Draw1DElements ()
  {
    shared_ptr<Mesh> mesh = GetMesh();

    if (element1dlist)
      glDeleteLists (element1dlist, 1);
 
    element1dlist = glGenLists (1);
    glNewList (element1dlist, GL_COMPILE);

    int npt = (1 << subdivisions) + 1;
    Array<double> pref(npt), values(npt);
    Array<Point<3> > points(npt);

    const SolData * sol = NULL;
    if (scalfunction != -1) sol = soldata[scalfunction];

    int ncomp = 0;
    if (sol) ncomp = sol->components;
    Array<double> mvalues(ncomp);


    for (int i = 0; i < npt; i++)
      pref[i] = double(i) / (npt-1);

    for (SegmentIndex i = 0; i < mesh -> GetNSeg(); i++)
      {
        // mesh->GetCurvedElements().
        // CalcMultiPointSegmentTransformation (&pref, i, &points, NULL);
        // const Segment & seg = mesh -> LineSegment(i);
        for (int j = 0; j < npt; j++)
          mesh->GetCurvedElements().
            CalcSegmentTransformation (pref[j], i, points[j]);
        if (sol)
          {
            for (int j = 0; j < npt; j++)
              {
                sol->solclass->GetSegmentValue (i, pref[j], &mvalues[0]);
                values[j] = ExtractValue (sol, scalcomp, &mvalues[0]);
                points[j](1) += scaledeform * values[j];
              }
          }

        glBegin (GL_LINE_STRIP);
        for (int i = 0; i < npt; i++)
          glVertex3dv (points[i]);
        glEnd();
      }

    glEndList ();
  }
  
  void  VisualSceneSolution :: DrawSurfaceElements ()
  {
    shared_ptr<Mesh> mesh = GetMesh();

    static int timer = NgProfiler::CreateTimer ("Solution::DrawSurfaceElements");
    NgProfiler::RegionTimer reg (timer);
  
    
#ifdef PARALLELGL

    if (id == 0 && ntasks > 1)
      {
	InitParallelGL();

	par_surfellists.SetSize (ntasks);

	MyMPI_SendCmd ("redraw");
	MyMPI_SendCmd ("solsurfellist");

	for ( int dest = 1; dest < ntasks; dest++ )
	  MyMPI_Recv (par_surfellists[dest], dest, MPI_TAG_VIS);

	if (surfellist)
	  glDeleteLists (surfellist, 1);

	surfellist = glGenLists (1);
	glNewList (surfellist, GL_COMPILE);
	
	for ( int dest = 1; dest < ntasks; dest++ )
	  glCallList (par_surfellists[dest]);
	
	glEndList();
	return;
      }
#endif



    if (surfellist)
      glDeleteLists (surfellist, 1);
    
    surfellist = glGenLists (1);
    glNewList (surfellist, GL_COMPILE);

    
    const SolData * sol = NULL;
    
    if (scalfunction != -1)
      sol = soldata[scalfunction];
    
    if (mesh->GetTimeStamp () > solutiontimestamp)
      sol = NULL;

    if (sol && sol->solclass) sol->solclass->SetMultiDimComponent (multidimcomponent);



    glLineWidth (1.0f);

    GLfloat col_grey[] = { 0.6f, 0.6f, 0.6f, 1.0f };
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, col_grey);
        

    int nse = mesh->GetNSE();

    SetTextureMode (usetexture);

    CurvedElements & curv = mesh->GetCurvedElements();

    int n = 1 << subdivisions;
    int npt = sqr(n+1);

    Array<Point<2> > pref (npt);
    Array<Point<3> > points (npt);
    Array<Mat<3,2> > dxdxis (npt);
    Array<Vec<3> > nvs(npt);
    Array<double> values(npt);

    Array<double> mvalues(npt);
    if (sol && sol->draw_surface) mvalues.SetSize (npt * sol->components);

    Array<complex<double> > valuesc(npt);

    for (SurfaceElementIndex sei = 0; sei < nse; sei++)
      {
        const Element2d & el = (*mesh)[sei];

        if (vispar.drawdomainsurf > 0)
          {
            if (mesh->GetDimension() == 3)
              {
                if (vispar.drawdomainsurf != mesh->GetFaceDescriptor(el.GetIndex()).DomainIn() &&
                    vispar.drawdomainsurf != mesh->GetFaceDescriptor(el.GetIndex()).DomainOut())
                  continue;
              }
            else
              {
                if (el.GetIndex() != vispar.drawdomainsurf) continue;
              }
          }



        if ( el.GetType() == QUAD || el.GetType() == QUAD6 )
          {
            bool curved = curv.IsSurfaceElementCurved (sei);


            for (int iy = 0, ii = 0; iy <= n; iy++)
              for (int ix = 0; ix <= n; ix++, ii++)
                pref[ii] = Point<2> (double(ix)/n, double(iy)/n);

            int npt = (n+1)*(n+1);
            if (curved)
              {
                for (int ii = 0; ii < npt; ii++)
                  {
                    Point<2> xref = pref[ii];
                    
                    mesh->GetCurvedElements().
                      CalcSurfaceTransformation (xref, sei, points[ii], dxdxis[ii]);
                    nvs[ii] = Cross (dxdxis[ii].Col(0), dxdxis[ii].Col(1));
                    nvs[ii].Normalize();
                  }
              }
            else
              {
		Point<3> lpi[4];
		Vec<3> vx, vy, vtwist;
		
		for (int k = 0; k < 4; k++)
		  GetPointDeformation (el[k]-1, lpi[k]);
		
		vx = lpi[1]-lpi[0];
		vy = lpi[3]-lpi[0];
		vtwist = (lpi[0]-lpi[1]) + (lpi[2]-lpi[3]);

                for (int ii = 0; ii < npt; ii++)
                  {
                    double x = pref[ii](0);
                    double y = pref[ii](1);
                    points[ii] = lpi[0] + x * vx + y * vy + x*y * vtwist;
                    for (int j = 0; j < 3; j++)
                      {
                        dxdxis[ii](j,0) = vx(j) + y*vtwist(j);
                        dxdxis[ii](j,1) = vy(j) + x*vtwist(j);
                      }
                  }

                Vec<3> nv = Cross (vx, vy);
                nv.Normalize();
                for (int ii = 0; ii < npt; ii++)
                  nvs[ii] = nv;
              }


            bool drawelem = false;
	    /*
            if (sol && sol->draw_surface) 
              {
                if (usetexture == 2)
                  for (int ii = 0; ii < npt; ii++)
                    drawelem = GetSurfValueComplex (sol, sei, -1, pref[ii](0), pref[ii](1), scalcomp, valuesc[ii]);
                else
                  for (int ii = 0; ii < npt; ii++)
                    drawelem = GetSurfValue (sol, sei, -1, pref[ii](0), pref[ii](1), scalcomp, values[ii]);
              }
	    */
            if (sol && sol->draw_surface) 
              {
		drawelem = GetMultiSurfValues (sol, sei, -1, npt, 
					       &pref[0](0), &pref[1](0)-&pref[0](0),
					       &points[0](0), &points[1](0)-&points[0](0),
					       &dxdxis[0](0), &dxdxis[1](0)-&dxdxis[0](0),
					       &mvalues[0], sol->components);
                if (usetexture == 2)
		  for (int ii = 0; ii < npt; ii++)
		    valuesc[ii] = ExtractValueComplex(sol, scalcomp, &mvalues[ii*sol->components]);
                else
		  for (int ii = 0; ii < npt; ii++)
		    values[ii] = ExtractValue(sol, scalcomp, &mvalues[ii*sol->components]);
              }

            
            if (deform)
              for (int ii = 0; ii < npt; ii++)
                points[ii] += GetSurfDeformation (sei, -1, pref[ii](0), pref[ii](1));


            int save_usetexture = usetexture;
            if (!drawelem)
              {
                usetexture = 0;
                SetTextureMode (0);
              }

            int ii = 0;

            glBegin (GL_QUADS);

            for (int iy = 0; iy < n; iy++, ii++)
              for (int ix = 0; ix < n; ix++, ii++)
                {
                  int index[] = { ii, ii+1, ii+n+2, ii+n+1 };
                  
                  for (int j = 0; j < 4; j++)
                    {
                      if (drawelem)
                        {
                          if (usetexture != 2)
                            SetOpenGlColor  (values[index[j]]);
                          else
                            glTexCoord2f ( valuesc[index[j]].real(),
                                           valuesc[index[j]].imag() );
                        }
                      else
                        glColor4fv (col_grey);
                      
                      glNormal3dv (nvs[index[j]]);
                      glVertex3dv (points[index[j]]);
                    }
                }
            glEnd();

            if (!drawelem && (usetexture != save_usetexture))
              {
                usetexture = save_usetexture;
                SetTextureMode (usetexture);
              }

          }
      }

    n = 1 << subdivisions;
    double invn = 1.0 / n;
    npt = (n+1)*(n+2)/2;

    for(SurfaceElementIndex sei = 0; sei < nse; sei++)
      {
        const Element2d & el = (*mesh)[sei];
	// if (el.GetIndex() <= 1) continue;

        if(vispar.drawdomainsurf > 0)
	  {
	    if (mesh->GetDimension() == 3)
	      {
		if (vispar.drawdomainsurf != mesh->GetFaceDescriptor(el.GetIndex()).DomainIn() &&
		    vispar.drawdomainsurf != mesh->GetFaceDescriptor(el.GetIndex()).DomainOut())
		  continue;
	      }
	    else
	      {
		if (el.GetIndex() != vispar.drawdomainsurf)
		  continue;
	      }
	  }
        
        if ( el.GetType() == TRIG || el.GetType() == TRIG6 )
          {
	    bool curved = curv.IsSurfaceElementCurved(sei);

            for (int iy = 0, ii = 0; iy <= n; iy++)
              for (int ix = 0; ix <= n-iy; ix++, ii++)
                pref[ii] = Point<2> (ix*invn, iy*invn);

            if (curved)
              {
                mesh->GetCurvedElements().
                  CalcMultiPointSurfaceTransformation (&pref, sei, &points, &dxdxis);

                for (int ii = 0; ii < npt; ii++)
                  nvs[ii] = Cross (dxdxis[ii].Col(0), dxdxis[ii].Col(1)).Normalize();
              }
            else
              {
		Point<3> p1 = mesh->Point (el[0]);
		Point<3> p2 = mesh->Point (el[1]);
		Point<3> p3 = mesh->Point (el[2]);

                Vec<3> vx = p1-p3;
                Vec<3> vy = p2-p3;
                for (int ii = 0; ii < npt; ii++)
                  {
                    points[ii] = p3 + pref[ii](0) * vx + pref[ii](1) * vy;
                    for (int j = 0; j < 3; j++)
                      {
                        dxdxis[ii](j,0) = vx(j);
                        dxdxis[ii](j,1) = vy(j);
                      }
                  }

                Vec<3> nv = Cross (vx, vy).Normalize();
                for (int ii = 0; ii < npt; ii++)
                  nvs[ii] = nv;
              }

            bool drawelem = false;
            if (sol && sol->draw_surface) 
              {
		drawelem = GetMultiSurfValues (sol, sei, -1, npt, 
					       &pref[0](0), &pref[1](0)-&pref[0](0),
					       &points[0](0), &points[1](0)-&points[0](0),
					       &dxdxis[0](0), &dxdxis[1](0)-&dxdxis[0](0),
					       &mvalues[0], sol->components);
                if (usetexture == 2)
		  for (int ii = 0; ii < npt; ii++)
		    valuesc[ii] = ExtractValueComplex(sol, scalcomp, &mvalues[ii*sol->components]);
                else
		  for (int ii = 0; ii < npt; ii++)
		    values[ii] = ExtractValue(sol, scalcomp, &mvalues[ii*sol->components]);
              }
            
            if (deform)
              for (int ii = 0; ii < npt; ii++)
                points[ii] += GetSurfDeformation (sei, -1, pref[ii](0), pref[ii](1));

            int save_usetexture = usetexture;
            if (!drawelem)
              {
                usetexture = 0;
                SetTextureMode (usetexture);
              }

            for (int iy = 0, ii = 0; iy < n; iy++)
              {
                glBegin (GL_TRIANGLE_STRIP);
                for (int ix = 0; ix <= n-iy; ix++, ii++)
                  for (int k = 0; k < 2; k++)
                    {
                      if (ix+iy+k > n) continue;
                      int hi = (k == 0) ? ii : ii+n-iy+1;
                      
                      if (drawelem)
                        {
                          if (usetexture != 2)
                            SetOpenGlColor (values[hi]); 
                          else
                            glTexCoord2f ( valuesc[hi].real(), valuesc[hi].imag() );
                        }
                      else
                        glColor4fv (col_grey);
                      
                      glNormal3dv (nvs[hi]);
                      glVertex3dv (points[hi]);
                    }
                glEnd();
              }


	    /*

 	    GLuint vboId[3];
 	    glGenBuffersARB (3, &vboId[0]);
// 	    cout << "vboId = " << vboId << endl;

 	    glBindBufferARB (GL_ARRAY_BUFFER_ARB, vboId[0]);
 	    glBufferDataARB (GL_ARRAY_BUFFER_ARB, points.Size()*sizeof(Point<3>), 
 			     &points[0][0], GL_STATIC_DRAW_ARB);


	    // not so fast as old-fashened style
	    glEnableClientState(GL_VERTEX_ARRAY);
	    // glVertexPointer(3, GL_DOUBLE, 0, &points[0][0]);
	    glVertexPointer(3, GL_DOUBLE, 0, 0);  //ARB

 	    glBindBufferARB (GL_ARRAY_BUFFER_ARB, vboId[1]);
 	    glBufferDataARB (GL_ARRAY_BUFFER_ARB, nvs.Size()*sizeof(Point<3>), 
 			     &nvs[0][0], GL_STATIC_DRAW_ARB);

	    glEnableClientState(GL_NORMAL_ARRAY);
	    // glNormalPointer(GL_DOUBLE, 0, &nvs[0][0]);
	    glNormalPointer(GL_DOUBLE, 0, 0);  // ARB

	    // if (drawelem && usetexture == 1)
	      {
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		// glTexCoordPointer(1, GL_DOUBLE, 0, &values[0]);

		glBindBufferARB (GL_ARRAY_BUFFER_ARB, vboId[2]);
		glBufferDataARB (GL_ARRAY_BUFFER_ARB, values.Size()*sizeof(double), 
				 &values[0], GL_STATIC_DRAW_ARB);
		glTexCoordPointer(1, GL_DOUBLE, 0, 0);
	      }

	    Array<int> gind;

            for (int iy = 0, ii = 0; iy < n; iy++,ii++)
              {
                for (int ix = 0; ix < n-iy; ix++, ii++)
                    {
                      int nv = (ix+iy+1 < n) ? 6 : 3;

		      int ind[] = { ii, ii+1, ii+n-iy+1,  
				    ii+n-iy+1, ii+1, ii+n-iy+2 };

// 		      if (ix == 0 && iy == 0)
// 		      for (int l = 0; l < 3; l++)
// 			{
// 			  if (drawelem)
// 			    {
// 			      if (usetexture != 2)
// 				// SetOpenGlColor (values[ind[l]]); 
// 				glTexCoord1f ( values[ind[l]] );
// 			      else
// 				glTexCoord2f ( valuesc[ind[l]].real(), valuesc[ind[l]].imag() );
// 			    }
// 			  else
// 			    glColor3fv (col_grey);
// 			}

		      for (int j = 0; j < nv; j++)
			gind.Append(ind[j]);
		      // glDrawElements(GL_TRIANGLES, nv, GL_UNSIGNED_INT, &ind[0]);
                    }
              }
	    glDrawElements(GL_TRIANGLES, gind.Size(), GL_UNSIGNED_INT, &gind[0]);

	    glDisableClientState(GL_VERTEX_ARRAY);
	    glDisableClientState(GL_NORMAL_ARRAY);
	    glDisableClientState(GL_TEXTURE_COORD_ARRAY);	    
	    
 	    glDeleteBuffersARB (3, &vboId[0]);
	    */
	    

            if (!drawelem && (usetexture != save_usetexture))
              {
                usetexture = save_usetexture;
                SetTextureMode (usetexture);
              }
	  }
      }
    glEndList ();


#ifdef PARALLELGL
    glFinish();
    if (id > 0)
      MyMPI_Send (surfellist, 0, MPI_TAG_VIS);
#endif
  }


  void  VisualSceneSolution :: DrawSurfaceElementLines ()
  {
    shared_ptr<Mesh> mesh = GetMesh();

#ifdef PARALLELGL
    if (id == 0 && ntasks > 1)
      {
	InitParallelGL();

	par_surfellists.SetSize (ntasks);

	MyMPI_SendCmd ("redraw");
	MyMPI_SendCmd ("solsurfellinelist");

	for ( int dest = 1; dest < ntasks; dest++ )
	  MyMPI_Recv (par_surfellists[dest], dest, MPI_TAG_VIS);

	if (linelist)
	  glDeleteLists (linelist, 1);

	linelist = glGenLists (1);
	glNewList (linelist, GL_COMPILE);
	
	for ( int dest = 1; dest < ntasks; dest++ )
	  glCallList (par_surfellists[dest]);
	
	glEndList();
	return;
      }
#endif

    if (linelist)
      glDeleteLists (linelist, 1);
    
    linelist = glGenLists (1);
    glNewList (linelist, GL_COMPILE);

    glLineWidth (1.0f);

    int nse = mesh->GetNSE();
    CurvedElements & curv = mesh->GetCurvedElements();

    int n = 1 << subdivisions;
    ArrayMem<Point<2>, 65> ptsloc(n+1);
    ArrayMem<Point<3>, 65> ptsglob(n+1);

    double trigpts[3][2]  = { { 0, 0 }, { 0, 1 }, { 1, 0} };
    double trigvecs[3][2] = { { 1, 0 }, { 0, -1 }, { -1, 1} };

    double quadpts[4][2]  = { { 0, 0 },  { 1, 1 }, { 0, 1}, { 1, 0 } };
    double quadvecs[4][2] = { { 1, 0 },  { -1, 0}, { 0, -1}, { 0, 1 } };

    for (SurfaceElementIndex sei = 0; sei < nse; sei++)
      {
        Element2d & el = (*mesh)[sei];

        int nv = (el.GetType() == TRIG || el.GetType() == TRIG6) ? 3 : 4;
        for (int k = 0; k < nv; k++)
          {
            Point<2> p0;
            Vec<2> vtau;
            if (nv == 3)
	      {
		p0 = Point<2>(trigpts[k][0], trigpts[k][1]);
		vtau = Vec<2>(trigvecs[k][0], trigvecs[k][1]);
	      }
            else
	      {
		p0 = Point<2>(quadpts[k][0], quadpts[k][1]);
		vtau = Vec<2>(quadvecs[k][0], quadvecs[k][1]);
	      }

            glBegin (GL_LINE_STRIP);

	    for (int ix = 0; ix <= n; ix++)
	      ptsloc[ix] = p0 + (double(ix) / n) * vtau;
	    
	    curv.CalcMultiPointSurfaceTransformation (&ptsloc, sei, &ptsglob, 0);
	    
	    for (int ix = 0; ix <= n; ix++)
	      {
		if (deform)
		  ptsglob[ix] += GetSurfDeformation (sei, k, ptsloc[ix](0), ptsloc[ix](1));
		glVertex3dv (ptsglob[ix]);
	      }
	    
            glEnd ();
          }
      }
    glEndList ();


#ifdef PARALLELGL
    glFinish();
    if (id > 0)
      MyMPI_Send (linelist, 0, MPI_TAG_VIS);
#endif
  }









  void VisualSceneSolution :: DrawIsoSurface(const SolData * sol, 
                                             const SolData * vsol,
                                             int comp)
  {
    shared_ptr<Mesh> mesh = GetMesh();

    if (!draw_isosurface) return;
    if (!sol) return;

   
    SetTextureMode (0);
    glColor3d (1.0, 0, 0);
    glEnable (GL_COLOR_MATERIAL);


    glBegin (GL_TRIANGLES);

    int ne = mesh->GetNE();

    const int edgei[6][2] =
      { { 0, 1 }, { 0, 2 }, { 0, 3 },
        { 1, 2 }, { 1, 3 }, { 2, 3 } };
    
    double edgelam[6];
    Point<3> edgep[6];
    Vec<3> normp[6];
    double nodevali[4];
    
    int cntce;
    int cpe1 = 0, cpe2 = 0, cpe3 = 0;
    
    int n = 1 << subdivisions;
    int n3 = (n+1)*(n+1)*(n+1);
    
    Array<Point<3> > grid(n3);
    Array<Point<3> > locgrid(n3);
    Array<Mat<3,3> > trans(n3);
    Array<double> val(n3);
    Array<Vec<3> > grads(n3);
    Array<int> compress(n3);
    
    MatrixFixWidth<3> pointmat(8);
    grads = Vec<3> (0.0);

    for (ElementIndex ei = 0; ei < ne; ei++)
      {
        // if(vispar.clipdomain > 0 && vispar.clipdomain != (*mesh)[ei].GetIndex()) continue;
        // if(vispar.donotclipdomain > 0 && vispar.donotclipdomain == (*mesh)[ei].GetIndex()) continue;

        ELEMENT_TYPE type = (*mesh)[ei].GetType();
        if (type == HEX || type == PRISM || type == TET || type == PYRAMID)
          {
            const Element & el = (*mesh)[ei];
            
            int ii = 0;
            int cnt_valid = 0;
            
            for (int ix = 0; ix <= n; ix++)
              for (int iy = 0; iy <= n; iy++)
                for (int iz = 0; iz <= n; iz++, ii++)
                  {
                    Point<3> ploc;
                    compress[ii] = ii;
                    
                    switch (type)
                      {
                      case PRISM:
                        if (ix+iy <= n)
                          {
                            ploc = Point<3> (double(ix) / n, double(iy) / n, double(iz) / n);
                            compress[ii] = cnt_valid;
                            cnt_valid++;
                          }
                        else
                          compress[ii] = -1;
                        break;
                      case TET:
                        if (ix+iy+iz <= n)
                          {
                            ploc = Point<3> (double(ix) / n, double(iy) / n, double(iz) / n);
                            compress[ii] = cnt_valid;
                            cnt_valid++;
                          }
                        else
                          compress[ii] = -1;
                        break;
                      case HEX:
                        ploc = Point<3> (double(ix) / n, double(iy) / n, double(iz) / n);
                        break;
                      case PYRAMID:
                        ploc = Point<3> (double(ix) / n * (1-double(iz)/n),
                                         double(iy) / n * (1-double(iz)/n),
                                         double(iz)/n);
                        break;
                      default:
                        cerr << "case not implementd 878234" << endl;
                        ploc = 0.0;
                      }
                    if (compress[ii] != -1)
                      locgrid[compress[ii]] = ploc;
                  }
            
            if (type != TET && type != PRISM) cnt_valid = n3;
            
            
            if (mesh->GetCurvedElements().IsHighOrder() || 1)
              {
                mesh->GetCurvedElements().
                  CalcMultiPointElementTransformation (&locgrid, ei, &grid, &trans);
              }
            else
              {
                Vector shape(el.GetNP());
                for (int k = 0; k < el.GetNP(); k++)
                  for (int j = 0; j < 3; j++)
                    pointmat(k,j) = (*mesh)[el[k]](j);
                
                for (int i = 0; i < cnt_valid; i++)
                  {
                    el.GetShapeNew (locgrid[i], shape);
                    Point<3> pglob;
                    for (int j = 0; j < 3; j++)
                      {
                        pglob(j) = 0;
                        for (int k = 0; k < el.GetNP(); k++)
                          pglob(j) += shape(k) * pointmat(k,j);
                      }
                    grid[i] = pglob;
                  }
              }

            bool has_pos = 0, has_neg = 0;
                
            for (int i = 0; i < cnt_valid; i++)
              {
                GetValue (sol, ei, &locgrid[i](0), &grid[i](0), &trans[i](0), comp, val[i]);
        
                val[i] -= minval;

                if (vsol)
                  GetValues (vsol, ei, &locgrid[i](0), &grid[i](0), &trans[i](0), &grads[i](0));
                grads[i] *= -1;


                if (val[i] > 0)
                  has_pos = 1;
                else
                  has_neg = 1;
              }

            if (!has_pos || !has_neg) continue;
            
            for (int ix = 0; ix < n; ix++)
              for (int iy = 0; iy < n; iy++)
                for (int iz = 0; iz < n; iz++)
                  {
                    int base = iz + (n+1)*iy + (n+1)*(n+1)*ix;
                    int pi[8] = 
                      { base, base+(n+1)*(n+1), base+(n+1)*(n+1)+(n+1), base+(n+1),
                        base+1, base+(n+1)*(n+1)+1, base+(n+1)*(n+1)+(n+1)+1, base+(n+1)+1 };
                    
                    for (int j = 0; j < 8; j++)
                      pi[j] = compress[pi[j]];
                    
                    int tets[6][4] = 
                      { { 1, 2, 4, 5 },
                        { 4, 5, 2, 8 },
                        { 2, 8, 5, 6 },
                        { 2, 3, 4, 8 },
                        { 2, 3, 8, 6 },
                        { 3, 8, 6, 7 } };
                    
                    for (int ii = 0; ii < 6; ii++)
                      {
                        int teti[4];
                        for (int k = 0; k < 4; k++)
                          teti[k] = pi[tets[ii][k]-1];
                        
                        bool is_valid = 1;
                        for (int j = 0; j < 4; j++)
                          if (teti[j] == -1) is_valid = 0;
                        
                        if (!is_valid) continue;
                        
                        for (int j = 0; j < 4; j++)
                          nodevali[j] = val[teti[j]];
                        
                        cntce = 0;
                        for (int j = 0; j < 6; j++)
                          {
                            int lpi1 = edgei[j][0];
                            int lpi2 = edgei[j][1];
                            if ( (nodevali[lpi1] > 0) !=
                                 (nodevali[lpi2] > 0) )
                              {
                                Point<3> p1 = grid[teti[lpi1]];
                                Point<3> p2 = grid[teti[lpi2]];

                                edgelam[j] = nodevali[lpi2] / (nodevali[lpi2] - nodevali[lpi1]);
                                edgep[j] = grid[teti[lpi1]] + (1-edgelam[j]) * (grid[teti[lpi2]]-grid[teti[lpi1]]);
                                normp[j] = grads[teti[lpi1]] + (1-edgelam[j]) * (grads[teti[lpi2]]-grads[teti[lpi1]]);
                                
                                cntce++;
                                cpe3 = cpe2;
                                cpe2 = cpe1;
                                cpe1 = j;
                                if (cntce >= 3)
                                  {
                                    if (!vsol)
                                      {
                                        Point<3> points[3];
                                        
                                        points[0] = edgep[cpe1];
                                        points[1] = edgep[cpe2];
                                        points[2] = edgep[cpe3];

                                        Vec<3> normal = Cross (points[2]-points[0], points[1]-points[0]);
                                        if ( ( (normal * (p2-p1)) > 0 ) == ( nodevali[lpi1] < 0) )
                                          normal *= -1;
                                        glNormal3dv (normal);

                                        glVertex3dv (points[0]);
                                        glVertex3dv (points[1]);
                                        glVertex3dv (points[2]);
                                      }
                                    else
                                      {
                                        glNormal3dv (normp[cpe1]);
                                        glVertex3dv (edgep[cpe1]);
                                        glNormal3dv (normp[cpe2]);
                                        glVertex3dv (edgep[cpe2]);
                                        glNormal3dv (normp[cpe3]);
                                        glVertex3dv (edgep[cpe3]);
                                      }
                                  }
                              }
                          }
                      }
                  }
          }
      }
    glEnd();
  }
    







  void  VisualSceneSolution :: DrawTrigSurfaceVectors(const Array< Point<3> > & lp, 
                                                      const Point<3> & pmin, const Point<3> & pmax,
                                                      const int sei, const SolData * vsol)
  {
    shared_ptr<Mesh> mesh = GetMesh();

    int dir,dir1,dir2;
    double s,t;

    Vec<3> n = Cross (lp[1]-lp[0], lp[2]-lp[0]);
    Vec<3> na (fabs (n(0)), fabs(n(1)), fabs(n(2)));
    if (na(0) > na(1) && na(0) > na(2))
      dir = 1;
    else if (na(1) > na(2))
      dir = 2;
    else 
      dir = 3;
    
    dir1 = (dir % 3) + 1;
    dir2 = (dir1 % 3) + 1;

    Point<2> p2d[3];

    int k;

    for (k = 0; k < 3; k++)
      {
        p2d[k] = Point<2> ((lp[k](dir1-1) - pmin(dir1-1)) / (2*rad),
                           (lp[k](dir2-1) - pmin(dir2-1)) / (2*rad));
      }

    
    double minx2d, maxx2d, miny2d, maxy2d;
    minx2d = maxx2d = p2d[0](0);
    miny2d = maxy2d = p2d[0](1);
    for (k = 1; k < 3; k++)
      {
        minx2d = min2 (minx2d, p2d[k](0));
        maxx2d = max2 (maxx2d, p2d[k](0));
        miny2d = min2 (miny2d, p2d[k](1));
        maxy2d = max2 (maxy2d, p2d[k](1));
      }
    
    double mat11 = p2d[1](0) - p2d[0](0);
    double mat21 = p2d[1](1) - p2d[0](1);
    double mat12 = p2d[2](0) - p2d[0](0);
    double mat22 = p2d[2](1) - p2d[0](1);
    
    double det = mat11*mat22-mat21*mat12;
    double inv11 = mat22/det;
    double inv21 = -mat21/det;
    double inv12 = -mat12/det;
    double inv22 = mat11/det;
          
    //    cout << "drawsurfacevectors. xoffset = " << xoffset << ", yoffset = ";
    //    cout << yoffset << endl;
    
    for (s = xoffset/gridsize; s <= 1+xoffset/gridsize; s += 1.0 / gridsize)
      if (s >= minx2d && s <= maxx2d)
        for (t = yoffset/gridsize; t <= 1+yoffset/gridsize; t += 1.0 / gridsize)
          if (t >= miny2d && t <= maxy2d)
            {
              double lam1 = inv11 * (s - p2d[0](0)) + inv12 * (t-p2d[0](1));
              double lam2 = inv21 * (s - p2d[0](0)) + inv22 * (t-p2d[0](1));
              
              if (lam1 >= 0 && lam2 >= 0 && lam1+lam2 <= 1)
                {
                  Point<3> cp;
                  for (k = 0; k < 3; k++)
                    cp(k) = lp[0](k) + 
                      lam1 * (lp[1](k)-lp[0](k)) + 
                      lam2 * (lp[2](k)-lp[0](k));

                  Point<2> xref(lam1, lam2);
                  if (mesh->GetCurvedElements().IsHighOrder())
                    mesh->GetCurvedElements().
                      CalcSurfaceTransformation (xref, sei, cp);
                  
                  Vec<3> v;
                  double values[6];
                  bool drawelem = 
                    GetSurfValues (vsol, sei, -1, lam1, lam2, values);
                  
                  if (!vsol->iscomplex)
                    for (k = 0; k < 3; k++)
                      v(k) = values[k];
                  else
                    {
                      if (!imag_part)
                        for (k = 0; k < 3; k++)
                          v(k) = values[2*k];
                      else
                        for (k = 0; k < 3; k++)
                          v(k) = values[2*k+1];
                    }
                  
                  if (mesh->GetDimension() == 2)
                    if ( (!vsol->iscomplex && vsol->components != 3) ||
                         (vsol->iscomplex && vsol->components != 6) )
                      v(2) = 0;
                  
                  double val = v.Length();

                  SetOpenGlColor  (val); // (val, minval, maxval, logscale);  // change JS

                  if (val > 1e-10 * maxval)
                    v *= (rad / val / gridsize * 0.5);
                  else 
                    drawelem = 0;

                  if ( drawelem ) 
                    DrawCone (cp, cp+4*v, 0.8*rad / gridsize);
                }
            }
    
  }



  void  VisualSceneSolution :: DrawSurfaceVectors ()
  {
    shared_ptr<Mesh> mesh = GetMesh();

    SurfaceElementIndex sei;

    const SolData * vsol = NULL;
    // bool drawelem;

    if (vecfunction != -1)
      vsol = soldata[vecfunction];

    if (mesh->GetTimeStamp () > solutiontimestamp)
      vsol = NULL;
     
    if (!vsol) return;


    Point<3> pmin = center - Vec3d (rad, rad, rad);
    Point<3> pmax = center - Vec3d (rad, rad, rad);


    // glColor3d (1.0, 1.0, 1.0);
    // glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    if (vsol->draw_surface && showsurfacesolution)
      {
        int nse = mesh->GetNSE();
        for (sei = 0; sei < nse; sei++)
          {
            const Element2d & el = (*mesh)[sei];
          
            if (el.GetType() == TRIG || el.GetType() == TRIG6)
              {
          
                Array< Point<3> > lp(3);

                lp[0] = mesh->Point(el[2]);
                lp[1] = mesh->Point(el[0]);
                lp[2] = mesh->Point(el[1]);

                DrawTrigSurfaceVectors(lp,pmin,pmax,sei,vsol);
                
                /*
                  Vec<3> n = Cross (lp[1]-lp[0], lp[2]-lp[0]);
                  Vec<3> na (fabs (n(0)), fabs(n(1)), fabs(n(2)));
                  if (na(0) > na(1) && na(0) > na(2))
                  dir = 1;
                  else if (na(1) > na(2))
                  dir = 2;
                  else 
                  dir = 3;

                  dir1 = (dir % 3) + 1;
                  dir2 = (dir1 % 3) + 1;
          
                  for (k = 0; k < 3; k++)
                  {
                  p2d[k] = Point<2> ((lp[k](dir1-1) - pmin(dir1-1)) / (2*rad),
                  (lp[k](dir2-1) - pmin(dir2-1)) / (2*rad));
                  }
          
                  double minx2d, maxx2d, miny2d, maxy2d;
                  minx2d = maxx2d = p2d[0](0);
                  miny2d = maxy2d = p2d[0](1);
                  for (k = 1; k < 3; k++)
                  {
                  minx2d = min2 (minx2d, p2d[k](0));
                  maxx2d = max2 (maxx2d, p2d[k](0));
                  miny2d = min2 (miny2d, p2d[k](1));
                  maxy2d = max2 (maxy2d, p2d[k](1));
                  }

                  double mat11 = p2d[1](0) - p2d[0](0);
                  double mat21 = p2d[1](1) - p2d[0](1);
                  double mat12 = p2d[2](0) - p2d[0](0);
                  double mat22 = p2d[2](1) - p2d[0](1);

                  double det = mat11*mat22-mat21*mat12;
                  double inv11 = mat22/det;
                  double inv21 = -mat21/det;
                  double inv12 = -mat12/det;
                  double inv22 = mat11/det;
          
                  //      cout << "drawsurfacevectors. xoffset = " << xoffset << ", yoffset = ";
                  //      cout << yoffset << endl;
          
                  for (s = xoffset/gridsize; s <= 1+xoffset/gridsize; s += 1.0 / gridsize)
                  if (s >= minx2d && s <= maxx2d)
                  for (t = yoffset/gridsize; t <= 1+yoffset/gridsize; t += 1.0 / gridsize)
                  if (t >= miny2d && t <= maxy2d)
                  {
                  double lam1 = inv11 * (s - p2d[0](0)) + inv12 * (t-p2d[0](1));
                  double lam2 = inv21 * (s - p2d[0](0)) + inv22 * (t-p2d[0](1));
                    
                  if (lam1 >= 0 && lam2 >= 0 && lam1+lam2 <= 1)
                  {
                  Point<3> cp;
                  for (k = 0; k < 3; k++)
                  cp(k) = lp[0](k) + 
                  lam1 * (lp[1](k)-lp[0](k)) + 
                  lam2 * (lp[2](k)-lp[0](k));

                  Vec<3> v;
                  double values[6];
                  drawelem = GetSurfValues (vsol, sei, lam1, lam2, values);

                  if (!vsol->iscomplex)
                  for (k = 0; k < 3; k++)
                  v(k) = values[k];
                  else
                  {
                  if (!imag_part)
                  for (k = 0; k < 3; k++)
                  v(k) = values[2*k];
                  else
                  for (k = 0; k < 3; k++)
                  v(k) = values[2*k+1];
                  }
                        
                  if (mesh->GetDimension() == 2)
                  if ( (!vsol->iscomplex && vsol->components != 3) ||
                  (vsol->iscomplex && vsol->components != 6) )
                  v(2) = 0;
                        
                  double val = v.Length();
                  SetOpenGlColor  (val, minval, maxval, logscale);

                  if (val > 1e-10 * maxval)
                  v *= (rad / val / gridsize * 0.5);
                  else drawelem = 0;
                  // "drawelem": added 07.04.2004 (FB)
                  if ( drawelem ) DrawCone (cp, cp+4*v, 0.8*rad / gridsize);


                  }
                  }
                */
              }
            else if (el.GetType() == QUAD)
              {
                /*
		  Array < Point<3> > lp(3);

		  lp[0] = mesh->Point(el[0]);
		  lp[1] = mesh->Point(el[1]);
		  lp[2] = mesh->Point(el[2]);

		  DrawTrigSurfaceVectors(lp,pmin,pmax,sei,vsol);

		  lp[0] = mesh->Point(el[0]);
		  lp[1] = mesh->Point(el[2]);
		  lp[2] = mesh->Point(el[3]);

		  DrawTrigSurfaceVectors(lp,pmin,pmax,sei,vsol);
                */
                
                Point<3> lp[4];
                Point<2> p2d[4];
                
                for (int k = 0; k < 4; k++)
                  lp[k] = mesh->Point (el[k]);
                
                
                Vec<3> n = Cross (lp[1]-lp[0], lp[2]-lp[0]);
                Vec<3> na (fabs (n(0)), fabs(n(1)), fabs(n(2)));
                int dir, dir1, dir2;
                if (na(0) > na(1) && na(0) > na(2))
                  dir = 1;
                else if (na(1) > na(2))
                  dir = 2;
                else 
                  dir = 3;
                
                dir1 = (dir % 3) + 1;
                dir2 = (dir1 % 3) + 1;
                
                for (int k = 0; k < 4; k++)
                  {
                    p2d[k] = Point<2> ((lp[k](dir1-1) - pmin(dir1-1)) / (2*rad),
                                       (lp[k](dir2-1) - pmin(dir2-1)) / (2*rad));
                  }
                
                double minx2d, maxx2d, miny2d, maxy2d;
                minx2d = maxx2d = p2d[0](0);
                miny2d = maxy2d = p2d[0](1);
                for (int k = 1; k < 4; k++)
                  {
                    minx2d = min2 (minx2d, p2d[k](0));
                    maxx2d = max2 (maxx2d, p2d[k](0));
                    miny2d = min2 (miny2d, p2d[k](1));
                    maxy2d = max2 (maxy2d, p2d[k](1));
                  }
                
                for (double s = xoffset/gridsize; s <= 1+xoffset/gridsize; s += 1.0 / gridsize)
                  if (s >= minx2d && s <= maxx2d)
                    for (double t = yoffset/gridsize; t <= 1+yoffset/gridsize; t += 1.0 / gridsize)
                      if (t >= miny2d && t <= maxy2d)
                        {
                          double lami[3];
                          Point3d p3d(2*rad*s+pmin(0), 2*rad*t+pmin(1),0);
                          
                          if (mesh->PointContainedIn2DElement (p3d, lami, sei+1))
                            {
                              Point<3> cp = p3d;
                              double lam1 = lami[0];
                              double lam2 = lami[1];
                              
                              //for (k = 0; k < 3; k++)
                              //cp(k) = lp[0](k) + 
                              //lam1 * (lp[1](k)-lp[0](k)) + 
                              //lam2 * (lp[2](k)-lp[0](k));
                              
                              
                              Vec<3> v;
                              double values[6];
                              bool drawelem = GetSurfValues (vsol, sei, -1, lam1, lam2, values);
                              (*testout) << "sei " << sei << " lam1 " << lam1 << " lam2 " << lam2 << " drawelem " << drawelem << endl;
                              
                              if (!vsol->iscomplex)
                                for (int k = 0; k < 3; k++)
                                  v(k) = values[k];
                              else
                                {
                                  if (!imag_part)
                                    for (int k = 0; k < 3; k++)
                                      v(k) = values[2*k];
                                  else
                                    for (int k = 0; k < 3; k++)
                                      v(k) = values[2*k+1];
                                }
                              
                              if (mesh->GetDimension() == 2)
                                if ( (!vsol->iscomplex && vsol->components != 3) ||
                                     (vsol->iscomplex && vsol->components != 6) )
                                  v(2) = 0;
                              
                              double val = v.Length();
                              SetOpenGlColor  (val); // , minval, maxval, logscale); july 09
                              
                              (*testout) << "v " << v << endl;
                              
                              if (val > 1e-10 * maxval)
                                v *= (rad / val / gridsize * 0.5);
                              
                              (*testout) << "v " << v << endl;
                              
                              if ( drawelem )
                                {
                                  DrawCone (cp, cp+4*v, 0.8*rad / gridsize);
                                  (*testout) << "cp " << cp << " rad " << rad << " gridsize " << gridsize << endl;
                                }
                              
                            }
                        }
              }
          }
      }
  }
  
  
  
  
  void VisualSceneSolution :: 
  DrawIsoLines (const Point<3> & p1, 
                const Point<3> & p2, 
                const Point<3> & p3,
                double val1, double val2, double val3)
  {
    DrawIsoLines2 (p1, p2, p1, p3, val1, val2, val1, val3); // , minval, maxval, n);
    DrawIsoLines2 (p2, p1, p2, p3, val2, val1, val2, val3); // , minval, maxval, n);
    DrawIsoLines2 (p3, p1, p3, p2, val3, val1, val3, val2); // , minval, maxval, n);
  }


  void VisualSceneSolution :: 
  DrawIsoLines2 (const Point<3> & hp1, 
                 const Point<3> & hp2, 
                 const Point<3> & hp3,
                 const Point<3> & hp4,
                 double val1, double val2, double val3, double val4)
  {
    int n = numisolines;
    Point<3> p1, p2, p3, p4;
    if (val1 < val2)
      {
        p1 = hp1; p2 = hp2;
      }
    else
      {
        p1 = hp2; p2 = hp1;
        swap (val1, val2);
      }

    if (val3 < val4)
      {
        p3 = hp3; p4 = hp4;
      }
    else
      {
        p3 = hp4; p4 = hp3;
        swap (val3, val4);
      }

    val2 += 1e-10;
    val4 += 1e-10;

    double fac = (maxval-minval) / n;
    double idelta1 = 1.0 / (val2 - val1);
    double idelta2 = 1.0 / (val4 - val3);

    int mini = int ((max2 (val1, val3) - minval) / fac);
    int maxi = int ((min2 (val2, val4) - minval) / fac);
    if (mini < 0) mini = 0;
    if (maxi > n-1) maxi = n-1;

    for (int i = mini; i <= maxi; i++)
      {
        double val = minval + i * fac;
        double lam1 = (val - val1) * idelta1;
        double lam2 = (val - val3) * idelta2;
        if (lam1 >= 0 && lam1 <= 1 && lam2 >= 0 && lam2 <= 1)
          {
            Point<3> lp1 = p1 + lam1 * (p2-p1);
            Point<3> lp2 = p3 + lam2 * (p4-p3);
            glVertex3dv (lp1 );
            glVertex3dv (lp2 );
            // glVertex3dv (lp2 );  // better ?
            // glVertex3dv (lp1 );  
          }
      }
  }



  void VisualSceneSolution :: 
  GetMinMax (int funcnr, int comp, double & minv, double & maxv) const
  {
    shared_ptr<Mesh> mesh = GetMesh();

    static int timer1 = NgProfiler::CreateTimer ("getminmax, vol");
    static int timer2 = NgProfiler::CreateTimer ("getminmax, surf");

#ifdef PARALLEL
    if (id == 0)
      {
	MyMPI_SendCmd ("redraw");
	MyMPI_SendCmd ("getminmax");
      }
    MyMPI_Bcast (funcnr);
    MyMPI_Bcast (comp);
#endif

    // double val;
    // bool considerElem;

    bool hasit = false;
#ifdef max
      #undef max
#endif
    minv = numeric_limits<double>::max();
    maxv = -numeric_limits<double>::max();

    if ((ntasks == 1) || (id > 0))
      if (funcnr != -1)
        {
          const SolData * sol = soldata[funcnr];

          if (sol->draw_volume)
            {
              NgProfiler::RegionTimer reg1 (timer1);
              
              int ne = mesh->GetNE();
              double hminv = numeric_limits<double>::max();
              double hmaxv = -numeric_limits<double>::max();
              bool hhasit = false;
#if defined _OPENMP && _OPENMP >= 201107
#pragma omp parallel for reduction (max : hmaxv) reduction (min : hminv) reduction (|| : hhasit)
#endif
              for (int i = 0; i < ne; i++) 
                {
                  double val;
                  bool considerElem = GetValue (sol, i, 0.333, 0.333, 0.333, comp, val);
                  if (considerElem)
                    {
                      if (val > hmaxv) hmaxv = val;
                      if (val < hminv) hminv = val;
                      hhasit = true;
                    }
                }

              minv = min(minv, hminv);
              maxv = max(maxv, hmaxv);
              hasit |= hhasit;
            }
          
          if (sol->draw_surface)
            {
              NgProfiler::RegionTimer reg2 (timer2);
              
              int nse = mesh->GetNSE();
              for (int i = 0; i < nse; i++)
                {
                  ELEMENT_TYPE type = mesh->SurfaceElement(i+1).GetType();
                  double val;
                  bool considerElem = (type == QUAD) 
                    ? GetSurfValue (sol, i, -1, 0.5, 0.5, comp, val)
                    : GetSurfValue (sol, i, -1, 0.3333333, 0.3333333, comp, val);
                  
                  if (considerElem)
                    {
                      if (val > maxv) maxv = val;
                      if (val < minv) minv = val;
                      hasit = true;
                    }
                }
            }
        }
    if (minv == maxv) maxv = minv+1e-6;
    if (!hasit) { minv = 0; maxv = 1; }
    
#ifdef PARALLEL
    if ((ntasks > 1) && (id == 0))
      {
	minv = 1e99;
	maxv = -1e99;
      }
    double hmin, hmax;
    MPI_Reduce (&minv, &hmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce (&maxv, &hmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    minv = hmin;
    maxv = hmax;
#endif
  }





  bool VisualSceneSolution :: 
  GetValues (const SolData * data, ElementIndex elnr, 
             double lam1, double lam2, double lam3,
             double * values) const
  {
    bool ok = false;
    switch (data->soltype)
      {
      case SOL_VIRTUALFUNCTION:
        {
          ok = data->solclass->GetValue (elnr, lam1, lam2, lam3, values);
          break;
        }
      default:
        {
          for (int i = 0; i < data->components; i++)
            ok = GetValue (data, elnr, lam1, lam2, lam3, i+1, values[i]);
        }
      }
    return ok;
  }

  bool VisualSceneSolution :: 
  GetValues (const SolData * data, ElementIndex elnr, 
             const double xref[], const double x[], const double dxdxref[], 
             double * values) const
  {
    bool ok = false;
    switch (data->soltype)
      {
      case SOL_VIRTUALFUNCTION:
        {
          ok = data->solclass->GetValue (elnr, xref, x, dxdxref, values);
          break;
        }
      default:
        {
          for (int i = 0; i < data->components; i++)
            ok = GetValue (data, elnr, xref[0], xref[1], xref[2], i+1, values[i]);
        }
      }
    return ok;
  }


  bool VisualSceneSolution :: 
  GetValue (const SolData * data, ElementIndex elnr, 
            const double xref[], const double x[], const double dxdxref[], 
            int comp, double & val) const
  {
    shared_ptr<Mesh> mesh = GetMesh();

    double lam1 = xref[0];
    double lam2 = xref[1];
    double lam3 = xref[2];
        
    val = 0;
    bool ok = 0;


    if (comp == 0)
      {
        ArrayMem<double,20> values(data->components);
        ok = GetValues (data, elnr, xref, x, dxdxref, &values[0]);

	val = ExtractValue (data, 0, &values[0]);
	return ok;
      }


    switch (data->soltype)
      {
      case SOL_VIRTUALFUNCTION:
        {
          double values[20];
          ok = data->solclass->GetValue (elnr, xref, x, dxdxref, values);

          val = values[comp-1];
          return ok;
        }
      case SOL_NODAL:
        {
          const Element & el = (*mesh)[elnr];

          double lami[8] = { 0.0 };
          int np = 0;
        
          switch (el.GetType())
            {
            case TET:
            case TET10:
              {
                lami[1] = lam1;
                lami[2] = lam2;
                lami[3] = lam3;
                lami[0] = 1-lam1-lam2-lam3;
                np = 4;
                break;
              }
            case PRISM:
            case PRISM12:
              {
                lami[0] = (1-lam3) * (1-lam1-lam2);
                lami[1] = (1-lam3) * lam1;
                lami[2] = (1-lam3) * lam2;
                lami[3] = (lam3) * (1-lam1-lam2);
                lami[4] = (lam3) * lam1;
                lami[5] = (lam3) * lam2;
                np = 6;
                break;
              }     
            default:
              cerr << "case not implementd 23523" << endl;
            }

          for (int i = 0; i < np; i++)
            val += lami[i] * data->data[(el[i]-1) * data->dist + comp-1];

          return 1;
        }

      case SOL_ELEMENT:
        {
          val = data->data[elnr * data->dist + comp-1];
          return 1;
        }

      case SOL_SURFACE_ELEMENT:
        return 0;

      case SOL_NONCONTINUOUS:
        {
          const Element & el = (*mesh)[elnr];

          double lami[8] = { 0.0 };
          int np = 0;

          switch (el.GetType())
            {
            case TET:
            case TET10:
              {
                lami[1] = lam1;
                lami[2] = lam2;
                lami[3] = lam3;
                lami[0] = 1-lam1-lam2-lam3;
                np = 4;
                break;
              }
            case PRISM:
            case PRISM12:
              {
                lami[0] = (1-lam3) * (1-lam1-lam2);
                lami[1] = (1-lam3) * lam1;
                lami[2] = (1-lam3) * lam2;
                lami[3] = (lam3) * (1-lam1-lam2);
                lami[4] = (lam3) * lam1;
                lami[5] = (lam3) * lam2;
                np = 6;
                break;
              }
            case PYRAMID:
              {
                if (lam3 > 1-1e-5)
                  {
                    lami[0] = lami[1] = lami[2] = lami[3] = 0;
                    lami[4] = 1;
                  }
                else
                  {
                    double x0 = lam1 / (1-lam3);
                    double y0 = lam2 / (1-lam3);
                    lami[0] = (1-x0) * (1-y0) * (1-lam3);
                    lami[1] = (  x0) * (1-y0) * (1-lam3);
                    lami[2] = (  x0) * (  y0) * (1-lam3);
                    lami[3] = (1-x0) * (  y0) * (1-lam3);
                    lami[4] = lam3;
                    np = 5;
                  }
                break;
              }
            default:
              np = 0;
            }

          int base;
          if (data->order == 1)
            base = 6 * elnr;
          else
            base = 10 * elnr;


          for (int i = 0; i < np; i++)
            val += lami[i] * data->data[(base+i) * data->dist + comp-1];

          return 1;
        }

      case SOL_MARKED_ELEMENTS:
        {
          val = (*mesh)[elnr].TestRefinementFlag();
          return 1;
        }
      
      case SOL_ELEMENT_ORDER:
        {
          val = (*mesh)[elnr].GetOrder();
          return 1;
        }

      default:
        cerr << "case not handled 7234" << endl;
      }
    return 0;
  }



  bool VisualSceneSolution :: 
  GetValue (const SolData * data, ElementIndex elnr, 
            double lam1, double lam2, double lam3,
            int comp, double & val) const
  {
    shared_ptr<Mesh> mesh = GetMesh();

    val = 0;
    bool ok = 0;

    if (comp == 0)
      {
        ArrayMem<double,20> values(data->components);
        ok = GetValues (data, elnr, lam1, lam2, lam3, &values[0]);
	val = ExtractValue (data, 0, &values[0]);
	return ok;
      }


    switch (data->soltype)
      {
      case SOL_VIRTUALFUNCTION:
        {
	  val = 0.0;
          double values[20];
	  ok = data->solclass->GetValue (elnr, lam1, lam2, lam3, values);

          val = values[comp-1];
          return ok;
        }
      case SOL_NODAL:
        {
          const Element & el = (*mesh)[elnr];

          double lami[8] = { 0.0 };
          int np = 0;
        
          switch (el.GetType())
            {
            case TET:
            case TET10:
              {
                lami[1] = lam1;
                lami[2] = lam2;
                lami[3] = lam3;
                lami[0] = 1-lam1-lam2-lam3;
                np = 4;
                break;
              }
            case PRISM:
            case PRISM12:
              {
                lami[0] = (1-lam3) * (1-lam1-lam2);
                lami[1] = (1-lam3) * lam1;
                lami[2] = (1-lam3) * lam2;
                lami[3] = (lam3) * (1-lam1-lam2);
                lami[4] = (lam3) * lam1;
                lami[5] = (lam3) * lam2;
                np = 6;
                break;
              }     
            default:
              cerr << "case not implemented 234324" << endl;
            }

          for (int i = 0; i < np; i++)
            val += lami[i] * data->data[(el[i]-1) * data->dist + comp-1];

          return 1;
        }

      case SOL_ELEMENT:
        {
          val = data->data[elnr * data->dist + comp-1];
          return 1;
        }

      case SOL_SURFACE_ELEMENT:
        return 0;

      case SOL_NONCONTINUOUS:
        {
          const Element & el = (*mesh)[elnr];

          double lami[8] = { 0.0 };
          int np = 0;

          switch (el.GetType())
            {
            case TET:
            case TET10:
              {
                lami[1] = lam1;
                lami[2] = lam2;
                lami[3] = lam3;
                lami[0] = 1-lam1-lam2-lam3;
                np = 4;
                break;
              }
            case PRISM:
            case PRISM12:
              {
                lami[0] = (1-lam3) * (1-lam1-lam2);
                lami[1] = (1-lam3) * lam1;
                lami[2] = (1-lam3) * lam2;
                lami[3] = (lam3) * (1-lam1-lam2);
                lami[4] = (lam3) * lam1;
                lami[5] = (lam3) * lam2;
                np = 6;
                break;
              }
            case PYRAMID:
              {
                if (lam3 > 1-1e-5)
                  {
                    lami[0] = lami[1] = lami[2] = lami[3] = 0;
                    lami[4] = 1;
                  }
                else
                  {
                    double x0 = lam1 / (1-lam3);
                    double y0 = lam2 / (1-lam3);
                    lami[0] = (1-x0) * (1-y0) * (1-lam3);
                    lami[1] = (  x0) * (1-y0) * (1-lam3);
                    lami[2] = (  x0) * (  y0) * (1-lam3);
                    lami[3] = (1-x0) * (  y0) * (1-lam3);
                    lami[4] = lam3;
                    np = 5;
                  }
                break;
              }
            default:
              np = 0;
            }

          int base;
          if (data->order == 1)
            base = 6 * elnr;
          else
            base = 10 * elnr;


          for (int i = 0; i < np; i++)
            val += lami[i] * data->data[(base+i) * data->dist + comp-1];

          return 1;
        }

      case SOL_MARKED_ELEMENTS:
        {
          val = (*mesh)[elnr].TestRefinementFlag();
          return 1;
        }
      
      case SOL_ELEMENT_ORDER:
        {
          val = (*mesh)[elnr].GetOrder();
          return 1;
        }
      default:
        cerr << "case not implemented 234234" << endl;
      }
    return 0;
  }







  bool VisualSceneSolution :: 
  GetValueComplex (const SolData * data, ElementIndex elnr, 
                   double lam1, double lam2, double lam3,
                   int comp, complex<double> & val) const
  {
    shared_ptr<Mesh> mesh = GetMesh();

    val = 0.0;
    bool ok = 0;

           
    switch (data->soltype)
      {
      case SOL_VIRTUALFUNCTION:
        {
          double values[20];
          ok = data->solclass->GetValue (elnr, lam1, lam2, lam3, values);
          val = complex<double> (values[comp-1], values[comp]);
          return ok;
        }
      default:
        cerr << "case not handled 234234" << endl;
      } 
    return 0;
  }
  

  bool VisualSceneSolution :: 
  GetMultiValues (const SolData * data, ElementIndex elnr, int facetnr, int npt,
		  const double * xref, int sxref,
		  const double * x, int sx,
		  const double * dxdxref, int sdxdxref,
		  double * val, int sval) const
  {
    bool drawelem = false;
    if (data->soltype == SOL_VIRTUALFUNCTION)
      drawelem = data->solclass->GetMultiValue(elnr, facetnr, npt, xref, sxref, x, sx, dxdxref, sdxdxref, val, sval);
    else
      for (int i = 0; i < npt; i++)
        drawelem = GetValues (data, elnr, xref+i*sxref, x+i*sx, dxdxref+i*sdxdxref, val+i*sval);
    return drawelem;
  }






  bool VisualSceneSolution :: 
  GetSurfValues (const SolData * data, SurfaceElementIndex selnr, int facetnr, 
                 double lam1, double lam2, 
                 double * values) const
  {
    bool ok = false;
    switch (data->soltype)
      {
      case SOL_VIRTUALFUNCTION:
        {
          ok = data->solclass->GetSurfValue (selnr, facetnr, lam1, lam2, values);
          // ok = 1;
          // values[0] = 1.0;
          break;
        }
      default:
        {
          for (int i = 0; i < data->components; i++)
            ok = GetSurfValue (data, selnr, facetnr, lam1, lam2, i+1, values[i]);
        }
      }
    return ok;
  }


  bool VisualSceneSolution :: 
  GetSurfValues (const SolData * data, SurfaceElementIndex selnr, int facetnr, 
                 const double xref[], const double x[], const double dxdxref[], 
                 double * values) const
  {
    bool ok = false;
    switch (data->soltype)
      {
      case SOL_VIRTUALFUNCTION:
        {
          ok = data->solclass->GetSurfValue (selnr, facetnr, xref, x, dxdxref, values);
          break;
        }
      default:
        {
          for (int i = 0; i < data->components; i++)
            ok = GetSurfValue (data, selnr, facetnr, xref[0], xref[1], i+1, values[i]);
        }
      }
    return ok;
  }

  bool VisualSceneSolution :: 
  GetMultiSurfValues (const SolData * data, SurfaceElementIndex elnr, int facetnr, int npt,
                      const double * xref, int sxref,
                      const double * x, int sx,
                      const double * dxdxref, int sdxdxref,
                      double * val, int sval) const
  {
    bool drawelem = false;
    if (data->soltype == SOL_VIRTUALFUNCTION)
      drawelem = data->solclass->GetMultiSurfValue(elnr, facetnr, npt, xref, sxref, x, sx, dxdxref, sdxdxref, val, sval);
    else
      for (int i = 0; i < npt; i++)
        drawelem = GetSurfValues (data, elnr, facetnr, xref+i*sxref, x+i*sx, dxdxref+i*sdxdxref, val+i*sval);
    return drawelem;
  }
  
  double VisualSceneSolution ::  ExtractValue (const SolData * data, int comp, double * values) const
  {
    double val = 0;
    if (comp == 0)
      {
        switch (evalfunc)
          {
          case FUNC_ABS:
            {
              for (int ci = 0; ci < data->components; ci++)
                val += sqr (values[ci]);
              val = sqrt (val);
              break;
            }
          case FUNC_ABS_TENSOR:
            {
              int d = 0;
              switch (data->components)
                {
                case 1: d = 1; break;
                case 3: d = 2; break;
                case 6: d = 3; break;
                }
              for (int ci = 0; ci < d; ci++)
                val += sqr (values[ci]);
              for (int ci = d; ci < data->components; ci++)
                val += 2*sqr (values[ci]);
              val = sqrt (val);
              break;
            }

          case FUNC_MISES:
            {
              int d = 0;
              switch(data->components)
                {
                case 1: d = 1; break;
                case 3: d = 2; break;
                case 6: d = 3; break;
                }
              int ci;
              double trace = 0.;
              for (ci = 0; ci < d; ci++)
                trace += 1./3.*(values[ci]);
              for (ci = 0; ci < d; ci++)
                val += sqr (values[ci]-trace);
              for (ci = d; ci < data->components; ci++)
                val += 2.*sqr (values[ci]);
              val = sqrt (val);
              break;
            }
          case FUNC_MAIN:
            {
              int d = 0;
              switch(data->components)
                {
                case 1: d = 1; break;
                case 3: d = 2; break;
                case 6: d = 3; break;
                }
              Mat<3,3> m ;
              Vec<3> ev;
              int ci;
              for (ci = 0; ci < d; ci++)
                m(ci,ci) = (values[ci]);
              m(0,1) = m(1,0) = values[3];
              m(0,2) = m(2,0) = values[4];
              m(1,2) = m(2,1) = values[5];

              EigenValues (m, ev);
              double help;
              for (int i=0; i<d; i++)
                {
                  for (int j=d-1; i<j; j--)
                    {
                      if ( abs(ev(j)) > abs(ev(j-1)) )
                        {
                          help = ev(j);
                          ev(j) = ev(j-1);
                          ev(j-1) = help;
                        }
                    }
                }
              val = (ev(0));
              break;
            }
          }
        return val;
      }

    return values[comp-1];
  }

  complex<double> VisualSceneSolution ::  ExtractValueComplex (const SolData * data, int comp, double * values) const
  {
    if (!data->iscomplex)
      return values[comp-1];
    else
      return complex<double> (values[comp-1], values[comp]);
  }




  bool VisualSceneSolution :: 
  GetSurfValueComplex (const SolData * data, SurfaceElementIndex selnr, int facetnr, 
                       double lam1, double lam2, 
                       int comp, complex<double> & val) const
  {
    switch (data->soltype)
      {
      case SOL_VIRTUALFUNCTION:
        {
          ArrayMem<double,20> values(data->components);
          bool ok;
          
          ok = data->solclass->GetSurfValue (selnr, facetnr, lam1, lam2, &values[0]);
          
          if (ok)
            {
              if (!data->iscomplex)
                val = values[comp-1];
              else
                val = complex<double> (values[comp-1], values[comp]);
            }
          
          return ok;
        }
      default:
        cerr << "case not implementd 6565" << endl;
      }
    return 0;
  }
  
  bool VisualSceneSolution :: 
  GetSurfValue (const SolData * data, SurfaceElementIndex selnr, int facetnr, 
                double lam1, double lam2, 
                int comp, double & val) const
  {
    shared_ptr<Mesh> mesh = GetMesh();

    bool ok;
    if (comp == 0)
      {
        val = 0;
        ArrayMem<double,20> values(data->components);
        ok = GetSurfValues (data, selnr, facetnr, lam1, lam2, &values[0]);
	val = ExtractValue (data, 0, &values[0]);
	return ok;
      }


    switch (data->soltype)
      {
      case SOL_VIRTUALFUNCTION:
        {
  
          ArrayMem<double,20> values(data->components);
          bool ok;

          ok = data->solclass->GetSurfValue (selnr, facetnr, lam1, lam2, &values[0]);

          if (ok)
            {
              if (!data->iscomplex)
                val =  values[comp-1];
              else
                {
                  // cout << "time = " << time << ", cos = " << cos(time) << endl;
     
                  // old version: val = values[comp-1]*cos(3*time) + values[comp]*sin(3*time);
                  // SZ: Sept 06 
                  if(comp%2==0) 
                    val =  values[comp-1]*cos(3*time) - values[comp-2]*sin(3*time);
                  else
                    val = values[comp-1]*cos(3*time) + values[comp]*sin(3*time);
         
         
         
                }
            }

          return ok;
        }


      case SOL_NODAL:
        {
          const Element2d & el = (*mesh)[selnr];

          double lami[8];
          int np, i;
          val = 0;
          double lam3 = 1-lam1-lam2;

          switch (el.GetType())
            {
            case TRIG:
              /*
                lami[0] = lam3;
                lami[1] = lam1;
                lami[2] = lam2;
              */
              lami[0] = lam1;
              lami[1] = lam2;
              lami[2] = lam3;
              np = 3;
              break;

            case TRIG6:
              /*
                lami[0] = lam3*(2*lam3-1);
                lami[1] = lam1*(2*lam1-1);
                lami[2] = lam2*(2*lam2-1);
              */
              // hierarchical basis:
              lami[0] = lam3;
              lami[1] = lam1;
              lami[2] = lam2;
              lami[3] = 4*lam1*lam2;
              lami[4] = 4*lam2*lam3;
              lami[5] = 4*lam1*lam3;
              np = 6;
              break;

            case QUAD:
            case QUAD6:
              lami[0] = (1-lam1)*(1-lam2);
              lami[1] = lam1 * (1-lam2);
              lami[2] = lam1 * lam2;
              lami[3] = (1-lam1) * lam2;
              np = 4;
              break;

            default:
              np = 0;
            }

          for (i = 0; i < np; i++)
            val += lami[i] * data->data[(el[i]-1) * data->dist + comp-1];

          return 1;
        }

      case SOL_ELEMENT:
        {
          int el1, el2;
          mesh->GetTopology().GetSurface2VolumeElement (selnr+1, el1, el2);
          el1--;

          val = data->data[el1 * data->dist+comp-1];
          return 1;
        }

      case SOL_NONCONTINUOUS:
        {
          val = 0;
          // ?????
          return 0;
        }

      case SOL_SURFACE_ELEMENT:
        {
          val = data->data[selnr * data->dist + comp-1];
          return 1;
        }

      case SOL_SURFACE_NONCONTINUOUS:
        {
          const Element2d & el = (*mesh)[selnr];

          double lami[8];
          int np = 0;
          val = 0;
          int order = data->order;

          switch (order)
            {
            case 0:
              return data->data[selnr * data->dist + comp-1];
            case 1:
              {
                switch (el.GetType())
                  {
                  case TRIG:
                  case TRIG6:
                    {
                      lami[1] = lam1;
                      lami[2] = lam2;
                      lami[0] = 1-lam1-lam2;
                      np = 3;
                      break;
                    }
                  default:
                    cerr << "case not implementd 2342" << endl;
                  }
                break;
              }
            case 2:
              {
                switch (el.GetType())
                  {
                  case TRIG:
                    {
                      lami[1] = lam1;
                      lami[2] = lam2;
                      lami[0] = 1-lam1-lam2;
                      np = 3;
                      break;
                    }
                  case TRIG6:
                    {
                      double lam3 = 1-lam1-lam2;
                      lami[1] = 2*lam1 * (lam1-0.5);
                      lami[2] = 2*lam2 * (lam2-0.5);
                      lami[0] = 2*lam3 * (lam3-0.5);
                      lami[3] = 4*lam1*lam2;
                      lami[4] = 4*lam2*lam3;
                      lami[5] = 4*lam1*lam3;
                      np = 6;
                      break;
                    }
                  default:
                    cerr << "case not implemented 8712" << endl;
                  }
                break;
              }
            }
        
          int base;
          if (order == 1)
            base = 4 * selnr;
          else 
            base = 9 * selnr;

          for (int i = 0; i < np; i++)
            val += lami[i] * data->data[(base+i) * data->dist + comp-1];

          return 1;
        }

      case SOL_MARKED_ELEMENTS:
        {
          val = (*mesh)[selnr].TestRefinementFlag();
          return 1;
        }
      
      case SOL_ELEMENT_ORDER:
        {       
          val = (*mesh)[selnr].GetOrder();
          return 1;
        }

      }
    return 0;
  }












  bool VisualSceneSolution :: 
  GetSurfValue (const SolData * data, SurfaceElementIndex selnr, int facetnr, 
                const double xref[], const double x[], const double dxdxref[], 
                int comp, double & val) const
  {
    shared_ptr<Mesh> mesh = GetMesh();

    double lam1 = xref[0], lam2 = xref[1];

    bool ok;
    if (comp == 0)
      {
        val = 0;
        ArrayMem<double,20> values(data->components);
        ok = GetSurfValues (data, selnr, facetnr, xref, x, dxdxref, &values[0]);
	val = ExtractValue (data, 0, &values[0]);
	return ok;
      }


    switch (data->soltype)
      {
      case SOL_VIRTUALFUNCTION:
        {
          ArrayMem<double,20> values(data->components);
          bool ok;

          // ok = data->solclass->GetSurfValue (selnr, lam1, lam2, &values[0]);
          // cout << "data->solclass = " << flush << data->solclass << endl;
          ok = data->solclass->GetSurfValue (selnr, facetnr, xref, x, dxdxref, &values[0]);
          // ok = 1;
          // values[0] = 1.0;

          if (ok)
            {
              if (!data->iscomplex)
                val =  values[comp-1];
              else
                {
                  // cout << "time = " << time << ", cos = " << cos(time) << endl;
     
                  // old version: val = values[comp-1]*cos(3*time) + values[comp]*sin(3*time);
                  // SZ: Sept 06 
                  if(comp%2==0) 
                    val =  values[comp-1]*cos(3*time) - values[comp-2]*sin(3*time);
                  else
                    val = values[comp-1]*cos(3*time) + values[comp]*sin(3*time);
                  
                }
            }

          return ok;
        }


      case SOL_NODAL:
        {
          const Element2d & el = (*mesh)[selnr];

          double lami[8];
          int np, i;
          val = 0;
          double lam3 = 1-lam1-lam2;

          switch (el.GetType())
            {
            case TRIG:
              /*
                lami[0] = lam3;
                lami[1] = lam1;
                lami[2] = lam2;
              */
              lami[0] = lam1;
              lami[1] = lam2;
              lami[2] = lam3;
              np = 3;
              break;

            case TRIG6:
              /*
                lami[0] = lam3*(2*lam3-1);
                lami[1] = lam1*(2*lam1-1);
                lami[2] = lam2*(2*lam2-1);
              */
              // hierarchical basis:
              lami[0] = lam3;
              lami[1] = lam1;
              lami[2] = lam2;
              lami[3] = 4*lam1*lam2;
              lami[4] = 4*lam2*lam3;
              lami[5] = 4*lam1*lam3;
              np = 6;
              break;

            case QUAD:
            case QUAD6:
              lami[0] = (1-lam1)*(1-lam2);
              lami[1] = lam1 * (1-lam2);
              lami[2] = lam1 * lam2;
              lami[3] = (1-lam1) * lam2;
              np = 4;
              break;

            default:
              np = 0;
            }

          for (i = 0; i < np; i++)
            val += lami[i] * data->data[(el[i]-1) * data->dist + comp-1];

          return 1;
        }

      case SOL_ELEMENT:
        {
          int el1, el2;
          mesh->GetTopology().GetSurface2VolumeElement (selnr+1, el1, el2);
          el1--;

          val = data->data[el1 * data->dist+comp-1];
          return 1;
        }

      case SOL_NONCONTINUOUS:
        {
          val = 0;
          // ?????
          return 0;
        }

      case SOL_SURFACE_ELEMENT:
        {
          val = data->data[selnr * data->dist + comp-1];
          return 1;
        }

      case SOL_SURFACE_NONCONTINUOUS:
        {
          const Element2d & el = (*mesh)[selnr];

          double lami[8] = { 0.0 };
          int np = 0;
          val = 0;
          int order = data->order;

          switch (order)
            {
            case 0:
              return data->data[selnr * data->dist + comp-1];
            case 1:
              {
                switch (el.GetType())
                  {
                  case TRIG:
                  case TRIG6:
                    {
                      lami[1] = lam1;
                      lami[2] = lam2;
                      lami[0] = 1-lam1-lam2;
                      np = 3;
                      break;
                    }
                  default:
                    cerr << "case not impl 234234" << endl;
                  }
                break;
              }
            case 2:
              {
                switch (el.GetType())
                  {
                  case TRIG:
                    {
                      lami[1] = lam1;
                      lami[2] = lam2;
                      lami[0] = 1-lam1-lam2;
                      np = 3;
                      break;
                    }
                  case TRIG6:
                    {
                      double lam3 = 1-lam1-lam2;
                      lami[1] = 2*lam1 * (lam1-0.5);
                      lami[2] = 2*lam2 * (lam2-0.5);
                      lami[0] = 2*lam3 * (lam3-0.5);
                      lami[3] = 4*lam1*lam2;
                      lami[4] = 4*lam2*lam3;
                      lami[5] = 4*lam1*lam3;
                      np = 6;
                      break;
                    }
                  default:
                    cerr << "case not implented 3234" << endl;
                  }
                break;
              }
            }
        
          int base;
          if (order == 1)
            base = 4 * selnr;
          else 
            base = 9 * selnr;

          for (int i = 0; i < np; i++)
            val += lami[i] * data->data[(base+i) * data->dist + comp-1];

          return 1;
        }

      case SOL_MARKED_ELEMENTS:
        {
          val = (*mesh)[selnr].TestRefinementFlag();
          return 1;
        }
      
      case SOL_ELEMENT_ORDER:
        {       
          val = (*mesh)[selnr].GetOrder();
          return 1;
        }

      }
    return 0;
  }









  Vec<3> VisualSceneSolution :: 
  GetDeformation (ElementIndex elnr, const Point<3> & p) const
  {
    Vec<3> def;
    if (deform && vecfunction != -1)
      {
        GetValues (soldata[vecfunction], elnr, p(0), p(1), p(2), &def(0));
        def *= scaledeform;

        if (soldata[vecfunction]->components == 2) def(2) = 0;
      }
    else
      def = 0;
    return def;
  }


  Vec<3> VisualSceneSolution :: 
  GetSurfDeformation (SurfaceElementIndex elnr, int facetnr, double lam1, double lam2) const
  {
    shared_ptr<Mesh> mesh = GetMesh();

    Vec<3> def;
    if (deform && vecfunction != -1)
      {
        // GetSurfValues (soldata[vecfunction], elnr, facetnr, lam1, lam2,  &def(0));
	double values[6];
	GetSurfValues (soldata[vecfunction], elnr, facetnr, lam1, lam2,  values);
	def = RealVec3d (values, soldata[vecfunction]->iscomplex, imag_part);
        def *= scaledeform;

        if (soldata[vecfunction]->components == 2) def(2) = 0;
      }
    else if (deform && scalfunction != -1 && mesh->GetDimension()==2)
      { // he: allow for 3d plots of 2d surfaces: usage: turn deformation on
        def = 0;
        GetSurfValue (soldata[scalfunction], elnr, facetnr, lam1, lam2, scalcomp, def(2));
        def *= scaledeform;
      }
    else
      def = 0;
    return def;
  }

  void VisualSceneSolution :: GetPointDeformation (int pnum, Point<3> & p, 
                                                   SurfaceElementIndex elnr) const
  {
    shared_ptr<Mesh> mesh = GetMesh();

    p = mesh->Point (pnum+1);
    if (deform && vecfunction != -1)
      {
        const SolData * vsol = soldata[vecfunction];
      
        Vec<3> v(0,0,0);
        if (vsol->soltype == SOL_NODAL)
          {
            v = Vec3d(vsol->data[pnum * vsol->dist],
                      vsol->data[pnum * vsol->dist+1],
                      vsol->data[pnum * vsol->dist+2]);
          }
        else if (vsol->soltype == SOL_SURFACE_NONCONTINUOUS)
          {
            const Element2d & el = (*mesh)[elnr];
            for (int j = 0; j < el.GetNP(); j++)
              if (el[j] == pnum+1)
                {
                  int base = (4*elnr+j-1) * vsol->dist;
                  v = Vec3d(vsol->data[base],
                            vsol->data[base+1],
                            vsol->data[base+2]);
                }
          }

        if (vsol->dist == 2) v(2) = 0;
      
        v *= scaledeform;
        p += v;
      }
  }




  void VisualSceneSolution :: GetClippingPlaneTrigs (Array<ClipPlaneTrig> & trigs,
                                                     Array<ClipPlanePoint> & pts)
  {
    shared_ptr<Mesh> mesh = GetMesh();

    static int timer_vals = NgProfiler::CreateTimer ("ClipPlaneTrigs - vertex values");
    static int timer1 = NgProfiler::CreateTimer ("ClipPlaneTrigs1");
    // static int timer1a = NgProfiler::CreateTimer ("ClipPlaneTrigs1a");
    // static int timer2 = NgProfiler::CreateTimer ("ClipPlaneTrigs2");
    static int timer3 = NgProfiler::CreateTimer ("ClipPlaneTrigs3");
    static int timer4 = NgProfiler::CreateTimer ("ClipPlaneTrigs4");
    static int timer4b = NgProfiler::CreateTimer ("ClipPlaneTrigs4b");


    NgProfiler::RegionTimer reg1 (timer1);

    
    int ne = mesh->GetNE();

    const int edgei[6][2] =
      { { 0, 1 }, { 0, 2 }, { 0, 3 },
        { 1, 2 }, { 1, 3 }, { 2, 3 } };

    double edgelam[6];
    // Point<3> edgep[6];
    double nodevali[4];

    int cntce;
    int cpe1 = 0, cpe2 = 0, cpe3 = 0;

    // Array<Element> loctets;
    // Array<Element> loctetsloc;
    // Array<Point<3> > pointsloc;

    int n = 1 << subdivisions;
    int n3 = (n+1)*(n+1)*(n+1);

    Array<Point<3> > grid(n3);
    Array<Point<3> > locgrid(n3);
    Array<Mat<3,3> > trans(n3);
    Array<double> val(n3);
    Array<bool> locposval(n3);
    Array<int> compress(n3);

    NgProfiler::StartTimer (timer_vals);
    Array<double,PointIndex::BASE> vertval(mesh->GetNV());
    Array<bool,PointIndex::BASE> posval(mesh->GetNV());
    for (PointIndex pi = vertval.Begin(); pi < vertval.End(); pi++)
      {
        Point<3> vert = (*mesh)[pi];
        vertval[pi] = 
          vert(0) * clipplane[0] + 
          vert(1) * clipplane[1] + 
          vert(2) * clipplane[2] + 
          clipplane[3]; 
        posval[pi] = vertval[pi] > 0;
      }
    NgProfiler::StopTimer (timer_vals);

    INDEX_2_CLOSED_HASHTABLE<int> edges(8*n3);  // point nr of edge
    

    for (ElementIndex ei = 0; ei < ne; ei++)
      {
        // NgProfiler::RegionTimer reg1a (timer1a);
        int first_point_of_element = pts.Size();

	locgrid.SetSize(n3);
        if(vispar.clipdomain > 0 && vispar.clipdomain != (*mesh)[ei].GetIndex()) continue;
        if(vispar.donotclipdomain > 0 && vispar.donotclipdomain == (*mesh)[ei].GetIndex()) continue;

        ELEMENT_TYPE type = (*mesh)[ei].GetType();
        if (type == HEX || type == PRISM || type == TET || type == TET10 || type == PYRAMID)
          {
            const Element & el = (*mesh)[ei];

            int ii = 0;
            int cnt_valid = 0;

            // NgProfiler::StartTimer (timer2);

            
            if (!mesh->GetCurvedElements().IsElementHighOrder(ei))
              {
                bool has_pos = 0, has_neg = 0;
                
                for (int i = 0; i < el.GetNP(); i++)
                  if (posval[el[i]])
                    has_pos = 1;
                  else
                    has_neg = 1;
                
                if (!has_pos || !has_neg) 
                  {
                    // NgProfiler::StopTimer (timer2);
                    continue;
                  }
              }
            


            if (type == TET || type == TET10)
              {
                for (int ix = 0; ix <= n; ix++)
                  for (int iy = 0; iy <= n; iy++)
                    for (int iz = 0; iz <= n; iz++, ii++)
                      {
                        if (ix+iy+iz <= n)
                          {
                            compress[ii] = cnt_valid;
                            locgrid[cnt_valid] = 
                              Point<3> (double(ix) / n, double(iy) / n, double(iz) / n);
                            cnt_valid++;
                          }
                        else
                          compress[ii] = -1;
                      }
              }
            
            else
              
              for (int ix = 0; ix <= n; ix++)
                for (int iy = 0; iy <= n; iy++)
                  for (int iz = 0; iz <= n; iz++, ii++)
                    {
                      Point<3> ploc;
                      compress[ii] = ii;
                      
                      switch (type)
                        {
                        case PRISM:
                          if (ix+iy <= n)
                            {
                              ploc = Point<3> (double(ix) / n, double(iy) / n, double(iz) / n);
                              compress[ii] = cnt_valid;
                              cnt_valid++;
                            }
                          else
                            compress[ii] = -1;
                          break;
                        case HEX:
                          ploc = Point<3> (double(ix) / n, double(iy) / n, double(iz) / n);
                          break;
                        case PYRAMID:
                          ploc = Point<3> (double(ix) / n * (1-double(iz)/n),
                                           double(iy) / n * (1-double(iz)/n),
                                           double(iz)/n);
                          if (iz == n) ploc = Point<3> (0,0,1-1e-8);
                          break;
                        default:
                          cerr << "clip plane trigs not implemented" << endl;
                          ploc = Point<3> (0,0,0);
                        }
                      if (compress[ii] != -1)
                        locgrid[compress[ii]] = ploc;
                    }

            if (type != TET && type != TET10 && type != PRISM) cnt_valid = n3;

	    locgrid.SetSize(cnt_valid);

            // NgProfiler::StopTimer (timer2);
            // NgProfiler::RegionTimer reg4(timer4);

            if (mesh->GetCurvedElements().IsHighOrder())
              {
                NgProfiler::RegionTimer reg4(timer4);
                mesh->GetCurvedElements().
                  CalcMultiPointElementTransformation (&locgrid, ei, &grid, 0);
              }
            else
              {
                NgProfiler::RegionTimer reg4(timer4b);
                Vector shape(el.GetNP());
                MatrixFixWidth<3> pointmat(el.GetNP());

                for (int k = 0; k < el.GetNP(); k++)
                  for (int j = 0; j < 3; j++)
                    pointmat(k,j) = (*mesh)[el[k]](j);
                
                for (int i = 0; i < cnt_valid; i++)
                  {
                    el.GetShapeNew (locgrid[i], shape);
                    Point<3> pglob;
                    for (int j = 0; j < 3; j++)
                      {
                        pglob(j) = 0;
                        for (int k = 0; k < el.GetNP(); k++)
                          pglob(j) += shape(k) * pointmat(k,j);
                      }
                    grid[i] = pglob;
                  }
              }

            NgProfiler::RegionTimer reg3(timer3);

            bool has_pos = false, all_pos = true;

            for (int i = 0; i < cnt_valid; i++)
              {
                val[i] = 
                  grid[i](0) * clipplane[0] + 
                  grid[i](1) * clipplane[1] + 
                  grid[i](2) * clipplane[2] + 
                  clipplane[3];
                
                locposval[i] = val[i] > 0;
                has_pos |= locposval[i];
                all_pos &= locposval[i];
                
                // if (val[i] > 0) has_pos = 1; else has_neg = 1;
              }
                
            // if (!has_pos || !has_neg) continue;
            if (!has_pos || all_pos) continue;
            
            edges.DeleteData();
            
            for (int ix = 0; ix < n; ix++)
              for (int iy = 0; iy < n; iy++)
                for (int iz = 0; iz < n; iz++)
                  {
                    int base = iz + (n+1)*iy + (n+1)*(n+1)*ix;
                    int pi[8] = 
                      { base, base+(n+1)*(n+1), base+(n+1)*(n+1)+(n+1), base+(n+1),
                        base+1, base+(n+1)*(n+1)+1, base+(n+1)*(n+1)+(n+1)+1, base+(n+1)+1 };

                    for (int j = 0; j < 8; j++)
                      pi[j] = compress[pi[j]];

                    bool has_pos = false, all_pos = true;
                    for (int j = 0; j < 8; j++)
                      if (pi[j] != -1)
                        {
                          has_pos |= locposval[pi[j]];
                          all_pos &= locposval[pi[j]];
                        }
                    if (!has_pos || all_pos) continue;


                    const int tets[6][4] = 
                      { { 1, 2, 4, 5 },
                        { 4, 5, 2, 8 },
                        { 2, 8, 5, 6 },
                        { 2, 3, 4, 8 },
                        { 2, 3, 8, 6 },
                        { 3, 8, 6, 7 } };
                    
                    for (int ii = 0; ii < 6; ii++)
                      {
                        int teti[4];
                        for (int k = 0; k < 4; k++)
                          teti[k] = pi[tets[ii][k]-1];

                        bool is_valid = true;
                        for (int j = 0; j < 4; j++)
                          is_valid &= (teti[j] != -1);
                        if (!is_valid) continue;
                        
                        bool has_pos = false, all_pos = true;
                        for (int j = 0; j < 4; j++)
                          {
                            has_pos |= locposval[teti[j]];
                            all_pos &= locposval[teti[j]];
                          }
                        if (!has_pos || all_pos) continue;

                        for (int j = 0; j < 4; j++)
                          nodevali[j] = val[teti[j]];
          
                        cntce = 0;
                        for (int j = 0; j < 6; j++)
                          {
                            int lpi1 = edgei[j][0];
                            int lpi2 = edgei[j][1];
                            if ( (nodevali[lpi1] > 0) != (nodevali[lpi2] > 0) )
                              {
                                cntce++;
                                cpe3 = cpe2;
                                cpe2 = cpe1;
                                cpe1 = j;
                                if (cntce >= 3)
                                  {
                                    ClipPlaneTrig cpt;
                                    cpt.elnr = ei;
                                  
                                    for (int k = 0; k < 3; k++)
                                      {
                                        int ednr;
                                        switch (k)
                                          {
                                          case 0: ednr = cpe1; break;
                                          case 1: ednr = cpe2; break;
                                          case 2: ednr = cpe3; break;
                                          }

                                        int pi1 = edgei[ednr][0];
                                        int pi2 = edgei[ednr][1];
                                        int pnr = -1;
                                        
                                        INDEX_2 pair (teti[pi1], teti[pi2]);
                                        pair.Sort();
                                        if (edges.Used(pair))
                                          pnr = edges.Get(pair);
                                        else
                                          {
                                            ClipPlanePoint cppt;
                                            cppt.elnr = ei;
                                            
                                            edgelam[ednr] = nodevali[pi2] / (nodevali[pi2] - nodevali[pi1]);

                                            Point<3> gp1 = grid[teti[pi1]];
                                            Point<3> gp2 = grid[teti[pi2]];
                                            cppt.p = gp2 + edgelam[ednr] * (gp1-gp2);
                                            
                                            Point<3> p1 = locgrid[teti[pi1]];
                                            Point<3> p2 = locgrid[teti[pi2]];
                                            cppt.lami =  p2 + edgelam[ednr] * (p1-p2);
                                           
 
                                            pnr = pts.Append (cppt)-1;
                                            edges.Set (pair, pnr);
                                          }

                                        cpt.points[k].pnr = pnr;
                                        cpt.points[k].locpnr = pnr-first_point_of_element;
                                      }
                                  
                                    trigs.Append (cpt);
                                  }
                              }
                          }
                      }
                  }
          }

        else
          {  // other elements not supported (JS, June 2007)
            continue; // return;
          }
      
      }
  }

  void VisualSceneSolution :: GetClippingPlaneGrid (Array<ClipPlanePoint> & pts)
  {
    shared_ptr<Mesh> mesh = GetMesh();

    Vec3d n(clipplane[0], clipplane[1], clipplane[2]);

    double mu = -clipplane[3] / n.Length2();
    Point3d p(mu*n.X(), mu * n.Y(), mu * n.Z());

    // n /= n.Length();
    n.Normalize();
    Vec3d t1, t2;
    n.GetNormal (t1);
    t2 = Cross (n, t1);

    double xi1, xi2;

    double xi1mid = (center - p) * t1;
    double xi2mid = (center - p) * t2;

    pts.SetSize(0);

    for (xi1 = xi1mid-rad+xoffset/gridsize; xi1 <= xi1mid+rad+xoffset/gridsize; xi1 += rad / gridsize)
      for (xi2 = xi2mid-rad+yoffset/gridsize; xi2 <= xi2mid+rad+yoffset/gridsize; xi2 += rad / gridsize)
        {
          Point3d hp = p + xi1 * t1 + xi2 * t2;
        
          int cindex(-1);
          bool allowindex(true);
          if(vispar.clipdomain > 0)
            {
              cindex = vispar.clipdomain;
            }
          else if(vispar.donotclipdomain > 0)
            {
              allowindex = false;
              cindex = vispar.donotclipdomain;
            }

          double lami[3];
          int elnr = mesh->GetElementOfPoint (hp, lami,0,cindex,allowindex)-1;

          if (elnr != -1)
            {
              ClipPlanePoint cpp;
              cpp.p = hp;
              cpp.elnr = elnr;
              cpp.lami(0) = lami[0];
              cpp.lami(1) = lami[1];
              cpp.lami(2) = lami[2];
              pts.Append (cpp);
            }
        }
  };




  void VisualSceneSolution :: DrawClipPlaneTrigs () 
  {
    shared_ptr<Mesh> mesh = GetMesh();

#ifdef PARALLELGL

    if (id == 0 && ntasks > 1)
      {
	InitParallelGL();

	Array<int> parlists (ntasks);

	MyMPI_SendCmd ("redraw");
	MyMPI_SendCmd ("clipplanetrigs");

	for ( int dest = 1; dest < ntasks; dest++ )
	  MyMPI_Recv (parlists[dest], dest, MPI_TAG_VIS);

	if (clipplanelist_scal)
	  glDeleteLists (clipplanelist_scal, 1);

	clipplanelist_scal = glGenLists (1);
	glNewList (clipplanelist_scal, GL_COMPILE);
	
	for ( int dest = 1; dest < ntasks; dest++ )
	  glCallList (parlists[dest]);
	
	glEndList();
	return;
      }
#endif





    if (clipplanelist_scal)
      glDeleteLists (clipplanelist_scal, 1);
    
    clipplanelist_scal = glGenLists (1);
    glNewList (clipplanelist_scal, GL_COMPILE);


    Array<ClipPlaneTrig> trigs;
    Array<ClipPlanePoint> points;
    GetClippingPlaneTrigs (trigs, points);
	    
    glNormal3d (-clipplane[0], -clipplane[1], -clipplane[2]);
    glColor3d (1.0, 1.0, 1.0);
    
    SetTextureMode (usetexture);

    SolData * sol = NULL;

    if (scalfunction != -1) 
      sol = soldata[scalfunction];

    if (sol -> draw_volume)
      {
	glBegin (GL_TRIANGLES);

    int maxlpnr = 0;
    for (int i = 0; i < trigs.Size(); i++)
      for (int j = 0; j < 3; j++)
        maxlpnr = max2 (maxlpnr, trigs[i].points[j].locpnr);

    Array<double> vals(maxlpnr+1);
    Array<complex<double> > valsc(maxlpnr+1);
    Array<int> elnrs(maxlpnr+1);
    Array<bool> trigok(maxlpnr+1);
    Array<Point<3> > locpoints(maxlpnr+1);
    Array<Point<3> > globpoints(maxlpnr+1);
    Array<Mat<3> > jacobi(maxlpnr+1);
    Array<double> mvalues( (maxlpnr+1) * sol->components);
    trigok = false;
    elnrs = -1;

    Point<3> p[3];
    // double val[3];
    // complex<double> valc[3];
    int lastelnr = -1;
    int nlp = -1;
    bool ok = false;

    for (int i = 0; i < trigs.Size(); i++)
      {
        const ClipPlaneTrig & trig = trigs[i];
	if (trig.elnr != lastelnr)
	  {
	    lastelnr = trig.elnr;
	    nlp = -1;

	    for (int ii = i; ii < trigs.Size(); ii++)
	      {
		if (trigs[ii].elnr != trig.elnr) break;
		for (int j = 0; j < 3; j++)
		  nlp = max (nlp, trigs[ii].points[j].locpnr);
	      }
	    nlp++;
	    locpoints.SetSize (nlp);

	    for (int ii = i; ii < trigs.Size(); ii++)
	      {
		if (trigs[ii].elnr != trig.elnr) break;
		for (int j = 0; j < 3; j++)
		  locpoints[trigs[ii].points[j].locpnr] = points[trigs[ii].points[j].pnr].lami;
	      }

	    mesh->GetCurvedElements().
	      CalcMultiPointElementTransformation (&locpoints, trig.elnr, 
						   &globpoints, &jacobi);

	    bool
	      drawelem = GetMultiValues (sol, trig.elnr, -1, nlp, 
					 &locpoints[0](0), &locpoints[1](0)-&locpoints[0](0),
					 &globpoints[0](0), &globpoints[1](0)-&globpoints[0](0),
					 &jacobi[0](0), &jacobi[1](0)-&jacobi[0](0),
					 &mvalues[0], sol->components);
	    
	    // cout << "have multivalues, comps = " << sol->components << endl;

	    // if (!drawelem) ok = false;
	    ok = drawelem;
	    if (usetexture != 2 || !sol->iscomplex)
	      for (int ii = 0; ii < nlp; ii++)
		vals[ii] = ExtractValue(sol, scalcomp, &mvalues[ii*sol->components]);
	    else
	      for (int ii = 0; ii < nlp; ii++)
		valsc[ii] = complex<double> (mvalues[ii*sol->components + scalcomp-1],
					     mvalues[ii*sol->components + scalcomp]);
	  }
	
	if(ok)
	  for(int j=0; j<3; j++)
	    {
	      if (usetexture != 2 || !sol->iscomplex)
		SetOpenGlColor (vals[trig.points[j].locpnr]);
	      else
		glTexCoord2f ( valsc[trig.points[j].locpnr].real(), 
			       valsc[trig.points[j].locpnr].imag() );

	      p[j] = points[trig.points[j].pnr].p;

	      if (deform)
		{
		  Point<3> ploc = points[trig.points[j].pnr].lami;
		  p[j] += GetDeformation (trig.elnr, ploc);
		}

	      glVertex3dv (p[j]);
	    }

      }
    glEnd();
      }
    glEndList ();


#ifdef PARALLELGL
    glFinish();
    if (id > 0)
      MyMPI_Send (clipplanelist_scal, 0, MPI_TAG_VIS);
#endif
  }










  void VisualSceneSolution ::
  SetOpenGlColor(double val)
  {
    if (usetexture == 1 && !logscale)
      {
        glTexCoord1f ( val );
        return;
      }

    double valmin = minval;
    double valmax = maxval;

    double value;

    if (!logscale)
      value = (val - valmin) / (valmax - valmin);
    else
      {
        if (valmax <= 0) valmax = 1;
        if (valmin <= 0) valmin = 1e-4 * valmax;
        value = (log(fabs(val)) - log(valmin)) / (log(valmax) - log(valmin));
      }

    if (!invcolor)
      value = 1 - value;


    if (value > 1) value = 1;
    if (value < 0) value = 0;

    value *= 4;

    static const double colp[][3] =
      {
        { 1, 0, 0 },
        { 1, 1, 0 },
        { 0, 1, 0 },
        { 0, 1, 1 },
        { 0, 0, 1 },
        { 1, 0, 1 },
        { 1, 0, 0 },
      };
  
    int i = int(value);
    double r = value - i;

    GLdouble col[3];
    for (int j = 0; j < 3; j++)
      col[j] = (1-r) * colp[i][j] + r * colp[i+1][j];
  
    glColor3dv (col);
  }



  void VisualSceneSolution ::
  SetTextureMode (int texturemode) const
  {
    switch (texturemode)
      {
      case 0:
        glDisable (GL_TEXTURE_1D);
        glDisable (GL_TEXTURE_2D);
        break;
      case 1:
        glEnable (GL_TEXTURE_1D);
        glDisable (GL_TEXTURE_2D);
        glColor3d (1.0, 1.0, 1.0);   
        break;
      case 2:
        glDisable (GL_TEXTURE_1D);
        glEnable (GL_TEXTURE_2D);
        glColor3d (1.0, 1.0, 1.0);   
        break;
      }
  }




  void VisualSceneSolution ::
  DrawCone (const Point<3> & p1, const Point<3> & p2, double r)
  {
    int n = 10, i;
    Vec<3> p1p2 = p2 - p1;

    p1p2.Normalize();
    Vec<3> p2p1 = -p1p2;

    Vec<3> t1 = p1p2.GetNormal();
    Vec<3> t2 = Cross (p1p2, t1);

    Point<3> oldp = p1 + r * t1;
    Vec<3> oldn = t1;

    Point<3> p;
    Vec<3> normal;

    Mat<2> rotmat;
    Vec<2> cs, newcs;
    cs(0) = 1;
    cs(1) = 0;
    rotmat(0,0) = rotmat(1,1) = cos(2*M_PI/n);
    rotmat(1,0) = sin(2*M_PI/n);
    rotmat(0,1) = -rotmat(1,0);

    glBegin (GL_TRIANGLES);
    for (i = 1; i <= n; i++)
      {
        /*
          phi = 2 * M_PI * i / n;
          normal = cos(phi) * t1 + sin(phi) * t2;
        */
        newcs = rotmat * cs;
        cs = newcs;
        normal = cs(0) * t1 + cs(1) * t2;

        p = p1 + r * normal;

        // cone
        glNormal3dv (normal);
        glVertex3dv (p);
        glVertex3dv (p2);
        glNormal3dv (oldn);
        glVertex3dv (oldp);

        // base-circle
        glNormal3dv (p2p1);
        glVertex3dv (p);
        glVertex3dv (p1);
        glVertex3dv (oldp);

        oldp = p;
        oldn = normal;
      }
    glEnd ();
  }



  void VisualSceneSolution ::
  DrawCylinder (const Point<3> & p1, const Point<3> & p2, double r)
  {
    int n = 10, i;
    Vec<3> p1p2 = p2 - p1;

    p1p2.Normalize();
    // Vec<3> p2p1 = -p1p2;

    Vec<3> t1 = p1p2.GetNormal();
    Vec<3> t2 = Cross (p1p2, t1);

    Point<3> oldhp1 = p1 + r * t1;
    Point<3> oldhp2 = p2 + r * t1;
    Vec<3> oldn = t1;

    Point<3> hp1, hp2;
    Vec<3> normal;

    Mat<2> rotmat;
    Vec<2> cs, newcs;
    cs(0) = 1;
    cs(1) = 0;
    rotmat(0,0) = rotmat(1,1) = cos(2*M_PI/n);
    rotmat(1,0) = sin(2*M_PI/n);
    rotmat(0,1) = -rotmat(1,0);

    glBegin (GL_QUADS);
    for (i = 1; i <= n; i++)
      {
        newcs = rotmat * cs;
        cs = newcs;
        normal = cs(0) * t1 + cs(1) * t2;

        hp1 = p1 + r * normal;
        hp2 = p2 + r * normal;

        // cylinder
        glNormal3dv (normal);

        glVertex3dv (hp1);
        glVertex3dv (hp2);
        glVertex3dv (oldhp2);
        glVertex3dv (oldhp1);

        oldhp1 = hp1;
        oldhp2 = hp2;
        oldn = normal;
      }
    glEnd ();
  }













  void VisualSceneSolution :: MouseDblClick (int px, int py)
  {
    vsmesh.SetClippingPlane();
    // vsmesh.BuildFilledList();
    vsmesh.MouseDblClick(px,py);
  }



#ifdef PARALLELGL

  void VisualSceneSolution :: Broadcast ()
  {
    MPI_Datatype type;
    int blocklen[] = 
      { 
	1, 1, 1, 1,
	1, 1, 1, 1, 
	1, 1, 1, 1, 
	1, 4, 1, 1, 
	1
      };
    MPI_Aint displ[] = { (char*)&usetexture - (char*)this,
			 (char*)&clipsolution - (char*)this,
			 (char*)&scalfunction - (char*)this,
			 (char*)&scalcomp - (char*)this,

			 (char*)&vecfunction - (char*)this,
			 (char*)&gridsize - (char*)this,
			 (char*)&autoscale - (char*)this,
			 (char*)&logscale - (char*)this,

			 (char*)&minval - (char*)this,
			 (char*)&maxval - (char*)this,
			 (char*)&numisolines - (char*)this,
			 (char*)&subdivisions - (char*)this,

			 (char*)&evalfunc - (char*)this,
			 (char*)&clipplane[0] - (char*)this,
			 (char*)&multidimcomponent - (char*)this, 
			 (char*)&deform - (char*)this,

			 (char*)&scaledeform - (char*)this 
    };


    MPI_Datatype types[] = { 
      MPI_INT, MPI_INT, MPI_INT, MPI_INT,
      MPI_INT, MPI_INT, MPI_INT, MPI_INT,
      MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT,
      MPI_INT, MPI_DOUBLE, MPI_INT, MPI_INT,
      MPI_DOUBLE
    };

    MPI_Type_create_struct (17, blocklen, displ, types, &type);
    MPI_Type_commit ( &type );

    MPI_Bcast (this, 1, type, 0, MPI_COMM_WORLD);
    MPI_Type_free (&type);
  }
  
#endif

}


#include "../include/nginterface.h"

void Ng_ClearSolutionData ()
{
#ifdef OPENGL
  // if (nodisplay) return;
  netgen::vssolution.ClearSolutionData();
#endif
}

void Ng_InitSolutionData (Ng_SolutionData * soldata)
{
  soldata -> name = NULL;
  soldata -> data = NULL;
  soldata -> components = 1;
  soldata -> dist = 1;
  soldata -> order = 1;
  soldata -> iscomplex = 0;
  soldata -> draw_surface = 1;
  soldata -> draw_volume = 1;
  soldata -> soltype = NG_SOLUTION_NODAL;
  soldata -> solclass = 0;
}

void Ng_SetSolutionData (Ng_SolutionData * soldata)
{
#ifdef OPENGL
  // if (nodisplay) return;
  //   vssolution.ClearSolutionData ();
  netgen::VisualSceneSolution::SolData * vss = new netgen::VisualSceneSolution::SolData;

  vss->name = new char[strlen (soldata->name)+1];
  strcpy (vss->name, soldata->name);

  vss->data = soldata->data;
  vss->components = soldata->components;
  vss->dist = soldata->dist;
  vss->order = soldata->order;
  vss->iscomplex = bool(soldata->iscomplex);
  vss->draw_surface = soldata->draw_surface;
  vss->draw_volume = soldata->draw_volume;
  vss->soltype = netgen::VisualSceneSolution::SolType (soldata->soltype);
  vss->solclass = soldata->solclass;
  netgen::vssolution.AddSolutionData (vss);
#endif
}



namespace netgen
{
  extern void Render (bool blocking);
}

void Ng_Redraw (bool blocking)
{
#ifdef OPENGL
  netgen::vssolution.UpdateSolutionTimeStamp();
  netgen::Render(blocking);
#endif
}

