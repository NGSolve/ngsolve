/*! \file writeOpenFOAM15x.cpp
*  \brief Export Netgen Mesh in the OpenFOAM 1.5+ File format
*  \author Philippose Rajan
*  \date 25 October 2009
*
*  This function extends the export capabilities of
*  Netgen to include the OpenFOAM 1.5+ File Format.
*
*  The OpenFOAM 1.5+ mesh format consists of a set of 5 files 
*  which together define the mesh points, faces, cells and 
*  boundary conditions. 
*
*  The files are:
*  1. points    -> A list of the point co-ordinates
*  2. faces     -> A list of the faces with format <n>(pnt_ind1 pnt_ind2 .... pnt_ind<n>)
*  3. owner     -> The owner cell of each face 
*  4. neighbour -> The neighbour cell of each face
*  5. boundary  -> The set of boundaries with name, start face, and num. of faces
*
*  For a detailed description of the format, refer to the following link:
*  http://openfoamwiki.net/index.php/Write_OpenFOAM_meshes
*
*/

#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>
#include <sys/stat.h>


namespace netgen
{
#include "writeuser.hpp"

   // Global Va
   Array<int> OF15x_owner_facelist;
   Array<int> OF15x_owner_celllist;
   Array<int> OF15x_neighbour_celllist;
   Array<int> OF15x_surfelem_bclist;
   Array<INDEX_2> OF15x_surfelem_lists;


   void WriteOpenFOAM15xBanner(ofstream & outfile)
   {
      static char FOAMversion[4] = "1.5";
      static char spaces[40];

      memset(spaces, ' ', 40);
      spaces[38 - strlen(FOAMversion)] = '\0';
      
      outfile << 
              "/*--------------------------------*- C++ -*----------------------------------*\\\n";

      outfile <<
              "| =========                 |                                                 |\n"
              "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
              "|  \\\\    /   O peration     | Version:  " << FOAMversion << spaces << "|\n"
              "|   \\\\  /    A nd           | Web:      http://www.OpenFOAM.org               |\n"
              "|    \\\\/     M anipulation  |                                                 |\n"
              "\\*---------------------------------------------------------------------------*/\n";

   }



   void WriteOpenFOAM15xDividerStart(ofstream & outfile)
   {
      outfile  <<
               "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";
   }



   void WriteOpenFOAM15xDividerEnd(ofstream & outfile)
   {
      outfile <<
              "// ************************************************************************* //\n";
   }



   bool BuildOpenFOAM15xLists (const Mesh & mesh)
   {
      // Clear all the arrays
      OF15x_owner_facelist.DeleteAll();
      OF15x_owner_celllist.DeleteAll();
      OF15x_neighbour_celllist.DeleteAll();
      OF15x_surfelem_bclist.DeleteAll();
      OF15x_surfelem_lists.DeleteAll();


      const_cast<Mesh&> (mesh).BuildElementSearchTree();
      const MeshTopology& meshtopo = mesh.GetTopology();
      
      // Update the mesh topology structures
      const_cast<Mesh&> (mesh).UpdateTopology();

      int ne = mesh.GetNE();
      int nse = mesh.GetNSE();
      int totfaces = meshtopo.GetNFaces();

      // Preset the size of the arrays to speed up future operations
      // Number of internal faces = total faces - num. of surface faces
      OF15x_owner_facelist.SetSize(totfaces - nse);
      OF15x_owner_celllist.SetSize(totfaces - nse);
      OF15x_neighbour_celllist.SetSize(totfaces - nse);
      OF15x_surfelem_bclist.SetSize(nse);
      OF15x_surfelem_lists.SetSize(nse);

      // Array used to keep track of Faces which have already been 
      // processed and added to the Owner list... In addition, also the 
      // location where the face appears in the Owner list is also stored 
      // to speed up creation of the Neighbour list
      Array<int> ownerfaces(totfaces);
      ownerfaces = 0;

      int owner_ind = 1;
      int bc_ind = 1;

      // Loop through all the volume elements
      for(int elind = 1; elind <= ne; elind++)
      {
         Element el = mesh.VolumeElement(elind);

         Array<int> locfaces;

         locfaces.SetSize(el.GetNFaces());

         // Get the face numbers of the faces of the current volume element
         // The values returned are given a sign depending on the orientation 
         // of the faces. This is used while writing the faces file, to 
         // determine whether or not to invert the face triangle before writing 
         // it to file
         meshtopo.GetElementFaces(elind,locfaces,true);

         // Loop through the faces
         for(int i = 1; i <= locfaces.Size(); i++)
         {
            int absfacenr = abs(locfaces.Elem(i));

            // Check if the face is a surface element (boundary face)
            // if not, add the current volume element and the corresponding face into 
            // the owner list
            if(!(meshtopo.GetFace2SurfaceElement(absfacenr)))
            {
               // Location of the current face in the owner list 
               // (if it exists... zero otherwise)
               int ownerface = ownerfaces.Elem(absfacenr);

               // If the face already exists in the owner list, add 
               // the current cell into the neighbour list, in the 
               // same location where the face appears in the owner list
               if(ownerface)
               {
                  OF15x_neighbour_celllist.Elem(ownerface) = elind;
               }
               // If it is a new face which has not been listed before, 
               // add the current cell into the owner list, and save 
               // the index location to be used later by the neighbour list
               else
               {
                  OF15x_owner_celllist.Elem(owner_ind) = elind;
                  OF15x_owner_facelist.Elem(owner_ind) = locfaces.Elem(i);
                  // Set the bit array to indicate that the face is already processed
                  ownerfaces.Elem(absfacenr) = owner_ind;

                  owner_ind++;
               }
            }
            // If the face is a boundary face, extract the boundary condition number of the 
            // face, and append that along with the face number and the current cell 
            // into the various surface elements lists
            else
            {
               Element2d sel = mesh.SurfaceElement(meshtopo.GetFace2SurfaceElement(absfacenr));
               OF15x_surfelem_bclist.Elem(bc_ind) = mesh.GetFaceDescriptor(sel.GetIndex()).BCProperty();
               OF15x_surfelem_lists.Elem(bc_ind) = INDEX_2(locfaces.Elem(i),elind);

               bc_ind++;
            }
         }
      }

      // This correction is required in cases where the mesh has been "uniform refined".... for 
      // some reason, the number of faces reported by Netgen is higher than the actual number 
      // of faces in the mesh
      OF15x_owner_facelist.SetSize(owner_ind-1);
      OF15x_owner_celllist.SetSize(owner_ind-1);
      OF15x_neighbour_celllist.SetSize(owner_ind-1);


      // Sort the list of surface elements in ascending order of boundary condition number
      // also sort the cell list in the same manner
      QuickSort(OF15x_surfelem_bclist,OF15x_surfelem_lists);

      int rng_start = 1;
      int rng_end = 1;

      for(int i = 1; i <= OF15x_owner_celllist.Size(); i++)
      {     
         // Order the face list of each owner cell in accordance to the cell list 
         // of the neighbours of that owner sorted in ascending order
         // This operation is performed by selecting the right range within the 
         // array (basically... all the neighbours of each owner cell set are 
         // extracted and sorted one set at a time)
         if((OF15x_owner_celllist.Elem(i) == OF15x_owner_celllist.Elem(rng_start)) && (i != OF15x_owner_celllist.Size()))
         {
            rng_end = i;
         }
         else
         {
            if(i == OF15x_owner_celllist.Size()) rng_end = i;

            FlatArray<int> neisort_celllist = OF15x_neighbour_celllist.Range(rng_start-1,rng_end);
            FlatArray<int> ownersort_facelist = OF15x_owner_facelist.Range(rng_start-1,rng_end);

            // Here, QuickSort seems to fail for some reason, but since the size of 
            // each array cannot be larger than 4 elements, the performance penalty 
            // of using the BubbleSort should be negligable
            BubbleSort(neisort_celllist,ownersort_facelist);
            //QuickSort(neisort_celllist,ownersort_facelist);
            
            // After sorting out the cell and face lists, replace the old list with the 
            // newly ordered lists
            for(int j = 1; j <= neisort_celllist.Size(); j++)
            {
               OF15x_neighbour_celllist.Elem(rng_start-1+j) = neisort_celllist.Elem(j);
               OF15x_owner_facelist.Elem(rng_start-1+j) = ownersort_facelist.Elem(j);
            }

            // initialise the range variables to the next set of owner cells
            rng_start = i;
            rng_end = i;
         }
      }
/*
      ofstream dbg("OpenFOAMDebug.log");

      dbg << " ------- Boundary List -------- " << endl;

      for(int i = 1; i <= OF15x_surfelem_bclist.Size(); i++)
      {
         dbg << "bc = " << OF15x_surfelem_bclist.Elem(i) 
              << " : face = " << OF15x_surfelem_lists.Elem(i).I1()
              << " : cell = " << OF15x_surfelem_lists.Elem(i).I2() << endl;
      }

      dbg << endl << " ------- Owner List ------- " << endl;

      for(int i = 1; i <= OF15x_owner_celllist.Size(); i++)
      {
         dbg << "Ind:" << i << " :: (" 
              << OF15x_owner_celllist.Elem(i) << " "
              << OF15x_owner_facelist.Elem(i) << ")" << endl;
      }

      dbg << endl << " ----- Neighbour List ----- " << endl;

      for(int i = 1; i <= OF15x_neighbour_celllist.Size(); i++)
      {
         dbg << "Ind:" << i << " :: (" 
              << OF15x_neighbour_celllist.Elem(i) << ")" << endl;
      }

      dbg.close();
*/
      return(false);
   }



   void WriteOpenFOAM15xNeighbour (ofstream & outfile)
   {
      // Write the OpenFOAM standard banner and dividers, etc...
      WriteOpenFOAM15xBanner(outfile);
      outfile << "FoamFile \n"
              << "{ \n"
              << "    version     2.0; \n"
              << "    format      ascii; \n"
              << "    class       labelList; \n"
              << "    location    \"constant\\polyMesh\"; \n"
              << "    object      neighbour; \n"
              << "} \n";
      WriteOpenFOAM15xDividerStart(outfile);

      outfile << endl << endl;

      int nneighbours = OF15x_neighbour_celllist.Size();

      outfile << nneighbours << endl;

      outfile << "(" << endl;

      // Write the neighbour cells to file
      for(int i = 1; i <= OF15x_neighbour_celllist.Size(); i++)
      {
         outfile << OF15x_neighbour_celllist.Elem(i) - 1 << endl;
      }
      outfile << ")" << endl << endl;
      WriteOpenFOAM15xDividerEnd(outfile);
   }



   void WriteOpenFOAM15xOwner (ofstream & outfile)
   {
      // Write the OpenFOAM standard banner and dividers, etc...
      WriteOpenFOAM15xBanner(outfile);
      outfile << "FoamFile \n"
              << "{ \n"
              << "    version     2.0; \n"
              << "    format      ascii; \n"
              << "    class       labelList; \n"
              << "    location    \"constant\\polyMesh\"; \n"
              << "    object      owner; \n"
              << "} \n";
      WriteOpenFOAM15xDividerStart(outfile);

      outfile << endl << endl;

      int nowners = OF15x_owner_celllist.Size() + OF15x_surfelem_lists.Size();

      outfile << nowners << endl;

      outfile << "(" << endl;

      // Write the owners of the internal cells to file
      for(int i = 1; i <= OF15x_owner_celllist.Size(); i++)
      {
         outfile << OF15x_owner_celllist.Elem(i) - 1 << endl;
      }

      // Write the owners of the boundary cells to file
      // (Written in order of ascending boundary condition numbers)
      for(int i = 1; i <= OF15x_surfelem_lists.Size(); i++)
      {
         outfile << OF15x_surfelem_lists.Elem(i).I2() - 1 << endl;
      }
      outfile << ")" << endl << endl;
      WriteOpenFOAM15xDividerEnd(outfile);
   }



   void WriteOpenFOAM15xFaces (ofstream & outfile, const Mesh & mesh)
   {
      const_cast<Mesh&> (mesh).BuildElementSearchTree();
      const MeshTopology& meshtopo = mesh.GetTopology();
      
      // Update the mesh topology structures
      const_cast<Mesh&> (mesh).UpdateTopology();


      // Write the OpenFOAM standard banner and dividers, etc...
      WriteOpenFOAM15xBanner(outfile);
      outfile << "FoamFile \n"
              << "{ \n"
              << "    version     2.0; \n"
              << "    format      ascii; \n"
              << "    class       faceList; \n"
              << "    location    \"constant\\polyMesh\"; \n"
              << "    object      faces; \n"
              << "} \n";
      WriteOpenFOAM15xDividerStart(outfile);

      outfile << endl << endl;

      int nfaces = OF15x_owner_facelist.Size() + OF15x_surfelem_lists.Size();

      outfile << nfaces << endl;

      outfile << "(" << endl;

      // Write the faces in the order specified in the owners lists of the 
      // internal cells and the boundary cells
      for(int i = 1; i <= OF15x_owner_facelist.Size(); i++)
      {
         int face_w_orientation = OF15x_owner_facelist.Elem(i);
         int facenr = abs(face_w_orientation);

         Array<int> facepnts;

         meshtopo.GetFaceVertices(facenr,facepnts);

         // Get the orientation of the face, and invert it if required
         // for a quad, inversion => swap 1 <=> 2 and 3 <=> 4
         // for a trig, inversion => swap 1 <=> 3
         if(face_w_orientation > 0)
         {
            int tmppnts = 0;

            if(facepnts.Size() == 4)
            {
               tmppnts = facepnts.Elem(1);
               facepnts.Elem(1) = facepnts.Elem(2);
               facepnts.Elem(2) = tmppnts;
               
               tmppnts = facepnts.Elem(3);
               facepnts.Elem(3) = facepnts.Elem(4);
               facepnts.Elem(4) = tmppnts;
            }
            else if(facepnts.Size() == 3)
            {
               tmppnts = facepnts.Elem(1);
               facepnts.Elem(1) = facepnts.Elem(3);
               facepnts.Elem(3) = tmppnts;
            }
         }

         outfile << facepnts.Size();
         outfile << "(";
         for(int j = 1; j <= facepnts.Size(); j++)
         {
            outfile << facepnts.Elem(j)-1;
            if(j != facepnts.Size()) outfile << " ";
         }
         outfile << ")" << endl;
      }

      for(int i = 1; i <= OF15x_surfelem_lists.Size(); i++)
      {
         int face_w_orientation = OF15x_surfelem_lists.Elem(i).I1();
         int facenr = abs(face_w_orientation);

         Array<int> facepnts;

         meshtopo.GetFaceVertices(facenr,facepnts);

         // Get the orientation of the face, and invert it if required
         // for a quad, inversion => swap 1 <=> 2 and 3 <=> 4
         // for a trig, inversion => swap 1 <=> 3
         if(face_w_orientation > 0)
         {
            int tmppnts = 0;

            if(facepnts.Size() == 4)
            {
               tmppnts = facepnts.Elem(1);
               facepnts.Elem(1) = facepnts.Elem(2);
               facepnts.Elem(2) = tmppnts;
               
               tmppnts = facepnts.Elem(3);
               facepnts.Elem(3) = facepnts.Elem(4);
               facepnts.Elem(4) = tmppnts;
            }
            else if(facepnts.Size() == 3)
            {
               tmppnts = facepnts.Elem(1);
               facepnts.Elem(1) = facepnts.Elem(3);
               facepnts.Elem(3) = tmppnts;
            }
         }

         outfile << facepnts.Size();
         outfile << "(";
         for(int j = 1; j <= facepnts.Size(); j++)
         {
            outfile << facepnts.Elem(j)-1;
            if(j != facepnts.Size()) outfile << " ";
         }
         outfile << ")" << endl;
      }

      outfile << ")" << endl << endl;
      WriteOpenFOAM15xDividerEnd(outfile);
   }


 
   void WriteOpenFOAM15xPoints (ofstream & outfile, const Mesh & mesh)
   {
      int np = mesh.GetNP();

      // Write the OpenFOAM standard banner and dividers, etc...
      WriteOpenFOAM15xBanner(outfile);
      outfile << "FoamFile \n"
              << "{ \n"
              << "    version     2.0; \n"
              << "    format      ascii; \n"
              << "    class       vectorField; \n"
              << "    location    \"constant\\polyMesh\"; \n"
              << "    object      points; \n"
              << "} \n";
      WriteOpenFOAM15xDividerStart(outfile);

      outfile << endl << endl;

      // Number of points in the following list
      outfile << np << endl;

      outfile.precision(6);
      outfile.setf (ios::fixed, ios::floatfield);
      outfile.setf (ios::showpoint);

      // Coordinate list starts here
      outfile << "(" << endl;

      for(int i = 1; i <= np; i++)
      {
         const Point3d & p = mesh.Point(i);

         // Write coordinates to file
         outfile << "(";
         outfile << p.X() << " ";
         outfile << p.Y() << " ";
         outfile << p.Z();
         outfile << ")" << endl;
      }
      outfile << ")" << endl << endl;
      WriteOpenFOAM15xDividerEnd(outfile);
   }



   void WriteOpenFOAM15xBoundary (ofstream & outfile)
   {
      // Write the OpenFOAM standard banner and dividers, etc...
      WriteOpenFOAM15xBanner(outfile);
      outfile << "FoamFile \n"
              << "{ \n"
              << "    version     2.0; \n"
              << "    format      ascii; \n"
              << "    class       polyBoundaryMesh; \n"
              << "    location    \"constant\\polyMesh\"; \n"
              << "    object      boundary; \n"
              << "} \n";
      WriteOpenFOAM15xDividerStart(outfile);

      outfile << endl;


      Array<INDEX_3> bcarray;
      int ind = 1;

      bcarray.Append(INDEX_3(OF15x_surfelem_bclist.Elem(1),1,0));
            
      for(int i = 2; i <= OF15x_surfelem_bclist.Size(); i++)
      {
         if(OF15x_surfelem_bclist.Elem(i) == bcarray.Elem(ind).I1())
         {
            bcarray.Elem(ind).I2() = bcarray.Elem(ind).I2()+1;
         }
         else
         {
            ind++;
            bcarray.Append(INDEX_3(OF15x_surfelem_bclist.Elem(i),1,i-1));
         }
      }

      outfile << bcarray.Size() << endl;
      outfile << "(" << endl;

      int startface = 0;

      for(int i = 1; i <= bcarray.Size(); i++)
      {
         startface = OF15x_owner_celllist.Size() + bcarray.Elem(i).I3();

         outfile << "    patch" << bcarray.Elem(i).I1() << endl
                 << "    {" << endl
                 << "        type            patch;" << endl
                 << "        physicalType    patch;" << endl
                 << "        nFaces          " << bcarray.Elem(i).I2() << ";" << endl
                 << "        startFace       " << startface << ";" << endl
                 << "    }" << endl;
      }

      outfile << ")" << endl << endl;
      WriteOpenFOAM15xDividerEnd(outfile);
   }



   void WriteOpenFOAM15xFormat (const Mesh & mesh, const string & casename)
   {
      int i,j;
      bool error = false;
      char casefiles[256];

      // Make sure that the mesh data has been updated
      const_cast<Mesh&> (mesh).Compress();
      const_cast<Mesh&> (mesh).RebuildSurfaceElementLists();
      const_cast<Mesh&> (mesh).ComputeNVertices();


      int np = mesh.GetNP();
      int nse = mesh.GetNSE();
      int ne = mesh.GetNE();

      cout << "Write OpenFOAM 1.5+ Mesh Files...." << endl;

      // Abort if there are no points, surface elements or volume elements
      if((np <= 0) || (ne <= 0) || (nse <= 0))
      {
         cout << "Export Error: Invalid mesh.... Aborting!" << endl;
         return;
      }

      // OpenFOAM only supports linear meshes!
      if(mparam.secondorder || mesh.GetCurvedElements().IsHighOrder())
      {
         cout << "Export Error: OpenFOAM 1.5+ does not support non-linear elements.... Aborting!" << endl;
         return;
      }

      if((mesh.SurfaceElement(nse/2).GetType() != TRIG) 
         && (mesh.SurfaceElement(nse/2).GetType() != QUAD)
         || (mesh.VolumeElement(ne/2).GetType() == TET10)
         || (mesh.VolumeElement(ne/2).GetType() == PRISM12))
      {
         cout << "Export Error: OpenFOAM 1.5+ does not support non-linear elements.... Aborting!" << endl;
         return;
      }


      cout << "Writing OpenFOAM 1.5+ Mesh files to case: " << casename << endl;

      // Create the case directory if it does not already exist
      // NOTE: This needs to be improved for the Linux variant....!!!
   #ifdef WIN32
      char casedir[256];
      sprintf(casedir, "mkdir %s\\constant\\polyMesh", casename.c_str());
      system(casedir);
   #else
      char casedir[256];
      mkdir(casename.c_str(), S_IRWXU|S_IRWXG);
      sprintf(casedir, "%s/constant", casename.c_str());
      mkdir(casedir, S_IRWXU|S_IRWXG);
      sprintf(casedir, "%s/constant/polyMesh", casename.c_str());
      mkdir(casedir, S_IRWXU|S_IRWXG);
   #endif

      // Open handles to the five required mesh files
      // points
      // faces
      // owner
      // neighbour
      // boundary
      sprintf(casefiles, "%s/constant/polyMesh/points", casename.c_str());
      ofstream outfile_pnts(casefiles);
      sprintf(casefiles, "%s/constant/polyMesh/faces", casename.c_str()); 
      ofstream outfile_faces(casefiles);
      sprintf(casefiles, "%s/constant/polyMesh/owner", casename.c_str()); 
      ofstream outfile_own(casefiles);
      sprintf(casefiles, "%s/constant/polyMesh/neighbour", casename.c_str()); 
      ofstream outfile_nei(casefiles);
      sprintf(casefiles, "%s/constant/polyMesh/boundary", casename.c_str()); 
      ofstream outfile_bnd(casefiles);

      ResetTime();

      // Build the owner, neighbour, faces and boundary lists 
      // from the Netgen mesh
      cout << endl << "Building Owner, Neighbour and Face Lists: ";

      error = BuildOpenFOAM15xLists(mesh);

      cout << "Done! (Time Elapsed = " << GetTime() << " sec)" << endl;


      // Write the "owner" file
      if(outfile_own.good() && !error)
      {
         cout << "Writing the owner file: ";
         WriteOpenFOAM15xOwner(outfile_own);
         outfile_own.close();
         cout << "Done! (Time Elapsed = " << GetTime() << " sec)" << endl;
      }
      else
      {
         cout << "Export Error: Error creating file: owner.... Aborting" << endl;
         error = true;
      }


      // Write the "neighbour" file
      if(outfile_nei.good() && !error)
      {
         cout << "Writing the neighbour file: ";
         WriteOpenFOAM15xNeighbour(outfile_nei);
         outfile_nei.close();
         cout << "Done! (Time Elapsed = " << GetTime() << " sec)" << endl;
      }
      else
      {
         cout << "Export Error: Error creating file: neighbour.... Aborting" << endl;
         error = true;
      }


      // Write the "faces" file
      if(outfile_faces.good() && !error)
      {
         cout << "Writing the faces file: ";
         WriteOpenFOAM15xFaces(outfile_faces, mesh);
         outfile_faces.close();
         cout << "Done! (Time Elapsed = " << GetTime() << " sec)" << endl;
      }
      else
      {
         cout << "Export Error: Error creating file: faces.... Aborting" << endl;
         error = true;
      }


      // Write the "points" file
      if(outfile_pnts.good() && !error)
      {
         cout << "Writing the points file: ";
         WriteOpenFOAM15xPoints(outfile_pnts,mesh);
         outfile_pnts.close();
         cout << "Done! (Time Elapsed = " << GetTime() << " sec)" << endl;
      }
      else
      {
         cout << "Export Error: Error creating file: points.... Aborting" << endl;
         error = true;
      }


      // Write the "boundary" file
      if(outfile_bnd.good() && !error)
      {
         cout << "Writing the boundary file: ";
         WriteOpenFOAM15xBoundary(outfile_bnd);
         outfile_bnd.close();
         cout << "Done! (Time Elapsed = " << GetTime() << " sec)" << endl;
      }
      else
      {
         cout << "Export Error: Error creating file: boundary.... Aborting" << endl;
         error = true;
      }

      if(!error)
      {
         cout << "OpenFOAM 1.5+ Export successfully completed (Time elapsed = " << GetTime() << " sec) !" << endl;
      }
      else
      {
         cout << "Error in OpenFOAM 1.5+ Export.... Aborted!" << endl;
      }
   }
}

