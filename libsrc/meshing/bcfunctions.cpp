
#include <mystdlib.h>
#include <meshing.hpp>


namespace netgen
{
   // Default colour to be used for boundary condition number "0"
   #define DEFAULT_R       0.0
   #define DEFAULT_G       1.0
   #define DEFAULT_B       0.0

   // Boundary condition number to use if a face does not have a 
   // colour assigned to it, or if the colour is the above defined 
   // default colour
   #define DEFAULT_BCNUM   1

   // Default tolerance for colour matching (using Euclidean distance)
   #define DEFAULT_EPS     2.5e-05




   /*! Philippose - 11/07/2009
       Function to check if two RGB colours are equal

       Note#1: Currently uses unweighted Euclidean Distance 
       for colour matching.

       Note#2: The tolerance used for deciding whether two 
       colours match is defined as "eps" and is currently 
       2.5e-5 (for square of distance)
   */
   bool ColourMatch(Vec3d col1, Vec3d col2, double eps = DEFAULT_EPS)
   {
      if(eps <= 0.0) eps = DEFAULT_EPS;
      
      bool colmatch = false;

      if(Dist2(col1,col2) < eps) colmatch = true;

      return colmatch;
   }
      




   /*! Philippose - 11/07/2009
       Function to create a list of all the unique colours 
       available in a given mesh
   */
   void GetFaceColours(Mesh & mesh, Array<Vec3d> & face_colours)
   {
      face_colours.SetSize(1);
      face_colours.Elem(1) = mesh.GetFaceDescriptor(1).SurfColour();
      
      for(int i = 1; i <= mesh.GetNFD(); i++)
      {
         Vec3d face_colour = mesh.GetFaceDescriptor(i).SurfColour();
         bool col_found = false;
         
         for(int j = 1; j <= face_colours.Size(); j++)
         {
            if(ColourMatch(face_colours.Elem(j),face_colour))
            {
               col_found = true;
               break;
            }
         }
         
         if(!col_found) face_colours.Append(face_colour);
      }

      if(printmessage_importance >= 3)
      {
         cout << endl << "-------- Face Colours --------" << endl;
         for( int i = 1; i <= face_colours.Size(); i++)
         {
            cout << face_colours.Elem(i) << endl;
         }
         cout << "------------------------------" << endl;
      }
   }






   /*! Philippose - 11/07/2009
       Assign boundary condition numbers based on a user defined 
       colour profile file.

       The default profile file is "netgen.ocf"

       If the mesh contains colours not defined in the profile,
       netgen automatically starts assigning each new colour a 
       new boundary condition number starting from the highest 
       boundary condition number specified in the profile file.
   */
   void AutoColourAlg_UserProfile(Mesh & mesh, ifstream & ocf)
   {
      char ocf_inp[100];
      bool header_found = false;

      // Number of colour specifications in the 
      // user profile file
      int numentries = 0;
      while((ocf.good()) && (!header_found))
      {
         ocf >> ocf_inp;
         if(strcmp(ocf_inp,"boundary_colours") == 0) header_found = true;
      }

      if(!header_found)
      {
         ocf.close();
         throw NgException("AutoColourAlg_UserProfile: Invalid or empty Boundary Colour Profile file\n");
         return;
      }

      // Read in the number of entries from file
      ocf >> numentries;
      if(numentries > 0)
      {
         if(!ocf.good())
         {
            ocf.close();
            throw NgException("AutoColourAlg_UserProfile: Invalid or empty Boundary Colour Profile file\n");
            return;
         }

         PrintMessage(3, "Number of colour entries: ", numentries);
      }
      else
      {
         ocf.close();
         PrintMessage(3, "AutoColourAlg_UserProfile: No Boundary Colour entries found.... no changes made!");
         return;
      }

      // Arrays to hold the specified RGB colour triplets as well 
      // as the associated boundary condition number
      Array<Vec3d> bc_colours(numentries);
      Array<int> bc_num(numentries);
      Array<bool> bc_used(numentries);
      
      // Actually read in the data from the file
      for(int i = 1; i <= numentries; i++)
      {
         int bcnum;
         // double col_red, col_green, col_blue;

         ocf >> bcnum;
         // Boundary condition number DEFAULT_BCNUM is reserved for 
         // faces which have the default colour Green (0.0,1.0,0.0)
         // To prevent confusion, no boundary numbery below this default 
         // are permitted
         if(bcnum < (DEFAULT_BCNUM + 1)) bcnum = DEFAULT_BCNUM+1;

         bc_num.Elem(i) = bcnum;
         bc_used.Elem(i) = false;
         ocf >> bc_colours.Elem(i).X() 
             >> bc_colours.Elem(i).Y() 
             >> bc_colours.Elem(i).Z();

         if(!ocf.good())
         {
            ocf.close();
            throw NgException("Boundary Colour file error: Number of entries do not match specified list size!!\n");
            return;
         }

         // Bound checking of the values
         // The RGB values should be between 0.0 and 1.0
         if(bc_colours.Elem(bcnum).X() < 0.0) bc_colours.Elem(bcnum).X() = 0.0;
         if(bc_colours.Elem(bcnum).X() > 1.0) bc_colours.Elem(bcnum).X() = 1.0;
         if(bc_colours.Elem(bcnum).Y() < 0.0) bc_colours.Elem(bcnum).X() = 0.0;
         if(bc_colours.Elem(bcnum).Y() > 1.0) bc_colours.Elem(bcnum).X() = 1.0;
         if(bc_colours.Elem(bcnum).Z() < 0.0) bc_colours.Elem(bcnum).X() = 0.0;
         if(bc_colours.Elem(bcnum).Z() > 1.0) bc_colours.Elem(bcnum).X() = 1.0;
      }

      PrintMessage(3, "Successfully loaded Boundary Colour Profile file....");
      ocf.close();

      // Find the highest boundary condition number in the list
      // All colours in the geometry which are not specified in the 
      // list will be given boundary condition numbers higher than this 
      // number
      int max_bcnum = DEFAULT_BCNUM;
      for(int i = 1; i <= bc_num.Size();i++)
      {
         if(bc_num.Elem(i) > max_bcnum) max_bcnum = bc_num.Elem(i);
      }

      PrintMessage(3, "Highest boundary number in list = ",max_bcnum);

      Array<Vec3d> all_colours;
      
      // Extract all the colours to see how many there are
      GetFaceColours(mesh,all_colours);
      PrintMessage(3,"\nNumber of colours defined in Mesh: ", all_colours.Size());

      if(all_colours.Size() == 0)
      {
         PrintMessage(3,"No colour data detected in Mesh... no changes made!");
         return;
      }

      int nfd = mesh.GetNFD();

      for(int face_index = 1; face_index <= nfd; face_index++)
      {
         // Temporary container for individual face colours
         Vec3d face_colour;

         // Get the colour of the face being currently processed
         face_colour = mesh.GetFaceDescriptor(face_index).SurfColour();
         if(!ColourMatch(face_colour,Vec3d(DEFAULT_R,DEFAULT_G,DEFAULT_B)))
         {
            // Boolean variable to check if the boundary condition was applied 
            // or not... not applied would imply that the colour of the face 
            // does not exist in the list of colours in the profile file
            bool bc_assigned = false;

            for(int col_index = 1; col_index <= bc_colours.Size(); col_index++)
            {
               if((ColourMatch(face_colour,bc_colours.Elem(col_index))) && (!bc_assigned))
               {
                  mesh.GetFaceDescriptor(face_index).SetBCProperty(bc_num.Elem(col_index));
                  bc_used.Elem(col_index) = true;
                  bc_assigned = true;
                  break;
               }
            }

            // If the colour was not found in the list, add it to the list, and assign 
            // the next free boundary condition number to it
            if(!bc_assigned)
            {
               max_bcnum++;
               bc_num.Append(max_bcnum);
               bc_colours.Append(face_colour);
               bc_used.Append(true);

               mesh.GetFaceDescriptor(face_index).SetBCProperty(max_bcnum);
            }
         }
         else
         {
            // Set the boundary condition number to the default one
            mesh.GetFaceDescriptor(face_index).SetBCProperty(DEFAULT_BCNUM);
         }
      }

      // User Information of the results of the operation
      Vec3d ref_colour(0.0,1.0,0.0);
      PrintMessage(3,"Colour based Boundary Condition Property details:");
      for(int bc_index = 0; bc_index <= bc_num.Size(); bc_index++)
      {
         if(bc_index > 0) ref_colour = bc_colours.Elem(bc_index);

         if(bc_index == 0) 
         {
            PrintMessage(3, "BC Property: ",DEFAULT_BCNUM);
            PrintMessage(3, "   RGB Face Colour = ",ref_colour,"","\n");
         }
         else if(bc_used.Elem(bc_index))
         {
            PrintMessage(3, "BC Property: ",bc_num.Elem(bc_index));
            PrintMessage(3, "   RGB Face Colour = ",ref_colour,"","\n");
         }
      }
   }




   
   /*! Philippose - 11/07/2009
       Assign boundary condition numbers based on the colours 
       assigned to each face in the mesh using an automated 
       algorithm.

       The particular algorithm used has been briefly explained 
       in the header file "occauxfunctions.hpp"
   */
   void AutoColourAlg_Sorted(Mesh & mesh)
   {
      Array<Vec3d> all_colours;
      Array<int> faces_sorted;
      Array<int> colours_sorted;

      // Extract all the colours to see how many there are
      GetFaceColours(mesh,all_colours);

      // Delete the default colour from the list since it will be accounted 
      // for automatically
      for(int i = 1; i <= all_colours.Size(); i++)
      {
         if(ColourMatch(all_colours.Elem(i),Vec3d(DEFAULT_R,DEFAULT_G,DEFAULT_B)))
         {
            all_colours.DeleteElement(i);
            break;
         }
      }
      PrintMessage(3,"\nNumber of colours defined in Mesh: ", all_colours.Size());

      if(all_colours.Size() == 0)
      {
         PrintMessage(3,"No colour data detected in Mesh... no changes made!");
         return;
      }

      // One more slot than the number of colours are required, to 
      // account for individual faces which have no colour data 
      // assigned to them in the CAD software
      faces_sorted.SetSize(all_colours.Size()+1);
      colours_sorted.SetSize(all_colours.Size()+1);
      faces_sorted = 0;
      
      // Slave Array to identify the colours the faces were assigned to, 
      // after the bubble sort routine to sort the automatic boundary 
      // identifiers according to the number of surface mesh elements 
      // of a given colour
      for(int i = 0; i <= all_colours.Size(); i++) colours_sorted[i] = i;

      // Used to hold the number of surface elements without any OCC 
      // colour definition
      int no_colour_faces = 0;

      // Index in the faces array assigned to faces without any 
      // or the default colour definition
      int no_colour_index = 0;

      int nfd = mesh.GetNFD();

      // Extract the number of surface elements having a given colour
      // And save this number into an array for later sorting
      for(int face_index = 1; face_index <= nfd; face_index++)
      {
         Array<SurfaceElementIndex> se_face;

         mesh.GetSurfaceElementsOfFace(face_index, se_face);

         Vec3d face_colour;

         face_colour = mesh.GetFaceDescriptor(face_index).SurfColour();
         if(!ColourMatch(face_colour,Vec3d(DEFAULT_R,DEFAULT_G,DEFAULT_B)))
         {
            for(int i = 1; i <= all_colours.Size(); i++)
            {
               if(ColourMatch(face_colour, all_colours.Elem(i)))
               {
                  faces_sorted[i] = faces_sorted[i] + se_face.Size();
               }
            }
         }
         else
         {
            // Add the number of surface elements without any colour 
            // definition separately
            no_colour_faces = no_colour_faces + se_face.Size();
         }
      }

      // Sort the face colour indices according to the number of surface 
      // mesh elements which have a specific colour
      BubbleSort(faces_sorted,colours_sorted);

      // Now update the array position assigned for surface elements 
      // without any colour definition with the number of elements
      faces_sorted[no_colour_index] = no_colour_faces;

      // Now actually assign the BC Property to the respective faces
      for(int face_index = 1; face_index <= nfd; face_index++)
      {
         Vec3d face_colour;

         face_colour = mesh.GetFaceDescriptor(face_index).SurfColour();
         if(!ColourMatch(face_colour,Vec3d(DEFAULT_R,DEFAULT_G,DEFAULT_B)))
         {
            for(int i = 0; i < colours_sorted.Size(); i++)
            {
               Vec3d ref_colour;
               if(i != no_colour_index) ref_colour = all_colours.Elem(colours_sorted[i]);

               if(ColourMatch(face_colour, ref_colour))
               {
                  mesh.GetFaceDescriptor(face_index).SetBCProperty(i + DEFAULT_BCNUM);
               }
            }
         }
         else
         {
            mesh.GetFaceDescriptor(face_index).SetBCProperty(DEFAULT_BCNUM);
         }

         PrintMessage(4,"Face number: ",face_index," ; BC Property = ",mesh.GetFaceDescriptor(face_index).BCProperty());
      }

      // User Information of the results of the operation
      Vec3d ref_colour(0.0,1.0,0.0);
      PrintMessage(3,"Colour based Boundary Condition Property details:");
      for(int i = 0; i < faces_sorted.Size(); i++)
      {
         if(colours_sorted[i] > 0) ref_colour = all_colours.Elem(colours_sorted[i]);

         PrintMessage(3, "BC Property: ",i + DEFAULT_BCNUM);
         PrintMessage(3, "   Nr. of Surface Elements = ", faces_sorted[i]);
         PrintMessage(3, "   Colour Index = ", colours_sorted[i]);
         PrintMessage(3, "   RGB Face Colour = ",ref_colour,"","\n");
      }
   }





   /*! Philippose - 13/07/2009
       Main function implementing automated assignment of 
       Boundary Condition numbers based on face colours

       This functionality is currently implemtented at the mesh 
       level, and hence allows colour based assignment of boundary 
       conditions for any geometry type within netgen which 
       supports face colours
   */
   void AutoColourBcProps(Mesh & mesh, const char * bccolourfile)
   {
      // Go directly to the alternate algorithm if no colour profile file was specified
      if(!bccolourfile)
      {
         PrintMessage(1,"AutoColourBcProps: Using Automatic Colour based boundary property assignment algorithm");
         AutoColourAlg_Sorted(mesh);
      }
      else
      {
         ifstream ocf(bccolourfile);

         // If there was an error opening the Colour profile file, jump to the alternate 
         // algorithm after printing a message
         if(!ocf)
         {
            PrintMessage(1,"AutoColourBcProps: Error loading Boundary Colour Profile file ", 
                         bccolourfile, " ....","Switching to Automatic Assignment algorithm!");

            AutoColourAlg_Sorted(mesh);
         }
         // If the file opens successfully, call the function which assigns boundary conditions 
         // based on the colour profile file
         else
         {
            PrintMessage(1, "AutoColourBcProps: Using Boundary Colour Profile file: ");
            PrintMessage(1, "  ", bccolourfile);
            AutoColourAlg_UserProfile(mesh, ocf);

            // Make sure the file is closed before exiting the function
            if(ocf.is_open())
            {
               ocf.close();
            }
         }
      }
   }
}
