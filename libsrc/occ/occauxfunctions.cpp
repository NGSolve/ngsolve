#ifdef OCCGEOMETRY

#include <mystdlib.h>
#include <meshing.hpp>

#include <occgeom.hpp>

namespace netgen
{
   void OCCAutoColourAlg_UserProfile(Mesh & mesh, OCCGeometry & occgeometry, ifstream & ocf)
   {
      int numentries = 0;

      ocf >> numentries;
      if(numentries > 0)
      {
         if(!ocf.good())
         {
            throw NgException("OCCAutoColourAlg_UserProfile: Invalid or empty Boundary Colour Profile file\n");
            return;
         }

         PrintMessage(3, "Number of colour entries: ", numentries);
      }
      else
      {
         PrintMessage(3, "OCCAutoColourAlg_UserProfile: No Boundary Colour entries found.... no changes made!");
         ocf.close();
         return;
      }

      Array<Vec3d> bc_colours(numentries);
      Array<int> bc_num(numentries);
      
      for(int i = 1; i <= numentries; i++)
      {
         int bcnum;
         double col_red, col_green, col_blue;

         ocf >> bcnum;
         if(bcnum < 1) bcnum = 1;

         bc_num.Elem(i) = bcnum;
         ocf >> bc_colours.Elem(i).X() 
             >> bc_colours.Elem(i).Y() 
             >> bc_colours.Elem(i).Z();

         if(!ocf.good())
            throw NgException("Boundary Colour file error: Number of entries do not match specified list size!!\n");

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
      int max_bcnum = 0;
      for(int i = 1; i <= bc_num.Size();i++)
      {
         if(bc_num.Elem(i) > max_bcnum) max_bcnum = bc_num.Elem(i);
      }

      PrintMessage(3, "Highest boundary number in list = ",max_bcnum);

      TDF_LabelSequence all_colours;
      
      // Extract all the colours to see how many there are
      occgeometry.face_colours->GetColors(all_colours);
      PrintMessage(3,"\nNumber of colours defined in OCC Mesh: ", all_colours.Length());

      if(all_colours.Length() == 0)
      {
         PrintMessage(3,"No colour data detected in OCC Mesh... no changes made!");
         return;
      }

      int nfd = mesh.GetNFD();

      for(int face_index = 1; face_index <= nfd; face_index++)
      {
         // Note: From the logic in file "occgenmesh.cpp", function "FindEdges" 
         // the Face Descriptor Number of an OCC face has a one-to-one mapping  
         // to the face number in the OCC face map (fmap)
         TopoDS_Face face = TopoDS::Face(occgeometry.fmap(face_index));
         Quantity_Color face_colour;

         if(occgeometry.face_colours->GetColor(face,XCAFDoc_ColorSurf,face_colour))
         {
            // Boolean variable to check if the boundary condition was applied 
            // or not... not applied would imply that the colour of the face 
            // does not exist in the list of colours in the profile file
            bool bc_assigned = false;

            for(int col_index = 1; col_index <= bc_colours.Size(); col_index++)
            {
               Quantity_Color bc_colour;

               double col_red = bc_colours.Elem(col_index).X();
               double col_green = bc_colours.Elem(col_index).Y();
               double col_blue = bc_colours.Elem(col_index).Z();

               bc_colour.SetValues(col_red,col_green,col_blue,Quantity_TOC_RGB);
               
               if((face_colour == bc_colour) && (!bc_assigned))
               {
                  mesh.GetFaceDescriptor(face_index).SetBCProperty(bc_num.Elem(col_index));
                  bc_assigned = true;
                  break;
               }
            }

            // If the colour was not found in the list, add it to the list, and assign 
            // the next free boundary condition number to it
            if(!bc_assigned)
            {
               double col_red = face_colour.Red();
               double col_green = face_colour.Green();
               double col_blue = face_colour.Blue();

               Vec3d new_colour(col_red,col_green,col_blue);

               max_bcnum++;
               bc_num.Append(max_bcnum);
               bc_colours.Append(new_colour);

               mesh.GetFaceDescriptor(face_index).SetBCProperty(max_bcnum);
            }
         }
         else
         {
            mesh.GetFaceDescriptor(face_index).SetBCProperty(0);
         }
      }
   }



   void OCCAutoColourAlg_Sorted(Mesh & mesh, OCCGeometry & occgeometry)
   {
      TDF_LabelSequence all_colours;
      Array<int> faces_sorted;
      Array<int> colours_sorted;

      // Extract all the colours to see how many there are
      occgeometry.face_colours->GetColors(all_colours);
      PrintMessage(3,"\nNumber of colours defined in OCC Mesh: ", all_colours.Length());

      if(all_colours.Length() == 0)
      {
         PrintMessage(3,"No colour data detected in OCC Mesh... no changes made!");
         return;
      }

      // One more slot than the number of colours are required, to 
      // account for individual faces which have no colour data 
      // assigned to them in the CAD software
      faces_sorted.SetSize(all_colours.Length()+1); 
      colours_sorted.SetSize(all_colours.Length()+1);
      faces_sorted = 0;
      
      // Slave Array to identify the colours the faces were assigned to, 
      // after the bubble sort routine to sort the automatic boundary 
      // identifiers according to the number of surface mesh elements 
      // of a given colour
      for(int i = 0; i <= all_colours.Length(); i++) colours_sorted[i] = i;

      // Used to hold the number of surface elements without any OCC 
      // colour definition
      int no_colour_faces = 0;

      // Index in the faces array assigned to faces without an 
      // OCC Colour definition
      int no_colour_index = 0;

      int nfd = mesh.GetNFD();

      // Extract the number of surface elements having a given OCC colour
      // And save this number into an array for later sorting
      for(int face_index = 1; face_index <= nfd; face_index++)
      {
         Array<SurfaceElementIndex> se_face;

         mesh.GetSurfaceElementsOfFace(face_index, se_face);

         // Note: From the logic in file "occgenmesh.cpp", function "FindEdges" 
         // the Face Descriptor Number of an OCC face has a one-to-one mapping  
         // to the face number in the OCC face map (fmap)
         TopoDS_Face face = TopoDS::Face(occgeometry.fmap(face_index));
         Quantity_Color face_colour;

         if(occgeometry.face_colours->GetColor(face,XCAFDoc_ColorSurf,face_colour))
         {
            for(int i = 1; i <= all_colours.Length(); i++)
            {
               Quantity_Color ref_colour;
               occgeometry.face_colours->GetColor(all_colours.Value(i),ref_colour);
               if(face_colour == ref_colour)
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
         TopoDS_Face face = TopoDS::Face(occgeometry.fmap(face_index));
         Quantity_Color face_colour;

         if(occgeometry.face_colours->GetColor(face,XCAFDoc_ColorSurf,face_colour))
         {
            for(int i = 0; i < colours_sorted.Size(); i++)
            {
               Quantity_Color ref_colour;
               if(i != no_colour_index)
                  occgeometry.face_colours->GetColor(all_colours.Value(colours_sorted[i]),ref_colour);

               if(face_colour == ref_colour)
               {
                  mesh.GetFaceDescriptor(face_index).SetBCProperty(i);
               }
            }
         }
         else
         {
            mesh.GetFaceDescriptor(face_index).SetBCProperty(no_colour_index);
         }

         PrintMessage(4,"Face number: ",face_index," ; BC Property = ",mesh.GetFaceDescriptor(face_index).BCProperty());
      }

      // User Information of the results of the operation
      PrintMessage(3,"OCC Colour based Boundary Condition Property details:");
      for(int i = 0; i < faces_sorted.Size(); i++)
      {
         Quantity_Color ref_colour;

         ref_colour.SetValues(1.0,1.0,0.0,Quantity_TOC_RGB);

         if(colours_sorted[i] > 0) occgeometry.face_colours->GetColor(all_colours.Value(colours_sorted[i]),ref_colour);

         PrintMessage(3, "BC Property: ",i);
         PrintMessage(3, "   Nr. of Surface Elements = ", faces_sorted[i]);
         PrintMessage(3, "   OCC Colour Index = ", colours_sorted[i]);
         PrintMessage(3, "   RGB Face Colour = (",ref_colour.Red(),","
                                                 ,ref_colour.Green(),","
                                                 ,ref_colour.Blue(),")","\n");
      }
   }



   void OCCAutoColourBcProps(Mesh & mesh, OCCGeometry & occgeometry, const char * occcolourfile)
   {
      // Go directly to the alternate algorithm if no colour profile file was specified
      if(!occcolourfile)
      {
         OCCAutoColourAlg_Sorted(mesh,occgeometry);
      }
      
      ifstream ocf(occcolourfile);

      // If there was an error opening the Colour profile file, jump to the alternate 
      // algorithm after printing a message
      if(!ocf)
      {
         PrintMessage(1,"OCCAutoColourBcProps: Error loading Boundary Colour Profile file ", 
                      occcolourfile, " ....","Switching to alternate algorithm!");

         OCCAutoColourAlg_Sorted(mesh,occgeometry);
      }
      // If the file opens successfully, call the function which assigns boundary conditions 
      // based on the colour profile file
      else
      {
         OCCAutoColourAlg_UserProfile(mesh,occgeometry,ocf);
      }
   }
}
#endif // #ifdef OCCGEOMETRY
