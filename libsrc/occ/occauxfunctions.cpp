#ifdef OCCGEOMETRY

#include <mystdlib.h>
#include <meshing.hpp>

#include <occgeom.hpp>

namespace netgen
{
   void OCCAutoColourBcProps(Mesh & mesh, OCCGeometry & occgeometry)
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
}
#endif // #ifdef OCCGEOMETRY
