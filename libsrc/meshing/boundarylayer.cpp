#include <mystdlib.h>
#include "meshing.hpp"

namespace netgen
{

   void InsertVirtualBoundaryLayer (Mesh & mesh)
   {
      cout << "Insert virt. b.l." << endl;

      int surfid;

      cout << "Boundary Nr:";
      cin >> surfid;

      int i, j;
      int np = mesh.GetNP();

      cout << "Old NP: " << mesh.GetNP() << endl;
      cout << "Trigs: " << mesh.GetNSE() << endl;

      BitArray bndnodes(np);
      Array<int> mapto(np);

      bndnodes.Clear();
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.LineSegment(i).edgenr;
         cout << "snr = " << snr << endl;
         if (snr == surfid)
         {
            bndnodes.Set (mesh.LineSegment(i)[0]);
            bndnodes.Set (mesh.LineSegment(i)[1]);
         }
      }
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.LineSegment(i).edgenr;
         if (snr != surfid)
         {
            bndnodes.Clear (mesh.LineSegment(i)[0]);
            bndnodes.Clear (mesh.LineSegment(i)[1]);
         }
      }

      for (i = 1; i <= np; i++)
      {
         if (bndnodes.Test(i))
            mapto.Elem(i) = mesh.AddPoint (mesh.Point (i));
         else
            mapto.Elem(i) = 0;
      }

      for (i = 1; i <= mesh.GetNSE(); i++)
      {
         Element2d & el = mesh.SurfaceElement(i);
         for (j = 1; j <= el.GetNP(); j++)
            if (mapto.Get(el.PNum(j)))
               el.PNum(j) = mapto.Get(el.PNum(j));
      }


      int nq = 0;
      for (i = 1; i <= mesh.GetNSeg(); i++)
      {
         int snr = mesh.LineSegment(i).edgenr;
         if (snr == surfid)
         {
            int p1 = mesh.LineSegment(i)[0];
            int p2 = mesh.LineSegment(i)[1];
            int p3 = mapto.Get (p1);
            if (!p3) p3 = p1;
            int p4 = mapto.Get (p2);
            if (!p4) p4 = p2;

            Element2d el(QUAD);
            el.PNum(1) = p1;
            el.PNum(2) = p2;
            el.PNum(3) = p3;
            el.PNum(4) = p4;
            el.SetIndex (2);
            mesh.AddSurfaceElement (el);
            nq++;
         }
      }

      cout << "New NP: " << mesh.GetNP() << endl;
      cout << "Quads: " << nq << endl;
   }





/*
   Philippose Rajan - 11 June 2009

   Function to calculate the surface normal at a given 
   vertex of a surface element, with respect to that 
   surface element.

   This function is used by the boundary layer generation 
   function, in order to calculate the effective direction 
   in which the prismatic layer should grow
*/
   void GetSurfaceNormal(Mesh & mesh, Element2d & el, int Vertex, Vec3d & SurfaceNormal)
   {
      int Vertex_A;
      int Vertex_B;

      Vertex_A = Vertex + 1;
      if(Vertex_A > el.GetNP()) Vertex_A = 1;

      Vertex_B = Vertex - 1;
      if(Vertex_B <= 0) Vertex_B = el.GetNP();

      Vec3d Vect_A,Vect_B;
      
      Vect_A = mesh.Point(el.PNum(Vertex_A)) - mesh.Point(el.PNum(Vertex));
      Vect_B = mesh.Point(el.PNum(Vertex_B)) - mesh.Point(el.PNum(Vertex));

      SurfaceNormal = Cross(Vect_A,Vect_B);
      SurfaceNormal.Normalize();
   }





/*
    Philippose Rajan - 11 June 2009
    
    Added an initial experimental function for 
    generating prismatic boundary layers on 
    a given set of surfaces.
    
    The number of layers, height of the first layer 
    and the growth / shrink factor can be specified 
    by the user

    Currently, the layer height is calculated using:
    height = h_first_layer * (growth_factor^(num_layers - 1))
*/
   void GenerateBoundaryLayer (Mesh & mesh, MeshingParameters & mp)
   {
      int i, j;

      ofstream dbg("BndLayerDebug.log");

      // Angle between a surface element and a growth-vector below which 
      // a prism is project onto that surface as a quad
      // (in degrees)
      double angleThreshold = 5.0;
      
      cout << "Generate Prismatic Boundary Layers (Experimental)...." << endl;

      // Use an array to support creation of boundary 
      // layers for multiple surfaces in the future...
      Array<int> surfid;
      int surfinp = 0;
      int prismlayers = 1;
      double hfirst = 0.01;
      double growthfactor = 1.0;

      // Monitor and print out the number of prism and quad elements 
      // added to the mesh
      int numprisms = 0;
      int numquads = 0;

      while(surfinp >= 0)
      {
         cout << "Enter Surface ID (-1 to end list): ";
         cin >> surfinp;
         if(surfinp >= 0) surfid.Append(surfinp);
      }

      cout << "Number of surfaces entered = " << surfid.Size() << endl; 
      cout << "Selected surfaces are:" << endl;

      for(i = 1; i <= surfid.Size(); i++)
      {
         cout << "Surface " << i << ": " << surfid.Elem(i) << endl;
      }
      
      cout << endl << "Enter number of prism layers: ";
      cin >> prismlayers;
      if(prismlayers < 1) prismlayers = 1;

      cout << "Enter height of first layer: ";
      cin >> hfirst;
      if(hfirst <= 0.0) hfirst = 0.01;

      cout << "Enter layer growth / shrink factor: ";
      cin >> growthfactor;
      if(growthfactor <= 0.0) growthfactor = 0.5;

      cout << "Old NP: " << mesh.GetNP() << endl;
      cout << "Old NSE: " << mesh.GetNSE() << endl;
      
      for(int layer = prismlayers; layer >= 1; layer--)
      {
         cout << "Generating layer: " << layer << endl;

         const MeshTopology& meshtopo = mesh.GetTopology();
         const_cast<MeshTopology &> (meshtopo).SetBuildEdges(true);
         const_cast<MeshTopology &> (meshtopo).SetBuildFaces(true);
         const_cast<MeshTopology &> (meshtopo).Update();

         double layerht = hfirst;

         if(growthfactor == 1)
         {
            layerht = layer * hfirst;
         }
         else
         {
            layerht = hfirst*(pow(growthfactor,(layer+1)) - 1)/(growthfactor - 1);
         }

         cout << "Layer Height = " << layerht << endl;

         // Need to store the old number of points and 
         // surface elements because there are new points and 
         // surface elements being added during the process
         int np = mesh.GetNP();
         int nse = mesh.GetNSE();

         // Safety measure to ensure no issues with mesh 
         // consistency
         int nseg = mesh.GetNSeg();

         // Indicate which points need to be remapped
         BitArray bndnodes(np);

         // Map of the old points to the new points
         Array<int> mapto(np);

         // Growth vectors for the prismatic layer based on 
         // the effective surface normal at a given point
         Array<Vec3d> growthvectors(np);

         // Bit array to identify all the points belonging 
         // to the surface of interest
         bndnodes.Clear();

         // Run through all the surface elements and mark the points 
         // belonging to those where a boundary layer has to be created.
         // In addition, also calculate the effective surface normal 
         // vectors at each of those points to determine the mesh motion 
         // direction
         cout << "Marking points for remapping...." << endl;

         for (i = 1; i <= nse; i++)
         {
            int snr = mesh.SurfaceElement(i).GetIndex();
            // cout << "snr = " << snr << endl;
            if (surfid.Contains(snr))
            {
               Element2d & sel = mesh.SurfaceElement(i);
               int selNP = sel.GetNP();
               for(j = 1; j <= selNP; j++)
               {
                  // Set the bitarray to indicate that the 
                  // point is part of the required set
                  bndnodes.Set(sel.PNum(j));
		  
                  // Vec3d& surfacenormal = Vec3d();   ????
                  Vec3d surfacenormal;

                  // Calculate the surface normal at the current point 
                  // with respect to the current surface element
                  GetSurfaceNormal(mesh,sel,j,surfacenormal);
                  
                  // Add the surface normal to the already existent one 
                  // (This gives the effective normal direction at corners 
                  //  and curved areas)
                  growthvectors.Elem(sel.PNum(j)) = growthvectors.Elem(sel.PNum(j)) 
                                                    + surfacenormal;
               }
            }
         }

         // Add additional points into the mesh structure in order to 
         // clone the surface elements.
         // Also invert the growth vectors so that they point inwards, 
         // and normalize them
         cout << "Cloning points and calculating growth vectors...." << endl;

         for (i = 1; i <= np; i++)
         {
            if (bndnodes.Test(i))
            {
               mapto.Elem(i) = mesh.AddPoint (mesh.Point (i));

               growthvectors.Elem(i).Normalize();
               growthvectors.Elem(i) *= -1.0;
            }
            else
            {
               mapto.Elem(i) = 0;
               growthvectors.Elem(i) = Vec3d(0,0,0);
            }
         }


         // Add quad surface elements at edges for surfaces which 
         // dont have boundary layers

         // Bit array to keep track of segments already processed
         BitArray segsel(nseg);

         // Set them all to "1" to initially activate all segments
         segsel.Set();

         cout << "Adding 2D Quad elements on required surfaces...." << endl;

         for (i = 1; i <= nseg; i++)
         {
            int seg_p1 = mesh.LineSegment(i)[0];
            int seg_p2 = mesh.LineSegment(i)[1];

            // Only go in if the segment is still active, and if both its 
            // surface index is part of the "hit-list"
            if(segsel.Test(i) && surfid.Contains(mesh.LineSegment(i).si))
            {
               // clear the bit to indicate that this segment has been processed
               segsel.Clear(i);

               // Find matching segment pair on other surface
               for(j = 1; j <= nseg; j++)
               {
                  int segpair_p1 = mesh.LineSegment(j)[1];
                  int segpair_p2 = mesh.LineSegment(j)[0];

                  // Find the segment pair on the neighbouring surface element
                  // Identified by: seg1[0] = seg_pair[1] and seg1[1] = seg_pair[0]
                  if(segsel.Test(j) && ((segpair_p1 == seg_p1) && (segpair_p2 == seg_p2)))
                  {
                     // clear bit to indicate that processing of this segment is done
                     segsel.Clear(j);

                     // Only worry about those surfaces which are not in the 
                     // boundary layer list
                     if(!surfid.Contains(mesh.LineSegment(j).si))
                     {
                        int pnt_commelem = 0;
                        int pnum_commelem = 0;
                        Array<int> pnt1_elems;
                        Array<int> pnt2_elems;
                       
                            
                        meshtopo.GetVertexSurfaceElements(segpair_p1,pnt1_elems);
                        meshtopo.GetVertexSurfaceElements(segpair_p2,pnt2_elems);
                        for(int k = 1; k <= pnt1_elems.Size(); k++)
                        {
                           Element2d pnt1_sel = mesh.SurfaceElement(pnt1_elems.Elem(k));
                           for(int l = 1; l <= pnt2_elems.Size(); l++)
                           {
                              Element2d pnt2_sel = mesh.SurfaceElement(pnt2_elems.Elem(l));
                              if((pnt1_sel.GetIndex() == mesh.LineSegment(j).si) 
                                 && (pnt2_sel.GetIndex() == mesh.LineSegment(j).si)
                                 && (pnt1_elems.Elem(k) == pnt2_elems.Elem(l)))
                              {
                                 pnt_commelem = pnt1_elems.Elem(k);
                              }
                           }
                        }

                        for(int k = 1; k <= mesh.SurfaceElement(pnt_commelem).GetNP(); k++)
                        {
                           if((mesh.SurfaceElement(pnt_commelem).PNum(k) != segpair_p1)
                              && (mesh.SurfaceElement(pnt_commelem).PNum(k) != segpair_p2))
                           {
                              pnum_commelem = mesh.SurfaceElement(pnt_commelem).PNum(k);
                           }
                        }

                        Vec3d surfelem_vect, surfelem_vect1;
                        
                        Element2d & commsel = mesh.SurfaceElement(pnt_commelem);

                        dbg << "NP= " << commsel.GetNP() << " : ";

                        for(int k = 1; k <= commsel.GetNP(); k++)
                        {
                           GetSurfaceNormal(mesh,commsel,k,surfelem_vect1);
                           surfelem_vect += surfelem_vect1;
                        }

                        surfelem_vect.Normalize();

                        double surfangle = Angle(growthvectors.Elem(segpair_p1),surfelem_vect);

                        dbg << "V1= " << surfelem_vect1 
                            << " : V2= " << surfelem_vect1
                            << " : V= " << surfelem_vect
                            << " : GV= " << growthvectors.Elem(segpair_p1)
                            << " : Angle= " << surfangle * 180 / 3.141592;

                  
                        // remap the segments to the new points
                        mesh.LineSegment(i)[0] = mapto.Get(seg_p1);
                        mesh.LineSegment(i)[1] = mapto.Get(seg_p2);
                        mesh.LineSegment(j)[1] = mapto.Get(seg_p1);
                        mesh.LineSegment(j)[0] = mapto.Get(seg_p2);

                        if((surfangle < (90 + angleThreshold) * 3.141592 / 180.0)
                           && (surfangle > (90 - angleThreshold) * 3.141592 / 180.0))
                        {
                           dbg << " : quad\n";
                           // Since the surface is lower than the threshold, change the effective 
                           // prism growth vector to match with the surface vector, so that 
                           // the Quad which is created lies on the original surface
                           //growthvectors.Elem(segpair_p1) = surfelem_vect;

                           // Add a quad element to account for the prism volume
                           // element which is going to be added 
                           Element2d sel(QUAD);
                           sel.PNum(4) = mapto.Get(seg_p1);
                           sel.PNum(3) = mapto.Get(seg_p2);
                           sel.PNum(2) = segpair_p2;
                           sel.PNum(1) = segpair_p1;
                           sel.SetIndex(mesh.LineSegment(j).si);
                           mesh.AddSurfaceElement(sel);
                           numquads++;
                        }
                        else
                        {
                           dbg << "\n";
                           for (int k = 1; k <= pnt1_elems.Size(); k++)
                           {
                              Element2d & pnt_sel = mesh.SurfaceElement(pnt1_elems.Elem(k));
                              if(pnt_sel.GetIndex() == mesh.LineSegment(j).si)
                              {
                                 for(int l = 1; l <= pnt_sel.GetNP(); l++)
                                 {
                                    if(pnt_sel.PNum(l) == segpair_p1)
                                    {
                                       pnt_sel.PNum(l) = mapto.Get(seg_p1);
                                    }
                                    else if(pnt_sel.PNum(l) == segpair_p2)
                                    {
                                       pnt_sel.PNum(l) = mapto.Get(seg_p2);
                                    }
                                 }
                              }
                           }

                           for (int k = 1; k <= pnt2_elems.Size(); k++)
                           {
                              Element2d & pnt_sel = mesh.SurfaceElement(pnt2_elems.Elem(k));
                              if(pnt_sel.GetIndex() == mesh.LineSegment(j).si)
                              {
                                 for(int l = 1; l <= pnt_sel.GetNP(); l++)
                                 {
                                    if(pnt_sel.PNum(l) == segpair_p1)
                                    {
                                       pnt_sel.PNum(l) = mapto.Get(seg_p1);
                                    }
                                    else if(pnt_sel.PNum(l) == segpair_p2)
                                    {
                                       pnt_sel.PNum(l) = mapto.Get(seg_p2);
                                    }
                                 }
                              }
                           }
                        }
                     }
                     else
                     {
                        // If the code comes here, it indicates that we are at 
                        // a line segment pair which is at the intersection 
                        // of two surfaces, both of which have to grow boundary 
                        // layers.... here too, remapping the segments to the 
                        // new points is required
                        mesh.LineSegment(i)[0] = mapto.Get(seg_p1);
                        mesh.LineSegment(i)[1] = mapto.Get(seg_p2);
                        mesh.LineSegment(j)[1] = mapto.Get(seg_p1);
                        mesh.LineSegment(j)[0] = mapto.Get(seg_p2);
                     }
                  }
               }
            }
         }

         // Add prismatic cells at the boundaries
         cout << "Generating prism boundary layer volume elements...." << endl;

         for (i = 1; i <= nse; i++)
         {
            Element2d & sel = mesh.SurfaceElement(i);
            if(surfid.Contains(sel.GetIndex()))
            {
               Element el(PRISM);
               for (j = 1; j <= sel.GetNP(); j++)
               {
                  // Check (Doublecheck) if the corresponding point has a 
                  // copy available for remapping
                  if (mapto.Get(sel.PNum(j)))
                  {
                     // Define the points of the newly added Prism cell
                     el.PNum(j+3) = mapto.Get(sel.PNum(j));
                     el.PNum(j) = sel.PNum(j);
                  }
               }

               el.SetIndex(1);
               el.Invert();
               mesh.AddVolumeElement(el);
               numprisms++;
            }
         }

         // Finally switch the point indices of the surface elements 
         // to the newly added ones
         cout << "Transferring boundary layer surface elements to new vertex references...." << endl;

         for (i = 1; i <= nse; i++)
         {
            Element2d & sel = mesh.SurfaceElement(i);
            if(surfid.Contains(sel.GetIndex()))
            {
               for (j = 1; j <= sel.GetNP(); j++)
               {
                  // Check (Doublecheck) if the corresponding point has a 
                  // copy available for remapping
                  if (mapto.Get(sel.PNum(j)))
                  {
                     // Map the surface elements to the new points
                     sel.PNum(j) = mapto.Get(sel.PNum(j));
                  }
               }
            }
         }

         // Lock all the prism points so that the rest of the mesh can be 
         // optimised without invalidating the entire mesh
         for (i = 1; i <= np; i++)
         {
            if(bndnodes.Test(i)) mesh.AddLockedPoint(i);
         }

         // Now, actually pull back the old surface points to create 
         // the actual boundary layers
         cout << "Moving and optimising boundary layer points...." << endl;
         
         for (i = 1; i <= np; i++)
         {
            Array<int> vertelems;

            if(bndnodes.Test(i))
            {
               MeshPoint pointtomove;

               pointtomove = mesh.Point(i);

               if(layer == prismlayers)
               {
                  mesh.Point(i).SetPoint(pointtomove + layerht * growthvectors.Elem(i));

                  meshtopo.GetVertexElements(i,vertelems);

                  for(j = 1; j <= vertelems.Size(); j++)
                  {
		    // double sfact = 0.9;
                     Element volel = mesh.VolumeElement(vertelems.Elem(j));
                     if(((volel.GetType() == TET) || (volel.GetType() == TET10)) && (!volel.IsDeleted()))
                     {
                        //while((volel.Volume(mesh.Points()) <= 0.0) && (sfact >= 0.0))
                        //{
                        //   mesh.Point(i).SetPoint(pointtomove + (sfact * layerht * growthvectors.Elem(i)));
                        //   mesh.ImproveMesh();

                        //   // Try to move the point back by one step but 
                        //   // if the volume drops to below zero, double back
                        //   mesh.Point(i).SetPoint(pointtomove + ((sfact + 0.1) * layerht * growthvectors.Elem(i)));
                        //   if(volel.Volume(mesh.Points()) <= 0.0)
                        //   {
                        //      mesh.Point(i).SetPoint(pointtomove + (sfact * layerht * growthvectors.Elem(i)));
                        //   }
                        //   sfact -= 0.1;
                        //}
                        volel.Delete();
                     }
                  }

                  mesh.Compress();
               }
               else
               {
                  mesh.Point(i).SetPoint(pointtomove + layerht * growthvectors.Elem(i));
               }
            }
         }
      }

      // Optimise the tet part of the volume mesh after all the modifications 
      // to the system are completed
      //OptimizeVolume(mparam,mesh);

      cout << "New NP: " << mesh.GetNP() << endl;
      cout << "Num of Quads: " << numquads << endl;
      cout << "Num of Prisms: " << numprisms << endl;
      cout << "Boundary Layer Generation....Done!" << endl;

      dbg.close();
   }

}

