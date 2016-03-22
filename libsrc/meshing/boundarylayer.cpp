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

      int i;
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
         for (int j = 1; j <= el.GetNP(); j++)
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
   void GetSurfaceNormal(Mesh & mesh, const Element2d & el, int Vertex, Vec3d & SurfaceNormal)
   {
      int Vertex_A;
      int Vertex_B;

      Vertex_A = Vertex + 1;
      if(Vertex_A > el.GetNP()) Vertex_A = 1;

      Vertex_B = Vertex - 1;
      if(Vertex_B <= 0) Vertex_B = el.GetNP();

      Vec3d Vect_A,Vect_B;
      
      Vect_A = mesh[el.PNum(Vertex_A)] - mesh[el.PNum(Vertex)];
      Vect_B = mesh[el.PNum(Vertex_B)] - mesh[el.PNum(Vertex)];

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
   void GenerateBoundaryLayer (Mesh & mesh, BoundaryLayerParameters & blp)
   {
      ofstream dbg("BndLayerDebug.log");

      // Angle between a surface element and a growth-vector below which 
      // a prism is project onto that surface as a quad
      // (in degrees)
      double angleThreshold = 5.0;


      Array<int> surfid (blp.surfid);
      int prismlayers = blp.prismlayers;
      double hfirst = blp.hfirst;
      double growthfactor = blp.growthfactor;
      Array<double> heights (blp.heights);

      bool grow_edges = false; // grow layer at edges
      

      // Monitor and print out the number of prism and quad elements 
      // added to the mesh
      int numprisms = 0;
      int numquads = 0;
      

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
          
          if(heights.Size()>0)
            {
              layerht = heights[layer-1];
            }
          else
            {
              if(growthfactor == 1)
                {
                  layerht = layer * hfirst;
                }
              else
                {
                  layerht = hfirst*(pow(growthfactor,(layer+1)) - 1)/(growthfactor - 1);
                }
            }
          
         cout << "Layer Height = " << layerht << endl;
         
         // Need to store the old number of points and 
         // surface elements because there are new points and 
         // surface elements being added during the process
         int np = mesh.GetNP();
         int nse = mesh.GetNSE();
         int ne = mesh.GetNE();

         // Safety measure to ensure no issues with mesh 
         // consistency
         int nseg = mesh.GetNSeg();

         // Indicate which points need to be remapped
         BitArray bndnodes(np+1);  // big enough for 1-based array

         // Map of the old points to the new points
         Array<PointIndex, PointIndex::BASE> mapto(np);

         // Growth vectors for the prismatic layer based on 
         // the effective surface normal at a given point
         Array<Vec3d, PointIndex::BASE> growthvectors(np);

         // Bit array to identify all the points belonging 
         // to the surface of interest
         bndnodes.Clear();

         // Run through all the surface elements and mark the points 
         // belonging to those where a boundary layer has to be created.
         // In addition, also calculate the effective surface normal 
         // vectors at each of those points to determine the mesh motion 
         // direction
         cout << "Marking points for remapping...." << endl;
         
         for (SurfaceElementIndex si = 0; si < nse; si++)
           if (surfid.Contains(mesh[si].GetIndex()))
             {
               const Element2d & sel = mesh[si];
               for(int j = 0; j < sel.GetNP(); j++)
                 {
                   // Set the bitarray to indicate that the 
                   // point is part of the required set
                   bndnodes.Set(sel[j]);
                   Vec3d surfacenormal;
                   
                   // Calculate the surface normal at the current point 
                   // with respect to the current surface element
                   GetSurfaceNormal(mesh,sel,j+1,surfacenormal);
                   
                   // Add the surface normal to the already existent one 
                   // (This gives the effective normal direction at corners 
                   //  and curved areas)
                   growthvectors[sel[j]] += surfacenormal;
                 }
             }

         if (!grow_edges)
           for (SegmentIndex sei = 0; sei <= nseg; sei++)
             {
               bndnodes.Clear (mesh[sei][0]);
               bndnodes.Clear (mesh[sei][1]);
             }

         // Add additional points into the mesh structure in order to 
         // clone the surface elements.
         // Also invert the growth vectors so that they point inwards, 
         // and normalize them
         cout << "Cloning points and calculating growth vectors...." << endl;

         for (PointIndex pi = 1; pi <= np; pi++)
           {
             if (bndnodes.Test(pi))
               {
                 mapto[pi] = mesh.AddPoint (mesh[pi]);
                 
                 growthvectors[pi].Normalize();
                 growthvectors[pi] *= -1.0;
               }
             else
               {
                 mapto[pi] = 0;
                 growthvectors[pi] = Vec3d(0,0,0);
               }
         }


         // Add quad surface elements at edges for surfaces which 
         // dont have boundary layers

         // Bit array to keep track of segments already processed
         BitArray segsel(nseg);

         // Set them all to "1" to initially activate all segments
         segsel.Set();

         cout << "Adding 2D Quad elements on required surfaces...." << endl;

         if (grow_edges)
         for (SegmentIndex sei = 0; sei <= nseg; sei++)
           {
             PointIndex seg_p1 = mesh[sei][0];
             PointIndex seg_p2 = mesh[sei][1];
             
             // Only go in if the segment is still active, and if both its 
             // surface index is part of the "hit-list"
             if(segsel.Test(sei) && surfid.Contains(mesh[sei].si))
               {
                 // clear the bit to indicate that this segment has been processed
                 segsel.Clear(sei);
                 
                 // Find matching segment pair on other surface
                 for (SegmentIndex sej = 0; sej < nseg; sej++)
                   {
                     PointIndex segpair_p1 = mesh[sej][1];
                     PointIndex segpair_p2 = mesh[sej][0];
                     
                     // Find the segment pair on the neighbouring surface element
                     // Identified by: seg1[0] = seg_pair[1] and seg1[1] = seg_pair[0]
                     if(segsel.Test(sej) && ((segpair_p1 == seg_p1) && (segpair_p2 == seg_p2)))
                       {
                         // clear bit to indicate that processing of this segment is done
                         segsel.Clear(sej);
                         
                         // Only worry about those surfaces which are not in the 
                         // boundary layer list
                         if(!surfid.Contains(mesh[sej].si))
                           {
                             SurfaceElementIndex pnt_commelem = 0;
                             Array<SurfaceElementIndex> pnt1_elems;
                             Array<SurfaceElementIndex> pnt2_elems;
                             
                            
                             meshtopo.GetVertexSurfaceElements(segpair_p1,pnt1_elems);
                             meshtopo.GetVertexSurfaceElements(segpair_p2,pnt2_elems);

                             for(int k = 0; k < pnt1_elems.Size(); k++)
                               {
                                 const Element2d & pnt1_sel = mesh.SurfaceElement(pnt1_elems[k]);
                                 for(int l = 0; l < pnt2_elems.Size(); l++)
                                   {
                                     const Element2d & pnt2_sel = mesh.SurfaceElement(pnt2_elems[l]);
                                     if((pnt1_sel.GetIndex() == mesh[sej].si) 
                                        && (pnt2_sel.GetIndex() == mesh[sej].si)
                                        && (pnt1_elems[k] == pnt2_elems[l]))
                                       {
                                         pnt_commelem = pnt1_elems[k];
                                       }
                                   }
                               }

                             /*
                               int pnum_commelem = 0;
                               for(int k = 1; k <= mesh.SurfaceElement(pnt_commelem).GetNP(); k++)
                               {
                               if((mesh.SurfaceElement(pnt_commelem).PNum(k) != segpair_p1)
                               && (mesh.SurfaceElement(pnt_commelem).PNum(k) != segpair_p2))
                               {
                               pnum_commelem = mesh.SurfaceElement(pnt_commelem).PNum(k);
                               }
                               }
                             */
                             
                             Vec3d surfelem_vect, surfelem_vect1;
                        
                             const Element2d & commsel = mesh.SurfaceElement(pnt_commelem);
                             
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
                             mesh[sei][0] = mapto[seg_p1];
                             mesh[sei][1] = mapto[seg_p2];
                             mesh[sej][1] = mapto[seg_p1];
                             mesh[sej][0] = mapto[seg_p2];
                             
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
                                 sel.PNum(4) = mapto[seg_p1];
                                 sel.PNum(3) = mapto[seg_p2];
                                 sel.PNum(2) = segpair_p2;
                                 sel.PNum(1) = segpair_p1;
                                 sel.SetIndex(mesh[sej].si);
                                 mesh.AddSurfaceElement(sel);
                                 numquads++;
                               }
                             else
                               {
                                 dbg << "\n";
                                 for (int k = 0; k < pnt1_elems.Size(); k++)
                                   {
                                     Element2d & pnt_sel = mesh.SurfaceElement(pnt1_elems[k]);
                                     if(pnt_sel.GetIndex() == mesh[sej].si)
                                       {
                                         for(int l = 0; l < pnt_sel.GetNP(); l++)
                                           {
                                             if(pnt_sel[l] == segpair_p1)
                                               pnt_sel[l] = mapto[seg_p1];
                                             else if (pnt_sel[l] == segpair_p2)
                                               pnt_sel[l] = mapto[seg_p2];
                                           }
                                       }
                                   }
                                 
                                 for (int k = 0; k < pnt2_elems.Size(); k++)
                                   {
                                     Element2d & pnt_sel = mesh.SurfaceElement(pnt2_elems[k]);
                                     if(pnt_sel.GetIndex() == mesh[sej].si)
                                       {
                                         for(int l = 0; l < pnt_sel.GetNP(); l++)
                                           {
                                             if(pnt_sel[l] == segpair_p1)
                                               pnt_sel[l] = mapto.Get(seg_p1);
                                             else if (pnt_sel[l] == segpair_p2)
                                               pnt_sel[l] = mapto.Get(seg_p2);
                                           }
                                       }
                                   }
                               }
                             // }
                           }
                         else
                           {
                             // If the code comes here, it indicates that we are at 
                             // a line segment pair which is at the intersection 
                             // of two surfaces, both of which have to grow boundary 
                             // layers.... here too, remapping the segments to the 
                             // new points is required
                             mesh[sei][0] = mapto.Get(seg_p1);
                             mesh[sei][1] = mapto.Get(seg_p2);
                             mesh[sej][1] = mapto.Get(seg_p1);
                             mesh[sej][0] = mapto.Get(seg_p2);
                           }
                       }
                   }
               }
           }
         
         // Add prismatic cells at the boundaries
         cout << "Generating prism boundary layer volume elements...." << endl;

         for (SurfaceElementIndex si = 0; si < nse; si++)
           {
             Element2d & sel = mesh.SurfaceElement(si);
             if(surfid.Contains(sel.GetIndex()))
               {
                 /*
                 Element el(PRISM);
                 for (int j = 0; j < sel.GetNP(); j++)
                   {
                     // Check (Doublecheck) if the corresponding point has a 
                     // copy available for remapping
                     if (mapto.Get(sel[j]))
                       {
                         // Define the points of the newly added Prism cell
                         el[j+3] = mapto[sel[j]];
                         el[j] = sel[j];
                       }
                     else
                       {
                         el[j+3] = sel[j];
                         el[j] = sel[j];
                       }
                   }
                 
                 el.SetIndex(1);
                 el.Invert();
                 mesh.AddVolumeElement(el);
                 numprisms++;
                 */
                 // cout << "add element: " << endl;
                 int classify = 0;
                 for (int j = 0; j < 3; j++)
                   if (mapto[sel[j]])
                     classify += (1 << j);

                 // cout << "classify = " << classify << endl;

                 ELEMENT_TYPE types[] = { PRISM, TET, TET, PYRAMID, 
                                          TET, PYRAMID, PYRAMID, PRISM };
                 int nums[] = { sel[0], sel[1], sel[2], mapto[sel[0]], mapto[sel[1]], mapto[sel[2]] };
                 int vertices[][6] = 
                   { 
                     { 0, 1, 2, 0, 1, 2 },   // should not occur
                     { 0, 2, 1, 3, 0, 0 },
                     { 0, 2, 1, 4, 0, 0 },
                     { 0, 1, 4, 3, 2, 0 },

                     { 0, 2, 1, 5, 0, 0 },
                     { 2, 0, 3, 5, 1, 0 }, 
                     { 1, 2, 5, 4, 0, 0 },
                     { 0, 2, 1, 3, 5, 4 }
                   };

                 Element el(types[classify]);
                 for (int i = 0; i < 6; i++)
                   el[i] = nums[vertices[classify][i]];
                   if(blp.new_matnrs.Size() > 0)
                      el.SetIndex(blp.new_matnrs[layer-1]);
                   else
                      el.SetIndex(blp.new_matnr);
                 // cout << "el = " << el << endl;
                 if (classify != 0)
                   mesh.AddVolumeElement(el);
               }
           }
         
         // Finally switch the point indices of the surface elements 
         // to the newly added ones
         cout << "Transferring boundary layer surface elements to new vertex references...." << endl;
         
         for (int i = 1; i <= nse; i++)
           {
             Element2d & sel = mesh.SurfaceElement(i);
             if(surfid.Contains(sel.GetIndex()))
              {
                for (int j = 1; j <= sel.GetNP(); j++)
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
         for (int i = 1; i <= ne; i++)
           {
             Element & el = mesh.VolumeElement(i);
             if(el.GetIndex() != blp.bulk_matnr)
              {
                for (int j = 1; j <= el.GetNP(); j++)
                  {
                    // Check (Doublecheck) if the corresponding point has a 
                    // copy available for remapping
                    if (mapto.Get(el.PNum(j)))
                      {
                        // Map the surface elements to the new points
                        el.PNum(j) = mapto.Get(el.PNum(j));
                      }
                  }
              }
           }



         
         // Lock all the prism points so that the rest of the mesh can be 
         // optimised without invalidating the entire mesh
         for (PointIndex pi = mesh.Points().Begin(); pi < mesh.Points().End(); pi++)
         {
           if(bndnodes.Test(pi)) mesh.AddLockedPoint(pi);
         }

         // Now, actually pull back the old surface points to create 
         // the actual boundary layers
         cout << "Moving and optimising boundary layer points...." << endl;
         
         for (int i = 1; i <= np; i++)
           {
            Array<ElementIndex> vertelems;

            if(bndnodes.Test(i))
              {
                MeshPoint pointtomove;
                
                pointtomove = mesh.Point(i);
                
                if(layer == prismlayers)
                  {
                    mesh.Point(i).SetPoint(pointtomove + layerht * growthvectors.Elem(i));
                    
                    meshtopo.GetVertexElements(i,vertelems);
                    
                    for(int j = 1; j <= vertelems.Size(); j++)
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
		  }
                else
                  {
                    mesh.Point(i).SetPoint(pointtomove + layerht * growthvectors.Elem(i));
                  }
              }
           }
	 mesh.Compress();
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

