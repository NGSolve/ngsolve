/*! \file writegmsh2.cpp
*  \brief Export Netgen Mesh in the GMSH v2.xx File format
*  \author Philippose Rajan
*  \date 02 November 2008
*
*  This function extends the export capabilities of
*  Netgen to include the GMSH v2.xx File Format.
*
*  Current features of this function include:
*
*  1. Exports Triangles, Quadrangles and Tetrahedra \n
*  2. Supports upto second order elements of each type
*
*/


#include <mystdlib.h>

#include <myadt.hpp>
#include <linalg.hpp>
#include <csg.hpp>
#include <meshing.hpp>

namespace netgen
{
#include "writeuser.hpp"

   // Mapping of entities from Netgen definitions to GMSH definitions
   enum GMSH_ELEMENTS {GMSH_TRIG = 2, GMSH_TRIG6 = 9,
      GMSH_QUAD = 3, GMSH_QUAD8 = 16,
      GMSH_TET = 4, GMSH_TET10 = 11};
   const int triGmsh[7] = {0,1,2,3,6,4,5};
   const int quadGmsh[9] = {0,1,2,3,4,5,8,6,7};
   const int tetGmsh[11] = {0,1,2,3,4,5,8,6,7,10,9};


   /*! GMSH v2.xx mesh format export function
   *
   *  This function extends the export capabilities of
   *  Netgen to include the GMSH v2.xx File Format.
   *
   *  Current features of this function include:
   *
   *  1. Exports Triangles, Quadrangles and Tetrahedra \n
   *  2. Supports upto second order elements of each type
   *
   */
   void WriteGmsh2Format (const Mesh & mesh,
      const CSGeometry & geom,
      const string & filename)
   {
      ofstream outfile (filename.c_str());
      outfile.precision(6);
      outfile.setf (ios::fixed, ios::floatfield);
      outfile.setf (ios::showpoint);

      int np = mesh.GetNP();  /// number of points in mesh
      int ne = mesh.GetNE();  /// number of 3D elements in mesh
      int nse = mesh.GetNSE();  /// number of surface elements (BC)
      int i, j, k, l;


      /*
      * 3D section : Volume elements (currently only tetrahedra)
      */

      if ((ne > 0)
         && (mesh.VolumeElement(1).GetNP() <= 10)
         && (mesh.SurfaceElement(1).GetNP() <= 6))
      {
         cout << "Write GMSH v2.xx Format \n";
         cout << "The GMSH v2.xx export is currently available for elements upto 2nd Order\n" << endl;

         int inverttets = mparam.inverttets;
         int invertsurf = mparam.inverttrigs;

         /// Prepare GMSH 2.0 file (See GMSH 2.0 Documentation)
         outfile << "$MeshFormat\n";
         outfile << (float)2.0 << " "
            << (int)0 << " "
            << (int)sizeof(double) << "\n";
         outfile << "$EndMeshFormat\n";

         /// Write nodes
         outfile << "$Nodes\n";
         outfile << np << "\n";

         for (i = 1; i <= np; i++)
         {
            const Point3d & p = mesh.Point(i);
            outfile << i << " "; /// node number
            outfile << p.X() << " ";
            outfile << p.Y() << " ";
            outfile << p.Z() << "\n";
         }

         outfile << "$EndNodes\n";

         /// write elements (both, surface elements and volume elements)
         outfile << "$Elements\n";
         outfile << ne + nse << "\n";  ////  number of elements + number of surfaces BC

         for (i = 1; i <= nse; i++)
         {
            int elType = 0;

            Element2d el = mesh.SurfaceElement(i);
            if(invertsurf) el.Invert();

            if(el.GetNP() == 3) elType = GMSH_TRIG;	//// GMSH Type for a 3 node triangle
            if(el.GetNP() == 6) elType = GMSH_TRIG6;  //// GMSH Type for a 6 node triangle
            if(elType == 0)
            {
               cout << " Invalid surface element type for Gmsh 2.0 3D-Mesh Export Format !\n";
               return;
            }

            outfile << i;
            outfile << " ";
            outfile << elType;
            outfile << " ";
            outfile << "2";                  //// Number of tags (2 => Physical and elementary entities)
            outfile << " ";
            outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << " ";
            /// that means that physical entity = elementary entity (arbitrary approach)
            outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << " ";
            for (j = 1; j <= el.GetNP(); j++)
            {
               outfile << " ";
               outfile << el.PNum(triGmsh[j]);
            }
            outfile << "\n";
         }


         for (i = 1; i <= ne; i++)
         {
            int elType = 0;

            Element el = mesh.VolumeElement(i);
            if (inverttets) el.Invert();

            if(el.GetNP() == 4) elType = GMSH_TET;    //// GMSH Element type for 4 node tetrahedron
            if(el.GetNP() == 10) elType = GMSH_TET10; //// GMSH Element type for 10 node tetrahedron
            if(elType == 0)
            {
               cout << " Invalid volume element type for Gmsh 2.0 3D-Mesh Export Format !\n";
               return;
            }

            outfile << nse + i;                       //// element number (Remember to add on surface elements)
            outfile << " ";
            outfile << elType;
            outfile << " ";
            outfile << "2";                   //// Number of tags (2 => Physical and elementary entities)
            outfile << " ";
            outfile << 100000 + el.GetIndex();
            /// that means that physical entity = elementary entity (arbitrary approach)
            outfile << " ";
            outfile << 100000 + el.GetIndex();   /// volume number
            outfile << " ";
            for (j = 1; j <= el.GetNP(); j++)
            {
               outfile << " ";
               outfile << el.PNum(tetGmsh[j]);
            }
            outfile << "\n";
         }
         outfile << "$EndElements\n";
      }
      /*
      * End of 3D section
      */


      /*
      * 2D section : available for triangles and quadrangles
      *              upto 2nd Order
      */
      else if(ne == 0)   /// means that there's no 3D element
      {
         cout << "\n Write Gmsh v2.xx Surface Mesh (triangle and/or quadrangles upto 2nd Order)" << endl;

         /// Prepare GMSH 2.0 file (See GMSH 2.0 Documentation)
         outfile << "$MeshFormat\n";
         outfile << (float)2.0 << " "
            << (int)0 << " "
            << (int)sizeof(double) << "\n";
         outfile << "$EndMeshFormat\n";

         /// Write nodes
         outfile << "$Nodes\n";
         outfile << np << "\n";

         for (i = 1; i <= np; i++)
         {
            const Point3d & p = mesh.Point(i);
            outfile << i << " "; /// node number
            outfile << p.X() << " ";
            outfile << p.Y() << " ";
            outfile << p.Z() << "\n";
         }
         outfile << "$EndNodes\n";

         /// write triangles & quadrangles
         outfile << "$Elements\n";
         outfile << nse << "\n";

         for (k = 1; k <= nse; k++)
         {
            int elType = 0;

            const Element2d & el = mesh.SurfaceElement(k);

            if(el.GetNP() == 3) elType = GMSH_TRIG;   //// GMSH Type for a 3 node triangle
            if(el.GetNP() == 6) elType = GMSH_TRIG6;  //// GMSH Type for a 6 node triangle
            if(el.GetNP() == 4) elType = GMSH_QUAD;   //// GMSH Type for a 4 node quadrangle
            if(el.GetNP() == 8) elType = GMSH_QUAD8;  //// GMSH Type for an 8 node quadrangle
            if(elType == 0)
            {
               cout << " Invalid surface element type for Gmsh 2.0 2D-Mesh Export Format !\n";
               return;
            }

            outfile << k;
            outfile << " ";
            outfile << elType;
            outfile << " ";
            outfile << "2";
            outfile << " ";
            outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << " ";
            /// that means that physical entity = elementary entity (arbitrary approach)
            outfile << mesh.GetFaceDescriptor (el.GetIndex()).BCProperty() << " ";
            for (l = 1; l <= el.GetNP(); l++)
            {
               outfile << " ";
               if((elType == GMSH_TRIG) || (elType == GMSH_TRIG6))
               {
                  outfile << el.PNum(triGmsh[l]);
               }
               else if((elType == GMSH_QUAD) || (elType == GMSH_QUAD8))
               {
                  outfile << el.PNum(quadGmsh[l]);
               }
            }
            outfile << "\n";
         }
         outfile << "$EndElements\n";
      }
      /*
      * End of 2D section
      */

      else
      {
         cout << " Invalid element type for Gmsh v2.xx Export Format !\n";
      }
   } // End: WriteGmsh2Format
} // End: namespace netgen


