#ifndef FILE_BCFUNCTIONS
#define FILE_BCFUNCTIONS

// Philippose - 14/03/2009
// Auxiliary functions for OCC Geometry
// Use this file and the corresponding ".cpp" 
// file to add miscellaneous functionality 
// to the OpenCascade Geometry support in Netgen
namespace netgen
{
   /*! \brief Automatically assign boundary conditions for meshes

       This function allows the boundary condition numbers of a 
       mesh created in Netgen to be automatically assigned based on 
       the colours of each face.

       Currently, two algorithms are utilised to assign the BC Properties:
       1. Automatic assignment using a user defined colour profile file 
          which defines which RGB colours are to be assigned to which 
          BC Property number
          - A default profile file exists in the Netgen folder called 
            "netgen.ocf"
       
       2. The second algorithm uses the following automated algorithm:
          - Extract all the colours present in the mesh
          - Use colour index 0 (zero) for all faces with no colour defined
          - Calculate the number of faces of the surface mesh for each colour
          - Sort the number of surface elements in ascending order, with the 
            colour indices as a slave
          - Use the indices of the sorted array as the BC property number

          Example: If there are 3 colours, present in the mesh and the number 
          of surface elements for each colour are:
          - Colour 0: 8500
          - Colour 1: 120
          - Colour 2: 2200
          - Colour 3: 575

          The above is sorted in ascending order and assigned as BC Properties:
          - BC Prop 0: 120  : Colour 1
          - BC Prop 1: 575  : Colour 3
          - BC Prop 2: 2200 : Colour 2
          - BC Prop 3: 8500 : Colour 0 (no colour defined)
   */
   //extern void OCCAutoColourBcProps(Mesh & mesh, OCCGeometry & occgeometry, const char *occcolourfile);
   extern void AutoColourBcProps(Mesh & mesh, const char *bccolourfile);

   extern void GetFaceColours(Mesh & mesh, Array<Vec3d> & face_colours);

   extern bool ColourMatch(Vec3d col1, Vec3d col2, double eps = 2.5e-05);
}
#endif

