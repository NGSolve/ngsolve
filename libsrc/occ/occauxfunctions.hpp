#ifndef FILE_OCCAUXFUNCTIONS
#define FILE_OCCAUXFUNCTIONS

// Philippose - 14/03/2009
// Auxiliary functions for OCC Geometry
// Use this file and the corresponding ".cpp" 
// file to add miscellaneous functionality 
// to the OpenCascade Geometry support in Netgen

namespace netgen
{
   /*! \brief Automatically assign boundary conditions for OCC meshes

       This function allows the boundary condition numbers of a 
       mesh created using an OpenCascade (STEP / IGES) geometry 
       to be assigned automatically.

       The boundary conditions are assigned based on the face 
       colour information (if any) contained in the geometry.

       Currently the following process is used to assign the BC Properties:
       - Extract all the colours present in the OCC Geometry
       - Use colour index 0 (zero) for all faces with no colour defined
       - Calculate the number of faces of the surface mesh for each colour
       - Sort the number of surface elements in ascending order, with the 
         colour indices as a slave
       - Use the indices of the sorted array as the BC property number

       Example: If there are 3 colours, present in the file and the number 
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
   extern void OCCAutoColourBcProps(Mesh & mesh, OCCGeometry & occgeometry, const char *occcolourfile);

}
#endif

