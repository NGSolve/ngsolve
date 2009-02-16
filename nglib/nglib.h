#ifndef NGLIB
#define NGLIB

/**************************************************************************/
/* File:   nglib.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   7. May. 2000                                                   */
/**************************************************************************/

/*!
   \file nglib.h
   \brief Library interface to the netgen meshing kernel
   \author Joachim Schoeberl
   \date 7. May 2000

   This header file provides access to the core functionality of the Netgen 
   Mesher via a library interface, without an interactive User Interface.

   The intention of providing these set of functions is to allow system 
   developers to integrate Netgen into top-level code, to act as the low 
   level mesh generation / optimisation kernel.  
*/

// Philippose - 14.02.2009
// Modifications for creating a DLL in Windows
#ifdef WIN32
   #ifdef NGLIB_EXPORTS
      #define DLL_HEADER   __declspec(dllexport)
   #else
      #define DLL_HEADER   __declspec(dllimport)
   #endif
#else
   #define DLL_HEADER 
#endif



// ** Constants used within Netgen *********************
/// Maximum allowed number of nodes per volume element
#define NG_VOLUME_ELEMENT_MAXPOINTS 10

/// Maximum allowed number of nodes per surface element
#define NG_SURFACE_ELEMENT_MAXPOINTS 6



// *** Data-types for accessing Netgen functionality ***
/// Data type for NETGEN mesh
typedef void * Ng_Mesh;

/// Data type for NETGEN CSG geometry
typedef void * Ng_CSG_Geometry;

/// Data type for NETGEN 2D geometry
typedef void * Ng_Geometry_2D;

/// Data type for NETGEN STL geometry
typedef void * Ng_STL_Geometry;



// *** Special Enum types used within Netgen ***********
/// Currently implemented surface element types
enum Ng_Surface_Element_Type 
   { NG_TRIG = 1, NG_QUAD = 2, NG_TRIG6 = 3 };

/// Currently implemented volume element types
enum Ng_Volume_Element_Type 
   { NG_TET = 1, NG_PYRAMID = 2, NG_PRISM = 3, NG_TET10 = 4 };

/// Values returned by Netgen functions
enum Ng_Result 
   { NG_OK                  = 0, 
     NG_SURFACE_INPUT_ERROR = 1,
     NG_VOLUME_FAILURE      = 2, 
     NG_STL_INPUT_ERROR     = 3,
     NG_SURFACE_FAILURE     = 4,
     NG_FILE_NOT_FOUND      = 5 
   };



// *** Classes required for use within Netgen **********
/// Netgen Meshing Parameters class
class Ng_Meshing_Parameters 
{
public:
   double maxh;                //!< Maximum global mesh size limit 
   double fineness;            //!< Mesh density: 0...1 (0 => coarse; 1 => fine)
   int secondorder;            //!< Generate second-order surface and volume elements
   char * meshsize_filename;   //!< Optional external mesh size file 
   int quad_dominated;         //!< Creates a Quad-dominated mesh 

   /*!
      Default constructor for the Mesh Parameters class

      Note: This constructor initialises the variables in the 
      class with the following default values
      - #maxh: 1000.0
      - #fineness: 0.5
      - #secondorder: 0.0
      - #meshsize_filename: 0
      - #quad_dominated: 0
   */
   DLL_HEADER Ng_Meshing_Parameters();
};




// *** Functions Exported by this Library *************

// General purpose initialisation / destruction functions

/*! \brief Initialise the Netgen library and prepare for use

    This function needs to be called by the third-party 
    program before beginning to use the other Netgen 
    specific functions.
*/
DLL_HEADER void Ng_Init ();


/*! \brief Exit the Netgen meshing kernel in a clean manner

    Use this function to exit the meshing sub-system in 
    a clean and orderly manner.
*/
DLL_HEADER void Ng_Exit ();
  

/*! \brief Create a new (and empty) Netgen Mesh Structure

    This function creates a new Netgen Mesh, initialises 
    it, and returns a pointer to the created mesh structure. 

    Use the returned pointer for subsequent operations 
    which involve mesh operations.

    \return Ng_Mesh* Pointer to a Netgen Mesh type #Ng_Mesh
*/
DLL_HEADER  Ng_Mesh * Ng_NewMesh ();


/*! \brief Delete an existing Netgen Mesh Structure

    Use this function to delete an existing Netgen mesh 
    structure and release the used memory. 

    \param mesh Pointer to an existing Netgen Mesh structure 
                of type #Ng_Mesh
*/
DLL_HEADER void Ng_DeleteMesh (Ng_Mesh * mesh);



// Common Mesh related utility functions

/*! \brief Add a point to a given Netgen Mesh Structure

    This function allows points to be directly added to a Netgen 
    mesh structure by providing the co-ordinates.

    Each call to the function allows only one point to be added.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param x    Pointer to an array of type double containing the co-ordinates 
                of the point to be added in the form: \n
                - x[0] = X co-ordinate
                - x[1] = Y co-ordinate
                - x[2] = Z co-ordinate
*/
DLL_HEADER void Ng_AddPoint (Ng_Mesh * mesh, double * x);


/*! \brief Add a surface element to a given Netgen Mesh Structure

    This function allows the top-level code to directly add individual 
    Surface Elements to a Netgen Mesh Structure by providing the type of 
    element to be added and the indices of the points which constitute the 
    element.

    <i>Note:</i> 
    - The points referred to by the surface elements must have been
      added prior to calling this function. 
    - Currently only triangular elements are supported, and the Surface Element 
      Type argument is not used.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param et   Surface Element type provided via the enumerated type 
                #Ng_Surface_Element_Type 
    \param pi   Pointer to an array of integers containing the indices of the 
                points which constitute the surface element being added
*/
DLL_HEADER void Ng_AddSurfaceElement (Ng_Mesh * mesh, Ng_Surface_Element_Type et, int * pi);


/*! \brief Add a volume element to a given Netgen Mesh Structure

    This function allows the top-level code to directly add individual 
    Volume Elements to a Netgen Mesh Structure by providing the type of 
    element to be added and the indices of the points which constitute the 
    element.

    <i>Note:</i> 
    - The points referred to by the volume elements must have been
      added prior to calling this function. 
    - Currently only tetrahedral elements are supported, and the Volume Element 
      Type argument is not used.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param et   Volume Element type provided via the enumerated type 
                #Ng_Volume_Element_Type 
    \param pi   Pointer to an array of integers containing the indices of the 
                points which constitute the volume element being added

*/
DLL_HEADER void Ng_AddVolumeElement (Ng_Mesh * mesh, Ng_Volume_Element_Type et, int * pi);
  


/*! \brief Returns the Number of Points present in the specified Mesh

    Given an already existent Netgen Mesh Structure, this function 
    returns the number of points currently present within the Mesh.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \return 
                Integer Data-type with the number of points in the Mesh
*/
DLL_HEADER int Ng_GetNP (Ng_Mesh * mesh);


/*! \brief Returns the Number of Surface Elements present in the specified Mesh

    Given an already existent Netgen Mesh Structure, this function 
    returns the number of surface elements currently present within 
    the Mesh.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \return 
                Integer Data-type with the number of surface elements in the Mesh
*/
DLL_HEADER int Ng_GetNSE (Ng_Mesh * mesh);


/*! \brief Returns the Number of Volume Elements present in the specified Mesh

    Given an already existent Netgen Mesh Structure, this function 
    returns the number of volume elements currently present within 
    the Mesh.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \return 
                Integer Data-type with the number of volume elements in the Mesh
*/
DLL_HEADER int Ng_GetNE (Ng_Mesh * mesh);


  
//  Return the Point Coordinates of a specified Point
// The x, y and z co-ordinates are returned in the array pointer as 
// x[0] = x ; x[1] = y ; x[2] = z
DLL_HEADER void Ng_GetPoint (Ng_Mesh * mesh, int num, double * x);



// return surface and volume element in pi
DLL_HEADER Ng_Surface_Element_Type 
Ng_GetSurfaceElement (Ng_Mesh * mesh, int num, int * pi);

DLL_HEADER Ng_Volume_Element_Type
Ng_GetVolumeElement (Ng_Mesh * mesh, int num, int * pi);


// Defines MeshSize Functions
DLL_HEADER void Ng_RestrictMeshSizeGlobal (Ng_Mesh * mesh, double h);
DLL_HEADER void Ng_RestrictMeshSizePoint (Ng_Mesh * mesh, double * p, double h);
DLL_HEADER void Ng_RestrictMeshSizeBox (Ng_Mesh * mesh, double * pmin, double * pmax, double h);
  
// generates volume mesh from surface mesh
DLL_HEADER Ng_Result Ng_GenerateVolumeMesh (Ng_Mesh * mesh, Ng_Meshing_Parameters * mp);

DLL_HEADER void Ng_SaveMesh(Ng_Mesh * mesh, const char* filename);
DLL_HEADER Ng_Mesh * Ng_LoadMesh(const char* filename);





// **********************************************************
// **   2D Meshing                                         **
// **********************************************************


// feeds points and boundary to mesh

DLL_HEADER void Ng_AddPoint_2D (Ng_Mesh * mesh, double * x);
DLL_HEADER void Ng_AddBoundarySeg_2D (Ng_Mesh * mesh, int pi1, int pi2);
  
// ask for number of points, elements and boundary segments
DLL_HEADER int Ng_GetNP_2D (Ng_Mesh * mesh);
DLL_HEADER int Ng_GetNE_2D (Ng_Mesh * mesh);
DLL_HEADER int Ng_GetNSeg_2D (Ng_Mesh * mesh);
  
//  return point coordinates
DLL_HEADER void Ng_GetPoint_2D (Ng_Mesh * mesh, int num, double * x);

// return 2d triangles
DLL_HEADER void Ng_GetElement_2D (Ng_Mesh * mesh, int num, int * pi, int * matnum = NULL);

// return 2d boundary segment
DLL_HEADER void Ng_GetSegment_2D (Ng_Mesh * mesh, int num, int * pi, int * matnum = NULL);


// load 2d netgen spline geometry
DLL_HEADER Ng_Geometry_2D * Ng_LoadGeometry_2D (const char * filename);

// generate 2d mesh, mesh is allocated by function
DLL_HEADER Ng_Result Ng_GenerateMesh_2D (Ng_Geometry_2D * geom,
                                         Ng_Mesh ** mesh,
                                         Ng_Meshing_Parameters * mp);
  
DLL_HEADER void Ng_HP_Refinement (Ng_Geometry_2D * geom,
                                  Ng_Mesh * mesh,
                                  int levels);
  




// **********************************************************
// **   STL Meshing                                        **
// **********************************************************


// loads geometry from STL file
DLL_HEADER Ng_STL_Geometry * Ng_STL_LoadGeometry (const char * filename, int binary = 0);


// generate new STL Geometry
DLL_HEADER Ng_STL_Geometry * Ng_STL_NewGeometry ();
  

// fills STL Geometry
// positive orientation
// normal vector may be null-pointer
DLL_HEADER void Ng_STL_AddTriangle (Ng_STL_Geometry * geom, 
                         double * p1, double * p2, double * p3, 
                         double * nv = NULL);

// add (optional) edges :
DLL_HEADER void Ng_STL_AddEdge (Ng_STL_Geometry * geom, 
                     double * p1, double * p2);

// after adding triangles (and edges) initialize
DLL_HEADER Ng_Result Ng_STL_InitSTLGeometry (Ng_STL_Geometry * geom);

// automatically generates edges:
DLL_HEADER Ng_Result Ng_STL_MakeEdges (Ng_STL_Geometry * geom,
                            Ng_Mesh* mesh,
                            Ng_Meshing_Parameters * mp);


// generates mesh, empty mesh must be already created.
DLL_HEADER Ng_Result Ng_STL_GenerateSurfaceMesh (Ng_STL_Geometry * geom,
                                      Ng_Mesh * mesh,
                                      Ng_Meshing_Parameters * mp);


#ifdef ACIS

// **********************************************************
// **   ACIS Meshing                                       **
// **********************************************************

/// Data type for NETGEN STL geomty
typedef void * Ng_ACIS_Geometry;

// loads geometry from STL file
DLL_HEADER Ng_ACIS_Geometry * Ng_ACIS_LoadGeometry (const char * filename);
  
// generates mesh, empty mesh must be already created.
DLL_HEADER Ng_Result Ng_ACIS_GenerateSurfaceMesh (Ng_ACIS_Geometry * geom,
                                                  Ng_Mesh * mesh,
                                                  Ng_Meshing_Parameters * mp);


#endif


#endif
