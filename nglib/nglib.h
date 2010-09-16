#ifndef NGLIB
#define NGLIB

/**************************************************************************/
/* File:   nglib.h                                                        */
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
   #ifdef NGLIB_EXPORTS || nglib_EXPORTS
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
#define NG_SURFACE_ELEMENT_MAXPOINTS 8



// *** Data-types for accessing Netgen functionality ***
/// Data type for NETGEN mesh
typedef void * Ng_Mesh;

/// Data type for NETGEN CSG geometry
typedef void * Ng_CSG_Geometry;

/// Data type for NETGEN 2D geometry
typedef void * Ng_Geometry_2D;

/// Data type for NETGEN STL geometry
typedef void * Ng_STL_Geometry;

#ifdef OCCGEOMETRY
/// Data type for NETGEN OpenCascade geometry
typedef void * Ng_OCC_Geometry;
typedef void * Ng_OCC_TopTools_IndexedMapOfShape;
#endif


// *** Special Enum types used within Netgen ***********
/// Currently implemented surface element types
enum Ng_Surface_Element_Type 
   { NG_TRIG = 1, NG_QUAD = 2, NG_TRIG6 = 3, NG_QUAD6 = 4, NG_QUAD8 = 5 };

/// Currently implemented volume element types
enum Ng_Volume_Element_Type 
   { NG_TET = 1, NG_PYRAMID = 2, NG_PRISM = 3, NG_TET10 = 4 };

/// Values returned by Netgen functions
enum Ng_Result 
   { 
     NG_ERROR               = -1,   
     NG_OK                  = 0, 
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
   int uselocalh;                      //!< Switch to enable / disable usage of local mesh size modifiers

   double maxh;                        //!< Maximum global mesh size allowed
   double minh;                        //!< Minimum global mesh size allowed

   double fineness;                    //!< Mesh density: 0...1 (0 => coarse; 1 => fine)
   double grading;                     //!< Mesh grading: 0...1 (0 => uniform mesh; 1 => aggressive local grading)

   double elementsperedge;             //!< Number of elements to generate per edge of the geometry
   double elementspercurve;            //!< Elements to generate per curvature radius

   int closeedgeenable;                //!< Enable / Disable mesh refinement at close edges
   double closeedgefact;               //!< Factor to use for refinement at close edges (larger => finer)

   int second_order;                   //!< Generate second-order surface and volume elements
   int quad_dominated;                 //!< Creates a Quad-dominated mesh 

   char * meshsize_filename;           //!< Optional external mesh size file 

   int optsurfmeshenable;              //!< Enable / Disable automatic surface mesh optimization
   int optvolmeshenable;               //!< Enable / Disable automatic volume mesh optimization

   int optsteps_3d;                     //!< Number of optimize steps to use for 3-D mesh optimization
   int optsteps_2d;                     //!< Number of optimize steps to use for 2-D mesh optimization

   // Philippose - 13/09/2010
   // Added a couple more parameters into the meshing parameters list 
   // from Netgen into Nglib
   int invert_tets;                    //!< Invert all the volume elements
   int invert_trigs;                   //!< Invert all the surface triangle elements

   int check_overlap;                  //!< Check for overlapping surfaces during Surface meshing
   int check_overlapping_boundary;     //!< Check for overlapping surface elements before volume meshing


   /*!
      Default constructor for the Mesh Parameters class

      Note: This constructor initialises the variables in the 
      class with the following default values
      - #uselocalh: 1
      - #maxh: 1000.0
      - #fineness: 0.5
      - #grading: 0.3
      - #elementsperedge: 2.0
      - #elementspercurve: 2.0
      - #closeedgeenable: 0
      - #closeedgefact: 2.0
      - #secondorder: 0
      - #meshsize_filename: null
      - #quad_dominated: 0
      - #optsurfmeshenable: 1
      - #optvolmeshenable: 1
      - #optsteps_2d: 3
      - #optsteps_3d: 3
      - #invert_tets: 0
      - #invert_trigs:0 
      - #check_overlap: 1
      - #check_overlapping_boundary: 1
   */
   DLL_HEADER Ng_Meshing_Parameters();



   /*!
       Reset the meshing parameters to their defaults

       This member function resets all the meshing parameters 
       of the object to the default values
   */
   DLL_HEADER void Reset_Parameters();



   /*!
       Transfer local meshing parameters to internal meshing parameters

       This member function transfers all the meshing parameters 
       defined in the local meshing parameters structure of nglib into 
       the internal meshing parameters structure used by the Netgen core
   */
   DLL_HEADER void Transfer_Parameters();
};




// *** Functions Exported by this Library *************

// ------------------------------------------------------------------
// Netgen library initialisation / destruction functions

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

    \return Ng_Mesh Pointer to a Netgen Mesh type #Ng_Mesh
*/
DLL_HEADER  Ng_Mesh * Ng_NewMesh ();


/*! \brief Delete an existing Netgen Mesh Structure

    Use this function to delete an existing Netgen mesh 
    structure and release the used memory. 

    \param mesh Pointer to an existing Netgen Mesh structure 
                of type #Ng_Mesh
*/
DLL_HEADER void Ng_DeleteMesh (Ng_Mesh * mesh);


/*! \brief Save a Netgen Mesh to disk

    This function allows a generated mesh structure to be saved 
    to disk.

    A Mesh saved using this function, will be written to disk 
    in the Netgen VOL file format.

    \param mesh    Pointer to an existing Netgen Mesh structure 
                   of type #Ng_Mesh
    \param filename Pointer to a character array containing the 
                    name of the file to which the mesh should 
                    be saved
*/
DLL_HEADER void Ng_SaveMesh(Ng_Mesh * mesh, const char* filename);


/*! \brief Load a Netgen VOL Mesh from disk into memory

    A Netgen mesh saved in the internal VOL format can be loaded 
    into a Netgen Mesh structure using this function. 

    \param filename Pointer to a character array containing the 
                    name of the file to load
    \return Ng_Mesh Pointer to a Netgen Mesh type #Ng_Mesh containing 
                    the mesh loaded from disk
*/
DLL_HEADER Ng_Mesh * Ng_LoadMesh(const char* filename);


/*! \brief Merge a Netgen VOL Mesh from disk into an existing mesh in memory

    A Netgen mesh saved in the internal VOL format can be merged 
    into an existing Netgen Mesh structure using this function. 

    \param mesh       Name of the Mesh structure already existent in memory
    \param filename   Pointer to a character array containing the 
                      name of the file to load
    \return Ng_Result Status of the merge operation
*/
DLL_HEADER Ng_Result Ng_MergeMesh(Ng_Mesh * mesh, const char* filename);


/*! \brief Merge one Netgen Mesh into another Netgen Mesh in the case 
    when both are already in memory

    (NOTE: FUNCTION STILL WORK IN PROGRESS!!!)

    This function can be used to merge two Netgen meshes already present 
    in memory.

    \param mesh1      Parent Mesh structure into which the second mesh 
                      will be merged
    \param mesh2      Child mesh structure which will get merged into 
                      the parent mesh
    \return Ng_Result Status of the merge operation
*/
DLL_HEADER Ng_Result Ng_MergeMesh(Ng_Mesh * mesh1, Ng_Mesh * mesh2);
// ------------------------------------------------------------------



// ------------------------------------------------------------------
// Basic Meshing functions for manually adding points, surface elements 
// and volume elements to a Netgen Mesh structure

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
  
// ------------------------------------------------------------------



// ------------------------------------------------------------------
// Local Mesh Size restriction / limiting utilities

/*! \brief Apply a global restriction on mesh element size

    This utility allows the user to apply a global mesh element 
    size limitation. 

    During mesh creation, in the absence of an explicit local 
    size restriction around the neighbourhood of a point within 
    the meshing domain, this global size restriction will be 
    utilised.

    <b>Note</b>: This function only limits the <b>Maximum</b> 
    size of an element within the mesh.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param h    Variable of type <i>double</i>, specifying the maximum
                allowable mesh size
*/
DLL_HEADER void Ng_RestrictMeshSizeGlobal (Ng_Mesh * mesh, double h);


/*! \brief Locally restrict the mesh element size at the given point

    Unlike the function #Ng_RestrictMeshSizeGlobal, this function 
    allows the user to locally restrict the maximum allowable mesh 
    size at a given point.

    The point is specified via its three cartesian co-ordinates.

    <b>Note</b>: This function only limits the <b>Maximum</b> size 
    of the elements around the specified point.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param p    Pointer to an Array of type <i>double</i>, containing 
                the three co-ordinates of the point in the form: \n
                - p[0] = X co-ordinate
                - p[1] = Y co-ordinate
                - p[2] = Z co-ordinate
    \param h    Variable of type <i>double</i>, specifying the maximum
                allowable mesh size at that point
*/
DLL_HEADER void Ng_RestrictMeshSizePoint (Ng_Mesh * mesh, double * p, double h);


/*! \brief Locally restrict the mesh element size within a specified box

    Similar to the function #Ng_RestrictMeshSizePoint, this function 
    allows the size of elements within a mesh to be locally limited.

    However, rather than limit the mesh size at a single point, this 
    utility restricts the local mesh size within a 3D Box region, specified 
    via the co-ordinates of the two diagonally opposite points of a cuboid.

    <b>Note</b>: This function only limits the <b>Maximum</b> size 
    of the elements within the specified region.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param pmin Pointer to an Array of type <i>double</i>, containing 
                the three co-ordinates of the first point of the cuboid: \n
                - pmin[0] = X co-ordinate
                - pmin[1] = Y co-ordinate
                - pmin[2] = Z co-ordinate
    \param pmax Pointer to an Array of type <i>double</i>, containing 
                the three co-ordinates of the opposite point of the 
                cuboid: \n
                - pmax[0] = X co-ordinate
                - pmax[1] = Y co-ordinate
                - pmax[2] = Z co-ordinate
    \param h    Variable of type <i>double</i>, specifying the maximum
                allowable mesh size at that point
*/
DLL_HEADER void Ng_RestrictMeshSizeBox (Ng_Mesh * mesh, double * pmin, double * pmax, double h);

// ------------------------------------------------------------------



// ------------------------------------------------------------------
// 3D Mesh Generation functions

/*! \brief Create a 3D Volume Mesh given a Surface Mesh

    After creating a surface mesh, this function can be utilised 
    to automatically generate the corresponding 3D Volume Mesh.

    Mesh generation parameters (such as grading, maximum element size, 
    etc.) are specified via the meshing parameters class which also 
    needs to be passed to this function.

    <b>Note</b>: Currently, Netgen generates pure tetrahedral volume 
    meshes.

    \param mesh Pointer to an existing Netgen Mesh structure of 
                type #Ng_Mesh
    \param mp   Pointer to a copy of the Meshing Parameters class
                (#Ng_Meshing_Parameters), filled up with the 
                required values

    \return Ng_Result Status of the Mesh Generation routine. More 
                      details regarding the return value can be 
                      found in the description of #Ng_Result
*/
DLL_HEADER Ng_Result Ng_GenerateVolumeMesh (Ng_Mesh * mesh, Ng_Meshing_Parameters * mp);

// ------------------------------------------------------------------



// ------------------------------------------------------------------
// Basic Mesh information functions

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

// ------------------------------------------------------------------



// ------------------------------------------------------------------
// Mesh Topology functions
// Use these functions to extract points, surface / volume elements, 
// perform topological searches, etc..etc...
  
//  Return the Point Coordinates of a specified Point
// The x, y and z co-ordinates are returned in the array pointer as 
// x[0] = x ; x[1] = y ; x[2] = z
DLL_HEADER void Ng_GetPoint (Ng_Mesh * mesh, int num, double * x);



// return surface and volume element in pi
DLL_HEADER Ng_Surface_Element_Type 
Ng_GetSurfaceElement (Ng_Mesh * mesh, int num, int * pi);

DLL_HEADER Ng_Volume_Element_Type
Ng_GetVolumeElement (Ng_Mesh * mesh, int num, int * pi);

// ------------------------------------------------------------------




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

// return 2d elements
DLL_HEADER Ng_Surface_Element_Type 
Ng_GetElement_2D (Ng_Mesh * mesh, int num, int * pi, int * matnum = NULL);

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



#ifdef OCCGEOMETRY

// **********************************************************
// **   OpenCascade Geometry / Meshing Utilities           **
// **********************************************************

// Create new OCC Geometry Object
DLL_HEADER Ng_OCC_Geometry * Ng_OCC_NewGeometry ();

// Delete an OCC Geometry Object
DLL_HEADER Ng_Result Ng_OCC_DeleteGeometry (Ng_OCC_Geometry * geom);

// Loads geometry from STEP file
DLL_HEADER Ng_OCC_Geometry * Ng_OCC_Load_STEP (const char * filename);

// Loads geometry from IGES file
DLL_HEADER Ng_OCC_Geometry * Ng_OCC_Load_IGES (const char * filename);

// Loads geometry from BREP file
DLL_HEADER Ng_OCC_Geometry * Ng_OCC_Load_BREP (const char * filename);

// Set the local mesh size based on geometry / topology
DLL_HEADER Ng_Result Ng_OCC_SetLocalMeshSize (Ng_OCC_Geometry * geom,
                                              Ng_Mesh * mesh,
                                              Ng_Meshing_Parameters * mp);

// Mesh the edges and add Face descriptors to prepare for surface meshing
DLL_HEADER Ng_Result Ng_OCC_GenerateEdgeMesh (Ng_OCC_Geometry * geom,
                                              Ng_Mesh * mesh,
                                              Ng_Meshing_Parameters * mp);

// Mesh the surfaces of an OCC geometry
DLL_HEADER Ng_Result Ng_OCC_GenerateSurfaceMesh (Ng_OCC_Geometry * geom,
                                                 Ng_Mesh * mesh,
                                                 Ng_Meshing_Parameters * mp); 

// Get the face map of an already loaded OCC geometry
DLL_HEADER Ng_Result Ng_OCC_GetFMap(Ng_OCC_Geometry * geom, 
                                    Ng_OCC_TopTools_IndexedMapOfShape * FMap);

#endif // OCCGEOMETRY



// **********************************************************
// **   Mesh refinement algorithms                         **
// **********************************************************

// uniform mesh refinement
DLL_HEADER void Ng_Uniform_Refinement (Ng_Mesh * mesh);


// uniform mesh refinement with geometry adaption:

DLL_HEADER void Ng_2D_Uniform_Refinement (Ng_Geometry_2D * geom,
					  Ng_Mesh * mesh);

DLL_HEADER void Ng_STL_Uniform_Refinement (Ng_STL_Geometry * geom,
					   Ng_Mesh * mesh);

DLL_HEADER void Ng_CSG_Uniform_Refinement (Ng_CSG_Geometry * geom,
					   Ng_Mesh * mesh);

#ifdef OCCGEOMETRY
DLL_HEADER void Ng_OCC_Uniform_Refinement (Ng_OCC_Geometry * geom,
					   Ng_Mesh * mesh);
#endif



// **********************************************************
// **   Second Order mesh algorithms                       **
// **********************************************************

// convert mesh to second order
DLL_HEADER void Ng_Generate_SecondOrder (Ng_Mesh * mesh);


// convert mesh to second order with geometry adaption:

DLL_HEADER void Ng_2D_Generate_SecondOrder (Ng_Geometry_2D * geom,
					  Ng_Mesh * mesh);

DLL_HEADER void Ng_STL_Generate_SecondOrder (Ng_STL_Geometry * geom,
					   Ng_Mesh * mesh);

DLL_HEADER void Ng_CSG_Generate_SecondOrder (Ng_CSG_Geometry * geom,
					   Ng_Mesh * mesh);

#ifdef OCCGEOMETRY
DLL_HEADER void Ng_OCC_Generate_SecondOrder (Ng_OCC_Geometry * geom,
					   Ng_Mesh * mesh);
#endif


#endif // NGLIB
