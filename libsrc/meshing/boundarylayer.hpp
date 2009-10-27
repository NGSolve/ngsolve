#ifndef FILE_BOUNDARYLAYER
#define FILE_BOUNDARYLAYER


///
extern void InsertVirtualBoundaryLayer (Mesh & mesh);

/// Create a typical prismatic boundary layer on the given 
/// surfaces
extern void GenerateBoundaryLayer (Mesh & mesh, MeshingParameters & mp);


#endif
