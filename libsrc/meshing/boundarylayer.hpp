#ifndef FILE_BOUNDARYLAYER
#define FILE_BOUNDARYLAYER


///
DLL_HEADER extern void InsertVirtualBoundaryLayer (Mesh & mesh);

/// Create a typical prismatic boundary layer on the given 
/// surfaces

class BoundaryLayerParameters
{
public:
  // parameters by Philippose ..
  Array<int> surfid;
  Array<double> heights;
  Array<double> new_matnrs;
  int prismlayers = 1;
  int bulk_matnr = 1;
  int new_matnr = 1;
  double hfirst = 0.01;
  double growthfactor = 1;
  bool optimize = true;
};

DLL_HEADER extern void GenerateBoundaryLayer (Mesh & mesh, BoundaryLayerParameters & blp);


#endif
