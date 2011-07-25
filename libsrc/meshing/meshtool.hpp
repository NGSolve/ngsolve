#ifndef FILE_MESHTOOL
#define FILE_MESHTOOL


///
extern void MeshQuality2d (const Mesh & mesh);

///
extern void MeshQuality3d (const Mesh & mesh,
			   Array<int> * inclass = NULL);

///
extern void SaveEdges (const Mesh & mesh, 
		       const char * geomfile, 
		       double h, 
		       char * filename);

///
extern void SaveSurfaceMesh (const Mesh & mesh,
			     double h,
			     char * filename);
/*
///
extern void Save2DMesh (
         const Mesh & mesh2d,
	 const Array<class SplineSegment*> * splines,
         ostream & outfile);
*/

class Surface;
///
extern void SaveVolumeMesh (
         const Array<Point3d> & points,
         const Array<Element> & elements,
         const Array<Element> & volelements,
         const Array<Surface*> & surfaces,
         char * filename);

///
void SaveVolumeMesh (const Mesh & mesh, 
		     const class CSGeometry & geometry,
		     char * filename);

///
extern int CheckCode ();


///
extern double CalcTetBadness (const Point3d & p1, const Point3d & p2,
			      const Point3d & p3, const Point3d & p4, 
			      double h,
			      const MeshingParameters & mp);
///
extern double CalcTetBadnessGrad (const Point3d & p1, const Point3d & p2,
				  const Point3d & p3, const Point3d & p4, 
				  double h, int pi,
				  Vec<3> & grad,
				  const MeshingParameters & mp);


/** Calculates volume of an element.
  The volume of the tetrahedron el is computed
 */
// extern double CalcVolume (const Array<Point3d> & points,
//        const Element & el);  

/** The total volume of all elements is computed.
  This function calculates the volume of the mesh */
extern double CalcVolume (const Array<Point3d> & points, 
	const Array<Element> & elements);

///
extern int CheckSurfaceMesh (const Mesh & mesh);

///
extern int CheckSurfaceMesh2 (const Mesh & mesh);
///
extern int CheckMesh3D (const Mesh & mesh);
///
extern void RemoveProblem (Mesh & mesh, int domainnr);
#endif
