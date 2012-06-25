#include <ngstd.hpp>
#include <bla.hpp>

#include <nginterface.h>
#include <nginterface_v2.hpp>

using namespace std;
#include <ngexception.hpp>


using namespace ngbla;

/*
int dim = 2;
int nv = 4;
int nseg = 0;
int ntrig = 2;

double points [][2] =
  { { 0, 0 },
    { 1, 0 },
    { 0, 1 },
    { 1, 1 } 
  };

int trigs[][3] = 
  { 
    { 1, 2, 3},
    { 2, 3, 4} 
  };

*/



// one-1 mesh example


int dim = 1;
int nv = 5;
int npoint = 2;  // geom points
int nseg = 4;
int ntrig = 0;

double points[][1] = { 0, 0.25, 0.5, 0.75, 1 };
int trigs[][3] = { { 0,0,0 } } ;

int segms[][2] = 
  {
    { 2, 1 },
    { 3, 2 },
    { 4, 3 },
    { 5, 4 } 
  };

int pnts[][1] =
  { { 1 }, { 5 } };





// here are the tested functions


// space dimension (2 or 3)
DLL_HEADER int Ng_GetDimension ()
{
  return dim;
}


DLL_HEADER int Ng_GetNNodes (int nt)
{
  switch (nt)
    {
    case 0: return nv;
    case 1: return 0;
    case 2: return 0;
    case 3: return 0;
    }
}















// just dummies, old style

extern "C"
{

  DLL_HEADER void Ng_LoadGeometry (const char * filename)
  {
    cout << "Ng-dummy ... cannot load geometry" << endl;
  }

  DLL_HEADER void Ng_LoadMesh (const char * filename)
  {
    cout << "Ng-dummy ... cannot load mesh" << endl;
  }

  DLL_HEADER void Ng_LoadMeshFromString (const char * mesh_as_string)
  {
    ;
  }


  // number of mesh points
  int Ng_GetNP ()
  {
    return nv;
  }





  // number of mesh elements
  DLL_HEADER int Ng_GetNE ()
  {
    return Ng_GetNElements (dim);
  }

  // number of surface triangles
  DLL_HEADER int Ng_GetNSE ()
  {
    return Ng_GetNElements (dim-1);
  }
  
  // Get Point coordintes, index from 1 .. np
  DLL_HEADER void Ng_GetPoint (int pi, double * p)
  {
    for (int i = 0; i < dim; i++)
      p[i] = points[pi][i];
  }
  
  // Get Element Points
  DLL_HEADER NG_ELEMENT_TYPE Ng_GetElement (int ei, int * epi, int * np)
  {
    for (int i = 0; i < dim+1; i++)
      epi[i] = trigs[ei][i];
    if (np) *np = dim+1;
    return NG_TRIG;
  }


  // Get Element Type
  DLL_HEADER NG_ELEMENT_TYPE Ng_GetElementType (int ei)
  {
    return NG_TRIG;
  }

  // Get sub-domain of element ei
  DLL_HEADER int Ng_GetElementIndex (int ei)
  {
    return 1;
  }

  DLL_HEADER void Ng_SetElementIndex(const int ei, const int index)
  {
    ;
  }

  // Get Material of element ei
  DLL_HEADER char * Ng_GetElementMaterial (int ei)
  {
    return (char*)"dummy";
  }

  // Get Material of domain dom
  DLL_HEADER char * Ng_GetDomainMaterial (int dom)
  {
    return (char*)"dummy";
  }
  
  // Get User Data
  DLL_HEADER int Ng_GetUserDataSize (char * id)
  {
    return 0;
  }

  DLL_HEADER void Ng_GetUserData (char * id, double * data)
  {
    ;
  }

  // Get Surface Element Points
  DLL_HEADER NG_ELEMENT_TYPE Ng_GetSurfaceElement (int ei, int * epi, int * np)
  {
    ;
  }

  // Get Surface Element Type
  DLL_HEADER NG_ELEMENT_TYPE Ng_GetSurfaceElementType (int ei)
  {
    ;
  }

  // Get Surface Element Index
  DLL_HEADER int Ng_GetSurfaceElementIndex (int ei)
  {
    return 1;
  }

  // Get Surface Element Surface Number
  DLL_HEADER int Ng_GetSurfaceElementSurfaceNumber (int ei)
  {
    return 0;
  }
  
  // Get Surface Element Number
  DLL_HEADER int Ng_GetSurfaceElementFDNumber (int ei)
  {
    return 0;
  }

  // Get BCName for Surface Element  
  DLL_HEADER char * Ng_GetSurfaceElementBCName (int ei)
  {
    return (char*)"bound";
  }


  DLL_HEADER char * Ng_GetBCNumBCName (int bcnr)
  {
    return (char*)"bound";
  }


  // Get normal vector of surface element node
  DLL_HEADER void Ng_GetNormalVector (int sei, int locpi, double * nv)
  {
    ;
  }
  

  DLL_HEADER void Ng_SetPointSearchStartElement(int el)
  {
    ;
  }
  
  // Find element of point, returns local coordinates
  DLL_HEADER int Ng_FindElementOfPoint (double * p, double * lami,
                                        int build_searchtrees, 
                                        const int * const indices, const int numind)
  {
    ;
  }
  
  // Find surface element of point, returns local coordinates
  DLL_HEADER int Ng_FindSurfaceElementOfPoint (double * p, double * lami,
                                               int build_searchtrees, 
                                               const int * const indices, const int numind)
  {
    ;
  }
  

  // is elment ei curved ?
  DLL_HEADER int Ng_IsElementCurved (int ei)
  {
    return false;
  }

  // is elment sei curved ?
  DLL_HEADER int Ng_IsSurfaceElementCurved (int sei)
  {
    return false;
  }

  /// Curved Elemens:
  /// xi..local coordinates
  /// x ..global coordinates
  /// dxdxi...D x D Jacobian matrix (row major storage)
  DLL_HEADER void Ng_GetElementTransformation (int ei, const double * xi, 
                                               double * x, double * dxdxi)
  {
    cout << "geteltrafo" << endl;
  }

  
  /// buffer must be at least 100 doubles, alignment of double
  DLL_HEADER void Ng_GetBufferedElementTransformation (int ei, const double * xi, 
                                                       double * x, double * dxdxi,
                                                       void * buffer, int buffervalid)
  {
    cout << "getbufferedtrafo" << endl;
  }
  


  /// Curved Elemens:
  /// xi..local coordinates
  /// x ..global coordinates
  /// dxdxi...D x D-1 Jacobian matrix (row major storage)
  /// curved ...is element curved ?
  DLL_HEADER void Ng_GetSurfaceElementTransformation (int sei, const double * xi, 
                                                      double * x, double * dxdxi)
  {
    cout << "get sel trafo" << endl;
  }
  
  /// Curved Elemens:
  /// xi..local coordinates
  /// sxi..step xi
  /// x ..global coordinates
  /// dxdxi...D x D Jacobian matrix (row major storage)
  DLL_HEADER void Ng_GetMultiElementTransformation (int ei, int n,
                                                    const double * xi, size_t sxi,
                                                    double * x, size_t sx,
                                                    double * dxdxi, size_t sdxdxi)
  {
    cout << "nget_getmultielementtrafo" << endl;
  }

   
  DLL_HEADER int Ng_GetSegmentIndex (int elnr)
  {
    return 0;
  }

  DLL_HEADER NG_ELEMENT_TYPE Ng_GetSegment (int elnr, int * epi, int * np)
  {
    ;
  }




  // Mark element for refinement
  DLL_HEADER void Ng_SetRefinementFlag (int ei, int flag) { ; }
  DLL_HEADER void Ng_SetSurfaceRefinementFlag (int sei, int flag) { ; }

  // Do local refinement

  DLL_HEADER void Ng_Refine (NG_REFINEMENT_TYPE reftype)  { ; }

  // Use second order elements
  DLL_HEADER void Ng_SecondOrder () { ; }
  DLL_HEADER void Ng_HighOrder (int order, bool rational) { ; }
  //void Ng_HPRefinement (int levels, double parameter = 0.125);
  DLL_HEADER void Ng_HPRefinement (int levels, double parameter,
				   bool setorders,bool ref_level) { ; }


  // Topology and coordinate information of master element:
  DLL_HEADER int Ng_ME_GetNVertices (NG_ELEMENT_TYPE et) { ; }
  DLL_HEADER int Ng_ME_GetNEdges (NG_ELEMENT_TYPE et) { ; }
  DLL_HEADER int Ng_ME_GetNFaces (NG_ELEMENT_TYPE et) { ; }

  DLL_HEADER const NG_POINT * Ng_ME_GetVertices (NG_ELEMENT_TYPE et);
  DLL_HEADER const NG_EDGE * Ng_ME_GetEdges (NG_ELEMENT_TYPE et);
  DLL_HEADER const NG_FACE * Ng_ME_GetFaces (NG_ELEMENT_TYPE et);

  DLL_HEADER void Ng_UpdateTopology() { ; }

  DLL_HEADER int Ng_GetNEdges() { return 0; }
  DLL_HEADER int Ng_GetNFaces() { return 0; }

  
  DLL_HEADER int Ng_GetElement_Edges (int elnr, int * edges, int * orient) { ; }
  DLL_HEADER int Ng_GetElement_Faces (int elnr, int * faces, int * orient) { ; }

  DLL_HEADER int Ng_GetSurfaceElement_Edges (int selnr, int * edges, int * orient) { ; }
  DLL_HEADER int Ng_GetSurfaceElement_Face (int selnr, int * orient) { ; }

  DLL_HEADER void Ng_GetSurfaceElementNeighbouringDomains(const int selnr, int & in, int & out) { ; }
       
  DLL_HEADER int Ng_GetFace_Vertices (int fnr, int * vert) { ; }
  DLL_HEADER void Ng_GetEdge_Vertices (int ednr, int * vert) { ; }
  DLL_HEADER int Ng_GetFace_Edges (int fnr, int * edge) { ; }

  DLL_HEADER int Ng_GetNVertexElements (int vnr) { ; }
  DLL_HEADER void Ng_GetVertexElements (int vnr, int * els) { ; }

  DLL_HEADER int Ng_GetElementOrder (int enr) { ; }
  DLL_HEADER void Ng_GetElementOrders (int enr, int * ox, int * oy, int * oz) { ; }

  DLL_HEADER void Ng_SetElementOrder (int enr, int order) { ; }
  DLL_HEADER void Ng_SetElementOrders (int enr, int ox, int oy, int oz) { ; }

  DLL_HEADER int Ng_GetSurfaceElementOrder (int enr) { ; }
  DLL_HEADER void Ng_GetSurfaceElementOrders (int enr, int * ox, int * oy) { ; }

  DLL_HEADER void Ng_SetSurfaceElementOrder (int enr, int order) { ; }
  DLL_HEADER void Ng_SetSurfaceElementOrders (int enr, int ox, int oy) { ; }

  // Multilevel functions:

  // number of levels:
  DLL_HEADER int Ng_GetNLevels () { return 1; }
  // get two parent nodes (indeed vertices !) of node ni
  DLL_HEADER void Ng_GetParentNodes (int ni, int * parents) { ; }

  // get parent element (first child has always same number)
  DLL_HEADER int Ng_GetParentElement (int ei) { ; }

  // get parent surface element (first child has always same number)
  DLL_HEADER int Ng_GetParentSElement (int ei);

  // representant of anisotropic cluster
  DLL_HEADER int Ng_GetClusterRepVertex (int vi) { ; }
  DLL_HEADER int Ng_GetClusterRepEdge (int edi) { ; }
  DLL_HEADER int Ng_GetClusterRepFace (int fai) { ; }
  DLL_HEADER int Ng_GetClusterRepElement (int eli) { ; }


  void Ng_SurfaceElementTransformation (int eli, double x, double y, 
					double * p3d, double * jacobian)
  {
    cout << "surfaceelementtrafo" << endl;
  }


  // the folling functions are 0-base  !!

  // number on distant processor 
  // returns pairs  (dist_proc, num_on_dist_proc)
  int NgPar_GetDistantNodeNums ( int nodetype, int locnum, int * pnums ) { ; } 
  int NgPar_GetNDistantNodeNums ( int nodetype, int locnum ) { ; }
  
  int NgPar_GetGlobalNodeNum (int nodetype, int locnum) { ; }

  
  
  // initialize solution data with default arguments
  DLL_HEADER void Ng_InitSolutionData (Ng_SolutionData * soldata) { ; }
  // set solution data
  DLL_HEADER void Ng_SetSolutionData (Ng_SolutionData * soldata) { ; }
  /// delete gridfunctions
  DLL_HEADER void Ng_ClearSolutionData() { ; }
  // redraw 
  DLL_HEADER void Ng_Redraw() { ; }
  //
  DLL_HEADER void Ng_SetVisualizationParameter (const char * name, 
                                                const char * value) { ; }
  

  // number of periodic vertices  
  DLL_HEADER int Ng_GetNPeriodicVertices (int idnr) { ; }
  // pairs should be an integer array of 2*npairs
  DLL_HEADER void Ng_GetPeriodicVertices (int idnr, int * pairs) { ; }

  // number of periodic edges  
  DLL_HEADER int Ng_GetNPeriodicEdges (int idnr) { ; }


  // pairs should be an integer array of 2*npairs
  DLL_HEADER void Ng_GetPeriodicEdges (int idnr, int * pairs) { ; }

  DLL_HEADER void RunParallel ( void * (*fun)(void *), void * in) { ; }

  DLL_HEADER void Ng_PushStatus (const char * str) { ; }
  DLL_HEADER void Ng_PopStatus () { ; }
  DLL_HEADER void Ng_SetThreadPercentage (double percent) { ; }
  DLL_HEADER void Ng_GetStatus (char ** str, double & percent) { ; }

  DLL_HEADER void Ng_SetTerminate(void) { ; }
  DLL_HEADER void Ng_UnSetTerminate(void) { ; }
  DLL_HEADER int Ng_ShouldTerminate(void) { ; }
  DLL_HEADER void Ng_SetRunning(int flag) { ; }
  DLL_HEADER int Ng_IsRunning() { return false; }
  
  //// added by Roman Stainko ....
  DLL_HEADER int Ng_GetVertex_Elements( int vnr, int* elems) { ; }
  DLL_HEADER int Ng_GetVertex_SurfaceElements( int vnr, int* elems ) { ; }
  DLL_HEADER int Ng_GetVertex_NElements( int vnr ) { ; }
  DLL_HEADER int Ng_GetVertex_NSurfaceElements( int vnr ) { ; }


  DLL_HEADER void Ng_InitPointCurve(double red, double green, double blue) { ; }
  DLL_HEADER void Ng_AddPointCurvePoint(const double * point) { ; }


  DLL_HEADER void Ng_SaveMesh ( const char * meshfile );
  DLL_HEADER void Ng_Bisect ( const char * refinementfile );

  // if qualityloss is not equal to NULL at input, a (1-based) list of qualitylosses (due to projection)
  // is saved in *qualityloss, its size is the return value
  DLL_HEADER int Ng_Bisect_WithInfo ( const char * refinementfile, double ** qualityloss);

  typedef void * Ng_Mesh;
  DLL_HEADER Ng_Mesh Ng_SelectMesh (Ng_Mesh mesh);


}









/*
  The new node interface ...
  it is 0-based !
*/

extern "C" {
  
  /*
    number of nodes of type nt
    nt = 0 is Vertex
    nt = 1 is Edge
    nt = 2 is Face
    nt = 3 is Cell
  */

  /*
    closure nodes of node (nt, nodenr):
    nodeset is bit-coded, bit 0 includes Vertices, bit 1 edges, etc
    E.g., nodeset = 6 includes edge and face nodes
    nodes consists of pairs of integers (nodetype, nodenr) 
    return value is number of nodes
  */
  DLL_HEADER int Ng_GetClosureNodes (int nt, int nodenr, int nodeset, int * nodes);

  
  /*
    number of dim-dimensional elements 
    dim = 3 ... volume elements
    dim = 2 ... surface elements
    dim = 1 ... segments
    dim = 0 ... not available
  */
  DLL_HEADER int Ng_GetNElements (int dim)
  {
    switch (dim)
      {
      case 0: return npoint;
      case 1: return nseg;
      case 2: return ntrig;
      case 3: return 0;
      }
  }

  /*
    closure nodes of dim-dimensional element elmentnr:
    nodeset is bit-coded, bit 0 includes Vertices, bit 1 edges, etc
    E.g., nodeset = 6 includes edge and face nodes
    nodes consists of pairs of integers (nodetype, nodenr) 
    return value is number of nodes
  */
  DLL_HEADER int Ng_GetElementClosureNodes (int dim, int elementnr, int nodeset, int * nodes) { ; }


  struct Ng_Tcl_Interp;
  typedef int (Ng_Tcl_CmdProc) (Ng_Tcl_Interp *interp, int argc, const char *argv[]);

  DLL_HEADER void Ng_Tcl_CreateCommand (Ng_Tcl_Interp * interp, 
                                        const char * cmdName, Ng_Tcl_CmdProc * proc);

  void Ng_Tcl_SetResult (Ng_Tcl_Interp * interp, const char * result);

}






// dummies, V2
namespace netgen
{
  ostream * testout = &cout;
  int printmessage_importance = 1;

  Ng_Point Ng_GetPoint (int nr)
  {
    ;
  }

  template<> Ng_Element Ng_GetElement<0> (int nr) 
  { 
    Ng_Element ret;
    ret.type = NG_PNT;
    ret.points.num = 1;
    ret.points.ptr  = &pnts[nr][0];
    ret.vertices.num = 1;
    ret.vertices.ptr  = &pnts[nr][0];
    ret.edges.num = 0;
    ret.faces.num = 0;
    return ret;
  }


  template<> Ng_Element Ng_GetElement<1> (int nr) 
  { 
    Ng_Element ret;
    ret.type = NG_SEGM; 
    ret.points.num = 2;
    ret.points.ptr  = &segms[nr][0];
    ret.vertices.num = 2;
    ret.vertices.ptr  = &segms[nr][0];

    ret.edges.num = 0;
    ret.edges.ptr = NULL;
    ret.faces.num = 0;
    ret.faces.ptr = NULL;
    return ret;
  }

  template<> Ng_Element Ng_GetElement<2> (int nr) 
  { 
    Ng_Element ret;
    ret.type = NG_TRIG; 
    ret.points.num = 3;
    ret.points.ptr  = &trigs[nr][0];

    ret.vertices.num = 3;
    ret.vertices.ptr = &trigs[nr][0];

    ret.edges.num = 0;
    ret.edges.ptr = NULL;

    ret.faces.num = 0;
    ret.faces.ptr = NULL;
    return ret;
  }

  template<> Ng_Element Ng_GetElement<3> (int nr) 
  { 
    cout << "calls ng_getelement<3>" << endl;
  }

  template<> Ng_Node<1> Ng_GetNode<1> (int nr)
  {
    ;
  }

  template<> Ng_Node<2> Ng_GetNode<2> (int nr)
  {
    ;
  }


  template <>
  void Ng_MultiElementTransformation<1,1> (int elnr, int npts,
					   const double * xi, size_t sxi,
					   double * x, size_t sx,
					   double * dxdxi, size_t sdxdxi)
  {
    double pts[2] = { points[segms[elnr][0]-1][0], 
		      points[segms[elnr][1]-1][0] };

    for (int i = 0; i < npts; i++)
      {
	double xi_ = xi[i*sxi];
	if (x)
	  {
	    x[i*sx] = pts[1] + xi_ * (pts[0]-pts[1]);
	  }
	if (dxdxi)
	  {
	    dxdxi[i*sdxdxi  ] = pts[0]-pts[1];
	  }
      }
  }



  template <>
  void Ng_MultiElementTransformation<1,2> (int elnr, int npts,
					   const double * xi, size_t sxi,
					   double * x, size_t sx,
					   double * dxdxi, size_t sdxdxi)
  {
    cout << "multi - trafo<1,2>" << endl;
  }

  template <>
  void Ng_MultiElementTransformation<2,2> (int elnr, int npts,
					   const double * xi, size_t sxi,
					   double * x, size_t sx,
					   double * dxdxi, size_t sdxdxi)
  {
    Vec<2> pts[3];
    for (int i = 0; i < 3; i++)
      pts[i] = Vec<2> (points[trigs[elnr][i]-1][0], 
		       points[trigs[elnr][i]-1][1]);

    for (int i = 0; i < npts; i++)
      {
	Vec<2> xi_ (xi[i*sxi], xi[i*sxi+1]);
	Vec<2> x_map = pts[2] + xi_(0) * (pts[0]-pts[2]) + xi_(1) * (pts[1]-pts[2]);
	Mat<2,2> grad; 
	grad.Col(0) = pts[0]-pts[2];
	grad.Col(1) = pts[1]-pts[2];
	if (x)
	  {
	    x[i*sx  ] = x_map(0);
	    x[i*sx+1] = x_map(1);
	  }
	if (dxdxi)
	  {
	    dxdxi[i*sdxdxi  ] = grad(0,0);
	    dxdxi[i*sdxdxi+1] = grad(0,1);
	    dxdxi[i*sdxdxi+2] = grad(1,0);
	    dxdxi[i*sdxdxi+3] = grad(1,1);
	  }
      }
  }
  
  template <>
  void Ng_MultiElementTransformation<2,3> (int elnr, int npts,
					   const double * xi, size_t sxi,
					   double * x, size_t sx,
					   double * dxdxi, size_t sdxdxi)
  {
    cout << "multi - trafo<2,3>" << endl;
  }
  
  template <>
  void Ng_MultiElementTransformation<3,3> (int elnr, int npts,
					   const double * xi, size_t sxi,
					   double * x, size_t sx,
					   double * dxdxi, size_t sdxdxi)
  {
    cout << "multi - trafo<3,3>" << endl;
  }
  


  
}




























namespace netgen
{

  BaseDynamicMem * BaseDynamicMem::first = 0;
  BaseDynamicMem * BaseDynamicMem::last = 0;


  BaseDynamicMem :: BaseDynamicMem ()
  {
    prev = last;
    next = 0;

    if (last) last->next = this;
    last = this;
    if (!first) first = this;

    size = 0;
    ptr = 0;
    name = 0;
  }
 
  BaseDynamicMem :: ~BaseDynamicMem ()
  {
    Free();

    if (next) next->prev = prev;
    else last = prev;
    if (prev) prev->next = next;
    else first = next;

    delete [] name;
  }

  void BaseDynamicMem :: SetName (const char * aname)
  {
    delete [] name;
    name = NULL;
    if (aname)
      {
	name = new char[strlen(aname)+1];
	strcpy (name, aname);
      }
  }


  void BaseDynamicMem :: Alloc (size_t s)
  {
    size = s;
    ptr = new char[s];

    if (!ptr)
      {
	cerr << "BaseynamicMem, cannot allocate " << s << " bytes" << endl;
	Print ();
	throw ("BaseDynamicMem::Alloc: out of memory");
      }
    // ptr = (char*)malloc (s);
    // ptr = (char*) _mm_malloc (s,16);
  }

  void BaseDynamicMem :: ReAlloc (size_t s)
  {
    if (size == s) return;

    char * old = ptr;
    ptr = new char[s];

    if (!ptr)
      {
	cerr << "BaseynamicMem, cannot Reallocate " << s << " bytes" << endl;
	Print ();
	throw ("BaseDynamicMem::Alloc: out of memory");
      }


    // ptr = (char*)malloc(s);
    // ptr = (char*) _mm_malloc (s,16);
    memmove (ptr, old, (s < size) ? s : size);
    delete [] old;
    // free (old);
    // _mm_free (old);
    size = s;
  }

  void BaseDynamicMem :: Free ()
  {
    delete [] ptr;
    // free (ptr);
    // _mm_free (ptr);
    ptr = 0;
  }

  void BaseDynamicMem :: Swap (BaseDynamicMem & m2)
  {
    size_t hi;
    char * cp;
    hi = size; size  = m2.size; m2.size = hi;
    cp = ptr; ptr = m2.ptr; m2.ptr = cp;
    cp = name; name = m2.name; m2.name = cp;
  }


  void BaseDynamicMem :: Print ()
  {
    cout << "****************** Dynamic Mem Report ****************" << endl;
    BaseDynamicMem * p = first;
    size_t mem = 0;
    int cnt = 0;
    while (p)
      {
	mem += p->size;
	cnt++;

	cout << setw(10) << p->size << " Bytes";
	cout << ", addr = " << (void*)p->ptr;
	if (p->name)
	  cout << " in block " << p->name;
	cout << endl;

	p = p->next;
      }

    if (mem > 100000000)
      cout << "memory in dynamic memory: " << mem/1048576 << " MB" << endl;
    else if (mem > 100000)
      cout << "memory in dynamic memory: " << mem/1024 << " kB" << endl;
    else
      cout << "memory in dynamic memory: " << mem << " Bytes" << endl;
    cout << "number of blocks:         " << cnt << endl;
    //  cout << "******************************************************" << endl;
  }


#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable:1684)
#endif

  void BaseDynamicMem :: GetUsed (int nr, char * ch)
  {
    BaseDynamicMem * p = first;

    for (int i = 0; i < nr; i++)
      ch[i] = '0';

    while (p)
      {
        long unsigned hptr = (long unsigned) (p->ptr);
	// uintptr_t hptr = reinterpret_cast<uintptr_t>(p->ptr); //??

	hptr /= (1024*1024);
	hptr /= (4096/nr);

	size_t blocks = p->size / (1024*1024);
	blocks /= (4096/nr);
	
	// cout << "ptr = " << (void*)(p->ptr) << ", size = " << p->size << ", hptr = " << hptr << " blocks = " << blocks << endl;

	for (size_t i = 0; i <= blocks; i++)
	  ch[hptr+i] = '1';

	p = p->next;
      }
    
    {

      /*
    BaseMoveableMem * pm = BaseMoveableMem::first;
    while (pm)
      {
        long unsigned hptr = (long unsigned) pm->ptr;
        // uintptr_t hptr = reinterpret_cast<uintptr_t>(pm->ptr);

	hptr /= (1024*1024);
	hptr /= (4096/nr);

	size_t blocks = pm->size / (1024*1024);
	blocks /= (4096/nr);
	
	// cout << "moveable, ptr = " << (void*)(pm->ptr) << ", size = " << pm->size << ", hptr = " << hptr << " blocks = " << blocks << endl;

	for (size_t i = 0; i <= blocks; i++)
	  ch[hptr+i] = '1';

	pm = pm->next;
      }
      */
    }



  }

#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

}











namespace netgen
{
  //using namespace netgen;



  NgException :: NgException (const string & s) 
    : what(s)
  {
    ; 
  }


  NgException :: ~NgException () 
  {
    ;
  }

  /// append string to description
  void NgException :: Append (const string & s)
  { 
    what += s; 
  }

}


