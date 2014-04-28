#ifndef NGINTERFACE_V2
#define NGINTERFACE_V2


/**************************************************************************/
/* File:   nginterface_v2.hpp                                             */
/* Author: Joachim Schoeberl                                              */
/* Date:   May  09                                                        */
/**************************************************************************/

/*
  C++ interface to Netgen
*/

namespace netgen
{
  struct T_EDGE2
  {
    int orient:1;
    int nr:31;    // 0-based
  };
  struct T_FACE2
  {
    int orient:3;
    int nr:29;    // 0-based
  };

  class Ng_Element
  {

    class Ng_Points
    {
    public:
      int num;
      const int * ptr;
  
      int Size() const { return num; }
      int operator[] (int i) const { return ptr[i]-1; }
    };


    class Ng_Vertices
    {
    public:
      int num;
      const int * ptr;
  
      int Size() const { return num; }
      int operator[] (int i) const { return ptr[i]-1; }
    };

    class Ng_Edges
    {
    public:
      int num;
      const T_EDGE2 * ptr;
  
      int Size() const { return num; }
      // int operator[] (int i) const { return abs (ptr[i])-1; }
      int operator[] (int i) const { return ptr[i].nr; }
    };

    class Ng_Faces
    {
    public:
      int num;
      const T_FACE2 * ptr;
  
      int Size() const { return num; }
      // int operator[] (int i) const { return (ptr[i]-1) / 8; }
      int operator[] (int i) const { return ptr[i].nr; }
    };

  public:
    NG_ELEMENT_TYPE type;
    NG_ELEMENT_TYPE GetType() const { return type; }
    
    Ng_Points points;      // all points
    Ng_Vertices vertices;
    Ng_Edges edges;
    Ng_Faces faces;
  };

  
  class Ng_Point
  {
    double * pt;
  public:
    Ng_Point (double * apt) : pt(apt) { ; }
    double operator[] (int i)
    { return pt[i]; }
    operator const double * () { return pt; }
  };




  template <int DIM> class Ng_Node;

  template <>
  class Ng_Node<1>
  {
    class Ng_Vertices
    {
    public:
      const int * ptr;
  
      int Size() const { return 2; }
      int operator[] (int i) const { return ptr[i]-1; }
    };


  public:
    Ng_Vertices vertices;
  };



  template <>
  class Ng_Node<2>
  {
    class Ng_Vertices
    {
    public:
      int nv;
      const int * ptr;
  
      int Size() const { return nv; }
      int operator[] (int i) const { return ptr[i]-1; }
    };

    class Ng_Edges
    {
    public:
      int ned;
      const int * ptr;
  
      int Size() const { return ned; }
      int operator[] (int i) const { return ptr[i]-1; }
    };


  public:
    Ng_Vertices vertices;
    Ng_Edges edges;
  };



    






  class DLL_HEADER Ngx_Mesh
  {
  private:
    class Mesh * mesh;
    
  public:
    // Ngx_Mesh () { ; }
    // Ngx_Mesh(class Mesh * amesh) : mesh(amesh) { ; }
    Ngx_Mesh(class Mesh * amesh = NULL);
    void LoadMesh (const string & filename);

    void LoadMesh (istream & str);
    void SaveMesh (ostream & str) const;
    void DoArchive (ngstd::Archive & archive);

    virtual ~Ngx_Mesh();

    bool Valid () { return mesh != NULL; }
    
    int GetDimension() const;
    int GetNLevels() const;

    int GetNElements (int dim) const;
    int GetNNodes (int nt) const;

    Ng_Point GetPoint (int nr) const;

    template <int DIM> 
    Ng_Element GetElement (int nr) const;

    template <int DIM> 
    int GetElementIndex (int nr) const;


    /// Curved Elements:
    /// elnr .. element nr
    /// xi..... DIM_EL local coordinates
    /// x ..... DIM_SPACE global coordinates 
    /// dxdxi...DIM_SPACE x DIM_EL Jacobian matrix (row major storage)
    template <int DIM_EL, int DIM_SPACE> 
    void ElementTransformation (int elnr,
                                const double * xi, 
                                double * x, 
                                double * dxdxi) const;
    
    
    /// Curved Elements:
    /// elnr .. element nr
    /// npts .. number of points
    /// xi..... DIM_EL local coordinates
    /// sxi ... step xi
    /// x ..... DIM_SPACE global coordinates
    /// dxdxi...DIM_SPACE x DIM_EL Jacobian matrix (row major storage)
    template <int DIM_EL, int DIM_SPACE> 
    void MultiElementTransformation (int elnr, int npts,
                                     const double * xi, size_t sxi,
                                     double * x, size_t sx,
                                     double * dxdxi, size_t sdxdxi) const;
    

    template <int DIM>
    Ng_Node<DIM> GetNode (int nr) const;
    
    
    template <int DIM>
    int GetNNodes ();

    // Find element of point, returns local coordinates
    template <int DIM>
    int FindElementOfPoint 
    (double * p, double * lami,
     bool build_searchtrees = false, 
     int * const indices = NULL, int numind = 0) const;
    
  };



  DLL_HEADER Ngx_Mesh * LoadMesh (const string & filename);
}


#ifdef HAVE_NETGEN_SOURCES
#include <meshing.hpp>

namespace netgen
{
#define NGX_INLINE inline
#include <nginterface_v2_impl.hpp>
}

#endif


#endif

