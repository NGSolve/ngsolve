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
      const int * ptr;
  
      int Size() const { return num; }
      int operator[] (int i) const { return abs (ptr[i])-1; }
    };

    class Ng_Faces
    {
    public:
      int num;
      const int * ptr;
  
      int Size() const { return num; }
      int operator[] (int i) const { return (ptr[i]-1) / 8; }
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
  public:
    double * pt;
    double operator[] (int i)
    { return pt[i]; }
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
    Ngx_Mesh(class Mesh * amesh);
    virtual ~Ngx_Mesh();
    
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
    Ng_Node<DIM> GetNode (int nr);
    
    
    template <int DIM>
    int GetNNodes ();

    // Find element of point, returns local coordinates
    template <int DIM>
    int FindElementOfPoint 
    (double * p, double * lami,
     bool build_searchtrees = false, 
     int * const indices = NULL, int numind = 0);
    
  };



  DLL_HEADER Ngx_Mesh * LoadMesh (const string & filename);



}
#endif

