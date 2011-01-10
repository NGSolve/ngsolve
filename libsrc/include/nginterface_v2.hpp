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




  template <int DIM> 
  DLL_HEADER int Ng_GetNElements ();

  template <int DIM> 
  DLL_HEADER Ng_Element Ng_GetElement (int nr);


  
  class Ng_Point
  {
  public:
    double * pt;
    double operator[] (int i)
    { return pt[i]; }
  };

  DLL_HEADER Ng_Point Ng_GetPoint (int nr);




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



    
  template <int DIM>
  DLL_HEADER Ng_Node<DIM> Ng_GetNode (int nr);
  

  template <int DIM>
  DLL_HEADER int Ng_GetNNodes ();





  /// Curved Elements:
  /// xi..... DIM_EL local coordinates
  /// sxi ... step xi
  /// x ..... DIM_SPACE global coordinates
  /// dxdxi...DIM_SPACE x DIM_EL Jacobian matrix (row major storage)
  template <int DIM_EL, int DIM_SPACE> 
  DLL_HEADER void Ng_MultiElementTransformation (int elnr, int npts,
                                                 const double * xi, size_t sxi,
                                                 double * x, size_t sx,
                                                 double * dxdxi, size_t sdxdxi);
  
  template <int DIM> 
  DLL_HEADER int Ng_GetElementIndex (int nr);
}
#endif

