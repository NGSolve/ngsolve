#ifndef CURVEDELEMS
#define CURVEDELEMS

/**************************************************************************/
/* File:   curvedelems.hpp                                                */
/* Author: Robert Gaisbauer (first version)                               */
/*         redesign by Joachim Schoeberl                                  */
/* Date:   27. Sep. 02, Feb 2006                                          */
/**************************************************************************/




class Refinement;


class CurvedElements
{
  const Mesh & mesh;

  Array<int> edgeorder;
  Array<int> faceorder;

  Array<int> edgecoeffsindex;
  Array<int> facecoeffsindex;

  Array< Vec<3> > edgecoeffs;
  Array< Vec<3> > facecoeffs;

  Array< double > edgeweight;  // for rational 2nd order splines

  int order;
  bool rational;

  bool ishighorder;

public:
  DLL_HEADER CurvedElements (const Mesh & amesh);
  DLL_HEADER ~CurvedElements();

  // bool IsHighOrder() const { return order > 1; }
  bool IsHighOrder() const { return ishighorder; }

  // void SetHighOrder (int aorder) { order=aorder; }
  void SetIsHighOrder (bool ho) { ishighorder = ho; }
  
  DLL_HEADER void BuildCurvedElements(const Refinement * ref, int aorder, bool arational = false);

  int GetOrder () { return order; }


  bool IsSegmentCurved (SegmentIndex segnr) const;
  bool IsSurfaceElementCurved (SurfaceElementIndex sei) const;
  bool IsElementCurved (ElementIndex ei) const;
  bool IsElementHighOrder (ElementIndex ei) const;


  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> & x)
  { CalcSegmentTransformation (xi, segnr, &x, NULL); };

  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Vec<3> & dxdxi)
  { CalcSegmentTransformation (xi, segnr, NULL, &dxdxi); };

  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> & x, Vec<3> & dxdxi)
  { CalcSegmentTransformation (xi, segnr, &x, &dxdxi, NULL); };

  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> & x, Vec<3> & dxdxi, bool & curved)
  { CalcSegmentTransformation (xi, segnr, &x, &dxdxi, &curved); };



  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Point<3> & x)
  { CalcSurfaceTransformation (xi, elnr, &x, NULL); };

  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Mat<3,2> & dxdxi)
  { CalcSurfaceTransformation (xi, elnr, NULL, &dxdxi); };

  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Point<3> & x, Mat<3,2> & dxdxi)
  { CalcSurfaceTransformation (xi, elnr, &x, &dxdxi, NULL); };

  void CalcSurfaceTransformation (const Point<2> & xi, SurfaceElementIndex elnr,
				  Point<3> & x, Mat<3,2> & dxdxi, bool & curved)
  { CalcSurfaceTransformation (xi, elnr, &x, &dxdxi, &curved); };





  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Point<3> & x)
  { CalcElementTransformation (xi, elnr, &x, NULL); };

  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Mat<3,3> & dxdxi)
  { CalcElementTransformation (xi, elnr, NULL, &dxdxi); };

  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Point<3> & x, Mat<3,3> & dxdxi)
  { CalcElementTransformation (xi, elnr, &x, &dxdxi /* , NULL */ ); };

  void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
				  Point<3> & x, Mat<3,3> & dxdxi,
                                  void * buffer, bool valid)
  { CalcElementTransformation (xi, elnr, &x, &dxdxi, /* NULL, */ buffer, valid ); };

  // void CalcElementTransformation (const Point<3> & xi, ElementIndex elnr,
  // 				  Point<3> & x, Mat<3,3> & dxdxi) // , bool & curved)
  //   { CalcElementTransformation (xi, elnr, &x, &dxdxi /* , &curved * ); }



  void CalcMultiPointSegmentTransformation (Array<double> * xi, SegmentIndex segnr,
					    Array<Point<3> > * x,
					    Array<Vec<3> > * dxdxi);

  template <int DIM_SPACE>
  void CalcMultiPointSegmentTransformation (SegmentIndex elnr, int n,
                                            const double * xi, size_t sxi,
                                            double * x, size_t sx,
                                            double * dxdxi, size_t sdxdxi);

  void CalcMultiPointSurfaceTransformation (Array< Point<2> > * xi, SurfaceElementIndex elnr,
					    Array< Point<3> > * x,
					    Array< Mat<3,2> > * dxdxi);

  template <int DIM_SPACE, typename T>
  void CalcMultiPointSurfaceTransformation (SurfaceElementIndex elnr, int n,
                                            const T * xi, size_t sxi,
                                            T * x, size_t sx,
                                            T * dxdxi, size_t sdxdxi);

  void CalcMultiPointElementTransformation (Array< Point<3> > * xi, ElementIndex elnr,
					    Array< Point<3> > * x,
					    Array< Mat<3,3> > * dxdxi);

  template <typename T>
  void CalcMultiPointElementTransformation (ElementIndex elnr, int n,
                                            const T * xi, size_t sxi,
                                            T * x, size_t sx,
                                            T * dxdxi, size_t sdxdxi);




private:
  
  void CalcSegmentTransformation (double xi, SegmentIndex segnr,
				  Point<3> * x = NULL, Vec<3> * dxdxi = NULL, bool * curved = NULL);

  void CalcSurfaceTransformation (Point<2> xi, SurfaceElementIndex elnr,
				  Point<3> * x = NULL, Mat<3,2> * dxdxi = NULL, bool * curved = NULL);

  void CalcElementTransformation (Point<3> xi, ElementIndex elnr,
				  Point<3> * x = NULL, Mat<3,3> * dxdxi = NULL, // bool * curved = NULL,
                                  void * buffer = NULL, bool valid = 0);






  class SegmentInfo
  {
  public:
    SegmentIndex elnr;
    int order;
    int nv;
    int ndof;
    int edgenr;
  };

  void CalcElementShapes (SegmentInfo &  elnr, double xi, Vector & shapes) const;
  void GetCoefficients (SegmentInfo & elnr, Array<Vec<3> > & coefs) const;
  void CalcElementDShapes (SegmentInfo & elnr, double xi, Vector & dshapes) const;


  class ElementInfo
  {
  public:
    ElementIndex elnr;
    int order;
    int nv;
    int ndof;
    int nedges;
    int nfaces;
    int edgenrs[12];
    int facenrs[6];
    Mat<3> hdxdxi;
    Vec<3> hcoefs[10]; // enough for second order tets
  };

  template <typename T>
  void CalcElementShapes (ElementInfo & info, Point<3,T> xi, TFlatVector<T> shapes) const;
  void GetCoefficients (ElementInfo & info, Vec<3> * coefs) const;
  template <typename T>  
  void CalcElementDShapes (ElementInfo & info, const Point<3,T> xi, MatrixFixWidth<3,T> dshapes) const;

  template <typename T>
  bool EvaluateMapping (ElementInfo & info, const Point<3,T> xi, Point<3,T> & x, Mat<3,3,T> & jac) const;  
  
  class SurfaceElementInfo
  {
  public:
    SurfaceElementIndex elnr;
    int order;
    int nv;
    int ndof;
    ArrayMem<int,4> edgenrs;
    int facenr;
  };

  template <typename T>
  void CalcElementShapes (SurfaceElementInfo & elinfo, const Point<2,T> xi, TFlatVector<T> shapes) const;
  template <int DIM_SPACE>
  void GetCoefficients (SurfaceElementInfo & elinfo, Array<Vec<DIM_SPACE> > & coefs) const;
  template <typename T>
  void CalcElementDShapes (SurfaceElementInfo & elinfo, const Point<2,T> xi, MatrixFixWidth<2,T> dshapes) const;

  template <int DIM_SPACE, typename T>
  bool EvaluateMapping (SurfaceElementInfo & info, const Point<2,T> xi, Point<DIM_SPACE,T> & x, Mat<DIM_SPACE,2,T> & jac) const;  
};



#endif
