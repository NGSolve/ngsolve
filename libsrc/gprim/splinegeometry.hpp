/*


JS, Nov 2007


The 2D/3D template-base classes should go into the libsrc/gprim directory

in geom2d only 2D - Geometry classes (with material properties etc.)


*/

#include "spline.hpp"


#ifndef _FILE_SPLINEGEOMETRY
#define _FILE_SPLINEGEOMETRY

namespace netgen
{


  template < int D >
  class SplineGeometry 
  {
    // protected:
  public:  
    Array < GeomPoint<D> > geompoints;
    Array < SplineSeg<D>* > splines;

    DLL_HEADER ~SplineGeometry();

    DLL_HEADER int Load (const Array<double> & raw_data, const int startpos = 0);

    DLL_HEADER void GetRawData (Array<double> & raw_data) const;


    const Array<SplineSeg<D>*> & GetSplines () const
    { return splines; }

    int GetNSplines (void) const { return splines.Size(); }
    string GetSplineType (const int i) const { return splines[i]->GetType(); }
    SplineSeg<D> & GetSpline (const int i) {return *splines[i];}
    const SplineSeg<D> & GetSpline (const int i) const {return *splines[i];}

    DLL_HEADER void GetBoundingBox (Box<D> & box) const;
    Box<D> GetBoundingBox () const 
    { Box<D> box; GetBoundingBox (box); return box; }

    int GetNP () const { return geompoints.Size(); }
    const GeomPoint<D> & GetPoint(int i) const { return geompoints[i]; }

    // void SetGrading (const double grading);
    DLL_HEADER void AppendPoint (const Point<D> & p, const double reffac = 1., const bool hpref = false);


    void AppendSegment(SplineSeg<D> * spline)
    {
      splines.Append (spline);
    }
  };



}

#endif // _FILE_SPLINEGEOMETRY
