/*

2d Spline curve for Mesh generator

*/


#include <mystdlib.h>
#include <linalg.hpp>
#include <gprim.hpp>
#include "splinegeometry.hpp"
 
namespace netgen
{


  template<int D>
  SplineGeometry<D> :: ~SplineGeometry()
  {
    for(int i = 0; i < splines.Size(); i++)
      delete splines[i];
  }


  template<int D>
  void SplineGeometry<D> :: GetRawData (Array<double> & raw_data) const
  {
    raw_data.Append(D);
    // raw_data.Append(elto0);

    raw_data.Append(splines.Size());
    for(int i=0; i<splines.Size(); i++)
      splines[i]->GetRawData(raw_data);
  }



  template<int D>
  int SplineGeometry<D> :: Load (const Array<double> & raw_data, const int startpos)
  {
    int pos = startpos;
    if(raw_data[pos] != D)
      throw NgException("wrong dimension of spline raw_data");

    pos++;

    // elto0 = raw_data[pos]; pos++;

    splines.SetSize(int(raw_data[pos]));
    pos++;

    Array< Point<D> > pts(3);

    for(int i=0; i<splines.Size(); i++)
      {
	int type = int(raw_data[pos]);
	pos++;
      
	for(int j=0; j<type; j++)
	  for(int k=0; k<D; k++)
	    {
	      pts[j](k) = raw_data[pos];
	      pos++;
	    }

	if (type == 2)
	  {
	    splines[i] = new LineSeg<D>(GeomPoint<D>(pts[0],1),
					GeomPoint<D>(pts[1],1));
	  }
	else if (type == 3)
	  {
	    splines[i] = new SplineSeg3<D>(GeomPoint<D>(pts[0],1),
					   GeomPoint<D>(pts[1],1),
					   GeomPoint<D>(pts[2],1));
	  }
	else
	  throw NgException("something wrong with spline raw data");

      }
    return pos;
  }







  

  template<int D>
  void SplineGeometry<D> :: GetBoundingBox (Box<D> & box) const
  {
    if (!splines.Size())
      {
	Point<D> auxp = 0.;
	box.Set (auxp);
	return;
      }

    Array<Point<D> > points;
    for (int i = 0; i < splines.Size(); i++)
      {
	splines[i]->GetPoints (20, points);

	if (i == 0) box.Set(points[0]);
	for (int j = 0; j < points.Size(); j++)
	  box.Add (points[j]);
      }
  }

  /*
  template<int D>
  void SplineGeometry<D> :: SetGrading (const double grading)
  { 
    elto0 = grading;
  }
  */

  template<int D>
  void SplineGeometry<D> :: AppendPoint (const Point<D> & p, const double reffac, const bool hpref)
  {
    geompoints.Append (GeomPoint<D>(p, reffac));
    geompoints.Last().hpref = hpref;
  }



  template class SplineGeometry<2>;
  template class SplineGeometry<3>;
}


