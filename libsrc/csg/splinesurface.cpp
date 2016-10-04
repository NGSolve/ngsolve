
#include <csg.hpp>

namespace netgen
{
SplineSurface :: SplineSurface() : OneSurfacePrimitive()
{ ; }

SplineSurface :: ~SplineSurface() { ; }

void SplineSurface :: AppendPoint(const Point<3> & p, const double reffac, const bool hpref)
{
  geompoints.Append(GeomPoint<3>(p,reffac));
  geompoints.Last().hpref = hpref;
}
  
  void SplineSurface :: AppendSegment(SplineSeg<3>* spline, string* bcname)
  {
      splines.Append(spline);
      bcnames.Append(bcname);
  }
    
double SplineSurface :: CalcFunctionValue (const Point<3> & point) const
{
  auto v1 = splines[0]->GetTangent(0);
  auto v2 = splines.Last()->GetTangent(1);
  auto p1 = splines[0]->GetPoint(0);
  auto n = Cross(v1,v2);
  n.Normalize();
  return n * (point-p1);
}

void SplineSurface :: CalcGradient (const Point<3> & point, Vec<3> & grad) const
{
  auto v1 = splines[0]->GetTangent(0);
  auto v2 = splines.Last()->GetTangent(1);
  auto p1 = splines[0]->GetPoint(0);
  grad = Cross(v1,v2);
  grad.Normalize();
}

  double SplineSurface :: HesseNorm() const
{
  return 0;
}

  Point<3> SplineSurface :: GetSurfacePoint() const
{
  return splines[0]->GetPoint(0);
}

  
  void SplineSurface :: CalcSpecialPoints(Array<Point<3>> & pts) const
  {
    for(auto pt : geompoints)
      {
	pts.Append(Point<3>(pt));
      }
  }

  Array<Plane*>* SplineSurface :: CreatePlanes() const
  {
    auto planes = new Array<Plane*>();
    auto sol = new Solid(new Plane(splines[0]->GetPoint(0),-(splines[0]->GetTangent(0))));
    for(auto spline : splines)
      {
	planes->Append(new Plane(spline->GetPoint(0),-spline->GetTangent(0)));
      }
    return planes;
  }


  void SplineSurface :: Print(ostream & str) const
{
  str << "SplineSurface " << endl;
}

INSOLID_TYPE SplineSurface :: BoxInSolid(const BoxSphere<3> & box) const
{
  
  int i;
  double val;
   
  val = CalcFunctionValue (box.Center());
  if (val > box.Diam() / 2) return IS_OUTSIDE;
  if (val < -box.Diam() / 2) return IS_INSIDE;
  
  Vec<3> n = GetNormalVector(box.Center());
  double cx = n[0];
  double cy = n[1];
  double cz = n[2];
  
  if (val > 0)
    {
      Vec<3> vdiag = box.PMax() - box.PMin();
      double modify = (vdiag(0) * fabs (cx) + 
		       vdiag(1) * fabs (cy) + 
		       vdiag(2) * fabs (cz) ) / 2;
      
      if (val - modify < 0)
	return DOES_INTERSECT;
      return IS_OUTSIDE;
    }
  else
    {
      Vec<3> vdiag = box.PMax() - box.PMin();
      double modify =  (vdiag(0) * fabs (cx) + 
			vdiag(1) * fabs (cy) + 
			vdiag(2) * fabs (cz) ) / 2;
      if (val + modify > 0)
	return DOES_INTERSECT;
      return IS_INSIDE;
    }
}
}
