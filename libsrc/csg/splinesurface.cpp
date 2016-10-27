
#include <csg.hpp>

namespace netgen
{
void SplineSurface :: AppendPoint(const Point<3> & p, const double reffac, const bool hpref)
{
  auto pp = Point<3>(p);
  geompoints.Append(GeomPoint<3>(pp,reffac));
  geompoints.Last().hpref = hpref;
}
  
  void SplineSurface :: AppendSegment(SplineSeg<3>* spline, string* bcname, double amaxh)
  {
    splines.Append(spline);
    bcnames.Append(bcname);
    maxh.Append(amaxh);
  }

  string* SplineSurface :: GetBCNameOf (Point<3> p1, Point<3> p2) const
  {
    
    double eps = 1e-5;
    for(int i=0; i<splines.Size(); i++)
      {
	auto pp1 = Point<3>(splines[i]->GetPoint(0));
	Project(pp1);
	auto pp2 = Point<3>(splines[i]->GetPoint(1));
	Project(pp2);
	if (((pp1-p1).Length()<eps && (pp2-p2).Length() < eps) || ((pp1-p2).Length() < eps && (pp2-p1).Length() < eps))
	  {
	    return bcnames[i];
	  }
      }
    return new string("default");
  }

  Array<OneSurfacePrimitive*>* SplineSurface :: CreateCuttingSurfaces() const
  {
    auto cuttings = new Array<OneSurfacePrimitive*>();
    for (auto cut : *cuts)
      cuttings->Append(cut);
    for(int i = 0; i<splines.Size(); i++)
      {
	auto spline = splines[i];
	auto lineseg = dynamic_cast<LineSeg<3>*>(spline);
	auto p1 = Point<3>(spline->GetPoint(0));
	Project(p1);
	auto p2 = Point<3>(spline->GetPoint(1));
	Project(p2);
	auto vec = Vec<3>(p2)-Vec<3>(p1);
	auto plane = new Plane(p1,-Cross(vec,baseprimitive->GetNormalVector(p1)));
	if(maxh[i]>0)
	  {
	  plane->SetMaxH(maxh[i]);
	  }
	cuttings->Append(plane);
      }
    return cuttings;
  }
  
  void SplineSurface :: Print(ostream & str) const
{
  str << "SplineSurface " << endl;
}

}
