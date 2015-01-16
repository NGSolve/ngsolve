#ifdef NG_PYTHON

#include <boost/python.hpp>
#include <../general/ngpython.hpp>

#include <meshing.hpp>
#include <geometry2d.hpp>

using namespace netgen;
namespace bp = boost::python;


DLL_HEADER void ExportGeom2d() 
{
  ModuleScope module("geom2d");

  bp::class_<SplineGeometry2d, boost::noncopyable>("SplineGeometry")
	.def("Load",&SplineGeometry2d::Load)
	.def("AppendPoint", FunctionPointer([](SplineGeometry2d &self, double px, double py)
	  {
		  Point<2> p;
		  p(0) = px;
		  p(1) = py;
		  self.geompoints.Append(GeomPoint<2>(p,1));
		  return self.geompoints.Size()-1;
	  }))
    .def("Append", FunctionPointer([](SplineGeometry2d &self, bp::list segment, int leftdomain, int rightdomain, int bc)
	  {
		  bp::extract<std::string> segtype(segment[0]);

		  SplineSegExt * seg;
		  if (segtype().compare("line") == 0)
		  {
			  bp::extract<int> point_index1(segment[1]);
			  bp::extract<int> point_index2(segment[2]);
			  //point_index1.check()

			  LineSeg<2> * l = new LineSeg<2>(self.GetPoint(point_index1()), self.GetPoint(point_index2()));
			  seg = new SplineSegExt(*l);
		  }
		  else if (segtype().compare("spline3") == 0)
		  {
			  bp::extract<int> point_index1(segment[1]);
			  bp::extract<int> point_index2(segment[2]);
			  bp::extract<int> point_index3(segment[3]);

			  SplineSeg3<2> * seg3 = new SplineSeg3<2>(self.GetPoint(point_index1()), self.GetPoint(point_index2()), self.GetPoint(point_index3()));
			  seg = new SplineSegExt(*seg3);
		  }
		  else
		  {
			  cout << "Appended segment is not a line or a spline3" << endl;
		  }
		  seg->leftdom = leftdomain;
		  seg->rightdom = rightdomain;
		  seg->hmax = 1e99;
		  seg->reffak = 1;
		  seg->copyfrom = -1;
		  seg->bc = (bc >= 0) ? bc : self.GetNSplines();
		  self.AppendSegment(seg);
	  }), (bp::arg("self"), bp::arg("point_indices"), bp::arg("leftdomain") = 1, bp::arg("rightdomain") = 0, bp::arg("bc")=-1))
	.def("AppendSegment", FunctionPointer([](SplineGeometry2d &self, bp::list point_indices, int leftdomain, int rightdomain)
		{
		  int npts = bp::len(point_indices);
		  SplineSegExt * seg;
		  //int a = bp::extract<int>(point_indices[0]);
		  if (npts == 2)
		  {
			  LineSeg<2> * l = new LineSeg<2>(self.GetPoint(bp::extract<int>(point_indices[0])), self.GetPoint(bp::extract<int>(point_indices[1])));
			  seg = new SplineSegExt(*l);
			  
		  }
		  else if (npts == 3)
		  {
			  SplineSeg3<2> * seg3 = new SplineSeg3<2>(self.GetPoint(bp::extract<int>(point_indices[0])), self.GetPoint(bp::extract<int>(point_indices[1])), self.GetPoint(bp::extract<int>(point_indices[2])));
			  seg = new SplineSegExt(*seg3);

		  }
		  seg->leftdom = leftdomain;
		  seg->rightdom = rightdomain;
		  seg->hmax = 1e99;
		  seg->reffak = 1;
		  seg->copyfrom = -1;
		  self.AppendSegment(seg);
		}), (bp::arg("self"), bp::arg("point_indices"), bp::arg("leftdomain") = 1, bp::arg("rightdomain") = 0) )
	//.def("AppendSegment", FunctionPointer([](SplineGeometry2d &self, int point_index1, int point_index2)//, int leftdomain, int rightdomain)
	//  {
	//	  LineSeg<2> * l = new LineSeg<2>(self.GetPoint(point_index1), self.GetPoint(point_index2));
	//	  SplineSegExt * seg = new SplineSegExt(*l);
	//	  seg->leftdom = 1;// leftdomain;
	//	  seg->rightdom = 0;// rightdomain;
	//	  seg->hmax = 1e99;
	//	  seg->reffak = 1;
	//	  seg->copyfrom = -1;

	//	  self.AppendSegment(seg);
	//  }))//, (bp::arg("self"), bp::arg("point_index1"), bp::arg("point_index2"), bp::arg("leftdomain") = 1, bp::arg("rightdomain") = 0) )
	//.def("AppendSegment", FunctionPointer([](SplineGeometry2d &self, int point_index1, int point_index2, int point_index3)//, int leftdomain, int rightdomain)
	//  {
	//	  SplineSeg3<2> * seg3 = new SplineSeg3<2>(self.GetPoint(point_index1), self.GetPoint(point_index2), self.GetPoint(point_index3));
	//	  SplineSegExt * seg = new SplineSegExt(*seg3);
	//	  seg->leftdom = 1;// leftdomain;
	//	  seg->rightdom = 0;// rightdomain;
	//	  seg->hmax = 1e99;
	//	  seg->reffak = 1;
	//	  seg->copyfrom = -1;
	//	  self.AppendSegment(seg);
	//  }))//, (bp::arg("self"), bp::arg("point_index1"), bp::arg("point_index2"), bp::arg("point_index3"), bp::arg("leftdomain") = 1, bp::arg("rightdomain") = 0 ) )
	.def("PlotData", FunctionPointer([](SplineGeometry2d &self)
	  {
		  Box<2> box(self.GetBoundingBox());
		  double xdist = box.PMax()(0) - box.PMin()(0);
		  double ydist = box.PMax()(1) - box.PMin()(1);
		  bp::tuple xlim = bp::make_tuple(box.PMin()(0) - 0.1*xdist, box.PMax()(0) + 0.1*xdist);
		  bp::tuple ylim = bp::make_tuple(box.PMin()(1) - 0.1*ydist, box.PMax()(1) + 0.1*ydist);

		  bp::list xpoints, ypoints;

		  for (int i = 0; i < self.splines.Size(); i++)
		  {
			  bp::list xp, yp;
			  if (self.splines[i]->GetType().compare("line")==0)
			  {
				  GeomPoint<2> p1 = self.splines[i]->StartPI();
				  GeomPoint<2> p2 = self.splines[i]->EndPI();
				  xp.append(p1(0));
				  xp.append(p2(0));
				  yp.append(p1(1));
				  yp.append(p2(1));
			  }
			  else if (self.splines[i]->GetType().compare("spline3")==0)
			  {
				  double len = self.splines[i]->Length();
				  int n = floor(len/(0.05*min(xdist,ydist)));
				  
				  for (int j = 0; j <= n; j++)
				  {
					  GeomPoint<2> point = self.splines[i]->GetPoint(j*1./n);
					  xp.append(point(0));
					  yp.append(point(1));
				  }
			  }
			  else
			  {
				  cout << "spline is neither line nor spline3" << endl;
			  }
			  xpoints.append(xp);
			  ypoints.append(yp);
				  
		  }
		  return bp::tuple(bp::make_tuple(xlim, ylim, xpoints, ypoints));

	  }))
	.def("PointData", FunctionPointer([](SplineGeometry2d &self)
	  {
		  bp::list xpoints, ypoints, pointindex;
		  
		  for (int i = 0; i < self.geompoints.Size(); i++)
		  {
			  pointindex.append(i);
			  xpoints.append(self.geompoints[i][0]);
			  ypoints.append(self.geompoints[i][1]);
		  }
		  return bp::tuple(bp::make_tuple(xpoints, ypoints, pointindex));
		  
	  }))
	.def("SegmentData", FunctionPointer([](SplineGeometry2d &self)
	  {
		  bp::list leftpoints, rightpoints, leftdom, rightdom;

		  for (int i = 0; i < self.splines.Size(); i++)
		  {
			  GeomPoint<2> point = self.splines[i]->GetPoint(0.5);
			  Vec<2> normal = self.GetSpline(i).GetTangent(0.5);
			  double temp = normal(0);
			  normal(0) = normal(1);
			  normal(1) = -temp;

			  leftdom.append(self.GetSpline(i).leftdom);
			  rightdom.append(self.GetSpline(i).rightdom);

			  rightpoints.append(bp::make_tuple(point(0), point(1), normal(0)<0, normal(1)<0));
			  leftpoints.append(bp::make_tuple(point(0), point(1), normal(0)<0, normal(1)<0));
		  }
		  return bp::tuple(bp::make_tuple(leftpoints, rightpoints, leftdom, rightdom));

	  }))
	.def("Print", FunctionPointer([](SplineGeometry2d &self)
	  {
		  for (int i = 0; i < self.geompoints.Size(); i++)
		  {
			  cout << i << " : " << self.geompoints[i][0] << " , " << self.geompoints[i][1] << endl;
		  }
		  //Box<2> box(self.GetBoundingBox());
		  //cout << box.PMin() << endl;
		  //cout << box.PMax() << endl;
		  cout << self.splines.Size() << endl;
		  for (int i = 0; i < self.splines.Size(); i++)
		  {
			  cout << self.splines[i]->GetType() << endl;
			  //cout << i << " : " << self.splines[i]->GetPoint(0.1) << " , " << self.splines[i]->GetPoint(0.5) << endl;
		  }
	  }))
	  .def("GenerateMesh", FunctionPointer([](SplineGeometry2d &self, MeshingParameters & mparam)
		{
		  shared_ptr<Mesh> mesh;
		  self.GenerateMesh(mesh, mparam, 0, 0);
		  return mesh;
	  }))
	  
	  ;
  
}

BOOST_PYTHON_MODULE(libgeom2d) {
	ExportGeom2d();
}

#endif

