#ifdef NG_PYTHON

#include <boost/python.hpp>
#include <meshing.hpp>
#include <geometry2d.hpp>

using namespace netgen;
namespace bp = boost::python;


//////////////////////////////////////////////////////////////////////
// Lambda to function pointer conversion
template <typename Function>
struct function_traits
  : public function_traits<decltype(&Function::operator())> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...) const> {
  typedef ReturnType (*pointer)(Args...);
  typedef ReturnType return_type;
};

template <typename Function>
typename function_traits<Function>::pointer
FunctionPointer (const Function& lambda) {
  return static_cast<typename function_traits<Function>::pointer>(lambda);
}


template <class T>
inline string ToString (const T& t)
{
  stringstream ss;
  ss << t;
  return ss.str();
}



void ExportGeom2d() 
{
  
  std::string nested_name = "geom2d";
  if( bp::scope() )
    nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".geom2d");
                                           
  bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
  
  cout << "exporting geom2d " << nested_name << endl;
  bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
  parent.attr("geom2d") = module ;
  
  bp::scope local_scope(module);

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
	.def("AppendSegment", FunctionPointer([](SplineGeometry2d &self, int point_index1, int point_index2)//, int leftdomain, int rightdomain)
	  {
		  LineSeg<2> * l = new LineSeg<2>(self.GetPoint(point_index1), self.GetPoint(point_index2));
		  SplineSegExt * seg = new SplineSegExt(*l);
		  seg->leftdom = 1;// leftdomain;
		  seg->rightdom = 0;// rightdomain;
		  seg->hmax = 1e99;
		  seg->reffak = 1;
		  seg->copyfrom = -1;

		  self.AppendSegment(seg);
	  }))//, (bp::arg("self"), bp::arg("point_index1"), bp::arg("point_index2"), bp::arg("leftdomain") = 1, bp::arg("rightdomain") = 0) )
	.def("AppendSegment", FunctionPointer([](SplineGeometry2d &self, int point_index1, int point_index2, int point_index3)//, int leftdomain, int rightdomain)
	  {
		  SplineSeg3<2> * seg3 = new SplineSeg3<2>(self.GetPoint(point_index1), self.GetPoint(point_index2), self.GetPoint(point_index3));
		  SplineSegExt * seg = new SplineSegExt(*seg3);
		  seg->leftdom = 1;// leftdomain;
		  seg->rightdom = 0;// rightdomain;
		  seg->hmax = 1e99;
		  seg->reffak = 1;
		  seg->copyfrom = -1;

		  self.AppendSegment(seg);
	  }))//, (bp::arg("self"), bp::arg("point_index1"), bp::arg("point_index2"), bp::arg("point_index3"), bp::arg("leftdomain") = 1, bp::arg("rightdomain") = 0 ) )
	.def("PlotData", FunctionPointer([](SplineGeometry2d &self)
	  {
		  Box<2> box(self.GetBoundingBox());
		  double xdist = box.PMax()(0) - box.PMin()(0);
		  double ydist = box.PMax()(1) - box.PMin()(1);
		  bp::tuple xlim = bp::make_tuple(box.PMin()(0) - 0.1*xdist, box.PMax()(0) + 0.1*xdist);
		  bp::tuple ylim = bp::make_tuple(box.PMin()(1) - 0.1*ydist, box.PMax()(1) + 0.1*ydist);

		  bp::list xpoints, ypoints;
		  GeomPoint<2> point = self.splines[0]->StartPI();
		  xpoints.append(point(0));
		  ypoints.append(point(1));

		  for (int i = 0; i < self.splines.Size(); i++)
		  {
			  if (self.splines[i]->GetType().compare("line")==0)
			  {
				  GeomPoint<2> point = self.splines[i]->EndPI();
				  xpoints.append(point(0));
				  ypoints.append(point(1));
			  }
			  else if (self.splines[i]->GetType().compare("spline3")==0)
			  {
				  double len = self.splines[i]->Length();
				  int n = floor(len/(0.05*min(xdist,ydist)));
				  //cout << n << endl;
				  for (int j = 1; j <= n; j++)
				  {
					  GeomPoint<2> point = self.splines[i]->GetPoint(j*1./n);
					  xpoints.append(point(0));
					  ypoints.append(point(1));
				  }
			  }
			  else
			  {
				  cout << "spline is neither line nor spline3" << endl;
			  }
				  
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
		  Mesh * mesh = NULL;
		  self.GenerateMesh(mesh, mparam, 0, 0);
		  return shared_ptr<Mesh>(mesh);
	  }))
	  
	  ;
  
}


#endif

