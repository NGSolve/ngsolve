#ifdef NG_PYTHON

#include <boost/python.hpp>
#include <boost/python/slice.hpp>

#include <mystdlib.h>
#include "meshing.hpp"


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



template <typename T, int BASE = 0, typename TIND = int>
void ExportArray ()
{
  string name = string("Array_") + typeid(T).name();
  bp::class_<Array<T,BASE,TIND>,boost::noncopyable>(name.c_str())
    .def ("__len__", &Array<T,BASE,TIND>::Size)
    .def ("__getitem__", 
          FunctionPointer ([](Array<T,BASE,TIND> & self, TIND i) -> T&
                           {
                             if (i < BASE || i >= BASE+self.Size())
                               bp::exec("raise IndexError()\n");
                             return self[i];
                           }),
          bp::return_value_policy<bp::reference_existing_object>())

    .def ("__iter__", 
          bp::range (FunctionPointer([](Array<T,BASE,TIND> & self) { return &self[BASE]; }),
                     FunctionPointer([](Array<T,BASE,TIND> & self) { return &self[BASE+self.Size()]; })))

    ;
}



void ExportNetgenMeshing() 
{
  
  std::string nested_name = "meshing";
  if( bp::scope() )
    nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".meshing");
                                           
  bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
  
  cout << "exporting meshing " << nested_name << endl;
  bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
  parent.attr("meshing") = module ;
  
  bp::scope local_scope(module);

  bp::class_<PointIndex>("PointId", bp::init<int>())
    .def("__repr__", &ToString<PointIndex>)
    .def("__str__", &ToString<PointIndex>)
    .add_property("nr", &PointIndex::operator int)
    ;
  
  bp::class_<Point<3>> ("Point")
    .def(bp::init<double,double,double>())
    ;
  
  bp::class_<MeshPoint,bp::bases<Point<3>>>("MeshPoint")
    .def(bp::init<Point<3>>())
    .add_property("p", FunctionPointer([](const MeshPoint & self)
                                       {
                                         bp::list l;
                                         l.append ( (self)[0] );
                                         l.append ( (self)[1] );
                                         l.append ( (self)[2] );
                                         return l;
                                       }))
    ;

  bp::class_<Element>("Element3D")
    .add_property("index", &Element::GetIndex, &Element::SetIndex)
    .add_property("vertices", 
                  FunctionPointer ([](const Element & self) -> bp::list
                                   {
                                     bp::list li;
                                     for (int i = 0; i < self.GetNV(); i++)
                                       li.append (self[i]);
                                     return li;
                                   }))
    ;
  
  ExportArray<Element>();
  ExportArray<Element2d>();
  ExportArray<MeshPoint,PointIndex::BASE,PointIndex>();
  ;
  
  
  bp::class_<Mesh,shared_ptr<Mesh>,boost::noncopyable>("Mesh")
    .def("__str__", &ToString<Mesh>)
    .def("Load", static_cast<void(Mesh::*)(const string & name)>(&Mesh::Load))
    .def("Save", static_cast<void(Mesh::*)(const string & name)const>(&Mesh::Save))

    .def("Elements3D", 
         static_cast<Array<Element>&(Mesh::*)()> (& &Mesh::VolumeElements),
         bp::return_value_policy<bp::reference_existing_object>())

    .def("Elements2D", 
         static_cast<Array<Element2d>&(Mesh::*)()> (& &Mesh::SurfaceElements),
         bp::return_value_policy<bp::reference_existing_object>())

    .def("Points", 
         static_cast<Mesh::T_POINTS&(Mesh::*)()> (& &Mesh::Points),
         bp::return_value_policy<bp::reference_existing_object>())


    .def("__getitem__", FunctionPointer ([](const Mesh & self, PointIndex pi)
                                         {
                                           return self[pi];
                                         }))

    .def ("Add", FunctionPointer ([](Mesh & self, MeshPoint p)
                                  {
                                    return self.AddPoint (Point3d(p));
                                  }))
    ;
  

  typedef MeshingParameters MP;
  bp::class_<MP> ("MeshingParameters")
    .def("__str__", &ToString<MP>)
    .add_property("maxh", 
                  FunctionPointer ([](const MP & mp ) { return mp.maxh; }),
                  FunctionPointer ([](MP & mp, double maxh) { return mp.maxh = maxh; }))
                  
    ;
}



BOOST_PYTHON_MODULE(libmesh) {
  ExportNetgenMeshing();
}



#endif




