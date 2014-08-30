#include <boost/python.hpp>

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



template <typename T>
void ExportArray ()
{
  string name = string("Array_") + typeid(T).name();
  bp::class_<Array<T>,boost::noncopyable>(name.c_str())
    .def ("__len__", &Array<T>::Size)
    .def ("__getitem__", 
          FunctionPointer ([](Array<T> & self, int i) -> T&
                           {
                             if (i < 0 || i >= self.Size())
                               bp::exec("raise IndexError()\n");
                             return self[i];
                           }),
          bp::return_value_policy<bp::reference_existing_object>())
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

  bp::class_<PointIndex>("PointId")
    .def("__repr__", &ToString<PointIndex>)
    .def("__str__", &ToString<PointIndex>)
    .add_property("nr", &PointIndex::operator int)
    ;
  
  bp::class_<Element>("Element3D")
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
  ;
  
  
  bp::class_<Mesh,shared_ptr<Mesh>,boost::noncopyable>("Mesh")
    .def("__str__", &ToString<Mesh>)
    .def("Load", static_cast<void(Mesh::*)(const string & name)>(&Mesh::Load))
    .def("Save", static_cast<void(Mesh::*)(const string & name)const>(&Mesh::Save))

    .def("Elements3D", 
         static_cast<Array<Element>&(Mesh::*)()> (& &Mesh::VolumeElements),
         bp::return_value_policy<bp::reference_existing_object>())
    
    /*
    .def("Elements2D", &Mesh::SurfaceElements,
         bp::return_value_policy<bp::reference_existing_object>())
    */

    ;


  bp::class_<MeshingParameters> ("MeshingParameters")
    ;
}



BOOST_PYTHON_MODULE(libmesh) {
  ExportNetgenMeshing();
}





