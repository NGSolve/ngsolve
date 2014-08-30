#include <boost/python.hpp>

#include <csg.hpp>


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



namespace netgen
{
  extern CSGeometry * ParseCSG (istream & istr);
}

void ExportCSG() 
{
  
  std::string nested_name = "csg";
  if( bp::scope() )
    nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".csg");
                                           
  bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));
  
  cout << "exporting csg " << nested_name << endl;
  bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
  parent.attr("csg") = module ;
  
  bp::scope local_scope(module);



  bp::class_<CSGeometry, boost::noncopyable> ("CSGeometry")
    .def("__init__", bp::make_constructor (FunctionPointer
                                          ([](const string & filename)
                                           {
                                             cout << "load geometry";
                                             ifstream ist(filename);
                                             shared_ptr<CSGeometry> geom(ParseCSG(ist));
                                             geom -> FindIdenticSurfaces(1e-8 * geom->MaxSize()); 
                                             return geom;
                                           })))
    .add_property ("ntlo", &CSGeometry::GetNTopLevelObjects)
    ;


  bp::def("GenerateMesh", FunctionPointer
          ([](CSGeometry & geo, MeshingParameters & param)
           {
             Mesh * dummy = NULL;
             cout << "Genrate Mesh, params = "; //  << param << endl;
             geo.GenerateMesh (dummy, param, 0, 6);
             return shared_ptr<Mesh> (dummy);
           }));
  
}





BOOST_PYTHON_MODULE(libcsg) {
  ExportCSG();
}



