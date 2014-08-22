#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <comp.hpp>
using namespace ngcomp;



struct PyExportNgComp {
  PyExportNgComp(BasePythonEnvironment & py_env);
};


enum color { red = 1, green = 2, blue = 4 };
color identity_ (color x) { return x; }

PyExportNgComp :: PyExportNgComp(BasePythonEnvironment & py_env) {

  /*
    // embedded py crasht hier, shared-libray klappt:
  //   py_env["color"] = 
  bp::enum_<color>("color")
    .value("red", red)
    .value("green", green)
    .value("blue", blue)
    ;
  bp::def ("identity", identity_);
  */
  /*
  // py_env["VorB"] = 
  bp::enum_<VorB>("VorB")
    .value("VOL", VOL)
    .value("BND", BND)
    .export_values()
    ;
    // liefert py-runtimer Error ?!?!?!?!

  bp::class_<ElementId> ("ElementId", bp::init<VorB,int>())
    .def(PyDefToString<ElementId>())
    .add_property("nr", &ElementId::Nr)    
    .def("IsVolume", &ElementId::IsVolume)
    .def("IsBoundary", &ElementId::IsBoundary)
    .def(PyDefToString<FESpace>())
    ;
  */


  bp::class_<Ngs_Element>("Ngs_Element", bp::no_init)
    .add_property("vertices", FunctionPointer([](Ngs_Element &el)->Array<int>{ return Array<int>(el.Vertices());} ))
    .add_property("edges", FunctionPointer([](Ngs_Element &el)->Array<int>{ return Array<int>(el.Edges());} ))
    .add_property("faces", FunctionPointer([](Ngs_Element &el)->Array<int>{ return Array<int>(el.Faces());} ))
    ;

  py_env["MeshAccess"] = 
    bp::class_<MeshAccess>("MeshAccess", "meshclass doc", bp::no_init)
    // .def("GetNV", &MeshAccess::GetNV)
    .def("GetElement", static_cast< Ngs_Element (MeshAccess::*)(int, bool)const> (&MeshAccess::GetElement),
         (bp::arg("arg1")=NULL,bp::arg("arg2")=0))
    .def("GetElementVertices", static_cast<void (MeshAccess::*)(int, Array<int> &) const>( &MeshAccess::GetElVertices))
    .add_property ("nv", &MeshAccess::GetNV, "number of vertices")
    .add_property ("ne", static_cast<int(MeshAccess::*)()const> (&MeshAccess::GetNE))
    ;



  bp::class_<FESpace,boost::noncopyable>("FESpace", bp::no_init)
    .def("GetDofNrs", FunctionPointer([](FESpace & self, int i) { Array<int> tmp; self.GetDofNrs(i,tmp); return tmp; }))
    .add_property ("ndof", FunctionPointer([](FESpace & self) { return self.GetNDof(); }))
    .def(PyDefToString<FESpace>())
    ;

  bp::class_<GridFunction,boost::noncopyable>("GridFunction", bp::no_init)
    .def(PyDefToString<GridFunction>())
    ;

  bp::class_<BilinearForm,boost::noncopyable>("BilinearForm", bp::no_init)
    .def(PyDefToString<BilinearForm>())
    ;



}


// Call constructor to export python classes
// static PyExportComp python_export_comp (PythonEnvironment::getInstance());

#endif // NGS_PYTHON
