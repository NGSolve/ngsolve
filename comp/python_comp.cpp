#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <comp.hpp>
using namespace ngcomp;



struct PyExportNgComp {
  PyExportNgComp(BasePythonEnvironment & py_env);
};


PyExportNgComp :: PyExportNgComp(BasePythonEnvironment & py_env) {

  void (*foo)(Ngs_Element &) = [](Ngs_Element &el){ cout << "hallo!" << endl; };

  bp::class_<Ngs_Element>("Ngs_Element", bp::no_init)
    .add_property("vertices", FunctionPointer([](Ngs_Element &el)->Array<int>{ return Array<int>(el.Vertices());} ))
    .add_property("edges", FunctionPointer([](Ngs_Element &el)->Array<int>{ return Array<int>(el.Edges());} ))
    .add_property("faces", FunctionPointer([](Ngs_Element &el)->Array<int>{ return Array<int>(el.Faces());} ))
    ;


  bp::class_<MeshAccess>("MeshAccess", bp::no_init)
    // .def("GetNV", &MeshAccess::GetNV)
    .def("GetElement", static_cast< Ngs_Element (MeshAccess::*)(int, bool)const> (&MeshAccess::GetElement),
         (bp::arg("arg1")=NULL,bp::arg("arg2")=0))
    .def("GetElementVertices", static_cast<void (MeshAccess::*)(int, Array<int> &) const>( &MeshAccess::GetElVertices))
    .add_property ("nv", &MeshAccess::GetNV)
    .add_property ("ne", static_cast<int(MeshAccess::*)()const> (&MeshAccess::GetNE))
    ;



  bp::class_<FESpace,boost::noncopyable>("FESpace", bp::no_init)
    .def("GetDofNrs", FunctionPointer([](FESpace & self, int i) { Array<int> tmp; self.GetDofNrs(i,tmp); return tmp; }))
    .add_property ("ndof", FunctionPointer([](FESpace & self) { return self.GetNDof(); }))
    .def(PyDefToString<FESpace>())
    ;

}


// Call constructor to export python classes
// static PyExportComp python_export_comp (PythonEnvironment::getInstance());

#endif // NGS_PYTHON
