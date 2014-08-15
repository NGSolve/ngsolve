#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <comp.hpp>
using namespace ngcomp;


struct PyExportComp {
    PyExportComp() {
        auto py_env = PythonEnvironment::getInstance();


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

    }
};


// Call constructor to export python classes
static PyExportComp python_export_comp;

#endif // NGS_PYTHON
