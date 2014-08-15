#ifdef NGS_PYTHON
#include <solve.hpp>
using namespace ngsolve;

#include <string>
#include <ostream>
#include <type_traits>

#include "ngs_python.hpp"
#include "../basiclinalg/python_bla.hpp"
#include "../basiclinalg/python_bla.cpp"  // a hack

using std::string;
using std::ostringstream;

PythonEnvironment PythonEnvironment::py_env;

PythonEnvironment::PythonEnvironment() {
    mainthread_id = std::this_thread::get_id();
    pythread_id = std::this_thread::get_id();
    cout << "Init Python environment." << endl;
    Py_Initialize();
    PyEval_InitThreads();

    try{
        main_module = bp::import("__main__");
        main_namespace = main_module.attr("__dict__");

        PyExportBla(*this);

        void (*foo)(Ngs_Element &) = [](Ngs_Element &el){ cout << "hallo!" << endl; };

        bp::class_<Ngs_Element>("Ngs_Element", bp::no_init)
            .add_property("vertices", FunctionPointer([](Ngs_Element &el)->Array<int>{ return Array<int>(el.Vertices());} ))
            ;

        bp::class_<MeshAccess>("MeshAccess", bp::no_init)
            .def("GetNV", &MeshAccess::GetNV)
            .def("GetElement", static_cast< Ngs_Element (MeshAccess::*)(int, bool)const> (&MeshAccess::GetElement),
                    (bp::arg("arg1")=NULL,bp::arg("arg2")=0))
            .def("GetElementVertices", static_cast<void (MeshAccess::*)(int, Array<int> &) const>( &MeshAccess::GetElVertices))
            .add_property ("nv", &MeshAccess::GetNV)
            ;
    }
    catch(bp::error_already_set const &) {
        PyErr_Print();
    }
    PyEval_ReleaseLock();
}

#endif // NGS_PYTHON
