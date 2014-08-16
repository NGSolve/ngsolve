#ifdef NGS_PYTHON

#include <string>
#include <ostream>
#include <type_traits>

#include "python_ngstd.hpp"

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

        exec("from sys import path");
        exec("from runpy import run_module");
        exec("from os import environ");

        exec("path.append(environ['NETGENDIR'])");
        exec("globals().update(run_module('expr',globals()))");
    }
    catch(bp::error_already_set const &) {
        PyErr_Print();
    }

    // Export ngstd classes
    py_env["FlatArrayD"] = bp::class_<FlatArray<double> >("FlatArrayD")
        .def(PyDefVector<FlatArray<double>, double>()) 
        .def(PyDefToString<FlatArray<double> >())
        .def(bp::init<int, double *>())
        ;

    py_env["ArrayD"] = bp::class_<Array<double>, bp::bases<FlatArray<double> > >("ArrayD")
        .def(bp::init<int>())
        ;

    py_env["FlatArrayI"] = bp::class_<FlatArray<int> >("FlatArrayI")
        .def(PyDefVector<FlatArray<int>, int>()) 
        .def(PyDefToString<FlatArray<int> >())
        .def(bp::init<int, int *>())
        ;

    py_env["ArrayI"] = bp::class_<Array<int>, bp::bases<FlatArray<int> > >("ArrayI")
        .def(bp::init<int>())
        ;
}

#endif // NGS_PYTHON
