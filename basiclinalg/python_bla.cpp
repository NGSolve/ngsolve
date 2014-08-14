#ifdef NGS_PYTHON
#include "../solve/ngs_python.hpp"
#include "python_bla.hpp"

void PyExportBla(PythonEnvironment &py_env) {
    py_env["FlatVector"] = bp::class_<FlatVector<double> >("FlatVector")
        .def(PyDefVector<FlatVector<double>, double>()) 
        .def(PyDefToString<FlatVector<double> >())
        .def(bp::init<int, double *>())
        .def(bp::self+=bp::self)
        .def(bp::self-=bp::self)
        .def(bp::self*=double())
        ;

    py_env["Vector"] = bp::class_<Vector<double>, bp::bases<FlatVector<double> > >("Vector")
        .def(bp::init<int>())
        ;

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

//             py_env["FlatMatrix"] = PyExportFlatMatrix("FlatMatrix");
//             py_env["Matrix"] = PyExportMatrix("Matrix");
}

#endif // NGS_PYTHON
