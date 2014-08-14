#ifdef NGS_PYTHON
#include "../solve/ngs_python.hpp"
#include "python_bla.hpp"

void PyExportBla(PythonEnvironment &py_env) {
    typedef FlatVector<double> FVD;
    py_env["FlatVector"] = bp::class_<FVD >("FlatVector")
        .def(PyDefVector<FVD, double>()) 
        .def(PyDefToString<FVD >())
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

    typedef FlatMatrix<double> FMD;
    py_env["FlatMatrix"] = bp::class_<FlatMatrix<double> >("FlatMatrix")
        .def(PyDefToString<FMD>())
        .def("Mult",        FunctionPointer( [](FMD &m, FVD &x, FVD &y, double s) { y  = s*m*x; }) )
        .def("MultAdd",     FunctionPointer( [](FMD &m, FVD &x, FVD &y, double s) { y += s*m*x; }) )
        .def("MultTrans",   FunctionPointer( [](FMD &m, FVD &x, FVD &y, double s) { y  = s*Trans(m)*x; }) )
        .def("MultTransAdd",FunctionPointer( [](FMD &m, FVD &x, FVD &y, double s) { y += s*Trans(m)*x; }) )
        .def("Get", FunctionPointer( [](FMD &m, int i, int j)->double { return m(i,j); } ) )
        .def("Set", FunctionPointer( [](FMD &m, int i, int j, double val) { m(i,j) = val; }) )
        .def(bp::self+=bp::self)
        .def(bp::self-=bp::self)
        .def(bp::self*=double())
        ;

    py_env["Matrix"] = bp::class_<Matrix<double>, bp::bases<FMD> >("Matrix")
        .def(bp::init<int, int>())
        ;

}

#endif // NGS_PYTHON
