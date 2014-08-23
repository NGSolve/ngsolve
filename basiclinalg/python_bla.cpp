#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <bla.hpp>
#include <memory>
using std::shared_ptr;
using namespace ngbla;

struct PyExportNgBla {
  PyExportNgBla(BasePythonEnvironment & py_env);
  PyExportNgBla();
};

static double test_double = 1;

// PyExportNgBla ::
//     PyExportNgBla(BasePythonEnvironment & py_env) {
//     }


BOOST_PYTHON_MODULE(Ngbla) {
    bp::def("Testfunc", FunctionPointer([]() { cout << "TEST!" <<endl; } ));
      // auto py_env = PythonEnvironment::getInstance();
        typedef FlatVector<double> FVD;
//         py_env["FlatVector"] = 
        bp::class_<FVD >("FlatVector")
            .def(PyDefVector<FVD, double>()) 
            .def(PyDefToString<FVD >())
            .def("Assign", FunctionPointer( [](FVD &self, FVD &v, double s) { self  = s*v; }) )
            .def("Add",    FunctionPointer( [](FVD &self, FVD &v, double s) { self += s*v; }) )
          .def("Range",    static_cast</* const */ FVD (FVD::*)(int,int) const> (&FVD::Range ) )
            .def(bp::init<int, double *>())
            .def(bp::self+=bp::self)
            .def(bp::self-=bp::self)
            .def(bp::self*=double())
            .def("GetTest", FunctionPointer( [] (FVD &self) -> double { return test_double; } ) )
            .def("GetTest1", FunctionPointer( [] (FVD &self){ return bp::ptr(test_double); } ) )

            .def("copy1", FunctionPointer( [] (FVD &self) -> FVD & { return self; } ),
            bp::return_value_policy<bp::reference_existing_object>())

            .def("copy2", FunctionPointer( [] (FVD &self) -> const FVD & { return self; } ),
            bp::return_value_policy<bp::copy_const_reference>())
            .def("copy3", FunctionPointer( [] (FVD &self) -> FVD & { return self; } ),
            bp::return_value_policy<bp::copy_non_const_reference>())
            .def("copy4", FunctionPointer( [] (FVD &self) -> FVD & { return self; } ),
            bp::return_internal_reference<>())
            .def("copy5", FunctionPointer( [] (FVD &self) -> const FVD & { return self; } ),
            bp::return_value_policy<bp::reference_existing_object>())

            ;
//         PyEnableVecExpr(py_env, "FlatVector");
//         py_env["Vector"] = 
        bp::class_<Vector<double>, bp::bases<FlatVector<double> > >("Vector")
            .def(bp::init<int>())
            ;

        typedef FlatMatrix<double> FMD;
//         py_env["FlatMatrix"] = 
        bp::class_<FlatMatrix<double> >("FlatMatrix")
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

//         PyRun_SimpleString("FlatVector.aaa = 1112");
//         PyEnableMatExpr(py_env, "FlatMatrix");

//         py_env["Matrix"] = 
        bp::class_<Matrix<double>, bp::bases<FMD> >("Matrix")
            .def(bp::init<int, int>())
            ;

        // execute bla.py
        // py_env.exec("globals().update(run_module('bla',globals()))");
    }

PyExportNgBla ::
    PyExportNgBla(BasePythonEnvironment & py_env) {
    cout << "appendInittab" << endl;
    PyImport_AppendInittab("Ngbla", PyInit_Ngbla);
    }

PyExportNgBla ::
    PyExportNgBla() {
    cout << "appendInittab" << endl;
    PyImport_AppendInittab("Ngbla", PyInit_Ngbla);
    }


// Call constructor to export python classes
// static PyExportNgBla python_export_bla ();

#endif // NGS_PYTHON
