#ifdef NGS_PYTHON
#include "../ngstd/python_ngstd.hpp"
#include <bla.hpp>

using namespace ngbla;

void ExportNgbla() {
    std::string nested_name = "ngbla";
    if( bp::scope() )
         nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".ngbla");
    
    bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

    cout << "exporting ngbla as " << nested_name << endl;
    bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
    parent.attr("ngbla") = module ;

    bp::scope ngbla_scope(module);

    bp::object expr_module = bp::import("expr");
    bp::object expr_namespace = expr_module.attr("__dict__");
    
    typedef FlatVector<double> FVD;
    bp::class_<FVD >("FlatVector")
        .def(PyDefVector<FVD, double>()) 
        .def(PyDefToString<FVD >())
        .def("Assign", FunctionPointer( [](FVD &self, FVD &v, double s) { self  = s*v; }) )
        .def("Add",    FunctionPointer( [](FVD &self, FVD &v, double s) { self += s*v; }) )
        .def("Range",    static_cast</* const */ FVD (FVD::*)(int,int) const> (&FVD::Range ) )
        .add_property("expr", bp::object(expr_namespace["VecExpr"]) )
        .add_property("data", bp::object(expr_namespace["VecExpr"]), bp::object(expr_namespace["expr_data"] ))
        .def(bp::init<int, double *>())
        .def(bp::self+=bp::self)
        .def(bp::self-=bp::self)
        .def(bp::self*=double())
        .def("__add__" , bp::object(expr_namespace["expr_add"]) )
        .def("__sub__" , bp::object(expr_namespace["expr_sub"]) )
        .def("__rmul__" , bp::object(expr_namespace["expr_rmul"]) )
        ;

    bp::class_<Vector<double>,  bp::bases<FlatVector<double> > >("Vector")
//         .def(bp::init<int>())
        .def("__init__", bp::make_constructor (FunctionPointer ([](bp::object const & x)
            {
                FlatVector<double> *res;
                bp::extract<int> en(x); 
                if(en.check()) {
                    // call ordinary constructor Vector(int)
                    res = new Vector<double>(en());
                }
                else {
                    // Assume x is of type VecExpr
                    try{
                        int s = bp::len(x);
                        res = new Vector<double>(s);
                        bp::object exp =  bp::object(res).attr("expr");
                        x.attr("AssignTo")(exp, 1.0);
                    }
                    catch(...) {
                        PyErr_Print();
                    }
                }
                return res;
            })))
        ;

    typedef FlatMatrix<double> FMD;
    bp::class_<FlatMatrix<double> >("FlatMatrix")
        .def(PyDefToString<FMD>())
        .def("Mult",        FunctionPointer( [](FMD &m, FVD &x, FVD &y, double s) { y  = s*m*x; }) )
        .def("MultAdd",     FunctionPointer( [](FMD &m, FVD &x, FVD &y, double s) { y += s*m*x; }) )
        .def("MultTrans",   FunctionPointer( [](FMD &m, FVD &x, FVD &y, double s) { y  = s*Trans(m)*x; }) )
        .def("MultTransAdd",FunctionPointer( [](FMD &m, FVD &x, FVD &y, double s) { y += s*Trans(m)*x; }) )
//         .def("Get", FunctionPointer( [](FMD &m, int i, int j)->double { return m(i,j); } ) )
        .def("__getitem__", FunctionPointer( [](FMD &m, bp::tuple ind )->double { return m(bp::extract<int>(ind[0]), bp::extract<int>(ind[1])); } ) )
//         .def("Set", FunctionPointer( [](FMD &m, int i, int j, double val) { m(i,j) = val; }) )
        .def("__setitem__", FunctionPointer( [](FMD &m, bp::tuple ind, double val ) { m(bp::extract<int>(ind[0]), bp::extract<int>(ind[1]))=val; } ) )
        .add_property("expr", bp::object(expr_namespace["MatExpr"]) )
        .add_property("data", bp::object(expr_namespace["MatExpr"]), bp::object(expr_namespace["expr_data"] ))
        .def(bp::self+=bp::self)
        .def(bp::self-=bp::self)
        .def(bp::self*=double())
        .def("__mul__" , bp::object(expr_namespace["expr_mul"]) )
        .def("__rmul__" , bp::object(expr_namespace["expr_rmul"]) )
        .def("Height", &FMD::Height )
        .def("Width", &FMD::Width )
        ;

    bp::class_<Matrix<double>, bp::bases<FMD> >("Matrix")
        .def(bp::init<int, int>())
        ;
}


BOOST_PYTHON_MODULE(libngbla) {
    ExportNgbla();
}

#endif // NGS_PYTHON
