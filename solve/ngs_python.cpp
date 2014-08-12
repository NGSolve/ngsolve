#include <solve.hpp>
using namespace ngsolve;



// #include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <string>
#include <ostream>

#include "ngs_python.hpp"

using namespace boost::python;
using std::string;
using std::ostringstream;


template<typename T>
struct PyVec : public T {
    PyVec( T &t ) : T(t) {}

    int Len() { return this->Size(); }

    double Get(int i) { 
      // if(i>=this->Size() || i<0) raiseIndexError();
        return (*this)[i];
    }

    void Set(int i, double x) {
        (*this)[i] = x;
        Ng_Redraw();
    }

    string toString() {
        std::ostringstream s;
        s << *this;
        return s.str();
    }
};

template <typename T>
PyVec<T> ToPyVec( T &t ) {
    return PyVec<T>(t);
}

////////////////////////////////////////////////////////////
template <typename T>
void PySelf(T &a) { } //return a; }

template <typename T>
void PyAssign(T &a, T &b) { a = b; }

template <typename T>
int PyVecLen(T &t) { return t.Size(); }

template <typename T>
double PyVecGet(T &t, int i) { 
  // if(i>=t.Size() || i<0) raiseIndexError();
    return t[i];
}

template <typename T>
void PyVecSet(T &t, int i, double x) {
    t[i] = x;
    Ng_Redraw();
}

template <typename T>
double PyMatGet(T &t, int i, int j) { 
 //     if(i>=t.Size() || i<0) raiseIndexError();
    return t(i,j);
}

template <typename T>
void PyMatSet(T &t, int i, int j, double x) {
    t(i,j) = x;
    Ng_Redraw();
}

template <typename T>
string PyToString(T &t ) {
    std::ostringstream s;
    s << t;
    return s.str();
}

template <typename T>
T *createObject( int size, double *data=0 ) {
    if(data)
        return new T(size);
    else
        return new T(size,data);
}

////////////////////////////////////////////////////////////
void Assign(FlatVector<double> &self, FlatVector<double> &v, double s) {
//     cout << "assign: " << s << "\t" << &self[0] << '\t' << &v[0] << endl;
    self = s*v;
}

void Add(FlatVector<double> &self, FlatVector<double> &v, double s) {
//     cout << "add   : " << s << "\t" << &self[0] << '\t' << &v[0] << endl;
    self += s*v;
}

void Mult(FlatMatrix<double> &self, FlatVector<double> &x, FlatVector<double> &y, double s) {
//     cout << "add   : " << s << "\t" << &self[0] << '\t' << &v[0] << endl;
    y = s*self*x;
}

void MultAdd(FlatMatrix<double> &self, FlatVector<double> &x, FlatVector<double> &y, double s) {
//     cout << "add   : " << s << "\t" << &self[0] << '\t' << &v[0] << endl;
    y += s*self*x;
}

void MultTrans(FlatMatrix<double> &self, FlatVector<double> &x, FlatVector<double> &y, double s) {
//     cout << "add   : " << s << "\t" << &self[0] << '\t' << &v[0] << endl;
    y = s*Trans(self)*x;
}

void MultTransAdd(FlatMatrix<double> &self, FlatVector<double> &x, FlatVector<double> &y, double s) {
//     cout << "add   : " << s << "\t" << &self[0] << '\t' << &v[0] << endl;
    y += s*Trans(self)*x;
}

class_<FlatVector<double> > &PyExportFlatVector(const char *name) {
    return class_<FlatVector<double> >(name)
        .def(init<int, double *>())
        .def("__str__", &PyToString<FlatVector<double> >)
        .def("__len__", &PyVecLen<FlatVector<double> >)
        .def("__getitem__", &PyVecGet<FlatVector<double> >)
        .def("__setitem__", &PyVecSet<FlatVector<double> >)
        .def("Assign", &Assign)
        .def("Add", &Add)
        .def(self+=self)
        .def(self-=self)
        .def(self*=double())
        ;
}

class_<Vector<double>, bases<FlatVector<double> > > PyExportVector(const char *name) {
    return class_<Vector<double> , bases<FlatVector<double> > >(name)
        .def(init<int>())
        ;
}


class_<FlatMatrix<double> > &PyExportFlatMatrix(const char *name) {
    double &(FlatMatrix<double>::*Get_p) (int, int) const;
    return class_<FlatMatrix<double> >(name)
        .def(init<int, double *>())
        .def("__str__", &PyToString<FlatMatrix<double> >)
//         .def("__len__", &PyVecLen<FlatMatrix<double> >)
//         .def("__getitem__", &PyVecGet<FlatMatrix<double> >)
//         .def("__setitem__", &PyVecSet<FlatMatrix<double> >)
//         .def("Assign", &Assign)
//         .def("Add", &Add)
        .def("Mult", &Mult)
        .def("MultAdd", &MultAdd)
        .def("MultTrans", &MultTrans)
        .def("MultTransAdd", &MultTransAdd)
        .def("Get", &PyMatGet<FlatMatrix<double> >)
        .def("Set", &PyMatSet<FlatMatrix<double> >)
        .def(self+=self)
        .def(self-=self)
        .def(self*=double())
        ;
}

class_<Matrix<double>, bases<FlatMatrix<double> > > PyExportMatrix(const char *name) {
    return class_<Matrix<double> , bases<FlatMatrix<double> > >(name)
        .def(init<int,int>())
        ;
}


// BOOST_PYTHON_MODULE(linhyp) {
// 
// 	class_<NumProcLinearHyperbolic<2> >("linhyp2", no_init )
// 		.def("Step", &NumProcLinearHyperbolic<2>::Step)
// 		.def("Redraw", &NumProcLinearHyperbolic<2>::Redraw)
//         .def_readwrite("dt", &NumProcLinearHyperbolic<2>::dt)
//         .def_readwrite("tend", &NumProcLinearHyperbolic<2>::tend)
//         .def_readwrite("t", &NumProcLinearHyperbolic<2>::t)
//         .def_readwrite("u", &NumProcLinearHyperbolic<2>::vecu)
//         .def_readwrite("hu", &NumProcLinearHyperbolic<2>::hu)
//         .def_readwrite("conv", &NumProcLinearHyperbolic<2>::conv)
//         .def_readwrite("w", &NumProcLinearHyperbolic<2>::w)
//         ;
// 
//     int (MeshAccess::*getNE_p)() const = &MeshAccess::GetNE;
//     class_<MeshAccess>("Meshaccess", no_init)
//         .def("GetNE", getNE_p)
//         .def("GetNP", &MeshAccess::GetNP)
//         .def("GetNV", &MeshAccess::GetNV)
//         .def("GetNEdges", &MeshAccess::GetNEdges)
//         .def("GetNFaces", &MeshAccess::GetNFaces)
//         .add_property("ne", getNE_p)
//         .add_property("np", &MeshAccess::GetNP)
//         .add_property("nv", &MeshAccess::GetNV)
//         .add_property("ned", &MeshAccess::GetNEdges)
//         .add_property("nf", &MeshAccess::GetNFaces)
//         ;
// };



// class_<FlatVector<double>, bases<FlatVector<double> > > PyFlatVectorD = PyExportFlatVector("FlatVector"); 
// class_<Vector<double>, bases<Vector<double> >  > PyVectorD = PyExportVector("Vector"); 

// auto PyFlatVectorD = PyExportFlatVector("FlatVector"); 
// auto PyVectorD = PyExportVector("Vector"); 

// void PythonEnvironment::Init() {
// 
//     Py_Initialize();
//     main_module = import("__main__");
//     main_namespace = main_module.attr("__dict__");
// 
//     PyRun_SimpleString("def raiseIndexError():\n\traise IndexError(\"that's enough!\")\n");
//     auto raiseIndexError = main_module.attr("raiseIndexError");
// 
//     main_namespace["FlatVector"] = PyExportFlatVector("FlatVector");
//     main_namespace["Vector"] = PyExportVector("Vector");
//     main_namespace["FlatMatrix"] = PyExportFlatMatrix("FlatMatrix");
//     main_namespace["Matrix"] = PyExportMatrix("Matrix");
// 
// }
// 
