#include <solve.hpp>
using namespace ngsolve;



// #include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <string>
#include <ostream>

#include "ngs_python.hpp"

using namespace boost::python;
using std::string;
using std::ostringstream;

PythonEnvironment PythonEnvironment::instance;

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



class_<FlatArray<int> > &PyExportFlatArray(const char *name) {
    return class_<FlatArray<int> >(name)
        .def(init<int, int *>())
        .def("__str__", &PyToString<FlatArray<int> >)
        .def("__len__", &PyVecLen<FlatArray<int> >)
        .def("__getitem__", &PyVecGet<FlatArray<int> >)
        .def("__setitem__", &PyVecSet<FlatArray<int> >)
      // .def("Assign", &Assign)
      // .def("Add", &Add)
      //  .def(self+=self)
      //  .def(self-=self)
      //  .def(self*=double())
        ;
}

class_<Array<int>, bases<FlatArray<int> > > PyExportArray(const char *name) {
    return class_<Array<int> , bases<FlatArray<int> > >(name)
        .def(init<int>())
        ;
}

