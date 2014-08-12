
#include <boost/python.hpp>

namespace bp = boost::python;

bp::class_<FlatVector<double> > &PyExportFlatVector(const char *name);
bp::class_<Vector<double>, bp::bases<FlatVector<double> > > PyExportVector(const char *name);
bp::class_<FlatMatrix<double> > &PyExportFlatMatrix(const char *name);
bp::class_<Matrix<double>, bp::bases<FlatMatrix<double> > > PyExportMatrix(const char *name);

bp::class_<FlatArray<int> > &PyExportFlatArray(const char *name);
bp::class_<Array<int>, bp::bases<FlatArray<int> > > PyExportArray(const char *name);



// using namespace boost::python;
// namespace bp = boost::python;
// 
// class PythonEnvironment {
//     // Private constructor
//     PythonEnvironment() {}
// 
//     bp::object main_module; // = import("__main__");
//     bp::object main_namespace; // = main_module.attr("__dict__");
//     public:
//     // Singleton pattern
//     static PythonEnvironment &getInstance() {
//         static PythonEnvironment instance;
//         return instance;
//     }
// 
//     void Init() {
// 
//         Py_Initialize();
//         main_module = bp::import("__main__");
//         main_namespace = main_module.attr("__dict__");
//         cout << "import readline..." << endl;
// //         main_module = bp::import("readline");
//         cout << "...finished" << endl;
// 
//         PyRun_SimpleString("def raiseIndexError():\n\traise IndexError(\"that's enough!\")\n");
//         auto raiseIndexError = main_module.attr("raiseIndexError");
// 
// //         main_namespace["FlatVector"] = PyExportFlatVector("FlatVector");
// //         main_namespace["Vector"] = PyExportVector("Vector");
// //         main_namespace["FlatMatrix"] = PyExportFlatMatrix("FlatMatrix");
// //         main_namespace["Matrix"] = PyExportMatrix("Matrix");
// 
//     }
// 
//     void exec_file(const char *file) {
//         bp::exec_file(file, main_namespace, main_namespace);
//     }
// 
// 
//     virtual ~PythonEnvironment() { }
// };
// 
// extern bp::object raiseIndexError;
// extern bp::class_<FlatVector<double> > *PyFlatVectorD;
// extern bp::class_<Vector<double>, bp::bases<FlatVector<double> >  > *PyVectorD;
// extern auto PyFlatVectorD -> decltype(PyExportFlatVector("FlatVector")); 
// extern auto PyVectorD -> decltype(PyExportVector("Vector")); 

// void InitPyNgs();





static Array<int> NgsElementGetVertices (Ngs_Element & el) { return Array<int> (el.Vertices()); }

