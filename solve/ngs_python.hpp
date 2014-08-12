
#include <boost/python.hpp>

// using namespace boost::python;
namespace bp = boost::python;

class PythonEnvironment {
    // Private constructor
    PythonEnvironment() {}

    bp::object main_module; // = import("__main__");
    bp::object main_namespace; // = main_module.attr("__dict__");
    public:
    // Singleton pattern
    static PythonEnvironment &getInstance() {
        static PythonEnvironment instance;
        return instance;
    }

    void Init();
    
    void exec_file(const char *file) {
        bp::exec_file(file, main_namespace, main_namespace);
    }


    virtual ~PythonEnvironment() { }
};

extern bp::object raiseIndexError;
extern bp::class_<FlatVector<double> > *PyFlatVectorD;
extern bp::class_<Vector<double>, bp::bases<FlatVector<double> >  > *PyVectorD;
// extern auto PyFlatVectorD -> decltype(PyExportFlatVector("FlatVector")); 
// extern auto PyVectorD -> decltype(PyExportVector("Vector")); 

// void InitPyNgs();


