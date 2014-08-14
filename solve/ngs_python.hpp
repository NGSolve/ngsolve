
#include <boost/python.hpp>
#include <solve.hpp>
#include <thread>

namespace bp = boost::python;

bp::class_<FlatVector<double> > &PyExportFlatVector(const char *name);
bp::class_<Vector<double>, bp::bases<FlatVector<double> > > PyExportVector(const char *name);
bp::class_<FlatMatrix<double> > &PyExportFlatMatrix(const char *name);
bp::class_<Matrix<double>, bp::bases<FlatMatrix<double> > > PyExportMatrix(const char *name);

bp::class_<FlatArray<int> > &PyExportFlatArray(const char *name);
bp::class_<Array<int>, bp::bases<FlatArray<int> > > PyExportArray(const char *name);


//////////////////////////////////////////////////////////////////////
// Lambda to function pointer conversion
template <typename Function>
struct function_traits
: public function_traits<decltype(&Function::operator())>
{};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...) const>
{
    typedef ReturnType (*pointer)(Args...);
};

template <typename Function>
    typename function_traits<Function>::pointer
FunctionPointer (const Function& lambda)
{
    return static_cast<typename function_traits<Function>::pointer>(lambda);
}
//////////////////////////////////////////////////////////////////////

using namespace ngsolve;
using namespace boost::python;
namespace bp = boost::python;

namespace netgen {
    extern string ngdir;
}

extern bp::object raiseIndexError;
extern AutoPtr<ngsolve::PDE> pde;
// extern bp::class_<FlatVector<double> > *PyFlatVectorD;
// extern bp::class_<Vector<double>, bp::bases<FlatVector<double> >  > *PyVectorD;
// extern auto PyFlatVectorD -> decltype(PyExportFlatVector("FlatVector")); 
// extern auto PyVectorD -> decltype(PyExportVector("Vector")); 

void InitPyNgs();
class AcquireGIL 
{
    public:
        inline AcquireGIL(){
            state = PyGILState_Ensure();
        }

        inline ~AcquireGIL(){
            PyGILState_Release(state);
        }
    private:
        PyGILState_STATE state;
};


class PythonEnvironment {
    public:
    static PythonEnvironment py_env;
    bp::object main_module; // = import("__main__");
    bp::object main_namespace; // = main_module.attr("__dict__");

    std::thread::id pythread_id;
    std::thread::id mainthread_id;

    // Private constructor
    PythonEnvironment() {
        mainthread_id = std::this_thread::get_id();
        pythread_id = std::this_thread::get_id();
        cout << "******************************" << endl;
        cout << "PythonEnvironment Constructor!" << endl;
        cout << "******************************" << endl;
        Py_Initialize();
        PyEval_InitThreads();

        try{
            main_module = bp::import("__main__");
            main_namespace = main_module.attr("__dict__");

            main_namespace["FlatVector"] = PyExportFlatVector("FlatVector");
            main_namespace["Vector"] = PyExportVector("Vector");
            main_namespace["FlatMatrix"] = PyExportFlatMatrix("FlatMatrix");
            main_namespace["Matrix"] = PyExportMatrix("Matrix");

            main_namespace["FlatArray"] = PyExportFlatArray("FlatArray");
            main_namespace["Array"] = PyExportArray("Array");

            void (*foo)(Ngs_Element &) = [](Ngs_Element &el){ cout << "hallo!" << endl; };

            bp::class_<Ngs_Element>("Ngs_Element", bp::no_init)
                .add_property("vertices2", FunctionPointer([](Ngs_Element &el)->Array<int>{ return Array<int>(el.Vertices());} ))
                .add_property("vertices", static_cast<Array<int> (*)(Ngs_Element &el)>([](Ngs_Element &el)->Array<int>{ return Array<int>(el.Vertices());} ));

            bp::class_<MeshAccess>("MeshAccess", bp::no_init)
                .def("GetNV", &MeshAccess::GetNV)
                .def("GetElement", static_cast< Ngs_Element (MeshAccess::*)(int, bool)const> (&MeshAccess::GetElement),
                        (bp::arg("arg1")=NULL,bp::arg("arg2")=0))
                .def("GetElementVertices", static_cast<void (MeshAccess::*)(int, Array<int> &) const>( &MeshAccess::GetElVertices))
                .add_property ("nv", &MeshAccess::GetNV);
        }
        catch(bp::error_already_set const &) {
            PyErr_Print();
        }
        PyEval_ReleaseLock();
        cout << "******************************" << endl;
        cout << "Constructor finished!" << endl;
        cout << "******************************" << endl;
    }

    public:
    // Singleton pattern
    static PythonEnvironment &getInstance() {
        return py_env;
    }

    void Spawn(string initfile) {
        if(pythread_id != mainthread_id) {
            cout << "Python thread already running!" << endl;
        } else {
            new std::thread([](string init_file_) {
                    try{
                    AcquireGIL gil_lock;
                    py_env.pythread_id = std::this_thread::get_id();

                    if (pde)
                      py_env.main_namespace["mesh"] = bp::ptr(&pde->GetMeshAccess(0));

                    py_env.exec_file(init_file_.c_str());
                    }
                    catch(bp::error_already_set const &) {
                    PyErr_Print();
                    }
                    cout << "Python shell finished." << endl;
                    py_env.pythread_id = py_env.mainthread_id;
                    }, initfile);
        }
        //         PyRun_SimpleString("def raiseIndexError():\n\traise IndexError(\"that's enough!\")\n");
        //         auto raiseIndexError = main_module.attr("raiseIndexError");
    }

    void exec_file(const char *file) {
        try{
            bp::exec_file(file, main_namespace, main_namespace);
        }
        catch(bp::error_already_set const &) {
            PyErr_Print();
        }
    }

    virtual ~PythonEnvironment() { }
    static Array<int> NgsElementGetVertices (Ngs_Element & el) { return Array<int> (el.Vertices()); }
};






