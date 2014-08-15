#ifndef PYTHON_NGSTD_HPP___
#define PYTHON_NGSTD_HPP___
#ifdef NGS_PYTHON

#include <boost/python.hpp>
// #include <solve.hpp>
#include <ngstd.hpp>
#include <thread>
namespace bp = boost::python;


// using namespace ngsolve;
namespace bp = boost::python;

/*
namespace netgen {
    extern string ngdir;
}

extern AutoPtr<ngsolve::PDE> pde;
*/

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
    bp::object main_module; 
    bp::object main_namespace; 

    std::thread::id pythread_id;
    std::thread::id mainthread_id;

    // Private constructor
    PythonEnvironment();

    public:
    // Singleton pattern
    static PythonEnvironment &getInstance() {
        return py_env;
    }
    
    auto operator[] ( const char *s ) -> decltype(main_namespace[s]) {
        return main_namespace[s];
    }


    void Spawn(string initfile) {
        if(pythread_id != mainthread_id) {
            cout << "Python thread already running!" << endl;
        } else {
            std::thread([](string init_file_) {
                    try{
                    AcquireGIL gil_lock;
                    py_env.pythread_id = std::this_thread::get_id();
                    py_env.exec_file(init_file_.c_str());
                    }
                    catch(bp::error_already_set const &) {
                    PyErr_Print();
                    }
                    cout << "Python shell finished." << endl;
                    py_env.pythread_id = py_env.mainthread_id;
                    }, initfile).detach();
        }
    }

    void exec(const char *s) {
        bp::exec(s);
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
  // static Array<int> NgsElementGetVertices (Ngs_Element & el) { return Array<int> (el.Vertices()); }
};


//////////////////////////////////////////////////////////////////////
// Lambda to function pointer conversion
template <typename Function>
struct function_traits
: public function_traits<decltype(&Function::operator())> {};

template <typename ClassType, typename ReturnType, typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...) const> {
    typedef ReturnType (*pointer)(Args...);
};

template <typename Function>
typename function_traits<Function>::pointer
FunctionPointer (const Function& lambda) {
    return static_cast<typename function_traits<Function>::pointer>(lambda);
}


//////////////////////////////////////////////////////////////////////
template< typename T>
struct PyDefToString : public boost::python::def_visitor<PyDefToString<T> > {
    template <class Tclass>
        void visit(Tclass& c) const {
            c.def("__str__", &PyDefToString<T>::ToString);
        }

    static string ToString(T &t ) {
        std::ostringstream s;
        s << t;
        return s.str();
    }
};


//////////////////////////////////////////////////////////////////////
template< typename T, typename TELEM = double >
struct PyDefBracketOperator : public boost::python::def_visitor<PyDefBracketOperator<T,TELEM> > {
    template <class Tclass>
        void visit(Tclass& c) const {
            c
                .def("__getitem__", &PyDefBracketOperator<T,TELEM>::Get)
                .def("__setitem__", &PyDefBracketOperator<T,TELEM>::Set)
                ;
        }

    static TELEM Get(T& self, int i) { 
        if( i<self.Size() && i>=0 )
            return self[i];
        RaiseIndexError();
        return TELEM();
    } 

    static void Set(T& self, int i, TELEM val) {
        if( i<self.Size() && i>=0 )
            self[i] = val;
        else
            RaiseIndexError();
    }

    static void RaiseIndexError() {
        PythonEnvironment::getInstance().exec("raise IndexError()\n");
    }

};

//////////////////////////////////////////////////////////////////////
// Python iterator protocoll
template <typename T>
class PyIterator {
    T &v;
    int size;
    int index;
    int startindex;
    typedef typename std::remove_reference<decltype(v[startindex])>::type TELEM;
    
    public:
    PyIterator(T &v_, int size_, int startindex_ = 0) : v(v_), size(size_), index(startindex_), startindex(startindex_) {}

    TELEM Next() { 
        if(index<startindex+size) return v[index++];
        else PythonEnvironment::getInstance().exec("raise StopIteration()\n");
        return TELEM();
    }

    static void Export () {
        bp::class_<PyIterator<T> >("PyIterator", bp::no_init).def("__next__", &PyIterator<T>::Next);
    }
};




//////////////////////////////////////////////////////////////////////
// Export len, bracket operator and iterator protocoll at once
template< typename T,  typename TELEM = double>
struct PyDefVector : public boost::python::def_visitor<PyDefVector<T,TELEM> > {
    template <class Tclass>
        void visit(Tclass& c) const {
            PyIterator<T>::Export();
            c
                .def("__len__", FunctionPointer( []( T& v) { return v.Size();} ) )
                .def(PyDefBracketOperator<T, TELEM>())
                .def("__iter__", FunctionPointer([](T &v) { return PyIterator<T>(v, v.Size(), 0 ); }))
                ;
        }
};


#endif // NGS_PYTHON
#endif // PYTHON_NGSTD_HPP___
