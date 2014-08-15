#ifdef NGS_PYTHON

#include <string>
#include <ostream>
#include <type_traits>

#include "python_ngstd.hpp"

using std::string;
using std::ostringstream;

PythonEnvironment PythonEnvironment::py_env;

PythonEnvironment::PythonEnvironment() {
    mainthread_id = std::this_thread::get_id();
    pythread_id = std::this_thread::get_id();
    cout << "Init Python environment." << endl;
    Py_Initialize();
    PyEval_InitThreads();

    try{
        main_module = bp::import("__main__");
        main_namespace = main_module.attr("__dict__");

//         PyExportBla(*this);

    }
    catch(bp::error_already_set const &) {
        PyErr_Print();
    }
    PyEval_ReleaseLock();
}

#endif // NGS_PYTHON
