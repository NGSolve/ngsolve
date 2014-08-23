#include <solve.hpp>
using namespace ngsolve;
#include "../ngstd/python_ngstd.hpp"

void func()
{
  cout << "hi" << endl;
  
  boost::python::exec("a = 1");
}


// NGS_DLL_HEADER AutoPtr<ngsolve::PDE> pde;


struct PyExportNgStd {
  PyExportNgStd(BasePythonEnvironment & py_env);
};

struct PyExportNgBla {
  PyExportNgBla(BasePythonEnvironment & py_env);
};

/*
struct PyExportNgComp {
  PyExportNgComp(BasePythonEnvironment & py_env);
};
*/
extern "C" PyObject * PyInit_Ngcomp();

void PyExportNgSolve (BasePythonEnvironment & py_env);


class PythonEnvironment : public BasePythonEnvironment
{
public:
  PythonEnvironment() 
  { 
    main_module = bp::import("__main__");
    main_namespace = main_module.attr("__dict__"); 

    exec("from sys import path");
    exec("from runpy import run_module");
    exec("from os import environ");
    
    exec("path.append(environ['NETGENDIR'])");
    exec("globals().update(run_module('expr',globals()))");
  }

  virtual void exec(const string s) {
    // cout << "exec from ngspy" << endl;
    bp::exec(s.c_str());
  }


  virtual void exec_file(const char *file) {
    // cout << "file exec from ngspy" << endl;
    try{
      bp::exec_file(file);
        }
    catch(bp::error_already_set const &) {
      PyErr_Print();
    }
  }
  
};


BOOST_PYTHON_MODULE(libngspy)
{
  namespace bp = boost::python;
  
  bp::def("func", func);

  PythonEnvironment py_env;
  PyExportNgStd ps(py_env);
  PyExportNgBla pbla(py_env);
  // PyExportNgComp pcomp(py_env);
  PyInit_Ngcomp();
  PyExportNgSolve (py_env);
}
