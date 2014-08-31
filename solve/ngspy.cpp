#include <solve.hpp>
using namespace ngsolve;
#include "../ngstd/python_ngstd.hpp"



BasePythonEnvironment py_env;

BasePythonEnvironment & GetPythonEnvironment () 
{
  return py_env;
}




int ExportNgstd();
int ExportNgbla();
int ExportNgfem();
int ExportNgla();
int ExportNgcomp();
int ExportNgsolve();


BOOST_PYTHON_MODULE(libngspy)
{
  ExportNgstd();
  ExportNgbla();
  ExportNgfem();
  ExportNgla();
  ExportNgcomp();
  ExportNgsolve();
}
