#include <solve.hpp>
using namespace ngsolve;
#include "../ngstd/python_ngstd.hpp"


/*
BasePythonEnvironment py_env;

BasePythonEnvironment & GetPythonEnvironment () 
{
  return py_env;
}
*/



void ExportNgstd();
void ExportNgbla();
void ExportNgfem();
void ExportNgla();
void ExportNgcomp();
void ExportNgsolve();


BOOST_PYTHON_MODULE(ngslib)
{
    try
    {
      ngsolve::MyMPI mympi(); // (argc, argv);
     
        ExportNgstd();
        ExportNgbla();
        ExportNgfem();
        ExportNgla();
        ExportNgcomp();
      
        ExportNgsolve();
    }
    catch (ngstd::Exception & e)
    {
        cerr << "\n\nCaught Exception:\n" << e.What() << endl;
    }
    catch (...)
    {
        cerr << "\n\nCaught Python Exception:\n" << endl;
        PyErr_Print();
    }
}


