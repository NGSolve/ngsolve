#include <solve.hpp>
//using namespace ngsolve;
//#include "../ngstd/python_ngstd.hpp"
#include<boost\python.hpp>
#include<iostream>

void __declspec(dllimport) ExportNgstd();
void __declspec(dllimport) ExportNgbla();
void __declspec(dllimport) ExportNgfem();
void __declspec(dllimport) ExportNgla();
void __declspec(dllimport) ExportNgcomp();
void __declspec(dllimport) ExportNgsolve();


// char * libargv = (char*) { "libngs" };
// MyMPI mpiinit (1, &libargv);



BOOST_PYTHON_MODULE(ngslib)
{

  /*
  cout << "do some lapack stuff" << endl;
  Matrix<> a(1000), b(1000), c(1000);
  a = 1; b = 2;
  c = a*b | Lapack;
  */

    try
    {
        ExportNgstd();
        ExportNgbla();
        ExportNgfem();
        ExportNgla();
        ExportNgcomp();      
        ExportNgsolve();
    }
    catch (ngstd::Exception & e)
    {
        std::cerr << "\n\nCaught Exception:\n" << e.What() << std::endl;
    }
    catch (...)
    {
        std::cerr << "\n\nCaught Python Exception:\n" << std::endl;
        PyErr_Print();
    }
}



