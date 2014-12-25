#include <solve.hpp>
using namespace ngsolve;
#include "../ngstd/python_ngstd.hpp"

#ifdef NGS_PYTHON

void ExportNgstd();
void ExportNgbla();
void ExportNgfem();
void ExportNgla();
void ExportNgcomp();
void ExportNgsolve();


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
        cerr << "\n\nCaught Exception:\n" << e.What() << endl;
    }
    catch (...)
    {
        cerr << "\n\nCaught Python Exception:\n" << endl;
        PyErr_Print();
    }
}



#endif
