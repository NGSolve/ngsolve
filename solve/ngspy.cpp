#undef NGS_EXPORTS

#include "../ngstd/python_ngstd.hpp"
#include <solve.hpp>
//using namespace ngsolve;
//#include "../ngstd/python_ngstd.hpp"
#include<boost/python.hpp>
// #include<iostream>

/*
void __declspec(dllimport) ExportNgstd();
void __declspec(dllimport) ExportNgbla();
void __declspec(dllimport) ExportNgfem();
void __declspec(dllimport) ExportNgla();
void __declspec(dllimport) ExportNgcomp();
void __declspec(dllimport) ExportNgsolve();
*/
void NGS_DLL_HEADER ExportNgstd();
void NGS_DLL_HEADER ExportNgbla();
void NGS_DLL_HEADER ExportNgfem();
void NGS_DLL_HEADER ExportNgla();
void NGS_DLL_HEADER ExportNgcomp();
void NGS_DLL_HEADER ExportNgsolve();

// char * libargv = (char*) { "libngs" };
// MyMPI mpiinit (1, &libargv);



BOOST_PYTHON_MODULE(ngslib)
{

  bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule("ngsolve"))));
  bp::object parent = bp::import("__main__");
  parent.attr("ngsolve") = module;
  bp::scope local_scope(module);

  
  /*
  cout << "do some lapack stuff" << endl;
  Matrix<> a(1000), b(1000), c(1000);
  a = 1; b = 2;
  c = a*b | Lapack;
  */
#ifdef PARALLEL

  int flag;
  MPI_Initialized (&flag);
  if (!flag)
    {
      cerr << "MPI not initialized !!!" << endl;
      cerr << "run 'from mpi4py import MPI' before importing ngsolve" << endl;
      throw Exception ("mpi not initialized as required from parallel-ngslib");
      /*
      const char * progname = "ngslib";
      typedef const char * pchar;
      pchar ptrs[2] = { progname, nullptr };
      pchar * pptr = &ptrs[0];
      int argc = 1;
      MPI_Init (&argc, (char***)&pptr);
      */
    }
  cout << "ngslib module, rank = " << MyMPI_GetId(MPI_COMM_WORLD) << "/" << MyMPI_GetNTasks(MPI_COMM_WORLD) << endl;
  
  NGSOStream::SetGlobalActive (MyMPI_GetId(MPI_COMM_WORLD) == 0);
#endif


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



