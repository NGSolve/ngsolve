#undef NGS_EXPORTS

#include "../ngstd/python_ngstd.hpp"
#include <solve.hpp>
//using namespace ngsolve;
// #include<iostream>

/*
void __declspec(dllimport) ExportNgstd();
void __declspec(dllimport) ExportNgbla();
void __declspec(dllimport) ExportNgfem();
void __declspec(dllimport) ExportNgla();
void __declspec(dllimport) ExportNgcomp();
void __declspec(dllimport) ExportNgsolve();
*/
void NGS_DLL_HEADER ExportNgstd(py::module &m);
void NGS_DLL_HEADER ExportNgbla(py::module &m);
void NGS_DLL_HEADER ExportNgfem(py::module &m);
void NGS_DLL_HEADER ExportNgsbem(py::module &m);
void NGS_DLL_HEADER ExportNgla(py::module &m);
void NGS_DLL_HEADER ExportNgcomp(py::module &m);
void NGS_DLL_HEADER ExportNgsolve(py::module &m);

// char * libargv = (char*) { "libngs" };
// MyMPI mpiinit (1, &libargv);



PYBIND11_MODULE(ngslib, m)
{
  py::module::import("pyngcore");
  m.attr("__name__") = "ngsolve";
  /*
  py::object module(py::handle<>(py::borrowed(PyImport_AddModule("ngsolve"))));
  py::object parent = py::import("__main__");
  parent.attr("ngsolve") = module;
  py::scope local_scope(module);
  */

  // oder doch einfach nur so:
//   py::scope().attr("__name__") = "ngsolve";

  
  /*
  cout << "do some lapack stuff" << endl;
  Matrix<> a(1000), b(1000), c(1000);
  a = 1; b = 2;
  c = a*b | Lapack;
  */
#ifdef PARALLEL

  // int flag;
  // MPI_Initialized (&flag);
  // if (!flag)
  //   {
  //     cerr << "MPI not initialized !!!" << endl;
  //     cerr << "run 'from mpi4py import MPI' before importing ngsolve" << endl;
  //     throw Exception ("mpi not initialized as required from parallel-ngslib");
  //     /*
  //     const char * progname = "ngslib";
  //     typedef const char * pchar;
  //     pchar ptrs[2] = { progname, nullptr };
  //     pchar * pptr = &ptrs[0];
  //     int argc = 1;
  //     MPI_Init (&argc, (char***)&pptr);
  //     */
  //   }
  // cout << "ngslib module, rank = " << MyMPI_GetId(MPI_COMM_WORLD) << "/" << MyMPI_GetNTasks(MPI_COMM_WORLD) << endl;  
  // NGSOStream::SetGlobalActive (MyMPI_GetId(MPI_COMM_WORLD) == 0);
  
#endif

    try
    {
        m.attr("__version__") = ngsolve_version;
        py::module ngstd = m.def_submodule("ngstd", "pybind ngstd");
        ExportNgstd(ngstd);
        py::module bla = m.def_submodule("bla", "pybind bla");
        ExportNgbla(bla);
        py::module la = m.def_submodule("la", "pybind la");
        ExportNgla(la);
        py::module fem = m.def_submodule("fem", "pybind fem");
        ExportNgfem(fem);
        
        py::module comp = m.def_submodule("comp", "pybind comp");
        ExportNgcomp(comp);
        py::module bem = m.def_submodule("bem", "pybind bem");
        ExportNgsbem(bem);
        
        py::module solve = m.def_submodule("solve", "pybind solve");
        ExportNgsolve(solve);
    }
    catch (ngstd::Exception & e)
    {
        std::cerr << "\n\nCaught Exception:\n" << e.What() << std::endl;
    }
    catch (std::exception & e)
    {
        std::cerr << "\n\nCaught exception:\n" << e.what() << std::endl;
    }
    catch (...)
    {
        std::cerr << "\n\nCaught Python Exception:\n" << std::endl;
        PyErr_Print();
    }
}



