/*********************************************************************/
/* File:   templates.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

#ifdef USE_MKL
#include <mkl.h>
#endif // USE_MKL

#include <ngstd.hpp>
#include <ngsolve_version.hpp>
#include <netgen_version.hpp>

namespace ngstd
{
  MPI_Comm ngs_comm;

  // int printmessage_importance = 5;
  // bool NGSOStream :: glob_active = true;
  const string ngsolve_version = NGSOLVE_VERSION;


#ifdef USE_MKL
  int mkl_max_threads = [] ()
  {
      auto mkl_max = mkl_get_max_threads();
      mkl_set_num_threads(1);
      return mkl_max;
  }();
#endif // USE_MKL

  static bool dummy = [] ()
  {
    ngcore::SetLibraryVersion("ngsolve", NGSOLVE_VERSION);
    auto ng_version_compiled = ngcore::VersionInfo(NETGEN_VERSION);
    auto ng_version=ngcore::GetLibraryVersion("netgen");
    if( ng_version != ng_version_compiled )
    {
      cerr << "================================================================" << endl;
      cerr << "WARNING: NGSolve was compiled with Netgen " << endl;
      cerr << "         version " << ng_version_compiled.to_string() << " but" << endl;
      cerr << "         version " << ng_version.to_string() << " is loaded at run-time!!!" << endl;
      cerr << "================================================================" << endl;
    }
    return true;
  }();
}
