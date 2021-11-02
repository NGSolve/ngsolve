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
#ifdef NETGEN_ENABLE_CHECK_RANGE
    const bool ngs_check_range = true;
#else  // NETGEN_ENABLE_CHECK_RANGE
    const bool ngs_check_range = false;
#endif // NETGEN_ENABLE_CHECK_RANGE
    const bool ng_check_range = IsRangeCheckEnabled();

    const int ngs_simd_width = GetDefaultSIMDSize();
    const int ng_simd_width = GetCompiledSIMDSize();

    if(ngs_check_range != ng_check_range || ngs_simd_width != ng_simd_width)
    {
        stringstream s;
        s << "Incompatible version of Netgen loaded!" << endl;
        s << "Range checks enabled (Negen, NGSolve): " << ng_check_range << "\t" << ngs_check_range << endl;
        s << "SIMD width (Negen, NGSolve):           " << ng_simd_width << "\t" << ngs_simd_width << endl;
        throw runtime_error(s.str());
    }

    return true;
  }();
}
