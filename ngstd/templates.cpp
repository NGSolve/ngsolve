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

namespace ngstd
{
  MPI_Comm ngs_comm;

  int printmessage_importance = 5;
  bool NGSOStream :: glob_active = true;
  const string ngsolve_version = NGSOLVE_VERSION;


#ifdef USE_MKL
  int mkl_max_threads = [] ()
  {
      auto mkl_max = mkl_get_max_threads();
      mkl_set_num_threads_local(1);
      return mkl_max;
  }();
#endif // USE_MKL

}
