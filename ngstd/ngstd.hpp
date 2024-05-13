#ifndef FILE_NGSTD
#define FILE_NGSTD

/*********************************************************************/
/* File:   ngstd.hpp                                                 */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   ng-standard classes
*/

#include <ngs_stdcpp_include.hpp>

namespace ngstd
{
  // NGS_DLL_HEADER extern int printmessage_importance;
  NGS_DLL_HEADER extern const std::string ngsolve_version;
}


/**
   namespace for standard data types and algorithms.

   Generic container classes: FlatArray, Array, ArrayMem

Specific data types Exception, BlockAllocator, EvalFunction, AutoDiff, AutoDiffDiff
*/


#include <ngs_defines.hpp>
#include <core/ngcore.hpp>
#include <core/autodiff.hpp>
#include <core/autodiffdiff.hpp>

namespace ngstd
{
  using namespace ngcore;
  // using ngcore::INT;
} // namespace ngstd

#include "ngs_utils.hpp"

#include "blockalloc.hpp"
#include "memusage.hpp"

#include "evalfunc.hpp"
#include "sample_sort.hpp"

// #include "polorder.hpp"
#include "stringops.hpp"
#include "statushandler.hpp"

namespace ngstd
{
#ifdef WIN32
  const char dirslash = '\\';
#else
  const char dirslash = '/';
#endif

  
  // using ngcore::NgMPI_Comm;
  enum { NG_MPI_TAG_CMD = 110 };
  enum { NG_MPI_TAG_SOLVE = 1110 };
}


inline void NOOP_Deleter(void *) { ; }

#endif
