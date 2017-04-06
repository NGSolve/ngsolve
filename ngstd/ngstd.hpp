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

#ifdef WIN32
   #ifdef NGINTERFACE_EXPORTS
      #define DLL_HEADER   __declspec(dllexport)
   #else
      #define DLL_HEADER   __declspec(dllimport)
   #endif

   #ifdef NGS_EXPORTS
      #define NGS_DLL_HEADER   __declspec(dllexport)
   #else
      #define NGS_DLL_HEADER   __declspec(dllimport)
   #endif


#else
   #define DLL_HEADER 
// #define NGS_DLL_HEADER 


/*
   #ifdef NGINTERFACE_EXPORTS
      #define DLL_HEADER   __declspec(dllexport)
   #else
      #define DLL_HEADER   __declspec(dllimport)
   #endif
*/

   #ifdef NGS_EXPORTS
      #define NGS_DLL_HEADER   __attribute__ ((visibility ("default")))
   #else
      #define NGS_DLL_HEADER
   #endif


#endif

/*
inline void * operator new (size_t cnt)
{
  static int cnt_new = 0;
  cnt_new++;
  std::cout << "private new called, cnt = " << cnt_new << ", bytes = " << cnt << std::endl;
  return operator new(cnt, std::nothrow);
}

inline void * operator new[] (size_t cnt)
{
  static int cnt_new = 0;
  cnt_new++;
  std::cout << "private new[] called, cnt = " << cnt_new << ", bytes = " << cnt << std::endl;
  return operator new[](cnt, std::nothrow);
}
*/





// #include "dynamicmem.hpp"


namespace ngstd
{
  NGS_DLL_HEADER extern ::std::ostream * testout;
  NGS_DLL_HEADER extern int printmessage_importance;
  NGS_DLL_HEADER extern const std::string ngsolve_version;
}



/*
namespace ngstd
{
  using netgen::DynamicMem;
}
*/


/**
   namespace for standard data types and algorithms.

   Generic container classes: FlatArray, Array, ArrayMem, Table, DynamicTable, HashTable, SymbolTable.

Specific data types Exception, BitArray, Flags, LocalHeap, BlockAllocator, NgProfiler, AutoPtr, EvalFunction, AutoDiff, AutoDiffDiff
*/
namespace ngstd
{
  using namespace std;
}

#include <ngs_defines.hpp>

// #include "mycomplex.hpp"  

#include "ngs_utils.hpp"
#include "archive_base.hpp"    
#include "ngsstream.hpp"  
#include "templates.hpp"
#include "exception.hpp"
#include "localheap.hpp"
#include "profiler.hpp"

#include "simd.hpp"
#include "simd_complex.hpp"
#include "tuple.hpp"
#include "array.hpp"
#include "table.hpp"
#include "symboltable.hpp"
#include "hashtable.hpp"
#include "bitarray.hpp"

#include "blockalloc.hpp"
#include "autoptr.hpp"
#include "memusage.hpp"
#include "flags.hpp"

#include "evalfunc.hpp"
// namespace ngstd { class EvalFunction; }

#include "autodiff.hpp"
#include "autodiffdiff.hpp"
#include "polorder.hpp"
#include "stringops.hpp"
#include "statushandler.hpp"

#include "taskmanager.hpp"
#include "mpiwrapper.hpp"
#ifndef WIN32
#include "sockets.hpp"
#endif
#include "archive.hpp"

namespace ngstd
{
#ifdef WIN32
  const char dirslash = '\\';
#else
  const char dirslash = '/';
#endif
}


inline void NOOP_Deleter(void *) { ; }


#ifdef GOLD
#include <ngstd_gold.hpp>
#endif



#include "cuda_ngstd.hpp"

#endif
