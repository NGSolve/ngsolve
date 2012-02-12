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
   #define NGS_DLL_HEADER 
#endif



#include "dynamicmem.hpp"


namespace ngstd
{
  DLL_HEADER extern ::std::ostream * testout;
}
using ngstd::testout;

namespace netgen
{
  // DLL_HEADER extern ::std::ostream * testout;
  DLL_HEADER extern int printmessage_importance;
}



#ifndef PARALLEL

enum { id = 0 };
enum { ntasks = 1 };

#else

namespace netgen {
  extern int id, ntasks;
}
using netgen::id;
using netgen::ntasks;

#endif



// using netgen::printmessage_importance;
using netgen::DynamicMem;
// using netgen::testout;



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




  
#include "ngsstream.hpp"  
#include "templates.hpp"
#include "exception.hpp"
#include "localheap.hpp"

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
#include "autodiff.hpp"
#include "autodiffdiff.hpp"
#include "stringops.hpp"
#include "profiler.hpp"
#include "statushandler.hpp"

namespace ngstd
{
#ifdef WIN32
  const char dirslash = '\\';
#else
  const char dirslash = '/';
#endif
}


#endif
