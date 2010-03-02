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



#ifdef WIN32 
// trick from http://social.msdn.microsoft.com/Forums/en/vclanguage/thread/ab642c88-2d2d-4f5d-9fd7-2341442d5a46
// all new/delete allocation is done from ngsolve heap

NGS_DLL_HEADER void * __cdecl my_operator_new_replacement(size_t _count);
NGS_DLL_HEADER void __cdecl my_operator_delete_replacement(void * _ptr);
NGS_DLL_HEADER void * __cdecl my_operator_new_array_replacement(size_t _count);
NGS_DLL_HEADER void __cdecl my_operator_delete_array_replacement(void * _ptr);

#ifndef NGS_EXPORTS
inline void * __cdecl operator new(size_t _count) {
    return my_operator_new_replacement(_count);
}
inline void __cdecl operator delete(void * _ptr) {
    my_operator_delete_replacement(_ptr);
}
inline void * __cdecl operator new[](size_t _count) {
    return my_operator_new_array_replacement(_count);
}
inline void __cdecl operator delete[](void * _ptr) {
    my_operator_delete_array_replacement(_ptr);
}

#endif

#endif


#include "dynamicmem.hpp"

namespace netgen
{
  DLL_HEADER extern ::std::ostream * testout;
  DLL_HEADER extern int printmessage_importance;
}

using netgen::printmessage_importance;
using netgen::DynamicMem;
using netgen::testout;



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

#endif
