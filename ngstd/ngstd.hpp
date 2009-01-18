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


namespace netgen
{
#include <parthreads.hpp>
#include <moveablemem.hpp>
#include <dynamicmem.hpp>

  extern ::std::ostream * testout;
  extern ::std::ostream * mycout; 
  extern int printmessage_importance;
}

using netgen::printmessage_importance;

using netgen::MoveableMem;
using netgen::DynamicMem;
using netgen::NgLock;  
using netgen::NgMutex;

using netgen::testout;
using netgen::mycout; 




/**
   namespace for standard data types and algorithms.

   Generic container classes: FlatArray, ARRAY, ArrayMem, Table, DynamicTable, HashTable, SymbolTable.

Specific data types Exception, BitArray, Flags, LocalHeap, BlockAllocator, NgProfiler, AutoPtr, EvalFunction, AutoDiff, AutoDiffDiff
*/
namespace ngstd
{
  using namespace std;

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
}

#endif
