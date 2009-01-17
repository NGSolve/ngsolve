/**************************************************************************/
/* File:   localheap.cpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2002                                                  */
/**************************************************************************/


#include <ngstd.hpp>

namespace ngstd
{
  using namespace ngstd;


#ifndef _OPENMP
  LocalHeap :: LocalHeap (unsigned int asize)
  {
    totsize = asize;
    data = new char[asize];
    p = data;
    owner = 1;
  }
#endif


  void LocalHeap :: ThrowException() throw (LocalHeapOverflow)
  {
    throw LocalHeapOverflow(totsize);
  }


  LocalHeapOverflow :: LocalHeapOverflow (int size) 
    : Exception("Local Heap overflow\n")
  {
    stringstream str;
    str << "Current heapsize is " << size << '\n';
    Append (str.str());
    Append ("please use 'define constant heapsize = xxx' with larger value\n");
  }
  
  LocalHeapOverflow :: ~LocalHeapOverflow ()
  {
    ; 
  }

}

