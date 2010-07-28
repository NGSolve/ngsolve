/**************************************************************************/
/* File:   localheap.cpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2002                                                  */
/**************************************************************************/


#include <ngstd.hpp>

namespace ngstd
{
  using namespace ngstd;


  LocalHeap :: LocalHeap (size_t asize, const char * aname)
  {
    totsize = asize;
    data = new char[asize];
    p = data;
    owner = 1;
    name = aname;
  }

  void LocalHeap :: ThrowException() throw (LocalHeapOverflow)
  {
    cout << "allocated: " << (p-data) << endl;
    cout << "throw LocalHeapOverflow, totsize = "<< totsize << endl;
    cout << "heap name = " << name << endl;
    throw LocalHeapOverflow(totsize);
  }


  LocalHeapOverflow :: LocalHeapOverflow (size_t size) 
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

