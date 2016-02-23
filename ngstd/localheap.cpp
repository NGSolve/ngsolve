/**************************************************************************/
/* File:   localheap.cpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2002                                                  */
/**************************************************************************/


#include <ngstd.hpp>

namespace ngstd
{

  LocalHeap :: LocalHeap (size_t asize, const char * aname, bool mult_by_threads)
  {
    if (mult_by_threads)
      asize *= TaskManager::GetMaxThreads();
    totsize = asize;
    try
      {
        data = new char[asize];
      }
    catch (exception & e)
      {
        throw Exception (ToString ("Could not allocate localheap, heapsize = ") + ToString(asize));
      }

    next = data + totsize;
    p = data;
    owner = true;
    name = aname;
    CleanUp();   // align pointer
  }

  LocalHeap LocalHeap :: Split() const
  {
    int pieces = TaskManager::GetNumThreads();
    int i = TaskManager::GetThreadId();
    size_t freemem = totsize - (p - data);
    size_t size_of_piece = freemem / pieces;
    return LocalHeap (p + i * size_of_piece, size_of_piece, name);
  }

  void LocalHeap :: ThrowException() // throw (LocalHeapOverflow)
  {
    /*
    cout << "allocated: " << (p-data) << endl;
    cout << "throw LocalHeapOverflow, totsize = "<< totsize << endl;
    cout << "heap name = " << name << endl;
    */
    throw LocalHeapOverflow(totsize);
  }


  LocalHeapOverflow :: LocalHeapOverflow (size_t size) 
    : Exception("Local Heap overflow\n")
  {
    stringstream str;
    str << "Current heapsize is " << size << '\n';
    Append (str.str());
    // Append ("please use 'define constant heapsize = xxx' with larger value\n");
  }
  
  LocalHeapOverflow :: ~LocalHeapOverflow ()
  {
    ; 
  }

}

