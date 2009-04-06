#ifndef FILE_LOCALHEAP
#define FILE_LOCALHEAP

/**************************************************************************/
/* File:   localheap.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2000                                                  */
/**************************************************************************/


/**
   Exception on heap overflow.
   Thrown by allocation on LocalHeap.
*/
class LocalHeapOverflow : public Exception
{
public:
  LocalHeapOverflow (int size);
  virtual ~LocalHeapOverflow ();
};
 


#ifndef _OPENMP

/**
   Optimized memory handler.
   One block of data is organized as stack memory. 
   One can allocate memory out of it. This increases the stack pointer.
   With \Ref{CleanUp}, the pointer is reset to the beginning or to a
   specific position. 
 */
class LocalHeap
{
  char * data;
  char * p;
  unsigned int totsize;
  bool owner;

public:
  /// Allocate one block of size asize.
  LocalHeap (unsigned int asize);

  /// Use provided memory for the LocalHeap
  LocalHeap (char * adata, unsigned int asize) throw ()
  {
    totsize = asize;
    data = adata;
    owner = 0;
    // p = data;
    CleanUp();
  }
  
  /// free memory
  ~LocalHeap ()
  {
    if (owner)
      delete [] data;
  }
  
  /// delete all memory on local heap
  void CleanUp() throw ()
  {
    // p = data;
    p = data;
    p += (16 - (long(p) & 15) );
  }

  /// returns heap-pointer
  void * GetPointer () throw ()
  {
    return p;
  }

  /// deletes memory back to heap-pointer
  void CleanUp (void * addr) throw ()
  {
    p = (char*)addr;
  }

  /// allocates size bytes of memory from local heap
  void * Alloc (unsigned int size) throw (LocalHeapOverflow)
  {
    char * oldp = p;
    
    // 16 byte allignment
    size += (16 - size % 16);
    p += size;

    if ( (p - data) >= int(totsize) )
      ThrowException();

    return oldp;
  }

  /// allocates size objects of type T on local heap
  template <typename T>
  T * Alloc (unsigned int size) throw (LocalHeapOverflow)
  {
    char * oldp = p;
    size *= sizeof (T);

    // 16 byte allignment
    size += (16 - size % 16);
    p += size;

    if ( (p - data) >= int(totsize) )
      ThrowException();

    return reinterpret_cast<T*> (oldp);
  }


  ///
  void ThrowException() throw (LocalHeapOverflow);

  /// free memory (dummy function)
  void Free (void * data) throw () 
  {
    ;
  }

  /// available memory on LocalHeap
  int Available () const throw () { return (totsize - (p-data)); }
};


/**
   Optimized memory handler.
   Provides static memory for the local heap. The template argument specifies the size in number of chars.
*/
template <int S>
class LocalHeapMem : public LocalHeap
{
  char mem[S];
public:
  LocalHeapMem () throw () : LocalHeap (mem, S) { ; }
};




#else



class LocalHeap
{
  enum { MAX_THREADS = 16 };
  char * data[MAX_THREADS];
  char * p[MAX_THREADS];
  unsigned int totsize;
  bool owner;

public:
  LocalHeap (unsigned int asize)
  {
    if (omp_in_parallel())
      {
        int maxth = omp_get_max_threads();
        for (int i = 0; i < maxth; i++)
          {
            data[i] = 0;
            p[i] = 0;
          }

        int tid = omp_get_thread_num();

        data[tid] = new char[asize];
        p[tid] = data[tid];
        totsize = asize;
        owner = 1;
      }
    else
      {
        int maxth = omp_get_max_threads();
        for (int i = 0; i < maxth; i++)
          {
            data[i] = new char[asize];
            p[i] = data[i];
          }
        totsize = asize;
        owner = 1;
      }
  }

  LocalHeap (char * adata, unsigned int asize) throw ()
  {
    int tid = omp_get_thread_num();
    totsize = asize;
    data[tid] = adata;
    owner = 0;
    CleanUp();
  }
  
  ~LocalHeap ()
  {
    if (owner)
      {
        int maxth = omp_get_max_threads();
        for (int i = 0; i < maxth; i++)
          delete [] data[i];
      }
  }
  
  void CleanUp() throw ()
  {
    int tid = omp_get_thread_num();
    p[tid] = data[tid];
    p[tid] += (16 - long(p[tid]) & 15);
  }

  void * GetPointer () throw ()
  {
    int tid = omp_get_thread_num();
    return p[tid];
  }

  void CleanUp (void * addr) throw ()
  {
    int tid = omp_get_thread_num();
    p[tid] = (char*)addr;
  }

  void * Alloc (unsigned int size) throw (LocalHeapOverflow)
  {
    int tid = omp_get_thread_num();

    char * oldp = p[tid];
    
    // 16 byte allignment
    size += (16 - size % 16);
    p[tid] += size;

    if ( (p[tid] - data[tid]) >= int(totsize) )
      ThrowException();

    return oldp;
  }


  template <typename T>
  T * Alloc (unsigned int size) throw (LocalHeapOverflow)
  {
    int tid = omp_get_thread_num();

    char * oldp = p[tid];
    size *= sizeof (T);

    // 16 byte allignment
    size += (16 - size % 16);
    p[tid] += size;

    if ( (p[tid] - data[tid]) >= int(totsize) )
      ThrowException();

    return reinterpret_cast<T*> (oldp);
  }


  ///
  void ThrowException() throw (LocalHeapOverflow);

  /// free memory (dummy function)
  void Free (void * data) throw ()  {;}

  /// available memory on LocalHeap
  int Available () const throw () 
  { 
    int tid = omp_get_thread_num();
    return (totsize - (p[tid]-data[tid])); 
  }
};


/**
   Optimized memory handler.
   Provides static memory for the local heap. The template argument specifies the size in number of chars.
*/
template <int S>
class LocalHeapMem : public LocalHeap
{
  char mem[S];
public:
  LocalHeapMem () throw () : LocalHeap (mem, S) { ; }
};


#endif








/**
   A reset for the heap-pointer of a LocalHeap..
   The constructor stores the heap-pointer, the constructor at the end of the regions resets the heap-pointer.
 */
class HeapReset
{
  LocalHeap & lh;
  void * pointer;
public:
  ///
  HeapReset (LocalHeap & alh) 
    : lh(alh), pointer (alh.GetPointer()) { ; }
    
  ///
  ~HeapReset () 
  {
    lh.CleanUp (pointer);
  }
};



#endif
