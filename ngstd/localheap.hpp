#ifndef FILE_LOCALHEAP
#define FILE_LOCALHEAP

/**************************************************************************/
/* File:   localheap.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2000                                                  */
/**************************************************************************/


namespace ngstd
{

  /**
     Exception on heap overflow.
     Thrown by allocation on LocalHeap.
  */
  class NGS_DLL_HEADER LocalHeapOverflow : public Exception
  {
  public:
    LocalHeapOverflow (size_t size);
    virtual ~LocalHeapOverflow ();
  };
 


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
    size_t totsize;
    bool owner;
    const char * name;

  public:
    /// Allocate one block of size asize.
    NGS_DLL_HEADER LocalHeap (size_t asize, const char * aname);

    /// Use provided memory for the LocalHeap
    LocalHeap (char * adata, size_t asize, const char  * aname) throw ()
    {
      totsize = asize;
      data = adata;
      owner = 0;
      // p = data;
      name = aname;
      CleanUp();
    }

    /// Use provided memory for the LocalHeap
    LocalHeap (const LocalHeap & lh2)
      : data(lh2.data), p(lh2.p), totsize(lh2.totsize), owner(false)
    { ; }

  
    /// free memory
    ~LocalHeap ()
    {
      if (owner)
	delete [] data;
    }
  
    /// delete all memory on local heap
    void CleanUp() throw ()
    {
      p = data;
      // p += (16 - (long(p) & 15) );
      p += (16 - (size_t(p) & 15) );
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
    void * Alloc (size_t size) throw (LocalHeapOverflow)
    {
      char * oldp = p;
    
      // 16 byte allignment
      size += (16 - size % 16);
      p += size;

      if ( size_t(p - data) >= totsize )
	ThrowException();

      return oldp;
    }

    /// allocates size objects of type T on local heap
    template <typename T>
    T * Alloc (size_t size) throw (LocalHeapOverflow)
    {
      char * oldp = p;
      size *= sizeof (T);

      // 16 byte allignment
      size += (16 - size % 16);
      p += size;

      if ( size_t(p - data) >= totsize )
	ThrowException();

      return reinterpret_cast<T*> (oldp);
    }

  private:
    ///
    NGS_DLL_HEADER void ThrowException() throw (LocalHeapOverflow);

  public:
    /// free memory (dummy function)
    void Free (void * data) throw () 
    {
      ;
    }

    /// available memory on LocalHeap
    size_t Available () const throw () { return (totsize - (p-data)); }

    /// Split free memory on heap into pieces for each openmp-thread
    LocalHeap Split () const
    {
#ifdef _OPENMP
      int pieces = omp_get_num_threads();
      int i = omp_get_thread_num();
#else
      int pieces = 1;
      int i = 0;
#endif
      size_t freemem = totsize - (p - data);
      size_t size_of_piece = freemem / pieces;
      return LocalHeap (p + i * size_of_piece, size_of_piece, name);
    }

    void ClearValues ()
    {
      for (size_t i = 0; i < totsize; i++) data[i] = 47;
    }

    size_t UsedSize ()
    {
      for (size_t i = totsize-1; i >= 0; i--)
        if (data[i] != 47) return i;
      return 0;
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
    LocalHeapMem (const char * aname) throw () : LocalHeap (mem, S, aname) { ; }
  };








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

}


inline void * operator new (size_t size, ngstd::LocalHeap & lh)  
{
  return lh.Alloc(size);
}

inline void operator delete (void * p, ngstd::LocalHeap & lh)  
{
  ; 
}



#endif
