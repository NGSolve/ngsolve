#ifndef FILE_LOCALHEAP
#define FILE_LOCALHEAP

/**************************************************************************/
/* File:   localheap.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2000                                                  */
/**************************************************************************/


namespace ngstd
{


  class Allocator
  {
  public:
    virtual ~Allocator() {}
    virtual void * Alloc (size_t size)
    {
      return new char[size];
    }
    virtual void Delete(void* p)
    {
      delete (char*) p;
    }
    virtual void ArrayDelete(void* p)
    {
        delete [] (char*) p;
    }
  };
  static Allocator global_alloc;

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
  class LocalHeap : public Allocator
  {
    char * data;
    char * next;
    char * p;
    size_t totsize;
  public:
    bool owner;
    const char * name;

#if defined(__MIC__) || defined (__AVX512F__)
    enum { ALIGN = 64 };
#else
    enum { ALIGN = 32 };
#endif  

  public:
    /// Allocate one block of size asize.
    NGS_DLL_HEADER LocalHeap (size_t asize, 
                              const char * aname = "noname",
                              bool mult_by_threads = false);

    /// Use provided memory for the LocalHeap
    INLINE LocalHeap (char * adata, size_t asize, 
                      const char  * aname = "noname") throw ()
    {
      totsize = asize;
      data = adata;
      next = data + totsize;
      owner = 0;
      // p = data;
      name = aname;
      CleanUp();
    }

    /*
    /// Use provided memory for the LocalHeap
    INLINE LocalHeap (const LocalHeap & lh2)
      : data(lh2.data), p(lh2.p), totsize(lh2.totsize), owner(false),
        name(lh2.name)
    {
      next = data + totsize;
    }
    */
    INLINE LocalHeap (const LocalHeap & lh2) = delete;

    INLINE LocalHeap (LocalHeap && lh2)
      : data(lh2.data), p(lh2.p), totsize(lh2.totsize), owner(lh2.owner),
        name(lh2.name)
    {
      next = data + totsize;
      lh2.owner = false;
    }
    
    INLINE LocalHeap Borrow() 
    {
      return LocalHeap (p, Available());
    }


    INLINE LocalHeap & operator= (LocalHeap && lh2)
    {
      if (owner)
        delete [] data;
      
      data = lh2.data;
      p = lh2.p;
      totsize = lh2.totsize;
      owner = lh2.owner;
      name = lh2.name;

      next = data + totsize;
      lh2.owner = false;
      return *this;
    }

    INLINE LocalHeap ()
      : data(nullptr), next(nullptr), p(nullptr), totsize(0), owner(false) { ; }

  
    /// free memory
    INLINE virtual ~LocalHeap ()
    {
      if (owner)
	delete [] data;
    }
  
    /// delete all memory on local heap
    INLINE void CleanUp() throw ()
    {
      p = data;
      // p += (16 - (long(p) & 15) );
      p += (ALIGN - (size_t(p) & (ALIGN-1) ) );
    }

    /// returns heap-pointer
    INLINE void * GetPointer () throw ()
    {
      return p;
    }

    /// deletes memory back to heap-pointer
    INLINE void CleanUp (void * addr) throw ()
    {
      p = (char*)addr;
    }

    /// allocates size bytes of memory from local heap
    INLINE void * Alloc (size_t size) final // throw (LocalHeapOverflow)
    {
      char * oldp = p;
    
      // 16 byte alignment
      size += (ALIGN - size % ALIGN);
      p += size;

      // if ( size_t(p - data) >= totsize )
#ifndef FULLSPEED
      if (likely(p >= next))
        ThrowException();
#endif
      return oldp;
    }

    /// allocates size objects of type T on local heap
    template <typename T>
    INLINE T * Alloc (size_t size) // throw (LocalHeapOverflow)
    {
      char * oldp = p;
      size *= sizeof (T);

      // 16 byte alignment
      size += (ALIGN - size % ALIGN);
      p += size;

#ifndef FULLSPEED
      if (likely(p >= next))
	ThrowException();
#endif

      return reinterpret_cast<T*> (oldp);
    }

    virtual void Delete(void* p) {}

    virtual void ArrayDelete(void* p) {}
  private: 
    ///
#ifndef __CUDA_ARCH__
    [[noreturn]] NGS_DLL_HEADER void ThrowException(); 
#else
    INLINE void ThrowException() { ; }
#endif

  public:
    /// free memory (dummy function)
    INLINE void Free (void * data) throw () 
    {
      ;
    }

    /// available memory on LocalHeap
    INLINE size_t Available () const throw () { return (totsize - (p-data)); }

    /// Split free memory on heap into pieces for each thread
    NGS_DLL_HEADER LocalHeap Split () const;

    /// Split free memory on heap into pieces
    INLINE LocalHeap Split (int partnr, int nparts) const
    {
      int pieces = nparts;
      int i = partnr;

      size_t freemem = totsize - (p - data);
      size_t size_of_piece = freemem / pieces;
      return LocalHeap (p + i * size_of_piece, size_of_piece, name);
    }


    INLINE void ClearValues ()
    {
      for (size_t i = 0; i < totsize; i++) data[i] = 47;
    }

    INLINE size_t UsedSize ()
    {
      for (size_t i = totsize-1; i != 0; i--)
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
    INLINE LocalHeapMem (const char * aname) throw () : LocalHeap (mem, S, aname) { ; }
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
    INLINE HeapReset (LocalHeap & alh) 
      : lh(alh), pointer (alh.GetPointer()) { ; }
    
    ///
    INLINE ~HeapReset () 
    {
      lh.CleanUp (pointer);
    }
  };

}



INLINE void * operator new (size_t size, ngstd::Allocator & alloc)  
{
  return alloc.Alloc(size);
}

INLINE void * operator new [] (size_t size, ngstd::Allocator & alloc)  
{
  return alloc.Alloc(size);
}


INLINE void operator delete (void * p, ngstd::Allocator & lh)  
{
  lh.Delete(p);
}

INLINE void operator delete [] (void * p, ngstd::Allocator & lh)  
{
  lh.ArrayDelete(p);
}



#endif
