#ifndef FILE_BLOCKALLOC
#define FILE_BLOCKALLOC

/**************************************************************************/
/* File:   blockalloc.hpp                                                 */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2000                                                  */
/**************************************************************************/

#include <ngstd.hpp>

namespace ngstd
{

/**
   Optimized memory handler.
   Memory handler allocates many objects at once.
   Maintains free list of deleted objects
 */
class BlockAllocator
{
  /// size of data
  unsigned int size;
  /// number of blocks allocated at once
  unsigned int blocks;
  /// single linked list of free elements
  void * freelist;
  /// pointers to blocks
  Array<char*> bablocks;
  ///
  int nels;
public:
  /// Create BlockAllocator for elements of size asize
  NGS_DLL_HEADER BlockAllocator (unsigned int asize, unsigned int ablocks = 100);
  /// Delete all memeory
  NGS_DLL_HEADER ~BlockAllocator ();

  /// Return pointer to new element
  NGS_DLL_HEADER void * Alloc ();
  /*
  {
    nels++;
    if (!freelist)
      Alloc2 ();

    void * p = freelist;
    freelist = *(void**)freelist;
    return p;
  }
  */


  /// Send memory to free-list
  NGS_DLL_HEADER void Free (void * p);
  /*
  {
    nels--;
    *(void**)p = freelist;
    freelist = p;
  }
  */

  /// number of allocated elements
  int NumElements () { return nels; }

  NGS_DLL_HEADER void Print (ostream * ost) const;

  const MemoryTracer& GetMemoryTracer() const { return mt; }

  void StartMemoryTracing() const
  {
    mt.Alloc(bablocks.Size() * size * blocks);
  }

private:
  NGS_DLL_HEADER void * Alloc2 ();
  MemoryTracer mt;
};


}

INLINE void * operator new (size_t /* size */, ngstd::BlockAllocator & ball)  
{
  return ball.Alloc();
}

INLINE void operator delete (void * p, ngstd::BlockAllocator & ball)  
{
  ball.Free (p);
}



#endif
