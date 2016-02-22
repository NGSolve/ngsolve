/**************************************************************************/
/* File:   optmem.cpp                                                     */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

/* 
   Abstract data type Array
*/


#include <mystdlib.h>
#include <myadt.hpp>

namespace netgen
{
  static mutex block_allocator_mutex;

  BlockAllocator :: BlockAllocator (unsigned asize, unsigned ablocks)
    : bablocks (0)
  {
    if (asize < sizeof(void*))
      asize = sizeof(void*);
    size = asize;
    blocks = ablocks;
    freelist = NULL;
  }

  BlockAllocator :: ~BlockAllocator ()
  {
    // cout << "****************** delete BlockAllocator " << endl;
    for (int i = 0; i < bablocks.Size(); i++)
      delete [] bablocks[i];
    bablocks.SetSize(0);
  }

  void * BlockAllocator :: Alloc ()
  {
    void * p;
    {
      lock_guard<mutex> guard(block_allocator_mutex); 
      //  return new char[size];
      if (!freelist)
        {
          // cout << "freelist = " << freelist << endl;
          // cout << "BlockAlloc: " << size*blocks << endl;
          char * hcp = new char [size * blocks];
          bablocks.Append (hcp);
          bablocks.Last() = hcp;
          for (unsigned i = 0; i < blocks-1; i++)
            *(void**)&(hcp[i * size]) = &(hcp[ (i+1) * size]);
          *(void**)&(hcp[(blocks-1)*size]) = NULL;
          freelist = hcp;
        }
      
      p = freelist;
      freelist = *(void**)freelist;
    }
    return p;
  }

  void BlockAllocator :: Free (void * p)
  {
    {
      lock_guard<mutex> guard(block_allocator_mutex); 
      if (bablocks.Size())
        {
          *(void**)p = freelist;
          freelist = p;
        }
    }
  }

}
