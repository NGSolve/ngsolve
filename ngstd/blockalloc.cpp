/**************************************************************************/
/* File:   blockalloc.cpp                                                 */
/* Author: Joachim Schoeberl                                              */
/* Date:   19. Apr. 2000                                                  */
/**************************************************************************/

/* 
   block allocator
*/


#include <ngstd.hpp>

namespace ngstd
{
  using namespace ngstd;


  BlockAllocator :: BlockAllocator (unsigned int asize, unsigned int ablocks)
    : bablocks (0)
  { 
    if (asize < sizeof(void*)) 
      asize = sizeof(void*);
    size = ((asize-1)/sizeof(void*) + 1)*sizeof(void*);
    blocks = ablocks;
    freelist = NULL;
    nels = 0;
  }

  BlockAllocator :: ~BlockAllocator ()
  {
    for (int i = 0; i < bablocks.Size(); i++)
      delete [] bablocks[i];
  }

  void * BlockAllocator :: Alloc2 ()
  {
    // cout << "blockallocator, alloc2 size = " << size << ", blocks = " << blocks << endl;
    //  return new char[size];
    //    if (!freelist)
      {
	char * hcp = new char [size * blocks];
	bablocks.Append (hcp);
	bablocks.Last() = hcp;

	for (unsigned int i = 0; i < blocks-1; i++)
	  *(void**)&(hcp[i * size]) = &(hcp[ (i+1) * size]);
	*(void**)&(hcp[(blocks-1)*size]) = NULL;
	freelist = hcp;
      }
      /*
    void * p = freelist;
    freelist = *(void**)freelist;
    return p;
      */
	  return (void*)freelist;
  }

  /*
  void BlockAllocator :: Free (void * p)
  {
    //  delete [] p;
    *(void**)p = freelist;
    freelist = p;
  }
  */
}

