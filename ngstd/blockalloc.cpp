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
    mt.Free(bablocks.Size() * size * blocks);
    for (int i = 0; i < bablocks.Size(); i++)
      delete [] bablocks[i];
  }

  void * BlockAllocator :: Alloc ()
  {
    nels++;
    if (!freelist)
      Alloc2 ();

    void * p = freelist;
    freelist = *(void**)freelist;
    return p;
  }


  void * BlockAllocator :: Alloc2 ()
  {
    // cout << "blockallocator, alloc2 size = " << size << ", blocks = " << blocks << endl;
    //  return new char[size];
    //    if (!freelist)
      {
	char * hcp = new char [size * blocks];
        mt.Alloc(size * blocks);
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


  void BlockAllocator :: Free (void * p)
  {
    nels--;
    *(void**)p = freelist;
    freelist = p;
  }


  void BlockAllocator :: Print (ostream * ost) const
  {
    int cnt = -1;
    void * p = freelist;

    while (p && cnt++ < 10)
      {
	*ost << "el " << cnt << " = " << p << endl;
	void * newp = *(void**)p;
	if (p == newp) 
	  {
	    cerr << "defect freelist, p = newp" << endl;
	    exit(1);
	  }
	p = newp;
      }
  }

}

