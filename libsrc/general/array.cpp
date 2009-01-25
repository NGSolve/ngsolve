#ifndef FILE_NGSTD_ArrayCPP
#define FILE_NGSTD_ArrayCPP
// necessary for SGI ????

/**************************************************************************/
/* File:   array.cpp                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type Array
*/

#include <mystdlib.h>
#include <myadt.hpp>
#include <assert.h>


namespace netgen
{
  //using namespace netgen;

#ifdef NONE  
  void BASE_Array :: ReSize (int minsize, int elementsize)
  {
    cout << "resize, minsize = " << minsize << endl;

    if (inc == -1)
      throw Exception ("Try to resize fixed size array");

    
    void * p;
    int nsize = (inc) ? allocsize + inc : 2 * allocsize;
    if (nsize < minsize) nsize = minsize;

    if (data)
      {
	p = new char [nsize * elementsize];
	
	int mins = (nsize < actsize) ? nsize : actsize; 
	memcpy (p, data, mins * elementsize);
	
	delete [] static_cast<char*> (data);
	data = p;
      }
    else
      {
	data = new char[nsize * elementsize];
      }
    
    allocsize = nsize;
    cout << "resize done" << endl;
  }
  
  
  
  void BASE_Array :: RangeCheck (int i) const
  {
    if (i < 0 || i >= actsize)
      throw ArrayRangeException ();
  }
  
  void BASE_Array :: CheckNonEmpty () const
  {
    if (!actsize)
      {
	throw Exception ("Array should not be empty");
	//      cerr << "Array souldn't be empty";
      }
  }
#endif
}
#endif

