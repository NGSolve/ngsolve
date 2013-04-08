/**************************************************************************/
/* File:   bitarray.cc                                                    */
/* Autho: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   data type BitArray
*/

#include <ngstd.hpp>

namespace ngstd
{
  using namespace ngstd;

  BitArray :: BitArray ()
  {
    size = 0;
    data = NULL;
  }

  BitArray :: BitArray (int asize)
  {
    size = 0;
    data = NULL;
    SetSize (asize);
  }

  BitArray :: BitArray (const BitArray & ba2)
  {
    size = 0;
    data = NULL;
    (*this) = ba2;
  }

  BitArray :: ~BitArray ()
  {
    if (data) delete [] data;
  }

  void BitArray :: SetSize (int asize)
  {
    if (size == asize) return;
    if (data) delete [] data;

    size = asize;
    data = new unsigned char [Addr (size)+1];
  }

  void BitArray :: Set ()
  {
    int i;
    if (!size) return;
    for (i = 0; i <= Addr (size); i++)
      data[i] = UCHAR_MAX;
  }

  void BitArray :: Clear ()
  {
    int i;
    if (!size) return;
    for (i = 0; i <= Addr (size); i++)
      data[i] = 0;
  }



  void BitArray :: Invert ()
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] ^= 255;
  }

  void BitArray :: And (const BitArray & ba2)
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] &= ba2.data[i];
  }


  void BitArray :: Or (const BitArray & ba2)
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] |= ba2.data[i];
  }


  BitArray & BitArray :: operator= (const BitArray & ba2)
  {
    SetSize (ba2.Size());
    if (!size) 
      return *this;
    for (int i = 0; i <= Addr (size); i++)
      data[i] = ba2.data[i];
    return *this;
  }

  ostream & operator<<(ostream & s, const BitArray & ba)
  {
    int n = ba.Size();
    for (int i = 0; i < n; i++)
      {
	if (i % 50 == 0) s << i << ": ";
	s << int(ba[i]);
	if (i % 50 == 49) s << "\n";
      }
    s << flush;
    return s;
  }


}
