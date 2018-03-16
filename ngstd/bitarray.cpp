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
  /*
  BitArray :: BitArray ()
  {
    size = 0;
    data = NULL;
  }
  */
  BitArray :: BitArray (size_t asize)
  {
    size = 0;
    data = NULL;
    SetSize (asize);
  }

  BitArray :: BitArray (size_t asize, LocalHeap & lh)
  {
    size = asize;
    data = new (lh) unsigned char [Addr (size)+1];
    owns_data = false;
  }
  
  BitArray :: BitArray (const BitArray & ba2)
  {
    size = 0;
    data = NULL;
    (*this) = ba2;
  }

  /*
  BitArray :: ~BitArray ()
  {
    if (owns_data)
      delete [] data;
  }
  */
  
  void BitArray :: SetSize (size_t asize)
  {
    if (size == asize) return;
    if (owns_data) delete [] data;

    size = asize;
    data = new unsigned char [Addr (size)+1];
  }

  BitArray & BitArray :: Set () throw()
  {
    if (!size) return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] = UCHAR_MAX;
    return *this;
  }

  BitArray & BitArray :: Clear () throw()
  {
    if (!size) return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] = 0;
    return *this;
  }



  BitArray & BitArray :: Invert ()
  {
    if (!size) return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] ^= 255;
    return *this;
  }

  BitArray & BitArray :: And (const BitArray & ba2)
  {
    if (!size) return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] &= ba2.data[i];
        return *this;
  }


  BitArray & BitArray :: Or (const BitArray & ba2)
  {
    if (!size) return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] |= ba2.data[i];
    return *this;
  }


  BitArray & BitArray :: operator= (const BitArray & ba2)
  {
    SetSize (ba2.Size());
    if (!size) 
      return *this;
    for (size_t i = 0; i <= Addr (size); i++)
      data[i] = ba2.data[i];
    return *this;
  }

  ostream & operator<<(ostream & s, const BitArray & ba)
  {
    size_t n = ba.Size();
    for (size_t i = 0; i < n; i++)
      {
	if (i % 50 == 0) s << i << ": ";
	s << int(ba[i]);
	if (i % 50 == 49) s << "\n";
      }
    s << flush;
    return s;
  }

  size_t BitArray :: NumSet () const
  {
    size_t cnt = 0;
    for (size_t i = 0; i < Size(); i++)
      if (Test(i)) cnt++;
    return cnt;
  }

  Archive & operator & (Archive & archive, BitArray & ba)
  {
    if (archive.Output())
      {
        archive << ba.Size();
        for (size_t i = 0; i < ba.Size(); i++)
          archive << ba[i];
      }
    else
      {
        int size;
        archive & size;
        ba.SetSize (size);
        ba.Clear();
        for (size_t i = 0; i < size; i++)
          {
            bool b;
            archive & b;
            if (b) ba.Set(i);
          }
      }
    return archive;
  }

  
}
