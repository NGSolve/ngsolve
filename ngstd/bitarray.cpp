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
    delete [] data;
  }

  void BitArray :: SetSize (int asize)
  {
    if (size == asize) return;
    if (data) delete [] data;

    size = asize;
    data = new unsigned char [Addr (size)+1];
  }

  void BitArray :: Set () throw()
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] = UCHAR_MAX;
  }

  void BitArray :: Clear () throw()
  {
    if (!size) return;
    for (int i = 0; i <= Addr (size); i++)
      data[i] = 0;
  }



  BitArray & BitArray :: Invert ()
  {
    if (!size) return *this;
    for (int i = 0; i <= Addr (size); i++)
      data[i] ^= 255;
    return *this;
  }

  BitArray & BitArray :: And (const BitArray & ba2)
  {
    if (!size) return *this;
    for (int i = 0; i <= Addr (size); i++)
      data[i] &= ba2.data[i];
        return *this;
  }


  BitArray & BitArray :: Or (const BitArray & ba2)
  {
    if (!size) return *this;
    for (int i = 0; i <= Addr (size); i++)
      data[i] |= ba2.data[i];
    return *this;
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

  int BitArray :: NumSet () const
  {
    int cnt = 0;
    for (int i = 0; i < Size(); i++)
      if (Test(i)) cnt++;
    return cnt;
  }

  Archive & operator & (Archive & archive, BitArray & ba)
  {
    if (archive.Output())
      {
        archive << ba.Size();
        for (int i = 0; i < ba.Size(); i++)
          archive << ba[i];
      }
    else
      {
        int size;
        archive & size;
        ba.SetSize (size);
        ba.Clear();
        for (int i = 0; i < size; i++)
          {
            bool b;
            archive & b;
            if (b) ba.Set(i);
          }
      }
    return archive;
  }

  
}
