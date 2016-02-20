#ifndef FILE_NGS_BITArray
#define FILE_NGS_BITArray

/**************************************************************************/
/* File:   bitarray.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/


namespace ngstd
{

/**
   A compressed array of bools.

   Provides bit-operations and whole array operations.
*/
class BitArray
{
  /// number of bits
  int size;

  /// the data
  unsigned char * data;
public:
  /// empty array
  NGS_DLL_HEADER BitArray ();
  /// array of asize bits
  NGS_DLL_HEADER BitArray (int asize);
  ///
  NGS_DLL_HEADER BitArray (const BitArray & ba2);

  template <typename T>
  INLINE BitArray (std::initializer_list<T> list) 
    : BitArray (list.size())
  {
    Clear();
    int cnt = 0;
    for (auto i = list.begin(); i < list.end(); i++, cnt++)
      if (*i) Set(cnt);
  }

  /// delete data
  NGS_DLL_HEADER ~BitArray ();

  /// Set size, loose values
  NGS_DLL_HEADER void SetSize (int asize);

  /// the size
  int Size () const
  { return size; }


  /// set all bits
  NGS_DLL_HEADER void Set () throw();

  /// clear all bits
  NGS_DLL_HEADER void Clear () throw();

  /// set bit i
  void Set (int i)
  {
#ifdef DEBUG
    if (i < 0 || i >= size)
      throw RangeException ("Bitarray::Set", i, 0, size-1);
#endif
    unsigned char * p = data+Addr(i);
    unsigned char mask = Mask(i);

    /*
#pragma omp atomic
    (*p) |= mask;
    */
    AsAtomic(*p) |= mask;
    
    // data[Addr(i)] |= Mask(i); 
  }

  /// clear bit i
  void Clear (int i)
  { 
#ifdef DEBUG
    if (i < 0 || i >= size)
      throw RangeException ("Bitarray::Clear", i, 0, size-1);
#endif
    data[Addr(i)] &= ~Mask(i); 
  }

  /// check bit i
  bool Test (int i) const
  {
    // return (data[i / CHAR_BIT] & (char(1) << (i % CHAR_BIT) ) ) ? 1 : 0;
    return (data[Addr(i)] & Mask(i)) ? true : false;
  }

  /// set all bits to b
  BitArray & operator= (bool b)
  {
    if (b) Set();
    else   Clear();
    return *this;
  }
  
  /// check bit i
  bool operator[] (int i) const
  {
    return Test(i);
  }

  
  /// invert all bits
  NGS_DLL_HEADER BitArray & Invert ();

  /// logical AND with ba2
  NGS_DLL_HEADER BitArray & And (const BitArray & ba2);

  /// logical OR with ba2
  NGS_DLL_HEADER BitArray & Or (const BitArray & ba2);

  /// copy from ba2
  NGS_DLL_HEADER BitArray & operator= (const BitArray & ba2);

  int NumSet () const;
private:
  ///
  unsigned char Mask (unsigned int i) const
  { return char(1) << (i % CHAR_BIT); }
    
  ///
  int Addr (unsigned int i) const
  { return (i / CHAR_BIT); }


};


/// prints BitArray
NGS_DLL_HEADER ostream & operator<<(ostream & s, const BitArray & ba);

  extern Archive & operator & (Archive & archive, BitArray & ba);

}

#endif
