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
class NGS_DLL_HEADER BitArray
{
  /// number of bits
  int size;

  /// the data
  unsigned char * data;
public:
  /// empty array
  BitArray ();
  /// array of asize bits
  BitArray (int asize);
  ///
  BitArray (const BitArray & ba2);
  /// delete data
  ~BitArray ();

  /// Set size, loose values
  void SetSize (int asize);

  /// the size
  int Size () const
  { return size; }


  /// set all bits
  void Set ();

  /// clear all bits
  void Clear ();

  /// set bit i
  void Set (int i)
  { 
    unsigned char * p = data+Addr(i);
    unsigned char mask = Mask(i);

#pragma omp atomic
    (*p) |= mask;
    // data[Addr(i)] |= Mask(i); 
  }

  /// clear bit i
  void Clear (int i)
  { data[Addr(i)] &= ~Mask(i); }

  /// check bit i
  inline bool Test (int i) const
  {
    return (data[i / CHAR_BIT] & (char(1) << (i % CHAR_BIT) ) ) ? 1 : 0;
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
    return (data[i / CHAR_BIT] & (char(1) << (i % CHAR_BIT) ) ) ? 1 : 0;
  }

  
  /// invert all bits
  void Invert ();

  /// logical AND with ba2
  void And (const BitArray & ba2);

  /// logical OR with ba2
  void Or (const BitArray & ba2);

  /// copy from ba2
  BitArray & operator= (const BitArray & ba2);

#ifdef PARALLEL
  const unsigned char * Data () const
  { return data; }

  unsigned char * Data ()
  { return data; }
#endif

private:
  ///
  unsigned char Mask (int i) const
  { return char(1) << (i % CHAR_BIT); }
    
  ///
  inline int Addr (int i) const
  { return (i / CHAR_BIT); }


};


/// prints BitArray
ostream & operator<<(ostream & s, const BitArray & ba);

}

#endif
