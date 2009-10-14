#ifndef FILE_NGSTD_HASHTABLE
#define FILE_NGSTD_HASHTABLE

/**************************************************************************/
/* File:   hashtable.hpp                                                  */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/


namespace ngstd
{


  /// N integers
  template <int N>
  class INT
  {
    /// data
    int i[N];

  public:
    ///
    INT () { }

    /// init i[0]
    INT (int ai1)
    { i[0] = ai1; }

    /// init i[0], i[1]
    INT (int ai1, int ai2)
    { i[0] = ai1; i[1] = ai2; }

    /// init i[0], i[1], i[2]
    INT (int ai1, int ai2, int ai3)
    { i[0] = ai1; i[1] = ai2; i[2] = ai3; }

    /// init i[0], i[1], i[2]
    INT (int ai1, int ai2, int ai3, int ai4)
    { i[0] = ai1; i[1] = ai2; i[2] = ai3; i[3] = ai4; }


    /// copy constructor
    INT (const INT & in2)
    { 
      for (int j = 0; j < N; j++) 
	i[j] = in2.i[j]; 
    }
  
    /// all ints equal ?
    bool operator== (const INT & in2) const
    { 
      for (int j = 0; j < N; j++) 
	if (i[j] != in2.i[j]) return 0;
      return 1; 
    }

    /// sort integers
    void Sort ()
    {
      for (int k = 0; k < N; k++)
	for (int l = k+1; l < N; l++)
	  if (i[k] > i[l]) 
	    Swap (i[k], i[l]);
    }

    /// access
    int & operator[] (int j)
    { return i[j]; }

    /// access
    const int & operator[] (int j) const
    { return i[j]; }

    void SetAll (int value)
    {
      for (int j = 0; j < N; j++)
	i[j] = value;
    }

    INT<N> & operator= (int value)
    {
      for (int j = 0; j < N; j++)
	i[j] = value;
      return *this;
    }

    INT<N> & operator= (INT<N> v2)
    {
      for (int j = 0; j < N; j++)
	i[j] = v2.i[j];
      return *this;
    }
  };

  /// sort 2 integers
  template <>
  inline void INT<2>::Sort ()
  {
    if (i[0] > i[1]) swap (i[0], i[1]);
  }

  /// sort 3 integers
  template <>
  inline void INT<3>::Sort ()
  {
    if (i[0] > i[1]) swap (i[0], i[1]);
    if (i[1] > i[2]) swap (i[1], i[2]);
    if (i[0] > i[1]) swap (i[0], i[1]);
  }

  /// Print integers
  template <int N>
  inline ostream & operator<<(ostream  & s, const INT<N> & i2)
  {
    for (int j = 0; j < N; j++)
      s << i2[j] << " ";
    return s;
  }


  template <int N>
  inline int HashValue (const INT<N> & ind, int size)
  {
    int sum = 0;
    for (int i = 0; i < N; i++)
      sum += ind[i];
    return sum % size;
  }

  /// hash value of 1 int
  inline int HashValue (const INT<1> & ind, int size) 
  {
    return ind[0] % size;
  }

  /// hash value of 2 int
  inline int HashValue (const INT<2> & ind, int size) 
  {
    return (113*ind[0]+ind[1]) % size;
  }

  /// hash value of 3 int
  inline int HashValue (const INT<3> & ind, int size) 
  {
    return (113*ind[0]+59*ind[1]+ind[2]) % size;
  }













  /**
     A hash-table.
     Generic identifiers are mapped to the generic type T.
     An open hashtable. The table is implemented by a \Ref{DynamicTable}.
     Identifiers must provide a \Ref{HashValue} function.
  */
  template <class T_HASH, class T>
  class HashTable
  {
    DynamicTable<T_HASH> hash;
    DynamicTable<T> cont;

  public:
    /// Constructs a hashtable of size bags.
    HashTable (int size)
      : hash(size), cont(size)
    {
      ;
    }

    /// Sets identifier ahash to value acont
    void Set (const T_HASH & ahash, const T & acont)
    {
      int bnr = HashValue (ahash, hash.Size());
      int pos = CheckPosition (bnr, ahash);
      if (pos != -1)
	cont.Set (bnr, pos, acont);
      else
	{
	  hash.Add (bnr, ahash);
	  cont.Add (bnr, acont);
	}        
    }

    /// get value of identifier ahash, exception if unused
    const T & Get (const T_HASH & ahash) const
    {
      int bnr = HashValue (ahash, hash.Size());
      int pos = Position (bnr, ahash);
      return cont.Get (bnr, pos);
    }

    /// is identifier used ?
    bool Used (const T_HASH & ahash) const
    {
      return (CheckPosition (HashValue (ahash, hash.Size()), ahash) != -1) 
	? 1 : 0;
    }

    /// number of hash entries
    int Size () const
    {
      return hash.Size();
    }

    /// size of hash entry
    int EntrySize (int bnr) const
    {
      return hash[bnr].Size();
    }

    /// get identifier and value of entry bnr, position colnr
    void GetData (int bnr, int colnr, T_HASH & ahash, T & acont)
    {
      ahash = hash[bnr][colnr];
      acont = cont[bnr][colnr];
    }

    /// set identifier and value of entry bnr, position colnr
    void SetData (int bnr, int colnr, const T_HASH & ahash, const T & acont)
    {
      hash[bnr][colnr] = ahash;
      cont[bnr][colnr] = acont;
    }    

    /// returns position of index. returns -1 on unused
    int CheckPosition (int bnr, const T_HASH & ind) const
    {
      for (int i = 0; i < hash[bnr].Size(); i++)
	if (hash[bnr][i] == ind)
	  return i;
      return -1;
    }

    /// returns position of index. exception on unused
    int Position (int bnr, const T_HASH & ind) const
    {
      for (int i = 0; i < hash[bnr].Size(); i++)
	if (hash[bnr][i] == ind)
	  return i;
      throw Exception ("Ask for unsused hash-value");
    }
  };






  /**
     A closed hash-table.
     All information is stored in one fixed array.
     The array should be allocated with the double size of the expected number of entries.
  */
  template <class T_HASH, class T>
  class ClosedHashTable
  {
  protected:
    ///
    int size;
    ///
    DynamicMem<T_HASH> hash;
    ///
    DynamicMem<T> cont;
    ///
    T_HASH invalid;
  public:
    ///
    ClosedHashTable (int asize)
      : size(asize), hash(asize), cont(asize)
    {
      hash.SetName ("i2-hashtable, hash");
      cont.SetName ("i2-hashtable, contents");
      invalid.SetAll (-1);
      for (int i = 0; i < size; i++)
	hash[i] = invalid;
    }

    /// 
    int Size() const
    {
      return size;
    }

    /// is position used
    bool UsedPos (int pos) const
    {
      return ! (hash[pos] == invalid); 
    }

    /// number of used elements
    int UsedElements () const
    {
      int cnt = 0;
      for (int i = 0; i < size; i++)
	if (hash[i] != invalid)
	  cnt++;
      return cnt;
    }

    int Position (const T_HASH ind) const
    {
      int i = HashValue(ind, size);
      while (1)
	{
	  if (hash[i] == ind) return i;
	  if (hash[i] == invalid) return -1;
	  i++;
	  if (i >= size) i = 0;
	}
    }
    // returns 1, if new postion is created
    int PositionCreate (const T_HASH ind, int & apos)
    {
      int i = HashValue (ind, size);

      while (1)
	{
	  if (hash[i] == invalid)
	    { 
	      hash[i] = ind; 
	      apos = i; 
	      return 1;
	    }
	  if (hash[i] == ind) 
	    { 
	      apos = i; 
	      return 0; 
	    }
	  i++;
	  if (i >= size) i = 0;
	}
    }


    ///
    void Set (const T_HASH & ahash, const T & acont)
    {
      int pos;
      PositionCreate (ahash, pos);
      hash[pos] = ahash;
      cont[pos] = acont;
    }
    ///
    const T & Get (const T_HASH & ahash) const
    {
      int pos = Position (ahash);
      return cont[pos];
    }

    ///
    bool Used (const T_HASH & ahash) const
    {
      int pos = Position (ahash);
      return (pos != -1);
    }

    void SetData (int pos, const T_HASH & ahash, const T & acont)
    {
      hash[pos] = ahash;
      cont[pos] = acont;
    }

    void GetData (int pos, T_HASH & ahash, T & acont) const
    {
      ahash = hash[pos];
      acont = cont[pos];
    }
  
    void SetData (int pos, const T & acont)
    {
      cont[pos] = acont;
    }

    void GetData (int pos, T & acont) const
    {
      acont = cont[pos];
    }

    void SetSize (int asize)
    {
      size = asize;
      hash.Alloc(size);
      cont.Alloc(size);
      for (int i = 0; i < size; i++)
	hash[i] = invalid;
    }

    void SetName (const char * aname)
    {
      cont.SetName(aname);
      hash.SetName(aname);
    }
  };

}


#endif
