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
  template <int N, typename T = int>
  class INT
  {
    /// data
    T i[(N>0)?N:1];

  public:
    ///
    INLINE INT () { }

    /// init all
    INLINE INT (T ai1)
    { 
      for (int j = 0; j < N; j++) { i[j] = ai1; }
    }

    /// init i[0], i[1]
    INLINE INT (T ai1, T ai2)
    { i[0] = ai1; i[1] = ai2; }

    /// init i[0], i[1], i[2]
    INLINE INT (T ai1, T ai2, T ai3)
    { i[0] = ai1; i[1] = ai2; i[2] = ai3; }

    /// init i[0], i[1], i[2]
    INLINE INT (T ai1, T ai2, T ai3, T ai4)
    { i[0] = ai1; i[1] = ai2; i[2] = ai3; i[3] = ai4; }

    /*
    /// copy constructor
    INT (const INT<N,T2> & in2)
    { 
      for (int j = 0; j < N; j++) 
	i[j] = in2.i[j]; 
    }
    */

    template <int N2, typename T2>
    INLINE INT (const INT<N2,T2> & in2)
    {
      if (N2 <= N)
        {
          for (int j = 0; j < N2; j++)
            i[j] = in2[j];
          for (int j = N2; j < N; j++)
            i[j] = 0;
        }
      else
        {
          for (int j = 0; j < N; j++)
            i[j] = in2[j];
        }
    }
  
    /// all ints equal ?
    INLINE bool operator== (const INT & in2) const
    { 
      for (int j = 0; j < N; j++) 
	if (i[j] != in2.i[j]) return 0;
      return 1; 
    }

    /// sort integers
    INLINE void Sort ()
    {
      for (int k = 0; k < N; k++)
	for (int l = k+1; l < N; l++)
	  if (i[k] > i[l]) 
	    Swap (i[k], i[l]);
    }

    /// access
    INLINE T & operator[] (int j)
    { return i[j]; }

    /// access
    INLINE const T & operator[] (int j) const
    { return i[j]; }

    INLINE void SetAll (T value)
    {
      for (int j = 0; j < N; j++)
	i[j] = value;
    }

    INLINE INT<N,T> & operator= (T value)
    {
      for (int j = 0; j < N; j++)
	i[j] = value;
      return *this;
    }

    template <typename T2>
    INLINE INT<N,T> & operator= (INT<N,T2> v2)
    {
      for (int j = 0; j < N; j++)
	i[j] = v2[j];
      return *this;
    }
  };

  /// sort 2 integers
  template <>
  INLINE void INT<2>::Sort ()
  {
    if (i[0] > i[1]) Swap (i[0], i[1]);
  }

  /// sort 3 integers
  template <>
  INLINE void INT<3>::Sort ()
  {
    if (i[0] > i[1]) Swap (i[0], i[1]);
    if (i[1] > i[2]) Swap (i[1], i[2]);
    if (i[0] > i[1]) Swap (i[0], i[1]);
  }

  /// Print integers
  template <int N, typename T>
  inline ostream & operator<<(ostream  & s, const INT<N,T> & i2)
  {
    for (int j = 0; j < N; j++)
      s << i2[j] << " ";
    return s;
  }


  template <int N>
  INLINE int HashValue (const INT<N> & ind, int size)
  {
    int sum = 0;
    for (int i = 0; i < N; i++)
      sum += ind[i];
    return sum % size;
  }

  /// hash value of 1 int
  INLINE int HashValue (const INT<1> & ind, int size) 
  {
    return ind[0] % size;
  }

  /// hash value of 2 int
  INLINE int HashValue (const INT<2> & ind, int size) 
  {
    return (113*ind[0]+ind[1]) % size;
  }

  /// hash value of 3 int
  INLINE int HashValue (const INT<3> & ind, int size) 
  {
    return (113*ind[0]+59*ind[1]+ind[2]) % size;
  }


  // using ngstd::max;

  template <int D, typename T>
  INLINE T Max (const INT<D,T> & i)
  {
    if (D == 0) return 0;
    T m = i[0];
    for (int j = 1; j < D; j++)
      if (i[j] > m) m = i[j];
    return m;
  }












  /**
     A hash-table.
     Generic identifiers are mapped to the generic type T.
     An open hashtable. The table is implemented by a DynamicTable.
     Identifiers must provide a HashValue method.
  */
  template <class T_HASH, class T>
  class HashTable
  {
    DynamicTable<T_HASH> hash;
    DynamicTable<T> cont;

  public:
    /// Constructs a hashtable of size bags.
    INLINE HashTable (int size)
      : hash(size), cont(size)
    { ; }
    INLINE ~HashTable () { ; }

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

    /// get value of identifier ahash, exception if unused
    const T & Get (int bnr, int pos) const
    {
      return cont.Get (bnr, pos);
    }

    /// is identifier used ?
    bool Used (const T_HASH & ahash) const
    {
      return (CheckPosition (HashValue (ahash, hash.Size()), ahash) != -1) 
	? 1 : 0;
    }

    /// is identifier used ?
    bool Used (const T_HASH & ahash, int & bnr, int & pos) const
    {
      bnr = HashValue (ahash, hash.Size());
      pos = CheckPosition (bnr, ahash);
      return (pos != -1);
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
    Array<T_HASH,size_t> hash;
    ///
    Array<T,size_t> cont;
    ///
    T_HASH invalid;
  public:
    ///
    ClosedHashTable (int asize)
      : size(asize), hash(asize), cont(asize)
    {
      // hash.SetName ("i2-hashtable, hash");
      // cont.SetName ("i2-hashtable, contents");
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




  template <int N, typename T>
  Archive & operator & (Archive & archive, INT<N,T> & mi)
  {
    for (int i = 0; i < N; i++)
      archive & mi[i];
    return archive;
  }
}

#endif
