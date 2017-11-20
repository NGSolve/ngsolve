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

    /// init i[0], i[1], i[2]
    INLINE INT (T ai1, T ai2, T ai3, T ai4, T ai5)
    { i[0] = ai1; i[1] = ai2; i[2] = ai3; i[3] = ai4; i[4] = ai5;}

    /// init i[0], i[1], i[2]
    INLINE INT (T ai1, T ai2, T ai3, T ai4, T ai5, T ai6, T ai7, T ai8, T ai9)
    { i[0] = ai1; i[1] = ai2; i[2] = ai3; i[3] = ai4; i[4] = ai5; i[5] = ai6; i[6] = ai7; i[7] = ai8; i[8] = ai9; }

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

    template <typename T2>
    INLINE INT (const BaseArrayObject<T2> & ao)
    {
      for (int j = 0; j < N; j++)
        i[j] = ao.Spec()[j];
    }
    
    INLINE size_t Size() const { return N; }
    /// all ints equal ?
    INLINE bool operator== (const INT & in2) const
    { 
      for (int j = 0; j < N; j++) 
	if (i[j] != in2.i[j]) return 0;
      return 1; 
    }

    /// any ints unequal ?
    INLINE bool operator!= (const INT & in2) const
    {
      for (int j = 0; j < N; j++)
        if (i[j] != in2.i[j]) return 1;
      return 0;
    }

    /// sort integers
    INLINE INT & Sort () & 
    {
      for (int k = 0; k < N; k++)
	for (int l = k+1; l < N; l++)
	  if (i[k] > i[l]) 
	    Swap (i[k], i[l]);
      return *this;
    }

    INLINE INT Sort () &&
    {
      for (int k = 0; k < N; k++)
	for (int l = k+1; l < N; l++)
	  if (i[k] > i[l]) 
	    Swap (i[k], i[l]);
      return *this;
    }

    /// access
    INLINE T & operator[] (int j)
    { return i[j]; }

    /// access
    INLINE const T & operator[] (int j) const
    { return i[j]; }

    /*
    INLINE void SetAll (T value)
    {
      for (int j = 0; j < N; j++)
	i[j] = value;
    }
    */

    operator FlatArray<T> () { return FlatArray<T> (N, &i[0]); } 

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
  INLINE INT<2> & INT<2>::Sort () & 
  {
    if (i[0] > i[1]) Swap (i[0], i[1]);
    return *this;
  }

  template <>
  INLINE INT<2> INT<2>::Sort () &&
  {
    if (i[0] > i[1]) Swap (i[0], i[1]);
    return *this;
  }

  /// sort 3 integers
  template <>
  INLINE INT<3> INT<3>::Sort () &&
  {
    if (i[0] > i[1]) Swap (i[0], i[1]);
    if (i[1] > i[2]) Swap (i[1], i[2]);
    if (i[0] > i[1]) Swap (i[0], i[1]);
    return *this;
  }

  /// Print integers
  template <int N, typename T>
  inline ostream & operator<<(ostream  & s, const INT<N,T> & i2)
  {
    for (int j = 0; j < N; j++)
      s << (int) i2[j] << " ";
    return s;
  }
  
  template <int N, typename T>
  auto begin(const INT<N,T> & ind)
  {
    return AOWrapperIterator<INT<N,T>> (ind, 0);
  }

  template <int N, typename T>
  auto end(const INT<N,T> & ind)
  {
    return AOWrapperIterator<INT<N,T>> (ind, N);    
  }





  
  template <int N, typename TI>
  INLINE size_t HashValue (const INT<N,TI> & ind, size_t size)
  {
    INT<N,size_t> lind = ind;    
    size_t sum = 0;
    for (int i = 0; i < N; i++)
      sum += lind[i];
    return sum % size;
  }

  /// hash value of 1 int
  template <typename TI>
  INLINE size_t HashValue (const INT<1,TI> & ind, size_t size) 
  {
    return ind[0] % size;
  }

  /// hash value of 2 int
  template <typename TI>  
  INLINE size_t HashValue (const INT<2,TI> & ind, size_t size) 
  {
    INT<2,size_t> lind = ind;
    return (113*lind[0]+lind[1]) % size;
  }

  /// hash value of 3 int
  template <typename TI>    
  INLINE size_t HashValue (const INT<3,TI> & ind, size_t size) 
  {
    INT<3,size_t> lind = ind;
    return (113*lind[0]+59*lind[1]+lind[2]) % size;
  }

  INLINE size_t HashValue (size_t ind, size_t size)
  {
    return ind%size;
  }
  INLINE size_t HashValue (int ind, size_t size)
  {
    return size_t(ind)%size;
  }
  





  
  template <int N, typename TI>
  INLINE size_t HashValue2 (const INT<N,TI> & ind, size_t mask)
  {
    INT<N,size_t> lind = ind;    
    size_t sum = 0;
    for (int i = 0; i < N; i++)
      sum += lind[i];
    return sum & mask;
  }

  /// hash value of 1 int
  template <typename TI>
  INLINE size_t HashValue2 (const INT<1,TI> & ind, size_t mask) 
  {
    return ind[0] & mask;
  }

  /// hash value of 2 int
  template <typename TI>  
  INLINE size_t HashValue2 (const INT<2,TI> & ind, size_t mask) 
  {
    INT<2,size_t> lind = ind;
    return (113*lind[0]+lind[1]) & mask;
  }

  /// hash value of 3 int
  template <typename TI>    
  INLINE size_t HashValue2 (const INT<3,TI> & ind, size_t mask) 
  {
    INT<3,size_t> lind = ind;
    return (113*lind[0]+59*lind[1]+lind[2]) & mask;
  }

  INLINE size_t HashValue2 (size_t ind, size_t mask)
  {
    return ind & mask;
  }
  INLINE size_t HashValue2 (int ind, size_t mask)
  {
    return size_t(ind) & mask;
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

  template <int D, typename T>
  INLINE T Min (const INT<D,T> & i)
  {
    if (D == 0) return 0;
    T m = i[0];
    for (int j = 1; j < D; j++)
      if (i[j] < m) m = i[j];
    return m;
  }

  template <int D, typename T>
  INLINE INT<D,T> Max (INT<D,T> i1, INT<D,T> i2)
  {
    INT<D,T> tmp;
    for (int i = 0; i < D; i++)
      tmp[i] = max2(i1[i], i2[i]);
    return tmp;
  }

  template <int D, typename T>
  INLINE INT<D,T> operator+ (INT<D,T> i1, INT<D,T> i2)
  {
    INT<D,T> tmp;
    for (int i = 0; i < D; i++)
      tmp[i] = i1[i]+i2[i];
    return tmp;
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
    /*
    DynamicTable<T_HASH> hash;
    DynamicTable<T> cont;
    */
    DynamicTable<pair<T_HASH,T>> table;
  public:
    /// Constructs a hashtable of size bags.
    INLINE HashTable (int size)
    // : hash(size), cont(size)
      : table(size)
    { ; }
    INLINE ~HashTable () { ; }

    /// Sets identifier ahash to value acont
    void Set (const T_HASH & ahash, const T & acont)
    {
      int bnr = HashValue (ahash, Size());
      int pos = CheckPosition (bnr, ahash);
      if (pos != -1)
	// cont.Set (bnr, pos, acont);
        table[bnr][pos].second = acont;
      else
	{
	  // hash.Add (bnr, ahash);
	  // cont.Add (bnr, acont);
          table.Add (bnr, make_pair(ahash, acont));
	}        
    }

    /// get value of identifier ahash, exception if unused
    const T & Get (const T_HASH & ahash) const
    {
      int bnr = HashValue (ahash, Size());
      int pos = Position (bnr, ahash);
      // return cont.Get (bnr, pos);
      return table.Get (bnr, pos).second;
    }

    /// get value of identifier ahash, exception if unused
    const T & Get (int bnr, int pos) const
    {
      // return cont.Get (bnr, pos);
      return table.Get (bnr, pos).second;
    }

    /// is identifier used ?
    bool Used (const T_HASH & ahash) const
    {
      // return (CheckPosition (HashValue (ahash, hash.Size()), ahash) != -1);
      return (CheckPosition (HashValue (ahash, table.Size()), ahash) != -1);
    }

    /// is identifier used ?
    bool Used (const T_HASH & ahash, int & bnr, int & pos) const
    {
      // bnr = HashValue (ahash, hash.Size());
      bnr = HashValue (ahash, Size());
      pos = CheckPosition (bnr, ahash);
      return (pos != -1);
    }


    /// number of hash entries
    size_t Size () const
    {
      // return hash.Size();
      return table.Size();
    }

    /// size of hash entry
    size_t EntrySize (int bnr) const
    {
      // return hash[bnr].Size();
      return table[bnr].Size();
    }

    /// get identifier and value of entry bnr, position colnr
    void GetData (int bnr, int colnr, T_HASH & ahash, T & acont) const
    {
      // ahash = hash[bnr][colnr];
      // acont = cont[bnr][colnr];
      ahash = table[bnr][colnr].first;
      acont = table[bnr][colnr].second;
    }

    /// set identifier and value of entry bnr, position colnr
    void SetData (int bnr, int colnr, const T_HASH & ahash, const T & acont)
    {
      // hash[bnr][colnr] = ahash;
      // cont[bnr][colnr] = acont;
      table[bnr][colnr] = make_pair(ahash, acont);
    }    

    /// returns position of index. returns -1 on unused
    int CheckPosition (int bnr, const T_HASH & ind) const
    {
      /*
      for (int i = 0; i < hash[bnr].Size(); i++)
	if (hash[bnr][i] == ind)
	  return i;
      */
      for (int i = 0; i < table[bnr].Size(); i++)
	if (table[bnr][i].first == ind)
	  return i;
      return -1;
    }

    /// returns position of index. exception on unused
    int Position (int bnr, const T_HASH & ind) const
    {
      for (int i = 0; i < table[bnr].Size(); i++)
	if (table[bnr][i].first == ind)
	  return i;
      throw Exception ("Ask for unsused hash-value");
    }

    T & operator[] (T_HASH ahash)
    {
      int bnr, pos;
      if (Used (ahash, bnr, pos))
        return table[bnr][pos].second;
      else
        {
	  // hash.Add (bnr, ahash);
	  // cont.Add (bnr, T(0));
          table.Add (bnr, make_pair(ahash, T(0)));
          // return cont[bnr][cont[bnr].Size()-1];
          return table[bnr][table[bnr].Size()-1].second;
        }
    }

    const T & operator[] (T_HASH ahash) const
    {
      return Get(ahash);
    }

    class Iterator
    {
      const HashTable & ht;
      int bnr;
      int pos;
    public:
      Iterator (const HashTable & aht, int abnr, int apos)
        : ht(aht), bnr(abnr), pos(apos) { ; }
      pair<T_HASH,T> operator* () const
      {
        T_HASH hash; 
        T data;
        ht.GetData (bnr, pos, hash, data);
        return pair<T_HASH,T> (hash, data);
      }

      Iterator & operator++() 
      {
        pos++;
        if (pos == ht.EntrySize(bnr))
          {
            pos = 0;
            bnr++;
            for ( ; bnr < ht.Size(); bnr++)
              if (ht.EntrySize(bnr) != 0) break;
          }
        return *this;
      }
      
      bool operator!= (const Iterator & it2) { return bnr != it2.bnr || pos != it2.pos; }
    };

    Iterator begin () const 
    {
      int i = 0;
      for ( ; i < Size(); i++)
        if (EntrySize(i) != 0) break;
      return Iterator(*this, i,0); 
    }
    Iterator end () const { return Iterator(*this, Size(),0); }
  };



  inline size_t RoundUp2 (size_t i)
  {
    size_t res = 1;
    while (res < i) res *= 2; // hope it will never be too large 
    return res; 
  }



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
    size_t size;
    size_t mask;
    ///
    size_t used;
    ///
    Array<T_HASH> hash;
    ///
    Array<T> cont;
    ///
    T_HASH invalid;
  public:
    ///
    ClosedHashTable (size_t asize = 128)
      : size(RoundUp2(asize)), used(0), hash(size), cont(size)
    {
      mask = size-1;
      invalid = -1; 
      hash = T_HASH(invalid);
    }

    ClosedHashTable (ClosedHashTable && ht2) = default;

      // who needs that ? 
    ClosedHashTable (FlatArray<T_HASH> _hash, FlatArray<T> _cont)
      : size(_hash.Size()), used(0), hash(_hash.Size(), _hash.Addr(0)), cont(_cont.Size(), _cont.Addr(0))
    {
      invalid = -1; 
      hash = T_HASH(invalid);
    }

    /// allocate on local heap
    ClosedHashTable (size_t asize, LocalHeap & lh)
      : size(asize), hash(asize, lh), cont(asize, lh)
    {
      invalid = -1; 
      hash = T_HASH(invalid);
    }

    ClosedHashTable & operator= (ClosedHashTable && ht2) = default;

    /// 
    size_t Size() const
    {
      return size;
    }

    /// is position used
    bool UsedPos (size_t pos) const
    {
      return ! (hash[pos] == invalid); 
    }

    /// number of used elements
    size_t UsedElements () const
    {
      return used;
      /*
      size_t cnt = 0;
      for (size_t i = 0; i < size; i++)
	if (hash[i] != invalid)
	  cnt++;
      return cnt;
      */
    }

    size_t Position (const T_HASH ind) const
    {
      size_t i = HashValue2(ind, mask);
      while (1)
	{
	  if (hash[i] == ind) return i;
	  if (hash[i] == invalid) return size_t(-1);
	  i++;
	  if (i >= size) i = 0;
	}
    }

    void DoubleSize()
    {
      ClosedHashTable tmp(2*Size());
      for (auto both : *this)
        tmp[both.first] = both.second;
      *this = move(tmp);
    }
    
    // returns true if new position is created
    bool PositionCreate (const T_HASH ind, size_t & apos)
    {
      if (UsedElements()*2 > Size()) DoubleSize();
      
      size_t i = HashValue2 (ind, mask);

      while (1)
	{
	  if (hash[i] == invalid)
	    { 
	      hash[i] = ind; 
	      apos = i;
              used++;
	      return true;
	    }
	  if (hash[i] == ind) 
	    { 
	      apos = i; 
	      return false; 
	    }
	  i++;
	  if (i >= size) i = 0;
	}
    }


    ///
    void Set (const T_HASH & ahash, const T & acont)
    {
      size_t pos;
      PositionCreate (ahash, pos);
      hash[pos] = ahash;
      cont[pos] = acont;
    }

    ///
    const T & Get (const T_HASH & ahash) const
    {
      size_t pos = Position (ahash);
      if (pos == size_t(-1))
        throw Exception (string("illegal key: ") + ToString(ahash) );
      return cont[pos];
    }

    ///
    bool Used (const T_HASH & ahash) const
    {
      return (Position (ahash) != size_t(-1));
    }

    void SetData (size_t pos, const T_HASH & ahash, const T & acont)
    {
      hash[pos] = ahash;
      cont[pos] = acont;
    }

    void GetData (size_t pos, T_HASH & ahash, T & acont) const
    {
      ahash = hash[pos];
      acont = cont[pos];
    }
  
    void SetData (size_t pos, const T & acont)
    {
      cont[pos] = acont;
    }

    void GetData (size_t pos, T & acont) const
    {
      acont = cont[pos];
    }

    pair<T_HASH,T> GetBoth (size_t pos) const
    {
      return pair<T_HASH,T> (hash[pos], cont[pos]);
    }

    const T & operator[] (T_HASH key) const { return Get(key); }
    T & operator[] (T_HASH key)
    {
      size_t pos;
      PositionCreate(key, pos);
      return cont[pos];
    }
    
    void SetSize (size_t asize)
    {
      size = asize;
      hash.Alloc(size);
      cont.Alloc(size);

      // for (size_t i = 0; i < size; i++)
      // hash[i] = invalid;
      hash = T_HASH(invalid);
    }

    void Delete (T_HASH key)
    {
      size_t pos = Position(key);
      if (pos == size_t(-1)) return;
      hash[pos] = invalid; used--;
      
      while (1)
        {
          size_t nextpos = pos+1;
          if (nextpos == size) nextpos = 0;
          if (hash[nextpos] == invalid) break;
          
          auto key = hash[nextpos];
          auto val = cont[nextpos];
          hash[pos] = invalid; used--;
          
          Set (key, val);
          pos = nextpos;
        }
    }
    
    class Iterator
    {
      const ClosedHashTable & tab;
      size_t nr;
    public:
      Iterator (const ClosedHashTable & _tab, size_t _nr)
        : tab(_tab), nr(_nr)
      {
        while (nr < tab.Size() && !tab.UsedPos(nr)) nr++;
      }
      Iterator & operator++()
      {
        nr++;
        while (nr < tab.Size() && !tab.UsedPos(nr)) nr++;
        return *this;
      }
      bool operator!= (const Iterator & it2) { return nr != it2.nr; }
      auto operator* () const
      {
        T_HASH hash;
        T val;
        tab.GetData(nr, hash,val);
        return std::make_pair(hash,val);
      }
    };

    Iterator begin() const { return Iterator(*this, 0); }
    Iterator end() const { return Iterator(*this, Size()); } 
  };

  template <class T_HASH, class T>  
  ostream & operator<< (ostream & ost,
                        const ClosedHashTable<T_HASH,T> & tab)
  {
    for (size_t i = 0; i < tab.Size(); i++)
      if (tab.UsedPos(i))
        {
          T_HASH key;
          T val;
          tab.GetData (i, key, val);
          ost << key << ": " << val << ", ";
        }
    return ost;
  }

  
  template <typename TI>  
  INLINE size_t HashValue (const INT<2,TI> ind)
  {
    INT<2,size_t> lind = ind;
    return 113*lind[0]+lind[1];
  }


  template <typename TKEY, typename T>
  class ParallelHashTable
  {
    class ClosedHT
    {
      Array<TKEY> keys;
      Array<T> values;
      size_t used;
      
    public:
      ClosedHT(size_t asize = 256) : keys(asize), values(asize), used(0)
      {
        keys = TKEY(-1);
      }
      
      size_t Used () const { return used; }

      ClosedHT & operator= (ClosedHT&&) = default;

      void Resize()
      {
        ClosedHT tmp(keys.Size()*2);
        for (size_t i = 0; i < keys.Size(); i++)
          if (keys[i] != TKEY(-1))
            {
              TKEY hkey = keys[i];
              T hval = values[i];
              size_t hhash = HashValue(hkey);
              size_t hhash2 = hhash / 256;
              tmp.DoSave(hkey, [hval] (T & v) { v = hval; }, hhash2);
            }
        (*this) = move(tmp);
      }
      
      template <typename TFUNC>
      auto Do (TKEY key, TFUNC func, size_t hash)
      {
        if (used > keys.Size()/2)
          Resize();
        return DoSave (key, func, hash);
      }
      
      template <typename TFUNC>
      auto DoSave (TKEY key, TFUNC func, size_t hash)
      {
        size_t pos = hash & (keys.Size()-1);
        while (1)
          {
            if (keys[pos] == key)
              break;
            if (keys[pos] == TKEY(-1))
              {
                keys[pos] = key;
                values[pos] = T(0);
                used++;
                break;
              }
            pos++;
            if (pos == keys.Size()) pos = 0;
          }
        return func(values[pos]);
      }
      
      T Get (TKEY key, size_t hash)
      {
        size_t pos = hash & (keys.Size()-1);
        while (1)
          {
            if (keys[pos] == key)
              return values[pos];
            if (keys[pos] == TKEY(-1))
              throw Exception ("ParallelHashTable::Get of unused key");
            pos++;
            if (pos == keys.Size()) pos = 0;
          }
      }

      size_t GetCosts (TKEY key, size_t hash)
      {
        size_t pos = hash & (keys.Size()-1);
        size_t costs = 1;
        while (1)
          {
            if (keys[pos] == key)
              return costs;
            if (keys[pos] == TKEY(-1))
              throw Exception ("ParallelHashTable::Get of unused key");
            costs++;
            pos++;
            if (pos == keys.Size()) pos = 0;
          }
      }


      template <typename TFUNC>
      void Iterate (TFUNC func) const
      {
        for (size_t i = 0; i < keys.Size(); i++)
          if (keys[i] != TKEY(-1))
            func(keys[i], values[i]);
      }
        
      void Print (ostream & ost) const
      {
        for (size_t i = 0; i < keys.Size(); i++)
          if (keys[i] != TKEY(-1))
            ost << keys[i] << ": " << values[i] << ", ";
      }
    };

    Array<ClosedHT> hts;
    class alignas(64) MyMutex64 : public MyMutex { };
    
    Array<MyMutex64> locks;

  public:
    ParallelHashTable() : hts(256), locks(256) { ; }
    size_t NumBuckets() const { return hts.Size(); }
    size_t Used() const
    {
      size_t used = 0;
      for (auto & ht : hts)
        used += ht.Used();
      return used;
    }  
    template <typename TFUNC>
    auto Do (TKEY key, TFUNC func)
    {
      size_t hash = HashValue(key);
      size_t hash1 = hash % 256;
      size_t hash2 = hash / 256;
      
      // locks[hash1].lock();
      // hts[hash1].Do (key, func, hash2);
      // locks[hash1].unlock();
      MyLock lock(locks[hash1]);
      return hts[hash1].Do (key, func, hash2);
    }
    
    T Get (TKEY key)
    {
      size_t hash = HashValue(key);
      size_t hash1 = hash % 256;
      size_t hash2 = hash / 256;
      
      return hts[hash1].Get (key, hash2);
    }

    auto GetCosts (TKEY key)
    {
      size_t hash = HashValue(key);
      size_t hash1 = hash % 256;
      size_t hash2 = hash / 256;
      
      return hts[hash1].GetCosts (key, hash2);
    }

    
    template <typename TFUNC>
    void Iterate(TFUNC func) const
    {
      for (auto & bucket : hts)
        bucket.Iterate(func);
    }

    template <typename TFUNC>
    void Iterate(size_t nr, TFUNC func) const
    {
      hts[nr].Iterate(func);
    }
    

    void Print (ostream & ost) const
    {
      for (size_t i : Range(hts))
        if (hts[i].Used() > 0)
          {
            ost << i << ": ";
            hts[i].Print(ost);
          }
    }
  };

  template <typename TKEY, typename T>
  inline ostream & operator<< (ostream & ost, const ParallelHashTable<TKEY,T> & ht)
  {
    ht.Print(ost);
    return ost;
  }


    
  template <int N, typename T>
  Archive & operator & (Archive & archive, INT<N,T> & mi)
  {
    for (int i = 0; i < N; i++)
      archive & mi[i];
    return archive;
  }
}

#endif
