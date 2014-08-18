#ifndef FILE_NGS_TABLE
#define FILE_NGS_TABLE

/**************************************************************************/
/* File:   table.hpp                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 2000                                                  */
/**************************************************************************/


namespace ngstd
{



#ifdef OLD

/**
   Base class of \Ref{Table} container.
   Provides index array pointing to first element in entry. 
 */
class BaseTable
{
protected:
  /// number of rows
  int size;
  /// pointer to first in row
  int * index;

  INLINE BaseTable ()
    : size(0), index(0) { ; }

public:  
  INLINE BaseTable (int asize, int entrysize)
  {
    size = asize;
    index = new int[size+1];
    for (int i = 0; i <= size; i++)
      index[i] = i*entrysize;
  }

  INLINE BaseTable (const FlatArray<int> & entrysize)
  {
    int i, cnt = 0;
    size  = entrysize.Size();
    
    index = new int[size+1];
    for (i = 0; i < size; i++)
      {
	index[i] = cnt;
	cnt += entrysize[i];
      }
    index[i] = cnt;
  }

  INLINE ~BaseTable ()
  {
    delete [] index;
  }
  INLINE int * Index() const { return index; }
};
#endif



template <class T>
class FlatTable 
{
protected:
  /// number of rows
  int size;
  /// pointer to first in row
  int * index;
  /// array of data 
  T * data;

public:
  INLINE FlatTable() { ; }

  INLINE FlatTable(int as, int * aindex, T * adata)
    : size(as), index(aindex), data(adata) { ; }

  /// Size of table
  INLINE int Size() const { return size; }

  /// Access entry
  INLINE const FlatArray<T> operator[] (int i) const 
  { 
    return FlatArray<T> (index[i+1]-index[i], data+index[i]); 
  }

  INLINE T * Data() const { return data; }

  INLINE FlatArray<T> AsArray() const
  {
    return FlatArray<T> (index[size], data);
  }

  INLINE FlatArray<int> IndexArray() const
  {
    return FlatArray<int> (size+1, index);
  }
};



/** 
    A compact Table container.
    A table contains size entries of variable size. 
    The entry sizes must be known at construction.
*/
template <class T>
class Table : public FlatTable<T>
{
protected:

  using FlatTable<T>::size;
  using FlatTable<T>::index;
  using FlatTable<T>::data;

  /*
  INLINE Table ()
    : data(NULL) { ; }
  */

public:
  /// Construct table of uniform entrysize
  INLINE Table (int asize, int entrysize)
  // : BaseTable (asize, entrysize) 
  { 
    size = asize;
    index = new int[size+1];
    for (int i = 0; i <= size; i++)
      index[i] = i*entrysize;
    data = new T [size*entrysize]; 
  }

  /// Construct table of variable entrysize
  INLINE Table (const FlatArray<int> & entrysize)
  // : BaseTable (entrysize)
  {
    int i, cnt = 0;
    size  = entrysize.Size();
    
    index = new int[size+1];
    for (i = 0; i < size; i++)
      {
	index[i] = cnt;
	cnt += entrysize[i];
      }
    index[i] = cnt;

    data = new T[cnt];
  }

  INLINE Table (const Table<int> & tab2) = delete;
  INLINE Table (Table<int> && tab2)
  {
    size = 0;
    index = NULL;
    data = NULL;

    Swap (size, tab2.size);
    Swap (index, tab2.index);
    Swap (data, tab2.data);
  }


  /// Delete data
  INLINE ~Table ()
  {
    delete [] data; 
    delete [] index;
  }

  /// Size of table
  // INLINE int Size() const { return size; }
  using FlatTable<T>::Size;
  
  /// number of elements in all rows
  INLINE int NElements() const { return index[size]; }

  /*
  /// Access entry.
  FlatArray<T> operator[] (int i) 
  { 
    return FlatArray<T> (index[i+1]-index[i], data+index[i]); 
  }
  */
  /*
  /// Access entry
  INLINE const FlatArray<T> operator[] (int i) const 
  { 
    return FlatArray<T> (index[i+1]-index[i], data+index[i]); 
  }
  */

  using FlatTable<T>::operator[];
};


/// Print table
template <class T>
inline ostream & operator<< (ostream & s, const Table<T> & table)
{
  for (int i = 0; i < table.Size(); i++)
    {
      s << i << ":";
      for (int j = 0; j < table[i].Size(); j++)
	s << " " << table[i][j];
      s << "\n";
    }
  s << flush;
  return s;
}




template <class T>
  class TableCreator
  {
  protected:  
    int mode;    // 1 .. cnt, 2 .. cnt entries, 3 .. fill table
    int nd;
    Array<int> cnt;
    Table<T> * table;
  public:
    TableCreator()
    { nd = 0; mode = 1; table = NULL; }
    TableCreator(int acnt)
    { nd = acnt; table = NULL; SetMode(2); }
    
    Table<T> * GetTable() { return table; }

    operator Table<T> && () 
    {
      Table<int> tmp (std::move(*table));
      delete table;
      table = NULL;
      return tmp;
    }

    bool Done () { return mode > 3; }
    void operator++(int) { SetMode (mode+1); }

    int GetMode () const { return mode; }
    void SetMode (int amode) 
    {
      mode = amode; 
      if (mode == 2) 
	{
	  cnt.SetSize(nd);
	  cnt = 0; 
	}
      if (mode == 3)
	{
	  table = new Table<T> (cnt);
	  int sum = 0;
	  for (int i = 0; i < cnt.Size(); i++)
	    {
	      int nsum = sum + cnt[i];
	      cnt[i] = sum;
	      sum = nsum;
	    }
	}
    }

    void Add (int blocknr, const T & data)
    {
      switch (mode)
	{
	case 1:
	  if (blocknr+1 > nd) nd = blocknr+1; 
	  break;
	case 2:
	  cnt[blocknr]++;
	  break;
	case 3:
	  // (*table)[blocknr][cnt[blocknr]++] = data;
	  table -> Data() [cnt[blocknr]++] = data;
	  break;
	}
    }


    void Add (int blocknr, IntRange range)
    {
      switch (mode)
	{
	case 1:
	  if (blocknr+1 > nd) nd = blocknr+1; 
	  break;
	case 2:
	  cnt[blocknr]+=range.Size();
	  break;
	case 3:
	  /*
	  for (int j = 0; j < range.Size(); j++)
	    (*table)[blocknr][cnt[blocknr]+j] = range.First()+j;
	  cnt[blocknr]+=range.Size();
	  */
	  for (int j = 0; j < range.Size(); j++)
	    table->Data()[cnt[blocknr]+j] = range.First()+j;
	  cnt[blocknr]+=range.Size();

	  break;
	}
    }

    void Add (int blocknr, const FlatArray<int> & dofs)
    {
      switch (mode)
	{
	case 1:
	  if (blocknr+1 > nd) nd = blocknr+1; 
	  break;
	case 2:
	  cnt[blocknr]+=dofs.Size();
	  break;
	case 3:
	  /*
	  for (int j = 0; j < dofs.Size(); j++)
	    (*table)[blocknr][cnt[blocknr]+j] = dofs[j];
	  cnt[blocknr]+=dofs.Size();
	  */
	  for (int j = 0; j < dofs.Size(); j++)
	    table->Data()[cnt[blocknr]+j] = dofs[j];
	  cnt[blocknr]+=dofs.Size();

	  break;
	}
    }
  };

  class BitArray;
  
  class NGS_DLL_HEADER FilteredTableCreator : public TableCreator<int>
  {
  protected:
    const BitArray* takedofs;  
  public:
    FilteredTableCreator(const BitArray* atakedofs) 
      : TableCreator<int>(), takedofs(atakedofs) { };
    FilteredTableCreator(int acnt, const BitArray* atakedofs)
      : TableCreator<int>(acnt),takedofs(atakedofs) { };
    void Add (int blocknr, int data);
    void Add (int blocknr, IntRange range);
    void Add (int blocknr, FlatArray<int> dofs);
  };








/// Base class to generic DynamicTable.
class BaseDynamicTable
{
protected:
  
  ///
  struct linestruct
  {
    ///
    int size;
    ///
    int maxsize;
    ///
    void * col;
  };
  
  ///
  Array<linestruct> data;
  ///
  char * oneblock;

public:
  ///
  NGS_DLL_HEADER BaseDynamicTable (int size);
  ///
  NGS_DLL_HEADER BaseDynamicTable (const Array<int> & entrysizes, int elemsize);
  ///
  NGS_DLL_HEADER ~BaseDynamicTable ();

  /// Changes Size of table to size, deletes data
  NGS_DLL_HEADER void SetSize (int size);
  ///
  NGS_DLL_HEADER void IncSize (int i, int elsize);

  NGS_DLL_HEADER void DecSize (int i);
};



/** 
    A dynamic table class.
   
    A DynamicTable contains entries of variable size. Entry sizes can
    be increased dynamically.
*/
template <class T>
class DynamicTable : public BaseDynamicTable
{
public:
  /// Creates table of size size
  DynamicTable (int size = 0)
    : BaseDynamicTable (size) { ; }

  /// Creates table with a priori fixed entry sizes.
  DynamicTable (const Array<int> & entrysizes)
    : BaseDynamicTable (entrysizes, sizeof(T)) { ; }

  /// Inserts element acont into row i. Does not test if already used.
  void Add (int i, const T & acont)
  {
    if (data[i].size == data[i].maxsize)
      IncSize (i, sizeof (T));
    else
      data[i].size++;
    static_cast<T*> (data[i].col) [data[i].size-1] = acont;
  }

  /// Inserts element acont into row i, iff not yet exists.
  void AddUnique (int i, const T & cont)
  {
    int es = EntrySize (i);
    int * line = const_cast<int*> (GetLine (i));
    for (int j = 0; j < es; j++)
      if (line[j] == cont)
	return;
    Add (i, cont);
  }


  /// Inserts element acont into row i. Does not test if already used.
  void AddEmpty (int i)
  {
    IncSize (i, sizeof (T));
  }

  /** Set the nr-th element in the i-th row to acont.
      Does not check for overflow. */
  void Set (int i, int nr, const T & acont)
  { static_cast<T*> (data[i].col)[nr] = acont; }
  

  /** Returns the nr-th element in the i-th row.
    Does not check for overflow. */
  const T & Get (int i, int nr) const
  { return static_cast<T*> (data[i].col)[nr]; }


  /** Returns pointer to the first element in row i. */
  const T * GetLine (int i) const
  { return static_cast<T*> (data[i].col); }
  
  
  /// Returns size of the table.
  int Size () const
  { return data.Size(); }

  /// Returns size of the i-th row.
  int EntrySize (int i) const
  { return data[i].size; }
  
  ///
  void DecEntrySize (int i)
  { DecSize(i); }

  /// Access entry i
  FlatArray<T> operator[] (int i) 
  { return FlatArray<T> (data[i].size, static_cast<T*> (data[i].col)); }

  typedef const FlatArray<T> ConstFlatArray;
  /// Access entry i
  ConstFlatArray operator[] (int i) const
  { return FlatArray<T> (data[i].size, static_cast<T*> (data[i].col)); }
};




/// Print table
template <class T>
inline ostream & operator<< (ostream & s, const DynamicTable<T> & table)
{
  for (int i = 0; i < table.Size(); i++)
    {
      s << i << ":";
      for (int j = 0; j < table[i].Size(); j++)
	s << " " << table[i][j];
      s << "\n";
    }
  s << flush;
  return s;
}

typedef DynamicTable<int> IntTable;

}

#endif

