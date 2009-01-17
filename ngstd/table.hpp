#ifndef FILE_NGS_TABLE
#define FILE_NGS_TABLE

/**************************************************************************/
/* File:   table.hh                                                       */
/* Author: Joachim Schoeberl                                              */
/* Date:   25. Mar. 2000                                                  */
/**************************************************************************/







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

public:  
  BaseTable (int asize, int entrysize);
  BaseTable (const FlatArray<int> & entrysize);
  ~BaseTable ();
  int * Index() const { return index; }
};


/** 
    A compact Table container.
    A table contains size entries of variable size. 
    The entry sizes must be known at construction.
*/
template <class T>
class Table : public BaseTable
{
protected:
  /// array of data 
  T * data;
  // MoveableMem<T> data;

public:
  /// Construct table of uniform entrysize
  Table (int asize, int entrysize)
    : BaseTable (asize, entrysize) 
  { 
    data = new T[index[size]]; 
    // data.Alloc (index[size]);
    // data.SetName ("NGS Table");
  }

  /// Construct table of variable entrysize
  Table (const FlatArray<int> & entrysize)
    : BaseTable (entrysize)
  {
    data = new T[index[size]]; 
    // data.Alloc (index[size]);
    // data.SetName ("NGS Table");
  }

  /// Delete data
  ~Table ()
  {
    delete [] data; 
    // data.Free();
  }

  /// Size of table
  int Size() const { return size; }
  
  /// number of elements in all rows
  int NElements() const { return index[size]; }

  /// Access entry.
  FlatArray<T> operator[] (int i) 
  { 
    return FlatArray<T> (index[i+1]-index[i], data+index[i]); 
  }

  /// Access entry
  FlatArray<T> operator[] (int i) const 
  { 
    return FlatArray<T> (index[i+1]-index[i], data+index[i]); 
  }

  T * Data() const { return data; }
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
  ARRAY<linestruct> data;
  ///
  char * oneblock;

public:
  ///
  BaseDynamicTable (int size);
  ///
  BaseDynamicTable (const ARRAY<int> & entrysizes, int elemsize);
  ///
  ~BaseDynamicTable ();

  /// Changes Size of table to size, deletes data
  void SetSize (int size);
  ///
  void IncSize (int i, int elsize);

  void DecSize (int i);
};



/** 
    A dynamic table class.
   
    A DynamicTable contains entries of variable size. Entry size can
    be dynamically increased.
*/
template <class T>
class DynamicTable : public BaseDynamicTable
{
public:
  /// Creates table of size size
  DynamicTable (int size = 0)
    : BaseDynamicTable (size) { ; }

  /// Creates table with a priori fixed entry sizes.
  DynamicTable (const ARRAY<int> & entrysizes)
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

#endif

