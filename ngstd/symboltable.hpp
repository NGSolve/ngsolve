#ifndef FILE_NGS_SYMBOLTABLE
#define FILE_NGS_SYMBOLTABLE

/**************************************************************************/
/* File:   symboltable.hpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

namespace ngstd
{

/**
  Base class for generic SymbolTable.
  Maintains the array of identifiers.
 */
class BaseSymbolTable
{
protected:
  /// identifiers
  Array <string> names;
  
public:
  /// 
  NGS_DLL_HEADER BaseSymbolTable () = default;
  ///
  NGS_DLL_HEADER BaseSymbolTable (const BaseSymbolTable & tab2);
  NGS_DLL_HEADER BaseSymbolTable (BaseSymbolTable && tab2) = default;

  /// deletes identifiers
  NGS_DLL_HEADER ~BaseSymbolTable ();
  /// delete all symbols
  NGS_DLL_HEADER void DelNames ();

  /// append new name (copy)
  NGS_DLL_HEADER void AppendName (const string & name);

  /// Index of symbol name, throws exception if unused
  NGS_DLL_HEADER int Index (const string & name) const;

  /// Index of symbol name, returns -1 if unused
  NGS_DLL_HEADER int CheckIndex (const string & name) const;
};






/** 
    A symbol table.
   
    The symboltable provides a mapping from string identifiers
    to the generic type T. The strings are copied.
*/
template <class T>
class SymbolTable : public BaseSymbolTable
{
  /// the data
  Array <T> data;
public:
  /// Creates a symboltable
  SymbolTable () = default;
  SymbolTable (const SymbolTable & tab2) = default;
  SymbolTable (SymbolTable && tab2) = default;

  /// number of identifiers
  int Size() const
  { 
    return data.Size(); 
  }

  /*
  /// Returns reference to element. exception for unused identifier
  T & operator[] (const char * name)
  {
    return data[Index (name)]; 
  }
  */

  /// Returns reference to element. exception for unused identifier
  T & operator[] (const string & name)
  {
    return data[Index (name)]; 
  }

  /*
  /// Returns element, error if not used
  const T & operator[] (const char * name) const
  {
    return data[Index (name)]; 
  }
  */
  const T & operator[] (const string & name) const
  {
    return data[Index (name)];
  }

  /// Returns reference to i-th element
  T & operator[] (int i)
  { return data[i]; } 

  /// Returns const reference to i-th element
  const T & operator[] (int i) const
  { return data[i]; } 

  /// Returns name of i-th element
  const string & GetName (int i) const
  { return names[i]; }


  /// Associates el to the string name, overrides if name is used
  void Set (const string & name, const T & el)
  {
    int i = CheckIndex (name);
    if (i >= 0) 
      data[i] = el;
    else
      {
	data.Append (el);
	AppendName (name);
      }
  }

  bool Used (const string & name) const
  {
    return (CheckIndex(name) >= 0) ? 1 : 0;
  }

  /// Deletes symboltable
  inline void DeleteAll ()
  {
    DelNames();
    data.DeleteAll();
  }

  SymbolTable<T> & operator= (const SymbolTable<T> & tab2)
  {
    for (int i = 0; i < tab2.Size(); i++)
      Set (tab2.GetName(i), tab2[i]);
    return *this;
  }
};

  template <typename T>
  ostream & operator<< (ostream & ost, const SymbolTable<T> & st)
  {
    for (int i = 0; i < st.Size(); i++)
      ost << st.GetName(i) << " : " << st[i] << endl;
    return ost;
  }

  template <typename T>
  INLINE void SymTabEl_DummyInit (T & data) { ; }
  INLINE void SymTabEl_DummyInit (string *& data) { data = nullptr; }

  template <typename T> 
  Archive & operator & (Archive & ar, SymbolTable<T> & table)
  {
    if (ar.Output())
      {
        ar << table.Size();
        for (int i = 0; i < table.Size(); i++)
          {
            ar << string (table.GetName(i));
            ar & table[i];
          }
      }
    else
      {
        int s;
        ar & s;
        for (int i = 0; i < s; i++)
          {
            string name;
            T entry;
            SymTabEl_DummyInit (entry);  // no un-initialized warning
            ar & name & entry;
            table.Set (name, entry);
          }
      }
    return ar;
  }


}

#endif
