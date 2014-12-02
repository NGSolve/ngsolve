/**************************************************************************/
/* File:   symboltable.cpp                                                */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Jun. 95                                                    */
/**************************************************************************/

/* 
   Abstract data type Symbol Table
*/
#ifndef FILE_SYMBOLTABLECC
#define FILE_SYMBOLTABLECC
// necessary for SGI ????




#include <ngstd.hpp>

namespace ngstd
{
  // BaseSymbolTable :: BaseSymbolTable () { ; }

  BaseSymbolTable :: BaseSymbolTable (const BaseSymbolTable & tab2)
    : names (tab2.names)
  { ; }


  BaseSymbolTable :: ~BaseSymbolTable() { ; }


  void BaseSymbolTable :: DelNames()
  {
    names.SetSize(0);
  }


  int BaseSymbolTable :: Index (const string & name) const
  {
    for (int i = 0; i < names.Size(); i++)
      if (names[i] == name) return i;

    stringstream str;
    str << "SymbolTable: unused name '" << name << "'" << endl;
    throw Exception (str.str());

  }


  int BaseSymbolTable :: CheckIndex (const string & name) const
  {
    for (int i = 0; i < names.Size(); i++)
      if (names[i] == name) return i;
    return -1;
  }


  void BaseSymbolTable :: AppendName (const string & name)
  {
    names.Append (name);
  }

}



#endif
