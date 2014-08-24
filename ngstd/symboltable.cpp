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
  BaseSymbolTable :: BaseSymbolTable ()
  {
    ;
  }


  BaseSymbolTable :: ~BaseSymbolTable()
  {
    DelNames();
  }


  void BaseSymbolTable :: DelNames()
  {
    /*
    for (int i = 0; i < names.Size(); i++)
      delete [] names[i];
    */
    names.SetSize (0);
  }

  int BaseSymbolTable :: Index (const string & name) const
  {
    /*
    if (name)
      for (int i = 0; i < names.Size(); i++)
	if (strcmp (names[i], name) == 0) return i;
    */
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
    /*
    char * hname = new char [strlen (name) + 1];
    strcpy (hname, name);
    names.Append (hname);
    */
    names.Append (name);
  }

}



#endif
