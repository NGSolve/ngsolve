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
  using namespace ngstd;


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
    for (int i = 0; i < names.Size(); i++)
      delete [] names[i];
    names.SetSize (0);
  }

  int BaseSymbolTable :: Index (const char * name) const
  {
    if (name)
      for (int i = 0; i < names.Size(); i++)
	if (strcmp (names[i], name) == 0) return i;
    
    stringstream str;
    if (name)
      str << "SymbolTable: unused name '" << name << "'" << endl;
    else
      str << "SymbolTable: name = 0" << endl;
      
    throw Exception (str.str());
  }

  int BaseSymbolTable :: CheckIndex (const char * name) const
  {
    if (name)
      for (int i = 0; i < names.Size(); i++)
	if (strcmp (names[i], name) == 0) return i;
    return -1;
  }


  void BaseSymbolTable :: AppendName (const char * name)
  {
    char * hname = new char [strlen (name) + 1];
    strcpy (hname, name);
    names.Append (hname);
  }

}



#endif
