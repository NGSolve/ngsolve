/**************************************************************************/
/* File:   archive.cpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   14. Feb. 2014                                                  */
/**************************************************************************/

#include <ngstd.hpp>

#include "archive.hpp"


namespace ngstd
{


  /* ******************* TextOutArchive ******************* */
  

  TextOutArchive :: TextOutArchive (string filename) 
    : fout(filename.c_str()) 
  { ; }

  
  bool TextOutArchive :: Output () { return true; }
  bool TextOutArchive :: Input () { return false; }
  
  Archive & TextOutArchive :: operator & (double & d) 
  {
    fout << d << '\n';
    return *this;
  }
  
  Archive & TextOutArchive :: operator & (int & i)
  {
    fout << i << '\n';
    return *this;
  }

  Archive & TextOutArchive :: operator & (short & i)
  {
    fout << i << '\n';
    return *this;
  }
  
  Archive & TextOutArchive :: operator & (unsigned char & i)
  {
    fout << int (i) << '\n';
    return *this;
  }
  
  Archive & TextOutArchive :: operator & (bool & b) 
  {
    fout << (b ? 't' : 'f') << '\n';
    return *this;
  }
  
  Archive & TextOutArchive :: operator & (string & str)
  {
    int len = str.length();
    fout << len << '\n';
    fout.write (&str[0], len);
    fout << '\n';
    return *this;
  }

  Archive & TextOutArchive :: operator & (char *& str)
  {
    int len = strlen (str);
    fout << len << '\n';
    fout.write (&str[0], len);
    fout << '\n';
    return *this;
  }





  /* ******************* TextInArchive ******************* */


  TextInArchive :: TextInArchive (string filename) 
    : fin(filename.c_str()) 
  { ; }

  bool TextInArchive :: Output () { return false; }
  bool TextInArchive :: Input () { return true; }


  Archive & TextInArchive :: operator & (double & d) 
  {
    fin >> d;
    return *this;
  }
  
  Archive & TextInArchive :: operator & (int & i)
  {
    fin >> i;
    return *this;
  }

  Archive & TextInArchive :: operator & (short & i)
  {
    fin >> i;
    return *this;
  }
  
  Archive & TextInArchive :: operator & (unsigned char & i)
  {
    int hi;
    fin >> hi;
    i = hi;
    return *this;
  }
  
  Archive & TextInArchive :: operator & (bool & b) 
  {
    char c;
    fin >> c;
    b = (c == 't');
    return *this;
  }
  
  Archive & TextInArchive :: operator & (string & str) 
  {
    int len;
    fin >> len;
    
    char ch;
    fin.get(ch); // '\n'
    
    str.resize(len);
    fin.get(&str[0], len+1, '\0');
    return *this;
  }

  Archive & TextInArchive :: operator & (char *& str) 
  {
    int len;
    fin >> len;
    
    char ch;
    fin.get(ch); // '\n'
    
    str = new char[len+1];
    fin.get(&str[0], len, '\0');
    str[len] = 0;
    return *this;
  }

  
}
