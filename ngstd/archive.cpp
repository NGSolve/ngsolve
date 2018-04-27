/**************************************************************************/
/* File:   archive.cpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   14. Feb. 2014                                                  */
/**************************************************************************/

#include <ngstd.hpp>


namespace ngstd
{

  /* ******************* BinaryOutArchive ******************* */
  

  BinaryOutArchive :: BinaryOutArchive (string filename) 
    : BinaryOutArchive(make_shared<ofstream>(filename.c_str()))
  { ; }

  Archive & BinaryOutArchive :: operator & (double & d) 
  {
    return Write(d);
    /*
    fout->write(reinterpret_cast<char*>(&d), sizeof(double));
    return *this;
    */
  }
  
  Archive & BinaryOutArchive :: operator & (int & i)
  {
    return Write(i);
    /*
    fout->write(reinterpret_cast<char*>(&i), sizeof(int));
    return *this;
    */
  }

  Archive & BinaryOutArchive :: operator & (short & i)
  {
    return Write(i);
    // fout->write(reinterpret_cast<char*>(&i), sizeof(short));
    // return *this;
  }
  
  Archive & BinaryOutArchive :: operator & (long & i)
  {
    return Write(i);
    // fout->write(reinterpret_cast<char*>(&i), sizeof(long));
    // return *this;
  }

  Archive & BinaryOutArchive :: operator & (size_t & i)
  {
    return Write(i);
    // fout->write(reinterpret_cast<char*>(&i), sizeof(size_t));
    // return *this;
  }

  Archive & BinaryOutArchive :: operator & (unsigned char & i)
  {
    return Write(i);
    // fout->write(reinterpret_cast<char *>(&i), sizeof(unsigned char));
    // return *this;
  }
  
  Archive & BinaryOutArchive :: operator & (bool & b) 
  {
    return Write(b);
    // fout->write(reinterpret_cast<char*>(&b), sizeof(bool));
    // return *this;
  }
  
  Archive & BinaryOutArchive :: operator & (string & str)
  {
    if (ptr > 0) FlushBuffer();    

    int len = str.length();
    fout->write (reinterpret_cast<char*>(&len), sizeof(int));
    fout->write (&str[0], len);
    return *this;
  }

  Archive & BinaryOutArchive :: operator & (char *& str)
  {
    if (ptr > 0) FlushBuffer();

    int len = strlen (str);
    fout->write (reinterpret_cast<char*>(&len), sizeof(int));
    fout->write (&str[0], len);
    return *this;
  }

  
  template <typename T>
  Archive & BinaryOutArchive :: Write (T x)
  {
    if (unlikely(ptr > BUFFERSIZE-sizeof(T)))
      {
        fout->write(&buffer[0], ptr);
        * (T*) (&buffer[0]) = x;        
        ptr = sizeof(T);
        return *this;
      }
    * (T*) (&buffer[ptr]) = x;
    ptr += sizeof(T);
    return *this;
    /*
    fout->write(reinterpret_cast<char*>(&x), sizeof(T));
    return *this;
    */
  }

  void BinaryOutArchive :: FlushBuffer()
  {
    if (ptr > 0)
      {
        fout->write(&buffer[0], ptr);
        ptr = 0;
      }
  }



  /* ******************* BinaryInArchive ******************* */


  BinaryInArchive :: BinaryInArchive (string filename) 
    : BinaryInArchive(make_shared<ifstream>(filename.c_str())) { ; }
  // { fin = make_shared<ifstream>(filename.c_str()); }

  Archive & BinaryInArchive :: operator & (double & d) 
  {
    fin->read(reinterpret_cast<char*>(&d), sizeof(double));
    return *this;
  }
  
  Archive & BinaryInArchive :: operator & (int & i)
  {
    fin->read(reinterpret_cast<char*>(&i), sizeof(int));
    return *this;
  }

  Archive & BinaryInArchive :: operator & (short & i)
  {
    fin->read(reinterpret_cast<char*>(&i), sizeof(short));
    return *this;
  }

  Archive & BinaryInArchive :: operator & (long & i)
  {
    fin->read(reinterpret_cast<char*>(&i), sizeof(long));
    return *this;
  }

  Archive & BinaryInArchive :: operator & (size_t & i)
  {
    fin->read(reinterpret_cast<char*>(&i), sizeof(size_t));
    return *this;
  }
  
  Archive & BinaryInArchive :: operator & (unsigned char & i)
  {
    fin->read(reinterpret_cast<char *>(&i), sizeof(unsigned char));
    return *this;
  }
  
  Archive & BinaryInArchive :: operator & (bool & b) 
  {
    fin->read(reinterpret_cast<char*>(&b), sizeof(bool));
    return *this;
  }
  
  Archive & BinaryInArchive :: operator & (string & str) 
  {
    int len;
    fin->read(reinterpret_cast<char*>(&len), sizeof(int));
    str.resize(len);
    fin->read(&str[0], len);
    return *this;
  }

  Archive & BinaryInArchive :: operator & (char *& str) 
  {
    int len;
    fin->read(reinterpret_cast<char*>(&len), sizeof(int));
    str = new char[len+1];
    fin->read(&str[0], len);
    str[len] = '\0';
    return *this;
  }


  Archive & BinaryInArchive :: Do (double * d, size_t n) 
  {
    cout << "load double array, size = " << n << endl;
    fin->read(reinterpret_cast<char*>(d), n*sizeof(double));
    return *this;
  }

  Archive & BinaryInArchive :: Do (int * i, size_t n) 
  {
    cout << "load int array, size = " << n << endl;
    fin->read(reinterpret_cast<char*>(i), n*sizeof(int));
    return *this;
  }

  Archive & BinaryInArchive :: Do (size_t * i, size_t n) 
  {
    cout << "load size_t array, size = " << n << endl;
    fin->read(reinterpret_cast<char*>(i), n*sizeof(size_t));
    return *this;
  }








  

  /* ******************* TextOutArchive ******************* */
  

  TextOutArchive :: TextOutArchive (string filename) 
    : TextOutArchive(make_shared<ofstream>(filename.c_str())) { ; }

  
  Archive & TextOutArchive :: operator & (double & d) 
  {
    *fout << d << '\n';
    return *this;
  }
  
  Archive & TextOutArchive :: operator & (int & i)
  {
    *fout << i << '\n';
    return *this;
  }

  Archive & TextOutArchive :: operator & (short & i)
  {
    *fout << i << '\n';
    return *this;
  }
  
  Archive & TextOutArchive :: operator & (long & i)
  {
    *fout << i << '\n';
    return *this;
  }

  Archive & TextOutArchive :: operator & (size_t & i)
  {
    *fout << i << '\n';
    return *this;
  }

  Archive & TextOutArchive :: operator & (unsigned char & i)
  {
    *fout << int (i) << '\n';
    return *this;
  }
  
  Archive & TextOutArchive :: operator & (bool & b) 
  {
    *fout << (b ? 't' : 'f') << '\n';
    return *this;
  }
  
  Archive & TextOutArchive :: operator & (string & str)
  {
    int len = str.length();
    *fout << len << '\n';
    fout->write (&str[0], len);
    *fout << '\n';
    return *this;
  }

  Archive & TextOutArchive :: operator & (char *& str)
  {
    int len = strlen (str);
    *fout << len << '\n';
    fout->write (&str[0], len);
    *fout << '\n';
    return *this;
  }





  /* ******************* TextInArchive ******************* */


  TextInArchive :: TextInArchive (string filename) 
    : TextInArchive(make_shared<ifstream>(filename.c_str())) { ; }


  Archive & TextInArchive :: operator & (double & d) 
  {
    *fin >> d;
    return *this;
  }
  
  Archive & TextInArchive :: operator & (int & i)
  {
    *fin >> i;
    return *this;
  }

  Archive & TextInArchive :: operator & (short & i)
  {
    *fin >> i;
    return *this;
  }

  Archive & TextInArchive :: operator & (long & i)
  {
    *fin >> i;
    return *this;
  }

  Archive & TextInArchive :: operator & (size_t & i)
  {
    *fin >> i;
    return *this;
  }
  
  Archive & TextInArchive :: operator & (unsigned char & i)
  {
    int hi;
    *fin >> hi;
    i = hi;
    return *this;
  }
  
  Archive & TextInArchive :: operator & (bool & b) 
  {
    char c;
    *fin >> c;
    b = (c == 't');
    return *this;
  }
  
  Archive & TextInArchive :: operator & (string & str) 
  {
    int len;
    *fin >> len;
    
    char ch;
    fin->get(ch); // '\n'
    
    str.resize(len);
    fin->get(&str[0], len+1, '\0');
    return *this;
  }

  Archive & TextInArchive :: operator & (char *& str) 
  {
    int len;
    *fin >> len;
    
    char ch;
    fin->get(ch); // '\n'
    
    str = new char[len+1];
    fin->get(&str[0], len, '\0');
    str[len] = 0;
    return *this;
  }

  
}
