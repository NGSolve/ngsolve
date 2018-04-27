#ifndef FILE_ARCHIVE
#define FILE_ARCHIVE


#include "archive_base.hpp" 


namespace ngstd
{
  class TextOutArchive : public Archive
  {
    shared_ptr<ostream> fout;
  public:
    TextOutArchive (string filename);
    TextOutArchive (shared_ptr<ostream> afout) : Archive(true), fout(afout) { ; }

    // virtual bool Output ();
    // virtual bool Input ();

    virtual Archive & operator & (double & d);
    virtual Archive & operator & (int & i);
    virtual Archive & operator & (short & i);
    virtual Archive & operator & (long & i);
    virtual Archive & operator & (size_t & i);
    virtual Archive & operator & (unsigned char & i);
    virtual Archive & operator & (bool & b);
    virtual Archive & operator & (string & str);
    virtual Archive & operator & (char *& str);
  };


  class TextInArchive : public Archive
  {
    shared_ptr<istream> fin;
  public:
    TextInArchive (string filename);
    TextInArchive (shared_ptr<istream> afin) : Archive(false), fin(afin) { ; }

    // virtual bool Output ();
    // virtual bool Input ();

    virtual Archive & operator & (double & d);
    virtual Archive & operator & (int & i);
    virtual Archive & operator & (short & i);
    virtual Archive & operator & (long & i);
    virtual Archive & operator & (size_t & i);
    virtual Archive & operator & (unsigned char & i);
    virtual Archive & operator & (bool & b);
    virtual Archive & operator & (string & str);
    virtual Archive & operator & (char *& str);
  };


  class BinaryOutArchive : public Archive
  {
    shared_ptr<ostream> fout;
    size_t ptr = 0;
    enum { BUFFERSIZE = 1024 };
    char buffer[BUFFERSIZE];
  public:
    BinaryOutArchive (string filename);
    BinaryOutArchive (shared_ptr<ostream> afout) : Archive(true), fout(afout) { ; }
    virtual ~BinaryOutArchive () { FlushBuffer(); }
    // virtual bool Output ();
    // virtual bool Input ();

    virtual Archive & operator & (double & d);
    virtual Archive & operator & (int & i);
    virtual Archive & operator & (short & i);
    virtual Archive & operator & (long & i);
    virtual Archive & operator & (size_t & i);
    virtual Archive & operator & (unsigned char & i);
    virtual Archive & operator & (bool & b);
    virtual Archive & operator & (string & str);
    virtual Archive & operator & (char *& str);

    template <typename T>
    Archive & Write (T x);
    void FlushBuffer();
  };


  class BinaryInArchive : public Archive
  {
    shared_ptr<istream> fin;
  public:
    BinaryInArchive (string filename);
    BinaryInArchive (shared_ptr<istream> afin) : Archive(false), fin(afin) { ; }

    // virtual bool Output ();
    // virtual bool Input ();

    virtual Archive & operator & (double & d);
    virtual Archive & operator & (int & i);
    virtual Archive & operator & (short & i);
    virtual Archive & operator & (long & i);
    virtual Archive & operator & (size_t & i);
    virtual Archive & operator & (unsigned char & i);
    virtual Archive & operator & (bool & b);
    virtual Archive & operator & (string & str);
    virtual Archive & operator & (char *& str);

    virtual Archive & Do (double * d, size_t n); 
    virtual Archive & Do (int * i, size_t n); 
    virtual Archive & Do (size_t * i, size_t n); 

  };










  // archive a pointer ...
  template <typename T>
  Archive & operator& (Archive & ar, T *& p)
  {
    if (ar.Output())
      ar & (*p);
    else
      {
        p = new T;
        ar & *p;
      }
    return ar;
  }

  // archive a shared pointer ...
  template <typename T>
  Archive & operator& (Archive & ar, shared_ptr<T> & p)
  {
    if (ar.Output())
      ar & (*p);
    else
      {
        p = make_shared<T>();
        ar & *p;
      }
    return ar;
  }



  // cannot archive type T
  template <typename T>
  Archive & operator& (Archive & ar, T & t)
  {
    if (ar.Output())
      {
        cout << string("Cannot_archive_object_of_type ") << string(typeid(T).name()) << endl;
        ar << string("Cannot_archive_object_of_type ") << string(typeid(T).name());
      }
    else
      {
        string dummy;
        ar & dummy;
        ar & dummy;
      }
    return ar;
  }
}




#endif
