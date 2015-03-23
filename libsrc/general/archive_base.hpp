#ifndef NGS_ARCHIVE_BASE
#define NGS_ARCHIVE_BASE

// copied from netgen

namespace ngstd
{

  class Archive
  {
  public:
    virtual bool Output () = 0;
    virtual bool Input () { return !Output(); }

    virtual Archive & operator & (double & d) = 0;
    virtual Archive & operator & (int & i) = 0;
    virtual Archive & operator & (long & i) = 0;
    virtual Archive & operator & (size_t & i) = 0;
    virtual Archive & operator & (short & i) = 0;
    virtual Archive & operator & (unsigned char & i) = 0;
    virtual Archive & operator & (bool & b) = 0;
    virtual Archive & operator & (string & str) = 0;
    virtual Archive & operator & (char *& str) = 0;

    template <typename T>
    Archive & Do (T * data, size_t n) 
    { for (size_t j = 0; j < n; j++) { (*this) & data[j]; }; return *this; };


    virtual Archive & Do (double * d, size_t n) 
    { for (size_t j = 0; j < n; j++) { (*this) & d[j]; }; return *this; };

    virtual Archive & Do (int * i, size_t n) 
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; };

    virtual Archive & Do (long * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; };

    virtual Archive & Do (size_t * i, size_t n) 
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; };

    virtual Archive & Do (short * i, size_t n) 
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; };

    virtual Archive & Do (unsigned char * i, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & i[j]; }; return *this; };

    virtual Archive & Do (bool * b, size_t n)
    { for (size_t j = 0; j < n; j++) { (*this) & b[j]; }; return *this; };


    // nvirtual Archive & Do (string * str, size_t n)
    // { for (size_t j = 0; j < n; j++) { (*this) & str[j]; }; return *this; };
    // virtual Archive & operator & (char *& str) = 0;

    template <typename T>
    Archive & operator << (const T & t)
    {
      T ht(t);
      (*this) & ht;
      return *this;
    }
  };


}


#endif
