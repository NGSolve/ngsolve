#ifndef NGS_ARCHIVE_BASE
#define NGS_ARCHIVE_BASE

namespace ngstd
{

  class Archive
  {
  public:
    virtual bool Output () = 0;
    virtual bool Input () { return !Output(); }

    virtual Archive & operator & (double & d) = 0;
    virtual Archive & operator & (int & i) = 0;
    virtual Archive & operator & (short & i) = 0;
    virtual Archive & operator & (unsigned char & i) = 0;
    virtual Archive & operator & (bool & b) = 0;
    virtual Archive & operator & (string & str) = 0;
    virtual Archive & operator & (char *& str) = 0;

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
