#ifndef NGS_ARCHIVE_BASE
#define NGS_ARCHIVE_BASE

// copied from netgen

#include <vector>
#include <map>

namespace ngstd
{

  class Archive
  {
    bool is_output;
  public:
    Archive (bool ais_output) : is_output(ais_output) { ; }
    virtual ~Archive() { ; }

    bool Output () { return is_output; }
    bool Input () { return !is_output; }

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


    // archive a pointer ...

    int cnt = 0;
    std::map<void*,int> ptr2nr;
    std::vector<void*> nr2ptr;
    
    template <typename T>
    Archive & operator& (T *& p)
    {
      if (Output())
        {
          if (!p)
            {
              int m2 = -2;
              (*this) & m2;
              return *this;
            }
          auto pos = ptr2nr.find( (void*) p);
          if (pos == ptr2nr.end())
            {
              ptr2nr[p] = cnt;
              int m1 = -1;
              (*this) & m1;
              cnt++;
              (*this) & (*p);
            }
          else
            {
              (*this) & pos->second;
            }
        }
      else
        {
          int nr;
          (*this) & nr;
          cout << "in, got nr " << nr << endl;
          if (nr == -2)
            {
              p = nullptr;
            }
          else if (nr == -1)
            {
              p = new T;
              cout << "create new ptr, p = " << p << endl;
              (*this) & *p;
              nr2ptr.push_back(p);
            }
          else
            {
              p = (T*)nr2ptr[nr];
              cout << "reuse ptr " << nr << ": " << p << endl;
            }
        }
      return *this;
    }


    
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
