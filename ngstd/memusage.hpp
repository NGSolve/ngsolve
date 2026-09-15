#ifndef FILE_MEMUSAGE
#define FILE_MEMUSAGE

/**************************************************************************/
/* File:   memusage.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   16. June 2002                                                  */
/**************************************************************************/

namespace ngstd
{

  struct MemorySize
  {
    size_t host = 0, device = 0;

    MemorySize () = default;
    template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    MemorySize (T ahost) : host(size_t(ahost)) { }
    MemorySize (size_t ahost, size_t adevice) : host(ahost), device(adevice) { }
    static MemorySize Device (size_t n) { return { size_t(0), n }; }

    MemorySize & operator+= (MemorySize b) { host += b.host; device += b.device; return *this; }
    friend MemorySize operator+ (MemorySize a, MemorySize b) { return a += b; }
  };

/**
   Reports amount of used memory
 */
class MemoryUsage
{
protected:
  string name;
  MemorySize nbytes;
  size_t nblocks;
  const void * owner = nullptr;
public:
  MemoryUsage () = default;
  MemoryUsage (const string & aname, MemorySize anbytes, const void * aowner = nullptr)
    : name(aname), nbytes(anbytes), nblocks(0), owner(aowner)
  { ; }
  MemoryUsage (const string & aname,
               MemorySize anbytes, size_t anblocks, const void * aowner = nullptr)
    : name(aname), nbytes(anbytes), nblocks(anblocks), owner(aowner)
  { ; }
  MemoryUsage (const MemoryUsage &) = default;
  MemoryUsage (MemoryUsage &&) = default;
  MemoryUsage & operator= (const MemoryUsage &) = default;
  MemoryUsage & operator= (MemoryUsage &&) = default;
  
  void AddName (const string & aname) { name += aname; }
  const string & Name() const { return name; }
  size_t NBytes () const { return nbytes.host; }
  size_t NDeviceBytes () const { return nbytes.device; }
  MemorySize Bytes () const { return nbytes; }
  size_t NBlocks () const { return nblocks; }
  const void * Owner () const { return owner; }
};

}

#endif
