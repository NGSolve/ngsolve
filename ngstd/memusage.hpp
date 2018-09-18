#ifndef FILE_MEMUSAGE
#define FILE_MEMUSAGE

/**************************************************************************/
/* File:   memusage.hpp                                                   */
/* Author: Joachim Schoeberl                                              */
/* Date:   16. June 2002                                                  */
/**************************************************************************/

namespace ngstd
{

/**
   Reports amount of used memory
 */
class MemoryUsage
{
protected:
  string name;
  size_t nbytes;
  size_t nblocks;
public:
  MemoryUsage () = default;
  MemoryUsage (const string & aname,
               size_t anbytes, size_t anblocks)
    : name(aname), nbytes(anbytes), nblocks(anblocks)
  { ; }
  MemoryUsage (const MemoryUsage &) = default;
  MemoryUsage (MemoryUsage &&) = default;
  MemoryUsage & operator= (const MemoryUsage &) = default;
  MemoryUsage & operator= (MemoryUsage &&) = default;
  
  void AddName (const string & aname) { name += aname; }
  const string & Name() const { return name; }
  size_t NBytes () const { return nbytes; }
  size_t NBlocks () const { return nblocks; }
};

}

#endif
