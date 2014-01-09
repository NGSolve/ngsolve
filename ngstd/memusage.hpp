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
class MemoryUsageStruct
{
protected:
  string name;
  int nbytes;
  int nblocks;
public:
  MemoryUsageStruct (const string & aname,
		     int anbytes,
		     int anblocks)
    : name(aname), nbytes(anbytes), nblocks(anblocks)
  { ; }

  void AddName (const string & aname) { name += aname; }
  const string & Name() const { return name; }
  int NBytes () const { return nbytes; }
  int NBlocks () const { return nblocks; }
};

}

#endif
