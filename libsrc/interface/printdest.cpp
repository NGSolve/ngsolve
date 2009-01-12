#include <mystdlib.h>
#include <myadt.hpp>

namespace netgen
{
  //Destination for messages, errors, ...
  void Ng_PrintDest(const char * s)
  {
    (*mycout) << s << flush;
  }
}
