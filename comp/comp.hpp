#ifndef FILE_COMP
#define FILE_COMP

/*********************************************************************/
/* File:   comp.hpp                                                  */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   NGS Components: Mesh, Bilinearform, ....
*/

#include <fem.hpp>
#include <la.hpp>


#include <soldata.hpp>


namespace ngmg
{
  class Prolongation;
  class TwoLevelMatrix;
  class MultigridPreconditioner;
}




namespace ngparallel
{
  class ParallelDofs;
}

/// namespace for ngs-components
namespace ngcomp
{
  using namespace std;
  using namespace ngstd;
  using ngstd::INT;
  using ngfem::ELEMENT_TYPE;

  using namespace ngla;
  using namespace ngfem;
  using namespace ngparallel;
}


#include <parallelngs.hpp>



namespace ngcomp
{

  using namespace std;
  
  // important message
  class IM
  {
    int value;
  public:
    IM (int val) : value(val) { ; }
    int Value () const { return value; }
  };

  
  class NGSOStream
  {
    ostream & ost;
    bool active;
  public:
    NGSOStream (ostream & aost, bool aactive)
      : ost(aost), active(aactive) { ; }
    bool Active () const { return active; }
    ostream & GetStream () { return ost; }
  };
  
  inline NGSOStream operator<< (ostream & ost, const IM & im)
  {
    return NGSOStream (ost, 
		       (im.Value() <= printmessage_importance)  &&
		       (id == 0) );
  }


  
  template <typename T>
  inline NGSOStream operator<< (NGSOStream ngsost, const T & id)
  {
    if (ngsost.Active())
      ngsost.GetStream() << id;
    return ngsost;
  }
  
  inline NGSOStream operator<< (NGSOStream ngsost, ostream& ( *pf )(ostream&))
  {
    if ( ngsost.Active() )
      ngsost.GetStream() << (*pf);

    return ngsost;
  }
  
  inline NGSOStream operator<< (NGSOStream ngsost, ios& ( *pf )(ios&))
  {
    if ( ngsost.Active() )
      ngsost.GetStream() << (*pf);
    
    return ngsost;
  }

  inline NGSOStream operator<< (NGSOStream ngsost, ios_base& ( *pf )(ios_base&))
  {
    if ( ngsost.Active() )
      ngsost.GetStream() << (*pf);
    
    return ngsost;
  }
}




#include "meshaccess.hpp"
#include "ngsobject.hpp"
#include "fespace.hpp"
#include "hcurlhdivfes.hpp"
#include "hdivfes.hpp"
#include "h1hofespace.hpp"
#include "l2hofespace.hpp"
#include "gridfunction.hpp"
#include "bilinearform.hpp"
#include "linearform.hpp"
#include "postproc.hpp"
#include "hdivhofespace.hpp" 
#include "hcurlhofespace.hpp" 
#include "facetfespace.hpp" 
#include "vectorfacetfespace.hpp"

#include "preconditioner.hpp"
#include "bddc.hpp"

#endif
