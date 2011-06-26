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
