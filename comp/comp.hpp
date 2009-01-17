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


#include <nginterface.h>
namespace netgen {
#include <soldata.hpp>
}

namespace ngmg
{
  class Prolongation;
  class TwoLevelMatrix;
  class MultigridPreconditioner;
}


namespace ngsolve
{
  class PDE;
}
using ngsolve::PDE;

namespace RS_AMG
{
  class OutBack;
}


#ifdef PARALLEL
namespace ngparallel
{
  class ParallelDofs;
  class ParallelMeshAccess;
}
#endif


/// namespace for ngs-components
namespace ngcomp
{
  using namespace std;
  using namespace ngstd;
  using ngstd::INT;

  using namespace ngla;
  using namespace ngfem;

  class MeshAccess;

#include "ngsobject.hpp"
#include "statushandler.hpp"
#include "meshaccess.hpp"
#include "fespace.hpp"
#include "hcurlhdivfes.hpp"
#include "hdivfes.hpp"
#include "h1hofespace.hpp"
#include "l2hofespace.hpp"
#include "gridfunction.hpp"
#include "bilinearform.hpp"
#include "linearform.hpp"
#include "preconditioner.hpp"
#include "postproc.hpp"
#include "hdivhofespace.hpp" 
#include "hdivhybridhofespace.hpp" 
#include "hcurlhofespace.hpp" 
#include "facetfespace.hpp" 
#include "vectorfacetfespace.hpp"
}


#endif
