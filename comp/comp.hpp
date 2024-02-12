#ifndef FILE_COMP
#define FILE_COMP

/*********************************************************************/
/* File:   comp.hpp                                                  */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

#include <fem.hpp>
#include <finiteelement.hpp>
#include <la.hpp>

/**
   namespace for NGS-components.
   It contains the access to the mesh topology MeshAccess, 
   finite element spaces derived from FESpace, and global objects such as
   Gridfunction, LinearForm, and BilinearForm, or a Preconditioner.
   Computational functions such as SetValues and CalcFlux.
 */
namespace ngcomp
{
  using namespace std;
  using namespace ngstd;

  using ngfem::ELEMENT_TYPE;

  using namespace ngla;
  using namespace ngfem;
}


namespace ngmg
{
  class Prolongation;
  class TwoLevelMatrix;
  class MultigridPreconditioner;
}


#include <parallelngs.hpp>
namespace ngcomp
{
  using namespace ngparallel;
}


// #include "pmltrafo.hpp"
#include "meshaccess.hpp"
#include "ngsobject.hpp"
#include "fespace.hpp"

#include "gridfunction.hpp"
#include "bilinearform.hpp"
#include "linearform.hpp"
#include "preconditioner.hpp"

#include "postproc.hpp"
#include "interpolate.hpp"

#include "tpfes.hpp"
#include "hcurlhdivfes.hpp"
#include "hdivfes.hpp"
#include "h1hofespace.hpp"
#include "l2hofespace.hpp"
#include "hdivhofespace.hpp"
#include "hdivhosurfacefespace.hpp" 
#include "hcurlhofespace.hpp" 
// #include "facetfespace.hpp" 
// #include "vectorfacetfespace.hpp"
// #include "h1lumping.hpp"
#include "periodic.hpp"
#include "discontinuous.hpp"
#include "hidden.hpp"
#include "reorderedfespace.hpp"

#include "facetsurffespace.hpp"
#include "normalfacetsurfacefespace.hpp"
#include "fesconvert.hpp"

// #include "bddc.hpp"

#endif
