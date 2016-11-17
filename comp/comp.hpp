#ifndef FILE_COMP
#define FILE_COMP

/*********************************************************************/
/* File:   comp.hpp                                                  */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

#include <fem.hpp>
#include <la.hpp>

#include <soldata.hpp>   // netgen visualization


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

  using ngstd::INT;
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



#include "meshaccess.hpp"
#include "ngsobject.hpp"
#include "fespace.hpp"

#include "gridfunction.hpp"
#include "bilinearform.hpp"
#include "linearform.hpp"
#include "preconditioner.hpp"
#include "numproc.hpp"
#include "pde.hpp"

#include "postproc.hpp"

#include "tpfes.hpp"
#include "hcurlhdivfes.hpp"
#include "hdivfes.hpp"
#include "h1hofespace.hpp"
#include "l2hofespace.hpp"
#include "hdivhofespace.hpp" 
#include "hcurlhofespace.hpp" 
#include "facetfespace.hpp" 
#include "vectorfacetfespace.hpp"

// #include "bddc.hpp"
#include "hypre_precond.hpp"
#include "vtkoutput.hpp"

#endif
