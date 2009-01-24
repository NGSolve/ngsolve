#ifndef FILE_FEM
#define FILE_FEM

/*********************************************************************/
/* File:   fem.hpp                                                   */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   Finite Element kernel functions
*/

#include <bla.hpp>

#ifdef NETGEN_ELTRANS
#include <nginterface.h>
#endif

/**
   namespace for finite elements.
   
   Definition of reference FiniteElement, mixed fe HDivFiniteElement and HCurlFiniteElementD 
   
   Element-matrix and element-vector assembling BilinearFormIntegrator, LinearFormIntegrator
*/
namespace ngfem
{
  using namespace std;
  using namespace ngstd;
  using ngstd::INT;
  using namespace ngbla;


#include "elementtopology.hpp"
#include "intrule.hpp"

#include "generic_recpol.hpp"
#include "recursive_pol.hpp"
#include "recursive_pol_trig.hpp"
#include "recursive_pol_tet.hpp"

#include "finiteelement.hpp"
#include "elementtransformation.hpp"
#include "coefficient.hpp"

#include "highorderfe.hpp"
#include "h1hofe.hpp"
#include "l2hofe.hpp"
#include "hdivfe.hpp"
#include "hcurlfe.hpp"
#include "specialelement.hpp"

#include "hdivhofe.hpp"
#include "hcurlhofe.hpp" 
#include "facetfe.hpp" 
#include "vectorfacetfe.hpp"


#include "integrator.hpp"
#include "bdbintegrator.hpp"
#include "bdbequations.hpp"
#include "hdiv_equations.hpp"

#include "pml.hpp" 
}


#endif
