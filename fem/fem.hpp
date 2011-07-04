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

/*
#ifdef NETGEN_ELTRANS
#include <nginterface.h>
#include <nginterface_v2.hpp>

// using namespace netgen;
#endif
*/


/**
   namespace for finite elements.
   
   Definition of reference FiniteElement, ScalarFiniteElement, and vector-valued elements HDivFiniteElement and HCurlFiniteElement 

   Definition of the geometry of the element, i.e., ElementTransformation

   IntegrationPoint on the reference element, and SpecificIntegrationPoint on the mapped element

   Element-matrix and element-vector calculation by BilinearFormIntegrator and LinearFormIntegrator
*/
namespace ngfem
{
  using namespace std;
  using namespace ngstd;
  using ngstd::INT;
  using namespace ngbla;
}


#include "elementtopology.hpp"
#include "intrule.hpp"

#include "generic_recpol.hpp"
#include "recursive_pol.hpp"
#include "recursive_pol_trig.hpp"
#include "recursive_pol_tet.hpp"

// #include "precomp.hpp"

#include "finiteelement.hpp"
#include "scalarfe.hpp"

#include "elementtransformation.hpp"
#include "coefficient.hpp"


#include "h1lofe.hpp"
#include "h1hofe.hpp"
#include "l2hofe.hpp"

#include "hdivfe.hpp"
#include "hcurlfe.hpp"

#include "hdivhofe.hpp"
#include "hcurlhofe.hpp" 


#include "specialelement.hpp"

#include "facetfe.hpp" 
#include "vectorfacetfe.hpp"


#include "integrator.hpp"
#include "diffop.hpp"
#include "bdbintegrator.hpp"
#include "bdbequations.hpp"
#include "hcurl_equations.hpp"
#include "hdiv_equations.hpp"
#include "elasticity_equations.hpp"

// #include "pml.hpp" 

using ngfem::ELEMENT_TYPE;

#endif
