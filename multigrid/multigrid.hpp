#ifndef FILE_MULTIGRID
#define FILE_MULTIGRID

/*********************************************************************/
/* File:   multigrid.hh                                              */
/* Author: Joachim Schoeberl                                         */
/* Date:   20. Apr. 2000                                             */
/*********************************************************************/

/* 
   Multigrid algorithms
*/


// #include <comp.hpp>

/// namespace for multigrid components

/*
namespace ngfem { };
namespace ngcomp { };

namespace ngmg
{
  using namespace std;
  using namespace ngstd;
  using ngcore::INT;
  using namespace ngla;
  using namespace ngfem;
  using namespace ngcomp;
}
*/


#include "mgpre.hpp"
#include "prolongation.hpp"
#include "smoother.hpp"
  // #include "vefc.hpp"
  // #include "evcoarse.hpp"


#endif
