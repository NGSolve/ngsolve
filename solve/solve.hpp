#ifndef FILE_SOLVE
#define FILE_SOLVE

/*********************************************************************/
/* File:   solve.hh                                                  */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/

/* 
   NGS Solvers: BVP, Instationary, ...
*/

#include <comp.hpp>
#include <multigrid.hpp>


#include <tcl.h>
#if TCL_MAJOR_VERSION==8 && TCL_MINOR_VERSION>=4
#define TCL_CONST_IS_CONST
#define tcl_const const
#else
#define tcl_const
#endif


/**
   A collection of solvers.
   Numerical procedures NumProc, e.g. for
   the solution of boundary value problems NumProcBVP or
   eigen value problems NumProcEVP, and many many more.
 */
namespace ngsolve
{
  using namespace std;

  using namespace ngstd;
  using namespace ngla;
  using namespace ngfem;
  using namespace ngcomp;
  using namespace ngmg;

// #ifdef SOCKETS
// #include "../markus/clientsocketaccess.hpp"
// #endif
#include "numproc.hpp"
#include "pde.hpp"


#ifdef WIN32
  const char dirslash = '\\';
#else
  const char dirslash = '/';
#endif

}


#endif
