/*********************************************************************/
/* File:   templates.cpp                                             */
/* Author: Joachim Schoeberl                                         */
/* Date:   25. Mar. 2000                                             */
/*********************************************************************/


#include <ngstd.hpp>
#include <ngsolve_version.hpp>

namespace ngstd
{
  MPI_Comm ngs_comm;
  MPI_Comm onlyme_comm;

  ostream * testout = &cout;
  int printmessage_importance = 5;
  bool NGSOStream :: glob_active = true;
  const string ngsolve_version = NGSOLVE_VERSION;

  template class Array<int>;
}

